"""One-shot end-to-end driver: cold start → scored variant DB.

Wraps the existing pipeline + recovery stack into a single command so a
new gene goes from "I have a symbol and an NCBI email" to "I have a
scored SQLite DB and a one-page report" without remembering which
script to run when.

Order of operations:

  1. doctor — validate .env keys + .venv tooling are reachable.
  2. extract — runs cli.automated_workflow.automated_variant_extraction_workflow
     (discovery → harvest → tier 1/2 filter → download → extract →
     migrate). Auto-resumes via GVF_RESUME_DIR when --resume-dir given.
  3. layers — runs scripts/recall_recovery/run_all_layers.py against
     the resulting DB. Auto-detects gold standard CSV under
     gene_variant_fetcher_gold_standard/normalized/ and the KCNH2 v12
     manual-recovery DB when --gene KCNH2.
  4. report — writes RUN_REPORT.md summarizing what happened, what
     scored, and what's gated on credentials.

Cold-start invocation (one command, no flags beyond gene + email):

    cli gvf-run --gene KCNQ1 --email me@uni.edu --output results/

Skip the LLM-bound extract step and just score an existing DB through
the recovery layers:

    cli gvf-run --gene KCNQ1 --email me@uni.edu --output results/ --skip extract
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[1]

logger = logging.getLogger("gvf_run")


# ---------------------------------------------------------------------------
# Step 1 — doctor
# ---------------------------------------------------------------------------


REQUIRED_ENV = ("NCBI_EMAIL",)
RECOMMENDED_ENV = (
    "NCBI_API_KEY",
    "OPENAI_API_KEY",  # or anthropic / azure_ai
    "ELSEVIER_API_KEY",
    "WILEY_API_KEY",
    "SPRINGER_API_KEY",
)
CREDENTIAL_UNLOCKS = (
    "ELSEVIER_INSTTOKEN",  # subscription full text
)


def doctor() -> dict:
    """Return a status dict describing the environment.

    Doesn't raise — just reports what's available so the caller can decide
    whether to proceed.
    """
    status: dict = {"required": {}, "recommended": {}, "unlocks": {}, "ok": True}
    for k in REQUIRED_ENV:
        present = bool(os.environ.get(k))
        status["required"][k] = present
        if not present:
            status["ok"] = False
    for k in RECOMMENDED_ENV + CREDENTIAL_UNLOCKS:
        bucket = "unlocks" if k in CREDENTIAL_UNLOCKS else "recommended"
        status[bucket][k] = bool(os.environ.get(k))
    # Quick reachability check — NCBI esearch
    try:
        import requests

        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi",
            params={
                "tool": "GVF",
                "email": os.environ.get("NCBI_EMAIL", "gvf@example.org"),
            },
            timeout=10,
        )
        status["ncbi_reachable"] = r.status_code == 200
    except Exception as e:
        status["ncbi_reachable"] = False
        status["ncbi_error"] = str(e)[:100]
    return status


# ---------------------------------------------------------------------------
# Step 2 — extract (wraps the existing workflow)
# ---------------------------------------------------------------------------


def step_extract(
    gene: str,
    email: str,
    output_dir: Path,
    pmid_file: Optional[Path],
    max_pmids: int,
    resume_dir: Optional[Path],
) -> Path:
    """Run the existing automated workflow. Returns the run dir it produced."""
    from cli.automated_workflow import automated_variant_extraction_workflow

    if resume_dir:
        os.environ["GVF_RESUME_DIR"] = str(resume_dir)
        logger.info("Resuming from %s", resume_dir)

    pmids = None
    if pmid_file:
        pmids = [
            line.split("#", 1)[0].strip()
            for line in pmid_file.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        pmids = [p for p in pmids if p.isdigit()]
        logger.info("Using %d PMIDs from %s (calibrated mode)", len(pmids), pmid_file)

    automated_variant_extraction_workflow(
        gene_symbol=gene,
        email=email,
        output_dir=str(output_dir),
        max_pmids=max_pmids,
        pmids=pmids,
        scout_first=True,
    )

    # Find the run dir the workflow just created (latest timestamped child).
    gene_root = output_dir / gene
    runs = sorted(
        (p for p in gene_root.iterdir() if p.is_dir()),
        key=lambda p: p.stat().st_mtime,
    )
    if not runs:
        raise RuntimeError(f"No run dir under {gene_root}")
    return runs[-1]


# ---------------------------------------------------------------------------
# Step 3 — layers
# ---------------------------------------------------------------------------


def _find_db(run_dir: Path, gene: str) -> Optional[Path]:
    db = run_dir / f"{gene}.db"
    if db.exists():
        return db
    # Fallback: any *.db file
    candidates = sorted(run_dir.glob("*.db"), key=lambda p: p.stat().st_mtime)
    return candidates[-1] if candidates else None


def _find_gold(gene: str) -> Optional[Path]:
    gold = (
        REPO_ROOT
        / "gene_variant_fetcher_gold_standard"
        / "normalized"
        / f"{gene}_recall_input.csv"
    )
    return gold if gold.exists() else None


def _find_v12_db(gene: str) -> Optional[Path]:
    """KCNH2 has a v12 manually-curated baseline DB on disk; other genes don't."""
    if gene.upper() != "KCNH2":
        return None
    candidates = list(
        REPO_ROOT.glob("results/KCNH2/*/end_to_end_*manual_recovery/KCNH2_v12*.db")
    )
    if not candidates:
        return None
    return sorted(candidates, key=lambda p: p.stat().st_mtime)[-1]


def step_layers(
    gene: str,
    db: Path,
    run_dir: Path,
    gold: Optional[Path],
    outdir: Path,
    with_v12: Optional[Path],
) -> Optional[Path]:
    """Run the recall-recovery driver. Returns the progression.csv path."""
    if not gold:
        logger.warning(
            "No gold standard CSV found for %s under gene_variant_fetcher_gold_standard/"
            "normalized/ — skipping recovery layers. The DB at %s is still usable; "
            "internal QC signals should be reviewed instead. See "
            "docs/NEW_GENE_RUNBOOK.md.",
            gene,
            db,
        )
        return None

    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "recall_recovery" / "run_all_layers.py"),
        "--gene",
        gene,
        "--db",
        str(db),
        "--gold",
        str(gold),
        "--pmc-dir",
        str(run_dir / "pmc_fulltext"),
        "--outdir",
        str(outdir),
    ]
    if with_v12:
        cmd.extend(["--with-v12", str(with_v12)])

    logger.info("layers → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("recovery layers failed: %s", result.stderr[-400:])
    return outdir / "progression.csv"


# ---------------------------------------------------------------------------
# Step 4 — report
# ---------------------------------------------------------------------------


def step_report(
    gene: str,
    db: Path,
    run_dir: Path,
    layer_outdir: Optional[Path],
    doctor_status: dict,
    started_at: float,
    duration_s: float,
    out_path: Path,
) -> None:
    """Write a RUN_REPORT.md summarizing the whole run."""
    lines: list[str] = []
    lines.append(f"# GVF Run Report — {gene}")
    lines.append("")
    lines.append(f"- Started: {datetime.fromtimestamp(started_at).isoformat()}")
    lines.append(f"- Duration: {duration_s/60:.1f} min")
    lines.append(f"- Run dir: `{run_dir}`")
    lines.append(f"- DB: `{db}`")
    lines.append("")
    lines.append("## Doctor")
    lines.append("")
    lines.append("Required:")
    for k, v in doctor_status.get("required", {}).items():
        lines.append(f"  - {k}: {'✓' if v else '✗'}")
    lines.append("Recommended:")
    for k, v in doctor_status.get("recommended", {}).items():
        lines.append(f"  - {k}: {'✓' if v else '–'}")
    lines.append("Subscription unlocks (not required, but lift recall):")
    for k, v in doctor_status.get("unlocks", {}).items():
        lines.append(f"  - {k}: {'✓' if v else '–'}")
    lines.append(f"- NCBI reachable: {doctor_status.get('ncbi_reachable')}")
    lines.append("")

    if layer_outdir and (layer_outdir / "progression.json").exists():
        progression = json.loads((layer_outdir / "progression.json").read_text())
        rows = progression.get("progression", [])
        lines.append("## Recall Progression (per layer)")
        lines.append("")
        lines.append(
            "| Layer | pmids | variant_rows | unique_variants | patients | affected | unaffected |"
        )
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for r in rows:
            lines.append(
                "| {layer} | {pmids} | {variant_rows} | {unique_variants} | {patients} | {affected} | {unaffected} |".format(
                    layer=r.get("layer", "?"),
                    pmids=_fmt(r.get("pmids")),
                    variant_rows=_fmt(r.get("variant_rows")),
                    unique_variants=_fmt(r.get("unique_variants")),
                    patients=_fmt(r.get("patients")),
                    affected=_fmt(r.get("affected")),
                    unaffected=_fmt(r.get("unaffected")),
                )
            )
        lines.append("")
        # Identify metrics still under 90%
        final = rows[-1] if rows else {}
        under = [
            k
            for k in (
                "pmids",
                "variant_rows",
                "unique_variants",
                "patients",
                "affected",
                "unaffected",
            )
            if final.get(k) is not None and final[k] < 90
        ]
        if under:
            lines.append("## Metrics under 90% target")
            for k in under:
                lines.append(f"- **{k}**: {final[k]:.1f}%")
            lines.append("")
        else:
            lines.append("## All six metrics ≥ 90%")
            lines.append("")
    else:
        lines.append("## Recall Progression")
        lines.append("")
        lines.append(
            "_Recovery layers skipped — no gold standard CSV found for this gene._"
        )
        lines.append(
            "Review internal QC signals instead (see docs/NEW_GENE_RUNBOOK.md)."
        )
        lines.append("")

    lines.append("## Next steps")
    if doctor_status.get("unlocks", {}).get("ELSEVIER_INSTTOKEN") is False:
        lines.append(
            "- ELSEVIER_INSTTOKEN is unset. Adding it (request from your library) "
            "typically lifts variant_rows by ~30-50 pp for cardiac channel genes "
            "because Heart Rhythm / JACC / Eur Heart J unlock."
        )
    if doctor_status.get("recommended", {}).get("OPENAI_API_KEY") is False:
        lines.append(
            "- OPENAI_API_KEY is unset. The extractor falls back to other "
            "providers when configured; without any LLM credential the "
            "extract step degrades to abstract-only."
        )
    if layer_outdir:
        lines.append(f"- Inspect per-layer outputs under `{layer_outdir}`.")
    lines.append(f"- The scored DB is at `{db}` for ad-hoc queries.")
    lines.append("")

    out_path.write_text("\n".join(lines))
    logger.info("Wrote %s", out_path)


def _fmt(v: Optional[float]) -> str:
    if v is None:
        return "—"
    return f"{v:.1f}%"


# ---------------------------------------------------------------------------
# Entry point (called from cli/__init__.py)
# ---------------------------------------------------------------------------


def run_gvf_pipeline(
    gene: str,
    email: str,
    output: Path,
    pmid_file: Optional[Path] = None,
    max_pmids: int = 1500,
    resume_dir: Optional[Path] = None,
    skip: Optional[list[str]] = None,
) -> int:
    """Execute the full pipeline. Returns exit code."""
    skip = {s.lower() for s in (skip or [])}
    started = time.time()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    output = Path(output).expanduser()
    output.mkdir(parents=True, exist_ok=True)

    # Step 1: doctor
    logger.info("🩺 Step 1: doctor")
    status = doctor()
    if not status["ok"] and "doctor" not in skip:
        logger.error(
            "Doctor failed required-env check; pass --skip doctor to override."
        )
        for k, v in status["required"].items():
            if not v:
                logger.error("  missing: %s", k)
        return 2

    # Step 2: extract (unless skipped)
    gene = gene.upper()
    run_dir: Optional[Path] = None
    if "extract" in skip:
        # Find existing run dir
        gene_root = output / gene
        if gene_root.is_dir():
            runs = sorted(
                (p for p in gene_root.iterdir() if p.is_dir()),
                key=lambda p: p.stat().st_mtime,
            )
            run_dir = runs[-1] if runs else None
        if not run_dir:
            logger.error(
                "--skip extract requires an existing run dir under %s", output / gene
            )
            return 2
        logger.info("⏭️  Step 2: extract — SKIPPED (using existing %s)", run_dir)
    else:
        logger.info("📚 Step 2: extract")
        try:
            run_dir = step_extract(
                gene=gene,
                email=email,
                output_dir=output,
                pmid_file=pmid_file,
                max_pmids=max_pmids,
                resume_dir=resume_dir,
            )
        except Exception as e:
            logger.exception("extract step failed: %s", e)
            return 3

    db = _find_db(run_dir, gene)
    if not db:
        logger.error("No DB produced in %s", run_dir)
        return 4

    # Step 3: layers
    layer_outdir: Optional[Path] = None
    if "layers" not in skip:
        logger.info("🧬 Step 3: recovery layers")
        gold = _find_gold(gene)
        v12 = _find_v12_db(gene)
        layer_outdir = run_dir / "layers"
        step_layers(
            gene=gene,
            db=db,
            run_dir=run_dir,
            gold=gold,
            outdir=layer_outdir,
            with_v12=v12,
        )
    else:
        logger.info("⏭️  Step 3: layers — SKIPPED")

    # Step 4: report
    if "report" not in skip:
        logger.info("📝 Step 4: report")
        report_path = run_dir / "RUN_REPORT.md"
        step_report(
            gene=gene,
            db=db,
            run_dir=run_dir,
            layer_outdir=layer_outdir,
            doctor_status=status,
            started_at=started,
            duration_s=time.time() - started,
            out_path=report_path,
        )
        # Copy to root output dir so users find it without spelunking
        try:
            shutil.copy2(report_path, output / f"{gene}_RUN_REPORT.md")
        except OSError:
            pass

    logger.info("✅ Done in %.1f min", (time.time() - started) / 60)
    return 0
