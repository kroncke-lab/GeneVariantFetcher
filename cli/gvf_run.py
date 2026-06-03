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
     gene_variant_fetcher_gold_standard/normalized/. The old KCNH2 v12
     manual-recovery DB is opt-in only because it is not cold-start behavior.
  4. source-qc — writes a gold-free source acquisition/replay worklist.
  5. source-recovery — optional: fetch source-QC queue, summarize actual
     usable full text, and staged-refresh accepted sources.
  6. report — writes RUN_REPORT.md summarizing what happened, what
     scored, and what's gated on credentials.

Cold-start invocation (one command, no flags beyond gene + email):

    gvf gvf-run KCNQ1 --email me@uni.edu --output results/

Skip the LLM-bound extract step and just score an existing DB through
the recovery layers:

    gvf gvf-run KCNQ1 --email me@uni.edu --output results/ --skip extract
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

from utils.bootstrap import (
    has_llm_provider_key,
    initialize_runtime,
    llm_provider_key_status,
)

REPO_ROOT = Path(__file__).resolve().parents[1]

logger = logging.getLogger("gvf_run")


# ---------------------------------------------------------------------------
# Step 1 — doctor
# ---------------------------------------------------------------------------


REQUIRED_ENV = ("NCBI_EMAIL",)
RECOMMENDED_ENV = (
    "NCBI_API_KEY",
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
    initialize_runtime()

    status: dict = {
        "required": {},
        "llm_providers": {},
        "recommended": {},
        "unlocks": {},
        "ok": True,
    }
    for k in REQUIRED_ENV:
        present = bool(os.environ.get(k))
        status["required"][k] = present
        if not present:
            status["ok"] = False
    provider_status = llm_provider_key_status()
    status["llm_providers"] = provider_status
    status["required"]["LLM provider key"] = has_llm_provider_key()
    if not has_llm_provider_key():
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
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "recall_recovery" / "run_all_layers.py"),
        "--gene",
        gene,
        "--db",
        str(db),
        "--pmc-dir",
        str(run_dir / "pmc_fulltext"),
        "--outdir",
        str(outdir),
        "--backup",
    ]
    if gold:
        cmd.extend(["--gold", str(gold)])
    else:
        logger.warning(
            "No gold standard CSV found for %s. Recovery layers will run in "
            "DB-PMID mode and per-layer recall scoring will be skipped.",
            gene,
        )
    if with_v12:
        cmd.extend(["--with-v12", str(with_v12)])

    logger.info("layers → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("recovery layers failed: %s", result.stderr[-400:])
    return outdir / "progression.csv"


# ---------------------------------------------------------------------------
# Step 4 — gold-free source acquisition QC
# ---------------------------------------------------------------------------


def step_source_qc(gene: str, run_dir: Path, outdir: Path) -> Optional[Path]:
    """Run the no-gold source acquisition audit. Returns summary JSON path."""
    outdir.mkdir(parents=True, exist_ok=True)
    summary = outdir / "source_acquisition_summary.json"
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "recall_audit" / "source_acquisition_audit.py"),
        "--gene",
        gene,
        "--run-dir",
        str(run_dir),
        "--out",
        str(outdir / "source_acquisition_worklist.csv"),
        "--summary-out",
        str(summary),
        "--fetch-input-out",
        str(outdir / "fetch_input.csv"),
        "--source-override-out",
        str(outdir / "source_override.csv"),
    ]
    logger.info("source-qc → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("source acquisition QC failed: %s", result.stderr[-400:])
        return None
    return summary


# ---------------------------------------------------------------------------
# Step 5 — optional no-gold source recovery loop
# ---------------------------------------------------------------------------


@dataclass
class SourceRecoveryResult:
    """Artifacts produced by the optional source-recovery loop."""

    fetch_dir: Path
    outcome_summary: Path
    fetched_source_override: Path
    refresh_summary: Optional[Path] = None
    active_db: Optional[Path] = None
    layer_outdir: Optional[Path] = None


def _csv_has_nonheader_rows(path: Path) -> bool:
    if not path.exists():
        return False
    rows = [
        line
        for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
        if line.strip()
    ]
    return len(rows) > 1


def _latest_refresh_summary_after(run_dir: Path, started_at: float) -> Optional[Path]:
    candidates = [
        p
        for p in run_dir.glob("refresh_*/refresh_summary.json")
        if p.stat().st_mtime >= started_at
    ]
    if not candidates:
        return None
    return sorted(candidates, key=lambda p: p.stat().st_mtime)[-1]


def step_source_recovery(
    *,
    gene: str,
    run_dir: Path,
    source_qc_dir: Path,
    gold: Optional[Path],
    run_recovery_layers: bool,
    timeout_s: int,
) -> Optional[SourceRecoveryResult]:
    """Run fetch → summarize → staged refresh from source-QC artifacts."""
    worklist = source_qc_dir / "source_acquisition_worklist.csv"
    fetch_input = source_qc_dir / "fetch_input.csv"
    source_override = source_qc_dir / "source_override.csv"
    fetch_dir = source_qc_dir / "fetch"
    outcome_summary = source_qc_dir / "acquisition_outcome_summary.json"
    fetched_source_override = source_qc_dir / "fetched_source_override.csv"

    if not worklist.exists():
        logger.warning("source recovery skipped: missing %s", worklist)
        return None

    if _csv_has_nonheader_rows(fetch_input):
        fetch_dir.mkdir(parents=True, exist_ok=True)
        fetch_cmd = [
            sys.executable,
            str(REPO_ROOT / "scripts" / "fetch_paywalled.py"),
            "--input",
            str(fetch_input),
            "--output",
            str(fetch_dir),
            "--headless",
            "--timeout-s",
            str(timeout_s),
        ]
        logger.info("source-recovery fetch → %s", " ".join(fetch_cmd))
        fetch_result = subprocess.run(fetch_cmd, capture_output=True, text=True)
        if fetch_result.returncode != 0:
            logger.warning(
                "source-recovery fetch returned %s; continuing with landed artifacts: %s",
                fetch_result.returncode,
                (fetch_result.stderr or fetch_result.stdout)[-400:],
            )
    else:
        logger.info("source-recovery fetch skipped: no fetch rows in %s", fetch_input)

    summary_cmd = [
        sys.executable,
        str(
            REPO_ROOT / "scripts" / "recall_audit" / "summarize_acquisition_outcome.py"
        ),
        "--gene",
        gene,
        "--worklist",
        str(worklist),
        "--fetch-output-dir",
        str(fetch_dir),
        "--out",
        str(outcome_summary),
        "--source-override-out",
        str(fetched_source_override),
    ]
    if gold:
        summary_cmd.extend(["--gold", str(gold)])
    logger.info("source-recovery summarize → %s", " ".join(summary_cmd))
    summary_result = subprocess.run(summary_cmd, capture_output=True, text=True)
    if summary_result.returncode != 0:
        logger.error("source-recovery summary failed: %s", summary_result.stderr[-400:])
        return None

    result = SourceRecoveryResult(
        fetch_dir=fetch_dir,
        outcome_summary=outcome_summary,
        fetched_source_override=fetched_source_override,
    )
    override_csvs = [
        path
        for path in (source_override, fetched_source_override)
        if _csv_has_nonheader_rows(path)
    ]
    if not override_csvs:
        logger.info("source-recovery refresh skipped: no usable source overrides")
        return result

    refresh_started_at = time.time()
    refresh_cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "refresh_run_db.py"),
        "--gene",
        gene,
        "--run-dir",
        str(run_dir),
        "--stage-extractions",
        "--only-forced-pmids",
    ]
    if not run_recovery_layers:
        refresh_cmd.append("--skip-recovery")
    for override_csv in override_csvs:
        refresh_cmd.extend(["--source-override-csv", str(override_csv)])

    logger.info("source-recovery refresh → %s", " ".join(refresh_cmd))
    refresh_result = subprocess.run(refresh_cmd, capture_output=True, text=True)
    if refresh_result.returncode != 0:
        logger.error("source-recovery refresh failed: %s", refresh_result.stderr[-800:])
        return result

    refresh_summary = _latest_refresh_summary_after(run_dir, refresh_started_at)
    if refresh_summary is None:
        logger.warning("source-recovery refresh summary not found under %s", run_dir)
        return result

    result.refresh_summary = refresh_summary
    try:
        refresh_data = json.loads(refresh_summary.read_text(encoding="utf-8"))
    except Exception:
        refresh_data = {}
    active_db = refresh_data.get("active_db")
    if active_db:
        result.active_db = Path(active_db)
    layers_dir = refresh_summary.parent / "layers"
    if layers_dir.exists():
        result.layer_outdir = layers_dir

    final_summary_cmd = [
        sys.executable,
        str(
            REPO_ROOT / "scripts" / "recall_audit" / "summarize_acquisition_outcome.py"
        ),
        "--gene",
        gene,
        "--worklist",
        str(worklist),
        "--fetch-output-dir",
        str(fetch_dir),
        "--refresh-summary",
        str(refresh_summary),
        "--out",
        str(outcome_summary),
        "--source-override-out",
        str(fetched_source_override),
    ]
    if gold:
        final_summary_cmd.extend(["--gold", str(gold)])
    logger.info("source-recovery final summarize → %s", " ".join(final_summary_cmd))
    final_summary = subprocess.run(final_summary_cmd, capture_output=True, text=True)
    if final_summary.returncode != 0:
        logger.warning(
            "source-recovery final summary failed: %s", final_summary.stderr[-400:]
        )
    return result


# ---------------------------------------------------------------------------
# Step 6 — report
# ---------------------------------------------------------------------------


def step_report(
    gene: str,
    db: Path,
    run_dir: Path,
    layer_outdir: Optional[Path],
    source_qc_summary: Optional[Path],
    source_recovery: Optional[SourceRecoveryResult],
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
    lines.append(f"- Duration: {duration_s / 60:.1f} min")
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
    lines.append("LLM provider keys:")
    for k, v in doctor_status.get("llm_providers", {}).items():
        lines.append(f"  - {k}: {'✓' if v else '–'}")
    lines.append("Subscription unlocks (not required, but lift recall):")
    for k, v in doctor_status.get("unlocks", {}).items():
        lines.append(f"  - {k}: {'✓' if v else '–'}")
    lines.append(f"- NCBI reachable: {doctor_status.get('ncbi_reachable')}")
    lines.append("")

    if layer_outdir and (layer_outdir / "progression.json").exists():
        progression = json.loads((layer_outdir / "progression.json").read_text())
        rows = progression.get("progression", [])
        scoring_enabled = bool(progression.get("scoring_enabled"))
        lines.append(
            "## Recall Progression (per layer)"
            if scoring_enabled
            else "## Recovery Progression"
        )
        lines.append("")
        if scoring_enabled:
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
            lines.append("| Layer | ClinVar added | PubTator added | Figures added |")
            lines.append("|---|---:|---:|---:|")
            for r in rows:
                lines.append(
                    "| {layer} | {clinvar} | {pubtator} | {figures} |".format(
                        layer=r.get("layer", "?"),
                        clinvar=r.get("clinvar_added", "—"),
                        pubtator=r.get("pubtator_added", "—"),
                        figures=r.get("figures_added", "—"),
                    )
                )
            lines.append("")
            lines.append(
                "_No gold standard CSV was available, so recovery ran without recall scoring._"
            )
            lines.append("")
    else:
        lines.append("## Recall Progression")
        lines.append("")
        lines.append(
            "_Recovery layers were skipped or did not produce progression output._"
        )
        lines.append(
            "Review internal QC signals instead (see docs/NEW_GENE_RUNBOOK.md)."
        )
        lines.append("")

    if source_qc_summary and source_qc_summary.exists():
        source_qc = json.loads(source_qc_summary.read_text())
        coverage = source_qc.get("pmid_coverage", {})
        lines.append("## Source Acquisition QC (Gold-Free)")
        lines.append("")
        lines.append("| Signal | PMIDs | Coverage |")
        lines.append("|---|---:|---:|")
        for label, key in (
            ("Usable full text now", "usable_fulltext_current"),
            ("Selected for fetch", "selected_for_fetch"),
            ("Selected for source refresh", "selected_for_source_refresh"),
            ("Manual or blocked", "selected_for_manual_or_blocked"),
            ("Zero-variant usable full text", "zero_variant_usable_fulltext"),
        ):
            item = coverage.get(key) or {}
            pmids = item.get("pmids", "—")
            total = item.get("total_pmids", "—")
            cov = item.get("coverage")
            cov_text = "—" if cov is None else f"{cov * 100:.1f}%"
            lines.append(f"| {label} | {pmids}/{total} | {cov_text} |")
        lines.append("")
        qc_dir = source_qc_summary.parent
        lines.append(f"- Worklist: `{qc_dir / 'source_acquisition_worklist.csv'}`")
        lines.append(f"- Fetch queue: `{qc_dir / 'fetch_input.csv'}`")
        lines.append(f"- Refresh source overrides: `{qc_dir / 'source_override.csv'}`")
        lines.append("")

    if source_recovery and source_recovery.outcome_summary.exists():
        outcome = json.loads(source_recovery.outcome_summary.read_text())
        coverage = outcome.get("pmid_coverage", {})
        recall = outcome.get("pmid_recall", {})
        lines.append("## Source Recovery")
        lines.append("")
        lines.append("| Signal | PMIDs | Coverage/Recall |")
        lines.append("|---|---:|---:|")
        for label, key in (
            ("Selected for fetch", "selected_for_fetch_download"),
            ("Fetch attempted", "fetch_attempted"),
            ("Usable full text downloaded", "usable_fulltext_downloaded"),
            ("Source refresh successful", "source_refresh_successful"),
        ):
            bucket = recall.get(key) or coverage.get(key) or {}
            pmids = bucket.get("pmids", "—")
            denominator = bucket.get("gold_pmids") or bucket.get("total_pmids") or "—"
            value = bucket.get("recall")
            if value is None:
                value = bucket.get("coverage")
            value_text = "—" if value is None else f"{value * 100:.1f}%"
            lines.append(f"| {label} | {pmids}/{denominator} | {value_text} |")
        lines.append("")
        lines.append(f"- Fetch output: `{source_recovery.fetch_dir}`")
        lines.append(f"- Outcome summary: `{source_recovery.outcome_summary}`")
        if source_recovery.refresh_summary:
            lines.append(f"- Refresh summary: `{source_recovery.refresh_summary}`")
        lines.append("")

    lines.append("## Next steps")
    if doctor_status.get("unlocks", {}).get("ELSEVIER_INSTTOKEN") is False:
        lines.append(
            "- ELSEVIER_INSTTOKEN is unset. Adding it (request from your library) "
            "typically lifts variant_rows by ~30-50 pp for cardiac channel genes "
            "because Heart Rhythm / JACC / Eur Heart J unlock."
        )
    if not any(doctor_status.get("llm_providers", {}).values()):
        lines.append(
            "- No LLM provider key is set. Add OPENAI_API_KEY, AZURE_AI_API_KEY, "
            "or ANTHROPIC_API_KEY before running extraction."
        )
    if layer_outdir:
        lines.append(f"- Inspect per-layer outputs under `{layer_outdir}`.")
    if source_qc_summary and source_qc_summary.exists() and not source_recovery:
        lines.append(
            "- For no-gold source recovery, run fetch_paywalled.py on the "
            "source QC fetch queue, then refresh_run_db.py with "
            "`--source-override-csv` and `--stage-extractions`, or rerun "
            "`gvf-run` with `--source-recovery`."
        )
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
    include_v12: bool = False,
    skip: Optional[list[str]] = None,
    source_recovery: bool = False,
    source_recovery_timeout_s: int = 120,
) -> int:
    """Execute the full pipeline. Returns exit code."""
    initialize_runtime()

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
        if not any(status.get("llm_providers", {}).values()):
            logger.error(
                "  missing: one of OPENAI_API_KEY, AZURE_AI_API_KEY, ANTHROPIC_API_KEY"
            )
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

    gold = _find_gold(gene)
    source_qc_summary: Optional[Path] = None
    source_recovery_result: Optional[SourceRecoveryResult] = None
    layer_outdir: Optional[Path] = None

    if source_recovery and "source-qc" in skip:
        logger.warning("source recovery requested but source-qc was skipped")

    if source_recovery and "source-qc" not in skip:
        logger.info("🔎 Step 3: source acquisition QC")
        source_qc_summary = step_source_qc(
            gene=gene,
            run_dir=run_dir,
            outdir=run_dir / "source_qc",
        )
        if "source-recovery" not in skip and source_qc_summary is not None:
            logger.info("🛟 Step 4: source recovery")
            source_recovery_result = step_source_recovery(
                gene=gene,
                run_dir=run_dir,
                source_qc_dir=source_qc_summary.parent,
                gold=gold,
                run_recovery_layers=("layers" not in skip),
                timeout_s=source_recovery_timeout_s,
            )
            if source_recovery_result and source_recovery_result.active_db:
                db = source_recovery_result.active_db
            if source_recovery_result and source_recovery_result.layer_outdir:
                layer_outdir = source_recovery_result.layer_outdir
        elif "source-recovery" in skip:
            logger.info("⏭️  Step 4: source recovery — SKIPPED")

    # Step 3: layers
    if source_recovery_result and source_recovery_result.layer_outdir:
        logger.info(
            "⏭️  recovery layers already ran during source recovery (%s)",
            source_recovery_result.layer_outdir,
        )
    elif "layers" not in skip:
        logger.info("🧬 Step 3: recovery layers")
        v12 = _find_v12_db(gene) if include_v12 else None
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

    # Step 4: gold-free source QC
    if source_qc_summary is None and "source-qc" not in skip:
        logger.info("🔎 Step 4: source acquisition QC")
        source_qc_summary = step_source_qc(
            gene=gene,
            run_dir=run_dir,
            outdir=run_dir / "source_qc",
        )
    elif "source-qc" in skip:
        logger.info("⏭️  Step 4: source acquisition QC — SKIPPED")
    else:
        logger.info("⏭️  Step 4: source acquisition QC — already ran")

    # Step 5: report
    if "report" not in skip:
        logger.info("📝 Step 5: report")
        report_path = run_dir / "RUN_REPORT.md"
        step_report(
            gene=gene,
            db=db,
            run_dir=run_dir,
            layer_outdir=layer_outdir,
            source_qc_summary=source_qc_summary,
            source_recovery=source_recovery_result,
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
