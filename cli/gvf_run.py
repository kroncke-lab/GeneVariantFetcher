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
    "GVF_EZPROXY_PREFIX",  # institutional proxy -> Wiley/AHA/Cloudflare publishers
)


def browser_recovery_status() -> dict:
    """Report whether the Playwright browser-recovery tier can actually run.

    Playwright is a declared dependency and Tier 3.5 / source-recovery relies on
    it to fetch free-after-embargo and Cloudflare-gated publisher HTML. When the
    package or its Chromium binary is missing the recovery tier silently
    degrades to abstract-only stubs, so surface it loudly here instead. Cheap:
    resolves the executable path without launching a browser.
    """
    result = {"available": False, "reason": ""}
    try:
        from playwright.sync_api import sync_playwright
    except Exception:
        result["reason"] = "playwright not installed (pip install playwright)"
        return result
    try:
        with sync_playwright() as p:
            exe = p.chromium.executable_path
        if exe and Path(exe).exists():
            result["available"] = True
        else:
            result["reason"] = (
                "chromium binary missing (run: python -m playwright install chromium)"
            )
    except Exception as e:  # driver start / path resolution failed
        result["reason"] = f"chromium unavailable: {str(e)[:120]}"
    return result


def institutional_auth_status() -> dict:
    """Report whether authenticated-publisher recovery is likely to work.

    The Tier 3.5 browser pool and EZproxy routing only get through Cloudflare /
    subscription walls when the run carries *some* institutional credential:
    either an EZproxy prefix/host env var, or logged-in publisher/SSO cookies in
    the local Chrome profile. Neither is required to run GVF, but when both are
    absent the source-recovery pass silently degrades to abstract-only for
    paywalled papers — which is the dominant remaining recall gap. Surface it here
    so it is a visible, fixable precondition rather than a silent ceiling.

    Cheap and prompt-free: it inspects env configuration and whether the cookie
    backend is importable. It does NOT read the Chrome cookie store (that can pop
    a macOS keychain prompt); actual cookie counts are written by the recovery
    pass to ``<run>/source_qc/fetch/auth_status.json``.
    """
    ezproxy_env = next(
        (
            k
            for k in (
                "GVF_EZPROXY_PREFIX",
                "GVF_EZPROXY_HOST",
                "PROXY_LOGIN_PREFIX",
                "PROXY_HOST",
            )
            if os.environ.get(k)
        ),
        None,
    )
    cookie_domains_env = bool(
        os.environ.get("GVF_COOKIE_DOMAINS")
        or os.environ.get("GVF_SSO_COOKIE_DOMAINS")
        or os.environ.get("COOKIE_DOMAIN")
    )
    try:
        import browser_cookie3  # noqa: F401

        cookie_backend = True
    except Exception:
        cookie_backend = False

    ezproxy = bool(ezproxy_env)
    ready = ezproxy or cookie_backend
    if ezproxy:
        reason = f"EZproxy configured via {ezproxy_env}"
    elif cookie_backend:
        reason = (
            "no EZproxy env; will reuse logged-in Chrome cookies — confirm you are "
            "signed in to your institution/publishers in Chrome"
        )
    else:
        reason = (
            "no EZproxy env and browser_cookie3 not installed — paywalled "
            "subscription content cannot be recovered"
        )
    return {
        "ready": ready,
        "ezproxy_configured": ezproxy,
        "ezproxy_env_var": ezproxy_env,
        "cookie_domains_env": cookie_domains_env,
        "cookie_backend_available": cookie_backend,
        "reason": reason,
        "hint": (
            "Set GVF_EZPROXY_PREFIX (e.g. "
            "https://login.proxy.library.<inst>.edu/login?url=) or sign in to your "
            "institution/publishers in Chrome so cookies can be reused."
        ),
    }


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
    # Browser-recovery tier: a hard dependency for source recovery, but easy to
    # leave un-provisioned (package installed, Chromium binary not). Surface it.
    status["browser_recovery"] = browser_recovery_status()
    # Institutional auth readiness: whether authenticated paywalled recovery can
    # actually reach subscription content. Advisory only — never flips status.ok.
    status["institutional_auth"] = institutional_auth_status()
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
    disease: Optional[str] = None,
    include_all_clinigen_phenotypes: bool = False,
    extraction_top_n: Optional[int] = None,
    extraction_priority_offset: int = 0,
    extraction_triage_mode: Optional[str] = None,
    extraction_triage_model: Optional[str] = None,
    extraction_triage_include_defer: bool = False,
    extraction_triage_max_llm: Optional[int] = None,
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

    if disease:
        logger.info("Scoping discovery + Tier-2 filter by disease: %s", disease)

    automated_variant_extraction_workflow(
        gene_symbol=gene,
        email=email,
        output_dir=str(output_dir),
        max_pmids=max_pmids,
        pmids=pmids,
        scout_first=True,
        disease=disease,
        include_all_clinigen_phenotypes=include_all_clinigen_phenotypes,
        extraction_top_n=extraction_top_n,
        extraction_priority_offset=extraction_priority_offset,
        extraction_triage_mode=extraction_triage_mode,
        extraction_triage_model=extraction_triage_model,
        extraction_triage_include_defer=extraction_triage_include_defer,
        extraction_triage_max_llm=extraction_triage_max_llm,
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
    stage_failures: Optional[list[str]] = None,
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
        if stage_failures is not None:
            stage_failures.append(
                f"recovery layers (run_all_layers.py) exited {result.returncode}"
            )
        # Don't hand back a progression.csv path the failed run never produced.
        return None
    return outdir / "progression.csv"


# ---------------------------------------------------------------------------
# Step 4 — gold-free source acquisition QC
# ---------------------------------------------------------------------------


def step_source_qc(
    gene: str,
    run_dir: Path,
    outdir: Path,
    stage_failures: Optional[list[str]] = None,
) -> Optional[Path]:
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
        if stage_failures is not None:
            stage_failures.append(
                f"source acquisition QC (source_acquisition_audit.py) exited "
                f"{result.returncode}"
            )
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
    stage_failures: Optional[list[str]] = None,
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
            if stage_failures is not None:
                stage_failures.append(
                    f"source-recovery fetch (fetch_paywalled.py) exited "
                    f"{fetch_result.returncode}"
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
        if stage_failures is not None:
            stage_failures.append(
                f"source-recovery summarize (summarize_acquisition_outcome.py) exited "
                f"{summary_result.returncode}"
            )
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
        if stage_failures is not None:
            stage_failures.append(
                f"source-recovery refresh (refresh_run_db.py) exited "
                f"{refresh_result.returncode}"
            )
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
# Step 5 — report
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
    stage_failures: Optional[list[str]] = None,
    stage_warnings: Optional[list[str]] = None,
    paper_final_check: Optional[dict] = None,
) -> None:
    """Write a RUN_REPORT.md summarizing the whole run."""
    stage_failures = stage_failures or []
    stage_warnings = stage_warnings or []
    lines: list[str] = []
    lines.append(f"# GVF Run Report — {gene}")
    lines.append("")
    lines.append(f"- Started: {datetime.fromtimestamp(started_at).isoformat()}")
    lines.append(f"- Duration: {duration_s / 60:.1f} min")
    lines.append(f"- Run dir: `{run_dir}`")
    lines.append(f"- DB: `{db}`")
    if stage_failures:
        stage_status = "⚠️ completed with warnings (stage failures recorded)"
    elif stage_warnings:
        stage_status = "✓ core stages ok; best-effort warnings recorded"
    else:
        stage_status = "✓ all stages ok"
    lines.append(f"- Stage status: {stage_status}")
    lines.append("")

    if stage_failures:
        lines.append("## Stage Warnings")
        lines.append("")
        lines.append(
            "Non-fatal stage failures (a subprocess exited nonzero but the run "
            "continued). Results may be incomplete for these stages:"
        )
        lines.append("")
        for failure in stage_failures:
            lines.append(f"  - ⚠️ {failure}")
        lines.append("")

    if stage_warnings:
        lines.append("## Best-effort Warnings")
        lines.append("")
        lines.append(
            "These quality/metadata stages did not remove core extracted evidence, "
            "so they do not change the process exit code:"
        )
        lines.append("")
        for warning in stage_warnings:
            lines.append(f"  - ⚠️ {warning}")
        lines.append("")

    if paper_final_check is not None:
        lines.append("## Paper Final Check")
        lines.append("")
        skip_reason = paper_final_check.get("skipped")
        if isinstance(skip_reason, str):
            lines.append(f"_Skipped: {skip_reason}_")
        else:
            lines.append(
                "| Papers | Checked | Cached/skipped | Source-grounded | Flags | Errors | Missing groups |"
            )
            lines.append("|---:|---:|---:|---:|---:|---:|---:|")
            lines.append(
                "| {papers} | {checked} | {skipped} | {grounded} | {flags} | "
                "{errors} | {missing} |".format(
                    papers=paper_final_check.get("papers", 0),
                    checked=paper_final_check.get("checked", 0),
                    skipped=paper_final_check.get("skipped", 0),
                    grounded=paper_final_check.get("source_grounded", 0),
                    flags=paper_final_check.get("flagged_facts", 0),
                    errors=paper_final_check.get("error", 0),
                    missing=paper_final_check.get("missing_carriers", 0),
                )
            )
            empty_skips = int(paper_final_check.get("skipped_empty_no_source", 0) or 0)
            if empty_skips:
                lines.append("")
                lines.append(
                    f"- {empty_skips} paper(s) had neither extracted count facts "
                    "nor usable source text and were explicitly skipped."
                )
            lines.append("")
            lines.append(
                "Soft results remain in `paper_final_check`; source-grounded "
                "reported/missing groups remain in `paper_carrier_groups`."
            )
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
    browser = doctor_status.get("browser_recovery") or {}
    if browser and not browser.get("available", True):
        lines.append(
            "- ⚠️ Browser recovery is UNAVAILABLE "
            f"({browser.get('reason') or 'unknown'}). Source recovery cannot fetch "
            "free-after-embargo or paywalled publisher HTML and will fall back to "
            "abstract-only stubs. Fix: `python -m playwright install chromium` "
            "(and `pip install playwright` if missing)."
        )
    auth = doctor_status.get("institutional_auth") or {}
    if auth and not auth.get("ready", True):
        lines.append(
            "- ⚠️ No institutional auth path detected "
            f"({auth.get('reason') or 'unknown'}). Paywalled subscription content "
            "(Wiley / AHA / Karger / Sage / Liebert) will NOT be recovered. "
            f"{auth.get('hint') or ''}"
        )
    if doctor_status.get("unlocks", {}).get("ELSEVIER_INSTTOKEN") is False:
        lines.append(
            "- ELSEVIER_INSTTOKEN is unset. Adding it (request from your library) "
            "typically lifts variant_rows by ~30-50 pp for cardiac channel genes "
            "because Heart Rhythm / JACC / Eur Heart J unlock."
        )
    if doctor_status.get("unlocks", {}).get("GVF_EZPROXY_PREFIX") is False:
        lines.append(
            "- GVF_EZPROXY_PREFIX/HOST is unset. Setting it (institutional proxy, "
            "plus a logged-in Chrome for its session cookie) routes Wiley / AHA / "
            "Cloudflare-gated publishers through your campus subscription so stub "
            "papers with a known DOI can be recovered instead of skipped."
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


def step_corpus_sync(run_dir: Path, stage_warnings: Optional[list[str]] = None) -> None:
    """Incrementally fold this run's fetched source into the consolidated corpus.

    Scoped to ``run_dir`` so it is fast; the builder is idempotent and only
    adds new (gene, PMID) papers or upgrades compromised categories, so a
    rerun never re-fetches what the corpus already holds.
    """
    builder = REPO_ROOT / "scripts" / "build_source_corpus.py"
    if not builder.exists():
        return
    logger.info("📦 Step 4.5: corpus sync (folding new source into corpus/)")
    try:
        result = subprocess.run(
            [sys.executable, str(builder), "--apply", "--roots", str(run_dir)],
            capture_output=True,
            text=True,
            timeout=1800,
        )
        for ln in result.stdout.splitlines():
            if ln.startswith(("actions:", "corpus ", "MODE")):
                logger.info("  %s", ln)
        if result.returncode != 0:
            logger.warning(
                "corpus sync returned %s: %s", result.returncode, result.stderr[-300:]
            )
            if stage_warnings is not None:
                stage_warnings.append(
                    f"corpus sync (build_source_corpus.py) exited {result.returncode}"
                )
    except Exception as e:  # noqa: BLE001 - corpus sync is best-effort
        logger.warning("corpus sync failed (non-fatal): %s", e)
        if stage_warnings is not None:
            stage_warnings.append(f"corpus sync failed: {e}")


# ---------------------------------------------------------------------------
# Step 6 — publish to the Variant_Browser review DB (opt-in)
# ---------------------------------------------------------------------------


def _find_review_repo(explicit: Optional[Path]) -> Optional[Path]:
    """Locate the sibling Variant_Browser checkout that owns gvf_publish.sh.

    Resolution order: an explicit ``--review-repo``, then the
    ``GVF_REVIEW_REPO`` / ``VARIANT_BROWSER_DIR`` env vars, then the
    conventional sibling ``<repo parent>/Variant_Browser``. Returns the repo
    path only if its ``scripts/gvf_publish.sh`` actually exists, else None so
    the caller can warn-and-skip rather than fail the run.
    """
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    for env_key in ("GVF_REVIEW_REPO", "VARIANT_BROWSER_DIR"):
        val = os.environ.get(env_key)
        if val:
            candidates.append(Path(val).expanduser())
    candidates.append(REPO_ROOT.parent / "Variant_Browser")
    for repo in candidates:
        if (repo / "scripts" / "gvf_publish.sh").exists():
            return repo
    return None


def step_backfill_metadata(
    *,
    db: Path,
    run_dir: Path,
    email: Optional[str],
    stage_warnings: Optional[list[str]] = None,
) -> None:
    """Fill papers.{first_author, journal, publication_date, doi, pmc_id}.

    Local abstract/artifact caches first; PubMed ESummary fills the rest (mainly
    doi / pmc_id). Runs before report + publish so the scored DB, the run report,
    and any Variant_Browser publish all carry real citations instead of bare PMIDs.
    Best-effort: a failure here never fails the run.
    """
    try:
        from scripts.backfill_paper_metadata import run_backfill

        results = run_backfill(
            [str(db)],
            roots=[run_dir],
            fetch_missing=bool(email),
            email=email,
        )
        filled = results.get(str(db), {})
        detail = ", ".join(f"{col}+{n}" for col, n in filled.items() if n)
        logger.info("📚 paper metadata backfill: %s", detail or "already complete")
    except Exception as e:  # noqa: BLE001
        logger.warning("paper metadata backfill failed (%s); continuing", e)
        if stage_warnings is not None:
            stage_warnings.append(f"metadata backfill failed: {e}")


def step_publish_review(
    *,
    gene: str,
    db: Path,
    disease: Optional[str],
    review_repo: Optional[Path],
    timeout_s: int = 600,
) -> bool:
    """Push this gene's scored DB into the Variant_Browser review DB.

    Opt-in final step. Shells out to ``Variant_Browser/scripts/gvf_publish.sh``,
    which owns the Azure ``vb-curation`` creds and the GVF→browser variant
    matching — GVF does not need DB creds or to duplicate that logic. The
    publish is idempotent on the browser side (re-running replaces the gene's
    carrier data on the current snapshot).

    Best-effort: a missing repo, a launch failure, a timeout, or a non-zero
    exit are all logged and warned, never fatal to the GVF run. Returns True
    only on a clean publish.
    """
    repo = _find_review_repo(review_repo)
    if repo is None:
        logger.warning(
            "publish-review skipped: could not locate "
            "Variant_Browser/scripts/gvf_publish.sh. Pass --review-repo or set "
            "GVF_REVIEW_REPO / VARIANT_BROWSER_DIR."
        )
        return False

    script = repo / "scripts" / "gvf_publish.sh"
    cmd = ["bash", str(script), gene, str(db)]
    if disease:
        cmd.append(disease)
    logger.info("📤 publish-review → %s", " ".join(cmd))
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
    except subprocess.TimeoutExpired:
        logger.warning(
            "publish-review timed out after %ss; GVF run not affected", timeout_s
        )
        return False
    except Exception as e:  # noqa: BLE001 - publish is best-effort
        logger.warning("publish-review failed to launch (%s); GVF run not affected", e)
        return False

    if result.returncode == 0:
        for ln in (result.stdout or "").strip().splitlines()[-3:]:
            logger.info("  %s", ln)
        logger.info("✅ Published %s to the review DB (%s)", gene, repo)
        return True

    logger.warning(
        "publish-review failed (exit %s); GVF run not affected. Output tail:\n%s",
        result.returncode,
        (result.stderr or result.stdout or "")[-800:],
    )
    return False


# ---------------------------------------------------------------------------
# Full-coverage steps (opt-in via --full-coverage). Thin wrappers over
# self-contained pipeline modules; they do not modify the extraction pipeline.
# ---------------------------------------------------------------------------


def step_full_coverage_walk(
    gene: str,
    run_dir: Path,
    *,
    model: str,
    max_workers: int,
    start_offset: int,
    min_new_variants: int,
) -> dict:
    """Walk priority extraction to variant-yield taper for full literature coverage."""
    from pipeline.full_coverage import run_walk_to_taper

    return run_walk_to_taper(
        gene,
        run_dir,
        model=model,
        max_workers=max_workers,
        start_offset=start_offset,
        min_new_variants=min_new_variants,
        logger=logger,
    )


def step_carrier_guard(db: Path) -> dict:
    """Neutralize implausible per-variant carrier counts (cohort/allele-number misreads)."""
    from pipeline.carrier_guard import apply_carrier_guard

    return apply_carrier_guard(db, logger=logger)


def step_vf_enrich(gene: str, db: Path) -> dict:
    """variantFeatures: canonical names + in silico scores + wrong-gene FP quarantine."""
    from pipeline.vf_enrichment import enrich_and_quarantine

    return enrich_and_quarantine(gene, db, logger=logger)


def step_trust_gate(db: Path) -> dict:
    """Soft-quarantine gold-free-implausible count facts into the trust tier."""
    from pipeline.trust_gate import apply_trust_gate

    return apply_trust_gate(db)


def step_paper_final_check(
    db: Path, run_dir: Optional[Path] = None, gene: Optional[str] = None
) -> dict:
    """Per-paper LLM review (default gpt-5.6-sol at xhigh). When the paper's
    on-disk source text is available it produces a carrier/phenotype summary + a
    missed-carrier completeness signal (paper_carrier_groups); otherwise it runs
    the DB-only sniff test over the extracted counts vs their provenance. Records
    soft results in paper_final_check; never mutates or deletes a count.

    Self-gating: returns a ``{"skipped": reason}`` dict (rather than raising) when
    settings can't load, the check is disabled, or the model provider has no
    credentials — so a keyless clone degrades cleanly instead of erroring."""
    try:
        from config.settings import get_settings

        settings = get_settings()
    except Exception as e:  # noqa: BLE001 - misconfigured settings must not fail run
        return {"skipped": f"settings unavailable: {e}"}

    if not settings.paper_final_check_enabled:
        return {"skipped": "disabled (paper_final_check_enabled=false)"}
    model = settings.get_paper_final_check_model()
    if not _paper_check_reachable(model):
        return {"skipped": f"model {model} unreachable (no provider credentials)"}

    from pipeline.paper_final_check import apply_paper_final_check

    return apply_paper_final_check(db, run_dir=run_dir, gene=gene)


def _paper_check_reachable(model: str) -> bool:
    """True when credentials for the final-check model's provider are present, so
    the step skips cleanly on a keyless clone instead of erroring every paper."""
    m = (model or "").lower()
    if m.startswith(("azure_ai/", "azure/")):
        return bool(
            os.environ.get("AZURE_AI_API_KEY") and os.environ.get("AZURE_AI_API_BASE")
        )
    if m.startswith("anthropic/"):
        return bool(os.environ.get("ANTHROPIC_API_KEY"))
    if m.startswith("openai/"):
        return bool(
            os.environ.get("OPENAI_API_KEY") or os.environ.get("AZURE_AI_API_KEY")
        )
    # Bare (unprefixed) Azure deployment names still need Azure credentials.
    base = m.rsplit("/", 1)[-1]
    if base.startswith(("gpt-5", "gpt5", "kimi", "grok", "deepseek")):
        return bool(
            os.environ.get("AZURE_AI_API_KEY") and os.environ.get("AZURE_AI_API_BASE")
        )
    return True  # unknown provider: attempt, let it warn per-paper if misconfigured


def _paper_final_check_error_warning(result: object) -> Optional[str]:
    """Describe per-paper reviewer errors without turning soft QC into data loss.

    The final check intentionally continues after individual LLM failures. Its
    aggregate result therefore needs explicit orchestration-level surfacing;
    otherwise an all-paper outage looks like a successful quality pass.
    """
    if not isinstance(result, dict):
        return None
    try:
        errors = int(result.get("error") or 0)
        checked = int(result.get("checked") or 0)
    except (TypeError, ValueError):
        return "paper final check returned invalid error statistics"
    if errors <= 0:
        return None
    if checked > 0 and errors >= checked:
        return f"paper final check failed for all {checked} checked paper(s)"
    return f"paper final check recorded {errors} error(s) across {checked} checked paper(s)"


# ---------------------------------------------------------------------------
# Entry point (called from cli/__init__.py)
# ---------------------------------------------------------------------------


EXIT_STAGE_WARNINGS = 3  # completed, but a completeness-affecting stage failed


def _write_run_status(
    run_dir: Path,
    gene: str,
    status: str,
    exit_code: int,
    stage_failures: list[str],
    stage_warnings: list[str],
    started_at: float,
) -> None:
    """Write a machine-readable RUN_STATUS.json so a fleet orchestrator keys off
    structured stage status (and the process exit code) instead of scraping
    RUN_REPORT.md. Best-effort: never raises into the run.
    """
    try:
        (run_dir / "RUN_STATUS.json").write_text(
            json.dumps(
                {
                    "gene": gene,
                    "status": status,
                    "exit_code": exit_code,
                    "duration_seconds": int(time.time() - started_at),
                    "stage_failures": list(stage_failures),
                    "stage_warnings": list(stage_warnings),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
    except OSError as e:
        logger.warning("could not write RUN_STATUS.json: %s", e)


def run_gvf_pipeline(
    gene: str,
    email: str,
    output: Path,
    pmid_file: Optional[Path] = None,
    max_pmids: int = 1500,
    resume_dir: Optional[Path] = None,
    include_v12: bool = False,
    skip: Optional[list[str]] = None,
    source_recovery: bool = True,
    source_recovery_timeout_s: int = 120,
    disease: Optional[str] = None,
    include_all_clinigen_phenotypes: bool = False,
    corpus_sync: bool = True,
    publish_review: bool = False,
    review_repo: Optional[Path] = None,
    publish_timeout_s: int = 600,
    extraction_top_n: Optional[int] = None,
    extraction_priority_offset: int = 0,
    extraction_triage_mode: Optional[str] = None,
    extraction_triage_model: Optional[str] = None,
    extraction_triage_include_defer: bool = False,
    extraction_triage_max_llm: Optional[int] = None,
    full_coverage: bool = False,
    extraction_model: Optional[str] = None,
    extraction_workers: Optional[int] = None,
    taper_min_variants: int = 8,
    carrier_guard: bool = True,
    vf_enrich: bool = True,
) -> int:
    """Execute the full pipeline. Returns exit code."""
    initialize_runtime()

    skip = {s.lower() for s in (skip or [])}
    started = time.time()
    # Accumulates non-fatal stage failures (recovery subprocesses that exit
    # nonzero but let the run continue) so the final status and RUN_REPORT.md
    # reflect them instead of "✅ Done" hiding a swallowed error.
    stage_failures: list[str] = []
    # Best-effort stages (corpus sync, metadata backfill, quality passes) whose
    # failure does not remove core evidence: recorded for visibility but they do
    # NOT flip the exit code.
    stage_warnings: list[str] = []

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    output = Path(output).expanduser()
    output.mkdir(parents=True, exist_ok=True)

    # Honor the --email flag as NCBI_EMAIL so the documented cold-start command
    # (gvf gvf-run <GENE> --email you@lab.edu) works on a fresh clone with no
    # pre-populated .env. The explicit flag wins over any env/.env value, which
    # also keeps the doctor reachability check consistent with the email the
    # extraction step actually uses.
    if email:
        os.environ["NCBI_EMAIL"] = email

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

    # Advisory: authenticated paywalled recovery needs an institutional credential.
    auth_status = status.get("institutional_auth") or {}
    if not auth_status.get("ready", True):
        logger.warning("🔒 institutional auth: %s", auth_status.get("reason"))
    elif not auth_status.get("ezproxy_configured"):
        logger.info("🔒 institutional auth: %s", auth_status.get("reason"))

    # Step 2: extract (unless skipped)
    gene = gene.upper()
    run_dir: Optional[Path] = None
    effective_extraction_top_n = extraction_top_n
    seeded_priority_count = 0
    if (
        full_coverage
        and "extract" not in skip
        and (effective_extraction_top_n is None or effective_extraction_top_n <= 0)
    ):
        effective_extraction_top_n = 1000
        logger.info(
            "Full coverage enabled: seeding priority extraction with top %d before walk-to-taper",
            effective_extraction_top_n,
        )
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
                disease=disease,
                include_all_clinigen_phenotypes=include_all_clinigen_phenotypes,
                extraction_top_n=effective_extraction_top_n,
                extraction_priority_offset=extraction_priority_offset,
                extraction_triage_mode=extraction_triage_mode,
                extraction_triage_model=extraction_triage_model,
                extraction_triage_include_defer=extraction_triage_include_defer,
                extraction_triage_max_llm=extraction_triage_max_llm,
            )
            if effective_extraction_top_n and effective_extraction_top_n > 0:
                seeded_priority_count = effective_extraction_top_n
        except Exception as e:
            logger.exception("extract step failed: %s", e)
            return 3

    db = _find_db(run_dir, gene)
    if not db:
        logger.error("No DB produced in %s", run_dir)
        return 4

    # Step 2.5: full-coverage walk-to-taper extraction (opt-in via --full-coverage).
    # Continues priority extraction in offset batches until variant yield tapers,
    # on a high-TPM model — vs the single bounded top-N pass above.
    if full_coverage and "walk" not in skip:
        logger.info("🚶 Step 2.5: full-coverage walk-to-taper extraction")
        try:
            walk_start_offset = extraction_priority_offset + seeded_priority_count
            stats = step_full_coverage_walk(
                gene=gene,
                run_dir=run_dir,
                model=extraction_model or "azure_ai/gpt-5.4",
                max_workers=extraction_workers or 10,
                start_offset=walk_start_offset,
                min_new_variants=taper_min_variants,
            )
            logger.info("full-coverage walk: %s", stats)
            db = _find_db(run_dir, gene) or db
        except Exception as e:  # noqa: BLE001 - best-effort; keep the bounded extraction
            logger.exception("full-coverage walk failed (%s); continuing", e)
            stage_warnings.append(f"full-coverage walk failed: {e}")
    elif full_coverage and "walk" in skip:
        logger.info("⏭️  Step 2.5: full-coverage walk — SKIPPED")

    gold = _find_gold(gene)
    source_qc_summary: Optional[Path] = None
    source_qc_attempted = False
    source_recovery_result: Optional[SourceRecoveryResult] = None
    layer_outdir: Optional[Path] = None
    paper_final_check_result: Optional[dict] = None

    if source_recovery and "source-qc" in skip:
        logger.warning("source recovery requested but source-qc was skipped")

    if source_recovery and "source-qc" not in skip:
        logger.info("🔎 Step 3: source acquisition QC")
        source_qc_summary = step_source_qc(
            gene=gene,
            run_dir=run_dir,
            outdir=run_dir / "source_qc",
            stage_failures=stage_failures,
        )
        source_qc_attempted = True
        if "source-recovery" not in skip and source_qc_summary is not None:
            logger.info("🛟 Step 4: source recovery")
            source_recovery_result = step_source_recovery(
                gene=gene,
                run_dir=run_dir,
                source_qc_dir=source_qc_summary.parent,
                gold=gold,
                run_recovery_layers=("layers" not in skip),
                timeout_s=source_recovery_timeout_s,
                stage_failures=stage_failures,
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
            stage_failures=stage_failures,
        )
    else:
        logger.info("⏭️  Step 3: layers — SKIPPED")

    # Step 3.5: full-coverage data-quality passes on the finalized DB (opt-in).
    # carrier-guard neutralizes cohort/allele-number-as-carrier misreads; vf-enrich
    # attaches canonical names + in silico scores and quarantines wrong-gene FPs.
    if full_coverage and carrier_guard and "carrier-guard" not in skip:
        logger.info("🛡️  Step 3.5: carrier-count guard")
        try:
            step_carrier_guard(db=db)
        except Exception as e:  # noqa: BLE001
            logger.exception("carrier guard failed (%s); continuing", e)
            stage_warnings.append(f"carrier guard failed: {e}")
    elif full_coverage and not carrier_guard:
        logger.info("⏭️  Step 3.5: carrier-count guard — SKIPPED")
    if full_coverage and vf_enrich and "vf-enrich" not in skip:
        logger.info("🔬 Step 3.6: variantFeatures enrich + wrong-gene FP quarantine")
        try:
            step_vf_enrich(gene=gene, db=db)
        except Exception as e:  # noqa: BLE001
            logger.exception("vf-enrich failed (%s); continuing", e)
            stage_warnings.append(f"vf-enrich failed: {e}")
    elif full_coverage and not vf_enrich:
        logger.info("⏭️  Step 3.6: variantFeatures enrich — SKIPPED")

    # Step 3.7: per-fact trust gate — soft-quarantine gold-free-implausible count
    # facts into a two-tier (trusted/quarantine) DB. Default ON: this is the
    # primary automated quality control for unattended operation. It only sets
    # trust_tier (never NULLs/deletes), so a failure degrades to the DDL default
    # ('trusted') rather than losing data — hence a warning, not a hard failure.
    if "trust-gate" not in skip:
        logger.info("🚦 Step 3.7: per-fact trust gate")
        try:
            logger.info("trust gate: %s", step_trust_gate(db=db))
        except Exception as e:  # noqa: BLE001
            logger.exception("trust gate failed (%s); continuing", e)
            stage_warnings.append(f"trust gate failed: {e}")
    else:
        logger.info("⏭️  Step 3.7: trust gate — SKIPPED")

    # Step 3.8: per-paper final check (sniff test). A strong reasoning model
    # (default gpt-5.6-sol at xhigh, per config.settings) reviews each paper's
    # extracted counts against their captured provenance and records a soft
    # verdict in the paper_final_check table — it never mutates or deletes a
    # count. Default ON. The step self-gates (disabled / unreachable / bad
    # settings → clean skip); a genuine failure warns, it does not fail the run.
    if "paper-final-check" not in skip:
        logger.info("🧪 Step 3.8: per-paper final check")
        try:
            paper_final_check_result = step_paper_final_check(
                db=db, run_dir=run_dir, gene=gene
            )
            if isinstance(paper_final_check_result, dict) and isinstance(
                paper_final_check_result.get("skipped"), str
            ):
                logger.info(
                    "⏭️  Step 3.8: paper final check — %s",
                    paper_final_check_result["skipped"],
                )
            else:
                logger.info("paper final check: %s", paper_final_check_result)
            pfc_warning = _paper_final_check_error_warning(paper_final_check_result)
            if pfc_warning:
                logger.warning("%s", pfc_warning)
                stage_warnings.append(pfc_warning)
        except Exception as e:  # noqa: BLE001
            logger.exception("paper final check failed (%s); continuing", e)
            stage_warnings.append(f"paper final check failed: {e}")
    else:
        logger.info("⏭️  Step 3.8: paper final check — SKIPPED")

    # Step 4: gold-free source QC. Only run if Step 3 did not already attempt it —
    # a failed Step-3 attempt also leaves source_qc_summary None, and re-running
    # would repeat the failing subprocess and double-record the stage failure.
    if not source_qc_attempted and "source-qc" not in skip:
        logger.info("🔎 Step 4: source acquisition QC")
        source_qc_summary = step_source_qc(
            gene=gene,
            run_dir=run_dir,
            outdir=run_dir / "source_qc",
            stage_failures=stage_failures,
        )
        source_qc_attempted = True
    elif "source-qc" in skip:
        logger.info("⏭️  Step 4: source acquisition QC — SKIPPED")
    else:
        logger.info("⏭️  Step 4: source acquisition QC — already ran")

    # Step 4.5: fold newly-fetched source into the consolidated corpus cache
    # (idempotent incremental merge — adds new papers / upgrades compromised
    # categories only, so the next run reuses them instead of re-fetching).
    if corpus_sync and "corpus-sync" not in skip:
        step_corpus_sync(run_dir=run_dir, stage_warnings=stage_warnings)
    elif "corpus-sync" in skip or not corpus_sync:
        logger.info("⏭️  Step 4.5: corpus sync — SKIPPED")

    # Step 4.6: backfill paper bibliographic metadata (journal/year/doi/pmc/author)
    # so the report and any review-DB publish render real citations, not bare PMIDs.
    if "metadata-backfill" not in skip:
        logger.info("📚 Step 4.6: paper metadata backfill")
        step_backfill_metadata(
            db=db, run_dir=run_dir, email=email, stage_warnings=stage_warnings
        )
    else:
        logger.info("⏭️  Step 4.6: paper metadata backfill — SKIPPED")

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
            stage_failures=stage_failures,
            stage_warnings=stage_warnings,
            paper_final_check=paper_final_check_result,
        )
        # Copy to root output dir so users find it without spelunking
        try:
            shutil.copy2(report_path, output / f"{gene}_RUN_REPORT.md")
        except OSError:
            pass

    # Step 6: publish to the Variant_Browser review DB (opt-in)
    if publish_review and "publish-review" not in skip:
        logger.info("📤 Step 6: publish to review DB")
        step_publish_review(
            gene=gene,
            db=db,
            disease=disease,
            review_repo=review_repo,
            timeout_s=publish_timeout_s,
        )
    elif publish_review and "publish-review" in skip:
        logger.info("⏭️  Step 6: publish to review DB — SKIPPED")

    exit_code = EXIT_STAGE_WARNINGS if stage_failures else 0
    status = "completed_with_warnings" if stage_failures else "completed"
    _write_run_status(
        run_dir, gene, status, exit_code, stage_failures, stage_warnings, started
    )

    if stage_failures:
        logger.warning(
            "⚠️  Done in %.1f min with %d stage failure(s) (exit %d) — see "
            "RUN_STATUS.json / RUN_REPORT.md:",
            (time.time() - started) / 60,
            len(stage_failures),
            exit_code,
        )
        for failure in stage_failures:
            logger.warning("   - %s", failure)
    else:
        logger.info("✅ Done in %.1f min", (time.time() - started) / 60)
    if stage_warnings:
        logger.info(
            "   (%d non-fatal stage warning(s) recorded in RUN_STATUS.json)",
            len(stage_warnings),
        )
    return exit_code
