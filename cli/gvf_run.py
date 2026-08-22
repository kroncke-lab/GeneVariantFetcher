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
     gene_variant_fetcher_gold_standard/normalized/. The canonical KCNH2
     baseline is opt-in only because merging it is not cold-start behavior.
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
import sqlite3
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from utils.bootstrap import (
    has_llm_provider_key,
    initialize_runtime,
    llm_provider_key_status,
)
from utils.env_utils import get_env_bool, local_data_discovery_disabled
from utils.local_storage import (
    EXTERNAL_PATHS,
    external_link_target,
    external_path_state,
    external_volume_name,
    local_storage_allowed,
)

REPO_ROOT = Path(__file__).resolve().parents[1]

logger = logging.getLogger("gvf_run")


def validate_prior_vf_enrichment(db: Path) -> dict:
    """Prove a resumed DB already has complete VariantFeatures coverage.

    A validated skip avoids re-reading the external ~38 GB warehouse when the
    active DB was already enriched and wrong-gene rows were quarantined. Schema
    presence alone is insufficient: every live variant must have an enrichment
    row, and the quarantine audit table must exist.
    """
    result = {
        "valid": False,
        "variants": 0,
        "enrichment_rows": 0,
        "uncovered": 0,
        "quarantined": 0,
    }
    try:
        connection = sqlite3.connect(f"file:{Path(db).resolve()}?mode=ro", uri=True)
        try:
            tables = {
                str(row[0])
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
            }
            required = {"variants", "vf_enrichment", "quarantined_variants"}
            if not required.issubset(tables):
                result["reason"] = (
                    "missing prior VariantFeatures audit tables: "
                    + ", ".join(sorted(required - tables))
                )
                return result
            result["variants"] = int(
                connection.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
            )
            result["enrichment_rows"] = int(
                connection.execute("SELECT COUNT(*) FROM vf_enrichment").fetchone()[0]
            )
            result["uncovered"] = int(
                connection.execute(
                    """
                    SELECT COUNT(*)
                    FROM variants AS v
                    LEFT JOIN vf_enrichment AS e USING (variant_id)
                    WHERE e.variant_id IS NULL
                    """
                ).fetchone()[0]
            )
            result["quarantined"] = int(
                connection.execute(
                    "SELECT COUNT(*) FROM quarantined_variants"
                ).fetchone()[0]
            )
        finally:
            connection.close()
    except (OSError, sqlite3.Error) as exc:
        result["reason"] = str(exc)
        return result

    if result["uncovered"]:
        result["reason"] = (
            f"{result['uncovered']} live variant(s) lack VariantFeatures audit rows"
        )
        return result
    result["valid"] = True
    result["reason"] = "all live variants retain prior VariantFeatures audit coverage"
    return result


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


def _apply_local_storage_checks(status: dict) -> None:
    """Fold external-storage link state into a doctor status dict.

    A *dangling* link means the volume was expected and is not mounted — the
    accident case (drive removed, wrong machine). Block at Step 1 rather than
    running cold and re-fetching the whole corpus over the network. A *local*
    real directory also blocks, because writing there builds a second corpus on
    the internal disk that no later run reads; `GVF_ALLOW_LOCAL_CORPUS=1` opts
    in on a machine with no external volume. An *absent* link is a legitimate
    fresh checkout or collaborator setup, so it does not block — the write-side
    guard catches that case later with mount instructions. `--skip doctor`
    overrides either way.

    Split out of ``doctor`` so it is testable without doctor's network probes.
    """
    status["local_storage"] = {
        name: external_path_state(name) for name in EXTERNAL_PATHS
    }
    corpus_override = os.environ.get("GVF_CORPUS_DIR", "").strip()
    for name, state in status["local_storage"].items():
        if name == "corpus" and corpus_override:
            # A deliberate path override bypasses the repo-local corpus link;
            # its target may legitimately start absent and be created by the
            # corpus builder.
            continue
        if state == "dangling":
            # Name the drive from the link, so this is right on any machine.
            volume = external_volume_name(name)
            remedy = (
                f"attach '{volume}' and mount it at /Volumes/{volume}"
                if volume
                else f"restore its target ({external_link_target(name)})"
            )
            status["required"][f"{name}/ external storage ({remedy})"] = False
            status["ok"] = False
        elif state == "local" and not local_storage_allowed():
            status["required"][
                f"{name}/ is local (set GVF_ALLOW_LOCAL_CORPUS=1 to opt in)"
            ] = False
            status["ok"] = False


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
    # EZproxy accepts four spellings (PREFIX/HOST × GVF_/plain); the BMPR2 run's
    # doctor printed "GVF_EZPROXY_PREFIX: –" while the proxy was fully configured
    # via GVF_EZPROXY_HOST. Report the resolved state, not one literal env var.
    try:
        from harvesting.browser_html import ezproxy as _ezproxy

        status["unlocks"]["GVF_EZPROXY_PREFIX"] = _ezproxy.is_configured()
    except Exception:  # noqa: BLE001 - doctor must never crash on an import
        pass
    # Browser-recovery tier: a hard dependency for source recovery, but easy to
    # leave un-provisioned (package installed, Chromium binary not). Surface it.
    status["browser_recovery"] = browser_recovery_status()
    _apply_local_storage_checks(status)
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
    status_path = run_dir / "RUN_STATUS.json"
    if status_path.is_file():
        try:
            status = json.loads(status_path.read_text(encoding="utf-8"))
            active_ref = status.get("active_db") if isinstance(status, dict) else None
            if isinstance(active_ref, str) and active_ref.strip():
                active_db = Path(active_ref.strip()).expanduser()
                if not active_db.is_absolute():
                    active_db = run_dir / active_db
                if active_db.is_file():
                    return active_db
                logger.warning(
                    "RUN_STATUS.json names missing active_db %s; falling back to discovery",
                    active_db,
                )
        except (OSError, ValueError, TypeError) as exc:
            logger.warning(
                "Could not read %s (%s); falling back to discovery", status_path, exc
            )
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
    """Return the consolidated KCNH2 baseline; other genes have no v12 input."""
    if gene.upper() != "KCNH2":
        return None
    baseline = REPO_ROOT / "validation_runs" / "canonical_baseline" / "KCNH2.db"
    return baseline if baseline.is_file() else None


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
        "--supplement-input-out",
        str(outdir / "supplement_input.csv"),
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
    supplement_input = source_qc_dir / "supplement_input.csv"
    supplement_source_override = source_qc_dir / "supplement_source_override.csv"
    linked_source_override = source_qc_dir / "linked_supplement_source_override.csv"
    linked_report = source_qc_dir / "linked_supplement_report.csv"
    source_override = source_qc_dir / "source_override.csv"
    fetch_dir = source_qc_dir / "fetch"
    outcome_summary = source_qc_dir / "acquisition_outcome_summary.json"
    fetched_source_override = source_qc_dir / "fetched_source_override.csv"

    if not worklist.exists():
        logger.warning("source recovery skipped: missing %s", worklist)
        return None

    # Linked-supplement recovery: papers whose markup advertised supplement
    # files that never landed on disk. Self-gating — it walks the harvest dir
    # and no-ops unless a paper has recorded links and an empty supplements
    # directory, so it costs one Europe PMC request per *affected* paper only.
    # Runs first because it is free (no credentials, no browser) and its
    # overrides join the same refresh below.
    linked_cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "fetch_linked_supplements.py"),
        "--gene",
        gene,
        "--harvest-dir",
        str(run_dir / "pmc_fulltext"),
        "--report-out",
        str(linked_report),
        "--source-override-out",
        str(linked_source_override),
    ]
    logger.info("source-recovery linked-supplement fetch → %s", " ".join(linked_cmd))
    linked_result = subprocess.run(linked_cmd, capture_output=True, text=True)
    if linked_result.returncode != 0:
        logger.warning(
            "source-recovery linked-supplement fetch returned %s: %s",
            linked_result.returncode,
            (linked_result.stderr or linked_result.stdout)[-400:],
        )
        if stage_failures is not None:
            stage_failures.append(
                "source-recovery linked-supplement fetch "
                f"(fetch_linked_supplements.py) exited {linked_result.returncode}"
            )

    if _csv_has_nonheader_rows(supplement_input):
        supplement_cmd = [
            sys.executable,
            str(REPO_ROOT / "scripts" / "fetch_elsevier_supplements.py"),
            "--gene",
            gene,
            "--harvest-dir",
            str(run_dir / "pmc_fulltext"),
            "--input",
            str(supplement_input),
            "--source-override-out",
            str(supplement_source_override),
        ]
        logger.info(
            "source-recovery supplement-only fetch → %s", " ".join(supplement_cmd)
        )
        supplement_result = subprocess.run(
            supplement_cmd, capture_output=True, text=True
        )
        if supplement_result.returncode != 0:
            logger.warning(
                "source-recovery supplement-only fetch returned %s: %s",
                supplement_result.returncode,
                (supplement_result.stderr or supplement_result.stdout)[-400:],
            )
            if stage_failures is not None:
                stage_failures.append(
                    "source-recovery supplement-only fetch "
                    f"(fetch_elsevier_supplements.py) exited "
                    f"{supplement_result.returncode}"
                )
    else:
        logger.info(
            "source-recovery supplement-only fetch skipped: no rows in %s",
            supplement_input,
        )

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
        for path in (
            source_override,
            supplement_source_override,
            linked_source_override,
            fetched_source_override,
        )
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
    # No explicit env=: refresh_run_db.py re-extracts, so it must inherit the
    # GVF_LLM_TRACE_DIR / GVF_LLM_TRACE_RUN_ID that open_trace_session exported
    # into os.environ, or its model calls land outside this run's trace tree.
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
        gate = paper_final_check.get("gate")
        if isinstance(skip_reason, str):
            lines.append(f"_Skipped: {skip_reason}_")
            lines.append("")
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
        if paper_final_check.get("review_skipped"):
            lines.append(
                "The live reviewer was skipped; Step 3.9 recomposed durable "
                "stored findings against the current extraction payload."
            )
            lines.append("")
        if isinstance(gate, dict) and not gate.get("skipped"):
            lines.append(
                "| Gate-applied facts | Applied bindings | Unresolved | "
                "Unverified quote | Stale papers | Advisory reason | Trusted | "
                "Quarantined |"
            )
            lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
            lines.append(
                "| {facts} | {bindings} | {unresolved} | {unverified} | "
                "{stale} | {advisory} | {trusted} | {quarantine} |".format(
                    facts=gate.get("applied_facts", 0),
                    bindings=gate.get("applied_bindings", 0),
                    unresolved=gate.get("unresolved_actionable", 0),
                    unverified=gate.get("unverified_actionable", 0),
                    stale=gate.get("stale_papers", 0),
                    advisory=gate.get("advisory_reason_flags", 0),
                    trusted=gate.get("trusted", 0),
                    quarantine=gate.get("quarantine", 0),
                )
            )
            lines.append("")
            lines.append(
                "Raw count values remain unchanged. Exact, source-quoted "
                "fact/field contradictions are composed into `field_trust`; default "
                "scoring projects those quarantined fields as missing while "
                "retaining variant identity. Reported/missing groups remain "
                "in `paper_carrier_groups` as a re-extraction signal."
            )
        else:
            lines.append(
                "Final-check findings remain recorded in "
                "`paper_final_check`; no actionable trust composition was "
                "reported for this run."
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
    storage = doctor_status.get("local_storage") or {}
    if storage:
        lines.append("External storage:")
        for name, state in storage.items():
            mark = {"linked": "✓", "dangling": "✗", "absent": "–"}.get(state, "?")
            lines.append(f"  - {name}/: {mark} ({state})")
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
            ("Selected for supplement-only fetch", "selected_for_supplement_fetch"),
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
        lines.append(f"- Supplement-only queue: `{qc_dir / 'supplement_input.csv'}`")
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


def step_corpus_sync(
    run_dir: Path,
    stage_warnings: Optional[list[str]] = None,
    gene: Optional[str] = None,
) -> None:
    """Incrementally fold this run's fetched source into the consolidated corpus.

    Scoped to ``run_dir`` so it is fast; the builder is idempotent and only
    adds new (gene, PMID) papers or upgrades compromised categories, so a
    rerun never re-fetches what the corpus already holds. ``gene`` is passed
    as both the inference hint and the corpus scope: a run-dir root has no gene
    path component below it, while an incremental run must not scan/hash the
    unrelated external corpus.
    """
    builder = REPO_ROOT / "scripts" / "build_source_corpus.py"
    if not builder.exists():
        return
    logger.info("📦 Step 4.5: corpus sync (folding new source into corpus/)")
    cmd = [sys.executable, str(builder), "--apply", "--roots", str(run_dir)]
    corpus_override = os.environ.get("GVF_CORPUS_DIR", "").strip()
    if corpus_override:
        cmd += ["--out", corpus_override]
    if gene:
        cmd += ["--assume-gene", gene, "--gene", gene]
    try:
        result = subprocess.run(
            cmd,
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

    The sibling guess is skipped when local-data discovery is disabled, so the
    offline suite cannot pick up a developer's real Variant_Browser checkout.
    """
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    for env_key in ("GVF_REVIEW_REPO", "VARIANT_BROWSER_DIR"):
        val = os.environ.get(env_key)
        if val:
            candidates.append(Path(val).expanduser())
    if not local_data_discovery_disabled():
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
    pmid_file: Optional[Path] = None,
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

    publish_env = os.environ.copy()
    publish_env.pop("GVF_PMID_FILE", None)
    if pmid_file is not None:
        manifest = Path(pmid_file).expanduser().resolve()
        if not manifest.is_file():
            logger.warning(
                "publish-review refused: pinned PMID manifest does not exist: %s",
                manifest,
            )
            return False
        publish_env["GVF_PMID_FILE"] = str(manifest)
    elif publish_env.get("GVF_FULL_DB_REPLACE") != "1":
        logger.warning(
            "publish-review refused: no pinned PMID manifest was supplied. Pass "
            "--pmid-file, or explicitly set GVF_FULL_DB_REPLACE=1 for a full "
            "replacement."
        )
        return False

    script = repo / "scripts" / "gvf_publish.sh"
    cmd = ["bash", str(script), gene, str(db)]
    if disease:
        cmd.append(disease)
    logger.info("📤 publish-review → %s", " ".join(cmd))
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            env=publish_env,
        )
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


def step_count_repair(db: Path, *, rules: Optional[set[str]] = None) -> dict:
    """Run selected deterministic count-repair rules."""
    from pipeline.count_repair import apply_count_repair

    return apply_count_repair(db, rules=rules, logger=logger)


def _count_recovery_enabled() -> bool:
    """Read the default-off feature flag independently of full settings.

    Prefers the validated ``Settings.count_recovery_enabled`` field so the
    setting is actually consumed (it was previously dead: nothing read it and
    this function went straight to the env var). Falls back to the raw env var
    only when Settings cannot be constructed — otherwise an unrelated invalid
    setting would turn explicit enablement into a silent skip, and
    ``step_count_recovery`` should be the thing that raises.
    """
    try:
        from config.settings import get_settings

        return bool(get_settings().count_recovery_enabled)
    except Exception:  # noqa: BLE001 - fall back to the raw flag
        return get_env_bool("COUNT_RECOVERY_ENABLED", False)


def step_count_recovery(gene: str, db: Path, run_dir: Optional[Path] = None) -> dict:
    """Fill NULL per-variant counts for variants already found in each paper.

    The additive counterpart to Step 3.5: the guards remove counts that should
    not be there, this fills counts the extractor never emitted. Never
    overwrites an existing count, and every value it writes is grounded in a
    verbatim quote from the source that proves a per-variant role.

    Imports the source resolver from ``pipeline.count_recovery`` (packaged), not
    ``scripts.recover_counts`` (not in the wheel) — the old import degraded this
    stage to a warning-and-continue in every installed deployment.
    """
    from config.settings import get_settings
    from pipeline.count_recovery import (
        default_source_roots,
        make_source_resolver,
        recover_counts,
    )
    from utils.llm_utils import litellm_completion

    settings = get_settings()
    model = settings.get_count_recovery_model()
    settings.require_provider_credentials_for_model(model, label="count recovery")
    fields = tuple(
        f.strip()
        for f in (settings.count_recovery_fields or "").split(",")
        if f.strip()
    )
    roots: list[Path] = []
    if run_dir is not None and (run_dir / "pmc_fulltext").is_dir():
        roots.append(run_dir / "pmc_fulltext")
    for root in default_source_roots(db):
        if root not in roots:
            roots.append(root)

    stats = recover_counts(
        db,
        gene,
        source_for_pmid=make_source_resolver(gene, roots),
        llm_caller=litellm_completion,
        model=model,
        reasoning_effort=settings.count_recovery_reasoning_effort,
        fields=fields or ("carriers",),
        max_variants_per_call=settings.count_recovery_max_variants_per_call,
        max_source_chars=settings.count_recovery_max_source_chars,
    )
    stats.pop("result_objects", None)
    logger.info(
        "count recovery: %d/%d gap(s) filled across %d paper(s) (%d rejected)",
        stats["counts_written"],
        stats["gaps_found"],
        stats["papers_with_gaps"],
        stats["counts_rejected"],
    )
    return stats


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


def step_paper_final_check_gate(db: Path) -> dict:
    """Compose exact source-grounded final-check flags into count trust.

    Raw count columns and variant identities remain unchanged. Only the named
    fact/field trusted projection is quarantined.
    """
    from config.settings import get_settings
    from pipeline.paper_final_check_gate import apply_paper_final_check_gate

    settings = get_settings()
    if not settings.paper_final_check_gate_enabled:
        return {"skipped": "disabled (paper_final_check_gate_enabled=false)"}
    return apply_paper_final_check_gate(
        db,
        min_severity=settings.paper_final_check_gate_min_severity,
        require_source_grounded=(
            settings.paper_final_check_gate_require_source_grounded
        ),
    )


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
EXIT_INSTITUTIONAL_BLOCK = (
    5  # refused to START a full run: institutional access degraded
)


def _attempt_ezproxy_self_heal(access) -> Optional[object]:
    """Try a headless EZproxy session refresh and re-probe; None when not healed.

    Fires only for the failure classes a relogin can fix — an expired session
    (login redirect) or a missing session cookie — and only when the operator
    has already bootstrapped the dedicated browser profile
    (``scripts/ezproxy_relogin.py --bootstrap``). Disable with
    ``GVF_EZPROXY_AUTOHEAL=0``. A Cloudflare wall or probe error is not
    healable by a new cookie, so those still block immediately.
    """
    if (os.environ.get("GVF_EZPROXY_AUTOHEAL") or "").strip().lower() in (
        "0",
        "false",
        "no",
    ):
        return None
    if not (getattr(access, "login_redirect", False) or access.ez_cookies == 0):
        return None
    try:
        from utils.env_utils import local_data_discovery_disabled

        # The profile check probes a home-directory path; the offline suite
        # must never discover developer state through it.
        if local_data_discovery_disabled():
            return None
        from scripts.ezproxy_relogin import (
            profile_is_bootstrapped,
            resolve_profile_dir,
        )
    except ImportError:
        return None  # installed package without scripts/ — manual refresh only
    profile_dir = resolve_profile_dir()
    if not profile_is_bootstrapped(profile_dir):
        logger.info(
            "institutional preflight: no bootstrapped EZproxy profile at %s — "
            "run scripts/ezproxy_relogin.py --bootstrap once to enable self-heal.",
            profile_dir,
        )
        return None
    logger.info(
        "🔁 institutional preflight: attempting EZproxy session self-heal "
        "(profile: %s)…",
        profile_dir,
    )
    try:
        result = subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "scripts" / "ezproxy_relogin.py"),
                "--skip-verify",
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )
    except Exception as e:  # noqa: BLE001 - self-heal must never brick the run
        logger.warning("EZproxy self-heal errored (%s); keeping original result.", e)
        return None
    if result.returncode != 0:
        logger.warning(
            "EZproxy self-heal did not complete (exit %d): %s",
            result.returncode,
            (result.stderr or result.stdout or "").strip()[-300:],
        )
        return None
    try:
        from cli.institutional_preflight import probe_institutional_access

        healed = probe_institutional_access()
    except Exception as e:  # noqa: BLE001
        logger.warning("post-heal re-probe errored (%s); keeping original result.", e)
        return None
    if healed.viable:
        logger.info("✅ EZproxy self-heal succeeded — session refreshed.")
    else:
        logger.warning(
            "EZproxy self-heal refreshed cookies but access is still blocked (%s).",
            healed.reason,
        )
    return healed


@dataclass
class TraceSession:
    """Where this gvf-run invocation writes its LLM traces, and under whose id.

    A gvf-run that skips extraction reuses a *previous* run's directory. Writing
    into that run's ``llm_traces/`` and rebuilding its ``trace_manifest.json``
    under the current run id destroyed the earlier manifest and presented stale
    traces as this run's artifact. So a reused directory gets its own
    ``llm_traces_<run id>/`` and its own report file, and the earlier run's
    artifacts are never touched.
    """

    run_id: str
    trace_root: Path
    report_path: Path
    enabled: bool
    owns_run_dir: bool


def _new_gvf_run_id() -> str:
    return f"gvfrun-{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}-{uuid.uuid4().hex[:8]}"


def _existing_run_id(run_dir: Path) -> Optional[str]:
    """The run id recorded by the extraction workflow that created ``run_dir``."""
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        return json.loads(manifest_path.read_text(encoding="utf-8")).get("run_id")
    except (OSError, json.JSONDecodeError):
        return None


def open_trace_session(
    run_dir: Path, gvf_run_id: str, *, extraction_ran: bool
) -> TraceSession:
    """Configure process-wide tracing for the post-extraction gvf-run stages.

    ``gvf-run`` never configured tracing itself: the only production
    ``configure_llm_tracing`` lived inside the extraction workflow, so a
    ``--skip extract`` run's count recovery and per-paper final check ran
    completely untraced unless an operator happened to export
    ``GVF_LLM_TRACE_DIR``.
    """
    from utils.llm_trace import (
        configure_llm_tracing,
        resolve_trace_location,
        tracing_enabled_by_environment,
    )
    from utils.llm_trace_html import TRACE_REPORT_NAME

    enabled = tracing_enabled_by_environment()
    if extraction_ran:
        # The extraction workflow already selected and EXPORTED a per-run child
        # plus its id, so resolve_trace_location returns them verbatim and this
        # run's extraction and post-extraction records share one tree and one id.
        default_root = run_dir / "llm_traces"
        default_run_id = _existing_run_id(run_dir) or gvf_run_id
        owns = True
    else:
        # Skip-extract reuses an OLDER run's directory. Its own child keeps this
        # pass from writing into — and later relabelling — that run's manifest.
        default_root = run_dir / f"llm_traces_{gvf_run_id}"
        default_run_id = gvf_run_id
        owns = False
    location = resolve_trace_location(default_run_id, default_root=default_root)
    trace_root, run_id = location.root, location.run_id
    report_path = (
        run_dir / TRACE_REPORT_NAME
        if owns
        else run_dir / f"{gvf_run_id}_{TRACE_REPORT_NAME}"
    )
    configure_llm_tracing(trace_root, run_id=run_id, enabled=enabled)
    # Subprocess stages that make model calls (scripts/refresh_run_db.py
    # re-extraction) inherit the resolved child and id, so their calls land in
    # this run's trace tree instead of nowhere — and, because the id is set, they
    # use the directory verbatim instead of appending another level.
    if enabled:
        os.environ["GVF_LLM_TRACE_DIR"] = str(trace_root)
        os.environ["GVF_LLM_TRACE_RUN_ID"] = run_id
        logger.info("🧾 LLM traces -> %s (run %s)", trace_root, run_id)
    else:
        logger.info("🧾 LLM tracing disabled by GVF_LLM_TRACE_ENABLED")
    return TraceSession(
        run_id=run_id,
        trace_root=trace_root,
        report_path=report_path,
        enabled=enabled,
        owns_run_dir=owns,
    )


def finalize_trace_artifacts(
    session: TraceSession,
    gene: str,
    stage_warnings: list[str],
) -> Optional[dict]:
    """Rebuild this run's manifest + HTML report, refusing to relabel another run."""
    from utils.llm_trace import (
        TRACE_MANIFEST_NAME,
        TraceRunMismatchError,
        build_trace_manifest,
    )
    from utils.llm_trace_html import build_trace_html_report

    if not session.enabled or not session.trace_root.is_dir():
        return None
    try:
        manifest = build_trace_manifest(
            session.trace_root,
            output_path=session.trace_root / TRACE_MANIFEST_NAME,
            run_id=session.run_id,
        )
    except TraceRunMismatchError as exc:
        warning = f"refused to relabel another run's trace manifest: {exc}"
        logger.warning("%s", warning)
        stage_warnings.append(warning)
        return None
    except Exception as exc:  # noqa: BLE001 - tracing must not change the outcome
        warning = f"could not build LLM trace manifest: {exc}"
        logger.warning("%s", warning)
        stage_warnings.append(warning)
        return None

    logger.info(
        "🧾 LLM trace manifest: %s calls, %s decision events, integrity=%s → %s",
        manifest["llm_call_count"],
        manifest["decision_event_count"],
        manifest["verification"]["level"],
        session.trace_root / TRACE_MANIFEST_NAME,
    )
    try:
        data = build_trace_html_report(
            session.trace_root,
            output_path=session.report_path,
            run_dir=session.report_path.parent,
            title=f"{gene} · LLM curation trace review",
            run_id=session.run_id,
        )
    except Exception as exc:  # noqa: BLE001
        warning = f"could not build LLM trace report: {exc}"
        logger.warning("%s", warning)
        stage_warnings.append(warning)
        return {"manifest": manifest}
    logger.info("🌐 LLM trace report → %s", session.report_path)
    for gap in data["coverage"]["missing_decision_links"]:
        logger.info("   route-coverage gap: %s (%s)", gap["stage"], gap["reason"])
    if data.get("omissions"):
        logger.info(
            "   %d trace body/record omission(s) listed in the report",
            len(data["omissions"]),
        )
    return {
        "run_id": session.run_id,
        "trace_root": str(session.trace_root),
        "report": str(session.report_path),
        "llm_call_count": manifest["llm_call_count"],
        "decision_event_count": manifest["decision_event_count"],
        "integrity_level": manifest["verification"]["level"],
        "integrity_errors": manifest["verification"]["errors"],
        "missing_decision_links": data["coverage"]["missing_decision_links"],
        "omissions": len(data.get("omissions") or []),
        "sharded": bool(data.get("sharded")),
    }


def _write_run_status(
    run_dir: Path,
    gene: str,
    status: str,
    exit_code: int,
    stage_failures: list[str],
    stage_warnings: list[str],
    started_at: float,
    active_db: Path,
    integrity: Optional[dict] = None,
    count_recovery: Optional[dict] = None,
    llm_trace: Optional[dict] = None,
    gold_free_run: bool = False,
) -> None:
    """Write a machine-readable RUN_STATUS.json so a fleet orchestrator keys off
    structured stage status (and the process exit code) instead of scraping
    RUN_REPORT.md.

    ``active_db`` names the finalized database consumed by the post-extraction
    stages. Paths inside ``run_dir`` are stored relative to the directory that
    contains RUN_STATUS.json; paths outside it are stored as absolute paths.
    Consumers can therefore resolve a relative value against
    ``RUN_STATUS.json.parent``. Best-effort: never raises into the run.
    """
    try:
        try:
            active_db_ref = str(active_db.relative_to(run_dir))
        except ValueError:
            resolved_run_dir = run_dir.resolve()
            resolved_active_db = active_db.resolve()
            try:
                active_db_ref = str(resolved_active_db.relative_to(resolved_run_dir))
            except ValueError:
                active_db_ref = str(resolved_active_db)

        payload = {
            "gene": gene,
            "status": status,
            "severity": "warning" if stage_failures else "ok",
            "exit_code": exit_code,
            "active_db": active_db_ref,
            "duration_seconds": int(time.time() - started_at),
            "stage_failures": list(stage_failures),
            "stage_warnings": list(stage_warnings),
            "gold_access": {
                "disabled": bool(gold_free_run),
                "mode": "disabled" if gold_free_run else "auto_discovery_allowed",
            },
        }
        if integrity is not None:
            payload["source_integrity"] = integrity
        if count_recovery is not None:
            payload["count_recovery"] = count_recovery
        if llm_trace is not None:
            payload["llm_trace"] = llm_trace
        (run_dir / "RUN_STATUS.json").write_text(
            json.dumps(payload, indent=2),
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
    gold_free_run: bool = False,
    skip: Optional[list[str]] = None,
    source_recovery: bool = True,
    source_recovery_timeout_s: int = 120,
    allow_degraded_institutional: bool = False,
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
    require_vf_enrich: bool = False,
) -> int:
    """Execute the full pipeline. Returns exit code.

    Wraps the pipeline so the trace identity this invocation publishes
    (``GVF_LLM_TRACE_DIR`` / ``GVF_LLM_TRACE_RUN_ID``, exported for nested and
    subprocess stages) is restored on every exit path. A leaked run id would make
    the NEXT in-process run adopt this run's identity and write into its trace
    directory — which then mixes runs and makes the manifest rebuild refuse.
    """
    previous_dir = os.environ.get("GVF_LLM_TRACE_DIR")
    previous_run = os.environ.get("GVF_LLM_TRACE_RUN_ID")
    try:
        return _run_gvf_pipeline(
            gene=gene,
            email=email,
            output=output,
            pmid_file=pmid_file,
            max_pmids=max_pmids,
            resume_dir=resume_dir,
            include_v12=include_v12,
            gold_free_run=gold_free_run,
            skip=skip,
            source_recovery=source_recovery,
            source_recovery_timeout_s=source_recovery_timeout_s,
            allow_degraded_institutional=allow_degraded_institutional,
            disease=disease,
            include_all_clinigen_phenotypes=include_all_clinigen_phenotypes,
            corpus_sync=corpus_sync,
            publish_review=publish_review,
            review_repo=review_repo,
            publish_timeout_s=publish_timeout_s,
            extraction_top_n=extraction_top_n,
            extraction_priority_offset=extraction_priority_offset,
            extraction_triage_mode=extraction_triage_mode,
            extraction_triage_model=extraction_triage_model,
            extraction_triage_include_defer=extraction_triage_include_defer,
            extraction_triage_max_llm=extraction_triage_max_llm,
            full_coverage=full_coverage,
            extraction_model=extraction_model,
            extraction_workers=extraction_workers,
            taper_min_variants=taper_min_variants,
            carrier_guard=carrier_guard,
            vf_enrich=vf_enrich,
            require_vf_enrich=require_vf_enrich,
        )
    finally:
        for key, value in (
            ("GVF_LLM_TRACE_DIR", previous_dir),
            ("GVF_LLM_TRACE_RUN_ID", previous_run),
        ):
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _run_gvf_pipeline(
    *,
    gene: str,
    email: str,
    output: Path,
    pmid_file: Optional[Path] = None,
    max_pmids: int = 1500,
    resume_dir: Optional[Path] = None,
    include_v12: bool = False,
    gold_free_run: bool = False,
    skip: Optional[list[str]] = None,
    source_recovery: bool = True,
    source_recovery_timeout_s: int = 120,
    allow_degraded_institutional: bool = False,
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
    require_vf_enrich: bool = False,
) -> int:
    initialize_runtime()

    if gold_free_run and include_v12:
        logger.error("--gold-free-run is incompatible with --with-v12")
        return EXIT_STAGE_WARNINGS

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

    # A "full-dataset" run is one that actually attempts paywalled/supplementary
    # recovery: source recovery ON, not a --pmid-file calibration run, and the
    # recovery steps not skipped. Only these are guarded/audited for institutional
    # access; fast (--no-source-recovery) and calibration (--pmid-file) runs are
    # expected to be abstract-limited and must pass through untouched.
    institutional_run = (
        source_recovery
        and pmid_file is None
        and "source-recovery" not in skip
        and "source-qc" not in skip
    )

    # Advisory (cheap, presence-only) — always shown.
    auth_status = status.get("institutional_auth") or {}
    if auth_status:
        if not auth_status.get("ready", True):
            logger.warning("🔒 institutional auth: %s", auth_status.get("reason"))
        elif not auth_status.get("ezproxy_configured"):
            logger.info("🔒 institutional auth: %s", auth_status.get("reason"))

    # Live institutional-access preflight guard. A full-dataset run must not START
    # if paywall access is dead, otherwise the orchestrator silently degrades every
    # unreachable paywalled paper to an abstract-only stub and still reports
    # success. Unlike the advisory above (presence-only), this routes a real Wiley
    # DOI through EZproxy and verifies licensed full text comes back.
    if not institutional_run:
        logger.info("🔓 institutional preflight skipped (not a full-recovery run)")
    elif os.environ.get("GVF_PREFLIGHT_SKIP"):
        logger.warning(
            "🔓 institutional preflight skipped (GVF_PREFLIGHT_SKIP set) — live "
            "paywall access was NOT verified."
        )
    else:
        logger.info("🔐 institutional preflight: probing live paywall access…")
        access = None
        try:
            from cli.institutional_preflight import (
                format_block_message,
                probe_institutional_access,
            )

            access = probe_institutional_access()
        except Exception as e:  # noqa: BLE001 - a guard bug must not brick every run
            logger.exception("institutional preflight probe errored: %s", e)
            stage_warnings.append(
                "institutional preflight probe could not run; access unverified"
            )
        if access is not None and access.should_block:
            healed = _attempt_ezproxy_self_heal(access)
            if healed is not None:
                access = healed
        if access is not None:
            for ln in access.lines:
                logger.info("   %s", ln)
            if access.should_block and not allow_degraded_institutional:
                print(format_block_message(access))
                logger.error(
                    "institutional preflight blocked the run (exit %d). Override "
                    "with --allow-degraded-institutional if this is intentional.",
                    EXIT_INSTITUTIONAL_BLOCK,
                )
                return EXIT_INSTITUTIONAL_BLOCK
            if access.should_block:
                logger.warning(
                    "⚠️  institutional access DEGRADED (%s) — proceeding anyway via "
                    "--allow-degraded-institutional; expect abstract-only coverage "
                    "for paywalled papers.",
                    access.reason,
                )
                stage_failures.append(
                    f"institutional access degraded (overridden): {access.reason}"
                )
            else:
                logger.info(
                    "✅ institutional preflight: live full-text access confirmed."
                )

    # Step 2: extract (unless skipped)
    gene = gene.upper()
    run_dir: Optional[Path] = None
    # Identifies THIS invocation, independent of the run dir it works in, so a
    # skip-extract pass can isolate its own traces without borrowing (and then
    # overwriting) the identity of the run that produced the directory.
    gvf_run_id = _new_gvf_run_id()
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

    # Trace the WHOLE remaining lifecycle, not just extraction. On the
    # skip-extract path this points at a run-scoped directory so the previous
    # run's manifest and report are never rewritten or relabelled.
    trace_session = open_trace_session(
        run_dir, gvf_run_id, extraction_ran="extract" not in skip
    )

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

    if gold_free_run:
        logger.info(
            "🔒 Gold access disabled for this run; recovery and reporting remain gold-free"
        )
        gold = None
    else:
        gold = _find_gold(gene)
    source_qc_summary: Optional[Path] = None
    source_qc_attempted = False
    source_recovery_result: Optional[SourceRecoveryResult] = None
    layer_outdir: Optional[Path] = None
    paper_final_check_result: Optional[dict] = None
    paper_final_check_gate_result: Optional[dict] = None

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

    # Step 3.45: adopt counts already read from patient-level figures. This must
    # precede every guard and the trust gate so they judge the adopted values,
    # and precede count recovery so it does not spend a model call on a filled
    # gap. Suspicious raw values remain intact and are handled by trust projection.
    if "count-repair" not in skip:
        logger.info("🧮 Step 3.45: adopt deterministic figure counts")
        try:
            summary = step_count_repair(db=db, rules={"adopt_figure_counts"})
            if summary.get("rows_changed"):
                logger.info(
                    "   adopted counts on %d rows (%d penetrance rows created)",
                    summary["rows_changed"],
                    summary.get("rows_created", 0),
                )
        except Exception as e:  # noqa: BLE001
            logger.exception("figure-count adoption failed (%s); continuing", e)
            stage_warnings.append(f"figure-count adoption failed: {e}")
    else:
        logger.info("⏭️  Step 3.45: adopt deterministic figure counts — SKIPPED")

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

    # Step 3.55: targeted count recovery — fill NULL per-variant counts for
    # variants the extractor already found. Default OFF (COUNT_RECOVERY_ENABLED)
    # pending a clean v2 measurement; unlike the guards above it spends LLM calls.
    count_recovery_status: Optional[dict] = None
    if _count_recovery_enabled() and "count-recovery" not in skip:
        if "trust-gate" in skip:
            # Recovered counts land as `quarantine` for Step 3.7 to promote.
            # Running 3.55 with 3.7 skipped would leave every recovered value
            # quarantined and unevaluated, so refuse rather than half-apply.
            failure = (
                "count recovery refused: it lands values as quarantine for the "
                "trust gate to promote, and --skip trust-gate would leave them "
                "unevaluated. Run without --skip trust-gate, or also skip "
                "count-recovery."
            )
            logger.error("%s", failure)
            stage_failures.append(failure)
            count_recovery_status = {"status": "refused", "reason": failure}
        else:
            logger.info("🔢 Step 3.55: targeted count recovery")
            try:
                count_recovery_stats = step_count_recovery(
                    gene=gene, db=db, run_dir=run_dir
                )
                count_recovery_status = {
                    key: value
                    for key, value in count_recovery_stats.items()
                    if key != "results"
                }
                count_recovery_status["status"] = "ran"
                if count_recovery_stats.get("all_batches_failed"):
                    # Total failure is a stage FAILURE (non-zero exit), not a
                    # warning: every paper of paid model calls produced nothing.
                    failure = (
                        "count recovery failed on every attempted paper "
                        f"({count_recovery_stats.get('batch_failures', 0)} batch "
                        "failure(s)); no counts recovered"
                    )
                    logger.error("%s", failure)
                    stage_failures.append(failure)
                    count_recovery_status["status"] = "failed"
                elif count_recovery_stats.get("papers_failed"):
                    warning = (
                        "count recovery had model/parse failures for "
                        f"{count_recovery_stats['papers_failed']} paper(s)"
                    )
                    logger.warning("%s", warning)
                    stage_warnings.append(warning)
            except Exception as e:  # noqa: BLE001
                # Explicit enablement that cannot run is a failure, not a
                # silently-successful run with zero recovered counts.
                logger.exception("count recovery failed (%s)", e)
                stage_failures.append(f"count recovery failed: {e}")
                count_recovery_status = {"status": "error", "reason": str(e)}

    # NOT gated on full_coverage. This is the only check that validates a
    # variant's gene against an independent reference, and it was reachable only
    # in full-coverage discovery mode — so every calibrated `--pmid-file` run,
    # which is how the collaborator-facing review sets are built, silently
    # skipped it. Even the "SKIPPED" line below sat inside the same branch, so
    # nothing said it had not run.
    # ...but only when the warehouse is actually reachable. It is a 36GB SQLite
    # on an external volume; forcing it unconditionally would hang a collaborator
    # laptop or CI, and this same DB previously wedged the offline test suite for
    # ~685s. Unavailable is a loud warning, never a silent no-op.
    vf_db_path = None
    try:
        from utils.gene_metadata import default_variantfeatures_db_path
        from utils.env_utils import local_data_discovery_disabled

        if not local_data_discovery_disabled():
            vf_db_path = default_variantfeatures_db_path()
    except Exception as e:  # noqa: BLE001 - never fail a run on a QC lookup
        logger.debug("vf-enrich: could not resolve variantFeatures DB (%s)", e)

    vf_enrich_required = require_vf_enrich or publish_review

    if vf_enrich and "vf-enrich" not in skip and (full_coverage or vf_db_path):
        logger.info("🔬 Step 3.6: variantFeatures enrich + wrong-gene FP quarantine")
        try:
            vf_stats = step_vf_enrich(gene=gene, db=db)
            logger.info("vf-enrich: %s", vf_stats)
            enriched = isinstance(vf_stats, dict) and vf_stats.get("enriched") is True
            quarantined = (
                isinstance(vf_stats, dict) and vf_stats.get("quarantined") is True
            )
            if not enriched or not quarantined:
                issue = (
                    "vf-enrich did not run "
                    f"({vf_stats.get('reason') or vf_stats.get('rc') or vf_stats}); "
                    "wrong-gene FP quarantine was NOT applied"
                )
                logger.warning("%s", issue)
                (stage_failures if vf_enrich_required else stage_warnings).append(issue)
        except Exception as e:  # noqa: BLE001
            logger.exception("vf-enrich failed (%s); continuing", e)
            issue = f"vf-enrich failed: {e}"
            (stage_failures if vf_enrich_required else stage_warnings).append(issue)
    elif vf_enrich and "vf-enrich" not in skip:
        issue = (
            "vf-enrich SKIPPED: no readable variantFeatures DB "
            "(set VARIANTFEATURES_DB). The wrong-gene FP quarantine — the only "
            "check that validates a variant's gene against an independent "
            "reference — was NOT applied to this run."
        )
        logger.warning("%s", issue)
        (stage_failures if vf_enrich_required else stage_warnings).append(issue)
    elif vf_enrich and "vf-enrich" in skip:
        prior_vf = validate_prior_vf_enrichment(db)
        if prior_vf.get("valid"):
            logger.info(
                "⏭️  Step 3.6: VariantFeatures warehouse read — SKIPPED; "
                "validated prior DB coverage: %s",
                prior_vf,
            )
        else:
            issue = (
                "vf-enrich skipped without complete prior DB coverage "
                f"({prior_vf.get('reason') or prior_vf}); wrong-gene FP quarantine "
                "cannot be verified"
            )
            logger.warning("%s", issue)
            (stage_failures if vf_enrich_required else stage_warnings).append(issue)
    else:
        logger.info("⏭️  Step 3.6: variantFeatures enrich — SKIPPED (disabled)")
        issue = "vf-enrich disabled; wrong-gene FP quarantine was NOT applied"
        (stage_failures if vf_enrich_required else stage_warnings).append(issue)

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

    # Step 3.8: per-paper final check (sniff test). PARKED (default off since
    # 2026-07-26): one @xhigh call per paper cost more time and money than its
    # measured effect justified for a step that only RECORDS findings. Set
    # PAPER_FINAL_CHECK_ENABLED=1 to revive it, together with
    # PAPER_FINAL_CHECK_GATE_ENABLED=1 -- Step 3.9 without a live reviewer can
    # only refuse stale findings and fail acceptance. When enabled, a strong
    # reasoning model reviews each paper's extracted counts against
    # source/provenance and records structured findings; it never mutates or
    # deletes a raw count, and Step 3.9 composes exact source-grounded fact/field
    # flags into the trusted projection.
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

    # Step 3.9: actionable final-check gate. Never blanket-quarantines a paper:
    # only at/above-threshold source-grounded findings carrying exact fact IDs and
    # named count fields can affect trust. Raw values and identities remain in
    # the DB. Quote-verified missing carrier groups fail run acceptance because
    # there is no safe count to synthesize; they require re-extraction.
    # Recompose from the durable current paper_final_check rows even when the
    # live reviewer is skipped/unreachable. A transient outage must not silently
    # lift a previously grounded quarantine after Step 3.7 resets structural
    # reasons. Operators can explicitly skip only this deterministic composer
    # with --skip paper-final-check-gate.
    if "paper-final-check-gate" not in skip:
        logger.info("🔒 Step 3.9: paper final-check trust composition")
        try:
            paper_final_check_gate_result = step_paper_final_check_gate(db=db)
            logger.info(
                "paper final-check trust composition: %s",
                paper_final_check_gate_result,
            )
            if isinstance(paper_final_check_result, dict):
                paper_final_check_result["gate"] = paper_final_check_gate_result
            else:
                # The composer intentionally runs from durable rows even when
                # the live reviewer is skipped. Keep that enforcement visible
                # in RUN_REPORT instead of dropping the whole section.
                paper_final_check_result = {
                    "papers": paper_final_check_gate_result.get("papers", 0),
                    "checked": 0,
                    "skipped": paper_final_check_gate_result.get("papers", 0),
                    "source_grounded": paper_final_check_gate_result.get(
                        "source_grounded_papers", 0
                    ),
                    "flagged_facts": paper_final_check_gate_result.get(
                        "flags_total", 0
                    ),
                    "missing_carriers": paper_final_check_gate_result.get(
                        "missing_carriers", 0
                    ),
                    "error": 0,
                    "review_skipped": True,
                    "gate": paper_final_check_gate_result,
                }
            unresolved = int(
                paper_final_check_gate_result.get("unresolved_actionable") or 0
            )
            missing = int(paper_final_check_gate_result.get("missing_carriers") or 0)
            ungrounded = int(
                paper_final_check_gate_result.get("ungrounded_actionable") or 0
            )
            unverified = int(
                paper_final_check_gate_result.get("unverified_actionable") or 0
            )
            advisory_reasons = int(
                paper_final_check_gate_result.get("advisory_reason_flags") or 0
            )
            stale_actionable = int(
                paper_final_check_gate_result.get("stale_actionable") or 0
            )
            stale_missing = int(
                paper_final_check_gate_result.get("stale_missing_carriers") or 0
            )
            stale_reviewer = int(
                paper_final_check_gate_result.get("stale_reviewer_version") or 0
            )
            if unresolved:
                stage_failures.append(
                    f"paper final-check gate left {unresolved} actionable "
                    "flag binding(s) unresolved"
                )
            if missing:
                stage_failures.append(
                    f"paper final check found {missing} quote-verified missing "
                    "carrier group(s); re-extraction required"
                )
            if stale_actionable or stale_missing:
                stage_failures.append(
                    "paper final-check gate refused stale findings "
                    f"({stale_actionable} actionable flag(s), {stale_missing} "
                    f"missing carrier group(s), {stale_reviewer} outdated "
                    "reviewer generation(s)); "
                    "reviewer replay required"
                )
            if ungrounded:
                stage_warnings.append(
                    f"paper final-check gate kept {ungrounded} DB-only "
                    "actionable flag(s) advisory because source was unavailable"
                )
            if advisory_reasons:
                stage_warnings.append(
                    f"paper final-check gate kept {advisory_reasons} weak/"
                    "unsupported finding(s) advisory; only objective contradiction "
                    "reason codes enforce field trust"
                )
            if unverified:
                stage_failures.append(
                    f"paper final-check gate refused {unverified} actionable "
                    "finding(s) without a verbatim source-verified quote; reviewer "
                    "replay required"
                )
        except Exception as e:  # noqa: BLE001
            logger.exception("paper final-check trust composition failed (%s)", e)
            stage_failures.append(f"paper final-check trust composition failed: {e}")
    else:
        logger.info("⏭️  Step 3.9: paper final-check trust composition — SKIPPED")

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

    # Step 4.4: source-integrity audit (ground truth). The preflight guard is
    # point-in-time; this catches a session that expired mid-run or a per-publisher
    # gap by measuring how many resolved papers actually landed as full text vs.
    # abstract-only stubs. Runs before corpus sync (while run_dir holds the run's
    # own FULL_CONTEXT files) and only escalates to a stage failure on a
    # full-recovery run at near-total failure (a healthy run is ~30-50% full text).
    integrity_status: Optional[dict] = None
    try:
        from cli.institutional_preflight import audit_source_integrity

        integ = audit_source_integrity(run_dir)
        integrity_status = {
            "full_text": integ.full_text,
            "abstract_only": integ.abstract_only,
            "total": integ.total,
            "abstract_only_ratio": round(integ.ratio, 3),
            "threshold": integ.threshold,
            "degraded": bool(integ.degraded and institutional_run),
        }
        if integ.total and integrity_status["degraded"]:
            # Full-recovery run that came back near-empty of full text: escalate.
            logger.warning("🧪 %s", integ.message)
            stage_failures.append(integ.message)
        elif integ.total and integ.degraded:
            # High stub ratio, but this run deliberately skipped recovery
            # (fast/calibration) — that's expected, so don't cry "fix access".
            logger.info(
                "🧪 source integrity: %d/%d papers abstract-only "
                "(expected for a fast/calibration run)",
                integ.abstract_only,
                integ.total,
            )
        elif integ.total:
            logger.info("🧪 %s", integ.message)
    except Exception as e:  # noqa: BLE001 - audit is best-effort, never fatal
        logger.warning("source-integrity audit failed (%s); continuing", e)

    # Step 4.5: fold newly-fetched source into the consolidated corpus cache
    # (idempotent incremental merge — adds new papers / upgrades compromised
    # categories only, so the next run reuses them instead of re-fetching).
    if corpus_sync and "corpus-sync" not in skip:
        step_corpus_sync(run_dir=run_dir, stage_warnings=stage_warnings, gene=gene)
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
    if publish_review and "publish-review" not in skip and stage_failures:
        logger.error(
            "🛑 Step 6: publish REFUSED — %d required stage failure(s)",
            len(stage_failures),
        )
    elif publish_review and "publish-review" not in skip:
        logger.info("📤 Step 6: publish to review DB")
        step_publish_review(
            gene=gene,
            db=db,
            disease=disease,
            review_repo=review_repo,
            pmid_file=pmid_file,
            timeout_s=publish_timeout_s,
        )
    elif publish_review and "publish-review" in skip:
        logger.info("⏭️  Step 6: publish to review DB — SKIPPED")

    # Rebuild after all gvf-run stages, because count recovery, claim checks, and
    # other post-extraction reviewers may add calls after the inner workflow
    # wrote its first manifest.
    trace_summary = finalize_trace_artifacts(trace_session, gene, stage_warnings)

    exit_code = EXIT_STAGE_WARNINGS if stage_failures else 0
    status = "completed_with_warnings" if stage_failures else "completed"
    _write_run_status(
        run_dir,
        gene,
        status,
        exit_code,
        stage_failures,
        stage_warnings,
        started,
        db,
        integrity_status,
        count_recovery=count_recovery_status,
        llm_trace=trace_summary,
        gold_free_run=gold_free_run,
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
