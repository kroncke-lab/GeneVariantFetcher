#!/usr/bin/env python3
"""Drive the post-insttoken re-extraction experiment end-to-end.

Steps per gene (KCNH2, KCNQ1, RYR2, SCN5A):

1. Read the insttoken unlock CSV and select PMIDs where outcome == "FULLTEXT".
2. Delete the matching ``extractions/{GENE}_PMID_{pmid}.json`` files so the
   extraction step actually re-runs them. Every other PMID's extraction JSON
   is preserved, so the workflow only spends LLM budget on the unlocked set.
3. Run ``gvf gvf-run <GENE> --resume-dir <existing>`` so harvest steps that
   already completed are reused. The new bodies on disk flow into the DB.
4. Score the freshly-migrated SQLite against the per-gene gold standard via
   ``scripts/run_recall_suite.py``.
5. Compute MAE versus the pre-experiment baseline summary.

Run this from the repo root with the project venv activated::

    .venv/bin/python scripts/run_insttoken_reextract_experiment.py \
        --email you@example.com

Pass ``--dry-run`` to print the plan without deleting or running anything.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]

logger = logging.getLogger("insttoken_experiment")


@dataclass(frozen=True)
class GenePlan:
    gene: str
    run_dir: Path  # resume_dir
    output_root: Path  # parent of <GENE>/<timestamp>/
    unlock_csv: Path
    db_path: Path  # gvf-run writes the DB here after migration

    @property
    def extractions_dir(self) -> Path:
        return self.run_dir / "extractions"


# Maps from RECALL_STATUS.md "Where the unlocked bodies live" + per-gene run dirs.
PLANS: tuple[GenePlan, ...] = (
    GenePlan(
        gene="KCNH2",
        run_dir=REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/KCNH2/20260518_213938",
        output_root=REPO_ROOT / "validation_runs/turnkey_e2e_20260518_213934/results",
        unlock_csv=REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/KCNH2/20260518_213938/pmc_fulltext/insttoken_unlock_results.csv",
        db_path=REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/KCNH2/20260518_213938/KCNH2.db",
    ),
    GenePlan(
        gene="KCNQ1",
        run_dir=REPO_ROOT
        / "validation_runs/20260517_203904/results/KCNQ1/20260517_204424",
        output_root=REPO_ROOT / "validation_runs/20260517_203904/results",
        unlock_csv=REPO_ROOT
        / "validation_runs/20260517_203904/results/KCNQ1/20260517_204424/pmc_fulltext/insttoken_unlock_results.csv",
        db_path=REPO_ROOT
        / "validation_runs/20260517_203904/results/KCNQ1/20260517_204424/KCNQ1.db",
    ),
    GenePlan(
        gene="RYR2",
        run_dir=REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938",
        output_root=REPO_ROOT / "validation_runs/turnkey_e2e_20260518_213934/results",
        unlock_csv=REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938/pmc_fulltext/insttoken_unlock_results.csv",
        db_path=REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938/RYR2.db",
    ),
    GenePlan(
        gene="SCN5A",
        run_dir=REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938",
        output_root=REPO_ROOT / "validation_runs/turnkey_e2e_20260518_213934/results",
        unlock_csv=REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938/pmc_fulltext/insttoken_unlock_results.csv",
        db_path=REPO_ROOT
        / "validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938/SCN5A.db",
    ),
)

DEFAULT_BASELINE_SUMMARY = (
    REPO_ROOT / "recall_metrics/post_insttoken_20260521/summary.json"
)


def _setup_logging(log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    fmt = logging.Formatter(
        "%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    root.handlers = [sh, fh]


def _unlocked_pmids(unlock_csv: Path) -> list[str]:
    if not unlock_csv.exists():
        raise FileNotFoundError(f"Missing unlock CSV: {unlock_csv}")
    pmids: list[str] = []
    with unlock_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            outcome = (row.get("outcome") or "").strip().upper()
            pmid = (row.get("PMID") or row.get("pmid") or "").strip()
            if outcome == "FULLTEXT" and pmid.isdigit():
                pmids.append(pmid)
    return pmids


def _delete_extraction_jsons(
    plan: GenePlan,
    pmids: Iterable[str] | None,
    *,
    dry_run: bool,
    full_reextract: bool,
) -> list[Path]:
    """Move extraction JSONs into a backup dir so they will be re-run.

    When ``full_reextract`` is True, *every* ``*.json`` file in the gene's
    ``extractions/`` directory is moved. Otherwise only the per-PMID files
    matching ``pmids`` are moved.
    """
    deleted: list[Path] = []
    if not plan.extractions_dir.exists():
        logger.warning(
            "No extractions dir for %s at %s; skipping deletes",
            plan.gene,
            plan.extractions_dir,
        )
        return deleted
    label = "full_reextract" if full_reextract else "pre_insttoken"
    backup_root = (
        plan.extractions_dir.parent
        / f"extractions_{label}_backup_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
    )

    if full_reextract:
        targets = sorted(plan.extractions_dir.glob("*.json"))
    else:
        targets = []
        for pmid in pmids or ():
            candidate = plan.extractions_dir / f"{plan.gene}_PMID_{pmid}.json"
            if candidate.exists():
                targets.append(candidate)

    for target in targets:
        if dry_run:
            deleted.append(target)
            continue
        backup_root.mkdir(parents=True, exist_ok=True)
        shutil.move(str(target), str(backup_root / target.name))
        deleted.append(target)

    if deleted and not dry_run:
        logger.info(
            "Moved %d extraction JSON(s) for %s to %s",
            len(deleted),
            plan.gene,
            backup_root,
        )
    elif deleted:
        logger.info(
            "[dry-run] would move %d extraction JSON(s) for %s",
            len(deleted),
            plan.gene,
        )
    return deleted


def _gvf_run_cmd(
    plan: GenePlan, email: str, *, skip_layers: bool, max_pmids: int
) -> list[str]:
    cmd = [
        sys.executable,
        "-m",
        "cli",
        "gvf-run",
        plan.gene,
        "--email",
        email,
        "--output",
        str(plan.output_root),
        "--resume-dir",
        str(plan.run_dir),
        "--max-pmids",
        str(max_pmids),
    ]
    if skip_layers:
        cmd.extend(["--skip", "layers"])
    return cmd


def _run_gvf_run(
    plan: GenePlan,
    email: str,
    *,
    skip_layers: bool,
    dry_run: bool,
    log_dir: Path,
    max_pmids: int,
) -> int:
    cmd = _gvf_run_cmd(plan, email, skip_layers=skip_layers, max_pmids=max_pmids)
    logger.info("$ %s", " ".join(cmd))
    if dry_run:
        return 0
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"gvf_run_{plan.gene}.log"
    started = time.time()
    with log_path.open("w") as logfile:
        proc = subprocess.run(
            cmd, cwd=str(REPO_ROOT), stdout=logfile, stderr=subprocess.STDOUT
        )
    elapsed = time.time() - started
    logger.info(
        "gvf-run %s exit=%d elapsed=%.1f min log=%s",
        plan.gene,
        proc.returncode,
        elapsed / 60.0,
        log_path,
    )
    return proc.returncode


def _run_gvf_runs_parallel(
    plans: list[GenePlan],
    email: str,
    *,
    skip_layers: bool,
    dry_run: bool,
    log_dir: Path,
    max_pmids: int,
) -> dict[str, int]:
    """Launch all gvf-runs concurrently; wait for all; return per-gene exit codes."""
    log_dir.mkdir(parents=True, exist_ok=True)
    procs: list[tuple[GenePlan, subprocess.Popen, Path, float]] = []
    for plan in plans:
        cmd = _gvf_run_cmd(plan, email, skip_layers=skip_layers, max_pmids=max_pmids)
        logger.info("$ [parallel] %s", " ".join(cmd))
        if dry_run:
            continue
        log_path = log_dir / f"gvf_run_{plan.gene}.log"
        logfile = log_path.open("w")
        started = time.time()
        # start_new_session=True puts each child in its own process group,
        # so a kill on the parent doesn't propagate. Combined with running
        # this entire script under `screen -dmS`, the children survive
        # parent-shell death.
        proc = subprocess.Popen(
            cmd,
            cwd=str(REPO_ROOT),
            stdout=logfile,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            start_new_session=True,
        )
        procs.append((plan, proc, log_path, started))
        logger.info("  launched %s as pid %d log=%s", plan.gene, proc.pid, log_path)

    if dry_run:
        return {p.gene: 0 for p in plans}

    exit_codes: dict[str, int] = {}
    # Poll loop so we can log per-gene completion as it happens.
    pending = list(procs)
    while pending:
        time.sleep(30)
        still_pending: list[tuple[GenePlan, subprocess.Popen, Path, float]] = []
        for plan, proc, log_path, started in pending:
            rc = proc.poll()
            if rc is None:
                still_pending.append((plan, proc, log_path, started))
                continue
            elapsed = time.time() - started
            logger.info(
                "gvf-run %s exit=%d elapsed=%.1f min log=%s",
                plan.gene,
                rc,
                elapsed / 60.0,
                log_path,
            )
            exit_codes[plan.gene] = rc
        pending = still_pending
    return exit_codes


def _score(plans: list[GenePlan], outdir: Path, *, dry_run: bool, log_dir: Path) -> int:
    db_args = []
    for plan in plans:
        db_args.extend(["--db", f"{plan.gene}={plan.db_path}"])
    cmd = [
        sys.executable,
        "scripts/run_recall_suite.py",
        "--score",
        "--genes",
        ",".join(p.gene for p in plans),
        "--outdir",
        str(outdir),
        *db_args,
    ]
    logger.info("$ %s", " ".join(cmd))
    if dry_run:
        return 0
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "score.log"
    with log_path.open("w") as logfile:
        proc = subprocess.run(
            cmd, cwd=str(REPO_ROOT), stdout=logfile, stderr=subprocess.STDOUT
        )
    logger.info("score exit=%d log=%s", proc.returncode, log_path)
    return proc.returncode


def _mae_rows(summary_json: Path, outdir: Path, *, dry_run: bool) -> int:
    cmd = [
        sys.executable,
        "scripts/recall_mae.py",
        "rows",
        "--summary",
        str(summary_json),
        "--outdir",
        str(outdir),
    ]
    logger.info("$ %s", " ".join(cmd))
    if dry_run:
        return 0
    proc = subprocess.run(cmd, cwd=str(REPO_ROOT))
    return proc.returncode


def _mae_runs(before: Path, after: Path, outdir: Path, *, dry_run: bool) -> int:
    cmd = [
        sys.executable,
        "scripts/recall_mae.py",
        "runs",
        "--before",
        str(before),
        "--after",
        str(after),
        "--outdir",
        str(outdir),
    ]
    logger.info("$ %s", " ".join(cmd))
    if dry_run:
        return 0
    proc = subprocess.run(cmd, cwd=str(REPO_ROOT))
    return proc.returncode


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--email",
        required=True,
        help="Email for NCBI E-utilities; passed to gvf gvf-run.",
    )
    parser.add_argument(
        "--genes",
        default="KCNH2,KCNQ1,RYR2,SCN5A",
        help="Comma-separated genes to run (default: all four with gold).",
    )
    parser.add_argument(
        "--baseline-summary",
        type=Path,
        default=DEFAULT_BASELINE_SUMMARY,
        help="summary.json to use as the 'before' baseline for runs-mode MAE.",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Run label suffix for output directories (default: timestamp).",
    )
    parser.add_argument(
        "--skip-layers",
        action="store_true",
        help="Pass --skip layers to gvf-run (skip recovery enrichment; measures pure extraction lift).",
    )
    parser.add_argument(
        "--full-reextract",
        action="store_true",
        help=(
            "Re-extract every PMID per gene, not just the insttoken-unlocked set. "
            "Backs up *all* extraction JSONs so gvf-run regenerates them from "
            "the (post-insttoken) full text on disk."
        ),
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help=(
            "Launch all per-gene gvf-runs concurrently (independent process "
            "groups). Use when wall-clock matters more than easy debugging."
        ),
    )
    parser.add_argument(
        "--max-pmids",
        type=int,
        default=1500,
        help=(
            "Per-source PMID discovery cap passed to gvf-run. Default 1500 is "
            "the gvf-run default. Raise this when gold PMIDs are being missed: "
            "an audit showed only ~40 percent of SCN5A gold is captured at 1500."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print plan; do not delete files or run any subcommands.",
    )
    args = parser.parse_args(argv)

    selected_genes = {g.strip().upper() for g in args.genes.split(",") if g.strip()}
    plans = [p for p in PLANS if p.gene in selected_genes]
    if not plans:
        print(f"No matching genes in {args.genes}", file=sys.stderr)
        return 2

    label = args.label or datetime.now(timezone.utc).strftime(
        "post_reextract_%Y%m%dT%H%M%SZ"
    )
    score_outdir = REPO_ROOT / "recall_metrics" / label
    mae_rows_outdir = REPO_ROOT / "recall_metrics" / "mae" / f"{label}_rows"
    mae_runs_outdir = REPO_ROOT / "recall_metrics" / "mae" / f"{label}_runs"
    log_dir = score_outdir / "logs"

    _setup_logging(log_dir / "experiment.log")
    logger.info("=" * 72)
    logger.info("Insttoken re-extraction experiment label=%s", label)
    logger.info("Genes: %s", ", ".join(p.gene for p in plans))
    logger.info("Baseline summary: %s", args.baseline_summary)
    logger.info("Score outdir: %s", score_outdir)
    logger.info("MAE rows outdir: %s", mae_rows_outdir)
    logger.info("MAE runs outdir: %s", mae_runs_outdir)
    logger.info("Dry run: %s", args.dry_run)
    logger.info("=" * 72)

    # Sanity-check inputs
    for plan in plans:
        for path in (plan.run_dir, plan.unlock_csv):
            if not path.exists():
                logger.error("Missing path for %s: %s", plan.gene, path)
                return 3

    # Step 1: identify unlocked PMIDs + delete (back up) extractions per gene
    total_deleted = 0
    per_gene_unlocked: dict[str, list[str]] = {}
    for plan in plans:
        unlocked = _unlocked_pmids(plan.unlock_csv)
        per_gene_unlocked[plan.gene] = unlocked
        deleted = _delete_extraction_jsons(
            plan,
            unlocked,
            dry_run=args.dry_run,
            full_reextract=args.full_reextract,
        )
        logger.info(
            "%s: %d unlocked PMIDs; %d extraction JSON(s) queued for re-extract%s",
            plan.gene,
            len(unlocked),
            len(deleted),
            " (FULL re-extract)" if args.full_reextract else "",
        )
        total_deleted += len(deleted)
    logger.info("Total extraction JSONs queued for re-extract: %d", total_deleted)

    # Step 2: per-gene gvf-run (sequential or parallel)
    if args.parallel:
        exit_codes = _run_gvf_runs_parallel(
            plans,
            email=args.email,
            skip_layers=args.skip_layers,
            dry_run=args.dry_run,
            log_dir=log_dir,
            max_pmids=args.max_pmids,
        )
        failed = {g: rc for g, rc in exit_codes.items() if rc != 0}
        if failed:
            logger.error(
                "Parallel gvf-runs had failures: %s. Continuing to scoring "
                "with whatever DBs landed so we still get a partial report.",
                failed,
            )
    else:
        for plan in plans:
            rc = _run_gvf_run(
                plan,
                email=args.email,
                skip_layers=args.skip_layers,
                dry_run=args.dry_run,
                log_dir=log_dir,
                max_pmids=args.max_pmids,
            )
            if rc != 0:
                logger.error(
                    "gvf-run failed for %s (exit=%d); aborting downstream steps.",
                    plan.gene,
                    rc,
                )
                return rc

    # Step 3: score against gold
    rc = _score(plans, score_outdir, dry_run=args.dry_run, log_dir=log_dir)
    if rc != 0:
        logger.error("Scoring failed (exit=%d)", rc)
        return rc

    new_summary = score_outdir / "summary.json"

    # Step 4: MAE (both modes)
    if not args.dry_run and not new_summary.exists():
        logger.error("Expected scored summary not found at %s", new_summary)
        return 4

    rc = _mae_rows(new_summary, mae_rows_outdir, dry_run=args.dry_run)
    if rc != 0:
        logger.error("MAE rows mode failed (exit=%d)", rc)
        return rc

    rc = _mae_runs(
        args.baseline_summary, new_summary, mae_runs_outdir, dry_run=args.dry_run
    )
    if rc != 0:
        logger.error("MAE runs mode failed (exit=%d)", rc)
        return rc

    # Final manifest
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "label": label,
        "genes": list(per_gene_unlocked.keys()),
        "per_gene_unlocked_pmids": per_gene_unlocked,
        "total_queued_re_extractions": total_deleted,
        "baseline_summary": str(args.baseline_summary),
        "new_summary": str(new_summary),
        "score_outdir": str(score_outdir),
        "mae_rows_outdir": str(mae_rows_outdir),
        "mae_runs_outdir": str(mae_runs_outdir),
        "skip_layers": args.skip_layers,
        "full_reextract": args.full_reextract,
        "max_pmids": args.max_pmids,
    }
    manifest_path = score_outdir / "experiment_manifest.json"
    if not args.dry_run:
        manifest_path.write_text(json.dumps(manifest, indent=2))
    logger.info("Wrote experiment manifest to %s", manifest_path)
    logger.info("DONE label=%s", label)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
