#!/usr/bin/env python3
"""Operational driver (NOT committed): walk run_priority_extraction in windows
until per-batch variant yield tapers.

For each gene it advances --priority-offset in fixed steps, migrating each batch
so the DB reflects cumulative coverage, and stops when consecutive batches stop
adding new distinct variants. Already-extracted PMIDs are skipped by
run_priority_extraction (idempotent), so re-covering early ranks is cheap and
closes any gaps from the prior partial walk.

Yield signal = growth in COUNT(DISTINCT variant_id) in <gene>.db, which:
  * correctly keeps BRCA1 going (dense VUS tables keep adding variants), and
  * stops APOE once high-score ranks are done (it is an allele gene: deep papers
    repeat the same epsilon variants and add ~0 new distinct variants).

Run in background; tail the walk log. Safe to kill -9 and re-run (resumes).
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PRIORITY = ROOT / "scripts" / "run_priority_extraction.py"
TRIAGE_MODEL = "anthropic/claude-haiku-4-5-20251001"


def db_counts(db: Path) -> dict:
    if not db.exists():
        return {"variants": 0, "rows": 0, "papers": 0, "carriers": 0}
    con = sqlite3.connect(str(db))
    try:
        cur = con.cursor()

        def one(q: str) -> int:
            try:
                return int(cur.execute(q).fetchone()[0] or 0)
            except Exception:
                return 0

        return {
            "variants": one("SELECT COUNT(DISTINCT variant_id) FROM variant_papers"),
            "rows": one("SELECT COUNT(*) FROM variant_papers"),
            "papers": one("SELECT COUNT(*) FROM papers"),
            "carriers": one(
                "SELECT COALESCE(SUM(total_carriers_observed),0) FROM penetrance_data"
            ),
        }
    finally:
        con.close()


def log(walk_log: Path, msg: str) -> None:
    line = f"{datetime.now().isoformat(timespec='seconds')} {msg}"
    print(line, flush=True)
    with walk_log.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def walk_gene(
    gene: str,
    run_dir: Path,
    step: int,
    max_cand: int,
    min_new_variants: int,
    patience: int,
    walk_log: Path,
) -> None:
    db = run_dir / f"{gene}.db"
    summary = run_dir / "extraction_priority" / "priority_extraction_summary.json"
    log(
        walk_log,
        f"=== {gene}: walk start; run_dir={run_dir} max_cand={max_cand} "
        f"step={step} min_new_variants={min_new_variants} patience={patience}",
    )
    prev = db_counts(db)
    log(walk_log, f"{gene} baseline: {prev}")

    seen_yield = False
    low_streak = 0
    offset = 0
    while offset < max_cand:
        t0 = time.time()
        cmd = [
            sys.executable,
            str(PRIORITY),
            "--gene",
            gene,
            "--run-dir",
            str(run_dir),
            "--top-n",
            str(step),
            "--priority-offset",
            str(offset),
            "--skip-preprocess",
            # Deterministic (rule-based) triage: instant vs ~40min/batch for the
            # serial hybrid-LLM triage. It is explicitly calibrated to avoid
            # APOE-like association noise; grok extraction still gates quality.
            "--triage-mode",
            "deterministic",
        ]
        rc = subprocess.run(cmd, cwd=str(ROOT)).returncode
        dt = time.time() - t0
        now = db_counts(db)
        dV = now["variants"] - prev["variants"]
        dRows = now["rows"] - prev["rows"]
        dCar = now["carriers"] - prev["carriers"]
        extract_now = None
        try:
            extract_now = json.loads(summary.read_text())["stats"].get(
                "extraction_triage_extract_now"
            )
        except Exception:
            pass
        log(
            walk_log,
            f"{gene} offset={offset} rc={rc} {dt / 60:.1f}min "
            f"extract_now={extract_now} dVariants={dV} dRows={dRows} dCarriers={dCar} "
            f"| cum variants={now['variants']} rows={now['rows']} "
            f"papers={now['papers']} carriers={now['carriers']}",
        )

        if rc != 0:
            log(
                walk_log,
                f"{gene} batch rc={rc}; stopping this gene to avoid burning spend on a broken state",
            )
            return

        if dV >= min_new_variants:
            seen_yield = True
            low_streak = 0
        elif seen_yield:
            low_streak += 1

        if seen_yield and low_streak >= patience:
            log(
                walk_log,
                f"{gene} TAPER reached at offset={offset} "
                f"(<{min_new_variants} new variants x{patience}); stopping.",
            )
            return
        prev = now
        offset += step

    log(
        walk_log,
        f"{gene} reached candidate cap offset={offset} >= {max_cand}; full coverage walked.",
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("genes", nargs="+")
    ap.add_argument("--run-timestamp", default="20260616_132646")
    ap.add_argument("--step", type=int, default=1000)
    ap.add_argument("--min-new-variants", type=int, default=8)
    ap.add_argument("--patience", type=int, default=2)
    ap.add_argument("--walk-log", default=str(ROOT / "results" / "_priority_walk.log"))
    args = ap.parse_args()

    walk_log = Path(args.walk_log)
    walk_log.parent.mkdir(parents=True, exist_ok=True)

    for gene in args.genes:
        run_dir = ROOT / "results" / gene / args.run_timestamp
        cand_tsv = run_dir / "extraction_priority" / "priority_candidates.tsv"
        try:
            max_cand = sum(1 for _ in cand_tsv.open()) - 1
        except Exception:
            max_cand = 0
        if max_cand <= 0:
            log(walk_log, f"{gene}: no candidate ranking found at {cand_tsv}; skipping")
            continue
        walk_gene(
            gene,
            run_dir,
            args.step,
            max_cand,
            args.min_new_variants,
            args.patience,
            walk_log,
        )
    log(walk_log, "=== walk driver complete for all genes ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
