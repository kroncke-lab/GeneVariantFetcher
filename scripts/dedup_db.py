#!/usr/bin/env python3
"""Collapse exact-duplicate child rows in a GVF SQLite DB (idempotent).

Exact-duplicate ``penetrance_data`` / ``individual_records`` / ``functional_data``
/ ``phenotypes`` / ``variant_metadata`` rows accumulate when one variant is
represented across several supplement-table cells (multi-cohort /
multi-classification tables) or when a paper is re-migrated on a refresh. Because
the recall scorer SUMS ``total_carriers_observed`` across the rows linked to a
(pmid, variant), these duplicates inflate the carrier/affected counts N-fold —
e.g. KCNQ1 PMID 32893267 V254M scored 4x gold (100 vs 25).

This script removes those exact duplicates (keeping one of each identical set);
genuinely different rows (real sub-cohort splits) are left alone. It backs up the
DB first and is safe to re-run — a clean DB yields zero removals.

The matching write-path fix is the exact-dup insert guards in
``harvesting.migrate_to_sqlite.insert_variant_data`` (new duplicates are no longer
written); this script back-fills DBs built before those guards existed.

Typical usage:

    .venv/bin/python scripts/dedup_db.py --db path/to/GENE.db

After dedup, re-score with scripts/run_recall_suite.py to measure the MAE delta.
"""

from __future__ import annotations

import argparse
import shutil
import sqlite3
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from harvesting.migrate_to_sqlite import dedup_existing_rows  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", required=True, type=Path, help="SQLite DB to dedup.")
    ap.add_argument(
        "--no-backup",
        action="store_true",
        help="Skip the .before_dedup.db backup (not recommended).",
    )
    args = ap.parse_args(argv)

    db = args.db.expanduser().resolve()
    if not db.exists():
        ap.error(f"DB not found: {db}")

    if not args.no_backup:
        backup = db.with_name(db.name + ".before_dedup.db")
        shutil.copy2(db, backup)
        print(f"backup: {backup}")

    conn = sqlite3.connect(str(db))
    try:
        removed = dedup_existing_rows(conn)
    finally:
        conn.close()

    total = sum(removed.values())
    for table, n in removed.items():
        print(f"  {table}: -{n}")
    print(f"removed {total} exact-duplicate row(s) from {db.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
