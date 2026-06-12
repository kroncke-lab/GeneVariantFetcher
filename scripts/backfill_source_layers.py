#!/usr/bin/env python3
"""Backfill ``variant_papers.source_layer`` without rebuilding a GVF DB."""

from __future__ import annotations

import argparse
import shutil
import sqlite3
import sys
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from harvesting.migrate_to_sqlite import create_database_schema  # noqa: E402
from utils.source_layers import source_layer_sql_case  # noqa: E402


def _columns(con: sqlite3.Connection, table: str) -> set[str]:
    return {row[1] for row in con.execute(f"PRAGMA table_info({table})")}


def _source_layer_expr(con: sqlite3.Connection) -> str:
    columns = _columns(con, "variant_papers")
    return source_layer_sql_case(
        "source_location" if "source_location" in columns else "NULL",
        "source_layer" if "source_layer" in columns else None,
        additional_notes_expr=(
            "additional_notes" if "additional_notes" in columns else None
        ),
    )


def summarize(db_path: Path) -> list[tuple[str, int]]:
    con = sqlite3.connect(str(db_path))
    try:
        expr = _source_layer_expr(con)
        rows = con.execute(
            f"""
            SELECT {expr} AS inferred_layer, COUNT(*) AS n
            FROM variant_papers
            GROUP BY inferred_layer
            ORDER BY inferred_layer
            """
        ).fetchall()
        return [(str(layer), int(n)) for layer, n in rows]
    finally:
        con.close()


def backfill(db_path: Path, *, backup_suffix: str) -> tuple[Path, int]:
    backup = db_path.with_suffix(
        db_path.suffix + f".before_source_layer_{backup_suffix}"
    )
    shutil.copy2(db_path, backup)

    con = create_database_schema(str(db_path))
    try:
        expr = _source_layer_expr(con)
        cur = con.execute(
            f"""
            UPDATE variant_papers
            SET source_layer = {expr}
            WHERE source_layer IS NULL
               OR TRIM(source_layer) = ''
               OR LOWER(TRIM(source_layer)) = 'manual_or_legacy'
            """
        )
        con.commit()
        return backup, int(cur.rowcount if cur.rowcount is not None else 0)
    finally:
        con.close()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("db", nargs="+", type=Path, help="SQLite DB path(s)")
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Mutate DBs. Default is dry-run summary only.",
    )
    parser.add_argument(
        "--backup-suffix",
        default=datetime.now().strftime("%Y%m%d_%H%M%S"),
        help="Suffix used for .before_source_layer_* backups.",
    )
    args = parser.parse_args(argv)

    for raw_path in args.db:
        db_path = raw_path.expanduser().resolve()
        if not db_path.exists():
            raise FileNotFoundError(db_path)
        before = summarize(db_path)
        print(f"\n{db_path}")
        print("  inferred layers before/apply target:")
        for layer, count in before:
            print(f"    {layer}: {count}")
        if not args.apply:
            continue
        backup, updated = backfill(db_path, backup_suffix=args.backup_suffix)
        after = summarize(db_path)
        print(f"  backup: {backup}")
        print(f"  rows updated: {updated}")
        print("  explicit layers after:")
        for layer, count in after:
            print(f"    {layer}: {count}")
    if not args.apply:
        print("\nDry run only. Re-run with --apply to write source_layer.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
