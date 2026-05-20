#!/usr/bin/env python3
"""Show gold-standard and pipeline rows for one PMID side by side."""

from __future__ import annotations

import argparse
from pathlib import Path

try:
    from scripts.recall_audit.common import (
        canonical_variant,
        display_variant,
        format_table,
        load_gold_rows,
        parse_int,
        query_pipeline_rows,
        resolve_gene_db,
    )
except ModuleNotFoundError:  # pragma: no cover - direct script execution
    from common import (
        canonical_variant,
        display_variant,
        format_table,
        load_gold_rows,
        parse_int,
        query_pipeline_rows,
        resolve_gene_db,
    )


def balance_flag(row: dict) -> str:
    carriers = parse_int(row.get("total_carriers_observed"))
    affected = parse_int(row.get("affected_count"))
    unaffected = parse_int(row.get("unaffected_count"))
    if carriers is None or affected is None or unaffected is None:
        return ""
    return "imbalance" if affected + unaffected != carriers else ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True, help="Gene symbol, e.g. KCNH2")
    parser.add_argument("--pmid", required=True, help="PMID to inspect")
    parser.add_argument("--run-dir", help="Validation run directory containing dbs/")
    parser.add_argument("--db", help="Explicit SQLite DB path")
    parser.add_argument("--gold-dir", help="Gold-standard normalized directory")
    args = parser.parse_args()

    gene = args.gene.upper()
    db_path = resolve_gene_db(gene, args.db, args.run_dir)
    gold_rows = load_gold_rows(gene, args.pmid, args.gold_dir)
    sqlite_rows = query_pipeline_rows(Path(db_path), args.pmid)

    print(f"# {gene} PMID {args.pmid}")
    print(f"Gold rows: {len(gold_rows)}")
    print(f"Pipeline penetrance rows: {len(sqlite_rows)}")
    print()

    gold_table = [
        [
            row.get("variant", ""),
            row.get("carriers", ""),
            row.get("affected", ""),
            row.get("unaffected", ""),
        ]
        for row in gold_rows
    ]
    print("## Gold")
    print(format_table(["variant", "carriers", "affected", "unaffected"], gold_table))
    print()

    sqlite_table = [
        [
            display_variant(row),
            row.get("cdna_notation") or "",
            row.get("total_carriers_observed") or "",
            row.get("affected_count") or "",
            row.get("unaffected_count") or "",
            balance_flag(row),
            row.get("source_location") or "",
        ]
        for row in sqlite_rows
    ]
    print("## Pipeline")
    print(
        format_table(
            ["variant", "cdna", "carriers", "affected", "unaffected", "flag", "source"],
            sqlite_table,
        )
    )
    print()

    gold_keys = {canonical_variant(row.get("variant")) for row in gold_rows}
    sqlite_keys = {canonical_variant(display_variant(row)) for row in sqlite_rows}
    gold_keys.discard("")
    sqlite_keys.discard("")

    only_gold = sorted(gold_keys - sqlite_keys)
    only_sqlite = sorted(sqlite_keys - gold_keys)
    imbalances = [row for row in sqlite_rows if balance_flag(row)]

    print("## Differences")
    print(f"Gold-only variants ({len(only_gold)}): {', '.join(only_gold) or '-'}")
    print(
        f"Pipeline-only variants ({len(only_sqlite)}): {', '.join(only_sqlite) or '-'}"
    )
    if imbalances:
        print(f"Rows where affected + unaffected != carriers: {len(imbalances)}")
    else:
        print("Rows where affected + unaffected != carriers: 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
