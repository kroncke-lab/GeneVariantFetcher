#!/usr/bin/env python3
"""List variants with many penetrance rows in the same PMID."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

try:
    from scripts.recall_audit.common import (
        canonical_variant,
        display_variant,
        query_pipeline_rows,
        resolve_gene_db,
        write_csv_rows,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import (
        canonical_variant,
        display_variant,
        query_pipeline_rows,
        resolve_gene_db,
        write_csv_rows,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--db", help="Explicit SQLite DB path")
    parser.add_argument("--run-dir", help="Validation run directory containing dbs/")
    parser.add_argument("--pmid", help="Restrict to one PMID")
    parser.add_argument("--min-rows", type=int, default=4)
    parser.add_argument("--out", help="Write CSV here instead of stdout")
    args = parser.parse_args()

    gene = args.gene.upper()
    rows = query_pipeline_rows(resolve_gene_db(gene, args.db, args.run_dir))
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        if args.pmid and str(row.get("pmid") or "") != str(args.pmid):
            continue
        variant = canonical_variant(display_variant(row))
        if variant:
            grouped[(str(row.get("pmid") or ""), variant)].append(row)

    output: list[dict[str, object]] = []
    for (pmid, variant), group_rows in grouped.items():
        if len(group_rows) < args.min_rows:
            continue
        for row in group_rows:
            output.append(
                {
                    "gene": gene,
                    "pmid": pmid,
                    "normalized_variant": variant,
                    "group_rows": len(group_rows),
                    "penetrance_id": row.get("penetrance_id"),
                    "protein_notation": row.get("protein_notation") or "",
                    "cdna_notation": row.get("cdna_notation") or "",
                    "genomic_position": row.get("genomic_position") or "",
                    "carriers": row.get("total_carriers_observed") or "",
                    "affected": row.get("affected_count") or "",
                    "unaffected": row.get("unaffected_count") or "",
                    "source_location": row.get("source_location") or "",
                }
            )

    output.sort(
        key=lambda row: (
            -int(row["group_rows"]),
            row["pmid"],
            row["normalized_variant"],
        )
    )
    write_csv_rows(
        output,
        [
            "gene",
            "pmid",
            "normalized_variant",
            "group_rows",
            "penetrance_id",
            "protein_notation",
            "cdna_notation",
            "genomic_position",
            "carriers",
            "affected",
            "unaffected",
            "source_location",
        ],
        Path(args.out) if args.out else None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
