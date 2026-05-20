#!/usr/bin/env python3
"""Detect PMIDs where one study-wide N was smeared across many variants."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from statistics import pstdev

try:
    from scripts.recall_audit.common import (
        display_variant,
        parse_int,
        query_pipeline_rows,
        resolve_gene_db,
        write_csv_rows,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import (
        display_variant,
        parse_int,
        query_pipeline_rows,
        resolve_gene_db,
        write_csv_rows,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--db", help="Explicit SQLite DB path")
    parser.add_argument("--run-dir", help="Validation run directory containing dbs/")
    parser.add_argument("--min-variants", type=int, default=3)
    parser.add_argument("--min-carriers", type=int, default=5)
    parser.add_argument("--out", help="Write CSV here instead of stdout")
    args = parser.parse_args()

    gene = args.gene.upper()
    rows = query_pipeline_rows(resolve_gene_db(gene, args.db, args.run_dir))
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[str(row.get("pmid") or "")].append(row)

    flagged: list[dict[str, object]] = []
    for pmid, pmid_rows in grouped.items():
        variants = {display_variant(row) for row in pmid_rows if display_variant(row)}
        carrier_values = [
            parse_int(row.get("total_carriers_observed"))
            for row in pmid_rows
            if parse_int(row.get("total_carriers_observed")) is not None
        ]
        if len(variants) < args.min_variants or len(carrier_values) < args.min_variants:
            continue
        if len(set(carrier_values)) == 1 and carrier_values[0] > args.min_carriers:
            affected_values = sorted(
                {
                    parse_int(row.get("affected_count"))
                    for row in pmid_rows
                    if parse_int(row.get("affected_count")) is not None
                }
            )
            unaffected_values = sorted(
                {
                    parse_int(row.get("unaffected_count"))
                    for row in pmid_rows
                    if parse_int(row.get("unaffected_count")) is not None
                }
            )
            flagged.append(
                {
                    "gene": gene,
                    "pmid": pmid,
                    "variant_count": len(variants),
                    "row_count": len(pmid_rows),
                    "shared_carriers": carrier_values[0],
                    "carrier_stddev": f"{pstdev(carrier_values):.3f}",
                    "affected_values": ";".join(str(v) for v in affected_values),
                    "unaffected_values": ";".join(str(v) for v in unaffected_values),
                    "sample_variants": ";".join(sorted(variants)[:12]),
                }
            )

    flagged.sort(key=lambda row: (-int(row["shared_carriers"]), row["pmid"]))
    write_csv_rows(
        flagged,
        [
            "gene",
            "pmid",
            "variant_count",
            "row_count",
            "shared_carriers",
            "carrier_stddev",
            "affected_values",
            "unaffected_values",
            "sample_variants",
        ],
        Path(args.out) if args.out else None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
