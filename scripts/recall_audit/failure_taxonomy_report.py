#!/usr/bin/env python3
"""Bucket scored recall discrepancies into actionable failure classes."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

try:
    from scripts.recall_audit.common import (
        parse_bool,
        parse_number,
        read_csv_rows,
        resolve_recall_score,
        write_csv_rows,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import (
        parse_bool,
        parse_number,
        read_csv_rows,
        resolve_recall_score,
        write_csv_rows,
    )

BUCKETS = [
    "acquisition_miss",
    "false_positive_pipeline_only",
    "notation_drift",
    "count_inflation",
    "count_deflation",
    "phenotype_inversion",
    "missing_phenotype",
    "multicohort_no_aggregation",
    "other_matched_count_mismatch",
    "other_matched",
    "other_unmatched",
]


def _matched(row: dict[str, str]) -> bool:
    return row.get("match_type") not in {"", "none", None}


def _ratio_close(
    excel: float | None, sqlite: float | None, tolerance: float = 0.10
) -> bool:
    if excel is None or sqlite is None:
        return False
    if excel == 0:
        return sqlite == 0
    return abs(excel - sqlite) / abs(excel) <= tolerance


def _multicohort_keys(rows: Iterable[dict[str, str]]) -> set[tuple[str, str]]:
    grouped: dict[tuple[str, str], list[tuple[str | None, str | None, str | None]]] = (
        defaultdict(list)
    )
    for row in rows:
        if parse_bool(row.get("missing_in_sqlite")):
            continue
        variant = row.get("sqlite_variant_norm") or row.get("sqlite_variant_raw")
        if not variant:
            continue
        key = (row.get("pmid", ""), variant)
        counts = (
            row.get("sqlite_carriers_total"),
            row.get("sqlite_affected"),
            row.get("sqlite_unaffected"),
        )
        grouped[key].append(counts)
    return {
        key
        for key, values in grouped.items()
        if len(values) >= 5 and len(set(values)) >= 2
    }


def classify(row: dict[str, str], multicohort: set[tuple[str, str]]) -> str:
    missing_sqlite = parse_bool(row.get("missing_in_sqlite"))
    missing_excel = parse_bool(row.get("missing_in_excel"))
    match_type = row.get("match_type", "")

    if match_type == "none" and missing_sqlite:
        return "acquisition_miss"
    if missing_excel:
        return "false_positive_pipeline_only"
    if match_type in {"fuzzy", "fuzzy_cdna_bridge"}:
        return "notation_drift"

    variant = row.get("sqlite_variant_norm") or row.get("sqlite_variant_raw")
    if (row.get("pmid", ""), variant) in multicohort:
        return "multicohort_no_aggregation"

    carriers_diff = parse_number(row.get("carriers_diff"))
    excel_carriers = parse_number(row.get("excel_carriers_total"))
    sqlite_carriers = parse_number(row.get("sqlite_carriers_total"))
    affected_diff = parse_number(row.get("affected_diff"))
    unaffected_diff = parse_number(row.get("unaffected_diff"))
    excel_affected = parse_number(row.get("excel_affected"))
    sqlite_affected = parse_number(row.get("sqlite_affected"))
    excel_unaffected = parse_number(row.get("excel_unaffected"))
    sqlite_unaffected = parse_number(row.get("sqlite_unaffected"))

    if (
        _matched(row)
        and carriers_diff is not None
        and carriers_diff < -5
        and excel_carriers
        and abs(carriers_diff) / excel_carriers > 1.0
    ):
        return "count_inflation"
    if (
        _matched(row)
        and carriers_diff is not None
        and carriers_diff > 5
        and excel_carriers
        and sqlite_carriers is not None
        and sqlite_carriers <= excel_carriers / 2
    ):
        return "count_deflation"
    if (
        _matched(row)
        and _ratio_close(excel_carriers, sqlite_carriers)
        and affected_diff is not None
        and unaffected_diff is not None
        and affected_diff * unaffected_diff < 0
        and abs(affected_diff) > 5
        and abs(unaffected_diff) > 5
    ):
        return "phenotype_inversion"
    if _matched(row) and (
        (excel_affected is not None and sqlite_affected is None)
        or (excel_unaffected is not None and sqlite_unaffected is None)
    ):
        return "missing_phenotype"
    if _matched(row) and parse_bool(row.get("count_mismatch")):
        return "other_matched_count_mismatch"
    if _matched(row):
        return "other_matched"
    return "other_unmatched"


def discrepancy_files(root: Path, genes: list[str] | None) -> list[Path]:
    wanted = {gene.upper() for gene in genes or []}
    files = sorted(root.glob("*/discrepancies.csv"))
    if wanted:
        files = [path for path in files if path.parent.name.upper() in wanted]
    return files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--recall-score", help="Directory containing per-gene recall_score outputs"
    )
    parser.add_argument(
        "--run-dir", help="Validation run directory containing recall_score/"
    )
    parser.add_argument(
        "--gene", action="append", help="Restrict to one gene; repeatable"
    )
    parser.add_argument("--out", help="Write CSV here instead of stdout")
    args = parser.parse_args()

    root = resolve_recall_score(args.recall_score, args.run_dir)
    output_rows: list[dict[str, object]] = []

    for path in discrepancy_files(root, args.gene):
        rows = read_csv_rows(path)
        multicohort = _multicohort_keys(rows)
        by_pmid: dict[str, Counter[str]] = defaultdict(Counter)
        for row in rows:
            bucket = classify(row, multicohort)
            by_pmid[row.get("pmid", "")][bucket] += 1

        for pmid, counts in sorted(
            by_pmid.items(), key=lambda item: (-sum(item[1].values()), item[0])
        ):
            out_row: dict[str, object] = {
                "gene": path.parent.name,
                "pmid": pmid,
                "total_rows": sum(counts.values()),
            }
            out_row.update({bucket: counts.get(bucket, 0) for bucket in BUCKETS})
            output_rows.append(out_row)

    fieldnames = ["gene", "pmid", "total_rows", *BUCKETS]
    write_csv_rows(output_rows, fieldnames, Path(args.out) if args.out else None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
