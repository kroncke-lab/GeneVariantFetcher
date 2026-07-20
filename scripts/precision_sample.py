#!/usr/bin/env python3
"""Precision calibration sampler: turn the FP-upper-bound proxy into a real number.

This is a periodic **calibration** instrument, not a per-run review step. The
recall scorer reports ``precision_vs_counted_gold_pmids`` -- a false-positive
UPPER BOUND, not clean precision, because "extra" DB rows on gold PMIDs may be
true positives the curator simply never recorded. A one-off adjudicated sample
resolves that ambiguity and anchors the true-precision number the automated
quality gates are calibrated against -- run occasionally to set/validate those
gates, not on every extraction run:

    sample  -> draw N random counted-extra rows (SQLite-only, count-bearing, on
               gold PMIDs -- the exact population compute_precision_summary
               counts) and write a review worksheet CSV.
    score   -> read the filled worksheet back and compute TRUE precision with a
               Wilson 95% confidence interval, excluding 'unsure' verdicts.

The counted-extra population is selected with the scorer's own comparison rows,
so its definition cannot drift from ``compute_precision_summary``.

Examples
--------
Draw a 150-row sample across the four cardiac gold genes::

    python scripts/precision_sample.py sample \\
        --gold-dir gene_variant_fetcher_gold_standard/normalized \\
        --db KCNH2=results/KCNH2/latest/KCNH2.db \\
        --db KCNQ1=results/KCNQ1/latest/KCNQ1.db \\
        --db SCN5A=results/SCN5A/latest/SCN5A.db \\
        --db RYR2=results/RYR2/latest/RYR2.db \\
        --n 150 --seed 20260709 \\
        --out precision_sample.csv

Fill the ``verdict`` column (real_tp | false_positive | unsure), then::

    python scripts/precision_sample.py score --worksheet precision_sample_FILLED.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sqlite3
import sys
from pathlib import Path
from typing import Any, Optional

BASE_DIR = Path(__file__).resolve().parents[1]
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from cli.compare_variants import (  # noqa: E402
    ComparisonRow,
    aggregate_excel_data,
    aggregate_sqlite_data,
    compare_data,
    extract_sqlite_data,
    introspect_sqlite,
    load_excel_data,
)

WORKSHEET_COLUMNS = [
    "gene",
    "pmid",
    "pubmed_url",
    "variant",
    "carriers",
    "affected",
    "unaffected",
    "source_layer",
    # Reviewer fills this in: real_tp | false_positive | unsure
    "verdict",
    "reviewer_notes",
]

TP_VERDICTS = {"real_tp", "tp", "real", "true_positive", "correct", "confirm"}
FP_VERDICTS = {"false_positive", "fp", "wrong", "junk", "not_in_paper"}


def build_comparison_rows(
    gold_path: Path,
    db_path: Path,
    match_mode: str = "fuzzy",
    fuzzy_threshold: float = 0.80,
) -> list[ComparisonRow]:
    """Build the scorer's comparison rows for one gene (no files written)."""
    excel_df, detected_columns = load_excel_data(gold_path, None, None)
    conn = sqlite3.connect(f"{db_path.resolve().as_uri()}?mode=ro", uri=True)
    try:
        table_info = introspect_sqlite(conn)
        # This sample calibrates the trust gate, so make its projection explicit.
        sqlite_df = extract_sqlite_data(conn, table_info, trust_mode="trusted")
    finally:
        conn.close()
    gold_aggregated = aggregate_excel_data(excel_df, detected_columns)
    sqlite_aggregated = aggregate_sqlite_data(sqlite_df)
    return compare_data(gold_aggregated, sqlite_aggregated, match_mode, fuzzy_threshold)


def select_counted_extras(rows: list[ComparisonRow]) -> list[ComparisonRow]:
    """DB-only rows on gold PMIDs that carry a count.

    Matches ``compute_precision_summary``'s ``counted_extra_on_gold_pmids``
    population exactly: ``missing_in_excel`` and PMID in the gold PMID set and at
    least one non-None extracted count.
    """
    gold_pmids = {r.pmid for r in rows if not r.missing_in_excel and r.pmid}
    out = []
    for r in rows:
        if not r.missing_in_excel or r.pmid not in gold_pmids:
            continue
        if any(
            v is not None
            for v in (r.sqlite_carriers_total, r.sqlite_affected, r.sqlite_unaffected)
        ):
            out.append(r)
    return out


def _row_to_worksheet(gene: str, r: ComparisonRow) -> dict[str, Any]:
    return {
        "gene": gene,
        "pmid": r.pmid,
        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{r.pmid}/" if r.pmid else "",
        "variant": r.sqlite_variant_raw or r.sqlite_variant_norm or "",
        "carriers": "" if r.sqlite_carriers_total is None else r.sqlite_carriers_total,
        "affected": "" if r.sqlite_affected is None else r.sqlite_affected,
        "unaffected": "" if r.sqlite_unaffected is None else r.sqlite_unaffected,
        "source_layer": r.sqlite_source_layer or "",
        "verdict": "",
        "reviewer_notes": "",
    }


def write_worksheet(sampled: list[dict[str, Any]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=WORKSHEET_COLUMNS)
        writer.writeheader()
        writer.writerows(sampled)


def wilson_interval(k: int, n: int, z: float = 1.96) -> tuple[Optional[float], ...]:
    """Wilson score interval for k successes in n trials. Returns (p, lo, hi)."""
    if n == 0:
        return (None, None, None)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    margin = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (p, max(0.0, center - margin), min(1.0, center + margin))


def tally_verdicts(worksheet: Path) -> dict[str, Any]:
    """Read a filled worksheet and compute true precision + Wilson CI."""
    tp = fp = unsure = blank = 0
    unknown: list[str] = []
    with open(worksheet, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            verdict = (row.get("verdict") or "").strip().lower()
            if not verdict:
                blank += 1
            elif verdict in TP_VERDICTS:
                tp += 1
            elif verdict in FP_VERDICTS:
                fp += 1
            elif verdict in {"unsure", "unknown", "?"}:
                unsure += 1
            else:
                unknown.append(verdict)
    adjudicated = tp + fp
    p, lo, hi = wilson_interval(tp, adjudicated)
    return {
        "true_positives": tp,
        "false_positives": fp,
        "unsure": unsure,
        "blank": blank,
        "unrecognized_verdicts": sorted(set(unknown)),
        "adjudicated": adjudicated,
        "true_precision": p,
        "wilson_95_ci": None if p is None else [lo, hi],
    }


def _cmd_sample(args: argparse.Namespace) -> int:
    rng = random.Random(args.seed)
    pool: list[dict[str, Any]] = []
    per_gene: dict[str, int] = {}
    for spec in args.db:
        if "=" not in spec:
            print(f"--db expects GENE=/path/to.db, got: {spec}", file=sys.stderr)
            return 2
        gene, db_path = spec.split("=", 1)
        gold_path = args.gold_dir / f"{gene}_recall_input.csv"
        if not gold_path.exists():
            print(f"[{gene}] gold not found: {gold_path}", file=sys.stderr)
            return 2
        rows = build_comparison_rows(gold_path, Path(db_path))
        extras = select_counted_extras(rows)
        per_gene[gene] = len(extras)
        pool.extend(_row_to_worksheet(gene, r) for r in extras)

    print(f"counted-extra population: {len(pool)} rows ({per_gene})")
    if not pool:
        print("no counted extras to sample", file=sys.stderr)
        return 1
    n = min(args.n, len(pool))
    if n < args.n:
        print(f"population smaller than --n; sampling all {n} rows")
    sampled = rng.sample(pool, n)
    sampled.sort(key=lambda r: (r["gene"], r["pmid"]))
    write_worksheet(sampled, args.out)
    print(f"wrote {n}-row worksheet -> {args.out}")
    print("Fill the 'verdict' column (real_tp | false_positive | unsure), then run:")
    print(f"  python scripts/precision_sample.py score --worksheet {args.out}")
    return 0


def _cmd_score(args: argparse.Namespace) -> int:
    result = tally_verdicts(args.worksheet)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(result, indent=2))
    p = result["true_precision"]
    print(json.dumps(result, indent=2))
    if p is not None:
        lo, hi = result["wilson_95_ci"]
        print(
            f"\nTRUE precision on counted extras: {p:.1%} "
            f"(95% CI {lo:.1%}-{hi:.1%}, n={result['adjudicated']}); "
            f"{result['unsure']} unsure excluded."
        )
    else:
        print("\nNo adjudicated rows (fill the verdict column first).")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("sample", help="draw a random counted-extra worksheet")
    s.add_argument(
        "--db",
        action="append",
        default=[],
        metavar="GENE=PATH",
        help="Gene to DB path (repeatable).",
    )
    s.add_argument(
        "--gold-dir",
        type=Path,
        default=BASE_DIR / "gene_variant_fetcher_gold_standard" / "normalized",
        help="Dir holding <GENE>_recall_input.csv gold files.",
    )
    s.add_argument("--n", type=int, default=150, help="Sample size (default 150).")
    s.add_argument(
        "--seed", type=int, default=20260709, help="RNG seed for reproducibility."
    )
    s.add_argument("--out", type=Path, default=Path("precision_sample.csv"))
    s.set_defaults(func=_cmd_sample)

    sc = sub.add_parser("score", help="score a filled worksheet -> true precision")
    sc.add_argument("--worksheet", type=Path, required=True)
    sc.add_argument("--json-out", type=Path, default=None)
    sc.set_defaults(func=_cmd_score)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
