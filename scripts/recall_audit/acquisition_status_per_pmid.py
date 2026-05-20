#!/usr/bin/env python3
"""Classify acquisition status for gold-standard PMIDs."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

try:
    from scripts.recall_audit.common import (
        DEFAULT_CONTEXT_MIN_BYTES,
        DEFAULT_GOLD_DIR,
        DEFAULT_RESULTS_DIR,
        context_status,
        find_full_contexts,
        load_gold_rows,
        normalize_pmid,
        repo_path,
        write_csv_rows,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import (
        DEFAULT_CONTEXT_MIN_BYTES,
        DEFAULT_GOLD_DIR,
        DEFAULT_RESULTS_DIR,
        context_status,
        find_full_contexts,
        load_gold_rows,
        normalize_pmid,
        repo_path,
        write_csv_rows,
    )


def infer_genes(gold_dir: Path) -> list[str]:
    return sorted(
        path.name.split("_recall_input.csv")[0]
        for path in gold_dir.glob("*_recall_input.csv")
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene", action="append", help="Restrict to one gene; repeatable"
    )
    parser.add_argument("--gold-dir", default=str(DEFAULT_GOLD_DIR))
    parser.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR))
    parser.add_argument("--min-bytes", type=int, default=DEFAULT_CONTEXT_MIN_BYTES)
    parser.add_argument("--out", help="Write CSV here instead of stdout")
    parser.add_argument(
        "--summary", action="store_true", help="Print status counts to stderr"
    )
    args = parser.parse_args()

    gold_dir = repo_path(args.gold_dir)
    results_dir = repo_path(args.results_dir)
    genes = [gene.upper() for gene in args.gene] if args.gene else infer_genes(gold_dir)

    rows: list[dict[str, object]] = []
    status_counts: Counter[str] = Counter()
    for gene in genes:
        pmids = sorted(
            {
                normalize_pmid(row.get("pmid"))
                for row in load_gold_rows(gene, None, gold_dir)
            }
        )
        for pmid in pmids:
            if not pmid:
                continue
            contexts = find_full_contexts(results_dir, gene, pmid)
            chosen = contexts[0] if contexts else None
            status = context_status(chosen, args.min_bytes)
            status_counts[status] += 1
            rows.append(
                {
                    "gene": gene,
                    "pmid": pmid,
                    "status": status,
                    "context_path": str(chosen) if chosen else "",
                    "context_bytes": chosen.stat().st_size if chosen else "",
                    "candidate_contexts": len(contexts),
                }
            )

    write_csv_rows(
        rows,
        [
            "gene",
            "pmid",
            "status",
            "context_path",
            "context_bytes",
            "candidate_contexts",
        ],
        Path(args.out) if args.out else None,
    )
    if args.summary:
        import sys

        print(dict(status_counts), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
