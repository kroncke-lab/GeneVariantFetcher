#!/usr/bin/env python3
"""Score a post-lock source overlay as an explicitly unblinded diagnostic."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from benchmarks.codex_paper_eval.apply_source_backed_overlay import (
    apply_overlay,
    digest,
    load_json,
)
from benchmarks.codex_paper_eval.run_eval import (
    CARDIAC_GENES,
    aggregate,
    gold_csv_path,
    load_gold,
    score_one,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--overlay", type=Path, required=True)
    parser.add_argument("--gold-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="omit per-gene and per-paper details for a compact checked artifact",
    )
    args = parser.parse_args()

    overlay = load_json(args.overlay)
    baseline_digest = digest(args.predictions)
    if baseline_digest != overlay.get("baseline_predictions_sha256"):
        raise SystemExit("baseline prediction digest does not match overlay")
    selection = load_json(args.selection)
    predictions = load_json(args.predictions)
    apply_overlay(predictions, overlay)
    prediction_map = {
        (paper["gene"], str(paper["pmid"])): paper for paper in predictions["papers"]
    }
    scores = []
    for paper in selection["papers"]:
        gene, pmid = paper["gene"], str(paper["pmid"])
        scores.append(
            score_one(
                gene,
                pmid,
                prediction_map[(gene, pmid)],
                load_gold(args.gold_root, gene, pmid),
            )
        )
    present_genes = [
        gene for gene in CARDIAC_GENES if any(row["gene"] == gene for row in scores)
    ]
    report = {
        "schema_version": 1,
        "status": "post_lock_source_backed_diagnostic",
        "blinded": False,
        "gold_used_only_after_overlay": True,
        "baseline_predictions_sha256": baseline_digest,
        "selection_sha256": digest(args.selection),
        "overlay_sha256": digest(args.overlay),
        "gold_sha256": {
            gene: digest(gold_csv_path(args.gold_root, gene)) for gene in present_genes
        },
        "overall": aggregate(scores),
    }
    if not args.summary_only:
        report["by_gene"] = {
            gene: aggregate([row for row in scores if row["gene"] == gene])
            for gene in present_genes
        }
        report["papers"] = scores
    args.out.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
