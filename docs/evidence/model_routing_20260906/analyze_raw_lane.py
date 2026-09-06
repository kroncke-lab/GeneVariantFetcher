"""Unvalidated clinical proposals on the same locked reference rows, after score."""

import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from benchmarks.codex_paper_eval.run_eval import merge_notation_twins
from docs.evidence.tranche_validation_20260905.summarize import count_summary, key

HERE = Path(__file__).resolve().parent
RUN = (
    ROOT
    / "benchmarks/codex_paper_eval/runs/20260906_model12_grok43_astra_clinical_verified"
)
FIELDS = ("carriers", "affected", "unaffected")


def main():
    lock = json.loads((RUN / "LOCK.json").read_text())
    assert (
        hashlib.sha256((RUN / "predictions.json").read_bytes()).hexdigest()
        == lock["predictions_sha256"]
    )
    predictions = json.loads((RUN / "predictions.json").read_text())
    report = json.loads((RUN / "report.json").read_text())
    lane = report["provenance_lane_scores"]["reader_raw_additive"]
    primary = {key(p): p for p in report["papers"]}
    extras = dict.fromkeys(FIELDS, 0)
    by_paper = {}
    for paper in predictions["papers"]:
        rows, _ = merge_notation_twins(
            paper["comparison_variants"]["reader_raw_additive"], paper["gene"]
        )
        by_paper[key(paper)] = {r["variant"]: r for r in rows}
    for scored in lane["papers"]:
        assert scored["matched_variants"] == primary[key(scored)]["matched_variants"]
        assert scored["extra_predictions"] == primary[key(scored)]["extra_predictions"]
        for name in scored["extra_predictions"]:
            row = by_paper[key(scored)][name]
            for f in FIELDS:
                extras[f] += row.get(f) is not None
    with (RUN / "figures/data/gold_difference.csv").open() as handle:
        rows = [r for r in csv.DictReader(handle) if r["analysis_run"] == RUN.name]
    for r in rows:
        name = r["predicted_variant"]
        raw = by_paper[key(r)][name].get(r["measure"]) if name else None
        evaluated = raw if raw is not None else 0
        difference = evaluated - float(r["gold_count"])
        r.update(
            automated_count_raw=raw,
            automated_count_evaluated=evaluated,
            difference=difference,
            absolute_difference=abs(difference),
            exact=int(difference == 0),
            status="identity_miss"
            if not name
            else "abstained"
            if raw is None
            else "supplied",
        )
    result = {
        "classification": "Unvalidated additive proposals; not accepted predictions, source truth or a promotion result",
        "identities_and_matching_unchanged_from_primary": True,
        "overall": lane["overall"],
        "counts": {
            f: count_summary([r for r in rows if r["measure"] == f]) for f in FIELDS
        },
        "non_null_count_fields_on_identity_extras": extras,
        "by_paper": {
            ":".join(k): {
                f: count_summary([r for r in rows if key(r) == k and r["measure"] == f])
                for f in FIELDS
            }
            for k in by_paper
        },
        "without_both_preflagged_disputes": {
            f: count_summary(
                [
                    r
                    for r in rows
                    if r["measure"] == f and r["pmid"] not in {"20129283", "25814417"}
                ]
            )
            for f in FIELDS
        },
    }
    (HERE / "raw_lane_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in result.items() if k not in ["by_paper", "overall"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
