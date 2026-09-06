"""Recompute diagnostic pooled/subset effects from completed, locked pairs.

The registered verdicts remain in each candidate's compare_with_secondary.json.
These additional two-sided 95% intervals are descriptive, use PMID clusters,
and do not introduce new acceptance rules. No answer-key file is read here.
"""

import argparse
import copy
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from benchmarks.codex_paper_eval.run_eval import merge_notation_twins

OUT = Path(__file__).resolve().parent
RUNS = ROOT / "benchmarks/codex_paper_eval/runs"
FIELDS = ("carriers", "affected", "unaffected")


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key(row):
    return str(row["gene"]), str(row["pmid"])


def rate(n, d):
    return n / d if d else None


def load_arm(number, arm):
    run = RUNS / f"20260905_protocol_cont120_{number}_{arm}"
    lock = read(run / "LOCK.json")
    for name in ("selection", "predictions"):
        assert sha(run / f"{name}.json") == lock[f"{name}_sha256"]
    report = read(run / "report.json")
    papers = {key(p): p for p in report["papers"]}
    assert len(papers) == 120
    extra_counts = {}
    for paper in read(run / "predictions.json")["papers"]:
        k = key(paper)
        merged, _ = merge_notation_twins(paper.get("variants", []), k[0])
        extra_names = set(papers[k]["extra_predictions"])
        extras = [row for row in merged if row["variant"] in extra_names]
        assert len(extras) == papers[k]["fp"]
        assert (
            sum(any(row.get(field) is not None for field in FIELDS) for row in extras)
            == papers[k]["counted_precision"]["counted_extra_rows"]
        )
        extra_counts[k] = {
            field: sum(row.get(field) is not None for row in extras) for field in FIELDS
        }
    with (run / "figures/data/gold_difference.csv").open(newline="") as handle:
        # Paired figures contain both arms. Never double-count their rows.
        rows = [r for r in csv.DictReader(handle) if r["analysis_run"] == run.name]
    unique = {(r["gene"], r["pmid"], r["gold_row_index"], r["measure"]) for r in rows}
    assert len(unique) == len(rows), "duplicate count rows"
    assert all(key(r) in papers for r in rows)
    return {
        "run_id": run.name,
        "papers": papers,
        "rows": rows,
        "extra_counts": extra_counts,
        "selection": {key(p): p for p in read(run / "selection.json")["papers"]},
        "snapshot": read(run / "frozen_corpus/source_snapshot.json"),
        "sha256": {
            name: sha(run / name)
            for name in (
                "LOCK.json",
                "selection.json",
                "predictions.json",
                "report.json",
                "figures/data/gold_difference.csv",
                "frozen_corpus/source_snapshot.json",
            )
        },
    }


def count_summary(rows):
    supplied = [r for r in rows if r["status"] == "supplied"]
    error = sum(float(r["absolute_difference"]) for r in rows)
    overcount = sum(max(0.0, float(r["difference"])) for r in rows)
    undercount = sum(max(0.0, -float(r["difference"])) for r in rows)
    assert overcount + undercount == error
    return {
        "asserted_rows": len(rows),
        "supplied_rows": len(supplied),
        "supplied_fraction": rate(len(supplied), len(rows)),
        "supplied_exact_rows": sum(int(r["exact"]) for r in supplied),
        "supplied_exact_fraction": rate(
            sum(int(r["exact"]) for r in supplied), len(supplied)
        ),
        "conditional_mae": rate(
            sum(float(r["absolute_difference"]) for r in supplied), len(supplied)
        ),
        "absolute_error_sum": error,
        "overcount_units": overcount,
        "undercount_units_including_omissions": undercount,
        "end_to_end_mae": rate(error, len(rows)),
        "identity_miss_rows": sum(r["status"] == "identity_miss" for r in rows),
        "abstained_rows": sum(r["status"] == "abstained" for r in rows),
        "nonzero_gold_rows": sum(float(r["gold_count"]) != 0 for r in rows),
        "supplied_on_nonzero_gold": sum(float(r["gold_count"]) != 0 for r in supplied),
        "omitted_nonzero_gold_fields": sum(
            r["status"] != "supplied" and float(r["gold_count"]) > 0 for r in rows
        ),
        "supplied_zero_on_positive_gold_fields": sum(
            float(r["automated_count_evaluated"]) == 0 and float(r["gold_count"]) > 0
            for r in supplied
        ),
        "supplied_positive_on_gold_zero_fields": sum(
            float(r["automated_count_evaluated"]) > 0 and float(r["gold_count"]) == 0
            for r in supplied
        ),
        "supplied_wrong_fields": sum(not int(r["exact"]) for r in supplied),
    }


def summary(arm, allowed):
    papers = [p for k, p in arm["papers"].items() if k in allowed]
    tp, fp, fn = (sum(p[f] for p in papers) for f in ("tp", "fp", "fn"))
    counts = {
        field: count_summary(
            [r for r in arm["rows"] if key(r) in allowed and r["measure"] == field]
        )
        for field in FIELDS
    }
    counts["affected_unaffected_combined"] = count_summary(
        [
            r
            for r in arm["rows"]
            if key(r) in allowed and r["measure"] in ("affected", "unaffected")
        ]
    )
    for field in FIELDS:
        counts[field]["extra_variant_rows_with_supplied_count"] = sum(
            arm["extra_counts"][k][field] for k in allowed
        )
    counted = {
        name: sum(p["counted_precision"][name] for p in papers)
        for name in ("count_bearing_matched_rows", "counted_extra_rows", "matched_rows")
    }
    counted["precision_among_count_bearing_predictions"] = rate(
        counted["count_bearing_matched_rows"],
        counted["count_bearing_matched_rows"] + counted["counted_extra_rows"],
    )
    return {
        "attempts": len(papers),
        "pmids": len({p["pmid"] for p in papers}),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "recall": rate(tp, tp + fn),
        "precision": rate(tp, tp + fp),
        "counts": counts,
        "counted_precision": counted,
    }


def diagnostic_intervals(base, candidate, allowed):
    ids = sorted({pmid for _, pmid in allowed})
    positions = {pmid: i for i, pmid in enumerate(ids)}
    # tp, fp, fn, then asserted rows and absolute error for each count field.
    terms = np.zeros((len(ids), 2, 9))
    for j, arm in enumerate((base, candidate)):
        for k, paper in arm["papers"].items():
            if k in allowed:
                terms[positions[k[1]], j, :3] += [paper[f] for f in ("tp", "fp", "fn")]
        for row in arm["rows"]:
            if key(row) in allowed:
                i = positions[row["pmid"]]
                offset = 3 + 2 * FIELDS.index(row["measure"])
                terms[i, j, offset] += 1
                terms[i, j, offset + 1] += float(row["absolute_difference"])
    rng = np.random.default_rng(2026090502)
    weights = rng.multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=10000)
    totals = np.einsum("bi,ijk->bjk", weights, terms)
    distributions = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        for name, numerator, denominator in (
            ("recall", totals[:, :, 0], totals[:, :, 0] + totals[:, :, 2]),
            ("precision", totals[:, :, 0], totals[:, :, 0] + totals[:, :, 1]),
            *[
                (
                    f"{field}_end_to_end_mae",
                    totals[:, :, 4 + 2 * i],
                    totals[:, :, 3 + 2 * i],
                )
                for i, field in enumerate(FIELDS)
            ],
            (
                "affected_unaffected_combined_end_to_end_mae",
                totals[:, :, 6] + totals[:, :, 8],
                totals[:, :, 5] + totals[:, :, 7],
            ),
        ):
            rates = numerator / denominator
            delta = rates[:, 1] - rates[:, 0]
            delta = delta[np.isfinite(delta)]
            distributions[name] = {
                "two_sided_95_percentile_interval": np.quantile(
                    delta, [0.025, 0.975], method="inverted_cdf"
                ).tolist(),
                "valid_resamples": len(delta),
            }
    return {
        "classification": "descriptive, no acceptance decision",
        "cluster_unit": "PMID",
        "clusters": len(ids),
        "seed": 2026090502,
        "resamples": 10000,
        "deltas": distributions,
    }


def compare(base, candidate, allowed):
    a, b = summary(base, allowed), summary(candidate, allowed)
    assert a["tp"] + a["fn"] == b["tp"] + b["fn"]
    changes = {
        field: b[field] - a[field]
        for field in ("tp", "fp", "fn", "recall", "precision")
    }
    changes["counts"] = {}
    for field in a["counts"]:
        old, new = a["counts"][field], b["counts"][field]
        assert old["asserted_rows"] == new["asserted_rows"]
        changes["counts"][field] = {
            "supplied_rows": new["supplied_rows"] - old["supplied_rows"],
            "end_to_end_mae": new["end_to_end_mae"] - old["end_to_end_mae"],
            "relative_error_reduction": rate(
                old["absolute_error_sum"] - new["absolute_error_sum"],
                old["absolute_error_sum"],
            ),
        }
    return {
        "baseline": a,
        "candidate": b,
        "delta_candidate_minus_baseline": changes,
        "diagnostic_bootstrap": diagnostic_intervals(base, candidate, allowed),
        "count_error_transitions": count_error_transitions(base, candidate, allowed),
    }


def count_error_transitions(base, candidate, allowed):
    def row_map(arm):
        return {
            (r["gene"], r["pmid"], r["gold_row_index"], r["measure"]): r
            for r in arm["rows"]
            if key(r) in allowed
        }

    old_rows, new_rows = row_map(base), row_map(candidate)
    assert old_rows.keys() == new_rows.keys()
    groups = {}
    for k, old in old_rows.items():
        new = new_rows[k]
        assert old["gold_count"] == new["gold_count"]
        group_key = (old["measure"], old["status"], new["status"])
        group = groups.setdefault(
            group_key,
            {
                "field": group_key[0],
                "baseline_status": group_key[1],
                "candidate_status": group_key[2],
                "rows": 0,
                "nonzero_gold_rows": 0,
                "baseline_absolute_error": 0.0,
                "candidate_absolute_error": 0.0,
            },
        )
        group["rows"] += 1
        group["nonzero_gold_rows"] += int(float(old["gold_count"]) != 0)
        group["baseline_absolute_error"] += float(old["absolute_difference"])
        group["candidate_absolute_error"] += float(new["absolute_difference"])
    for group in groups.values():
        group["error_reduction"] = (
            group["baseline_absolute_error"] - group["candidate_absolute_error"]
        )
    return [groups[k] for k in sorted(groups)]


def paper_count_changes(base, candidate):
    output = []
    for k in sorted(base["papers"]):
        row = {"gene": k[0], "pmid": k[1]}
        for name, arm in (("baseline", base), ("candidate", candidate)):
            selected = [r for r in arm["rows"] if key(r) == k]
            row[name] = {
                field: sum(
                    float(r["absolute_difference"])
                    for r in selected
                    if r["measure"] == field
                )
                for field in FIELDS
            }
            for field in ("affected", "unaffected"):
                row[name][f"supplied_{field}"] = sum(
                    r["status"] == "supplied" for r in selected if r["measure"] == field
                )
        row["delta"] = {
            field: row["candidate"][field] - value
            for field, value in row["baseline"].items()
        }
        output.append(row)
    return output


def failed_paper_sensitivity(base, candidate, allowed):
    """Descriptive only: treat the archived operational failure as no output."""
    retry = read(OUT / "operational_retry_03_baseline.json")
    assert retry["retry_exit_code"] == 0
    affected = {k for k in allowed if k[0] == retry["gene"]}
    assert len(affected) == 1
    counterfactual = copy.deepcopy(base)
    for k in affected:
        paper = counterfactual["papers"][k]
        paper["fn"] += paper["tp"]
        paper["tp"] = paper["fp"] = 0
        for field in (
            "count_bearing_matched_rows",
            "counted_extra_rows",
            "matched_rows",
        ):
            paper["counted_precision"][field] = 0
        counterfactual["extra_counts"][k] = dict.fromkeys(FIELDS, 0)
    for row in counterfactual["rows"]:
        if key(row) in affected:
            gold = float(row["gold_count"])
            row.update(
                status="identity_miss",
                automated_count_raw="",
                automated_count_evaluated="0",
                difference=str(-gold),
                absolute_difference=str(abs(gold)),
                exact=str(int(gold == 0)),
            )
    result = compare(counterfactual, candidate, allowed)
    expected_lost_tp = sum(base["papers"][k]["tp"] for k in affected)
    assert result["baseline"]["tp"] == summary(base, allowed)["tp"] - expected_lost_tp
    return {
        "classification": "descriptive operational sensitivity; not a registered verdict",
        "scenario": "Baseline BRCA1 operational failure treated as empty output instead of its first successful unchanged-runtime retry. Candidate remains locked output.",
        "affected_attempts": [dict(gene=k[0], pmid=k[1]) for k in sorted(affected)],
        "comparison": result,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tranches", nargs="+", choices=["02", "03"])
    args = parser.parse_args()
    exposure = read(OUT / "prior_exposure.json")
    result = {"tranches": {}, "artifacts": {}}
    pooled = [{"papers": {}, "rows": [], "extra_counts": {}} for _ in range(2)]
    pooled_novel = set()
    for number in args.tranches:
        base, candidate = (load_arm(number, arm) for arm in ("baseline", "candidate"))
        allowed = set(base["papers"])
        assert allowed == set(candidate["papers"])
        assert base["snapshot"] == candidate["snapshot"], "initial source files differ"
        novel = {
            key(row)
            for row in exposure[f"tranche_{number}"]["previously_unscored_attempts"]
        }
        assert novel <= allowed
        item = {
            "registered_comparison": read(
                RUNS / candidate["run_id"] / "compare_with_secondary.json"
            ),
            "all_registered_attempts": compare(base, candidate, allowed),
            "previously_unscored_attempts": compare(base, candidate, novel),
            "by_gene": {
                gene: {
                    "baseline": summary(base, {k for k in allowed if k[0] == gene}),
                    "candidate": summary(
                        candidate, {k for k in allowed if k[0] == gene}
                    ),
                }
                for gene in sorted({k[0] for k in allowed})
            },
            "sources": {
                "initial_snapshot_equal": True,
                "actual_rendering_equal_attempts": sum(
                    base["selection"][k]["source_sha256"]
                    == candidate["selection"][k]["source_sha256"]
                    for k in allowed
                ),
                "attempts": len(allowed),
            },
            "paper_deltas": [
                {
                    "gene": k[0],
                    "pmid": k[1],
                    "previously_unscored": k in novel,
                    "baseline": {
                        f: base["papers"][k][f] for f in ("tp", "fp", "fn", "count")
                    },
                    "candidate": {
                        f: candidate["papers"][k][f]
                        for f in ("tp", "fp", "fn", "count")
                    },
                }
                for k in sorted(allowed)
            ],
        }
        result["tranches"][number] = item
        if number == "03":
            item["operational_failure_sensitivity"] = failed_paper_sensitivity(
                base, candidate, allowed
            )
        (OUT / f"paper_count_changes_{number}.json").write_text(
            json.dumps(paper_count_changes(base, candidate), indent=2) + "\n"
        )
        for index, arm in enumerate((base, candidate)):
            result["artifacts"][arm["run_id"]] = arm["sha256"]
            assert not set(pooled[index]["papers"]) & allowed, (
                "pooled duplicate attempts"
            )
            pooled[index]["papers"].update(arm["papers"])
            pooled[index]["rows"].extend(arm["rows"])
            pooled[index]["extra_counts"].update(arm["extra_counts"])
        pooled_novel |= novel
    if len(args.tranches) > 1:
        result["pooled"] = {
            "all_registered_attempts": compare(*pooled, set(pooled[0]["papers"])),
            "previously_unscored_attempts": compare(*pooled, pooled_novel),
        }
    prior = read(OUT.parent / "phenotype_failure_panel_20260905/manifest.json")
    for path, expected in prior["input_sha256"].items():
        assert sha(ROOT / path) == expected, path
    result["prior_locked_inputs_verified"] = len(prior["input_sha256"])
    (OUT / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                n: r["all_registered_attempts"]["delta_candidate_minus_baseline"]
                for n, r in result["tranches"].items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
