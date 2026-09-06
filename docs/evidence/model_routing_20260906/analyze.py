"""Summarize already-locked and scored calibration arms; no answer-key reads."""

import csv
from collections import Counter
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from docs.evidence.tranche_validation_20260905.summarize import count_summary, sha, key
from benchmarks.codex_paper_eval.run_eval import merge_notation_twins

HERE = Path(__file__).resolve().parent
RUNS = ROOT / "benchmarks/codex_paper_eval/runs"
FIELDS = ("carriers", "affected", "unaffected")


def load(name):
    run = RUNS / name
    lock = json.loads((run / "LOCK.json").read_text())
    for what in ("selection", "predictions"):
        assert sha(run / (what + ".json")) == lock[what + "_sha256"]
    report = json.loads((run / "report.json").read_text())
    with (run / "figures/data/gold_difference.csv").open() as file:
        rows = [r for r in csv.DictReader(file) if r["analysis_run"] == name]
    assert len(
        {(r["gene"], r["pmid"], r["gold_row_index"], r["measure"]) for r in rows}
    ) == len(rows)
    predictions = json.loads((run / "predictions.json").read_text())
    report_by_key = {key(p): p for p in report["papers"]}
    extra_counts = {}
    for paper in predictions["papers"]:
        merged, _ = merge_notation_twins(paper.get("variants", []), paper["gene"])
        scored = report_by_key[key(paper)]
        extra_names = set(scored["extra_predictions"])
        extras = [row for row in merged if row["variant"] in extra_names]
        assert len(extras) == scored["fp"]
        extra_counts[key(paper)] = {
            f: sum(row.get(f) is not None for row in extras) for f in FIELDS
        }

    def metrics(papers, count_rows):
        tp, fp, fn = (sum(p[k] for p in papers) for k in ("tp", "fp", "fn"))
        return {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "recall": tp / (tp + fn) if tp + fn else None,
            "precision": tp / (tp + fp) if tp + fp else None,
            "non_null_count_fields_on_identity_extras": {
                f: sum(extra_counts[key(p)][f] for p in papers) for f in FIELDS
            },
            "counted_extra_rows": sum(
                p["counted_precision"]["counted_extra_rows"] for p in papers
            ),
            "counts": {
                f: count_summary([r for r in count_rows if r["measure"] == f])
                for f in FIELDS
            },
        }

    all_metrics = metrics(report["papers"], rows)
    by_paper = {
        ":".join(key(p)): metrics([p], [r for r in rows if key(r) == key(p)])
        for p in report["papers"]
    }
    return {
        "run_id": name,
        "overall": all_metrics,
        "without_20129283": metrics(
            [p for p in report["papers"] if str(p["pmid"]) != "20129283"],
            [r for r in rows if str(r["pmid"]) != "20129283"],
        ),
        "without_both_preflagged_disputes": metrics(
            [
                p
                for p in report["papers"]
                if str(p["pmid"]) not in {"20129283", "25814417"}
            ],
            [r for r in rows if str(r["pmid"]) not in {"20129283", "25814417"}],
        ),
        "without_deterministic_shortcuts_posthoc": metrics(
            [
                p
                for p in report["papers"]
                if str(p["pmid"]) not in {"20129283", "30059973"}
            ],
            [r for r in rows if str(r["pmid"]) not in {"20129283", "30059973"}],
        ),
        "by_paper": by_paper,
        "by_gold_provenance": report.get("by_gold_provenance"),
        "secondary_lanes": report.get("provenance_lane_scores"),
        "run_lock_sha256": sha(run / "LOCK.json"),
        "report_sha256": sha(run / "report.json"),
        "primary_model_calls": {},
        "incremental_clinical": predictions.get("clinical_reader_audits"),
        "clinical_usage_status": {
            ":".join(key(p)): p["token_usage"].get("status")
            for p in predictions["papers"]
        }
        if predictions.get("clinical_reader_audits") is not None
        else None,
    }


def main():
    names = [
        "20260906_model12_grok43",
        "20260906_model12_astra_medium_verified",
        "20260906_model12_grok43_astra_clinical_verified",
    ]
    result = {
        "classification": "Descriptive opened calibration, not a corpus/human ceiling estimate or promotion gate",
        "arms": {name: load(name) for name in names},
    }
    budget = json.loads((HERE / "budget.json").read_text())
    assert not [r for r in budget["calls"] if r["status"] == "reserved"], (
        "Calls still running"
    )
    models = {}
    for r in budget["calls"]:
        group = models.setdefault(
            r["arm"],
            {
                "calls": 0,
                "returned_proxy_usd": 0,
                "unknown_usage_reserve_usd": 0,
                "models": {},
            },
        )
        group["calls"] += 1
        if r["status"] == "returned":
            group["returned_proxy_usd"] += r["api_proxy_usd"]
            group["unknown_usage_reserve_usd"] += r.get("unknown_retry_reserve_usd", 0)
        else:
            group["unknown_usage_reserve_usd"] += r["reserved_usd"]
        entry = group["models"].setdefault(
            r["model"], {"calls": 0, "returned_proxy_usd": 0, "papers": []}
        )
        entry["calls"] += 1
        entry["returned_proxy_usd"] += r.get("api_proxy_usd", 0)
        if r["pmid"] and r["pmid"] not in entry["papers"]:
            entry["papers"].append(r["pmid"])
    result["cost_by_arm"] = models
    result["returned_api_proxy_usd"] = sum(
        g["returned_proxy_usd"] for g in models.values()
    )
    result["unknown_usage_reserve_usd"] = sum(
        g["unknown_usage_reserve_usd"] for g in models.values()
    )
    result["additional_uncertainty_margin_usd"] = budget[
        "smoke_uncertainty_reserve_usd"
    ]
    result["conservative_envelope_usd"] = (
        result["returned_api_proxy_usd"]
        + result["unknown_usage_reserve_usd"]
        + result["additional_uncertainty_margin_usd"]
    )
    assert result["conservative_envelope_usd"] <= budget["api_ceiling_usd"]
    result["clinical_trigger_summary"] = {
        name: {
            "triggered": sum(
                a["trigger"]["triggered"] for a in arm["incremental_clinical"]
            ),
            "not_triggered": sum(
                not a["trigger"]["triggered"] for a in arm["incremental_clinical"]
            ),
            "completed_without_error": sum(
                a["trigger"]["triggered"] and not a.get("error")
                for a in arm["incremental_clinical"]
            ),
            "budget_not_dispatched": sum(
                status == "budget_not_dispatched"
                for status in arm["clinical_usage_status"].values()
            ),
            "errors": [
                {"gene": a["gene"], "pmid": a["pmid"], "error": a["error"]}
                for a in arm["incremental_clinical"]
                if a.get("error")
            ],
            "accepted_field_fills": sum(
                len(a["accepted"]) for a in arm["incremental_clinical"]
            ),
            "raw_field_fills": sum(
                len(a["raw_fills"]) for a in arm["incremental_clinical"]
            ),
            "contradictions": sum(
                len(a["contradictions"]) for a in arm["incremental_clinical"]
            ),
            "rejection_reasons": dict(
                Counter(
                    item.get("reason")
                    for a in arm["incremental_clinical"]
                    for item in a["rejected"]
                )
            ),
        }
        for name, arm in result["arms"].items()
        if arm["incremental_clinical"] is not None
    }
    for name, arm in result["arms"].items():
        short = name.removeprefix("20260906_model12_")
        primary = {
            "grok43": "azure_ai/grok-4.3",
            "grok46_verified": "azure_ai/grok-4.6",
            "astra_medium_verified": "azure_ai/gpt-6-astra",
        }.get(short)
        arm["primary_model_calls"] = (
            models.get(short, {}).get("models", {}).get(primary) if primary else None
        )
    (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    for name, arm in result["arms"].items():
        m = arm["overall"]
        print(
            name,
            m["tp"],
            m["fp"],
            m["fn"],
            "recall",
            round(m["recall"] * 100, 2),
            "precision",
            round(m["precision"] * 100, 2),
        )
        for field, c in m["counts"].items():
            print(
                field,
                "exact",
                c["supplied_exact_rows"],
                "of supplied",
                c["supplied_rows"],
                "asserted",
                c["asserted_rows"],
                "positive reference covered",
                c["supplied_on_nonzero_gold"],
                "/",
                c["nonzero_gold_rows"],
                "AE",
                c["absolute_error_sum"],
                "wrong",
                c["supplied_wrong_fields"],
            )
    print("cost", result["returned_api_proxy_usd"])


if __name__ == "__main__":
    main()
