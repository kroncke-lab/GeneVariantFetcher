"""Verify locked local evidence and emit a compact receipt; no API calls."""

import json
import math
import sys
from pathlib import Path

from summarize import OUT, ROOT, RUNS, read, sha

sys.path.insert(0, str(ROOT))
from benchmarks.codex_paper_eval.setup_production_eval import runtime_fingerprint


def main():
    runtime = read(OUT / "candidate_runtime.json")
    assert runtime_fingerprint() == runtime["runtime"]
    for name, expected in runtime["files"].items():
        assert sha(ROOT / name) == expected, name
    old = read(OUT.parent / "phenotype_failure_panel_20260905/manifest.json")
    for name, expected in old["input_sha256"].items():
        assert sha(ROOT / name) == expected, name
    results = read(OUT / "results.json")
    assert results["prior_locked_inputs_verified"] == 34
    runs = {}
    for number in ("02", "03"):
        pair = results["tranches"][number]
        assert not pair["registered_comparison"]["passed"]
        assert not pair["registered_comparison"]["secondary_count_endpoint"]["passed"]
        assert pair["sources"]["actual_rendering_equal_attempts"] == (
            120 if number == "02" else 119
        )
        for arm in ("baseline", "candidate"):
            run = RUNS / f"20260905_protocol_cont120_{number}_{arm}"
            lock = read(run / "LOCK.json")
            for name in ("selection", "predictions", "setup"):
                assert sha(run / f"{name}.json") == lock[f"{name}_sha256"]
            for name, expected in results["artifacts"][run.name].items():
                assert sha(run / name) == expected
            selection = read(run / "selection.json")["papers"]
            assert len(selection) == 120
            for paper in selection:
                assert sha(Path(paper["source"])) == paper["source_sha256"]
            states = []
            for status in lock["production_run_statuses"]:
                path = ROOT / status["status"]
                assert sha(path) == status["sha256"]
                data = read(path)
                assert data["exit_code"] == 0 and data["status"] == "completed"
                assert data["gold_access"]["disabled"]
                assert data["gold_access"]["gold_derived_alias_files_disabled"]
                states.append(
                    {
                        "gene": status["gene"],
                        "severity": data["severity"],
                        "stage_warnings": data.get("stage_warnings", []),
                        "stage_failures": data.get("stage_failures", []),
                        "trace_omissions": data.get("llm_trace", {}).get("omissions"),
                        "missing_decision_links": data.get("llm_trace", {}).get(
                            "missing_decision_links", []
                        ),
                        "source_integrity": data.get("source_integrity", {}),
                    }
                )
            assert len(states) == 7
            runs[run.name] = {
                "selected_attempts": len(selection),
                "completed_gold_disabled_gene_jobs": len(states),
                "status_details": states,
            }
    supplement = read(OUT / "supplement_check_results.json")
    run = RUNS / supplement["run_id"]
    for name, expected in supplement["artifact_sha256"].items():
        assert sha(run / name) == expected
    assert [supplement["overall"][f] for f in ("tp", "fp", "fn")] == [20, 0, 0]
    assert all(
        supplement["overall"]["count"][f]["predicted"] == 0
        for f in ("carriers", "affected", "unaffected")
    )
    corpus = ROOT / "corpus"
    assert corpus.is_symlink() and corpus.is_dir()
    for name, expected in read(OUT / "supplement_corpus_merge.json")[
        "verified_source_files"
    ].items():
        assert sha(ROOT / name) == expected, name
    for item in read(OUT / "figure_qa.json")["checked"]:
        assert sha(ROOT / item["path"]) == item["sha256"]
    for name, expected in read(OUT / "updated_forecast.json")["input_sha256"].items():
        assert sha(OUT / name) == expected
    budget = read(OUT / "budget.json")
    extra = sum(
        row["api_proxy_usd"]
        for group in (
            "completed_arms",
            "operational_failed_attempts",
            "additional_opened_checks",
        )
        for row in budget[group]
    )
    assert math.isclose(extra, 40.04289465)
    assert math.isclose(
        budget["prior_api_test_proxy_usd"] + extra, budget["campaign_api_proxy_usd"]
    )
    assert math.isclose(
        budget["campaign_api_proxy_usd"] + budget["remaining_api_envelope_usd"], 100
    )
    sensitivity = results["tranches"]["03"]["operational_failure_sensitivity"][
        "comparison"
    ]
    main_comparison = results["tranches"]["03"]["all_registered_attempts"]
    assert sensitivity["baseline"] == main_comparison["baseline"]
    assert sensitivity["candidate"] == main_comparison["candidate"]
    registry = (
        ROOT / "benchmarks/evaluation_tiers/mixed_gold_continuation_120/registry.json"
    )
    assert (
        sha(registry)
        == "9a6bbbf6b386aeff7c6cb500e46c0b4a76cd2c5ef7a54f58dc980dc7959c1f54"
    )
    receipt = {
        "classification": "Local integrity and arithmetic verification; not proof of clinical correctness",
        "runtime": runtime["runtime"],
        "runtime_files_verified": len(runtime["files"]),
        "prior_locked_input_files_verified": len(old["input_sha256"]),
        "registry_sha256": sha(registry),
        "runs": runs,
        "initial_snapshot_equal_attempts": 240,
        "actual_primary_text_equal_attempts": 239,
        "successful_attempt_executions": 480,
        "additional_supplement_check_attempts": 1,
        "preserved_operational_failed_attempts": 1,
        "new_api_proxy_usd": extra,
        "campaign_api_proxy_usd": budget["campaign_api_proxy_usd"],
        "supplement_source_hashes_verified": 6,
        "figure_qa_hashes_verified": 12,
        "operational_empty_failure_sensitivity_equal": True,
        "caveats": [
            "Source classes are availability flags, not proof of count-bearing evidence.",
            "RUN_STATUS completion includes the recorded stage warnings and missing trace-decision links below.",
            "Frozen availability matches, but one primary input differs due to protocol validation.",
            "Gold agreement does not adjudicate H558R unit or control-study attribution.",
            "Both per-tranche identity and carrier gates fail; no promotion.",
        ],
        "analysis_files_sha256": {
            name: sha(OUT / name)
            for name in (
                "results.json",
                "outlier_sensitivity.json",
                "updated_forecast.json",
                "source_review_03_receipt.json",
                "budget.json",
            )
        },
    }
    (OUT / "verification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        {
            k: v
            for k, v in receipt.items()
            if k not in ("runs", "analysis_files_sha256", "caveats")
        }
    )


if __name__ == "__main__":
    main()
