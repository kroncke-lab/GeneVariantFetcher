from pathlib import Path

import pytest

from benchmarks.tier2_relevance_eval.run_shadow import (
    _aggregate_telemetry,
    _load_pinned_selections,
    _prepare_outdir,
    _quality_summary,
    _stable_sample,
)


def test_prepare_outdir_refuses_to_mix_runs(tmp_path: Path):
    outdir = tmp_path / "shadow"
    _prepare_outdir(outdir)
    (outdir / "trace.json").write_text("{}", encoding="utf-8")

    with pytest.raises(FileExistsError, match="choose a fresh --outdir"):
        _prepare_outdir(outdir)


def test_load_pinned_selections_rejects_duplicate_pmids(tmp_path: Path):
    manifest = tmp_path / "cohort.tsv"
    manifest.write_text(
        "group\tpmid\nproductive_positive\t1\nhigh_confidence_negative\t1\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate PMID"):
        _load_pinned_selections(manifest)


def test_stable_sample_is_order_independent():
    rows = [{"pmid": str(value)} for value in range(10)]

    forward = _stable_sample(rows, 4, seed="fixed")
    reverse = _stable_sample(reversed(rows), 4, seed="fixed")

    assert forward == reverse
    assert len({row["pmid"] for row in forward}) == 4


def test_quality_summary_separates_final_and_raw_agreement():
    rows = [
        {
            "group": "productive_positive",
            "historical_final_decision": "pass",
            "historical_raw_model_decision": "pass",
            "luna_final_decision": "pass",
            "luna_raw_model_decision": "pass",
            "luna_fail_open_target_gene_signal": False,
            "luna_call_error": None,
        },
        {
            "group": "fail_open_boundary",
            "historical_final_decision": "pass",
            "historical_raw_model_decision": "fail",
            "luna_final_decision": "pass",
            "luna_raw_model_decision": "pass",
            "luna_fail_open_target_gene_signal": False,
            "luna_call_error": None,
        },
    ]

    summary = _quality_summary(rows)

    assert summary["final_decision_agreement"] == 1.0
    assert summary["raw_model_decision_agreement"] == 0.5
    assert summary["groups"]["productive_positive"]["luna_expected_rate"] == 1.0


def test_aggregate_telemetry_sums_usage_and_reports_models():
    summary = _aggregate_telemetry(
        {
            "1": {
                "calls": 1,
                "input_tokens": 100,
                "output_tokens": 20,
                "reasoning_tokens": 5,
                "total_tokens": 120,
                "duration_seconds": 2.0,
                "requested_models": ["azure_ai/gpt-5.6-luna"],
                "response_models": ["gpt-5.6-luna-snapshot"],
                "reasoning_efforts": ["xhigh"],
                "failures": 0,
            },
            "2": {
                "calls": 2,
                "input_tokens": 200,
                "output_tokens": 40,
                "reasoning_tokens": 10,
                "total_tokens": 240,
                "duration_seconds": 4.0,
                "requested_models": ["azure_ai/gpt-5.6-luna"],
                "response_models": ["gpt-5.6-luna-snapshot"],
                "reasoning_efforts": ["xhigh"],
                "failures": 1,
            },
        }
    )

    assert summary["calls"] == 3
    assert summary["total_tokens"] == 360
    assert summary["reasoning_tokens"] == 15
    assert summary["provider_failures"] == 1
    assert summary["requested_models"] == ["azure_ai/gpt-5.6-luna"]
    assert summary["reasoning_efforts"] == ["xhigh"]
