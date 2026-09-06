"""Decision-score regressions: failures, acceptance independence and additive fills."""

import importlib.util
import json
from pathlib import Path

import pytest


@pytest.fixture
def scorer(monkeypatch, tmp_path):
    path = (
        Path(__file__).resolve().parents[2]
        / "docs/evidence/astra_value_20260906/score.py"
    )
    spec = importlib.util.spec_from_file_location("astra_value_score", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "HERE", tmp_path)
    return module


@pytest.fixture
def example():
    query = {"id": "q1", "variant": "p.R123H", "field": "carriers"}
    text = "The R123H variant was observed in 3 carriers."
    packet = {
        "name": "synthetic",
        "gene": "RYR2",
        "pmid": "0",
        "source_units": {"L1": text},
        "queries": [query],
    }
    properties = {
        "id": {"type": "string", "enum": ["q1"]},
        "value": {"type": ["integer", "null"]},
        "basis": {"type": "string", "enum": ["explicit", "derived", "unknown"]},
        "sources": {
            "type": "array",
            "items": {"type": "string", "enum": ["L1"]},
        },
        "quote": {"type": ["string", "null"]},
        "reason": {"type": "string"},
    }
    schema = {
        "type": "object",
        "required": ["claims", "limitations"],
        "additionalProperties": False,
        "properties": {
            "claims": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": list(properties),
                    "additionalProperties": False,
                    "properties": properties,
                },
            },
            "limitations": {"type": "array", "items": {"type": "string"}},
        },
    }
    plan = {
        "name": "synthetic__test",
        "arm": "grok_first",
        "body": {"response_format": {"json_schema": {"schema": schema}}},
    }
    claim = {
        "id": "q1",
        "value": 3,
        "basis": "explicit",
        "sources": ["L1"],
        "quote": text,
        "reason": "Literal count.",
    }
    return packet, plan, claim


def grade(scorer, example, claims, expected=3, finish="stop", dispatched=True):
    packet, plan, _ = example
    response = {
        "status": "returned",
        "http_status": 200,
        "response": {
            "choices": [
                {
                    "finish_reason": finish,
                    "message": {
                        "content": json.dumps({"claims": claims, "limitations": []})
                    },
                }
            ]
        },
        "seconds": 1,
        "api_proxy_usd": 0.1,
        "accounted_usd": 0.12,
    }
    return scorer.grade_cell(
        plan, packet, {"q1": expected}, response if dispatched else None
    )["queries"][0]


def abstention(claim):
    return {**claim, "value": None, "basis": "unknown", "quote": None, "sources": []}


def test_reference_correctness_never_controls_literal_acceptance(scorer, example):
    claim = example[2]
    correct = grade(scorer, example, [claim])
    changed_reference = grade(scorer, example, [claim], expected=9)
    assert correct["accepted_correct"]
    assert changed_reference["accepted_wrong"]
    assert correct["literal_validator"] == changed_reference["literal_validator"]


@pytest.mark.parametrize("failure", ["length", "missing", "duplicate", "undispatched"])
def test_failed_negative_query_does_not_earn_abstention_credit(
    scorer, example, failure
):
    claim = abstention(example[2])
    rows = (
        []
        if failure == "missing"
        else [claim, claim]
        if failure == "duplicate"
        else [claim]
    )
    result = grade(
        scorer,
        example,
        rows,
        expected=None,
        finish="length" if failure == "length" else "stop",
        dispatched=failure != "undispatched",
    )
    metrics = scorer.metrics([result])
    assert metrics["planned_queries"] == metrics["reference_unknown"] == 1
    assert metrics["outcomes"]["response_failure_or_missing_claim"] == 1
    assert metrics["outcomes"]["correct_abstention"] == 0
    assert not result["numeric_exact"]


@pytest.mark.parametrize("invalid_value", [True, -1, 3.5])
def test_invalid_count_cannot_be_accepted(scorer, example, invalid_value):
    result = grade(scorer, example, [{**example[2], "value": invalid_value}])
    assert result["numeric_outcome"] == "response_failure_or_missing_claim"
    assert not result["current_validator_lane_accepted"]


def test_numeric_agreement_does_not_turn_derived_or_invented_quote_into_acceptance(
    scorer, example
):
    derived = grade(scorer, example, [{**example[2], "basis": "derived"}])
    invented_quote = grade(
        scorer,
        example,
        [
            {
                **example[2],
                "quote": "The R123H variant was observed in 3 affected people.",
            }
        ],
    )
    assert derived["numeric_exact"] and invented_quote["numeric_exact"]
    assert not derived["literal_validator"]["attempted"]
    assert not derived["current_validator_lane_accepted"]
    assert "quote_not_verbatim_in_packet" in invented_quote["mechanical_errors"]
    assert not invented_quote["current_validator_lane_accepted"]


def test_valid_abstention_is_not_a_missing_response(scorer, example):
    claim = abstention(example[2])
    correct = grade(scorer, example, [claim], expected=None)
    missed = grade(scorer, example, [claim], expected=3)
    assert correct["numeric_outcome"] == "correct_abstention"
    assert missed["numeric_outcome"] == "null_miss"
    assert correct["numeric_exact"] and not missed["numeric_exact"]


def test_missing_only_overlay_preserves_wrong_nonnull_and_surfaces_conflict(
    scorer, example
):
    first = grade(scorer, example, [{**example[2], "value": 2}])
    added = grade(scorer, example, [example[2]])
    out = scorer.compare_overlay([first], [added], "raw_candidates")
    assert out["queries"][0]["selected_value"] == 2
    assert out["conflict_preserving_first"] == 1
    assert out["new_correct_nonnull"] == 0
    assert out["final_wrong_nonnull"] == 1


def test_validator_overlay_accepts_grounded_null_fill(scorer, example):
    first = grade(scorer, example, [abstention(example[2])])
    added = grade(scorer, example, [example[2]])
    out = scorer.compare_overlay([first], [added], "current_validator")
    assert out["new_correct_nonnull"] == 1
    assert out["new_wrong_nonnull"] == 0
    assert out["queries"][0]["selected_from"] == "overlay"


def test_actual_responses_cannot_be_read_before_lock(scorer):
    with pytest.raises(RuntimeError, match="locked before scoring"):
        scorer.verify_inputs()


@pytest.fixture
def remediation_scorer(monkeypatch, tmp_path):
    path = (
        Path(__file__).resolve().parents[2]
        / "docs/evidence/astra_value_20260906/score_remediation.py"
    )
    spec = importlib.util.spec_from_file_location("astra_remediation_score", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "HERE", tmp_path)
    return module


@pytest.mark.parametrize("altered_file", ["amendment", "boundary_receipt"])
def test_remediation_rejects_altered_concurrency_evidence(
    remediation_scorer, altered_file
):
    scorer = remediation_scorer
    receipt = scorer.HERE / "remediation_boundary_stop.json"
    receipt.write_text(json.dumps({"pending_calls": 0}))
    amendment = scorer.HERE / "remediation_concurrency_amendment.json"
    amendment.write_text(
        json.dumps(
            {
                "plan_sha256": "fixture_plan",
                "recorded_unix": 2,
                "stop_receipt_sha256": scorer.primary.sha(receipt),
            }
        )
    )
    lock = {
        "plan_sha256": "fixture_plan",
        "locked_unix": 3,
        "concurrency_amendment_sha256": scorer.primary.sha(amendment),
    }
    plan = {"prepared_unix": 1}
    assert scorer.verify_concurrency_amendment(plan, lock)
    target = amendment if altered_file == "amendment" else receipt
    target.write_text(target.read_text() + "\n")
    with pytest.raises(AssertionError):
        scorer.verify_concurrency_amendment(plan, lock)


@pytest.mark.parametrize("unknown_arm", ["grok", "astra"])
def test_remediation_cost_ratio_requires_matching_known_usage(
    scorer, remediation_scorer, example, unknown_arm
):
    row = grade(scorer, example, [example[2]])
    common = {
        "packet": "synthetic",
        "queries": [row],
        "dispatched": True,
        "request_successful_stop": True,
        "seconds": 5,
    }
    grok = {**common, "api_proxy_usd": 0.02, "accounted_usd": 0.025}
    astra = {**common, "api_proxy_usd": 0.1, "accounted_usd": 0.12}
    out = remediation_scorer.subset_comparison(
        "same_paper", ["synthetic"], [grok], [astra]
    )
    assert out["grok_to_astra_returned_proxy_ratio"] == pytest.approx(0.2)
    (grok if unknown_arm == "grok" else astra)["api_proxy_usd"] = None
    out = remediation_scorer.subset_comparison(
        "same_paper_unknown_usage", ["synthetic"], [grok], [astra]
    )
    assert not out["all_planned_paired_requests_dispatched_with_known_usage"]
    assert out["grok_to_astra_returned_proxy_ratio"] is None


def test_42_query_sensitivity_removes_only_disputed_affected_inference(
    scorer, remediation_scorer, example
):
    exact = grade(scorer, example, [example[2]])
    wrong = grade(scorer, example, [example[2]], expected=9)
    rows = [{**exact, "id": f"q{i}"} for i in range(40)]
    rows += [
        {**wrong, "packet": "ryr193_abstract", "id": f"q1_{field}", "field": field}
        for field in ("carriers", "affected", "unaffected")
    ]
    full, sensitivity = remediation_scorer.scoped_metrics({"synthetic_arm": rows})
    assert full["synthetic_arm"]["planned_queries"] == 43
    assert full["synthetic_arm"]["outcomes"]["wrong_nonnull"] == 3
    assert sensitivity["synthetic_arm"]["planned_queries"] == 42
    assert sensitivity["synthetic_arm"]["outcomes"]["wrong_nonnull"] == 2
    assert len(rows) == 43  # No mutation or replacement of primary observations.
