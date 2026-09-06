"""Native locks retain traced API failures without fabricating zero usage."""

import hashlib
import json

import pytest

from benchmarks.codex_paper_eval.run_eval import validate_predictions


def fixture(tmp_path):
    record = {
        "trace_id": "failed-request",
        "record_type": "llm_call",
        "context": {"gene": "SCN5A", "pmid": "123"},
        "response": {"success": False, "error": {"type": "Timeout"}, "usage": None},
    }
    path = tmp_path / "failed.json"
    path.write_text(json.dumps(record))
    ref = {
        "trace_id": "failed-request",
        "record_type": "llm_call",
        "path": path.name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "context": {"gene": "SCN5A", "pmid": "123", "stage": "paper_curation"},
    }
    record["context"]["stage"] = "paper_curation"
    path.write_text(json.dumps(record))
    ref["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    (tmp_path / "trace_index.jsonl").write_text(json.dumps(ref) + "\n")
    paper = {
        "gene": "SCN5A",
        "pmid": "123",
        "tool": "text",
        "tool_rationale": "source-only attempt",
        "elapsed_seconds": 1200,
        "source_completeness": "full_text",
        "curation_rationale": "Traced timeout; retain previous result.",
        "variants": [],
        "token_usage": {
            "telemetry_available": False,
            "input_tokens": None,
            "output_tokens": None,
            "total_tokens": None,
            "status": "unknown_failed_call",
            "unknown_failed_call_trace_ids": ["failed-request"],
            "known_usage": {"input_tokens": 0, "output_tokens": 0, "total_tokens": 0},
        },
        "llm_trace_refs": [ref]
        + [
            {"context": {"stage": s}}
            for s in (
                "representation_route",
                "representation_route_decision",
                "paper_curation_decision",
            )
        ],
    }
    return (
        {"papers": [{"gene": "SCN5A", "pmid": "123"}]},
        {
            "schema_version": 2,
            "papers": [paper],
            "token_usage": {
                "telemetry_available": False,
                "input_tokens": None,
                "output_tokens": None,
                "total_tokens": None,
                "known_usage": {
                    "input_tokens": 0,
                    "output_tokens": 0,
                    "total_tokens": 0,
                },
            },
        },
        record,
        path,
    )


def test_native_traced_timeout_accepts_explicit_null_usage(tmp_path):
    selection, predictions, _, _ = fixture(tmp_path)
    assert validate_predictions(selection, predictions, trace_root=tmp_path) == []
    assert validate_predictions(selection, predictions) == [
        "SCN5A:123: missing exact token telemetry"
    ]


@pytest.mark.parametrize(
    "case",
    [
        "zero",
        "missing_id",
        "bad_hash",
        "wrong_paper",
        "success",
        "known_usage",
        "non_object",
    ],
)
def test_native_unknown_usage_requires_authentic_failed_trace(tmp_path, case):
    selection, predictions, record, path = fixture(tmp_path)
    p = predictions["papers"][0]
    if case == "zero":
        p["token_usage"]["total_tokens"] = 0
    elif case == "missing_id":
        p["token_usage"]["unknown_failed_call_trace_ids"] = []
    elif case == "bad_hash":
        p["llm_trace_refs"][0]["sha256"] = "bad"
    else:
        if case == "wrong_paper":
            record["context"]["pmid"] = "another-paper"
        elif case == "success":
            record["response"]["success"] = True
        elif case == "known_usage":
            record["response"]["usage"] = {"total_tokens": 20}
        elif case == "non_object":
            record = []
        path.write_text(json.dumps(record))
        p["llm_trace_refs"][0]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    assert "SCN5A:123: missing exact token telemetry" in validate_predictions(
        selection, predictions, trace_root=tmp_path
    )


@pytest.mark.parametrize(
    "case", ["missing_index", "missing_ref", "duplicate_id", "path_escape"]
)
def test_unknown_failure_cannot_bypass_write_time_binding(tmp_path, case):
    selection, predictions, record, path = fixture(tmp_path)
    paper = predictions["papers"][0]
    if case == "missing_index":
        (tmp_path / "trace_index.jsonl").unlink()
    elif case == "missing_ref":
        paper["llm_trace_refs"].pop(0)
    elif case == "duplicate_id":
        paper["token_usage"]["unknown_failed_call_trace_ids"].append("failed-request")
    elif case == "path_escape":
        paper["llm_trace_refs"][0]["path"] = "../failed.json"
    assert "SCN5A:123: missing exact token telemetry" in validate_predictions(
        selection, predictions, trace_root=tmp_path
    )


def test_known_sibling_usage_is_preserved_alongside_unknown_call(tmp_path):
    selection, predictions, record, path = fixture(tmp_path)
    record["trace_id"] = "successful-sibling"
    record["response"] = {
        "success": True,
        "usage": {"prompt_tokens": 11, "completion_tokens": 4, "total_tokens": 15},
    }
    path = tmp_path / "success.json"
    path.write_text(json.dumps(record))
    entry = {
        "trace_id": record["trace_id"],
        "record_type": "llm_call",
        "context": record["context"],
        "path": path.name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    with (tmp_path / "trace_index.jsonl").open("a") as handle:
        handle.write(json.dumps(entry) + "\n")
    assert "SCN5A:123: missing exact token telemetry" in validate_predictions(
        selection, predictions, trace_root=tmp_path
    )
    predictions["papers"][0]["token_usage"]["known_usage"] = {
        "input_tokens": 11,
        "output_tokens": 4,
        "total_tokens": 15,
    }
    assert (
        "run: failed-call usage requires null totals and the complete known subset"
        in validate_predictions(selection, predictions, trace_root=tmp_path)
    )
    predictions["token_usage"]["known_usage"] = {
        "input_tokens": 11,
        "output_tokens": 4,
        "total_tokens": 15,
    }
    assert validate_predictions(selection, predictions, trace_root=tmp_path) == []


def test_unknown_run_usage_cannot_claim_exact_total(tmp_path):
    selection, predictions, _, _ = fixture(tmp_path)
    predictions["token_usage"].update(telemetry_available=True, total_tokens=0)
    assert (
        "run: failed-call usage requires null totals and the complete known subset"
        in validate_predictions(selection, predictions, trace_root=tmp_path)
    )
