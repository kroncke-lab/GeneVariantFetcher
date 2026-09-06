"""A real pre-dispatch refusal is zero usage; an uncertain API call is not."""

import hashlib
import json

import pytest

from docs.evidence.model_routing_20260906.bind_failed_usage import (
    documented_budget_non_dispatch,
)


@pytest.fixture
def case(tmp_path):
    root = tmp_path / "llm_traces"
    root.mkdir()
    error = "RuntimeError: Campaign API budget prevents call to azure_ai/gpt-6-astra; reservation $7.123"
    record = {
        "record_type": "decision_event",
        "trace_id": "refusal-1",
        "context": {"gene": "SCN5A", "pmid": "32533946", "stage": "paper_curation"},
        "event": {
            "type": "paper_curation",
            "data": {"status": "failed", "error": error},
        },
    }
    path = root / "refusal.json"
    path.write_text(json.dumps(record))
    entry = {
        "record_type": "decision_event",
        "trace_id": record["trace_id"],
        "context": dict(record["context"]),
        "path": path.name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    index = root / "trace_index.jsonl"
    index.write_text(json.dumps(entry) + "\n")
    paper = {"gene": "SCN5A", "pmid": "32533946", "llm_trace_refs": [dict(entry)]}
    audit = {"error": error}
    budget = {"calls": []}
    return tmp_path, paper, audit, budget, record, entry, path, index


def test_genuine_guard_refusal_has_no_dispatched_usage(case):
    folder, paper, audit, budget, _, entry, _, _ = case
    assert documented_budget_non_dispatch(paper, audit, folder, budget) == entry


def test_timeout_never_qualifies_as_zero(case):
    folder, paper, audit, budget, *_ = case
    audit["error"] = "Timeout: APITimeoutError - Request timed out."
    assert documented_budget_non_dispatch(paper, audit, folder, budget) is None


def test_existing_reservation_prevents_zero_even_if_failed(case):
    folder, paper, audit, budget, *_ = case
    budget["calls"].append(
        {
            "arm": "grok43_astra_clinical_verified",
            "gene": "SCN5A",
            "pmid": "32533946",
            "status": "failed_usage_unknown_reservation_retained",
        }
    )
    assert documented_budget_non_dispatch(paper, audit, folder, budget) is None


def test_tampered_refusal_is_rejected(case):
    folder, paper, audit, budget, _, _, path, _ = case
    path.write_text("{}")
    assert documented_budget_non_dispatch(paper, audit, folder, budget) is None


def test_missing_paper_reference_is_rejected(case):
    folder, paper, audit, budget, *_ = case
    paper["llm_trace_refs"] = []
    assert documented_budget_non_dispatch(paper, audit, folder, budget) is None


def test_indexed_api_call_prevents_zero(case):
    folder, paper, audit, budget, _, _, _, index = case
    root = folder / "llm_traces"
    path = root / "call.json"
    record = {
        "record_type": "llm_call",
        "trace_id": "api-1",
        "context": {"gene": "SCN5A", "pmid": "32533946", "stage": "paper_curation"},
        "response": {"success": False, "usage": None},
    }
    path.write_text(json.dumps(record))
    entry = {
        **record,
        "path": path.name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    with index.open("a") as handle:
        handle.write(json.dumps(entry) + "\n")
    assert documented_budget_non_dispatch(paper, audit, folder, budget) is None


def test_forged_context_fails_even_with_valid_digest(case):
    folder, paper, audit, budget, record, entry, path, index = case
    record["context"]["pmid"] = "123"
    path.write_text(json.dumps(record))
    entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    index.write_text(json.dumps(entry) + "\n")
    paper["llm_trace_refs"] = [entry]
    assert documented_budget_non_dispatch(paper, audit, folder, budget) is None


def test_refusal_with_wrong_error_cannot_establish_non_dispatch(case):
    folder, paper, audit, budget, *_ = case
    audit["error"] += " forged"
    assert documented_budget_non_dispatch(paper, audit, folder, budget) is None


def test_mislabelled_index_cannot_hide_an_actual_paper_call(case):
    folder, paper, audit, budget, _, _, _, index = case
    path = folder / "llm_traces/call.json"
    record = {
        "record_type": "llm_call",
        "trace_id": "api-1",
        "context": {"gene": "SCN5A", "pmid": "32533946", "stage": "paper_curation"},
        "response": {"success": False, "usage": None},
    }
    path.write_text(json.dumps(record))
    entry = {
        "record_type": "llm_call",
        "trace_id": "api-1",
        "context": {"gene": "SCN5A", "pmid": "different"},
        "path": path.name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    with index.open("a") as handle:
        handle.write(json.dumps(entry) + "\n")
    assert documented_budget_non_dispatch(paper, audit, folder, budget) is None
