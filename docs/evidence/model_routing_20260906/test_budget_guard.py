"""Offline campaign accounting checks; no provider requests."""

import json
import pytest
import budget_guard


@pytest.fixture
def guard(tmp_path, monkeypatch):
    import utils.llm_trace as trace
    import utils.llm_utils as llm

    path = tmp_path / "budget.json"
    path.write_text(
        json.dumps(
            {
                "calls": [],
                "api_ceiling_usd": 150,
                "per_call_reserve_limit_usd": 35,
                "smoke_uncertainty_reserve_usd": 5,
            }
        )
    )
    monkeypatch.setenv("GVF_EXPERIMENT_BUDGET", str(path))
    original = trace.capture_llm_call
    monkeypatch.setattr(trace, "capture_llm_call", original)
    monkeypatch.setattr(llm, "capture_llm_call", original)
    budget_guard.install()
    assert llm.capture_llm_call is trace.capture_llm_call
    return path, trace.capture_llm_call


def invoke(fn, call, request=None):
    return fn(
        provider="mock",
        requested_model="azure_ai/gpt-6-astra",
        resolved_model="openai/gpt-6-astra",
        request=request or {"max_completion_tokens": 1000, "messages": []},
        call=call,
    )


def test_release_and_keep_usage(guard):
    path, fn = guard
    response = {
        "usage": {"prompt_tokens": 100, "completion_tokens": 20},
        "model": "gpt-6-astra",
    }
    assert invoke(fn, lambda: response)[0] == response
    row = json.loads(path.read_text())["calls"][0]
    assert row["accounted_usd"] == pytest.approx(0.002)
    assert row["reserved_usd"] > row["accounted_usd"]


def test_unknown_failure_retains_reservation(guard):
    path, fn = guard

    def fail():
        raise RuntimeError("network")

    with pytest.raises(RuntimeError, match="network"):
        invoke(fn, fail)
    row = json.loads(path.read_text())["calls"][0]
    assert row["status"] == "failed_usage_unknown_reservation_retained"
    assert "accounted_usd" not in row


def test_cap_and_ceiling_refuse_before_dispatch(guard):
    path, fn = guard

    def forbidden():
        raise AssertionError("must not dispatch")

    with pytest.raises(RuntimeError, match="explicit output"):
        invoke(fn, forbidden, {"messages": []})
    data = json.loads(path.read_text())
    data["api_ceiling_usd"] = 5
    path.write_text(json.dumps(data))
    with pytest.raises(RuntimeError, match="budget prevents"):
        invoke(fn, forbidden)
    assert json.loads(path.read_text())["calls"] == []


def test_http_body_cap(guard):
    path, fn = guard
    invoke(
        fn,
        lambda: {"usage": {"input_tokens": 10, "output_tokens": 5}},
        {"body": {"max_output_tokens": 1000}},
    )
    assert json.loads(path.read_text())["calls"][0]["status"] == "returned"


def test_sdk_retry_keeps_unknown_attempt_cost_after_success(guard):
    import logging

    path, fn = guard

    def retry_then_return():
        logging.getLogger("openai._base_client").info(
            "Retrying request to /chat/completions in 0.4 seconds"
        )
        return {
            "usage": {
                "prompt_tokens": 100,
                "completion_tokens": 20,
                "total_tokens": 120,
            }
        }

    invoke(fn, retry_then_return)
    row = json.loads(path.read_text())["calls"][0]
    assert row["sdk_retry_count"] == 1
    assert row["unknown_retry_reserve_usd"] == pytest.approx(row["reserved_usd"] / 3)
    assert row["accounted_usd"] > row["api_proxy_usd"]
