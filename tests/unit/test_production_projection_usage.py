"""A failed production call must not become a zero-cost exact total."""

import pytest
from benchmarks.codex_paper_eval.db_to_predictions import (
    add_usage,
    empty_usage,
    finalized_usage,
)


def record(identity, usage, *, success=True, error=None):
    return {
        "trace_id": identity,
        "context": {"model": "azure_ai/gpt-6-astra"},
        "response": {"usage": usage, "success": success, "error": error},
    }


def test_mixed_returned_and_failed_usage_retains_known_subset_and_unknown_total():
    bucket = empty_usage()
    add_usage(
        bucket,
        record("ok", {"input_tokens": 100, "output_tokens": 2, "total_tokens": 132}),
    )
    add_usage(bucket, record("timeout", None, success=False, error="timeout"))
    final = finalized_usage(bucket)
    assert final["telemetry_available"] is False
    assert (
        final["input_tokens"] is final["output_tokens"] is final["total_tokens"] is None
    )
    assert final["known_usage"] == {
        "input_tokens": 100,
        "output_tokens": 32,
        "total_tokens": 132,
    }
    assert final["unknown_failed_call_trace_ids"] == ["timeout"]
    model = final["models"]["azure_ai/gpt-6-astra"]
    assert model["total_tokens"] is None
    assert model["known_usage"] == final["known_usage"]
    assert final["llm_calls"] == 2 and final["successful_calls"] == 1
    assert bucket["total_tokens"] == 132  # Export does not mutate the accumulator.


def test_returned_only_usage_and_no_call_usage_are_exact():
    bucket = empty_usage()
    assert finalized_usage(bucket)["total_tokens"] == 0
    assert finalized_usage(bucket)["telemetry_available"] is True
    add_usage(
        bucket,
        record("ok", {"prompt_tokens": 3, "completion_tokens": 4, "total_tokens": 7}),
    )
    final = finalized_usage(bucket)
    assert final["telemetry_available"] is True
    assert final["total_tokens"] == 7
    assert "known_usage" not in final


@pytest.mark.parametrize("success,error", [(True, None), (False, None)])
def test_unexplained_missing_usage_cannot_claim_exact_zero(success, error):
    with pytest.raises(ValueError, match="traceable failed API call"):
        add_usage(empty_usage(), record("invalid", None, success=success, error=error))
