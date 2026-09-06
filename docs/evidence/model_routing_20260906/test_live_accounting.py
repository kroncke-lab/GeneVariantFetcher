from reconcile_live_budget import reconcile


def test_standard_context_astra_keeps_three_full_attempts_without_premium():
    bound = 100000
    original = 6 * (bound * 10 + 32000 * 50) / 1e6
    row = {
        "model": "azure_ai/gpt-6-astra",
        "status": "reserved",
        "reserved_usd": original,
        "output_cap": 32000,
    }
    data = {"calls": [row]}
    reconcile(data)
    reconcile(data)
    assert row["original_reserved_usd"] == original
    assert row["reserved_usd"] == 3 * (bound * 10 + 32000 * 50) / 1e6


def test_possible_long_context_keeps_premium():
    original = 6 * (300000 * 10 + 32000 * 50) / 1e6
    row = {
        "model": "azure_ai/gpt-6-astra",
        "status": "reserved",
        "reserved_usd": original,
        "output_cap": 32000,
    }
    reconcile({"calls": [row]})
    assert row["reserved_usd"] == original


def test_grok_separate_reasoning_is_billed_and_never_reduced():
    row = {
        "model": "azure_ai/grok-4.6",
        "status": "returned",
        "reserved_usd": 1,
        "api_proxy_usd": 0.00008,
        "accounted_usd": 0.00008,
        "usage": {"prompt_tokens": 10, "completion_tokens": 10, "total_tokens": 1010},
    }
    reconcile({"calls": [row]})
    assert row["billable_output_tokens_conservative"] == 1000
    assert row["accounted_usd"] == (10 * 2 + 1000 * 6) / 1e6
