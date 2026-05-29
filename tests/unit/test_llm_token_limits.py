from utils.llm_utils import (
    build_reasoning_effort_kwargs,
    build_responses_reasoning_param,
    clamp_max_tokens,
)


def test_new_azure_deployments_have_large_output_budget():
    assert clamp_max_tokens("azure_ai/gpt-5.4", 60000, warn=False) == 15000
    assert clamp_max_tokens("azure_ai/DeepSeek-V4-Pro", 60000, warn=False) == 15000


def test_gpt5_and_grok43_families_have_large_output_budget():
    # gpt-5.5 was previously unrecognized and fell to the 4000 default.
    assert clamp_max_tokens("azure_ai/gpt-5.5", 60000, warn=False) == 15000
    # Generic gpt-5 fallback catches future gpt-5.x ids we have not enumerated.
    assert clamp_max_tokens("gpt-5.9", 60000, warn=False) == 15000
    # grok-4.3 is covered by both the explicit entry and the grok-4 substring.
    assert clamp_max_tokens("azure_ai/grok-4.3", 60000, warn=False) == 15000


def test_reasoning_effort_kwargs_gate():
    # OpenAI-style reasoning models get the param passed through.
    assert build_reasoning_effort_kwargs("azure_ai/gpt-5.5", "high") == {
        "reasoning_effort": "high"
    }
    assert build_reasoning_effort_kwargs("o3-mini", "low") == {
        "reasoning_effort": "low"
    }
    # Anthropic (extended thinking, not wired yet) and Grok (no knob) -> no-op.
    assert build_reasoning_effort_kwargs("anthropic/claude-sonnet-4-6", "high") == {}
    assert build_reasoning_effort_kwargs("azure_ai/grok-4.3", "high") == {}
    # Falsy effort -> provider default (empty kwargs).
    assert build_reasoning_effort_kwargs("gpt-5.5", None) == {}


def test_responses_reasoning_param_gate():
    # Responses-API vision models (gpt-5 family) get a nested reasoning fragment.
    assert build_responses_reasoning_param("azure_ai/gpt-5.3-codex-1", "high") == {
        "reasoning": {"effort": "high"}
    }
    assert build_responses_reasoning_param("gpt-5.3-codex-1", "minimal") == {
        "reasoning": {"effort": "minimal"}
    }
    # Non-OpenAI model or no effort -> no fragment (body unchanged).
    assert build_responses_reasoning_param("azure_ai/grok-4.3", "high") == {}
    assert build_responses_reasoning_param("gpt-5.5", None) == {}
