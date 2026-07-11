"""Unit tests for the gpt-5 / Responses API routing in pedigree_extractor."""

import json
from unittest.mock import patch, MagicMock

from pipeline.pedigree_extractor import (
    _call_azure_responses_api_vision,
    _strip_provider_prefix,
    _uses_responses_api,
)


def test_uses_responses_api_detects_gpt5_family():
    assert _uses_responses_api("gpt-5") is True
    assert _uses_responses_api("gpt-5-codex") is True
    assert _uses_responses_api("gpt-5.3-codex-1") is True
    assert _uses_responses_api("gpt-5.6-sol") is True
    assert _uses_responses_api("azure_ai/gpt-5") is True
    assert _uses_responses_api("azure_ai/gpt-5.3-codex-1") is True
    assert _uses_responses_api("azure_ai/gpt-5.6-sol") is True


def test_uses_responses_api_rejects_other_models():
    assert _uses_responses_api("gpt-4o") is False
    assert _uses_responses_api("gpt-4.1") is False
    assert _uses_responses_api("azure_ai/Kimi-K2.6-1") is False
    assert _uses_responses_api("claude-3-5-sonnet") is False
    assert _uses_responses_api("gemini-1.5-pro") is False


def test_strip_provider_prefix():
    assert _strip_provider_prefix("azure_ai/gpt-5.3-codex-1") == "gpt-5.3-codex-1"
    assert _strip_provider_prefix("gpt-5.3-codex-1") == "gpt-5.3-codex-1"
    assert _strip_provider_prefix("azure_ai/Kimi-K2.6-1") == "Kimi-K2.6-1"


def _make_fake_response(status_code: int, body: dict | str) -> MagicMock:
    fake = MagicMock()
    fake.status_code = status_code
    if isinstance(body, dict):
        fake.json.return_value = body
        fake.text = json.dumps(body)
    else:
        fake.text = body
        fake.json.side_effect = json.JSONDecodeError("not json", body, 0)
    return fake


def test_responses_api_parses_output_text():
    expected_json = {"is_pedigree": True, "confidence": 0.92}
    api_payload = {
        "id": "resp_abc",
        "object": "response",
        "output": [
            {
                "type": "message",
                "content": [
                    {"type": "output_text", "text": json.dumps(expected_json)},
                ],
            }
        ],
        "usage": {"total_tokens": 42},
    }
    env = {"AZURE_AI_API_BASE": "https://h.example", "AZURE_AI_API_KEY": "k"}
    with patch.dict("os.environ", env, clear=False):
        with patch(
            "pipeline.pedigree_extractor.requests.post",
            return_value=_make_fake_response(200, api_payload),
        ) as post:
            result = _call_azure_responses_api_vision(
                deployment="gpt-5.3-codex-1",
                prompt="Is this a pedigree?",
                image_data_url="data:image/png;base64,AAAA",
                max_output_tokens=1024,
            )
    assert result == expected_json
    call_args, call_kwargs = post.call_args
    url = call_kwargs.get("url") or (call_args[0] if call_args else "")
    assert "v1/responses" in url
    body = call_kwargs["json"]
    assert body["model"] == "gpt-5.3-codex-1"
    # Input must use the Responses-API content shape
    user_content = body["input"][0]["content"]
    types = [c["type"] for c in user_content]
    assert "input_text" in types and "input_image" in types


def test_responses_api_returns_none_on_500():
    env = {"AZURE_AI_API_BASE": "https://h.example", "AZURE_AI_API_KEY": "k"}
    with patch.dict("os.environ", env, clear=False):
        with patch(
            "pipeline.pedigree_extractor.requests.post",
            return_value=_make_fake_response(500, "boom"),
        ):
            result = _call_azure_responses_api_vision(
                deployment="gpt-5.3-codex-1",
                prompt="?",
                image_data_url="data:image/png;base64,AAAA",
                max_output_tokens=1024,
            )
    assert result is None


def test_responses_api_returns_none_when_no_output_text():
    api_payload = {"output": [{"type": "reasoning"}]}
    env = {"AZURE_AI_API_BASE": "https://h.example", "AZURE_AI_API_KEY": "k"}
    with patch.dict("os.environ", env, clear=False):
        with patch(
            "pipeline.pedigree_extractor.requests.post",
            return_value=_make_fake_response(200, api_payload),
        ):
            result = _call_azure_responses_api_vision(
                deployment="gpt-5.3-codex-1",
                prompt="?",
                image_data_url="data:image/png;base64,AAAA",
                max_output_tokens=1024,
            )
    assert result is None


def test_responses_api_requires_env_vars():
    with patch.dict(
        "os.environ", {"AZURE_AI_API_BASE": "", "AZURE_AI_API_KEY": ""}, clear=False
    ):
        result = _call_azure_responses_api_vision(
            deployment="gpt-5.3-codex-1",
            prompt="?",
            image_data_url="data:image/png;base64,AAAA",
            max_output_tokens=1024,
        )
    assert result is None
