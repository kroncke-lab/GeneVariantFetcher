"""Tests for provider edge cases in utils.llm_utils."""

from utils import llm_utils
from utils.llm_utils import BaseLLMCaller


def _response(content: str):
    class _Choice:
        def __init__(self, value: str) -> None:
            self.message = type("M", (), {"content": value})
            self.finish_reason = "stop"

    class _Response:
        def __init__(self, value: str) -> None:
            self.choices = [_Choice(value)]

    return _Response(content)


def test_call_llm_json_retries_empty_content_once(monkeypatch):
    calls: list[str] = []

    def fake_completion(**_kwargs):
        calls.append("call")
        if len(calls) == 1:
            return _response("")
        return _response('{"ok": true}')

    monkeypatch.setattr(llm_utils, "completion", fake_completion)
    monkeypatch.setattr(llm_utils, "wait_for_llm_rate_limit", lambda _model: None)

    caller = BaseLLMCaller(model="azure_ai/Kimi-K2.6-1", max_tokens=8192)

    assert caller.call_llm_json("return json") == {"ok": True}
    assert len(calls) == 2
