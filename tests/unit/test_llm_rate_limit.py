"""Tests for model-aware LLM rate-limit selection."""

from config.settings import get_settings
from utils.llm_utils import _resolve_rpm_limit


def test_rpm_resolution_uses_actual_model_provider(monkeypatch):
    monkeypatch.setenv("NCBI_EMAIL", "user@example.org")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-anthropic")
    monkeypatch.setenv("MODEL_PROVIDER", "anthropic")
    monkeypatch.setenv("AZURE_RPM", "7")
    monkeypatch.setenv("ANTHROPIC_RPM", "211")
    monkeypatch.delenv("LLM_REQUESTS_PER_MINUTE", raising=False)
    get_settings.cache_clear()

    try:
        assert _resolve_rpm_limit("azure_ai/grok-4-20-reasoning") == 7
        assert _resolve_rpm_limit("anthropic/claude-haiku-4-5-20251001") == 211
    finally:
        get_settings.cache_clear()
