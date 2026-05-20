"""Tests for model-provider-aware settings resolution."""

from config.settings import Settings


def _settings(**overrides):
    base = {
        "ncbi_email": "user@test.org",
        "anthropic_api_key": "anthropic-key",
        "model_provider": "anthropic",
    }
    base.update(overrides)
    return Settings(_env_file=None, **base)


def test_explicit_tier3_models_win_over_anthropic_provider_default():
    settings = _settings(tier3_models="azure_ai/grok-4-20-reasoning")

    assert settings.get_tier3_models() == ["azure_ai/grok-4-20-reasoning"]


def test_anthropic_provider_uses_default_when_tier3_unset(monkeypatch):
    monkeypatch.delenv("TIER3_MODELS", raising=False)
    settings = _settings()

    assert settings.get_tier3_models() == [
        "anthropic/claude-sonnet-4-6",
        "anthropic/claude-opus-4-7",
    ]
