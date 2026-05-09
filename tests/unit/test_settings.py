"""Tests for configuration validation."""

import pytest

from config.settings import get_settings


@pytest.fixture(autouse=True)
def clear_settings_cache():
    """Ensure the settings cache is cleared before and after each test."""
    get_settings.cache_clear()
    yield
    get_settings.cache_clear()


def test_missing_required_settings(monkeypatch):
    """Missing required values should raise a clear validation error."""

    for env_var in ["OPENAI_API_KEY", "NCBI_EMAIL"]:
        monkeypatch.setenv(env_var, "")

    with pytest.raises(ValueError) as exc:
        get_settings()

    assert "Missing required settings" in str(exc.value)


def test_placeholder_values_are_rejected(monkeypatch):
    """Placeholder values should be caught by validation."""

    monkeypatch.setenv("OPENAI_API_KEY", "your_api_key")
    monkeypatch.setenv("NCBI_EMAIL", "your_email@example.com")

    with pytest.raises(ValueError) as exc:
        get_settings()

    assert "placeholder values" in str(exc.value)


def test_valid_settings_are_loaded(monkeypatch):
    """Valid settings should be returned with defaults applied."""

    monkeypatch.setenv("OPENAI_API_KEY", "sk-openai")
    monkeypatch.setenv("NCBI_EMAIL", "user@example.org")
    monkeypatch.setenv("INTERN_MODEL", "gpt-4o-mini")

    settings = get_settings()

    assert settings.openai_api_key == "sk-openai"
    assert settings.ncbi_email == "user@example.org"
    assert settings.intern_model == "gpt-4o-mini"


def _baseline_env(monkeypatch, provider: str):
    """Set the minimum env vars required for Settings to load.

    Settings reads `.env` directly via SettingsConfigDict(env_file=...),
    so `monkeypatch.delenv` cannot mask values that live in the file.
    To make the provider-default tests independent of the developer's
    .env, we explicitly setenv the per-provider knobs to the documented
    defaults — the tests are about the resolution logic, not the values.
    """
    monkeypatch.setenv("NCBI_EMAIL", "user@example.org")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-anthropic")
    monkeypatch.setenv("MODEL_PROVIDER", provider)
    # Clear explicit overrides so provider-aware paths run.
    monkeypatch.delenv("LLM_REQUESTS_PER_MINUTE", raising=False)
    monkeypatch.delenv("MAX_WORKERS", raising=False)
    # Pin provider knobs to their settings.py defaults so the test result
    # is independent of the loaded .env (which can override these).
    monkeypatch.setenv("AZURE_RPM", "50")
    monkeypatch.setenv("ANTHROPIC_RPM", "200")
    monkeypatch.setenv("AZURE_MAX_WORKERS", "3")
    monkeypatch.setenv("ANTHROPIC_MAX_WORKERS", "10")


def test_provider_aware_rpm_anthropic(monkeypatch):
    _baseline_env(monkeypatch, "anthropic")
    settings = get_settings()
    assert settings.get_requests_per_minute() == 200
    assert settings.get_max_workers() == 10


def test_provider_aware_rpm_azure(monkeypatch):
    _baseline_env(monkeypatch, "azure")
    settings = get_settings()
    assert settings.get_requests_per_minute() == 50
    assert settings.get_max_workers() == 3


def test_explicit_overrides_win_over_provider_defaults(monkeypatch):
    _baseline_env(monkeypatch, "anthropic")
    monkeypatch.setenv("LLM_REQUESTS_PER_MINUTE", "75")
    monkeypatch.setenv("MAX_WORKERS", "5")
    settings = get_settings()
    assert settings.get_requests_per_minute() == 75
    assert settings.get_max_workers() == 5
