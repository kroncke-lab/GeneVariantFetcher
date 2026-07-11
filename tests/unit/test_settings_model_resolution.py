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


def test_azure_provider_uses_deployment_defaults_when_tiers_unset(monkeypatch):
    for env_var in (
        "TIER2_MODEL",
        "TIER3_MODELS",
        "TIER3_ADJUDICATOR_MODELS",
        "TABLE_ROUTER_MODEL",
    ):
        monkeypatch.delenv(env_var, raising=False)

    settings = _settings(
        model_provider="azure",
        azure_deployment_kimi="Kimi-K2.6-1",
        azure_deployment_grok="grok-4.3",
    )

    assert settings.get_tier2_model() == "azure_ai/gpt-5.4"
    assert settings.get_table_router_model() == "azure_ai/Kimi-K2.6-1"
    assert settings.get_tier3_models() == ["azure_ai/grok-4.3"]
    assert settings.get_tier3_adjudicator_models() == ["azure_ai/gpt-5.4"]


def test_experimental_azure_models_are_prefixed():
    settings = _settings(
        azure_deployment_gpt54="gpt-5.4",
        azure_deployment_deepseek="DeepSeek-V4-Pro",
    )

    assert settings.get_experimental_azure_models() == [
        "azure_ai/gpt-5.4",
        "azure_ai/DeepSeek-V4-Pro",
    ]


def test_experimental_azure_models_accept_existing_prefix():
    settings = _settings(
        azure_deployment_gpt54="azure_ai/gpt-5.4",
        azure_deployment_deepseek="",
    )

    assert settings.get_experimental_azure_models() == ["azure_ai/gpt-5.4"]


def test_default_early_debate_models_include_gpt54_and_deepseek(monkeypatch):
    monkeypatch.delenv("EARLY_DEBATE_MODELS", raising=False)
    settings = _settings(
        early_debate_models="",
        azure_deployment_gpt54="gpt-5.4",
        azure_deployment_deepseek="DeepSeek-V4-Pro",
    )

    assert settings.get_early_debate_models() == [
        "azure_ai/gpt-5.4",
        "azure_ai/DeepSeek-V4-Pro",
        "azure_ai/Kimi-K2.6-1",
    ]


def test_explicit_early_debate_models_win():
    settings = _settings(
        early_debate_models=("azure_ai/grok-4-20-reasoning,anthropic/claude-sonnet-4-6")
    )

    assert settings.get_early_debate_models() == [
        "azure_ai/grok-4-20-reasoning",
        "anthropic/claude-sonnet-4-6",
    ]


def test_final_adjudication_models_are_anthropic_reserved_defaults(monkeypatch):
    # Isolate from a developer's live .env / exported FINAL_* vars.
    for env_var in (
        "FINAL_ADJUDICATOR_MODELS",
        "FINAL_ARBITER_MODEL",
        "FINAL_ADJUDICATOR_REASONING_EFFORT",
        "FINAL_ARBITER_REASONING_EFFORT",
    ):
        monkeypatch.delenv(env_var, raising=False)

    settings = _settings()

    assert settings.get_final_adjudicator_models() == ["anthropic/claude-sonnet-5"]
    assert settings.get_final_arbiter_model() == "anthropic/claude-opus-4-8"


def test_explicit_final_adjudication_models_win():
    settings = _settings(
        final_adjudicator_models=(
            "anthropic/claude-sonnet-5,anthropic/claude-opus-4-8"
        ),
        final_arbiter_model="anthropic/claude-opus-4-8-max",
    )

    assert settings.get_final_adjudicator_models() == [
        "anthropic/claude-sonnet-5",
        "anthropic/claude-opus-4-8",
    ]
    assert settings.get_final_arbiter_model() == "anthropic/claude-opus-4-8-max"


def test_final_adjudication_accepts_gpt56_sol_with_max_effort_alias():
    settings = _settings(
        final_adjudicator_models="azure_ai/gpt-5.6-sol",
        final_adjudicator_reasoning_effort="max",
        final_arbiter_model="azure_ai/gpt-5.6-sol",
        final_arbiter_reasoning_effort="xhigh",
    )

    assert settings.get_final_adjudicator_models() == ["azure_ai/gpt-5.6-sol"]
    assert settings.get_final_arbiter_model() == "azure_ai/gpt-5.6-sol"
    # "max" normalizes to Azure-supported xhigh
    assert settings.final_adjudicator_reasoning_effort == "xhigh"
    assert settings.final_arbiter_reasoning_effort == "xhigh"


def test_gpt56_deployment_alias_is_isolated_to_paper_final_check(monkeypatch):
    for env_var in (
        "PAPER_FINAL_CHECK_MODEL",
        "TIER2_MODEL",
        "FINAL_ADJUDICATOR_MODELS",
        "FINAL_ARBITER_MODEL",
    ):
        monkeypatch.delenv(env_var, raising=False)
    settings = _settings(
        model_provider="azure",
        azure_deployment_gpt54="gpt-5.4",
        azure_deployment_gpt56_sol="gpt-5.6-sol-regional",
    )

    assert settings.get_tier2_model() == "azure_ai/gpt-5.4"
    assert settings.get_final_adjudicator_models() == ["anthropic/claude-sonnet-5"]
    assert settings.get_final_arbiter_model() == "anthropic/claude-opus-4-8"
    assert settings.get_paper_final_check_model() == "azure_ai/gpt-5.6-sol-regional"


def test_explicit_paper_final_check_model_wins_over_deployment_alias():
    settings = _settings(
        azure_deployment_gpt56_sol="gpt-5.6-sol-regional",
        paper_final_check_model="anthropic/claude-sonnet-5",
    )

    assert settings.get_paper_final_check_model() == "anthropic/claude-sonnet-5"


def test_gpt56_deployment_environment_alias_configures_paper_check(monkeypatch):
    monkeypatch.delenv("PAPER_FINAL_CHECK_MODEL", raising=False)
    monkeypatch.setenv("AZURE_DEPLOYMENT_GPT56_SOL", "gpt-5.6-sol-environment")

    settings = _settings()

    assert settings.get_paper_final_check_model() == "azure_ai/gpt-5.6-sol-environment"
