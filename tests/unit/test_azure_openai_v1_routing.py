"""Regression coverage for Azure Foundry's OpenAI-compatible v1 endpoint."""

from utils.llm_utils import resolve_litellm_model_and_kwargs


def test_openai_v1_base_rewrites_azure_model_and_preserves_credentials(
    monkeypatch,
):
    monkeypatch.setenv(
        "AZURE_AI_API_BASE",
        "https://example.services.ai.azure.com/openai/v1/",
    )
    monkeypatch.setenv("AZURE_AI_API_KEY", "test-azure-key")

    model, kwargs = resolve_litellm_model_and_kwargs(
        "azure_ai/gpt-5.6-sol",
        temperature=0,
        max_tokens=8192,
    )

    assert model == "openai/gpt-5.6-sol"
    assert kwargs["api_base"] == ("https://example.services.ai.azure.com/openai/v1")
    assert kwargs["api_key"] == "test-azure-key"
    assert "temperature" not in kwargs
    assert kwargs["max_tokens"] == 8192


def test_resource_root_keeps_standard_azure_ai_route(monkeypatch):
    monkeypatch.setenv("AZURE_AI_API_BASE", "https://example.services.ai.azure.com")
    monkeypatch.setenv("AZURE_AI_API_KEY", "test-azure-key")

    model, kwargs = resolve_litellm_model_and_kwargs(
        "azure_ai/gpt-5.4",
        temperature=0,
    )

    assert model == "azure_ai/gpt-5.4"
    assert kwargs == {"temperature": 0}
