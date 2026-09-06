"""Regression coverage for Azure Foundry's OpenAI-compatible v1 endpoint."""

from utils.llm_utils import resolve_litellm_model_and_kwargs
import pytest
import json


def test_new_models_can_use_a_separate_resource_without_moving_old_stages(monkeypatch):
    monkeypatch.setenv(
        "AZURE_AI_API_BASE", "https://old.services.ai.azure.com/openai/v1"
    )
    monkeypatch.setenv("AZURE_AI_API_KEY", "old-key")
    monkeypatch.setenv("NEW_RESOURCE_KEY", "new-key")
    monkeypatch.setenv(
        "AZURE_AI_MODEL_ROUTES",
        json.dumps(
            {
                "gpt-6-astra": {
                    "api_base": "https://new.services.ai.azure.com/openai/v1",
                    "api_key_env": "NEW_RESOURCE_KEY",
                }
            }
        ),
    )
    _, old = resolve_litellm_model_and_kwargs("azure_ai/grok-4.3")
    _, new = resolve_litellm_model_and_kwargs("azure_ai/gpt-6-astra")
    assert old["api_key"] == "old-key"
    assert new["api_key"] == "new-key"
    assert old["api_base"] != new["api_base"]
    monkeypatch.delenv("NEW_RESOURCE_KEY")
    with pytest.raises(ValueError, match="requires NEW_RESOURCE_KEY"):
        resolve_litellm_model_and_kwargs("azure_ai/gpt-6-astra")


@pytest.mark.parametrize("temperature", [0, 1, None])
def test_astra_omits_all_unsupported_sampling_parameters(monkeypatch, temperature):
    monkeypatch.setenv(
        "AZURE_AI_API_BASE", "https://example.services.ai.azure.com/openai/v1"
    )
    model, kwargs = resolve_litellm_model_and_kwargs(
        "azure_ai/gpt-6-astra",
        temperature=temperature,
        top_p=0.9,
        logprobs=False,
        top_logprobs=2,
        max_tokens=32000,
        reasoning_effort="high",
    )
    assert model == "openai/gpt-6-astra"
    assert not {"temperature", "top_p", "logprobs", "top_logprobs"}.intersection(kwargs)
    assert "max_tokens" not in kwargs
    assert kwargs["max_completion_tokens"] == 32000
    assert kwargs["reasoning_effort"] == "high"


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


@pytest.mark.parametrize("model", ["gpt-6-astra", "grok-4.6"])
@pytest.mark.parametrize("effort", ["medium", "high"])
def test_new_model_effort_reaches_actual_sdk_request(monkeypatch, model, effort):
    import httpx
    from openai import OpenAI
    from utils.llm_utils import litellm_completion

    sent = []

    def handle(request):
        sent.append(json.loads(request.content))
        return httpx.Response(
            200,
            json={
                "id": "wire-test",
                "object": "chat.completion",
                "created": 1,
                "model": model,
                "choices": [
                    {
                        "index": 0,
                        "message": {"role": "assistant", "content": "{}"},
                        "finish_reason": "stop",
                    }
                ],
                "usage": {
                    "prompt_tokens": 10,
                    "completion_tokens": 2,
                    "total_tokens": 12,
                },
            },
        )

    client = OpenAI(
        api_key="test-only-key",
        base_url="https://mock.invalid/v1",
        http_client=httpx.Client(transport=httpx.MockTransport(handle)),
        max_retries=0,
    )
    litellm_completion(
        model="openai/" + model,
        client=client,
        messages=[{"role": "user", "content": "Return JSON"}],
        max_tokens=32000,
        reasoning_effort=effort,
        temperature=0,
        response_format={"type": "json_object"},
    )
    assert len(sent) == 1
    assert sent[0]["reasoning_effort"] == effort
    cap = "max_completion_tokens" if model == "gpt-6-astra" else "max_tokens"
    assert sent[0][cap] == 32000
    if model == "gpt-6-astra":
        assert "temperature" not in sent[0]


@pytest.mark.parametrize(
    "route",
    [
        None,
        {
            "api_base": "https://new.services.ai.azure.com",
            "api_key_env": "NEW_RESOURCE_KEY",
        },
    ],
)
def test_selected_malformed_route_never_falls_back_to_old_resource(monkeypatch, route):
    monkeypatch.setenv("AZURE_AI_MODEL_ROUTES", json.dumps({"gpt-6-astra": route}))
    monkeypatch.setenv("NEW_RESOURCE_KEY", "test-key")
    with pytest.raises(ValueError):
        resolve_litellm_model_and_kwargs("azure_ai/gpt-6-astra")


@pytest.mark.parametrize("entry", ["figure_text", "figure_variant", "pedigree"])
def test_astra_responses_paths_use_routed_resource_and_response_fields(
    monkeypatch, entry
):
    from types import SimpleNamespace
    import harvesting.figure_text_extractor as text_reader
    import harvesting.figure_variant_reader as variant_reader
    import pipeline.pedigree_extractor as pedigree

    monkeypatch.setenv(
        "AZURE_AI_API_BASE", "https://old.services.ai.azure.com/openai/v1"
    )
    monkeypatch.setenv("AZURE_AI_API_KEY", "old-key")
    monkeypatch.setenv("NEW_RESOURCE_KEY", "new-key")
    monkeypatch.setenv(
        "AZURE_AI_MODEL_ROUTES",
        json.dumps(
            {
                "gpt-6-astra": {
                    "api_base": "https://new.services.ai.azure.com/openai/v1",
                    "api_key_env": "NEW_RESOURCE_KEY",
                }
            }
        ),
    )
    settings = SimpleNamespace(vision_reasoning_effort="medium")
    monkeypatch.setattr(text_reader, "get_settings", lambda: settings)
    monkeypatch.setattr(pedigree, "get_settings", lambda: settings)
    sent = []

    def post(url, **kwargs):
        sent.append((url, kwargs))
        return SimpleNamespace(
            status_code=200,
            json=lambda: {
                "status": "completed",
                "output": [
                    {
                        "type": "message",
                        "content": [{"type": "output_text", "text": '{"ok":true}'}],
                    }
                ],
                "usage": {"input_tokens": 10, "output_tokens": 4, "total_tokens": 14},
            },
        )

    monkeypatch.setattr(text_reader.requests, "post", post)
    if entry == "pedigree":
        result = pedigree._call_azure_responses_api_vision(
            deployment="gpt-6-astra",
            prompt="Read this",
            image_data_url="data:image/png;base64,AAAA",
            max_output_tokens=4096,
        )
        assert result == {"ok": True}
    elif entry == "figure_variant":
        assert (
            variant_reader._call_responses_api(
                "data:image/png;base64,AAAA", "azure_ai/gpt-6-astra", "Read this"
            )
            == '{"ok":true}'
        )
    else:
        assert (
            text_reader.call_responses_api_vision(
                "Read this",
                "data:image/png;base64,AAAA",
                "azure_ai/gpt-6-astra",
                attempt_role="test",
            )
            == '{"ok":true}'
        )
    assert len(sent) == 1
    url, request = sent[0]
    assert url == "https://new.services.ai.azure.com/openai/v1/responses?api-version=v1"
    assert request["headers"]["api-key"] == "new-key"
    body = request["json"]
    assert body["model"] == "gpt-6-astra"
    assert body["max_output_tokens"] == 4096
    assert body["reasoning"] == {"effort": "medium"}
    assert not {"temperature", "max_tokens", "max_completion_tokens"}.intersection(body)
