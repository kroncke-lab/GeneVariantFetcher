"""Vision-LLM text extraction from supplement/figure images.

Disabled by default. Enable by setting the ``GVF_EXTRACT_FIGURE_TEXT``
environment variable to a truthy value (``1``, ``true``, ``yes``) before
running the pipeline, or pass an injectable ``image_text_extractor`` callable
to :func:`harvesting.paywall_context_enrichment.enrich_paywall_full_context`
for testing.

No network calls are made unless the feature is explicitly enabled.
"""

from __future__ import annotations

import base64
import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List

import requests

from config.settings import get_settings
from utils.llm_trace import (
    OUTCOME_DISCARDED,
    OUTCOME_ERROR,
    OUTCOME_PARSED,
    attempt_link_summary,
    capture_llm_call,
    last_llm_trace,
    llm_attempt_ledger,
    llm_trace_scope,
    note_llm_attempt,
    note_llm_outcome,
    record_trace_event,
)
from utils.llm_utils import (
    azure_responses_api_url,
    build_reasoning_effort_kwargs,
    build_responses_reasoning_param,
    litellm_completion,
    normalize_azure_ai_api_base,
)

logger = logging.getLogger(__name__)

_IMAGE_SUFFIXES = frozenset(
    {".png", ".jpg", ".jpeg", ".gif", ".tiff", ".tif", ".webp", ".bmp"}
)

_RESPONSES_API_PREFIXES = (
    "gpt-5",
    "azure_ai/gpt-5",
)

_EXTRACT_PROMPT = """\
This image is from a biomedical journal article (supplement or figure).

Extract ALL text and tabular data you can see. If there is a table, reproduce
it in Markdown table format. If there is free text or labels, quote them
verbatim. If there is no useful text (e.g. a purely decorative figure), reply
with an empty string.

Return ONLY the extracted content — no preamble, no commentary."""


def is_image_path(path: "str | Path") -> bool:
    """Return True when *path* has a recognised image file extension."""
    return Path(path).suffix.lower() in _IMAGE_SUFFIXES


def _uses_responses_api(model: str) -> bool:
    name = (model or "").lower()
    return any(name.startswith(prefix) for prefix in _RESPONSES_API_PREFIXES)


def _strip_provider_prefix(model: str) -> str:
    if model.startswith("azure_ai/"):
        return model[len("azure_ai/") :]
    return model


def extract_images_to_markdown(
    image_paths: List[Path],
    model: str,
    *,
    gene: "str | None" = None,
    pmid: "str | None" = None,
) -> str:
    """Call a vision LLM on each image and return a combined markdown section.

    Returns an empty string when ``image_paths`` is empty — no API call is
    made. Individual image failures are logged and skipped; they do not abort
    the batch.

    ``gene``/``pmid`` come from the caller because a figure OCR prompt contains
    neither: without them every figure read for every paper collapsed into one
    ``_unscoped/`` group in the trace tree.
    """
    if not image_paths:
        return ""

    parts: List[str] = []
    for idx, img_path in enumerate(image_paths, start=1):
        try:
            text = _extract_one(img_path, model, gene=gene, pmid=pmid)
        except Exception as exc:
            logger.warning(
                "Figure text extraction failed for %s: %s", img_path.name, exc
            )
            text = ""
        if text and text.strip():
            parts.append(f"### Figure Image {idx}: {img_path.name}\n\n{text.strip()}\n")

    if not parts:
        return ""

    return "\n\n## FIGURE IMAGE TEXT\n\n" + "\n".join(parts)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _image_to_data_url(image_path: Path) -> str:
    ext = image_path.suffix.lower()
    mime = {
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".png": "image/png",
        ".gif": "image/gif",
        ".webp": "image/webp",
        ".tiff": "image/tiff",
        ".tif": "image/tiff",
        ".bmp": "image/bmp",
    }.get(ext, "image/png")
    b64 = base64.b64encode(image_path.read_bytes()).decode()
    return f"data:{mime};base64,{b64}"


def _extract_one(
    image_path: Path,
    model: str,
    *,
    gene: "str | None" = None,
    pmid: "str | None" = None,
) -> str:
    data_url = _image_to_data_url(image_path)
    responses_path = _uses_responses_api(model)
    with (
        llm_trace_scope(
            gene=gene,
            pmid=pmid,
            stage="figure_text_extraction",
            component="figure_text_extractor",
            operation="responses_api" if responses_path else "chat_completions",
            image_path=str(image_path),
        ),
        llm_attempt_ledger(),
    ):
        trace_id: "str | None" = None
        try:
            if responses_path:
                text = _extract_one_responses_api(data_url, model)
                trace_id = note_llm_attempt(last_llm_trace(), role="figure_ocr")
            else:
                response = litellm_completion(
                    model=model,
                    messages=[
                        {
                            "role": "user",
                            "content": [
                                {"type": "text", "text": _EXTRACT_PROMPT},
                                {"type": "image_url", "image_url": {"url": data_url}},
                            ],
                        }
                    ],
                    temperature=0,
                    max_tokens=2048,
                    **build_reasoning_effort_kwargs(
                        model, get_settings().vision_reasoning_effort
                    ),
                )
                trace_id = note_llm_attempt(last_llm_trace(), role="figure_ocr")
                text = (response.choices[0].message.content or "").strip()
        except BaseException as exc:
            if trace_id is None:
                trace_id = note_llm_attempt(last_llm_trace(), role="figure_ocr")
            note_llm_outcome(trace_id, OUTCOME_ERROR)
            _record_figure_text_decision(
                image_path,
                model,
                characters=0,
                decision_source="call_failed",
                error=f"{type(exc).__name__}: {exc}",
            )
            raise
        # An empty extraction is a legitimate outcome for a decorative figure —
        # it must stay distinguishable from a provider failure and from real text.
        note_llm_outcome(trace_id, OUTCOME_PARSED if text else OUTCOME_DISCARDED)
        _record_figure_text_decision(
            image_path,
            model,
            characters=len(text),
            decision_source="text_extracted" if text else "empty_output",
        )
        return text


def _record_figure_text_decision(
    image_path: Path,
    model: str,
    *,
    characters: int,
    decision_source: str,
    error: "str | None" = None,
) -> None:
    record_trace_event(
        "figure_text_extraction_decision",
        {
            "image": image_path.name,
            "image_path": str(image_path),
            "model": model,
            "characters_extracted": characters,
            "text_extracted": bool(characters),
            "decision_source": decision_source,
            "error": error,
            "selection_policy": (
                "OCR one figure image. An empty result means the figure carried "
                "no usable text; it is not a failure and adds nothing to the "
                "unified source."
            ),
            **attempt_link_summary(),
        },
    )


def _extract_one_responses_api(image_data_url: str, model: str) -> str:
    """Call Azure AI Foundry Responses API for GPT-5-family vision models."""
    return call_responses_api_vision(
        _EXTRACT_PROMPT, image_data_url, model, attempt_role="figure_ocr"
    )


def normalize_responses_usage(data: Any) -> Any:
    """Alias Responses-API usage counters to the chat-completions names in-place.

    ``capture_llm_call`` stores ``response.usage`` verbatim, and cost tooling
    sums the litellm path's ``prompt_tokens``/``completion_tokens`` — the
    Responses API reports ``input_tokens``/``output_tokens`` instead, so every
    raw-HTTP vision call traced as zero cost. Counters the provider did not
    report stay absent: a missing usage block must never become zeros.
    """
    if not isinstance(data, dict):
        return data
    usage = data.get("usage")
    if not isinstance(usage, dict):
        return data
    tokens_in = usage.get("input_tokens")
    tokens_out = usage.get("output_tokens")
    if isinstance(tokens_in, int) and "prompt_tokens" not in usage:
        usage["prompt_tokens"] = tokens_in
    if isinstance(tokens_out, int) and "completion_tokens" not in usage:
        usage["completion_tokens"] = tokens_out
    if (
        "total_tokens" not in usage
        and isinstance(tokens_in, int)
        and isinstance(tokens_out, int)
    ):
        usage["total_tokens"] = tokens_in + tokens_out
    details = usage.get("output_tokens_details")
    reasoning = details.get("reasoning_tokens") if isinstance(details, dict) else None
    if isinstance(reasoning, int) and "completion_tokens_details" not in usage:
        usage["completion_tokens_details"] = {"reasoning_tokens": reasoning}
    return data


def responses_output_text(data: Any) -> str:
    """First ``output_text`` block in a Responses-API envelope, stripped."""
    if not isinstance(data, dict):
        return ""
    for item in data.get("output", []) or []:
        if not isinstance(item, dict) or item.get("type") != "message":
            continue
        for content in item.get("content", []) or []:
            if isinstance(content, dict) and content.get("type") == "output_text":
                return (content.get("text") or "").strip()
    return ""


def responses_incomplete_reason(data: Any) -> "str | None":
    """The provider's reason when the envelope stopped early, else None."""
    if not isinstance(data, dict) or data.get("status") != "incomplete":
        return None
    details = data.get("incomplete_details")
    reason = details.get("reason") if isinstance(details, dict) else None
    return str(reason) if reason else "incomplete"


def _post_responses_api(
    url: str,
    key: str,
    model: str,
    prompt: str,
    image_data_url: str,
    reasoning_param: Dict[str, Any],
) -> dict:
    body = {
        "model": _strip_provider_prefix(model),
        "input": [
            {
                "role": "user",
                "content": [
                    {"type": "input_text", "text": prompt},
                    {"type": "input_image", "image_url": image_data_url},
                ],
            }
        ],
        # GPT-5 family spends part of this budget on reasoning before visible
        # output, so keep it higher than the chat-completions max_tokens.
        "max_output_tokens": 4096,
        **reasoning_param,
    }

    def send_request() -> dict:
        response = requests.post(
            url,
            headers={"api-key": key, "Content-Type": "application/json"},
            json=body,
            timeout=120,
        )
        if response.status_code != 200:
            raise RuntimeError(
                f"Responses API returned {response.status_code}: {response.text[:300]}"
            )
        try:
            payload = response.json()
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                f"Responses API returned invalid JSON: {response.text[:300]}"
            ) from exc
        # Aliased before capture so the trace's usage carries the counter
        # names cost tooling sums.
        return normalize_responses_usage(payload)

    data, _trace = capture_llm_call(
        provider="azure_openai_responses_http",
        requested_model=model,
        resolved_model=_strip_provider_prefix(model),
        request={"endpoint": url, "body": body},
        call=send_request,
    )
    return data


def call_responses_api_vision(
    prompt: str,
    image_data_url: str,
    model: str,
    *,
    attempt_role: str,
) -> str:
    """POST one vision prompt to the Azure Responses API and return its text.

    gpt-5-family models spend part of ``max_output_tokens`` on reasoning before
    any visible output; when the whole budget goes to reasoning the envelope
    comes back ``status=incomplete`` with no output text. That is a token-cap
    burn, not a clean "nothing to report": retry once at reasoning effort
    ``low``, and if the retry is also empty accept it — with the outcome
    recorded in the trace — rather than re-burning the cap at the same effort.
    """
    base = normalize_azure_ai_api_base()
    key = os.environ.get("AZURE_AI_API_KEY", "")
    if not base or not key:
        raise RuntimeError(
            "Responses API figure extraction requires AZURE_AI_API_BASE and "
            "AZURE_AI_API_KEY"
        )
    url = azure_responses_api_url(base)
    configured_effort = get_settings().vision_reasoning_effort
    reasoning_param = build_responses_reasoning_param(model, configured_effort)

    data = _post_responses_api(url, key, model, prompt, image_data_url, reasoning_param)
    text = responses_output_text(data)
    incomplete = responses_incomplete_reason(data)
    if text or not incomplete:
        return text

    retry_param = build_responses_reasoning_param(model, "low")
    if not retry_param or retry_param == reasoning_param:
        # Already at low, or the model has no effort knob: a retry would
        # re-burn the same cap for the same nothing.
        record_trace_event(
            "figure_vision_token_cap",
            {
                "outcome": "accepted_empty_no_retry",
                "incomplete_reason": incomplete,
                "reasoning_effort": configured_effort,
                "model": model,
            },
        )
        return ""
    first_trace_id = note_llm_attempt(last_llm_trace(), role=attempt_role)
    note_llm_outcome(first_trace_id, OUTCOME_DISCARDED)
    record_trace_event(
        "figure_vision_token_cap",
        {
            "outcome": "retry_low_effort",
            "incomplete_reason": incomplete,
            "reasoning_effort": configured_effort,
            "retry_reasoning_effort": "low",
            "model": model,
        },
    )
    data = _post_responses_api(url, key, model, prompt, image_data_url, retry_param)
    text = responses_output_text(data)
    if not text:
        record_trace_event(
            "figure_vision_token_cap",
            {
                "outcome": "accepted_empty_after_retry",
                "incomplete_reason": responses_incomplete_reason(data),
                "retry_reasoning_effort": "low",
                "model": model,
            },
        )
    return text
