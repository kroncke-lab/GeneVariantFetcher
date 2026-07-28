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
from typing import List

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
    base = normalize_azure_ai_api_base()
    key = os.environ.get("AZURE_AI_API_KEY", "")
    if not base or not key:
        raise RuntimeError(
            "Responses API figure extraction requires AZURE_AI_API_BASE and "
            "AZURE_AI_API_KEY"
        )

    url = azure_responses_api_url(base)
    body = {
        "model": _strip_provider_prefix(model),
        "input": [
            {
                "role": "user",
                "content": [
                    {"type": "input_text", "text": _EXTRACT_PROMPT},
                    {"type": "input_image", "image_url": image_data_url},
                ],
            }
        ],
        # GPT-5 family spends part of this budget on reasoning before visible
        # output, so keep it higher than the chat-completions max_tokens.
        "max_output_tokens": 4096,
        **build_responses_reasoning_param(
            model, get_settings().vision_reasoning_effort
        ),
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
            return response.json()
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                f"Responses API returned invalid JSON: {response.text[:300]}"
            ) from exc

    data, _trace = capture_llm_call(
        provider="azure_openai_responses_http",
        requested_model=model,
        resolved_model=_strip_provider_prefix(model),
        request={"endpoint": url, "body": body},
        call=send_request,
    )

    for item in data.get("output", []) or []:
        if item.get("type") != "message":
            continue
        for content in item.get("content", []) or []:
            if content.get("type") == "output_text":
                return (content.get("text") or "").strip()
    return ""
