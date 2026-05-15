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
from litellm import completion

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


def extract_images_to_markdown(image_paths: List[Path], model: str) -> str:
    """Call a vision LLM on each image and return a combined markdown section.

    Returns an empty string when ``image_paths`` is empty — no API call is
    made. Individual image failures are logged and skipped; they do not abort
    the batch.
    """
    if not image_paths:
        return ""

    parts: List[str] = []
    for idx, img_path in enumerate(image_paths, start=1):
        try:
            text = _extract_one(img_path, model)
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


def _extract_one(image_path: Path, model: str) -> str:
    data_url = _image_to_data_url(image_path)
    if _uses_responses_api(model):
        return _extract_one_responses_api(data_url, model)

    response = completion(
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
    )
    return (response.choices[0].message.content or "").strip()


def _extract_one_responses_api(image_data_url: str, model: str) -> str:
    """Call Azure AI Foundry Responses API for GPT-5-family vision models."""
    base = os.environ.get("AZURE_AI_API_BASE", "").rstrip("/")
    key = os.environ.get("AZURE_AI_API_KEY", "")
    if not base or not key:
        raise RuntimeError(
            "Responses API figure extraction requires AZURE_AI_API_BASE and "
            "AZURE_AI_API_KEY"
        )

    url = f"{base}/openai/v1/responses?api-version=v1"
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
    }
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
        data = response.json()
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            f"Responses API returned invalid JSON: {response.text[:300]}"
        ) from exc

    for item in data.get("output", []) or []:
        if item.get("type") != "message":
            continue
        for content in item.get("content", []) or []:
            if content.get("type") == "output_text":
                return (content.get("text") or "").strip()
    return ""
