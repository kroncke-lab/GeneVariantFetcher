"""Vision-LLM API for extracting variants from figure images.

The existing :mod:`harvesting.figure_text_extractor` does a generic OCR pass.
This module is the focused counterpart: a *variant-aware* reader that
prompts a multimodal model to find HGVS-style mutations, optionally with
their reported carrier/affected/unaffected counts, and returns structured
JSON.

Typical use cases:

  * A figure is a *mutation map* (linear topology diagram with circles at
    each variant position labeled ``R176W``, ``L552S`` etc.).
  * A figure or supplement is an *image-only table* listing variants and
    cohort counts the publisher's HTML stripped.
  * A pedigree figure whose caption names a single variant
    (``Family carrying KCNH2 G604S``).

The reader is gene-scoped — every call must pass the gene symbol so the
prompt anchors on it and avoids false positives from other genes mentioned
in the same image.

Wire-up notes:

  * Reads paths off the per-PMID ``_artifacts.json`` and per-figures
    directory the harvester already creates.
  * Vision model is resolved from ``config.settings.Settings.get_vision_model``
    unless the caller passes ``--model`` / a model argument.
  * Returns a :class:`FigureReadResult` per image plus an aggregated list
    of distinct variants for convenience.
"""

from __future__ import annotations

import base64
import json
import logging
import os
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from config.settings import get_settings
from utils.llm_trace import (
    OUTCOME_DISCARDED,
    OUTCOME_ERROR,
    OUTCOME_PARSE_FAILED,
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

_RESPONSES_API_PREFIXES = ("gpt-5", "azure_ai/gpt-5")


_VARIANT_PROMPT = """\
You are reading a figure from a biomedical journal article about variants in \
the {gene} gene. Extract every {gene} variant the image shows.

For each variant return one JSON object with these keys (omit any you cannot \
read):

  "protein"         protein-level HGVS or short form, e.g. "p.Arg176Trp",
                    "R176W", "G601S", "K897fs", "R1014*"
  "cdna"            cDNA-level HGVS if visible, e.g. "c.526C>T", "c.2515del"
  "carriers"        total carriers reported for this variant (integer)
  "affected"        affected carriers (integer)
  "unaffected"      unaffected carriers (integer)
  "context"         short string saying where in the image it appears
                    (e.g. "Table 2 row 4", "transmembrane S6", "Figure 1 black circle")

Rules:
  * ONLY return {gene} variants. If a variant clearly belongs to a different
    gene shown in the image (e.g. KCNQ1 in a multi-gene compendium), skip it.
  * If the image is a mutation map / topology diagram, treat every labeled
    residue as a variant even if no counts are given.
  * Do NOT copy the carrier total onto affected. affected and unaffected are
    only people the image marks as such (filled vs empty pedigree symbols,
    or explicit Aff/Unaff columns). If you cannot count those groups
    separately, omit affected and unaffected.
  * Never emit unaffected=0 just because you did not see an unaffected count.
  * Single-letter and three-letter forms are both acceptable.
  * If the image has no {gene} variants, return an empty array.

Return JSON ONLY. Top-level: {{"variants": [ ... ]}}. No commentary, no
markdown fences, no preamble."""


@dataclass
class FigureReadResult:
    """One figure's worth of variant data."""

    image_path: str
    variants: List[Dict[str, Any]] = field(default_factory=list)
    error: Optional[str] = None

    @property
    def ok(self) -> bool:
        return self.error is None


@dataclass
class PMIDFigureReport:
    """Aggregated result for one PMID."""

    pmid: str
    gene: str
    per_figure: List[FigureReadResult] = field(default_factory=list)

    @property
    def distinct_variants(self) -> List[Dict[str, Any]]:
        """Deduplicate variants by (protein, cdna) keeping the first sighting."""
        seen: Dict[tuple, Dict[str, Any]] = {}
        for fig in self.per_figure:
            for v in fig.variants:
                key = (
                    (v.get("protein") or "").strip().lower(),
                    (v.get("cdna") or "").strip().lower(),
                )
                if key == ("", ""):
                    continue
                seen.setdefault(key, v)
        return list(seen.values())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pmid": self.pmid,
            "gene": self.gene,
            "per_figure": [asdict(f) for f in self.per_figure],
            "distinct_variants": self.distinct_variants,
        }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def is_image_path(path: "str | Path") -> bool:
    return Path(path).suffix.lower() in _IMAGE_SUFFIXES


def find_pmid_figures(pmc_fulltext_dir: Path, pmid: str) -> List[Path]:
    """Return images from a flat run cache or nested canonical corpus."""
    fig_dir = pmc_fulltext_dir / f"{pmid}_figures"
    if not fig_dir.is_dir():
        fig_dir = pmc_fulltext_dir / pmid / f"{pmid}_figures"
    if not fig_dir.is_dir():
        return []
    imgs = []
    for suffix in _IMAGE_SUFFIXES:
        imgs.extend(sorted(fig_dir.glob(f"*{suffix}")))
    return imgs


def read_figures_for_pmid(
    pmid: str,
    gene: str,
    pmc_fulltext_dir: Path,
    *,
    model: Optional[str] = None,
    max_images: Optional[int] = None,
) -> PMIDFigureReport:
    """Run the variant reader over every figure on disk for *pmid*."""
    images = find_pmid_figures(pmc_fulltext_dir, pmid)
    if max_images is not None:
        images = images[:max_images]
    return read_images(images, gene, pmid=pmid, model=model)


def read_images(
    image_paths: List[Path],
    gene: str,
    *,
    pmid: str = "",
    model: Optional[str] = None,
) -> PMIDFigureReport:
    """Run the reader over an arbitrary list of image paths."""
    model = model or _default_vision_model()
    report = PMIDFigureReport(pmid=pmid, gene=gene)
    for img in image_paths:
        # PMID is in hand here and must be forwarded: _trace_parent needs BOTH
        # gene and pmid, so dropping it collapsed every figure read for every
        # paper into a single gene-level group.
        result = _read_one(img, gene, model, pmid=pmid)
        report.per_figure.append(result)
    return report


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------


def _default_vision_model() -> str:
    try:
        from config.settings import get_settings

        return get_settings().get_vision_model()
    except Exception:
        return os.environ.get("VISION_MODEL", "anthropic/claude-sonnet-4-6")


def _uses_responses_api(model: str) -> bool:
    name = (model or "").lower()
    return any(name.startswith(prefix) for prefix in _RESPONSES_API_PREFIXES)


def _strip_provider_prefix(model: str) -> str:
    if model.startswith("azure_ai/"):
        return model[len("azure_ai/") :]
    return model


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


def _parse_response(text: str) -> List[Dict[str, Any]]:
    """Coerce the model output to a list of variant dicts.

    Accepts JSON object with ``variants`` key, raw JSON list, or fenced JSON.
    Returns ``[]`` on parse failure rather than raising.
    """
    return _parse_response_status(text)[0]


def _parse_response_status(text: str) -> tuple[List[Dict[str, Any]], str]:
    """``(variants, status)`` where status is ``ok``/``empty_output``/``unparseable``.

    The biomedical outcome is unchanged — an unparseable response still yields no
    variants — but "the model returned nothing about this figure" and "the model
    returned something we could not read" are different facts for a curator, and
    a bare ``[]`` conflated them.
    """
    if not text:
        return [], "empty_output"
    cleaned = text.strip()
    if not cleaned:
        return [], "empty_output"
    # Strip code fences if present
    if cleaned.startswith("```"):
        cleaned = re.sub(r"^```(?:json)?\s*", "", cleaned)
        cleaned = re.sub(r"\s*```$", "", cleaned)
    try:
        data = json.loads(cleaned)
    except json.JSONDecodeError:
        # Try to find a JSON object inside the text
        m = re.search(r"\{[^{}]*\"variants\"\s*:\s*\[.*?\]\s*\}", cleaned, re.DOTALL)
        if not m:
            logger.warning("Could not parse figure-reader response: %s", cleaned[:200])
            return [], "unparseable"
        try:
            data = json.loads(m.group(0))
        except json.JSONDecodeError:
            return [], "unparseable"
    if isinstance(data, list):
        variants = [v for v in data if isinstance(v, dict)]
        return _sanitize_figure_phenotype_counts(variants), "ok"
    if isinstance(data, dict):
        vars_ = data.get("variants", [])
        if isinstance(vars_, list):
            variants = [v for v in vars_ if isinstance(v, dict)]
            return _sanitize_figure_phenotype_counts(variants), "ok"
    return [], "unparseable"


def _sanitize_figure_phenotype_counts(
    variants: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Drop figure-reader phenotype copies before they reach extraction."""
    from pipeline.phenotype_count_guard import sanitize_copied_phenotype

    for variant in variants:
        variant.setdefault("source_layer", "figure")
        sanitize_copied_phenotype(variant)
    return variants


def _read_one(
    image_path: Path, gene: str, model: str, *, pmid: str = ""
) -> FigureReadResult:
    try:
        data_url = _image_to_data_url(image_path)
    except Exception as exc:
        return FigureReadResult(image_path=str(image_path), error=f"read_image: {exc}")

    prompt = _VARIANT_PROMPT.format(gene=gene)
    responses_path = _uses_responses_api(model)

    with (
        llm_trace_scope(
            gene=gene,
            pmid=pmid or None,
            stage="figure_variant_read",
            component="figure_variant_reader",
            operation="responses_api" if responses_path else "chat_completions",
            image_path=str(image_path),
        ),
        llm_attempt_ledger(),
    ):
        trace_id: Optional[str] = None
        try:
            if responses_path:
                raw = _call_responses_api(data_url, model, prompt)
            else:
                raw = _call_chat_completions(data_url, model, prompt)
            trace_id = note_llm_attempt(last_llm_trace(), role="figure_variant_read")
        except Exception as exc:
            trace_id = note_llm_attempt(last_llm_trace(), role="figure_variant_read")
            note_llm_outcome(trace_id, OUTCOME_ERROR)
            _record_figure_variant_decision(
                image_path,
                gene,
                model,
                variant_count=0,
                decision_source="call_failed",
                error=str(exc),
            )
            return FigureReadResult(image_path=str(image_path), error=str(exc))

        variants, status = _parse_response_status(raw)
        if status == "unparseable":
            note_llm_outcome(trace_id, OUTCOME_PARSE_FAILED)
            decision_source = "unparseable_response"
        elif status == "empty_output":
            note_llm_outcome(trace_id, OUTCOME_DISCARDED)
            decision_source = "empty_output"
        elif not variants:
            # Parsed cleanly and reported no variants: a real, usable answer.
            note_llm_outcome(trace_id, OUTCOME_PARSED)
            decision_source = "no_variants_reported"
        else:
            note_llm_outcome(trace_id, OUTCOME_PARSED)
            decision_source = "variants_read"
        _record_figure_variant_decision(
            image_path,
            gene,
            model,
            variant_count=len(variants),
            decision_source=decision_source,
            variants=[
                str(v.get("variant") or v.get("notation") or "") for v in variants
            ],
        )
        return FigureReadResult(image_path=str(image_path), variants=variants)


def _record_figure_variant_decision(
    image_path: Path,
    gene: str,
    model: str,
    *,
    variant_count: int,
    decision_source: str,
    error: Optional[str] = None,
    variants: Optional[List[str]] = None,
) -> None:
    record_trace_event(
        "figure_variant_read_decision",
        {
            "image": image_path.name,
            "image_path": str(image_path),
            "gene": gene,
            "model": model,
            "variant_count": variant_count,
            "variants": [v for v in (variants or []) if v][:40],
            "decision_source": decision_source,
            "error": error,
            "selection_policy": (
                "Read one figure image for target-gene variants. An unparseable "
                "response and a clean 'no variants' answer both yield zero "
                "variants but are recorded distinctly."
            ),
            **attempt_link_summary(),
        },
    )


def _call_chat_completions(data_url: str, model: str, prompt: str) -> str:
    response = litellm_completion(
        model=model,
        messages=[
            {
                "role": "user",
                "content": [
                    {"type": "text", "text": prompt},
                    {"type": "image_url", "image_url": {"url": data_url}},
                ],
            }
        ],
        temperature=0,
        max_tokens=2048,
        **build_reasoning_effort_kwargs(model, get_settings().vision_reasoning_effort),
    )
    return (response.choices[0].message.content or "").strip()


def _call_responses_api(data_url: str, model: str, prompt: str) -> str:
    base = normalize_azure_ai_api_base()
    key = os.environ.get("AZURE_AI_API_KEY", "")
    if not base or not key:
        raise RuntimeError(
            "Responses API requires AZURE_AI_API_BASE and AZURE_AI_API_KEY"
        )

    url = azure_responses_api_url(base)
    body = {
        "model": _strip_provider_prefix(model),
        "input": [
            {
                "role": "user",
                "content": [
                    {"type": "input_text", "text": prompt},
                    {"type": "input_image", "image_url": data_url},
                ],
            }
        ],
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
