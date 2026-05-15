"""Enrich a recovered paywall full-text context with captions and supplements.

When ``scripts/fetch_paywalled.py`` rescues a paywalled article (either via a
Tier 3.5 publisher strategy or the Europe PMC fallback), the basic
``FULL_CONTEXT.md`` we used to write was just the body markdown produced by
the DOM extractor. That dropped two evidence sources that the canonical PMC
harvest path preserves:

  1. Figure / table / supplementary-material *captions* (extracted by
     ``harvesting.figure_extractor``). Variants frequently appear here
     ("Table 2. 86 KCNH2 missense variants identified...", "Pedigree of
     family carrying KCNH2 G604S").
  2. Downloaded *supplement files* (PDF, Excel, Word, HTML), which often
     contain the largest variant tables.

This module gives ``fetch_paywalled.py`` (and any future paywall recovery
caller) one place to assemble those signals into the same unified layout the
main orchestrator emits — main text, caption block, supplement markdown —
so downstream extraction sees the same shape regardless of how the
``FULL_CONTEXT.md`` was harvested.

The enricher never modifies the *body* content that the quality gate
validated; captions and supplements are appended after it, so the gate's
verdict on the rescued article is unaffected.
"""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional
from urllib.parse import urljoin, urlparse

from .figure_extractor import (
    CaptionExtractionResult,
    extract_from_html,
    merge_results,
    render_captions_markdown,
    save_captions_json,
)
from .figure_text_extractor import extract_images_to_markdown, is_image_path
from .supplement_processing_service import (
    SupplementFileResult,
    process_supplement_files,
)

logger = logging.getLogger(__name__)


# Extensions whose body content we can convert to markdown — anything else is
# pointless to download for variant extraction (icons, fonts, video etc.).
_CONVERTIBLE_EXTS = (
    ".pdf",
    ".xlsx",
    ".xls",
    ".docx",
    ".doc",
    ".csv",
    ".tsv",
    ".txt",
    ".html",
    ".htm",
    ".xml",
    ".zip",
)

# Don't pull down monster files — large supplement zips occasionally contain
# multi-gigabyte raw sequencing data that won't help us and would blow the
# disk. Body-text supplements (variant lists) are virtually always <25 MB.
_DEFAULT_SUPP_SIZE_LIMIT_BYTES = 25 * 1024 * 1024
_IMAGE_EXTS = (".jpg", ".jpeg", ".png", ".gif", ".webp", ".tif", ".tiff", ".bmp")


@dataclass
class EnrichmentResult:
    """What the enricher produced for one PMID."""

    unified_markdown: str
    captions: CaptionExtractionResult = field(default_factory=CaptionExtractionResult)
    supplement_markdown: str = ""
    supplement_results: List[SupplementFileResult] = field(default_factory=list)
    captions_path: Optional[Path] = None

    @property
    def figure_caption_count(self) -> int:
        return len(self.captions.figures)

    @property
    def table_caption_count(self) -> int:
        return len(self.captions.tables)

    @property
    def supplement_count(self) -> int:
        return sum(1 for r in self.supplement_results if r.downloaded)


def _supplement_filename(entry: Dict[str, Any], idx: int) -> str:
    """Pick a safe filename from a supplement link record."""
    raw = (entry.get("name") or "").strip() or f"supplement_{idx}"
    # Strip path separators and querystring residue so we never escape the
    # supplements directory.
    safe = raw.replace("/", "_").replace("\\", "_").split("?")[0].split("#")[0]
    return safe or f"supplement_{idx}"


def _looks_convertible(filename: str) -> bool:
    lower = filename.lower()
    return lower.endswith(_CONVERTIBLE_EXTS)


def _safe_image_label(value: str, fallback: str) -> str:
    raw = (value or "").strip() or fallback
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", raw).strip("._-")
    return safe[:80] or fallback


def _caption_image_entries(captions: CaptionExtractionResult) -> List[tuple[str, str]]:
    """Return ``(label, image_url)`` pairs from figure and image-table captions."""
    entries: List[tuple[str, str]] = []
    for fig in captions.figures:
        if fig.image_url:
            entries.append((fig.label or fig.figure_id or "figure", fig.image_url))
    for tbl in captions.tables:
        image_url = getattr(tbl, "image_url", None)
        if image_url:
            entries.append((tbl.label or tbl.table_id or "table", image_url))
    return entries


def _download_caption_images(
    *,
    captions: CaptionExtractionResult,
    output_dir: Path,
    pmid: str,
    session: Any,
    source_url: Optional[str],
    size_limit_bytes: int,
    timeout_s: int,
    max_images: int = 12,
) -> List[Path]:
    """Download caption-associated article/table images for optional OCR."""
    if session is None or captions.is_empty():
        return []

    image_paths: List[Path] = []
    images_dir = output_dir / f"{pmid}_figures" / "html_images"
    seen_urls: set[str] = set()

    for idx, (label, raw_url) in enumerate(
        _caption_image_entries(captions)[:max_images], start=1
    ):
        resolved = urljoin(source_url or "", raw_url)
        parsed = urlparse(resolved)
        if parsed.scheme not in {"http", "https"}:
            continue
        if resolved in seen_urls:
            continue
        seen_urls.add(resolved)

        url_name = Path(parsed.path).name
        suffix = Path(url_name).suffix.lower()
        if suffix not in _IMAGE_EXTS:
            suffix = ".jpg"
        safe_label = _safe_image_label(label, f"image_{idx}")
        filename = f"{idx:02d}_{safe_label}{suffix}"
        file_path = images_dir / filename

        if _default_download(
            url=resolved,
            file_path=file_path,
            session=session,
            size_limit_bytes=size_limit_bytes,
            timeout_s=timeout_s,
        ):
            image_paths.append(file_path)

    return image_paths


def _default_download(
    *,
    url: str,
    file_path: Path,
    session: Any,
    size_limit_bytes: int,
    timeout_s: int,
) -> bool:
    """Stream a supplement to disk via the shared requests session.

    Returns True on a successful, size-bounded write. Failures (HTTP errors,
    oversized payloads, write errors) return False without raising — the
    caller logs the failure as a per-file ``download_failed`` row so the rest
    of enrichment can continue.
    """
    if session is None:
        return False
    try:
        with session.get(url, stream=True, timeout=timeout_s) as resp:
            if resp.status_code != 200:
                logger.info(
                    "Supplement download non-200: %s -> %s", url, resp.status_code
                )
                return False
            total = 0
            file_path.parent.mkdir(parents=True, exist_ok=True)
            with open(file_path, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=64 * 1024):
                    if not chunk:
                        continue
                    total += len(chunk)
                    if total > size_limit_bytes:
                        logger.info(
                            "Supplement exceeded size limit %d bytes: %s",
                            size_limit_bytes,
                            url,
                        )
                        try:
                            file_path.unlink(missing_ok=True)
                        except OSError:
                            pass
                        return False
                    fh.write(chunk)
            return total > 0
    except Exception as exc:
        logger.info("Supplement download error for %s: %s", url, exc)
        try:
            file_path.unlink(missing_ok=True)
        except OSError:
            pass
        return False


def enrich_paywall_full_context(
    *,
    body_markdown: str,
    html: Optional[str],
    supp_files: Optional[List[Dict[str, Any]]],
    pmid: str,
    output_dir: Path,
    converter: Any,
    session: Any = None,
    extra_captions: Optional[CaptionExtractionResult] = None,
    download_supplements: bool = True,
    max_supplements: int = 12,
    supplement_size_limit_bytes: int = _DEFAULT_SUPP_SIZE_LIMIT_BYTES,
    download_timeout_s: int = 60,
    image_text_extractor: Optional[Callable[[List[Path]], str]] = None,
    source_url: Optional[str] = None,
) -> EnrichmentResult:
    """Append caption block and supplement markdown to a rescued body.

    Args:
        body_markdown: The body content the quality gate already validated.
        html: Raw HTML of the rescued page, used for caption extraction. May
            be None if only body markdown is available.
        supp_files: List of ``{url, name, ...}`` dicts produced by a Tier 3.5
            strategy or supplement scraper.
        pmid: PubMed ID — used to scope sidecar directories.
        output_dir: Run output directory; supplements land under
            ``<output_dir>/<pmid>_supplements/`` and the captions index lands
            under ``<output_dir>/<pmid>_figures/captions_index.json``.
        converter: Instance of :class:`harvesting.format_converters.FormatConverter`.
        session: ``requests.Session`` used to download supplement files.
            If None, supplement download is skipped and only the link list is
            preserved in the caption block (via ``extra_captions``).
        extra_captions: Pre-computed captions to merge with anything we
            extract from ``html`` (e.g. JATS XML captions from an earlier
            harvest path).
        download_supplements: Set to False to skip the download attempt
            entirely (useful in tests).
        max_supplements: Cap on how many supplements we try to download.
        image_text_extractor: Optional callable ``(image_paths: List[Path]) ->
            str`` that extracts text from figure images. When *None* (default),
            extraction runs only if the ``GVF_EXTRACT_FIGURE_TEXT`` environment
            variable is set to a truthy value (``1``, ``true``, or ``yes``), in
            which case :func:`harvesting.figure_text_extractor.extract_images_to_markdown`
            is called with the configured vision model. Supply a stub here in
            tests to avoid real API calls.

    Returns:
        :class:`EnrichmentResult` with the assembled markdown and audit
        metadata for the caller.
    """
    body = body_markdown or ""

    captions = CaptionExtractionResult()
    if html:
        try:
            captions = extract_from_html(html)
        except Exception as exc:
            logger.warning("Caption extraction failed for PMID %s: %s", pmid, exc)
            captions = CaptionExtractionResult()

    if extra_captions is not None and not extra_captions.is_empty():
        captions = merge_results(captions, extra_captions)

    # Persist a captions sidecar — even if the unified markdown ends up tiny,
    # the captions JSON is useful for audit and downstream tooling.
    captions_path: Optional[Path] = None
    if not captions.is_empty():
        try:
            figures_dir = output_dir / f"{pmid}_figures"
            figures_dir.mkdir(parents=True, exist_ok=True)
            captions_path = figures_dir / "captions_index.json"
            save_captions_json(captions, captions_path)
        except Exception as exc:
            logger.warning(
                "Failed to write captions_index.json for PMID %s: %s", pmid, exc
            )
            captions_path = None

    supplement_markdown = ""
    supplement_results: List[SupplementFileResult] = []
    if download_supplements and supp_files:
        # Skip records that don't carry a URL or whose filename clearly isn't
        # a body-text supplement (image-only previews, SVG icons, etc).
        usable: List[Dict[str, Any]] = []
        for entry in supp_files[:max_supplements]:
            url = (entry.get("url") or "").strip()
            if not url:
                continue
            name = _supplement_filename(entry, len(usable) + 1)
            if not _looks_convertible(name):
                # Some publishers expose supplements with no explicit
                # extension — keep them only if the URL itself has one.
                if not _looks_convertible(url):
                    continue
            normalized = dict(entry)
            normalized["url"] = url
            normalized["name"] = name
            usable.append(normalized)

        if usable:
            supplements_dir = output_dir / f"{pmid}_supplements"

            def _cb(
                url: str,
                file_path: Path,
                _pmid: str,
                _filename: str,
                _supp: Dict[str, Any],
            ) -> bool:
                return _default_download(
                    url=url,
                    file_path=file_path,
                    session=session,
                    size_limit_bytes=supplement_size_limit_bytes,
                    timeout_s=download_timeout_s,
                )

            try:
                result = process_supplement_files(
                    supp_files=usable,
                    supplements_dir=supplements_dir,
                    pmid=pmid,
                    converter=converter,
                    download_callback=_cb,
                    extract_figures=True,
                    figures_dir=output_dir / f"{pmid}_figures",
                    logger=logger,
                    sleep_seconds=0.0,
                )
                supplement_markdown = result.supplement_markdown
                supplement_results = result.file_results
            except Exception as exc:
                logger.warning(
                    "Supplement processing crashed for PMID %s: %s", pmid, exc
                )

    # -----------------------------------------------------------------------
    # Optional vision pass: extract text from article/table images and from
    # figure images embedded in supplement files (e.g. image tables extracted
    # from PDFs or ZIPs).
    # Runs only when GVF_EXTRACT_FIGURE_TEXT is truthy OR an injectable
    # extractor is supplied (for tests).  Default is off — no API calls.
    # -----------------------------------------------------------------------
    figure_image_text = ""
    _env_enabled = os.environ.get("GVF_EXTRACT_FIGURE_TEXT", "").lower() in (
        "1",
        "true",
        "yes",
    )
    if image_text_extractor is not None or _env_enabled:
        image_paths: List[Path] = _download_caption_images(
            captions=captions,
            output_dir=output_dir,
            pmid=pmid,
            session=session,
            source_url=source_url,
            size_limit_bytes=supplement_size_limit_bytes,
            timeout_s=download_timeout_s,
        )
        for sfr in supplement_results:
            for nf in sfr.nested_files:
                p = Path(nf)
                if is_image_path(p) and p.exists():
                    image_paths.append(p)

        if image_paths:
            try:
                if image_text_extractor is not None:
                    figure_image_text = image_text_extractor(image_paths)
                else:
                    from config.settings import get_settings

                    model = get_settings().get_vision_model()
                    figure_image_text = extract_images_to_markdown(image_paths, model)
            except Exception as exc:
                logger.warning(
                    "Figure image text extraction failed for PMID %s: %s", pmid, exc
                )

    captions_md = render_captions_markdown(captions)
    unified = (
        body
        + (captions_md or "")
        + (supplement_markdown or "")
        + (figure_image_text or "")
    )

    return EnrichmentResult(
        unified_markdown=unified,
        captions=captions,
        supplement_markdown=supplement_markdown,
        supplement_results=supplement_results,
        captions_path=captions_path,
    )
