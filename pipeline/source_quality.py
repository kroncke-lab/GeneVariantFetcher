"""Helpers for deciding whether harvested markdown is usable as full text."""

from __future__ import annotations

from pathlib import Path

from config.constants import MIN_EXTRACTION_INPUT_SIZE


SUPPLEMENT_SURFACE_STATUS_MARKER = "GVF_SUPPLEMENT_SURFACE_STATUS:"
INCOMPLETE_SUPPLEMENT_SURFACE_STATUSES = frozenset({"unavailable", "scrape_error"})


def is_abstract_only_fallback_text(text: str) -> bool:
    """Return True for GVF's generated abstract-only fallback markdown."""
    head = text[:8192].lower()
    return "# abstract-only fallback" in head or (
        "full text could not be retrieved" in head
        and "contains only the pubmed abstract" in head
    )


def is_usable_fulltext_source(path: Path) -> bool:
    """Return True when a markdown source is eligible for full-text extraction."""
    try:
        if not path.exists() or path.stat().st_size == 0:
            return False
        with path.open("r", encoding="utf-8", errors="replace") as f:
            head = f.read(8192)
        if len(head) < MIN_EXTRACTION_INPUT_SIZE:
            return False
        return not is_abstract_only_fallback_text(head)
    except OSError:
        return False


def has_incomplete_supplement_surface(text: str) -> bool:
    """Return whether GVF explicitly recorded that supplements were uninspected."""

    head = text[:8192]
    for line in head.splitlines():
        if SUPPLEMENT_SURFACE_STATUS_MARKER not in line:
            continue
        status = line.split(SUPPLEMENT_SURFACE_STATUS_MARKER, 1)[1]
        status = status.strip(" -!><\t\r\n").casefold()
        return status in INCOMPLETE_SUPPLEMENT_SURFACE_STATUSES
    return False


def is_reusable_fulltext_source(path: Path) -> bool:
    """Return True only when source is usable and not known body-only/incomplete."""

    if not is_usable_fulltext_source(path):
        return False
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            return not has_incomplete_supplement_surface(handle.read(8192))
    except OSError:
        return False
