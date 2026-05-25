"""Helpers for deciding whether harvested markdown is usable as full text."""

from __future__ import annotations

from pathlib import Path


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
            return not is_abstract_only_fallback_text(f.read(8192))
    except OSError:
        return False
