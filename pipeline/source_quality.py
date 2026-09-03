"""Helpers for deciding whether harvested markdown is usable as full text."""

from __future__ import annotations

import logging
from pathlib import Path

from config.constants import MIN_EXTRACTION_INPUT_SIZE

logger = logging.getLogger(__name__)


SUPPLEMENT_SURFACE_STATUS_MARKER = "GVF_SUPPLEMENT_SURFACE_STATUS:"
INCOMPLETE_SUPPLEMENT_SURFACE_STATUSES = frozenset({"unavailable", "scrape_error"})

# A FULL_CONTEXT below this floor is an empty/truncated write, not a paper.
# Deliberately far below MIN_EXTRACTION_INPUT_SIZE: this floor only marks files
# that are outright empty shells so a populated sibling rendering can be
# preferred; it does not change what counts as usable full text.
EMPTY_FULL_CONTEXT_FLOOR_BYTES = 64

# Sibling renderings that can stand in for an empty FULL_CONTEXT, best first.
_FULL_CONTEXT_SIBLING_SUFFIXES = ("_CLEANED.md", "_DATA_ZONES.md")


def is_abstract_only_fallback_text(text: str) -> bool:
    """Return True for GVF's generated abstract-only fallback markdown."""
    head = text[:8192].lower()
    return "# abstract-only fallback" in head or (
        "full text could not be retrieved" in head
        and "contains only the pubmed abstract" in head
    )


def _is_article_shaped_sibling(path: Path) -> tuple[bool, str]:
    """Run the harvest-time content-quality validator over a sibling rendering.

    The length floor plus the abstract-only marker accept anything big enough:
    a cached reference list or an access-denied page above the floor passes
    both and would then be preferred over the FULL_CONTEXT it replaces. The
    canonical junk/paper-structure check is the one harvesting already applies
    to every rescued body, so a staged stand-in must clear the same bar.
    """
    from harvesting.content_validation import validate_content_quality

    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            head = handle.read(8192)
    except OSError as exc:
        return False, f"unreadable ({exc})"
    return validate_content_quality(head)


def is_usable_fulltext_source(path: Path) -> bool:
    """Return True when a markdown source is eligible for full-text extraction."""
    try:
        if not path.exists() or path.stat().st_size == 0:
            return False
        with path.open("r", encoding="utf-8", errors="replace") as f:
            head = f.read(8192)
        if len(head) < MIN_EXTRACTION_INPUT_SIZE:
            return False
        if is_abstract_only_fallback_text(head):
            return False
        from harvesting.content_validation import is_abstract_reference_shell

        return not is_abstract_reference_shell(head)
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


def demote_empty_full_context(
    full_context_path: Path,
    *,
    floor_bytes: int = EMPTY_FULL_CONTEXT_FLOOR_BYTES,
) -> Path:
    """Never let an empty FULL_CONTEXT win over a populated sibling rendering.

    corpus/KCNQ1/27114410 carries a 0-byte ``27114410_FULL_CONTEXT.md`` next to
    a 17.7 KB ``27114410_CLEANED.md``; any selector that prefers FULL_CONTEXT by
    name fed extraction nothing. When ``full_context_path`` is missing or below
    ``floor_bytes``, return the first usable sibling (``_CLEANED.md``, then
    ``_DATA_ZONES.md``) that also passes the harvest-time content-quality
    validator, and log a warning naming both files; otherwise return
    ``full_context_path`` unchanged. Neither file is ever modified.
    """
    name = full_context_path.name
    if not name.endswith("_FULL_CONTEXT.md"):
        return full_context_path
    try:
        size = full_context_path.stat().st_size
    except OSError:
        size = -1
    if size >= floor_bytes:
        return full_context_path
    for suffix in _FULL_CONTEXT_SIBLING_SUFFIXES:
        sibling = full_context_path.with_name(name.replace("_FULL_CONTEXT.md", suffix))
        if is_usable_fulltext_source(sibling):
            article_shaped, quality_reason = _is_article_shaped_sibling(sibling)
            if not article_shaped:
                # An access-denied page or a cached reference list clears the
                # floor; promoting it would feed extraction junk in place of
                # nothing.
                logger.warning(
                    "empty FULL_CONTEXT %s: sibling %s not staged (%s)",
                    full_context_path,
                    sibling,
                    quality_reason,
                )
                continue
            logger.warning(
                "empty FULL_CONTEXT demoted: %s is %d bytes; using sibling %s "
                "(%d bytes) instead",
                full_context_path,
                max(size, 0),
                sibling,
                sibling.stat().st_size,
            )
            return sibling
    return full_context_path


def is_reusable_fulltext_source(path: Path) -> bool:
    """Return True only when source is usable and not known body-only/incomplete."""

    if not is_usable_fulltext_source(path):
        return False
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            return not has_incomplete_supplement_surface(handle.read(8192))
    except OSError:
        return False
