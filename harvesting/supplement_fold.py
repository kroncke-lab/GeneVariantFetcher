"""Fold already-on-disk ``{pmid}_supplements/`` files into a source FULL_CONTEXT.

Supplement files are downloaded and converted to markdown only at initial
harvest time (``harvesting/orchestrator.py``), folded once into
``{pmid}_FULL_CONTEXT.md``. The re-extraction / replay path
(``scripts/refresh_run_db.py``) then discovers sources only by globbing
``*_DATA_ZONES.md`` / ``*_CLEANED.md`` / ``*_FULL_CONTEXT.md`` and never re-reads
``{pmid}_supplements/``. So a paper whose supplement tables were downloaded but
not folded into FULL_CONTEXT (stale binding, a thin harvest-time fold, or a
side-directory recovery run) silently loses those variants on re-extraction.

This module re-folds the on-disk supplement files into the PMID's FULL_CONTEXT
so the standard discovery path sees them, with three safety properties:

* **non-destructive** — the original FULL_CONTEXT is backed up once to
  ``{pmid}_FULL_CONTEXT.md.pre_fold_bak`` before the first fold;
* **idempotent** — the folded text is delimited by a sentinel and regenerated
  on each run, so re-folding never double-appends;
* **gene-safe** — it only assembles text. Gene scoping happens at parse time
  (the markdown and fixed-width table parsers scope rows by caption gene), and
  the downstream ``refresh_run_db`` explosion gate guards against a garbage
  blow-up reaching the DB.

``.zip`` supplements are intentionally skipped: ``_convert_supplement`` extracts
them onto disk inside the supplements dir, which a second fold pass would then
re-discover and double-count, breaking idempotency. Zips are handled at harvest
time.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Optional

from harvesting.supplement_processing_service import _convert_supplement

logger = logging.getLogger(__name__)

FOLD_BEGIN = "<!-- GVF_FOLDED_SUPPLEMENTS_BEGIN -->"
FOLD_END = "<!-- GVF_FOLDED_SUPPLEMENTS_END -->"

# Harvest-time assembly predates the sentinel block and appended converted
# supplements directly with this generated heading.  When a later supplement
# refresh introduces the sentinel form, strip that legacy tail first so the
# same files are not represented twice.
_LEGACY_FOLD_BEGIN_RE = re.compile(r"(?m)^# SUPPLEMENTAL FILE 1:\s*.*$")
_SUPPLEMENT_HEADING_RE = re.compile(
    r"(?m)^# SUPPLEMENTAL FILE \d+:\s*(?P<label>.+?)\s*$"
)

# Supplement extensions we re-fold. ``.zip`` is deliberately excluded (see the
# module docstring); image/binary types are skipped (no usable markdown).
_CONVERTIBLE_SUFFIXES = {
    ".xlsx",
    ".xls",
    ".docx",
    ".doc",
    ".pdf",
    ".csv",
    ".tsv",
    ".txt",
    ".html",
    ".htm",
    ".xml",
}


def _strip_folded_block(text: str) -> str:
    """Remove any previously-folded supplement block so re-folds don't stack."""
    begin = text.find(FOLD_BEGIN)
    if begin == -1:
        return text
    end = text.find(FOLD_END, begin)
    if end == -1:
        # Truncated/corrupt end marker: drop everything from the begin marker.
        return text[:begin].rstrip() + "\n"
    end += len(FOLD_END)
    return (text[:begin].rstrip() + "\n" + text[end:].lstrip()).rstrip() + "\n"


def _supplement_labels(text: str) -> set[str]:
    """Normalized file labels represented by generated supplement headings."""
    return {
        Path(match.group("label").strip()).name.casefold()
        for match in _SUPPLEMENT_HEADING_RE.finditer(text)
    }


def _strip_existing_supplement_blocks(text: str, replacement_markdown: str) -> str:
    """Remove replaceable folds while retaining richer legacy-only content."""
    base = _strip_folded_block(text)
    legacy = _LEGACY_FOLD_BEGIN_RE.search(base)
    if legacy is not None:
        # A heading match alone does not prove that the rest of the paper is a
        # generated tail. Strip only an exact copy of the replacement block;
        # otherwise retain the legacy text and append the sentinel block.
        legacy_tail = base[legacy.start() :].strip()
        if legacy_tail == replacement_markdown.strip():
            base = base[: legacy.start()].rstrip() + "\n"
    return base


def _convertible_files(supplements_dir: Path) -> list[Path]:
    if not supplements_dir.is_dir():
        return []
    return sorted(
        (
            p
            for p in supplements_dir.rglob("*")
            if p.is_file()
            and p.suffix.lower() in _CONVERTIBLE_SUFFIXES
            and "__MACOSX" not in p.parts
            and not any(
                part.startswith(".") for part in p.relative_to(supplements_dir).parts
            )
        ),
        key=lambda p: p.relative_to(supplements_dir).as_posix(),
    )


def _build_supplement_markdown_result(
    supplements_dir: Path,
    *,
    converter: Any = None,
    logger_obj: Any = None,
) -> tuple[str, int, int]:
    """Convert every convertible file in ``supplements_dir`` to combined markdown.

    Returns ``(markdown, files_converted, conversion_failures)``. ``converter``
    defaults to a fresh
    :class:`~harvesting.format_converters.FormatConverter` (only constructed when
    a file needs it; plain ``.csv``/``.txt`` are read directly).
    """
    if not supplements_dir.is_dir():
        return "", 0, 0
    # Recursive walk (not top-level iterdir): convertible files extracted from a
    # ``.zip`` supplement land in a subdirectory and would otherwise be missed.
    # The ``.zip`` itself stays excluded (not a convertible suffix), so we fold
    # the extracted files exactly once; the sentinel-delimited rebuild keeps this
    # idempotent. Skip macOS zip cruft and hidden/AppleDouble files.
    files = _convertible_files(supplements_dir)
    if not files:
        return "", 0, 0
    if converter is None:
        from harvesting.format_converters import FormatConverter

        converter = FormatConverter()

    parts: list[str] = []
    converted = 0
    failures = 0
    for idx, file_path in enumerate(files, 1):
        # Label by path relative to the supplements dir: a top-level file shows
        # just its name; a file extracted from a zip shows ``subdir/name`` so the
        # provenance (and the nested-zip recovery) is visible in FULL_CONTEXT.
        rel = file_path.relative_to(supplements_dir).as_posix()
        try:
            md, _figs, _nested = _convert_supplement(
                file_path=file_path,
                converter=converter,
                extract_figures=False,
                figures_dir=None,
                logger=logger_obj,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("supplement convert failed for %s: %s", rel, exc)
            failures += 1
            continue
        if md and md.strip():
            parts.append(f"\n\n# SUPPLEMENTAL FILE {idx}: {rel}\n\n{md}")
            converted += 1
    return "".join(parts).strip(), converted, failures


def build_supplement_markdown(
    supplements_dir: Path,
    *,
    converter: Any = None,
    logger_obj: Any = None,
) -> tuple[str, int]:
    """Convert available supplement text, preserving the public two-value API."""
    markdown, converted, _failures = _build_supplement_markdown_result(
        supplements_dir,
        converter=converter,
        logger_obj=logger_obj,
    )
    return markdown, converted


def fold_supplements_into_full_context(
    pmid: str,
    harvest_dir: Path,
    *,
    supplements_dir: Optional[Path] = None,
    converter: Any = None,
) -> Optional[Path]:
    """Fold on-disk supplements for ``pmid`` into its FULL_CONTEXT (idempotent).

    Returns the FULL_CONTEXT path when it was (re)written, otherwise ``None``
    (no FULL_CONTEXT present, or no convertible supplements found).
    """
    full_context = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
    if not full_context.is_file():
        return None
    supp_dir = supplements_dir or (harvest_dir / f"{pmid}_supplements")
    md, converted, failures = _build_supplement_markdown_result(
        supp_dir, converter=converter
    )
    if converted == 0 or not md:
        return None
    if failures:
        logger.warning(
            "refusing partial supplement re-fold for PMID %s: %d conversion failure(s)",
            pmid,
            failures,
        )
        return None

    original = full_context.read_text(encoding="utf-8", errors="replace")
    replacement_labels = _supplement_labels(md)
    begin = original.find(FOLD_BEGIN)
    if begin != -1:
        end = original.find(FOLD_END, begin)
        old_block = (
            original[begin:] if end == -1 else original[begin : end + len(FOLD_END)]
        )
        missing_old_labels = _supplement_labels(old_block) - replacement_labels
        if missing_old_labels:
            logger.warning(
                "refusing supplement re-fold for PMID %s: replacement omits %s",
                pmid,
                ", ".join(sorted(missing_old_labels)),
            )
            return None
    base = _strip_existing_supplement_blocks(original, md)

    folded = (
        base.rstrip()
        + f"\n\n{FOLD_BEGIN}\n\n"
        + "# FOLDED SUPPLEMENTS (re-extraction aid)\n"
        + md
        + f"\n\n{FOLD_END}\n"
    )
    if folded == original:
        return None
    backup = full_context.parent / (full_context.name + ".pre_fold_bak")
    if not backup.exists():
        backup.write_text(original, encoding="utf-8")
    full_context.write_text(folded, encoding="utf-8")
    logger.info(
        "folded %d supplement file(s) into %s (+%d chars over base)",
        converted,
        full_context.name,
        len(folded) - len(base),
    )
    return full_context
