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

import json
import logging
import re
from pathlib import Path
from typing import Any, Optional

from harvesting.supplement_processing_service import (
    SUPPLEMENT_IDENTITY_UNVERIFIED,
    _convert_supplement,
    pdf_supplement_identity,
)

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
_NESTED_FILE_HEADING_RE = re.compile(
    r"(?mi)^#{1,6}\s+Nested file:\s*(?P<label>.+?)\s*$"
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


def _archive_nested_labels(text: str) -> set[str]:
    """Basenames already expanded inside a harvest-time archive block."""
    return {
        Path(match.group("label").strip()).name.casefold()
        for match in _NESTED_FILE_HEADING_RE.finditer(text)
    }


def _exclude_supplement_blocks(markdown: str, labels: set[str]) -> str:
    """Remove converted-file blocks already present via an expanded archive."""
    if not labels:
        return markdown
    matches = list(_SUPPLEMENT_HEADING_RE.finditer(markdown))
    if not matches:
        return markdown
    kept: list[str] = []
    for idx, match in enumerate(matches):
        end = matches[idx + 1].start() if idx + 1 < len(matches) else len(markdown)
        label = Path(match.group("label").strip()).name.casefold()
        if label not in labels:
            kept.append(markdown[match.start() : end].strip())
    return "\n\n".join(kept).strip()


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


_CONVERSION_PLACEHOLDER_RE = re.compile(
    r"^\[(?:Invalid|Error (?:reading|converting)|PDF file available|"
    r"Legacy \.doc file available)",
    re.IGNORECASE,
)


def _is_conversion_placeholder(markdown: str) -> bool:
    """Return true when a converter emitted an error marker, not source text."""

    return bool(_CONVERSION_PLACEHOLDER_RE.match(str(markdown or "").strip()))


def _paper_identity_context(pmid: str, harvest_dir: Path) -> dict[str, Any]:
    """Best-effort article identity for the PDF supplement gate.

    Reads the ``{pmid}_artifacts.json`` sidecar when present for the DOI and
    the per-filename source URLs the harvester recorded. Everything is
    optional — the gate can still verify via the PMID or supplement-family
    markers when the manifest is missing or predates these fields.
    """
    identity: dict[str, Any] = {"pmid": str(pmid), "doi": None, "urls": {}}
    manifest_path = harvest_dir / f"{pmid}_artifacts.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return identity
    if isinstance(manifest, dict):
        identity["doi"] = manifest.get("doi") or None
        for entry in manifest.get("supplements") or []:
            if not isinstance(entry, dict):
                continue
            filename = str(entry.get("filename") or "").strip()
            url = str(entry.get("url") or "").strip()
            if filename and url:
                identity["urls"][filename.casefold()] = url
    return identity


def _build_supplement_markdown_result(
    supplements_dir: Path,
    *,
    converter: Any = None,
    logger_obj: Any = None,
    identity: Optional[dict[str, Any]] = None,
) -> tuple[str, int, int, set[str]]:
    """Convert every convertible file in ``supplements_dir`` to combined markdown.

    Returns ``(markdown, files_converted, conversion_failures,
    omitted_labels)`` where ``omitted_labels`` are file labels that are safe to
    drop from an existing fold block: conversion placeholders plus, when an
    ``identity`` context is supplied, PDFs quarantined by
    :func:`~harvesting.supplement_processing_service.pdf_supplement_identity`
    (kept on disk, never folded). ``converter``
    defaults to a fresh
    :class:`~harvesting.format_converters.FormatConverter` (only constructed when
    a file needs it; plain ``.csv``/``.txt`` are read directly).
    """
    if not supplements_dir.is_dir():
        return "", 0, 0, set()
    # Recursive walk (not top-level iterdir): convertible files extracted from a
    # ``.zip`` supplement land in a subdirectory and would otherwise be missed.
    # The ``.zip`` itself stays excluded (not a convertible suffix), so we fold
    # the extracted files exactly once; the sentinel-delimited rebuild keeps this
    # idempotent. Skip macOS zip cruft and hidden/AppleDouble files.
    files = _convertible_files(supplements_dir)
    if not files:
        return "", 0, 0, set()
    if converter is None:
        from harvesting.format_converters import FormatConverter

        converter = FormatConverter()

    parts: list[str] = []
    converted = 0
    failures = 0
    omitted_labels: set[str] = set()
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
        if _is_conversion_placeholder(md):
            logger.warning("supplement convert produced placeholder for %s", rel)
            failures += 1
            omitted_labels.add(Path(rel).name.casefold())
            continue
        if file_path.suffix.lower() == ".pdf" and identity is not None:
            identity_ok, identity_reason = pdf_supplement_identity(
                text_head=md,
                filename=file_path.name,
                source_url=str(
                    (identity.get("urls") or {}).get(file_path.name.casefold(), "")
                ),
                pmid=str(identity.get("pmid") or ""),
                doi=identity.get("doi"),
                title=identity.get("title"),
            )
            if not identity_ok:
                logger.warning(
                    "%s: excluding %s from fold (%s); file kept on disk",
                    SUPPLEMENT_IDENTITY_UNVERIFIED,
                    rel,
                    identity_reason,
                )
                omitted_labels.add(Path(rel).name.casefold())
                continue
        if md and md.strip():
            parts.append(f"\n\n# SUPPLEMENTAL FILE {idx}: {rel}\n\n{md}")
            converted += 1
    return "".join(parts).strip(), converted, failures, omitted_labels


def build_supplement_markdown(
    supplements_dir: Path,
    *,
    converter: Any = None,
    logger_obj: Any = None,
) -> tuple[str, int]:
    """Convert available supplement text, preserving the public two-value API."""
    markdown, converted, _failures, _omitted_labels = _build_supplement_markdown_result(
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
    md, converted, failures, omitted_labels = _build_supplement_markdown_result(
        supp_dir,
        converter=converter,
        identity=_paper_identity_context(pmid, harvest_dir),
    )
    original = full_context.read_text(encoding="utf-8", errors="replace")
    if converted == 0 or not md:
        # No foldable text. Still rebuild when a prior fold carries text from
        # files the identity gate has since quarantined (KCNQ1 31520628's only
        # on-disk supplements are the two mis-bound CDC reports), so that
        # content does not survive a refresh; otherwise leave the file alone.
        if not omitted_labels or not (_supplement_labels(original) & omitted_labels):
            return None
        md = ""
    if failures:
        logger.warning(
            "PMID %s had %d supplement conversion failure(s); proceeding with "
            "successfully converted files",
            pmid,
            failures,
        )

    replacement_labels = _supplement_labels(md)
    begin = original.find(FOLD_BEGIN)
    if begin != -1:
        end = original.find(FOLD_END, begin)
        old_block = (
            original[begin:] if end == -1 else original[begin : end + len(FOLD_END)]
        )
        # A prior pipeline version folded converter error markers — and, before
        # the identity gate, mis-bound PDFs (KCNQ1 31520628's CDC reports) — as
        # if they were source. Those labels are safe to discard; otherwise they
        # can permanently block a newly recovered valid supplement from
        # refolding.
        missing_old_labels = (
            _supplement_labels(old_block) - replacement_labels - omitted_labels
        )
        if missing_old_labels:
            logger.warning(
                "refusing supplement re-fold for PMID %s: replacement omits %s",
                pmid,
                ", ".join(sorted(missing_old_labels)),
            )
            return None
    base = _strip_existing_supplement_blocks(original, md)
    # A harvest-time ZIP conversion writes each member under a ``Nested file``
    # heading inside the archive's SUPPLEMENTAL FILE block. The on-disk fold
    # later sees those extracted members as ordinary files. Appending them again
    # doubled large tables (and, for PMID 30702160, expanded a 2.8 MB source to
    # 43 MB). Keep only members not already represented in the archive block.
    md = _exclude_supplement_blocks(md, _archive_nested_labels(base))

    if md:
        folded = (
            base.rstrip()
            + f"\n\n{FOLD_BEGIN}\n\n"
            + "# FOLDED SUPPLEMENTS (re-extraction aid)\n"
            + md
            + f"\n\n{FOLD_END}\n"
        )
    else:
        # This can intentionally remove a redundant sentinel block written by
        # an older fold. The harvest-time archive copy remains in ``base``.
        folded = base.rstrip() + "\n"
    if folded == original:
        return None
    backup = full_context.parent / (full_context.name + ".pre_fold_bak")
    if not backup.exists():
        backup.write_text(original, encoding="utf-8")
    full_context.write_text(folded, encoding="utf-8")
    logger.info(
        "folded %d supplement file(s) into %s (+%d chars over base)",
        len(_supplement_labels(md)),
        full_context.name,
        len(folded) - len(base),
    )
    return full_context
