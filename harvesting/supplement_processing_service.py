"""Shared supplement download and conversion processing helpers."""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, List, Optional

logger = logging.getLogger(__name__)

# Recorded on a downloaded PDF whose content identity could not be tied to the
# parent article (see :func:`pdf_supplement_identity`). Such a file stays on
# disk but is excluded from context folding — KCNQ1 31520628 folded two cited
# CDC vital-statistics PDFs into FULL_CONTEXT because nothing checked identity
# after download.
SUPPLEMENT_IDENTITY_UNVERIFIED = "supplement_identity_unverified"

# How much of a PDF's converted text stands in for its "first page" during the
# identity check. Converters emit pages in order, so the head of the markdown
# is the front matter; scanning much further would let a citing document's own
# reference list satisfy the check (same trap scholar_pdf_fallback guards).
SUPPLEMENT_IDENTITY_HEAD_CHARS = 4_000

# 'Supplement'-family wording in the material's own text.
_SUPPLEMENT_TEXT_MARKER_RE = re.compile(
    r"(?i)\b(?:"
    r"supplement(?:s|al|ary)?|supporting\s+information|appendix|annex|"
    r"e-?appendix|extended\s+data|additional\s+file|esm|moesm\d*|"
    r"(?:table|figure|fig\.?|movie|video|data(?:set)?)\s*s\d+"
    r")\b"
)

# 'Supplement'-family markers in how the object itself is minted — publisher
# supplement filenames/paths (Elsevier mmc1, Springer MOESM/ESM, /doi/suppl/,
# JATS media ids like pone.0012345.s001). Structural identity of the link, not
# a host list.
_SUPPLEMENT_PATH_MARKER_RE = re.compile(
    r"(?i)(?:"
    r"suppl|supplement|supporting|appendix|mmc\d+|moesm\d*|(?:^|[._-])esm(?:[._-]|$)|"
    r"downloadsupplement|mediaobjects"
    r")"
)
_SUPPLEMENT_STEM_SUFFIX_RE = re.compile(r"(?i)(?:^|[._-])s[a-z]?\d{1,4}$")


@dataclass
class SupplementFileResult:
    """Per-file outcome of supplement processing.

    Used by the artifacts audit log so we can answer "did this paper's
    supplement N actually convert to text?" downstream.
    """

    filename: str
    path: str
    url: str = ""
    source: str = ""  # source label from upstream fetcher (e.g. "elsevier_api")
    description: str = ""
    extension: str = ""
    downloaded: bool = False
    converted_chars: int = 0
    figures_extracted: int = 0
    nested_files: List[str] = field(default_factory=list)
    error: Optional[str] = None
    size_bytes: int = 0
    identity_unverified: bool = False


@dataclass
class SupplementProcessingResult:
    """Summary of supplement processing outputs."""

    supplement_markdown: str
    downloaded_count: int
    total_figures_extracted: int
    file_results: List[SupplementFileResult] = field(default_factory=list)


def pdf_supplement_identity(
    *,
    text_head: str,
    filename: str = "",
    source_url: str = "",
    pmid: str = "",
    doi: Optional[str] = None,
    title: Optional[str] = None,
) -> tuple[bool, str]:
    """Fail closed unless a downloaded PDF identifies as the paper's material.

    A scraped file-shaped link can bind a completely unrelated document to a
    paper (KCNQ1 31520628 captured two CDC vital-statistics PDFs its authors
    merely cited). Verified when ANY of these holds, cheapest first:

    * the filename or recorded URL is minted as supplementary material
      (``mmc1``, ``MOESM``, ``/doi/suppl/``, JATS ``.s001`` media ids, …);
    * the extractable front matter carries 'supplement'-family wording;
    * the front matter carries the article DOI, PMID, or title.

    Same fail-closed posture as ``scholar_pdf_fallback._source_identity_matches``
    (whose normalizers this reuses): absence of every marker means the file's
    identity contradicts — or at best never confirms — the parent article, so
    it must not be folded into context. The file itself always stays on disk.
    """
    from harvesting.scholar_pdf_fallback import (
        _compact_identity_text,
        _normalized_doi,
    )

    stem = Path(str(filename or "")).stem
    for candidate in (filename, source_url):
        if candidate and _SUPPLEMENT_PATH_MARKER_RE.search(str(candidate)):
            return True, "supplement marker in filename/URL"
    if stem and _SUPPLEMENT_STEM_SUFFIX_RE.search(stem):
        return True, "supplement-style media id in filename"

    head = str(text_head or "")[:SUPPLEMENT_IDENTITY_HEAD_CHARS]
    if _SUPPLEMENT_TEXT_MARKER_RE.search(head):
        return True, "supplement marker in extractable text"

    compact_head = _compact_identity_text(head)
    compact_doi = _compact_identity_text(_normalized_doi(doi))
    if compact_doi and compact_doi in compact_head:
        return True, "article DOI in extractable text"

    clean_pmid = str(pmid or "").strip()
    if re.fullmatch(r"\d{5,9}", clean_pmid) and re.search(
        rf"(?<!\d){clean_pmid}(?!\d)", head
    ):
        return True, "article PMID in extractable text"

    compact_title = _compact_identity_text(title)
    if len(compact_title) >= 16 and compact_title in compact_head:
        return True, "article title in extractable text"

    return (
        False,
        "no supplement marker, article DOI, PMID, or title in filename/URL "
        "or extractable text",
    )


def _convert_supplement(
    *,
    file_path: Path,
    converter: Any,
    extract_figures: bool,
    figures_dir: Optional[Path],
    logger: Any,
) -> tuple[str, int, List[str]]:
    """Dispatch a supplement file to the right converter.

    Returns ``(markdown, figures_extracted, nested_files)``.
    """
    ext = file_path.suffix.lower()
    nested_files: List[str] = []

    if ext in {".xlsx", ".xls"}:
        return converter.excel_to_markdown(file_path), 0, nested_files
    if ext == ".docx":
        return converter.docx_to_markdown(file_path), 0, nested_files
    if ext == ".doc":
        return converter.doc_to_markdown(file_path), 0, nested_files
    if ext == ".pdf":
        if extract_figures and figures_dir is not None:
            text, images = converter.pdf_to_markdown_with_images(
                file_path, output_dir=figures_dir
            )
            count = len(images) if images else 0
            if count and logger is not None:
                logger.info("Extracted %s figures from %s", count, file_path.name)
            return text, count, nested_files
        return converter.pdf_to_markdown(file_path), 0, nested_files
    if ext == ".csv":
        try:
            text = file_path.read_text(encoding="utf-8", errors="ignore")
            return text + "\n\n", 0, nested_files
        except Exception as exc:
            return f"[Error reading CSV file: {exc}]\n\n", 0, nested_files
    if ext == ".tsv":
        return converter.tsv_to_markdown(file_path), 0, nested_files
    if ext == ".txt":
        try:
            text = file_path.read_text(encoding="utf-8", errors="ignore")
            return text + "\n\n", 0, nested_files
        except Exception as exc:
            return f"[Error reading text file: {exc}]\n\n", 0, nested_files
    if ext in {".html", ".htm"}:
        return converter.html_supplement_to_markdown(file_path), 0, nested_files
    if ext == ".xml":
        return converter.xml_supplement_to_markdown(file_path), 0, nested_files
    if ext == ".zip":
        zip_dest = file_path.parent / file_path.stem
        kwargs: dict[str, Any] = {"dest_dir": zip_dest}
        if extract_figures and figures_dir is not None:
            kwargs["figures_dir"] = figures_dir
            kwargs["extract_images"] = True
        result = converter.extract_zip_supplement(file_path, **kwargs)
        # extract_zip_supplement may return (paths, md) or (paths, md, figures)
        if len(result) == 3:
            extracted_paths, md, zip_figs = result
            figs_count = len(zip_figs) if zip_figs else 0
        else:
            extracted_paths, md = result
            figs_count = 0
        nested_files = [str(p) for p in extracted_paths]
        return md, figs_count, nested_files

    # Unknown / binary extension — record path but no markdown content.
    return f"[File available at: {file_path}]\n\n", 0, nested_files


def process_supplement_files(
    *,
    supp_files: list[dict[str, Any]],
    supplements_dir: Path,
    pmid: str,
    converter: Any,
    download_callback: Callable[[str, Path, str, str, dict[str, Any]], bool],
    extract_figures: bool = False,
    figures_dir: Optional[Path] = None,
    logger: Any = None,
    sleep_seconds: float = 0.5,
    sleep_fn: Callable[[float], None] = time.sleep,
    doi: Optional[str] = None,
    title: Optional[str] = None,
) -> SupplementProcessingResult:
    """Download and convert supplement files into markdown.

    Returns a :class:`SupplementProcessingResult` whose ``file_results`` list
    has one entry per supplement processed (whether download succeeded or
    not). The combined ``supplement_markdown`` preserves the prior format so
    downstream consumers (LLM extraction, scout) keep working.

    ``doi``/``title`` (optional) strengthen the PDF identity gate: a PDF whose
    identity cannot be tied to the parent article (see
    :func:`pdf_supplement_identity`) is kept on disk and recorded with
    ``identity_unverified`` but its text is excluded from the combined
    markdown.
    """
    supplements_dir.mkdir(exist_ok=True)
    if extract_figures and figures_dir is not None:
        figures_dir.mkdir(exist_ok=True)

    supplement_markdown = ""
    downloaded_count = 0
    total_figures_extracted = 0
    file_results: List[SupplementFileResult] = []

    for idx, supp in enumerate(supp_files, 1):
        url = supp.get("url", "")
        filename = supp.get("name", f"supplement_{idx}")
        if not url:
            continue

        file_path = supplements_dir / filename
        print(f"    Downloading: {filename}")

        per_file = SupplementFileResult(
            filename=filename,
            path=str(file_path),
            url=url,
            source=str(supp.get("source") or ""),
            description=str(supp.get("description") or ""),
            extension=file_path.suffix.lower(),
        )

        ok = False
        try:
            ok = bool(download_callback(url, file_path, pmid, filename, supp))
        except Exception as exc:
            per_file.error = f"download_callback raised: {exc}"
            if logger is not None:
                logger.warning(
                    "Supplement download callback raised for %s: %s", filename, exc
                )

        if ok:
            per_file.downloaded = True
            downloaded_count += 1
            try:
                per_file.size_bytes = (
                    file_path.stat().st_size if file_path.exists() else 0
                )
            except OSError:
                per_file.size_bytes = 0

            md, figs_count, nested = _convert_supplement(
                file_path=file_path,
                converter=converter,
                extract_figures=extract_figures,
                figures_dir=figures_dir,
                logger=logger,
            )
            per_file.converted_chars = len(md)
            per_file.figures_extracted = figs_count
            per_file.nested_files = nested
            total_figures_extracted += figs_count

            if file_path.suffix.lower() == ".pdf":
                identity_ok, identity_reason = pdf_supplement_identity(
                    text_head=md,
                    filename=filename,
                    source_url=url,
                    pmid=pmid,
                    doi=doi,
                    title=title,
                )
                if not identity_ok:
                    per_file.identity_unverified = True
                    per_file.error = SUPPLEMENT_IDENTITY_UNVERIFIED
                    # The ``logger`` parameter shadows the module logger and
                    # may be None; fall back so the exclusion is never silent.
                    (logger or logging.getLogger(__name__)).warning(
                        "%s: excluding %s from PMID %s context (%s); "
                        "file kept on disk at %s",
                        SUPPLEMENT_IDENTITY_UNVERIFIED,
                        filename,
                        pmid,
                        identity_reason,
                        file_path,
                    )
                    file_results.append(per_file)
                    sleep_fn(sleep_seconds)
                    continue

            supplement_markdown += f"\n\n# SUPPLEMENTAL FILE {idx}: {filename}\n\n"
            supplement_markdown += md
        else:
            per_file.error = per_file.error or "download_failed"

        file_results.append(per_file)
        sleep_fn(sleep_seconds)

    return SupplementProcessingResult(
        supplement_markdown=supplement_markdown,
        downloaded_count=downloaded_count,
        total_figures_extracted=total_figures_extracted,
        file_results=file_results,
    )
