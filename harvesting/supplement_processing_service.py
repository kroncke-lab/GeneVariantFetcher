"""Shared supplement download and conversion processing helpers."""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, List, Optional

from config.constants import MIN_EXTRACTION_INPUT_SIZE

logger = logging.getLogger(__name__)

# Recorded on a downloaded PDF whose content identity CONTRADICTS the parent
# article (see :func:`pdf_supplement_verdict`). Such a file stays on disk but
# is excluded from context folding — KCNQ1 31520628 folded two cited CDC
# vital-statistics PDFs into FULL_CONTEXT because nothing checked identity
# after download.
SUPPLEMENT_IDENTITY_UNVERIFIED = "supplement_identity_unverified"

# Three-state identity outcome. UNKNOWN is not QUARANTINE: a scanned mutation
# table named ``download.pdf`` carries no marker at all, so treating absence of
# evidence as evidence of a foreign document dropped real supplements.
SUPPLEMENT_IDENTITY_VERIFIED = "verified"
SUPPLEMENT_IDENTITY_QUARANTINE = "quarantine"
SUPPLEMENT_IDENTITY_UNKNOWN = "unknown"

# Where a scraped supplement link sat in the article markup, carried on the
# scraped entry's ``source`` label (the artifacts record's provenance field) so
# both the download-time and the fold-time gate can read it back.
# ``harvesting.supplement_scraper`` writes these.
SUPPLEMENT_PROVENANCE_REFERENCE_SECTION = "scraper_reference_section"
SUPPLEMENT_PROVENANCE_SUPPLEMENT_CONTAINER = "scraper_supplementary_container"

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

# A DOI printed in the front matter. Reaching the contradiction checks means
# the article's own DOI is not in the head, so any DOI here is a foreign one.
_DOI_TOKEN_RE = re.compile(r"(?i)\b10\.\d{4,9}/[^\s\"'<>,;)\]]+")

# A serial's masthead: a standing title above a volume/number/date line. This
# is what identifies KCNQ1 31520628's two cited CDC vital-statistics reports
# ("National Vital / Statistics Reports / Volume 60, Number 7 June 20, 2012"),
# which carry no DOI anywhere in their text and no reference-section
# provenance in a manifest written before provenance was recorded.
_SERIAL_MASTHEAD_RE = re.compile(
    r"(?im)^[^\n]{0,120}?\bvol(?:ume)?\.?\s*\d{1,4}\b[^\n]{0,60}?"
    r"\b(?:no|nos|number|issue)s?\.?\s*\d{1,4}\b[^\n]{0,60}?\b(?:19|20)\d{2}\b"
)
#: A serial masthead is repeated as a running header on every page, so where
#: its FIRST occurrence lands depends on how the front page happens to
#: extract. Scanning only the opening block caught KCNQ1 31520628's
#: ``nvsr60_07.pdf`` (header at offset 992) and missed its twin
#: ``nvsr60_08.pdf`` (offset 4,006) — the same publication, the same defect.
#: The window covers several pages so the signal is stable.
_SERIAL_MASTHEAD_SCAN_CHARS = 24_000

#: Within this opening block the standing title may sit on preceding lines;
#: past it, only an inline running header counts.
_SERIAL_MASTHEAD_FRONT_PAGE_CHARS = 1_200


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


def _has_serial_masthead(head: str) -> bool:
    """True when a PDF's first page opens like a numbered serial issue.

    Requires the standing title as well as the volume/number/date line — a
    bare ``vol 3 no 2 2011`` string inside a data table is not a masthead.
    """
    window = str(head or "")[:_SERIAL_MASTHEAD_SCAN_CHARS]
    match = _SERIAL_MASTHEAD_RE.search(window)
    if not match:
        return False
    # On a front page the standing title sits on the lines ABOVE the
    # volume/date line, so everything before the match is fair evidence.
    # Deeper in the document only a running header should count, and those
    # repeat the serial name inline — requiring it there stops a bare
    # "vol. 60, no. 7, 2012" inside a citation from passing as a masthead.
    if match.start() <= _SERIAL_MASTHEAD_FRONT_PAGE_CHARS:
        lead = window[: match.end()]
    else:
        line_start = window.rfind("\n", 0, match.start()) + 1
        lead = window[line_start : match.end()]
    volume = re.search(r"(?i)\bvol", lead)
    banner = lead[: volume.start()] if volume else ""
    return len(re.findall(r"[A-Za-z]{3,}", banner)) >= 2


def pdf_supplement_verdict(
    *,
    text_head: str,
    filename: str = "",
    source_url: str = "",
    pmid: str = "",
    doi: Optional[str] = None,
    title: Optional[str] = None,
    provenance: str = "",
) -> tuple[str, str]:
    """Decide whether a downloaded PDF may be folded into the paper's context.

    Returns ``(verdict, reason)`` with verdict one of
    ``SUPPLEMENT_IDENTITY_VERIFIED`` / ``_QUARANTINE`` / ``_UNKNOWN``.
    The file itself always stays on disk.

    Verified when ANY of these confirms the parent article, cheapest first:

    * the filename or recorded URL is minted as supplementary material
      (``mmc1``, ``MOESM``, ``/doi/suppl/``, JATS ``.s001`` media ids, …);
    * the link was harvested from a supplementary-material container;
    * the extractable front matter carries 'supplement'-family wording;
    * the front matter carries the article DOI, PMID, or title.

    Only a CONTRADICTION quarantines: reference-section provenance, a foreign
    DOI in the front matter, or a serial-publication masthead (KCNQ1
    31520628's two cited CDC vital-statistics PDFs). Anything else that merely
    fails to confirm is UNKNOWN and gets folded with an ``identity_unverified``
    flag — a genuine scanned mutation table named ``download.pdf`` carries no
    marker at all, so quarantining on absence dropped real supplements.
    """
    from harvesting.scholar_pdf_fallback import (
        _compact_identity_text,
        _normalized_doi,
    )

    stem = Path(str(filename or "")).stem
    for candidate in (filename, source_url):
        if candidate and _SUPPLEMENT_PATH_MARKER_RE.search(str(candidate)):
            return SUPPLEMENT_IDENTITY_VERIFIED, "supplement marker in filename/URL"
    if stem and _SUPPLEMENT_STEM_SUFFIX_RE.search(stem):
        return SUPPLEMENT_IDENTITY_VERIFIED, "supplement-style media id in filename"

    provenance_label = str(provenance or "")
    if SUPPLEMENT_PROVENANCE_SUPPLEMENT_CONTAINER in provenance_label:
        return (
            SUPPLEMENT_IDENTITY_VERIFIED,
            "link harvested from a supplementary-material container",
        )

    head = str(text_head or "")[:SUPPLEMENT_IDENTITY_HEAD_CHARS]
    if _SUPPLEMENT_TEXT_MARKER_RE.search(head):
        return SUPPLEMENT_IDENTITY_VERIFIED, "supplement marker in extractable text"

    compact_head = _compact_identity_text(head)
    compact_doi = _compact_identity_text(_normalized_doi(doi))
    if compact_doi and compact_doi in compact_head:
        return SUPPLEMENT_IDENTITY_VERIFIED, "article DOI in extractable text"

    clean_pmid = str(pmid or "").strip()
    if re.fullmatch(r"\d{5,9}", clean_pmid) and re.search(
        rf"(?<!\d){clean_pmid}(?!\d)", head
    ):
        return SUPPLEMENT_IDENTITY_VERIFIED, "article PMID in extractable text"

    compact_title = _compact_identity_text(title)
    if len(compact_title) >= 16 and compact_title in compact_head:
        return SUPPLEMENT_IDENTITY_VERIFIED, "article title in extractable text"

    if SUPPLEMENT_PROVENANCE_REFERENCE_SECTION in provenance_label:
        return (
            SUPPLEMENT_IDENTITY_QUARANTINE,
            "link harvested from reference/bibliography markup; no supplement "
            "marker, article DOI, PMID, or title",
        )
    head_dois = [token for token in _DOI_TOKEN_RE.findall(head) if token]
    # A token that prefixes the article DOI (or is prefixed by it) is the same
    # identifier rendered short — KCNQ1 18174212's own PDF prints its DOI
    # across the head boundary as "10.1113/jphysi". Not a contradiction.
    related = any(
        compact.startswith(compact_doi) or compact_doi.startswith(compact)
        for compact in (
            _compact_identity_text(_normalized_doi(token)) for token in head_dois
        )
        if compact
    )
    if compact_doi and head_dois and not related:
        return (
            SUPPLEMENT_IDENTITY_QUARANTINE,
            f"front matter carries a foreign DOI ({head_dois[0]}), not the "
            "article's; no supplement marker, article DOI, PMID, or title",
        )
    # A serial masthead is a WEAK signal and deliberately does not quarantine.
    # It cannot be separated from a paper's own article PDF using the anchors
    # this path actually has: LMNA 17334235's own full text opens
    # "EXPERIMENTAL and MOLECULAR MEDICINE, Vol. 39, No. 1, ... February 2007"
    # and contains neither its DOI nor its PMID anywhere, while the artifacts
    # record carries no title — so condemning on "masthead and no witness"
    # drops the paper's own text. The causal signal for KCNQ1 31520628's cited
    # CDC reports is reference-section provenance, checked above; a masthead
    # only annotates the fold so QC can see it.
    wide = str(text_head or "")[:_SERIAL_MASTHEAD_SCAN_CHARS]
    if _has_serial_masthead(wide):
        compact_wide = _compact_identity_text(wide)
        compact_title = _compact_identity_text(title)
        witness = (
            (compact_doi and compact_doi in compact_wide)
            or (
                re.fullmatch(r"\d{5,9}", clean_pmid)
                and re.search(rf"(?<!\d){clean_pmid}(?!\d)", wide)
            )
            or (len(compact_title) >= 24 and compact_title in compact_wide)
        )
        if not witness:
            return (
                SUPPLEMENT_IDENTITY_UNKNOWN,
                "reads like a numbered serial issue and names no supplement "
                "marker, article DOI, PMID, or title; folded unverified",
            )

    if len(head.strip()) < MIN_EXTRACTION_INPUT_SIZE:
        return (
            SUPPLEMENT_IDENTITY_UNKNOWN,
            f"too little extractable text to identify ({len(head.strip())} "
            "chars); folded unverified",
        )
    return (
        SUPPLEMENT_IDENTITY_UNKNOWN,
        "no supplement marker, article DOI, PMID, or title in filename/URL "
        "or extractable text, and nothing contradicting the article; folded "
        "unverified",
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

    ``doi``/``title`` (optional) strengthen the PDF identity gate (see
    :func:`pdf_supplement_verdict`). A PDF that CONTRADICTS the parent article
    is kept on disk and excluded from the combined markdown; one that merely
    fails to confirm is folded and recorded with ``identity_unverified``.
    Reference-section provenance rides in on each entry's ``source`` label.
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
                verdict, identity_reason = pdf_supplement_verdict(
                    text_head=md,
                    filename=filename,
                    source_url=url,
                    pmid=pmid,
                    doi=doi,
                    title=title,
                    provenance=per_file.source,
                )
                if verdict != SUPPLEMENT_IDENTITY_VERIFIED:
                    per_file.identity_unverified = True
                # The ``logger`` parameter shadows the module logger and may be
                # None; fall back so neither outcome is silent.
                effective_logger = logger or logging.getLogger(__name__)
                if verdict == SUPPLEMENT_IDENTITY_QUARANTINE:
                    per_file.error = SUPPLEMENT_IDENTITY_UNVERIFIED
                    effective_logger.warning(
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
                if verdict == SUPPLEMENT_IDENTITY_UNKNOWN:
                    effective_logger.warning(
                        "%s: folding %s into PMID %s context unverified (%s)",
                        SUPPLEMENT_IDENTITY_UNVERIFIED,
                        filename,
                        pmid,
                        identity_reason,
                    )

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
