"""Shared supplement download and conversion processing helpers."""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, List, Optional


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


@dataclass
class SupplementProcessingResult:
    """Summary of supplement processing outputs."""

    supplement_markdown: str
    downloaded_count: int
    total_figures_extracted: int
    file_results: List[SupplementFileResult] = field(default_factory=list)


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
) -> SupplementProcessingResult:
    """Download and convert supplement files into markdown.

    Returns a :class:`SupplementProcessingResult` whose ``file_results`` list
    has one entry per supplement processed (whether download succeeded or
    not). The combined ``supplement_markdown`` preserves the prior format so
    downstream consumers (LLM extraction, scout) keep working.
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

            supplement_markdown += f"\n\n# SUPPLEMENTAL FILE {idx}: {filename}\n\n"
            md, figs_count, nested = _convert_supplement(
                file_path=file_path,
                converter=converter,
                extract_figures=extract_figures,
                figures_dir=figures_dir,
                logger=logger,
            )
            supplement_markdown += md
            per_file.converted_chars = len(md)
            per_file.figures_extracted = figs_count
            per_file.nested_files = nested
            total_figures_extracted += figs_count
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
