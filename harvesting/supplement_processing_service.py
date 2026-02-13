"""Shared supplement download and conversion processing helpers."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional


@dataclass
class SupplementProcessingResult:
    """Summary of supplement processing outputs."""

    supplement_markdown: str
    downloaded_count: int
    total_figures_extracted: int


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
    """Download and convert supplement files into markdown."""
    supplements_dir.mkdir(exist_ok=True)
    if extract_figures and figures_dir is not None:
        figures_dir.mkdir(exist_ok=True)

    supplement_markdown = ""
    downloaded_count = 0
    total_figures_extracted = 0

    for idx, supp in enumerate(supp_files, 1):
        url = supp.get("url", "")
        filename = supp.get("name", f"supplement_{idx}")
        if not url:
            continue

        file_path = supplements_dir / filename
        print(f"    Downloading: {filename}")

        if download_callback(url, file_path, pmid, filename, supp):
            downloaded_count += 1
            ext = file_path.suffix.lower()
            supplement_markdown += f"\n\n# SUPPLEMENTAL FILE {idx}: {filename}\n\n"

            if ext in [".xlsx", ".xls"]:
                supplement_markdown += converter.excel_to_markdown(file_path)
            elif ext == ".docx":
                supplement_markdown += converter.docx_to_markdown(file_path)
            elif ext == ".doc":
                supplement_markdown += converter.doc_to_markdown(file_path)
            elif ext == ".pdf":
                if extract_figures and figures_dir:
                    text, images = converter.pdf_to_markdown_with_images(
                        file_path, output_dir=figures_dir
                    )
                    supplement_markdown += text
                    if images:
                        total_figures_extracted += len(images)
                        if logger is not None:
                            logger.info(
                                "Extracted %s figures from %s", len(images), filename
                            )
                else:
                    supplement_markdown += converter.pdf_to_markdown(file_path)
            elif ext in [".txt", ".csv"]:
                try:
                    text = file_path.read_text(encoding="utf-8", errors="ignore")
                    supplement_markdown += text + "\n\n"
                except Exception as exc:
                    supplement_markdown += f"[Error reading text file: {exc}]\n\n"
            else:
                supplement_markdown += f"[File available at: {file_path}]\n\n"

        sleep_fn(sleep_seconds)

    return SupplementProcessingResult(
        supplement_markdown=supplement_markdown,
        downloaded_count=downloaded_count,
        total_figures_extracted=total_figures_extracted,
    )
