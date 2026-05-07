"""Tests for the additive supplement-processing surface.

Covers the richer per-file metadata (``SupplementFileResult``) and the new
format dispatch (HTML / TSV / XML / ZIP) added during the figure +
supplement harvesting overhaul. The pre-existing test_supplement_processing_service.py
covers the original happy-path dispatch; these tests pin the new behavior.
"""

from __future__ import annotations

import zipfile
from pathlib import Path

import pytest

from harvesting.supplement_processing_service import (
    SupplementFileResult,
    SupplementProcessingResult,
    process_supplement_files,
)


class _StubConverter:
    """Records which converter method was called for each file."""

    def __init__(self):
        self.calls: list[str] = []

    def excel_to_markdown(self, p: Path) -> str:
        self.calls.append(f"excel:{p.name}")
        return f"# excel {p.name}\n"

    def docx_to_markdown(self, p: Path) -> str:
        self.calls.append(f"docx:{p.name}")
        return f"# docx {p.name}\n"

    def doc_to_markdown(self, p: Path) -> str:
        self.calls.append(f"doc:{p.name}")
        return f"# doc {p.name}\n"

    def pdf_to_markdown(self, p: Path) -> str:
        self.calls.append(f"pdf:{p.name}")
        return f"# pdf {p.name}\n"

    def pdf_to_markdown_with_images(self, p: Path, output_dir: Path):
        self.calls.append(f"pdf+img:{p.name}")
        return f"# pdf+img {p.name}\n", []

    def tsv_to_markdown(self, p: Path) -> str:
        self.calls.append(f"tsv:{p.name}")
        return f"# tsv {p.name}\n"

    def html_supplement_to_markdown(self, p: Path) -> str:
        self.calls.append(f"html:{p.name}")
        return f"# html {p.name}\n"

    def xml_supplement_to_markdown(self, p: Path) -> str:
        self.calls.append(f"xml:{p.name}")
        return f"# xml {p.name}\n"

    def extract_zip_supplement(self, p: Path, dest_dir: Path):
        self.calls.append(f"zip:{p.name}")
        dest_dir.mkdir(parents=True, exist_ok=True)
        nested = dest_dir / "inner.txt"
        nested.write_text("inner")
        return [nested], "# zip body\n"


def _touch_callback(
    url: str, file_path: Path, pmid: str, filename: str, supp: dict
) -> bool:
    file_path.write_text(f"data:{filename}")
    return True


def test_dispatch_routes_new_extensions(tmp_path: Path):
    converter = _StubConverter()
    supplements_dir = tmp_path / "supps"
    files = [
        {"url": "http://x/a.xlsx", "name": "a.xlsx", "source": "elsevier_api"},
        {"url": "http://x/b.html", "name": "b.html"},
        {"url": "http://x/c.tsv", "name": "c.tsv"},
        {"url": "http://x/d.xml", "name": "d.xml"},
    ]

    result = process_supplement_files(
        supp_files=files,
        supplements_dir=supplements_dir,
        pmid="1",
        converter=converter,
        download_callback=_touch_callback,
        sleep_fn=lambda _s: None,
    )

    assert isinstance(result, SupplementProcessingResult)
    assert result.downloaded_count == 4
    assert len(result.file_results) == 4
    extensions = sorted(r.extension for r in result.file_results)
    assert extensions == [".html", ".tsv", ".xlsx", ".xml"]
    assert all(r.converted_chars > 0 for r in result.file_results)
    # Source metadata threaded through.
    by_name = {r.filename: r for r in result.file_results}
    assert by_name["a.xlsx"].source == "elsevier_api"
    assert "html:b.html" in converter.calls
    assert "tsv:c.tsv" in converter.calls
    assert "xml:d.xml" in converter.calls


def test_zip_supplement_records_nested_files(tmp_path: Path):
    converter = _StubConverter()
    supplements_dir = tmp_path / "supps"
    files = [{"url": "http://x/bundle.zip", "name": "bundle.zip"}]

    result = process_supplement_files(
        supp_files=files,
        supplements_dir=supplements_dir,
        pmid="1",
        converter=converter,
        download_callback=_touch_callback,
        sleep_fn=lambda _s: None,
    )

    assert result.downloaded_count == 1
    [zip_result] = result.file_results
    assert zip_result.extension == ".zip"
    assert zip_result.nested_files
    assert zip_result.converted_chars > 0


def test_download_failure_is_recorded(tmp_path: Path):
    converter = _StubConverter()
    supplements_dir = tmp_path / "supps"
    files = [{"url": "http://x/a.xlsx", "name": "a.xlsx"}]

    def failing_cb(*a, **kw):
        return False

    result = process_supplement_files(
        supp_files=files,
        supplements_dir=supplements_dir,
        pmid="1",
        converter=converter,
        download_callback=failing_cb,
        sleep_fn=lambda _s: None,
    )

    assert result.downloaded_count == 0
    assert len(result.file_results) == 1
    fr = result.file_results[0]
    assert fr.downloaded is False
    assert fr.error == "download_failed"
    assert fr.converted_chars == 0


def test_pdf_with_image_extraction_increments_counter(tmp_path: Path):
    class _PDFConverter(_StubConverter):
        def pdf_to_markdown_with_images(self, p, output_dir):
            self.calls.append(f"pdf+img:{p.name}")
            return f"# pdf+img {p.name}", [{"path": str(output_dir / "f.png")}]

    converter = _PDFConverter()
    supplements_dir = tmp_path / "supps"
    figures_dir = tmp_path / "figs"
    files = [{"url": "http://x/a.pdf", "name": "a.pdf"}]
    result = process_supplement_files(
        supp_files=files,
        supplements_dir=supplements_dir,
        pmid="1",
        converter=converter,
        download_callback=_touch_callback,
        extract_figures=True,
        figures_dir=figures_dir,
        sleep_fn=lambda _s: None,
    )
    assert result.total_figures_extracted == 1
    assert result.file_results[0].figures_extracted == 1


def test_real_zip_extracted_through_format_converter(tmp_path: Path):
    """Smoke test: feed a real ZIP through the real FormatConverter."""
    from harvesting.format_converters import FormatConverter

    converter = FormatConverter()
    bundle = tmp_path / "bundle.zip"
    with zipfile.ZipFile(bundle, "w") as zf:
        zf.writestr("table.csv", "variant,carriers\np.G604S,5\np.A561V,3\n")

    supplements_dir = tmp_path / "supps"
    files = [{"url": "http://x/bundle.zip", "name": "bundle.zip"}]

    def cb(url, file_path, pmid, filename, supp):
        file_path.write_bytes(bundle.read_bytes())
        return True

    result = process_supplement_files(
        supp_files=files,
        supplements_dir=supplements_dir,
        pmid="1",
        converter=converter,
        download_callback=cb,
        sleep_fn=lambda _s: None,
    )

    assert result.downloaded_count == 1
    [zip_result] = result.file_results
    assert zip_result.nested_files
    # Real conversion should surface the variant text from the inner CSV.
    assert "G604S" in result.supplement_markdown
