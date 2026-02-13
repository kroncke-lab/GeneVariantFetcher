"""Tests for supplement processing service."""

from harvesting.supplement_processing_service import process_supplement_files


class _Converter:
    def excel_to_markdown(self, path):
        return "[excel]"

    def docx_to_markdown(self, path):
        return "[docx]"

    def doc_to_markdown(self, path):
        return "[doc]"

    def pdf_to_markdown(self, path):
        return "[pdf]"

    def pdf_to_markdown_with_images(self, path, output_dir):
        return "[pdf+img]", ["a.png", "b.png"]


class _Logger:
    def info(self, *args, **kwargs):
        return None


def test_process_supplement_files_routes_extensions(tmp_path):
    converter = _Converter()

    supp_files = [
        {"url": "u1", "name": "a.xlsx"},
        {"url": "u2", "name": "b.docx"},
        {"url": "u3", "name": "c.doc"},
        {"url": "u4", "name": "d.pdf"},
        {"url": "u5", "name": "e.txt"},
        {"url": "u6", "name": "f.bin"},
    ]

    def download_callback(url, file_path, pmid, filename, supp):
        if file_path.suffix == ".txt":
            file_path.write_text("text-file", encoding="utf-8")
        else:
            file_path.write_text("x", encoding="utf-8")
        return True

    result = process_supplement_files(
        supp_files=supp_files,
        supplements_dir=tmp_path / "supp",
        pmid="123",
        converter=converter,
        download_callback=download_callback,
        extract_figures=False,
        figures_dir=None,
        logger=_Logger(),
        sleep_seconds=0.0,
        sleep_fn=lambda _: None,
    )

    assert result.downloaded_count == 6
    assert "[excel]" in result.supplement_markdown
    assert "[docx]" in result.supplement_markdown
    assert "[doc]" in result.supplement_markdown
    assert "[pdf]" in result.supplement_markdown
    assert "text-file" in result.supplement_markdown
    assert "[File available at:" in result.supplement_markdown


def test_process_supplement_files_extracts_pdf_images(tmp_path):
    converter = _Converter()
    figures_dir = tmp_path / "figs"

    def download_callback(url, file_path, pmid, filename, supp):
        file_path.write_text("pdf", encoding="utf-8")
        return True

    result = process_supplement_files(
        supp_files=[{"url": "u", "name": "paper.pdf"}],
        supplements_dir=tmp_path / "supp",
        pmid="123",
        converter=converter,
        download_callback=download_callback,
        extract_figures=True,
        figures_dir=figures_dir,
        logger=_Logger(),
        sleep_seconds=0.0,
        sleep_fn=lambda _: None,
    )

    assert result.downloaded_count == 1
    assert result.total_figures_extracted == 2
    assert "[pdf+img]" in result.supplement_markdown
