"""Tests for supplement processing service."""

from harvesting.supplement_processing_service import (
    SUPPLEMENT_IDENTITY_UNVERIFIED,
    pdf_supplement_identity,
    process_supplement_files,
)


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

    def warning(self, *args, **kwargs):
        return None


def test_process_supplement_files_routes_extensions(tmp_path):
    converter = _Converter()

    # The .pdf entry uses the Elsevier mmc naming so the identity gate accepts
    # it (a bare "d.pdf" with unrelated stub text is correctly quarantined).
    supp_files = [
        {"url": "u1", "name": "a.xlsx"},
        {"url": "u2", "name": "b.docx"},
        {"url": "u3", "name": "c.doc"},
        {"url": "u4", "name": "mmc1.pdf"},
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
        supp_files=[{"url": "u", "name": "supplementary_figures.pdf"}],
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


# ---------------------------------------------------------------------------
# PDF supplement identity gate (KCNQ1 31520628: cited CDC vital-statistics
# PDFs were captured as supplements and folded into context)
# ---------------------------------------------------------------------------

_CDC_LIKE_TEXT = (
    "National Vital Statistics Reports Volume 60, Number 7. "
    "Deaths: preliminary statistics for 2011. Division of Vital Statistics. "
    "Mortality patterns are described by age, sex, race, and state of residence."
)


def test_pdf_identity_rejects_unrelated_document():
    ok, reason = pdf_supplement_identity(
        text_head=_CDC_LIKE_TEXT,
        filename="nvsr60_07.pdf",
        source_url="https://stats.example.gov/reports/nvsr60/nvsr60_07.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
    )
    assert ok is False
    assert "no supplement marker" in reason


def test_pdf_identity_accepts_article_doi_in_text():
    ok, reason = pdf_supplement_identity(
        text_head="Extra cohort tables. doi: 10.1016/j.ajog.2019.09.004\n" + "x " * 50,
        filename="cohort_tables.pdf",
        source_url="https://journal.example.org/files/cohort_tables.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
    )
    assert ok is True
    assert "DOI" in reason


def test_pdf_identity_accepts_supplement_marker_in_text():
    ok, _ = pdf_supplement_identity(
        text_head="Supplementary Table S3. Variant carrier counts by family.",
        filename="download.pdf",
        pmid="31520628",
    )
    assert ok is True


def test_pdf_identity_accepts_publisher_supplement_naming():
    # Elsevier mmc naming and /doi/suppl/ paths identify the object itself as
    # supplementary material, whatever its first page says.
    ok, _ = pdf_supplement_identity(
        text_head="Table of odds ratios only.", filename="mmc1.pdf", pmid="1"
    )
    assert ok is True
    ok, _ = pdf_supplement_identity(
        text_head="Table of odds ratios only.",
        filename="x.pdf",
        source_url="https://pub.example.org/doi/suppl/10.1056/NEJMoa042786/x.pdf",
        pmid="1",
    )
    assert ok is True
    # JATS media ids (pone.0012345.s001.pdf) are supplement-minted too.
    ok, _ = pdf_supplement_identity(
        text_head="Table of odds ratios only.",
        filename="pone.0012345.s001.pdf",
        pmid="1",
    )
    assert ok is True


def test_pdf_identity_accepts_pmid_in_text():
    ok, reason = pdf_supplement_identity(
        text_head="Companion data for PMID 31520628 cohort analysis.",
        filename="companion.pdf",
        pmid="31520628",
    )
    assert ok is True
    assert "PMID" in reason


def test_process_supplement_files_quarantines_unverified_pdf(tmp_path):
    """The unrelated PDF stays on disk but its text is excluded from context."""

    class _PdfConverter(_Converter):
        def pdf_to_markdown(self, path):
            return _CDC_LIKE_TEXT

    def download_callback(url, file_path, pmid, filename, supp):
        file_path.write_bytes(b"%PDF-1.4 fake")
        return True

    result = process_supplement_files(
        supp_files=[
            {
                "url": "https://stats.example.gov/reports/nvsr60_07.pdf",
                "name": "nvsr60_07.pdf",
            }
        ],
        supplements_dir=tmp_path / "supp",
        pmid="31520628",
        converter=_PdfConverter(),
        download_callback=download_callback,
        logger=_Logger(),
        sleep_seconds=0.0,
        sleep_fn=lambda _: None,
    )

    assert result.downloaded_count == 1
    assert _CDC_LIKE_TEXT not in result.supplement_markdown
    assert "nvsr60_07.pdf" not in result.supplement_markdown
    per_file = result.file_results[0]
    assert per_file.identity_unverified is True
    assert per_file.error == SUPPLEMENT_IDENTITY_UNVERIFIED
    # Never delete downloaded source.
    assert (tmp_path / "supp" / "nvsr60_07.pdf").exists()


def test_process_supplement_files_folds_pdf_that_carries_doi(tmp_path):
    class _PdfConverter(_Converter):
        def pdf_to_markdown(self, path):
            return "Cohort tables. doi:10.1016/j.ajog.2019.09.004"

    def download_callback(url, file_path, pmid, filename, supp):
        file_path.write_bytes(b"%PDF-1.4 fake")
        return True

    result = process_supplement_files(
        supp_files=[
            {"url": "https://journal.example.org/f/cohort.pdf", "name": "cohort.pdf"}
        ],
        supplements_dir=tmp_path / "supp",
        pmid="31520628",
        converter=_PdfConverter(),
        download_callback=download_callback,
        logger=_Logger(),
        sleep_seconds=0.0,
        sleep_fn=lambda _: None,
        doi="10.1016/j.ajog.2019.09.004",
    )

    per_file = result.file_results[0]
    assert per_file.identity_unverified is False
    assert per_file.error is None
    assert "cohort.pdf" in result.supplement_markdown
    assert "10.1016/j.ajog.2019.09.004" in result.supplement_markdown
