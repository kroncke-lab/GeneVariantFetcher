"""Tests for supplement processing service."""

from harvesting.supplement_processing_service import (
    SUPPLEMENT_IDENTITY_QUARANTINE,
    SUPPLEMENT_IDENTITY_UNKNOWN,
    SUPPLEMENT_IDENTITY_UNVERIFIED,
    SUPPLEMENT_IDENTITY_VERIFIED,
    SUPPLEMENT_PROVENANCE_REFERENCE_SECTION,
    SUPPLEMENT_PROVENANCE_SUPPLEMENT_CONTAINER,
    pdf_supplement_verdict,
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

#: A cited foreign report: a numbered serial that also prints its OWN doi,
#: which contradicts the parent article's. The masthead alone only annotates
#: (a paper's own PDF looks the same); the foreign DOI is what condemns it.
_CDC_LIKE_TEXT = (
    "National Vital Statistics Reports Volume 60, Number 7. "
    "Deaths: preliminary statistics for 2011. Division of Vital Statistics. "
    "doi:10.15620/cdc.nvsr6007. "
    "Mortality patterns are described by age, sex, race, and state of residence."
)


def test_pdf_identity_rejects_unrelated_document():
    verdict, reason = pdf_supplement_verdict(
        text_head=_CDC_LIKE_TEXT,
        filename="nvsr60_07.pdf",
        source_url="https://stats.example.gov/reports/nvsr60/nvsr60_07.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
    )
    assert verdict == SUPPLEMENT_IDENTITY_QUARANTINE
    assert "doi" in reason.lower()


def test_pdf_identity_accepts_article_doi_in_text():
    verdict, reason = pdf_supplement_verdict(
        text_head="Extra cohort tables. doi: 10.1016/j.ajog.2019.09.004\n" + "x " * 50,
        filename="cohort_tables.pdf",
        source_url="https://journal.example.org/files/cohort_tables.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
    )
    assert verdict == SUPPLEMENT_IDENTITY_VERIFIED
    assert "DOI" in reason


def test_pdf_identity_accepts_supplement_marker_in_text():
    verdict, _ = pdf_supplement_verdict(
        text_head="Supplementary Table S3. Variant carrier counts by family.",
        filename="download.pdf",
        pmid="31520628",
    )
    assert verdict == SUPPLEMENT_IDENTITY_VERIFIED


def test_pdf_identity_accepts_publisher_supplement_naming():
    # Elsevier mmc naming and /doi/suppl/ paths identify the object itself as
    # supplementary material, whatever its first page says.
    verdict, _ = pdf_supplement_verdict(
        text_head="Table of odds ratios only.", filename="mmc1.pdf", pmid="1"
    )
    assert verdict == SUPPLEMENT_IDENTITY_VERIFIED
    verdict, _ = pdf_supplement_verdict(
        text_head="Table of odds ratios only.",
        filename="x.pdf",
        source_url="https://pub.example.org/doi/suppl/10.1056/NEJMoa042786/x.pdf",
        pmid="1",
    )
    assert verdict == SUPPLEMENT_IDENTITY_VERIFIED
    # JATS media ids (pone.0012345.s001.pdf) are supplement-minted too.
    verdict, _ = pdf_supplement_verdict(
        text_head="Table of odds ratios only.",
        filename="pone.0012345.s001.pdf",
        pmid="1",
    )
    assert verdict == SUPPLEMENT_IDENTITY_VERIFIED


def test_pdf_identity_accepts_pmid_in_text():
    verdict, reason = pdf_supplement_verdict(
        text_head="Companion data for PMID 31520628 cohort analysis.",
        filename="companion.pdf",
        pmid="31520628",
    )
    assert verdict == SUPPLEMENT_IDENTITY_VERIFIED
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
        doi="10.1016/j.ajog.2019.09.004",
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


# ---------------------------------------------------------------------------
# Fail-open-on-absence: only a contradiction quarantines
# ---------------------------------------------------------------------------

_SCANNED_TABLE_TEXT = "Table 1\n\n" + "c.1A>G 3 carriers\n" * 4


def test_pdf_identity_folds_unmarked_scan_unverified():
    """A scanned mutation table named download.pdf carries no marker at all."""
    verdict, reason = pdf_supplement_verdict(
        text_head=_SCANNED_TABLE_TEXT,
        filename="download.pdf",
        source_url="https://journal.example.org/files/download.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
    )
    assert verdict == SUPPLEMENT_IDENTITY_UNKNOWN
    assert "too little extractable text" in reason


def test_pdf_identity_folds_boilerplate_first_page_unverified():
    verdict, _ = pdf_supplement_verdict(
        text_head="Downloaded from journals.example.org by guest. " * 40,
        filename="nihms-123456.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
    )
    assert verdict == SUPPLEMENT_IDENTITY_UNKNOWN


def test_pdf_identity_quarantines_foreign_doi():
    verdict, reason = pdf_supplement_verdict(
        text_head="Odds ratios by cohort. doi:10.1038/ng.3021\n" + "x " * 400,
        filename="table1.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
    )
    assert verdict == SUPPLEMENT_IDENTITY_QUARANTINE
    assert "foreign DOI" in reason


def test_pdf_identity_tolerates_doi_truncated_by_the_head_boundary():
    """KCNQ1 18174212's own PDF prints its DOI across the head cut.

    ``10.1113/jphysi`` is the article DOI rendered short, not a foreign one.
    """
    verdict, _ = pdf_supplement_verdict(
        text_head="Kv7.1 (KCNQ1) properties. DOI:10.1113/jphysi",
        filename="tjp0586-1785.pdf",
        pmid="18174212",
        doi="10.1113/jphysiol.2007.148254",
    )
    assert verdict == SUPPLEMENT_IDENTITY_UNKNOWN


def test_masthead_is_detected_below_the_first_page():
    """The masthead repeats as a running header, so its first occurrence
    depends on how page one extracts. Scanning only the opening block saw
    KCNQ1 31520628's nvsr60_07.pdf (offset 992) and missed its identical twin
    nvsr60_08.pdf (offset 4,006) — the same publication, the same defect."""
    text = (
        "Fetal and Perinatal Mortality, United States, 2006\n"
        + "Abstract text about maternal age and race. " * 90
        + "\nNational Vital Statistics Reports, Vol. 60, No. 8, August 28, 2012\n"
        + "Table 1. Fetal deaths by period of gestation.\n"
    )
    verdict, reason = pdf_supplement_verdict(
        text_head=text,
        filename="nvsr60_08.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
    )
    assert verdict == SUPPLEMENT_IDENTITY_UNKNOWN
    assert "numbered serial issue" in reason


def test_a_papers_own_article_pdf_is_never_quarantined():
    """LMNA 17334235's own full text opens with a journal masthead and prints
    neither its DOI nor its PMID; the artifacts record carries no title. A
    masthead rule that condemned on absent evidence would drop it."""
    text = (
        "EXPERIMENTAL and MOLECULAR MEDICINE, Vol. 39, No. 1, "
        "114-120, February 2007\n"
        "Lamin A/C mutations associated with familial and sporadic cases of "
        "dilated cardiomyopathy in Koreans\n" + "Body prose. " * 200
    )
    verdict, _ = pdf_supplement_verdict(
        text_head=text,
        filename="emm200713.pdf",
        pmid="17334235",
        doi="10.1038/emm.2007.13",
    )
    assert verdict != SUPPLEMENT_IDENTITY_QUARANTINE


def test_masthead_does_not_condemn_a_document_naming_this_article():
    """A journal's own article PDF carries a running header too, but it also
    prints the article DOI somewhere; a cited foreign serial never does."""
    text = (
        "The Journal of Physiology\n"
        + "Methods section prose. " * 120
        + "\nJ Physiol, Vol. 586, No. 7, 2008\n"
        + "https://doi.org/10.1113/jphysiol.2007.148254\n"
    )
    verdict, _ = pdf_supplement_verdict(
        text_head=text,
        filename="tjp0586-1785.pdf",
        pmid="18174212",
        doi="10.1113/jphysiol.2007.148254",
    )
    assert verdict != SUPPLEMENT_IDENTITY_QUARANTINE


def test_pdf_identity_quarantines_serial_masthead():
    """KCNQ1 31520628's CDC reports: no DOI anywhere, masthead front matter."""
    verdict, reason = pdf_supplement_verdict(
        text_head=(
            "National Vital\nStatistics Reports\n"
            "Volume 60, Number 7 June 20, 2012\n"
            "Estimated Pregnancy Rates for the United States, 1990-2008\n"
        ),
        filename="nvsr60_07.pdf",
        source_url="https://www.cdc.gov/nchs/data/nvsr/nvsr60/nvsr60_07.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
    )
    # A masthead only annotates: it cannot be told apart from a paper's own
    # article PDF using the anchors this path has (see LMNA 17334235).
    assert verdict == SUPPLEMENT_IDENTITY_UNKNOWN
    assert "numbered serial issue" in reason


def test_pdf_identity_keeps_article_pdf_without_an_issue_number():
    """A journal article PDF's own header is a volume and date, not a masthead.

    KCNQ1 12770875's on-disk ``3679.pdf`` opens exactly like this; requiring
    the issue number is what keeps it out of the masthead rule.
    """
    verdict, _ = pdf_supplement_verdict(
        text_head=(
            "Biophysical Journal Volume 84 June 2003 3679-3689\n"
            "Molecular determinants of KCNQ1 channel gating. " + "text " * 200
        ),
        filename="3679.pdf",
        pmid="12770875",
        doi="10.1016/S0006-3495(03)75097-8",
    )
    assert verdict == SUPPLEMENT_IDENTITY_UNKNOWN


def test_pdf_identity_quarantines_reference_section_provenance():
    verdict, reason = pdf_supplement_verdict(
        text_head="Deaths: preliminary statistics. " + "text " * 200,
        filename="report_2011.pdf",
        pmid="31520628",
        doi="10.1016/j.ajog.2019.09.004",
        provenance=SUPPLEMENT_PROVENANCE_REFERENCE_SECTION,
    )
    assert verdict == SUPPLEMENT_IDENTITY_QUARANTINE
    assert "reference/bibliography markup" in reason


def test_pdf_identity_accepts_supplementary_container_provenance():
    verdict, reason = pdf_supplement_verdict(
        text_head="Table of odds ratios only.",
        filename="download.pdf",
        pmid="31520628",
        provenance=SUPPLEMENT_PROVENANCE_SUPPLEMENT_CONTAINER,
    )
    assert verdict == SUPPLEMENT_IDENTITY_VERIFIED
    assert "supplementary-material container" in reason


def test_process_supplement_files_folds_unverified_pdf_with_flag(tmp_path):
    """The inverse of the quarantine path: unmarked source is kept, not dropped."""

    class _PdfConverter(_Converter):
        def pdf_to_markdown(self, path):
            return _SCANNED_TABLE_TEXT

    def download_callback(url, file_path, pmid, filename, supp):
        file_path.write_bytes(b"%PDF-1.4 fake")
        return True

    result = process_supplement_files(
        supp_files=[
            {
                "url": "https://journal.example.org/f/download.pdf",
                "name": "download.pdf",
            }
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
    assert per_file.identity_unverified is True
    assert per_file.error is None
    assert "download.pdf" in result.supplement_markdown
    assert "c.1A>G" in result.supplement_markdown


def test_process_supplement_files_quarantines_reference_section_pdf(tmp_path):
    """Reference-list provenance rides in on the scraped entry's source label."""

    class _PdfConverter(_Converter):
        def pdf_to_markdown(self, path):
            return "Deaths: preliminary statistics. " + "text " * 200

    def download_callback(url, file_path, pmid, filename, supp):
        file_path.write_bytes(b"%PDF-1.4 fake")
        return True

    result = process_supplement_files(
        supp_files=[
            {
                "url": "https://stats.example.gov/report_2011.pdf",
                "name": "report_2011.pdf",
                "source": SUPPLEMENT_PROVENANCE_REFERENCE_SECTION,
            }
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
    assert per_file.identity_unverified is True
    assert per_file.error == SUPPLEMENT_IDENTITY_UNVERIFIED
    assert "report_2011.pdf" not in result.supplement_markdown
    assert (tmp_path / "supp" / "report_2011.pdf").exists()
