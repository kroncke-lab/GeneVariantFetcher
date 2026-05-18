"""Unit tests for SupplementScraper.scrape_generic_supplements().

Tests the expanded keyword, URL pattern, and file extension matching
without requiring network access.
"""

import pytest

from harvesting.orchestrator import PMCHarvester
from harvesting.supplement_scraper import SupplementScraper


@pytest.fixture
def scraper():
    return SupplementScraper()


def _make_html(links: list[tuple[str, str]]) -> str:
    """Build minimal HTML with <a> tags from (href, text) pairs."""
    anchors = "\n".join(f'<a href="{href}">{text}</a>' for href, text in links)
    return f"<html><body>{anchors}</body></html>"


class TestKeywordMatching:
    """Links matched via text keywords."""

    def test_classic_supplement_keyword(self, scraper):
        html = _make_html([("/files/table.xlsx", "Supplementary Table 1")])
        result = scraper.scrape_generic_supplements(html, "https://example.com/article")
        assert len(result) == 1
        assert result[0]["name"] == "table.xlsx"

    def test_esm_keyword(self, scraper):
        html = _make_html([("/files/esm.pdf", "ESM Figure 1")])
        result = scraper.scrape_generic_supplements(html, "https://example.com/article")
        assert len(result) == 1

    def test_moesm_keyword(self, scraper):
        html = _make_html([("/files/data.csv", "MOESM1 of this article")])
        result = scraper.scrape_generic_supplements(html, "https://example.com/article")
        assert len(result) == 1

    def test_electronic_supplementary_keyword(self, scraper):
        html = _make_html([("/files/data.xlsx", "Electronic Supplementary Material")])
        result = scraper.scrape_generic_supplements(html, "https://example.com/article")
        assert len(result) == 1

    def test_extended_data_keyword(self, scraper):
        html = _make_html([("/files/ext.pdf", "Extended Data Table 3")])
        result = scraper.scrape_generic_supplements(html, "https://example.com/article")
        assert len(result) == 1

    def test_supporting_information_keyword(self, scraper):
        html = _make_html([("/files/si.docx", "Supporting Information S1")])
        result = scraper.scrape_generic_supplements(html, "https://example.com/article")
        assert len(result) == 1


class TestURLPatternMatching:
    """Links matched via URL patterns (the key improvement)."""

    def test_download_supplement_endpoint(self, scraper):
        """Oxford/Wiley downloadSupplement pattern — previously missed."""
        html = _make_html(
            [
                ("/downloadSupplement?file=data.xlsx", "Download"),
            ]
        )
        result = scraper.scrape_generic_supplements(
            html, "https://academic.oup.com/article"
        )
        assert len(result) == 1

    def test_media_objects_url(self, scraper):
        """Springer /MediaObjects/ pattern."""
        html = _make_html(
            [
                (
                    "/article/10.1038/s001/MediaObjects/12345_2020_1_MOESM1_ESM.pdf",
                    "Download",
                ),
            ]
        )
        result = scraper.scrape_generic_supplements(
            html, "https://link.springer.com/article"
        )
        assert len(result) == 1
        assert "MOESM1" in result[0]["name"]

    def test_mmc_pattern_in_url(self, scraper):
        """Elsevier mmc pattern detected via URL even without keyword text."""
        html = _make_html(
            [
                ("/science/article/pii/S00029297/mmc1.xlsx", "Table S1"),
            ]
        )
        result = scraper.scrape_generic_supplements(
            html, "https://www.sciencedirect.com"
        )
        assert len(result) == 1

    def test_suppl_path_segment(self, scraper):
        html = _make_html(
            [
                ("/doi/suppl/10.1161/CIR.123/suppl_file/data.pdf", "Click here"),
            ]
        )
        result = scraper.scrape_generic_supplements(
            html, "https://ahajournals.org/article"
        )
        assert len(result) == 1

    def test_supplementary_suffix_in_url(self, scraper):
        html = _make_html(
            [
                ("/articles/supplementary_data.xlsx", "Data"),
            ]
        )
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert len(result) == 1

    def test_pmc_relative_pdf_link_is_normalized(self, scraper):
        base = "https://pmc.ncbi.nlm.nih.gov/articles/PMC3566559/"

        result = scraper._normalize_pmc_url("pdf/nihms372211.pdf", base)

        assert result == (
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3566559/pdf/nihms372211.pdf"
        )

    def test_ncbi_ftp_supplement_urls_use_https(self, tmp_path):
        class Response:
            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=8192):
                yield b"%PDF-1.4\n%%EOF\n"

        class Session:
            def __init__(self):
                self.urls = []

            def get(self, url, **kwargs):
                self.urls.append(url)
                return Response()

        harvester = PMCHarvester(output_dir=tmp_path, gene_symbol="KCNH2")
        harvester.session = Session()

        ok = harvester.download_supplement(
            "ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_pdf/64/67/test.pdf",
            tmp_path / "test.pdf",
            "24596401",
            "test.pdf",
        )

        assert ok is True
        assert harvester.session.urls == [
            "https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_pdf/64/67/test.pdf"
        ]

    def test_pmc_relative_supplement_variants_are_normalized_for_download(
        self, tmp_path
    ):
        class Response:
            def __init__(self, body):
                self.body = body

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=8192):
                yield self.body

        class Session:
            def __init__(self):
                self.urls = []

            def get(self, url, **kwargs):
                self.urls.append(url)
                assert url.startswith("https://")
                if len(self.urls) == 1:
                    return Response(b"<html>expired PMC redirect</html>")
                return Response(b"%PDF-1.4\n%%EOF\n")

        base_url = "https://pmc.ncbi.nlm.nih.gov/articles/PMC5970051/"
        harvester = PMCHarvester(output_dir=tmp_path, gene_symbol="SCN5A")
        harvester.session = Session()
        harvester.scraper.get_pmc_supplement_url_variants = lambda url, base: [
            "/articles/instance/5970051/bin/test.pdf"
        ]

        ok = harvester.download_supplement(
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5970051/bin/test.pdf",
            tmp_path / "test.pdf",
            "29514831",
            "test.pdf",
            base_url=base_url,
            original_url="/articles/instance/5970051/bin/test.pdf",
        )

        assert ok is True
        assert harvester.session.urls == [
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5970051/bin/test.pdf",
            "https://pmc.ncbi.nlm.nih.gov/articles/instance/5970051/bin/test.pdf",
        ]


class TestFileExtensionMatching:
    """Links matched via file extension (with query-param stripping)."""

    def test_plain_extension(self, scraper):
        html = _make_html([("/files/table.xlsx", "Click to download")])
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert len(result) == 1

    def test_extension_with_query_params(self, scraper):
        """Previously missed: extension hidden behind query params."""
        html = _make_html(
            [
                ("/files/table.xlsx?token=abc123", "Download table"),
            ]
        )
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert len(result) == 1
        assert result[0]["name"] == "table.xlsx"

    def test_tsv_extension(self, scraper):
        html = _make_html([("/data/variants.tsv", "Variant list")])
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert len(result) == 1

    def test_xls_extension(self, scraper):
        html = _make_html([("/data/old_table.xls", "Legacy table")])
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert len(result) == 1


class TestFiltering:
    """Deduplication and ID-like filename filtering."""

    def test_deduplicates_by_filename(self, scraper):
        html = _make_html(
            [
                ("/path1/data.xlsx", "Supplement 1"),
                ("/path2/data.xlsx", "Supplement 2"),
            ]
        )
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert len(result) == 1

    def test_skips_pmcid_filenames(self, scraper):
        html = _make_html([("/articles/PMC3049907", "Supplementary info")])
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert len(result) == 0

    def test_skips_id_like_filenames_for_keyword_only_matches(self, scraper):
        html = _make_html([("/articles/ABC123", "Supplementary data")])
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert len(result) == 0

    def test_no_match_for_irrelevant_link(self, scraper):
        html = _make_html([("/about", "About this journal")])
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert len(result) == 0

    def test_skips_publisher_catalog_and_media_pack_pdfs(self, scraper):
        html = _make_html(
            [
                ("/images/pdf/Media-Pack-2024.pdf", "Media Pack"),
                ("/images/pdf/JournalCatalog2026.pdf", "Journal catalog"),
                ("/files/supplementary-table-1.xlsx", "Supplementary Table 1"),
            ]
        )

        result = scraper.scrape_generic_supplements(
            html, "https://www.eurekaselect.com/article"
        )

        assert [item["name"] for item in result] == ["supplementary-table-1.xlsx"]

    def test_skips_javascript_pseudo_links(self, scraper):
        html = _make_html([("javascript:;", "Supplementary material")])
        result = scraper.scrape_generic_supplements(html, "https://example.com")
        assert result == []


class TestFigshareExpansion:
    """Figshare landing pages should become direct downloadable files."""

    def test_generic_figshare_landing_expands_to_files(self, scraper, monkeypatch):
        monkeypatch.setattr(
            scraper,
            "_fetch_figshare_article",
            lambda article_id: {
                "files": [
                    {
                        "id": 6568608,
                        "name": "Extanded tables.doc",
                        "size": 162304,
                        "download_url": "https://ndownloader.figshare.com/files/6568608",
                        "mimetype": "application/msword",
                    }
                ]
            },
        )
        html = _make_html(
            [
                (
                    "https://figshare.com/articles/dataset/Supplementary_Material/4056510",
                    "Supplementary Material",
                )
            ]
        )

        result = scraper.scrape_generic_supplements(html, "https://karger.com/article")

        assert result == [
            {
                "url": "https://ndownloader.figshare.com/files/6568608",
                "name": "Extanded tables.doc",
                "figshare_article_id": "4056510",
                "figshare_file_id": 6568608,
                "size": 162304,
                "mimetype": "application/msword",
            }
        ]

    def test_karger_figshare_landing_expands_and_skips_javascript(
        self, scraper, monkeypatch
    ):
        monkeypatch.setattr(
            scraper,
            "_fetch_figshare_article",
            lambda article_id: {
                "files": [
                    {
                        "id": 1,
                        "name": "variant-table.doc",
                        "download_url": "https://ndownloader.figshare.com/files/1",
                    },
                    {
                        "id": 2,
                        "name": "participants.doc",
                        "download_url": "https://ndownloader.figshare.com/files/2",
                    },
                ]
            },
        )
        html = _make_html(
            [
                ("javascript:;", "Supplementary Material"),
                (
                    "https://figshare.com/articles/dataset/Supplementary_Material/4056510",
                    "Supplementary Material",
                ),
            ]
        )

        result = scraper.scrape_karger_supplements(html, "https://karger.com/article")

        assert [r["name"] for r in result] == [
            "variant-table.doc",
            "participants.doc",
        ]


class TestPMCPages:
    """PMC-specific behavior: original_url and base_url fields."""

    def test_pmc_page_adds_extra_fields(self, scraper):
        html = _make_html([("/articles/instance/123/bin/supp.pdf", "Supplement")])
        base = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC123456/"
        result = scraper.scrape_generic_supplements(html, base)
        assert len(result) == 1
        assert "original_url" in result[0]
        assert "base_url" in result[0]
        assert result[0]["base_url"] == base

    def test_non_pmc_page_no_extra_fields(self, scraper):
        html = _make_html([("/files/supp.pdf", "Supplement")])
        result = scraper.scrape_generic_supplements(html, "https://example.com/article")
        assert len(result) == 1
        assert "original_url" not in result[0]
        assert "base_url" not in result[0]


class TestHelperMethods:
    """Tests for _is_supplement_url and _has_file_extension."""

    def test_is_supplement_url_positive(self, scraper):
        assert scraper._is_supplement_url("/downloadSupplement?file=x")
        assert scraper._is_supplement_url("/article/MediaObjects/file.pdf")
        assert scraper._is_supplement_url("/path/mmc1.xlsx")
        assert scraper._is_supplement_url("/path/moesm3_data.pdf")

    def test_is_supplement_url_negative(self, scraper):
        assert not scraper._is_supplement_url("/about")
        assert not scraper._is_supplement_url("/article/main.html")

    def test_has_file_extension_plain(self, scraper):
        assert scraper._has_file_extension("/files/data.xlsx")
        assert scraper._has_file_extension("/files/report.pdf")

    def test_has_file_extension_with_query(self, scraper):
        assert scraper._has_file_extension("/files/data.xlsx?token=abc")

    def test_has_file_extension_negative(self, scraper):
        assert not scraper._has_file_extension("/articles/PMC123")
        assert not scraper._has_file_extension("/downloadSupplement?file=x")
