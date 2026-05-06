"""Fixture-based test for AHAStrategy supplement scraping.

We don't need a real Playwright browser to verify the AHA-specific
supplement scraper — feed the raw HTML and assert the parser extracts
the expected supplement URLs.
"""

from __future__ import annotations

from harvesting.browser_html.strategies.aha import AHAStrategy


AHA_HTML_FIXTURE = """
<html>
<head><title>Test article</title></head>
<body>
<article class="hlFld-Fulltext">
  <h1>Title</h1>
  <p>Body text...</p>
  <section id="supplementary-materials">
    <h2>Supplementary Materials</h2>
    <ul>
      <li>
        <a href="/doi/suppl/10.1161/CIRCRESAHA.118.123456/suppl_file/123456_supplement.pdf">
          Supplemental Material PDF
        </a>
      </li>
      <li>
        <a href="/doi/suppl/10.1161/CIRCRESAHA.118.123456/suppl_file/123456_table_s1.xlsx">
          Table S1
        </a>
      </li>
      <li>
        <a href="https://example.com/cite-this">Cite this article</a>
      </li>
    </ul>
  </section>
  <a href="#fig1">Figure 1</a>
</article>
</body>
</html>
"""


def test_aha_supplement_scraper_finds_pdf_and_xlsx():
    strategy = AHAStrategy()
    base_url = "https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.118.123456"
    supps = strategy._scrape_aha_supplements(AHA_HTML_FIXTURE, base_url)

    assert len(supps) == 2
    urls = [s["url"] for s in supps]
    assert any(u.endswith("_supplement.pdf") for u in urls)
    assert any(u.endswith("_table_s1.xlsx") for u in urls)


def test_aha_supplement_scraper_resolves_relative_urls():
    strategy = AHAStrategy()
    base_url = "https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.118.123456"
    supps = strategy._scrape_aha_supplements(AHA_HTML_FIXTURE, base_url)

    for supp in supps:
        assert supp["url"].startswith("https://www.ahajournals.org/doi/suppl/")


def test_aha_supplement_scraper_skips_non_supplement_links():
    strategy = AHAStrategy()
    base_url = "https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.118.123456"
    supps = strategy._scrape_aha_supplements(AHA_HTML_FIXTURE, base_url)

    # The "Cite this article" link should not be picked up — it isn't under
    # /doi/suppl/ and doesn't have a known file extension.
    for s in supps:
        assert "cite-this" not in s["url"]


def test_aha_filename_extracted_from_url():
    strategy = AHAStrategy()
    base_url = "https://www.ahajournals.org/doi/10.1161/test"
    supps = strategy._scrape_aha_supplements(AHA_HTML_FIXTURE, base_url)

    for supp in supps:
        # Filename should come from the URL path, not the link text.
        assert supp["name"].endswith((".pdf", ".xlsx"))


def test_aha_strategy_metadata():
    strategy = AHAStrategy()
    assert strategy.NAME == "aha"
    assert "10.1161" in strategy.DOI_PREFIXES
    assert "ahajournals.org" in strategy.DOMAINS
    assert strategy.EMBARGO_MONTHS == 12


def test_aha_strategy_matches_doi():
    strategy = AHAStrategy()
    assert strategy.matches("10.1161/CIRCRESAHA.118.123456")
    assert not strategy.matches("10.1093/europace/euaa067")


def test_aha_strategy_matches_url():
    strategy = AHAStrategy()
    assert strategy.matches("", "https://www.ahajournals.org/doi/abs/10.1161/...")
    assert not strategy.matches("", "https://onlinelibrary.wiley.com/doi/...")
