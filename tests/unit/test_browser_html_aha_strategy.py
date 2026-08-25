"""Fixture-based test for AHAStrategy supplement scraping.

We don't need a real Playwright browser to verify the AHA-specific
supplement scraper — feed the raw HTML and assert the parser extracts
the expected supplement URLs.
"""

from __future__ import annotations

from harvesting.browser_html.dom_extract import (
    extract_body_markdown,
    pick_better_markdown,
)
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


def test_strategy_encodes_legacy_doi_for_path_urls():
    strategy = AHAStrategy()
    doi = "10.1002/(SICI)1098-1004(1999)13:4<301::AID-HUMU7>3.0.CO;2-V"

    encoded = strategy.encode_doi_for_path(doi)

    assert "<" not in encoded
    assert ">" not in encoded
    assert "%3C301::AID-HUMU7%3E" in encoded


def test_readable_html_fallback_preserves_article_tables():
    strategy = AHAStrategy()
    html = """
    <html><body>
      <div class="hlFld-Fulltext">
        <h1>Mutation spectrum</h1>
        <p>KCNH2 variants were found in affected carriers across a large
        long QT syndrome cohort, and the article body includes a variant
        carrier table that must survive HTML conversion for downstream
        extraction.</p>
        <table>
          <caption>Summary of HERG mutations</caption>
          <tr><th>Variant</th><th>Carriers</th></tr>
          <tr><td>R176W</td><td>12</td></tr>
        </table>
      </div>
    </body></html>
    """

    markdown = strategy.extract_readable_html(html, selectors=[".hlFld-Fulltext"])

    assert markdown is not None
    assert "Summary of HERG mutations" in markdown
    assert "| Variant | Carriers |" in markdown
    assert "| R176W | 12 |" in markdown


def test_dom_extractor_preserves_collapsed_table_wrapped_in_figure():
    html = """
    <html><body>
      <section id="bodymatter">
        <div role="paragraph">
          This clinical cohort reports patient-level variants in a table.
          The article body contains enough methodological and clinical context
          for the body selector to qualify as a substantive source. The table
          remains the authoritative patient-level result and must survive even
          when the publisher initially hides most rows behind an expand control.
          Preserving those hidden rows is necessary for complete extraction.
        </div>
        <figure class="table">
          <figcaption>Patient Data</figcaption>
          <div class="collapsible-wrapper collapsed">
            <table>
              <thead><tr><th>Case</th><th>Gene</th><th>Variant</th></tr></thead>
              <tbody>
                <tr><td>31</td><td>RYR2</td><td>P164S</td></tr>
                <tr hidden><td>32</td><td>RYR2</td><td>R414L</td></tr>
              </tbody>
            </table>
          </div>
        </figure>
      </section>
    </body></html>
    """

    markdown = extract_body_markdown(html, ["#bodymatter"])

    assert markdown is not None
    assert "*Figure: Patient Data*" in markdown
    assert "| Case | Gene | Variant |" in markdown
    assert "| 31 | RYR2 | P164S |" in markdown
    assert "| 32 | RYR2 | R414L |" in markdown


def test_dom_extractor_does_not_duplicate_nested_table_rows():
    html = """
    <html><body><section id="bodymatter">
      <p>This article has sufficient body text to qualify for extraction. The
      clinical cohort details remain in the table below, with a nested layout
      table in one cell that must not become a second top-level result table.</p>
      <figure class="table"><figcaption>Patient Data</figcaption>
        <table>
          <tr><th>Case</th><th>Variant</th><th>Note</th></tr>
          <tr><td>1</td><td>P164S</td><td>
            <table><tr><td>nested layout note</td></tr></table>
          </td></tr>
        </table>
      </figure>
    </section></body></html>
    """

    markdown = extract_body_markdown(html, ["#bodymatter"])

    assert markdown is not None
    assert markdown.count("nested layout note") == 1
    assert markdown.count("| Case | Variant | Note |") == 1


def test_markdown_picker_prefers_complete_table_over_slightly_longer_prose():
    primary = "# MAIN TEXT\n\n" + ("Clinical prose without the cohort table. " * 80)
    dom = (
        "# MAIN TEXT\n\n"
        + ("Clinical prose. " * 130)
        + "\n\n| Case | Gene | Variant |\n"
        + "| --- | --- | --- |\n"
        + "| 31 | RYR2 | P164S |\n"
        + "| 32 | RYR2 | R414L |\n"
    )

    assert len(dom) < len(primary)
    assert pick_better_markdown(primary, dom, prefer_more_tables=True) == dom
    assert pick_better_markdown(primary, dom) == primary
