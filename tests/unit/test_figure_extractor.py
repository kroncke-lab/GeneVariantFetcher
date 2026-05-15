"""Tests for harvesting.figure_extractor.

Covers PMC JATS XML extraction (figures, table-wraps, supplementary-material)
and publisher-style HTML extraction (figure / figcaption / table caption /
supplementary-information sections). Also exercises the markdown renderer.
"""

from __future__ import annotations

import textwrap

import pytest

from harvesting.figure_extractor import (
    CaptionExtractionResult,
    extract_from_html,
    extract_from_jats_xml,
    merge_results,
    render_captions_markdown,
)


# ---------------------------------------------------------------------------
# JATS XML
# ---------------------------------------------------------------------------


JATS_FIXTURE = textwrap.dedent("""\
    <article xmlns:xlink="http://www.w3.org/1999/xlink">
      <body>
        <sec>
          <title>Results</title>
          <p>Body text.</p>
          <fig id="fig1">
            <label>Figure 1</label>
            <caption>
              <title>Pedigree of family carrying KCNH2 G604S</title>
              <p>Filled circles indicate affected probands. The proband is II-1.</p>
            </caption>
            <graphic xlink:href="figs/fig1.jpg"/>
          </fig>
          <table-wrap id="t1">
            <label>Table 1</label>
            <caption>
              <p>List of variants identified in the screening cohort.</p>
            </caption>
          </table-wrap>
        </sec>
        <sec sec-type="supplementary-material">
          <supplementary-material id="S1">
            <label>Supplementary Table S1</label>
            <caption>
              <p>Extended table of 86 KCNH2 variants identified.</p>
            </caption>
            <media xlink:href="supp/S1.xlsx" mimetype="application"
                   mime-subtype="vnd.openxmlformats-officedocument.spreadsheetml.sheet"/>
          </supplementary-material>
        </sec>
      </body>
    </article>
    """)


def test_extract_from_jats_xml_pulls_figure_caption_and_image():
    res = extract_from_jats_xml(JATS_FIXTURE)
    assert len(res.figures) == 1
    fig = res.figures[0]
    assert fig.label == "Figure 1"
    assert "G604S" in fig.title
    assert "II-1" in fig.text
    assert fig.image_url == "figs/fig1.jpg"
    assert fig.figure_id == "fig1"


def test_extract_from_jats_xml_pulls_table_caption():
    res = extract_from_jats_xml(JATS_FIXTURE)
    assert len(res.tables) == 1
    tbl = res.tables[0]
    assert tbl.label == "Table 1"
    assert "screening cohort" in tbl.text


def test_extract_from_jats_xml_pulls_supplementary_material():
    res = extract_from_jats_xml(JATS_FIXTURE)
    assert len(res.supplements) == 1
    supp = res.supplements[0]
    assert supp.label == "Supplementary Table S1"
    assert "86 KCNH2 variants" in supp.text
    assert supp.href == "supp/S1.xlsx"
    assert supp.media_type and "vnd.openxmlformats" in supp.media_type


def test_extract_from_jats_xml_handles_empty_input():
    assert extract_from_jats_xml("").is_empty()
    assert extract_from_jats_xml("not xml at all").is_empty()


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------


HTML_FIXTURE = textwrap.dedent("""\
    <html><body>
      <h1>Article title</h1>
      <p>Some article body.</p>
      <figure id="fig2">
        <img src="https://cdn.example.com/figs/fig2.png" alt=""/>
        <figcaption>
          <strong>Fig. 2.</strong> Pedigree showing the SCN5A R222Q variant
          segregating with disease in three generations.
        </figcaption>
      </figure>
      <table>
        <caption>Table 2. Variants identified in the cohort.</caption>
        <tr><th>Variant</th><th>Carriers</th></tr>
        <tr><td>p.G604S</td><td>5</td></tr>
      </table>
      <section class="tw xbox" id="t3-6124">
        <div class="caption p">
          <p><strong>Table 3</strong> Mutations visible only in a PMC table image.</p>
        </div>
        <p class="img-box">
          <img class="graphic" src="https://cdn.example.com/T3-6124.jpg"/>
        </p>
      </section>
      <section id="supplementary-information">
        <h2>Supplementary Information</h2>
        <p><strong>Supplementary Table S2</strong>: Extended list of 90 missense variants.
          <a href="supp/S2.xlsx">download</a>
        </p>
        <p>Supplementary Methods: detailed protocols.</p>
      </section>
    </body></html>
    """)


def test_extract_from_html_finds_figure():
    res = extract_from_html(HTML_FIXTURE)
    assert len(res.figures) == 1
    fig = res.figures[0]
    assert "R222Q" in fig.text
    assert fig.image_url == "https://cdn.example.com/figs/fig2.png"
    assert fig.figure_id == "fig2"
    assert "Fig" in fig.label or fig.label == "Fig. 2"


def test_extract_from_html_finds_table_caption():
    res = extract_from_html(HTML_FIXTURE)
    assert len(res.tables) == 2
    tbl = res.tables[0]
    assert "G604S" not in tbl.text  # caption only, not body cells
    assert "cohort" in tbl.text
    assert tbl.label.startswith("Table") or tbl.label.lower().startswith("table")


def test_extract_from_html_finds_pmc_image_table_caption_and_image():
    res = extract_from_html(HTML_FIXTURE)
    image_table = next(t for t in res.tables if t.table_id == "t3-6124")
    assert "PMC table image" in image_table.text
    assert image_table.label == "Table 3"
    assert image_table.image_url == "https://cdn.example.com/T3-6124.jpg"


def test_extract_from_html_finds_supplement_descriptions():
    res = extract_from_html(HTML_FIXTURE)
    assert len(res.supplements) >= 1
    labels = [s.label for s in res.supplements]
    assert any("S2" in lbl or "Supplementary Table" in lbl for lbl in labels)
    # The href should be picked up from the link inside.
    hrefs = [s.href for s in res.supplements if s.href]
    assert any(h and "S2.xlsx" in h for h in hrefs)


def test_extract_from_html_handles_empty():
    assert extract_from_html("").is_empty()
    assert extract_from_html("<html></html>").is_empty()


# ---------------------------------------------------------------------------
# Rendering & merging
# ---------------------------------------------------------------------------


def test_render_captions_markdown_includes_section_headings():
    xml_res = extract_from_jats_xml(JATS_FIXTURE)
    md = render_captions_markdown(xml_res)
    assert "## FIGURE CAPTIONS" in md
    assert "## TABLE CAPTIONS" in md
    assert "## SUPPLEMENT DESCRIPTIONS" in md
    assert "G604S" in md
    assert "86 KCNH2 variants" in md


def test_render_captions_markdown_empty_returns_empty_string():
    empty = CaptionExtractionResult()
    assert render_captions_markdown(empty) == ""


def test_merge_results_dedupes_by_text_fingerprint():
    a = extract_from_jats_xml(JATS_FIXTURE)
    b = extract_from_jats_xml(JATS_FIXTURE)
    merged = merge_results(a, b)
    # Same source twice should not double count.
    assert len(merged.figures) == 1
    assert len(merged.tables) == 1
    assert len(merged.supplements) == 1


def test_merge_results_preserves_distinct_items():
    xml_res = extract_from_jats_xml(JATS_FIXTURE)
    html_res = extract_from_html(HTML_FIXTURE)
    merged = merge_results(xml_res, html_res)
    # Distinct sources should accumulate.
    assert len(merged.figures) >= 2
    assert len(merged.tables) >= 3
    assert len(merged.supplements) >= 2
