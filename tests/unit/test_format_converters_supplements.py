"""Tests for new supplement-format handlers in format_converters.

Covers the additive surface added for the figure/supplement overhaul:
TSV, HTML supplement, XML supplement, and ZIP nested-extraction.
"""

from __future__ import annotations

import textwrap
import zipfile
from pathlib import Path
from unittest.mock import patch

import pytest

from harvesting.format_converters import FormatConverter


@pytest.fixture()
def converter() -> FormatConverter:
    return FormatConverter()


def test_tsv_to_markdown_outputs_pipe_table(converter, tmp_path: Path):
    p = tmp_path / "variants.tsv"
    p.write_text("variant\tcarriers\nG604S\t5\nA561V\t3\n", encoding="utf-8")
    md = converter.tsv_to_markdown(p)
    assert "G604S" in md
    assert "carriers" in md
    assert "|" in md  # markdown table


def test_tsv_to_markdown_falls_back_to_raw_text_for_ragged(tmp_path: Path):
    p = tmp_path / "weird.tsv"
    # Single-column "TSV" — pandas will still parse, but content survives.
    p.write_text("note\nfoo\nbar\n", encoding="utf-8")
    md = FormatConverter().tsv_to_markdown(p)
    assert "foo" in md and "bar" in md


def test_pdf_to_markdown_accepts_string_path(tmp_path: Path):
    p = tmp_path / "not_a_pdf.pdf"
    p.write_text("not a pdf", encoding="utf-8")

    md = FormatConverter().pdf_to_markdown(str(p))

    assert "[Invalid PDF file: not_a_pdf.pdf]" in md


def test_html_supplement_to_markdown_extracts_table_rows(converter, tmp_path: Path):
    html = textwrap.dedent("""\
        <html><head><title>Supplementary Table S2</title></head>
        <body>
        <table>
          <caption>Variants identified in cohort</caption>
          <tr><th>Variant</th><th>Carriers</th></tr>
          <tr><td>p.G604S</td><td>5</td></tr>
          <tr><td>p.A561V</td><td>3</td></tr>
        </table>
        <p>Notes: All variants are heterozygous.</p>
        <script>alert('drop me')</script>
        </body></html>
        """)
    p = tmp_path / "supp.html"
    p.write_text(html, encoding="utf-8")
    md = converter.html_supplement_to_markdown(p)
    assert "G604S" in md
    assert "A561V" in md
    assert "drop me" not in md  # script removed
    assert "heterozygous" in md
    assert "|" in md  # contains markdown table syntax


def test_xml_to_markdown_preserves_jats_table_wrap_rows(converter):
    xml = textwrap.dedent("""\
        <article>
          <front>
            <article-meta>
              <title-group><article-title>KCNH2 cohort</article-title></title-group>
              <abstract><p>Clinical cohort with KCNH2 variants.</p></abstract>
            </article-meta>
          </front>
          <body>
            <sec>
              <title>Results</title>
              <p>Variants are listed in Table 1.</p>
              <table-wrap id="t1">
                <label>Table 1</label>
                <caption><title>KCNH2 variants in affected probands</title></caption>
                <table>
                  <thead>
                    <tr><th>Variant</th><th>Affected</th><th>Unaffected</th></tr>
                  </thead>
                  <tbody>
                    <tr><td>p.Arg176Trp</td><td>4</td><td>1</td></tr>
                    <tr><td>p.Asn629Ser</td><td>2</td><td>0</td></tr>
                  </tbody>
                </table>
              </table-wrap>
            </sec>
          </body>
        </article>
        """)
    md = converter.xml_to_markdown(xml)
    assert "Table 1" in md
    assert "KCNH2 variants in affected probands" in md
    assert "| Variant | Affected | Unaffected |" in md
    assert "p.Arg176Trp" in md
    assert "p.Asn629Ser" in md


def test_pmc_html_to_markdown_preserves_article_table_rows(converter):
    html = textwrap.dedent("""\
        <html><body>
          <h1>KCNH2 cohort</h1>
          <section class="body main-article-body">
            <section>
              <h2>Results</h2>
              <p>Variants are listed in the table.</p>
              <table>
                <tr><th>Variant</th><th>Carriers</th></tr>
                <tr><td>p.Ala561Val</td><td>6</td></tr>
              </table>
            </section>
          </section>
        </body></html>
        """)
    md = converter.pmc_html_to_markdown(html)
    assert "p.Ala561Val" in md
    assert "| Variant | Carriers |" in md


def test_xml_supplement_to_markdown_strips_tags(converter, tmp_path: Path):
    xml = "<root><variant>p.G604S</variant><variant>p.A561V</variant></root>"
    p = tmp_path / "supp.xml"
    p.write_text(xml, encoding="utf-8")
    md = converter.xml_supplement_to_markdown(p)
    assert "G604S" in md
    assert "A561V" in md
    assert "<variant>" not in md


def test_xml_supplement_falls_back_for_malformed(converter, tmp_path: Path):
    p = tmp_path / "broken.xml"
    p.write_text("<root><not-closed>", encoding="utf-8")
    md = converter.xml_supplement_to_markdown(p)
    # Even when XML parsing fails, raw content should be preserved so
    # the variant scanner has a chance.
    assert "not-closed" in md or "root" in md


def test_extract_zip_supplement_processes_nested_files(converter, tmp_path: Path):
    zip_path = tmp_path / "bundle.zip"
    nested_csv = "variant,carriers\np.G604S,5\np.A561V,3\n"
    nested_txt = "Supplementary methods: detailed protocols."
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("table.csv", nested_csv)
        zf.writestr("methods.txt", nested_txt)

    extracted, md = converter.extract_zip_supplement(
        zip_path, dest_dir=tmp_path / "bundle"
    )
    assert len(extracted) == 2
    assert "G604S" in md
    assert "Supplementary methods" in md or "detailed protocols" in md
    # Both files should have been actually written to disk.
    assert all(p.exists() for p in extracted)


def test_extract_zip_supplement_rejects_path_traversal(converter, tmp_path: Path):
    zip_path = tmp_path / "evil.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("../escape.txt", "should not be written")
        zf.writestr("ok.txt", "fine")

    extracted, md = converter.extract_zip_supplement(
        zip_path, dest_dir=tmp_path / "out"
    )
    # Only the safe member should be extracted.
    names = [p.name for p in extracted]
    assert "ok.txt" in names
    assert not any(p.name == "escape.txt" for p in extracted)
    assert "fine" in md


def test_extract_zip_supplement_handles_bad_zip(converter, tmp_path: Path):
    p = tmp_path / "not-a-zip.zip"
    p.write_bytes(b"definitely not a zip file")
    extracted, md = converter.extract_zip_supplement(p, dest_dir=tmp_path / "out")
    assert extracted == []
    assert "Bad ZIP" in md or "ZIP" in md


def test_extract_zip_supplement_nested_pdf_figures_in_extracted(
    converter, tmp_path: Path
):
    """Nested PDF figure paths must appear in the returned extracted list."""
    zip_path = tmp_path / "bundle.zip"
    # Write a minimal valid PDF stub — real extraction is mocked below.
    fake_pdf_bytes = b"%PDF-1.4 stub"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("supp1.pdf", fake_pdf_bytes)

    figures_dir = tmp_path / "figs"
    fake_img_path = figures_dir / "fig_p1_0.png"

    fake_md = "### Page 1\n\nSome text.\n\n[Figure 1 from page 1: fig_p1_0.png (200x200px)]\n\n"
    fake_images = [
        {
            "page": 1,
            "index": 0,
            "width": 200,
            "height": 200,
            "size_bytes": 5000,
            "ext": "png",
            "filename": "fig_p1_0.png",
            "path": fake_img_path,
        }
    ]

    with patch.object(
        converter, "pdf_to_markdown_with_images", return_value=(fake_md, fake_images)
    ) as mock_pdf:
        extracted, md = converter.extract_zip_supplement(
            zip_path,
            dest_dir=tmp_path / "bundle",
            figures_dir=figures_dir,
            extract_images=True,
        )

    mock_pdf.assert_called_once()
    # The PDF itself and the extracted image path should both be in extracted.
    assert fake_img_path in extracted
    # Markdown should carry the figure placeholder from the nested PDF.
    assert "[Figure 1 from page 1" in md
    assert "supp1.pdf" in md or "fig_p1_0.png" in md


def test_pdf_to_markdown_with_images_creates_nested_output_dir(
    converter, tmp_path: Path
):
    """output_dir whose parent does not yet exist must be created without error."""
    deep_dir = tmp_path / "run" / "pmid_123" / "supplements"
    # deep_dir and its parents do NOT exist yet.
    assert not deep_dir.exists()

    # Provide a minimal valid PDF so the header check passes, then mock fitz.
    pdf_path = tmp_path / "test.pdf"
    pdf_path.write_bytes(b"%PDF-1.4\n%%EOF\n")

    fake_md = "### Page 1\n\nContent.\n\n"

    class FakePage:
        def get_text(self):
            return "Content."

        def get_images(self, full=True):
            return []

    class FakeDoc:
        def __iter__(self):
            return iter([FakePage()])

        def close(self):
            pass

    def fake_open(path):
        return FakeDoc()

    try:
        import fitz

        with patch.object(fitz, "open", side_effect=fake_open):
            md, images = converter.pdf_to_markdown_with_images(
                pdf_path, output_dir=deep_dir, extract_images=True
            )
        assert deep_dir.exists() or (deep_dir / "figures").exists()
    except ImportError:
        pytest.skip("PyMuPDF not installed")
