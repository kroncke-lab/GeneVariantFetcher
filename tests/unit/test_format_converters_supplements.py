"""Tests for new supplement-format handlers in format_converters.

Covers the additive surface added for the figure/supplement overhaul:
TSV, HTML supplement, XML supplement, and ZIP nested-extraction.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
import types
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


def test_dense_pdf_bypasses_markitdown_layout_sort(tmp_path: Path, monkeypatch):
    p = tmp_path / "dense_table.pdf"
    p.write_bytes(b"%PDF-1.7\n")

    class FakePage:
        def get_text(self, sort=False):
            return "p.Arg1Gly\tBMPR2\n" * 4_000

    class FakeDoc(list):
        def close(self):
            return None

    fake_fitz = types.SimpleNamespace(open=lambda _path: FakeDoc([FakePage()]))
    monkeypatch.setitem(sys.modules, "fitz", fake_fitz)
    converter = FormatConverter()
    converter.markitdown = types.SimpleNamespace(
        convert=lambda _path: pytest.fail("dense PDF reached MarkItDown")
    )

    md = converter.pdf_to_markdown(p)

    assert "p.Arg1Gly" in md
    assert md.startswith("### Page 1")


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


def test_xml_to_markdown_preserves_table_wrap_nested_inside_paragraph(converter):
    xml = textwrap.dedent("""\
        <article>
          <body>
            <sec>
              <title>Results</title>
              <p>Before the table.
                <table-wrap id="t2">
                  <label>Table 2</label>
                  <caption><title>Variants screened by gene</title></caption>
                  <table>
                    <thead><tr><th>Gene</th><th>Variant</th><th>Method</th></tr></thead>
                    <tbody>
                      <tr><td>BRCA1</td><td>c.5503C&gt;T</td><td>DHPLC</td></tr>
                      <tr><td>BRCA2</td><td>c.9117G&gt;A</td><td>DHPLC</td></tr>
                    </tbody>
                  </table>
                </table-wrap>
                After the table.
              </p>
            </sec>
          </body>
        </article>
        """)

    md = converter.xml_to_markdown(xml)

    assert "Before the table. After the table." in md
    assert "| Gene | Variant | Method |" in md
    assert "| BRCA1 | c.5503C>T | DHPLC |" in md
    assert "| BRCA2 | c.9117G>A | DHPLC |" in md
    assert md.count("### Table 2") == 1


def test_xml_to_markdown_does_not_duplicate_nested_table_in_body_container(converter):
    xml = textwrap.dedent("""\
        <article>
          <body>
            <statement>
              <p>Nested body container.
                <table-wrap id="t3">
                  <label>Table 3</label>
                  <table>
                    <tr><th>Gene</th><th>Variant</th></tr>
                    <tr><td>BMPR2</td><td>c.2695C&gt;T</td></tr>
                  </table>
                </table-wrap>
              </p>
            </statement>
          </body>
        </article>
        """)

    md = converter.xml_to_markdown(xml)

    assert md.count("### Table 3") == 1
    assert md.count("| BMPR2 | c.2695C>T |") == 1


def test_xml_to_markdown_uses_namespaced_main_body_not_subarticle(converter):
    xml = textwrap.dedent("""\
        <article xmlns="http://jats.nlm.nih.gov">
          <front><article-meta><title-group><article-title>Main cohort</article-title>
          </title-group></article-meta></front>
          <sub-article article-type="decision-letter">
            <body><sec><p>Reviewer decision text must not become source.</p></sec></body>
          </sub-article>
          <body><sec><title>Results</title>
            <p>BMPR2 c.2695C&gt;T was identified in the study cohort.</p>
          </sec></body>
        </article>
        """)

    md = converter.xml_to_markdown(xml)

    assert "Main cohort" in md
    assert "BMPR2 c.2695C>T was identified" in md
    assert "Reviewer decision text" not in md


def test_xml_to_markdown_preserves_back_matter_table_wrap_rows(converter):
    xml = textwrap.dedent("""\
        <article>
          <front>
            <article-meta>
              <title-group><article-title>SCN5A compendium</article-title></title-group>
            </article-meta>
          </front>
          <body>
            <sec>
              <title>Results</title>
              <p>The compendium is listed in Table 4.</p>
            </sec>
          </body>
          <back>
            <sec>
              <table-wrap id="T4">
                <label>Table 4</label>
                <caption><title>Compendium of SCN5A mutations</title></caption>
                <table>
                  <thead>
                    <tr><th>Exon</th><th>Nucleotide</th><th>Mutation</th><th>N</th></tr>
                  </thead>
                  <tbody>
                    <tr><td>Exon 28</td><td>5350 G&gt;A</td><td>E1784K</td><td>14</td></tr>
                    <tr><td>Exon 12</td><td>2582delT</td><td>F861WfsX90</td><td>11</td></tr>
                  </tbody>
                </table>
              </table-wrap>
            </sec>
          </back>
        </article>
        """)

    md = converter.xml_to_markdown(xml)

    assert "Table 4" in md
    assert "Compendium of SCN5A mutations" in md
    assert "| Exon | Nucleotide | Mutation | N |" in md
    assert "E1784K" in md
    assert "F861WfsX90" in md


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


def test_doc_to_markdown_uses_textutil_before_antiword_for_word_tables(
    converter, tmp_path: Path
):
    doc_path = tmp_path / "supp.doc"
    doc_path.write_bytes(b"legacy word placeholder")
    converter.markitdown = None

    textutil_text = (
        "Patient #\x07Gene\x07Aa change (p.)\x07Bases pair change (c.)"
        "\x07\x0713\x07SCN5A\x07Asn406Ser\x071217 A>G"
        "\x07\x0719\x07SCN5A\x07Glu1784Lys\x075350 G>A"
    )

    def fake_run(cmd, **kwargs):
        if cmd[0] == "soffice":
            raise FileNotFoundError
        if cmd[0] == "textutil":
            # The HTML leg fails here so the plain-text leg (the behavior this
            # test pins) is what produces the output.
            if "html" in cmd:
                return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="")
            return subprocess.CompletedProcess(cmd, 0, stdout=textutil_text, stderr="")
        raise AssertionError(f"unexpected converter command: {cmd}")

    with patch("subprocess.run", side_effect=fake_run):
        md = converter.doc_to_markdown(doc_path)

    assert "| Patient # | Gene | Aa change (p.) | Bases pair change (c.) |" in md
    assert "| 13 | SCN5A | Asn406Ser | 1217 A>G |" in md
    assert "| 19 | SCN5A | Glu1784Lys | 5350 G>A |" in md


def test_doc_to_markdown_prefers_textutil_html_for_table_cells(
    converter, tmp_path: Path
):
    """The HTML route must be attempted first and yield pipe-delimited cells.

    ``textutil -convert txt`` flattened KCNQ1 24667783's mutation table into
    delimiter-less runs ("…C>TMissensehetCM990761…"); ``-convert html``
    preserves <table> markup, so each cell survives as its own column.
    """
    doc_path = tmp_path / "supp.doc"
    doc_path.write_bytes(b"legacy word placeholder")
    converter.markitdown = None

    textutil_html = textwrap.dedent("""\
        <html><body>
          <p>Supplementary Table 1. KCNQ1 mutations.</p>
          <table>
            <tr><td>Nucleotide</td><td>Effect</td><td>Zygosity</td><td>HGMD</td></tr>
            <tr><td>477+5G>A</td><td>Splice</td><td>het</td><td>CS990751</td></tr>
            <tr><td>1032G>A</td><td>Missense</td><td>het</td><td>CM990761</td></tr>
          </table>
        </body></html>
        """)
    textutil_calls = []

    def fake_run(cmd, **kwargs):
        if cmd[0] == "soffice":
            raise FileNotFoundError
        if cmd[0] == "textutil":
            textutil_calls.append(list(cmd))
            if "html" in cmd:
                return subprocess.CompletedProcess(
                    cmd, 0, stdout=textutil_html, stderr=""
                )
            # Every route is scored, so the plain-text leg still runs; here it
            # recovers nothing, leaving the HTML rendering the only candidate.
            return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="")
        if cmd[0] == "antiword":
            return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="")
        raise AssertionError(f"unexpected converter command: {cmd}")

    with patch("subprocess.run", side_effect=fake_run):
        md = converter.doc_to_markdown(doc_path)

    assert textutil_calls[0][:3] == ["textutil", "-convert", "html"]
    assert "| Nucleotide | Effect | Zygosity | HGMD |" in md
    assert "| 477+5G>A | Splice | het | CS990751 |" in md
    assert "| 1032G>A | Missense | het | CM990761 |" in md
    # Non-table prose survives alongside the tables.
    assert "Supplementary Table 1. KCNQ1 mutations." in md
    # No flattened delimiter-less cell runs.
    assert "1032G>AMissensehet" not in md.replace(" ", "")  # sanity on join


def test_doc_html_route_does_not_touch_pdf_or_excel_paths(converter, tmp_path: Path):
    """Guard: the HTML-first conversion applies to .doc only."""
    pdf = tmp_path / "not_a_pdf.pdf"
    pdf.write_text("not a pdf", encoding="utf-8")

    def fail_run(cmd, **kwargs):
        raise AssertionError(f"pdf/xlsx conversion must not shell out: {cmd}")

    with patch("subprocess.run", side_effect=fail_run):
        md = converter.pdf_to_markdown(pdf)
    assert "[Invalid PDF file: not_a_pdf.pdf]" in md


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


def test_doc_to_markdown_prefers_antiword_when_textutil_flattens_tables(
    converter, tmp_path: Path
):
    """KCNQ1 24667783's real shape: textutil flattens the mutation table in
    EVERY output format (HTML keeps only a one-cell stub row), while
    ``antiword -t`` preserves the cells. A textutil rendering without table
    structure must not shadow an antiword rendering that has it.
    """
    doc_path = tmp_path / "supp.doc"
    doc_path.write_bytes(b"legacy word placeholder")
    converter.markitdown = None

    flattened_html = textwrap.dedent("""\
        <html><body>
          <table><tr><td></td><td>Clinical data</td></tr></table>
          <p>Gly306Arg916 G&gt;AMissensehetCM960900rs120074181N/ADELDEL</p>
        </body></html>
        """)
    antiword_tsv = (
        "SUPPLEMENTAL MATERIAL\n"
        "Patient\tGene\tAa change\tBases pair change\n"
        "37\tKCNQ1\tGly306Arg\t916 G>A\n"
        "41\tKCNQ1\tArg366Trp\t1096 C>T\n"
    )

    def fake_run(cmd, **kwargs):
        if cmd[0] == "soffice":
            raise FileNotFoundError
        if cmd[0] == "textutil":
            if "html" in cmd:
                return subprocess.CompletedProcess(
                    cmd, 0, stdout=flattened_html, stderr=""
                )
            # -convert txt: the flattened delimiter-less rendering.
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout="Gly306Arg916 G>AMissensehetCM960900rs120074181N/ADELDEL",
                stderr="",
            )
        if cmd[0] == "antiword":
            return subprocess.CompletedProcess(cmd, 0, stdout=antiword_tsv, stderr="")
        raise AssertionError(f"unexpected converter command: {cmd}")

    with patch("subprocess.run", side_effect=fake_run):
        md = converter.doc_to_markdown(doc_path)

    assert "| 37 | KCNQ1 | Gly306Arg | 916 G>A |" in md
    assert "| 41 | KCNQ1 | Arg366Trp | 1096 C>T |" in md


def test_doc_to_markdown_ignores_two_row_html_layout_table(converter, tmp_path: Path):
    """A Word HTML export's author/affiliation layout table is not a win.

    Two populated rows are what a title block looks like; the mutation table
    stays flattened in body text. Scoring must hand the file to the antiword
    rendering that carries the real rows instead of stopping at the first
    route that shows any table structure.
    """
    doc_path = tmp_path / "supp.doc"
    doc_path.write_bytes(b"legacy word placeholder")
    converter.markitdown = None

    layout_html = textwrap.dedent("""\
        <html><body>
          <table>
            <tr><td>Corresponding author</td><td>J. Doe, MD</td></tr>
            <tr><td>Affiliation</td><td>Dept. of Cardiology</td></tr>
          </table>
          <p>Supplementary Table 1. KCNQ1 variants.</p>
          <p>Gly306Arg916 G&gt;AMissensehetCM960900</p>
        </body></html>
        """)
    antiword_tsv = "Patient\tGene\tAa change\tBases pair change\n" + "".join(
        f"{n}\tKCNQ1\tArg{n}Trp\t{n}00 C>T\n" for n in range(1, 40)
    )

    def fake_run(cmd, **kwargs):
        if cmd[0] == "soffice":
            raise FileNotFoundError
        if cmd[0] == "textutil":
            if "html" in cmd:
                return subprocess.CompletedProcess(cmd, 0, stdout=layout_html)
            return subprocess.CompletedProcess(cmd, 1, stdout="")
        if cmd[0] == "antiword":
            return subprocess.CompletedProcess(cmd, 0, stdout=antiword_tsv)
        raise AssertionError(f"unexpected converter command: {cmd}")

    with patch("subprocess.run", side_effect=fake_run):
        md = converter.doc_to_markdown(doc_path)

    assert "| 39 | KCNQ1 | Arg39Trp | 3900 C>T |" in md
    assert "Corresponding author" not in md


def test_doc_to_markdown_keeps_narrative_over_two_row_antiword_table(
    converter, tmp_path: Path
):
    """A two-row antiword junk table must not cost the richer narrative.

    The inverse of the layout-table case: a thin table lead over a rendering
    with no tables at all is not decisive, so the longest usable text wins.
    """
    doc_path = tmp_path / "supp.doc"
    doc_path.write_bytes(b"legacy word placeholder")
    converter.markitdown = None

    narrative = (
        "Supplementary methods. Genomic DNA was extracted from peripheral "
        "blood leukocytes and all KCNQ1 exons were amplified. Carriers of "
        "Gly306Arg and Arg366Trp were genotyped in 42 relatives. " * 6
    )

    def fake_run(cmd, **kwargs):
        if cmd[0] == "soffice":
            raise FileNotFoundError
        if cmd[0] == "textutil":
            if "html" in cmd:
                return subprocess.CompletedProcess(cmd, 1, stdout="")
            return subprocess.CompletedProcess(cmd, 0, stdout=narrative)
        if cmd[0] == "antiword":
            return subprocess.CompletedProcess(
                cmd, 0, stdout="Page\tof\n1\t7\n", stderr=""
            )
        raise AssertionError(f"unexpected converter command: {cmd}")

    with patch("subprocess.run", side_effect=fake_run):
        md = converter.doc_to_markdown(doc_path)

    assert "Carriers of Gly306Arg and Arg366Trp were genotyped" in md
    assert "| 1 | 7 |" not in md


def test_doc_to_markdown_keeps_textutil_text_when_no_route_has_tables(
    converter, tmp_path: Path
):
    """A table-less .doc must keep the (richer) textutil rendering, not fall
    through to antiword's plain output."""
    doc_path = tmp_path / "notes.doc"
    doc_path.write_bytes(b"legacy word placeholder")
    converter.markitdown = None

    def fake_run(cmd, **kwargs):
        if cmd[0] == "soffice":
            raise FileNotFoundError
        if cmd[0] == "textutil":
            if "html" in cmd:
                return subprocess.CompletedProcess(
                    cmd,
                    0,
                    stdout="<html><body><p>Supplementary methods prose only.</p></body></html>",
                    stderr="",
                )
            return subprocess.CompletedProcess(
                cmd, 0, stdout="Supplementary methods prose only.", stderr=""
            )
        if cmd[0] == "antiword":
            return subprocess.CompletedProcess(
                cmd, 0, stdout="antiword prose rendering", stderr=""
            )
        raise AssertionError(f"unexpected converter command: {cmd}")

    with patch("subprocess.run", side_effect=fake_run):
        md = converter.doc_to_markdown(doc_path)

    assert "Supplementary methods prose only." in md
    assert "antiword prose rendering" not in md


# ---------------------------------------------------------------------------
# Binary-garbage detection on .doc conversion routes (RYR2 PMID 32508047:
# MarkItDown "converted" a 13.6 MB legacy .doc by decoding the raw OLE
# container as text; textutil does the same on that file).
# ---------------------------------------------------------------------------

from harvesting.format_converters import looks_like_binary_garbage  # noqa: E402


def _soup(n: int) -> str:
    unit = "\x00\x01g\x00\x00\x0b|ˇˇ\x02\x03\x04 –œ\x11‡°±\x1a·"
    return (unit * (n // len(unit) + 1))[:n]


def test_binary_garbage_detector_flags_soup_and_spares_real_text():
    assert looks_like_binary_garbage(_soup(50_000)) is True
    # Heterogeneous leak: a text-like first chunk must not launder a
    # binary-dense later chunk (the real soup's head 300 KB looked like text).
    hetero = ("real words about RYR2 carriers. " * 40_000) + _soup(1_200_000)
    assert looks_like_binary_garbage(hetero) is True

    prose = (
        "The proband carried RYR2 p.Arg420Trp. Family members were "
        "genotyped and phenotyped for CPVT.\n"
    ) * 100
    assert looks_like_binary_garbage(prose) is False
    pipe_table = "| c.1259G>A | 5 | 3 |\n|---|---|---|\n" * 200
    assert looks_like_binary_garbage(pipe_table) is False
    dna = "ACGTGGCTAAGCTTGACGTA" * 500
    assert looks_like_binary_garbage(dna) is False
    fixed_width = "R420W      12     3      0.87\n" * 200
    assert looks_like_binary_garbage(fixed_width) is False
    chinese = "基因变异与心律失常的相关性研究。患者携带罕见变异。" * 100
    assert looks_like_binary_garbage(chinese) is False
    # Sparse legitimate control chars (form feeds from PDF page breaks).
    paged = ("page text " * 200 + "\x0c") * 20
    assert looks_like_binary_garbage(paged) is False
    # Too short for a ratio verdict.
    assert looks_like_binary_garbage("\x00\x01\x02") is False
    assert looks_like_binary_garbage("") is False


def test_doc_to_markdown_rejects_markitdown_binary_leak(tmp_path, monkeypatch):
    import types

    # A junk .doc that no route can parse: markitdown "succeeds" with soup,
    # every subprocess tool is absent, OLE parsing fails, and the heuristic
    # byte scrape yields too little text — the chain must end at the
    # placeholder rather than return the soup.
    doc = tmp_path / "leaky.doc"
    doc.write_bytes(b"\x00\x01\x02junk" * 20)

    soup = _soup(20_000)
    converter = FormatConverter()
    converter.markitdown = types.SimpleNamespace(
        convert=lambda path: types.SimpleNamespace(text_content=soup)
    )

    def _no_tools(*args, **kwargs):
        raise FileNotFoundError("tool not installed")

    monkeypatch.setattr(subprocess, "run", _no_tools)

    md = converter.doc_to_markdown(doc)

    assert "\x00" not in md
    assert soup[:50] not in md
    assert "text extraction failed" in md  # the honest placeholder


def test_doc_to_markdown_keeps_good_markitdown_text(tmp_path):
    import types

    doc = tmp_path / "fine.doc"
    doc.write_bytes(b"irrelevant")
    good = "Table S5: RYR2 rare variants\n| c.1259G>A | 5 |\n" * 50
    converter = FormatConverter()
    converter.markitdown = types.SimpleNamespace(
        convert=lambda path: types.SimpleNamespace(text_content=good)
    )

    assert converter.doc_to_markdown(doc) == good
