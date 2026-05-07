"""Tests for new supplement-format handlers in format_converters.

Covers the additive surface added for the figure/supplement overhaul:
TSV, HTML supplement, XML supplement, and ZIP nested-extraction.
"""

from __future__ import annotations

import textwrap
import zipfile
from pathlib import Path

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
