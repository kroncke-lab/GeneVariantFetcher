"""Regression tests migrated from the supplement parser's demo runner."""

from harvesting.supplement_reference_parser import (
    check_supplement_gap,
    extract_supplement_urls_from_text,
    parse_supplement_references,
)


KARGER_TEXT = """
We identified 54 KCNH2 mutations (online suppl. table 1,
www.karger.com/doi/10.1159/000440608 for all online material).
The mutation spectrum is shown in supplementary figure S1.
See also Table S2 for detailed phenotype data.
"""


def test_parse_supplement_references_migrated_karger_fixture():
    parsed = parse_supplement_references(KARGER_TEXT)

    assert parsed["expected_tables"] == 2
    assert parsed["expected_figures"] == 1
    assert set(parsed["table_refs"]) == {"Table S1", "Table S2"}
    assert parsed["figure_refs"] == ["Figure S1"]
    assert parsed["mentioned_variants"] == 54


def test_extract_supplement_urls_migrated_karger_fixture():
    assert extract_supplement_urls_from_text(KARGER_TEXT) == [
        "https://www.karger.com/doi/10.1159/000440608"
    ]


def test_check_supplement_gap_migrated_karger_fixture():
    gap = check_supplement_gap(
        KARGER_TEXT, downloaded_count=0, extracted_variant_count=4
    )

    assert gap["has_gap"] is True
    assert gap["expected_tables"] == 2
    assert gap["downloaded_count"] == 0
    assert gap["mentioned_variants"] == 54
    assert gap["extracted_variants"] == 4
    assert len(gap["warnings"]) >= 2
