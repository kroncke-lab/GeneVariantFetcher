"""Tests for cheap paper-level extraction census estimates."""

from pipeline.extraction import ExpertExtractor


def _table(rows: int = 8) -> str:
    body = "\n".join(
        f"| p.Arg{idx}Trp | {idx + 1} | {idx} | 1 |" for idx in range(1, rows + 1)
    )
    return f"""
Table 1. KCNH2 variant carriers

| Protein | Carriers | Affected | Unaffected |
|---|---:|---:|---:|
{body}
"""


def test_paper_census_estimates_variant_and_count_ranges():
    extractor = ExpertExtractor(models=["noop"], tier_threshold=0)

    census = extractor._estimate_paper_census(
        _table(rows=4),
        "KCNH2",
        table_row_hint=6,
        scanner_variant_count=4,
        table_hint_variant_count=4,
    )

    assert census["version"] == "deterministic_v1"
    assert census["estimated_variant_rows"]["point"] == 6
    assert census["estimated_unique_variants"]["point"] >= 4
    assert census["estimated_carriers"]["point"] == 14
    assert census["estimated_affected"]["point"] == 10
    assert census["estimated_unaffected"]["point"] == 4
    assert "count_columns" in census["basis"]


def test_paper_census_drives_low_yield_risk_without_replacing_facts():
    extractor = ExpertExtractor(models=["noop"], tier_threshold=0)
    census = extractor._estimate_paper_census(
        _table(rows=10),
        "KCNH2",
        table_row_hint=12,
        scanner_variant_count=10,
        table_hint_variant_count=10,
    )
    extracted = {
        "variants": [],
        "extraction_metadata": {"paper_census": census},
    }

    risk = extractor._assess_extraction_risk(
        extracted_data=extracted,
        source_text=_table(rows=10),
        estimated_variants=12,
        scanner_variant_count=10,
        paper_census=census,
    )

    assert risk["requires_adjudication"] is True
    assert any(
        reason.startswith("low_variant_yield_vs_census") for reason in risk["reasons"]
    )
    assert extracted["variants"] == []
