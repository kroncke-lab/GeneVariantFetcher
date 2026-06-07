"""Tests for the adjudication-sheet row assembly (pipeline values + provenance)."""

from scripts.build_adjudication_sheet import (
    PIPELINE_COLUMNS,
    REVIEW_COLUMNS,
    _count_provenance,
    build_rows,
)


def _variant():
    return {
        "cdna_notation": "c.1A>T",
        "protein_notation": "p.Met1Leu",
        "clinical_significance": "pathogenic",
        "penetrance_data": {
            "total_carriers_observed": 5,
            "affected_count": 3,
            "unaffected_count": 2,
        },
        "patients": {"count": None},
        "count_provenance": {
            "carriers_column_label": "Carriers",
            "carriers_count_type": "explicit_total",
            "affected_count_type": "unknown",
            "unaffected_count_type": "unknown",
        },
        "source_location": "Table 2",
        "key_quotes": ["5 carriers, 3 affected by breast cancer"],
        "additional_notes": "from table 2",
    }


def test_build_rows_emits_variant_with_values_and_provenance():
    records = [
        {
            "pmid": "1",
            "title": "T",
            "url": "u",
            "fulltext_file": "papers/1.md",
            "abstract": "germline carriers were screened in families",
            "variants": [_variant()],
        }
    ]
    (row,) = build_rows(records)
    assert row["pipeline_variant"] == "c.1A>T / p.Met1Leu"
    assert row["carriers"] == 5 and row["affected"] == 3 and row["unaffected"] == 2
    assert "explicit_total" in row["count_provenance"]
    assert row["source_location"] == "Table 2"
    assert "5 carriers" in row["key_quote"]
    assert row["germline_or_somatic"] == "germline"  # from the abstract context
    # review columns exist but are blank for the assistant to fill
    assert all(col in row or col in REVIEW_COLUMNS for col in PIPELINE_COLUMNS)


def test_build_rows_marks_papers_with_no_extractions():
    records = [
        {
            "pmid": "2",
            "title": "T2",
            "url": "u2",
            "fulltext_file": "papers/2.md",
            "variants": [],
        }
    ]
    (row,) = build_rows(records)
    assert row["pipeline_variant"] == "(none extracted)"
    assert row["pmid"] == "2"


def test_count_provenance_skips_unknown_and_null():
    # all-unknown -> empty; a real count-type/label -> surfaced
    assert _count_provenance({"count_provenance": {}}) == ""
    v = {
        "count_provenance": {
            "carriers_count_type": "study_wide_total",
            "carriers_column_label": "N screened",
            "affected_count_type": "unknown",
        }
    }
    out = _count_provenance(v)
    assert "carriers: study_wide_total" in out and "N screened" in out
    assert "affected" not in out  # unknown is skipped
