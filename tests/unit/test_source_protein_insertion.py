"""Keep explicit insertion identity/counts without inventing a missing residue."""

import copy

import pytest

from harvesting.migrate_to_sqlite import sanitize_variant_notation
from pipeline.extraction import ExpertExtractor
from utils.protein_notation import normalize_partial_protein_insertion


@pytest.mark.parametrize(
    "source,expected",
    [
        ("p.1570-Phe1571insIle", "p.1570_F1571insI"),
        ("p.100_Glu101insLysGly", "p.100_E101insKG"),
        ("p.1_A2insC", "p.1_A2insC"),
    ],
)
def test_extraction_and_migration_preserve_source_and_observations(source, expected):
    row = {
        "gene_symbol": "SCN5A",
        "protein_notation": source,
        "source_notation": source,
        "patients": {"count": 11},
        "penetrance_data": {
            "affected_count": 6,
            "unaffected_count": 4,
            "uncertain_count": 1,
        },
    }
    x = ExpertExtractor.__new__(ExpertExtractor)
    filtered = x._filter_extraction_artifacts(
        {"variants": [copy.deepcopy(row)]}, "SCN5A"
    )
    migrated = copy.deepcopy(row)
    assert sanitize_variant_notation(migrated)
    assert filtered["variants"] == [migrated]
    assert migrated["protein_notation"] == expected
    assert migrated["source_notation"] == source
    assert migrated["patients"] == row["patients"]
    assert migrated["penetrance_data"] == row["penetrance_data"]
    assert normalize_partial_protein_insertion(expected) == expected
    assert not migrated.get("cdna_notation")


@pytest.mark.parametrize(
    "source",
    [
        "p.1570-Phe1572insIle",
        "p.1571-Phe1570insIle",
        "p.0-Ala1insIle",
        "p.1570-Phe1571insIle-junk",
        "p.1570-Phe1571delIle",
        "p.1570-1571insIle",
        "p.1570-Phe1571insStop",
        "p.1570-Xaa1571insIle",
        "c.1570-Phe1571insIle",
        "p.1570-Phe1571ins",
    ],
)
def test_ambiguous_or_malformed_source_is_not_promoted(source):
    assert normalize_partial_protein_insertion(source) is None
    row = {"gene_symbol": "SCN5A", "protein_notation": source}
    assert (
        ExpertExtractor.__new__(ExpertExtractor)._filter_extraction_artifacts(
            {"variants": [copy.deepcopy(row)]}, "SCN5A"
        )["variants"]
        == []
    )
    assert not sanitize_variant_notation(row)


def test_position_guard_still_rejects_impossible_insertion():
    row = {"gene_symbol": "SCN5A", "protein_notation": "p.9000-Phe9001insIle"}
    assert (
        ExpertExtractor.__new__(ExpertExtractor)._filter_extraction_artifacts(
            {"variants": [row]}, "SCN5A"
        )["variants"]
        == []
    )
