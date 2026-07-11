"""Regression tests for notation validation and variant identity fallbacks."""

from __future__ import annotations

import pytest

from harvesting.migrate_to_sqlite import sanitize_variant_notation
from pipeline.extraction import ExpertExtractor


VALID_PROTEIN_NOTATIONS = (
    "p.Ala100_Glu101insLys",
    "p.Gly262Alafs*98",
    "p.Gly24fsTer58",
    "p.Arg176delinsLysGly",
    "p.Ala100_Glu101delinsLysTrp",
)

TRAILING_JUNK_PROTEIN_NOTATIONS = (
    "p.Ala100_Glu101insLys-junk",
    "p.Gly262Alafs*98junk",
    "p.Gly24fsTer58more",
    "p.Arg176delinsLysGlyoops",
)


@pytest.mark.parametrize("notation", VALID_PROTEIN_NOTATIONS)
def test_extractor_and_migration_accept_valid_protein_forms(notation: str):
    extractor = ExpertExtractor(models=["gpt-4"])
    row = {"gene_symbol": "KCNH2", "protein_notation": notation}

    filtered = extractor._filter_extraction_artifacts(
        {"extraction_metadata": {}, "variants": [row]},
        "KCNH2",
    )

    assert filtered["variants"] == [row]
    assert sanitize_variant_notation(dict(row)) is True


@pytest.mark.parametrize("notation", TRAILING_JUNK_PROTEIN_NOTATIONS)
def test_extractor_and_migration_reject_trailing_protein_junk(notation: str):
    extractor = ExpertExtractor(models=["gpt-4"])
    row = {"gene_symbol": "KCNH2", "protein_notation": notation}

    filtered = extractor._filter_extraction_artifacts(
        {"extraction_metadata": {}, "variants": [dict(row)]},
        "KCNH2",
    )
    migrated = dict(row)

    assert filtered["variants"] == []
    assert sanitize_variant_notation(migrated) is False
    assert migrated["protein_notation"] is None


def test_extractor_clears_bad_protein_but_keeps_valid_cdna_identity():
    extractor = ExpertExtractor(models=["gpt-4"])
    source = {
        "gene_symbol": "KCNH2",
        "protein_notation": "p.Arg176Trp-trailing-junk",
        "cdna_notation": "c.526C>T",
    }

    filtered = extractor._filter_extraction_artifacts(
        {"extraction_metadata": {}, "variants": [dict(source)]},
        "KCNH2",
    )
    migrated = dict(source)

    assert len(filtered["variants"]) == 1
    assert filtered["variants"][0]["protein_notation"] is None
    assert filtered["variants"][0]["cdna_notation"] == "c.526C>T"
    assert sanitize_variant_notation(migrated) is True
    assert migrated["protein_notation"] is None
    assert migrated["cdna_notation"] == "c.526C>T"


def test_table_merge_does_not_dedupe_point_variants_by_shared_description():
    extractor = ExpertExtractor(models=["gpt-4"])
    variants = [
        {
            "gene_symbol": "KCNH2",
            "cdna_notation": "c.100delAinsT",
            "structural_description": "complex indel",
        },
        {
            "gene_symbol": "KCNH2",
            "cdna_notation": "c.200delGinsC",
            "structural_description": "complex indel",
        },
    ]

    merged = extractor._merge_table_variants(
        {"extraction_metadata": {}, "variants": []}, variants
    )

    assert [row["cdna_notation"] for row in merged["variants"]] == [
        "c.100delAinsT",
        "c.200delGinsC",
    ]
