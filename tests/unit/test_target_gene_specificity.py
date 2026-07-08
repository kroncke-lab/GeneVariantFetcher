import sqlite3

import pipeline.target_gene_specificity as target_gene_specificity
from utils.gene_metadata import clear_gene_metadata_cache

from pipeline.target_gene_specificity import (
    apply_target_gene_specificity,
    assess_variant_specificity,
)


def test_apoe_isoform_notation_is_supported_without_hgvs():
    variant = {"gene": "APOE", "protein_notation": "APOE ε4"}

    assessment = assess_variant_specificity(variant, gene_symbol="APOE")

    assert assessment["status"] == "supported"
    assert assessment["risk"] == "low"


def test_protein_position_beyond_gene_length_is_high_risk():
    variant = {"gene": "APOE", "protein_notation": "V717F"}

    assessment = assess_variant_specificity(variant, gene_symbol="APOE")

    assert assessment["status"] == "off_target_risk"
    assert assessment["risk"] == "high"


def test_variant_gene_alias_is_not_treated_as_off_target():
    variant = {"gene": "MYPBC3", "protein_notation": "p.Arg10Trp"}

    assessment = assess_variant_specificity(variant, gene_symbol="MYBPC3")

    assert assessment["status"] != "off_target_gene"


def test_variantfeatures_reference_mismatch_is_weak_qc_signal(tmp_path, monkeypatch):
    db_path = tmp_path / "variants.db"
    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        CREATE TABLE variant_consequences (
            gene_symbol TEXT,
            transcript_id TEXT,
            hgvs_p TEXT,
            hgvs_c TEXT,
            aa_pos INTEGER,
            aa_ref TEXT,
            aa_alt TEXT
        );
        INSERT INTO variant_consequences VALUES
            ('TESTGENE', 'ENST000TEST', 'ENSP000TEST:p.Gly10Val', 'ENST000TEST:c.29G>T', 10, 'G', 'V'),
            ('TESTGENE', 'ENST000TEST', 'ENSP000TEST:p.Ter101Ter', NULL, 101, '*', '*');
        """
    )
    conn.close()
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    clear_gene_metadata_cache()

    variant = {"gene": "TESTGENE", "protein_notation": "p.Ala10Val"}
    assessment = assess_variant_specificity(variant, gene_symbol="TESTGENE")

    assert assessment["status"] == "weak"
    assert assessment["risk"] == "medium"
    assert "VariantFeatures reference residue mismatch" in assessment["reasons"][0]
    assert assessment["variantfeatures"][0]["reference_residues"] == ["G"]
    assert assessment["protein_length"] == 100

    clear_gene_metadata_cache()


def test_apply_specificity_flags_but_keeps_rows_by_default():
    data = {
        "variants": [
            {"gene": "APOE", "protein_notation": "V717F"},
            {
                "gene": "APOE",
                "cdna_notation": "c.334T>C",
                "protein_notation": "p.Cys112Arg",
            },
        ]
    }

    stats = apply_target_gene_specificity(data, gene_symbol="APOE", policy="flag")

    assert stats["checked"] == 2
    assert stats["off_target_risk"] == 1
    assert len(data["variants"]) == 2
    assert data["variants"][0]["target_gene_specificity"]["risk"] == "high"
    assert data["extraction_metadata"]["target_gene_specificity"]["policy"] == "flag"


def test_apply_specificity_skips_context_scan_for_dense_sources(monkeypatch):
    monkeypatch.setattr(target_gene_specificity, "CONTEXT_SCAN_MAX_SOURCE_CHARS", 20)
    monkeypatch.setattr(target_gene_specificity, "CONTEXT_SCAN_MAX_VARIANTS", 2)
    monkeypatch.setattr(target_gene_specificity, "VARIANTFEATURES_MAX_VARIANTS", 2)
    data = {
        "variants": [
            {"gene": "BRCA1", "protein_notation": f"p.Arg{idx}Gly"}
            for idx in range(1, 4)
        ]
    }

    stats = apply_target_gene_specificity(
        data,
        gene_symbol="BRCA1",
        source_text="BRCA1 p.Arg1Gly " * 10,
        policy="flag",
    )

    assert stats["checked"] == 3
    metadata = data["extraction_metadata"]
    assert metadata["target_gene_specificity_context_scan"] == {
        "skipped": True,
        "reason": "source_or_variant_count_exceeds_context_scan_limit",
        "source_chars": 160,
        "variant_count": 3,
        "max_source_chars": 20,
        "max_variants": 2,
    }
    assert metadata["target_gene_specificity_variantfeatures"] == {
        "skipped": True,
        "reason": "variant_count_exceeds_variantfeatures_lookup_limit",
        "variant_count": 3,
        "max_variants": 2,
    }


def test_apply_specificity_clear_drops_high_risk_rows():
    data = {
        "variants": [
            {"gene": "APOE", "protein_notation": "V717F"},
            {"gene": "APOE", "protein_notation": "APOE ε4"},
        ]
    }

    stats = apply_target_gene_specificity(data, gene_symbol="APOE", policy="clear")

    assert stats["cleared"] == 1
    assert len(data["variants"]) == 1
    assert data["variants"][0]["protein_notation"] == "APOE ε4"
