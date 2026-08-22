"""Tests for SQLite migration notation guards."""

import json

from harvesting.migrate_to_sqlite import (
    _first_quote,
    create_database_schema,
    find_extraction_json_files,
    insert_paper_metadata,
    insert_variant_data,
    migrate_extraction_directory,
    normalize_affected_status,
    repair_extraction_data,
    validate_extraction_data,
)
from utils.pmid_utils import extract_pmid_from_filename


def test_first_quote_empty_list_is_none_not_bracket_literal():
    # Regression: an empty key_quotes list used to serialize to the literal "[]"
    # and land in fact_provenance.evidence_quote for every quote-less fact.
    assert _first_quote({"key_quotes": []}) is None
    assert _first_quote({"key_quotes": [""]}) is None
    assert _first_quote({}) is None
    assert _first_quote({"key_quotes": ["real quote"]}) == "real quote"


def test_migration_never_stores_bracket_literal_quote(tmp_path):
    db = tmp_path / "G.db"
    conn = create_database_schema(str(db))
    cur = conn.cursor()
    extraction = {
        "paper_metadata": {"pmid": "1", "title": "t"},
        "variants": [
            {"gene_symbol": "KCNH2", "protein_notation": "A561V", "key_quotes": []}
        ],
    }
    insert_paper_metadata(cur, extraction, replace_existing=False)
    insert_variant_data(cur, "1", extraction["variants"][0])
    conn.commit()
    junk = cur.execute(
        "SELECT COUNT(*) FROM fact_provenance WHERE evidence_quote = '[]'"
    ).fetchone()[0]
    assert junk == 0
    conn.close()


def test_migration_skips_malformed_protein_only_variant(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()

    extraction_data = {
        "paper_metadata": {"pmid": "32386560", "title": "GWAS table"},
        "variants": [{"gene_symbol": "KCNH2"}],
    }
    insert_paper_metadata(cursor, extraction_data)

    variant_id = insert_variant_data(
        cursor,
        "32386560",
        {
            "gene_symbol": "KCNH2",
            "protein_notation": "A",
            "penetrance_data": {
                "total_carriers_observed": 29762,
                "affected_count": 29762,
            },
        },
    )
    conn.commit()

    assert variant_id is None
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 0
    assert cursor.execute("SELECT COUNT(*) FROM penetrance_data").fetchone()[0] == 0
    conn.close()


def test_migration_keeps_valid_variant_notation(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()

    extraction_data = {
        "paper_metadata": {"pmid": "12345", "title": "Clinical table"},
        "variants": [{"gene_symbol": "KCNH2"}],
    }
    insert_paper_metadata(cursor, extraction_data)

    variant_id = insert_variant_data(
        cursor,
        "12345",
        {
            "gene_symbol": "KCNH2",
            "protein_notation": "p.Lys897Thr",
            "penetrance_data": {
                "total_carriers_observed": 7,
                "affected_count": 3,
                "unaffected_count": 4,
            },
        },
    )
    conn.commit()

    assert variant_id is not None
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    assert cursor.execute("SELECT COUNT(*) FROM penetrance_data").fetchone()[0] == 1
    conn.close()


def test_migration_folds_sparse_cdna_alias_into_richer_identity_and_count(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    extraction_data = {
        "paper_metadata": {"pmid": "27767231", "title": "BRCA1 tables"},
        "variants": [],
    }
    insert_paper_metadata(cursor, extraction_data)

    sparse_id = insert_variant_data(
        cursor,
        "27767231",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.121del",
            "penetrance_data": {"total_carriers_observed": 1},
            "source_table": "Supplementary validation table",
        },
    )
    rich_id = insert_variant_data(
        cursor,
        "27767231",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.121del",
            "protein_notation": "p.His41fs",
            "penetrance_data": {
                "total_carriers_observed": 1,
                "affected_count": 1,
            },
            "source_table": "Table 2",
        },
    )
    conn.commit()

    assert rich_id == sparse_id
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    assert cursor.execute(
        "SELECT cdna_notation, protein_notation FROM variants"
    ).fetchone() == ("c.121del", "p.His41fs")
    assert cursor.execute(
        """
        SELECT total_carriers_observed, affected_count, unaffected_count
        FROM penetrance_data
        """
    ).fetchall() == [(1, 1, None)]
    conn.close()


def test_migration_folds_sparse_alias_after_richer_identity(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "29470806", "title": "BRCA1 tables"}},
    )

    rich_id = insert_variant_data(
        cursor,
        "29470806",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.1002delC",
            "protein_notation": "p.Ser335AlafsTer6",
        },
    )
    sparse_id = insert_variant_data(
        cursor,
        "29470806",
        {"gene_symbol": "BRCA1", "cdna_notation": "c.1002delC"},
    )
    conn.commit()

    assert sparse_id == rich_id
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_migration_keeps_conflicting_protein_annotations_for_same_cdna(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "Conflicting annotations"}},
    )

    first_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.100C>T",
            "protein_notation": "p.Arg34Trp",
        },
    )
    second_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.100C>T",
            "protein_notation": "p.Arg34Gln",
        },
    )
    conn.commit()

    assert second_id != first_id
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 2
    conn.close()


def test_migration_folds_truncated_substitution_into_matching_frameshift(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "BRCA1 frameshift table"}},
    )

    short_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.1016dupA",
            "protein_notation": "V340G",
        },
    )
    frameshift_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.1016dupA",
            "protein_notation": "p.Val340Glyfs*6",
        },
    )
    conn.commit()

    assert frameshift_id == short_id
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    assert cursor.execute("SELECT protein_notation FROM variants").fetchone()[0] == (
        "p.Val340Glyfs*6"
    )
    conn.close()


def test_migration_folds_matching_frameshift_alias_in_reverse_order(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "BRCA1 frameshift table"}},
    )

    frameshift_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.1016dupA",
            "protein_notation": "p.Val340Glyfs*6",
        },
    )
    short_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.1016dupA",
            "protein_notation": "V340G",
        },
    )
    conn.commit()

    assert short_id == frameshift_id
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    assert cursor.execute("SELECT protein_notation FROM variants").fetchone()[0] == (
        "p.Val340Glyfs*6"
    )
    conn.close()


def test_migration_does_not_fold_unproven_frameshift_aliases(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "Conflicting annotations"}},
    )

    ids = []
    for cdna, proteins in (
        ("c.100dupA", ("V34G", "p.Val34Alafs*6")),
        ("c.101dupA", ("V34G", "p.Val34fs*6")),
    ):
        for protein in proteins:
            ids.append(
                insert_variant_data(
                    cursor,
                    "1",
                    {
                        "gene_symbol": "BRCA1",
                        "cdna_notation": cdna,
                        "protein_notation": protein,
                    },
                )
            )
    conn.commit()

    assert len(set(ids)) == 4
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 4
    conn.close()


def test_migration_keeps_protein_only_then_rich_identity_separate(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "Protein-only alias"}},
    )

    sparse_id = insert_variant_data(
        cursor,
        "1",
        {"gene_symbol": "BRCA1", "protein_notation": "p.Arg34Trp"},
    )
    rich_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.100C>T",
            "protein_notation": "R34W",
        },
    )
    conn.commit()

    assert rich_id != sparse_id
    assert cursor.execute(
        "SELECT cdna_notation, protein_notation FROM variants ORDER BY variant_id"
    ).fetchall() == [(None, "p.Arg34Trp"), ("c.100C>T", "R34W")]
    conn.close()


def test_migration_keeps_rich_then_protein_only_identity_separate(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "Protein-only alias"}},
    )

    rich_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.100C>T",
            "protein_notation": "R34W",
        },
    )
    sparse_id = insert_variant_data(
        cursor,
        "1",
        {"gene_symbol": "BRCA1", "protein_notation": "p.Arg34Trp"},
    )
    conn.commit()

    assert sparse_id != rich_id
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 2
    conn.close()


def test_migration_unwraps_predicted_frameshift_parentheses(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "Predicted frameshift"}},
    )

    variant_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.1016dupA",
            "protein_notation": "p.(Val340Glyfs*6)",
        },
    )
    conn.commit()

    assert variant_id is not None
    assert cursor.execute("SELECT protein_notation FROM variants").fetchone()[0] == (
        "p.Val340Glyfs*6"
    )
    conn.close()


def test_migration_folds_equivalent_protein_spellings_for_same_cdna(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "Equivalent annotations"}},
    )

    first_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.220C>T",
            "protein_notation": "p.Gln74*",
        },
    )
    second_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.220C>T",
            "protein_notation": "p.Gln74Ter",
        },
    )
    conn.commit()

    assert second_id == first_id
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_migration_keeps_distinct_cdna_variants_with_same_protein(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "Codon variants"}},
    )

    first_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.100C>T",
            "protein_notation": "p.Arg34Trp",
        },
    )
    second_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.101G>A",
            "protein_notation": "p.Arg34Trp",
        },
    )
    conn.commit()

    assert second_id != first_id
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 2
    conn.close()


def test_migration_keeps_conflicting_count_strata_on_folded_alias(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "Distinct cohorts"}},
    )

    rich_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.100C>T",
            "protein_notation": "p.Arg34Trp",
            "penetrance_data": {"total_carriers_observed": 2},
        },
    )
    sparse_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.100C>T",
            "penetrance_data": {"total_carriers_observed": 3},
        },
    )
    conn.commit()

    assert sparse_id == rich_id
    assert cursor.execute(
        "SELECT total_carriers_observed FROM penetrance_data ORDER BY penetrance_id"
    ).fetchall() == [(2,), (3,)]
    conn.close()


def test_migration_does_not_merge_arithmetically_impossible_partial_strata(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()
    insert_paper_metadata(
        cursor,
        {"paper_metadata": {"pmid": "1", "title": "Distinct count strata"}},
    )

    first_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.100C>T",
            "penetrance_data": {"total_carriers_observed": 2},
        },
    )
    second_id = insert_variant_data(
        cursor,
        "1",
        {
            "gene_symbol": "BRCA1",
            "cdna_notation": "c.100C>T",
            "protein_notation": "p.Arg34Trp",
            "penetrance_data": {
                "total_carriers_observed": 2,
                "affected_count": 50,
            },
        },
    )
    conn.commit()

    assert second_id == first_id
    assert cursor.execute(
        "SELECT total_carriers_observed, affected_count FROM penetrance_data "
        "ORDER BY penetrance_id"
    ).fetchall() == [(2, None), (2, 50)]
    conn.close()


def test_migration_keeps_intronic_range_cdna_indel(tmp_path):
    """A range endpoint may itself carry an intronic offset."""
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()

    extraction_data = {
        "paper_metadata": {"pmid": "19949876", "title": "BRCA2 table"},
        "variants": [{"gene_symbol": "BRCA2"}],
    }
    insert_paper_metadata(cursor, extraction_data)

    variant_id = insert_variant_data(
        cursor,
        "19949876",
        {
            "gene_symbol": "BRCA2",
            "cdna_notation": "c.7436_7437-2delAGAT",
        },
    )
    conn.commit()

    assert variant_id is not None
    assert (
        cursor.execute(
            "SELECT cdna_notation FROM variants WHERE variant_id = ?", (variant_id,)
        ).fetchone()[0]
        == "c.7436_7437-2delAGAT"
    )
    conn.close()


def test_validation_accepts_structural_variant_identity_without_notation():
    payload = {
        "paper_metadata": {"pmid": "16429403", "title": "BMPR2 deletions"},
        "variants": [
            {
                "gene_symbol": "BMPR2",
                "variant_class": "exon_deletion",
                "structural_description": "deletion of exons 5-7",
            }
        ],
    }

    valid, errors, warnings = validate_extraction_data(
        payload, "BMPR2_PMID_16429403.json"
    )

    assert valid
    assert errors == []
    assert not any("identity" in warning for warning in warnings)


def test_migration_keeps_uncertain_start_loss_notation(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()

    extraction_data = {
        "paper_metadata": {"pmid": "32893267", "title": "Supplemental table"},
        "variants": [{"gene_symbol": "KCNQ1"}],
    }
    insert_paper_metadata(cursor, extraction_data)

    variant_id = insert_variant_data(
        cursor,
        "32893267",
        {
            "gene_symbol": "KCNQ1",
            "protein_notation": "p.Met1?",
        },
    )
    conn.commit()

    assert variant_id is not None
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_migration_keeps_multi_residue_deletion(tmp_path):
    """Range HGVS like p.Asp2_Arg135del must survive sanitization.

    Regression test: the original PROTEIN_NOTATION_RE rejected these because
    it required the second AA-or-tag immediately after the first position
    digits, with no slot for `_<AAA><digits>`. Observed in PMID 28122216
    during 2026-05-09 cascade run.
    """
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()

    extraction_data = {
        "paper_metadata": {"pmid": "28122216", "title": "Engineered deletion"},
        "variants": [{"gene_symbol": "KCNH2"}],
    }
    insert_paper_metadata(cursor, extraction_data)

    variant_id = insert_variant_data(
        cursor,
        "28122216",
        {
            "gene_symbol": "KCNH2",
            "protein_notation": "p.Asp2_Arg135del",
        },
    )
    conn.commit()

    assert variant_id is not None
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_migration_keeps_hyphen_range_deletions_from_lqt_tables(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    cursor = conn.cursor()

    extraction_data = {
        "paper_metadata": {"pmid": "30758498", "title": "LQT1 supplement"},
        "variants": [{"gene_symbol": "KCNQ1"}],
    }
    insert_paper_metadata(cursor, extraction_data)

    first_id = insert_variant_data(
        cursor,
        "30758498",
        {"gene_symbol": "KCNQ1", "protein_notation": "p.A178-G189del"},
    )
    second_id = insert_variant_data(
        cursor,
        "30758498",
        {"gene_symbol": "KCNQ1", "protein_notation": "p.Q521-Y522delT"},
    )
    conn.commit()

    assert first_id is not None
    assert second_id is not None
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 2
    conn.close()


def test_migration_discovers_extract_folder_json_names(tmp_path):
    extraction_dir = tmp_path / "extractions"
    extraction_dir.mkdir()
    payload = {
        "paper_metadata": {"pmid": "12345678", "title": "Clinical table"},
        "variants": [
            {
                "gene_symbol": "KCNH2",
                "protein_notation": "p.Lys897Thr",
            }
        ],
    }
    (extraction_dir / "12345678_extraction.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )

    assert [p.name for p in find_extraction_json_files(extraction_dir)] == [
        "12345678_extraction.json"
    ]

    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    stats = migrate_extraction_directory(conn, extraction_dir)

    assert stats["successful"] == 1
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_migration_discovery_ignores_timestamped_backup_json(tmp_path):
    extraction_dir = tmp_path / "extractions"
    extraction_dir.mkdir()
    payload = {
        "paper_metadata": {"pmid": "26496715", "title": "Current extraction"},
        "variants": [{"gene_symbol": "KCNH2", "protein_notation": "p.Ala561Val"}],
    }
    for name in [
        "KCNH2_PMID_26496715.json",
        "KCNH2_PMID_26496715.20260514_figshare_pre.json",
        "KCNH2_PMID_26496715.backup.json",
        "notes.json",
    ]:
        (extraction_dir / name).write_text(json.dumps(payload), encoding="utf-8")

    assert [p.name for p in find_extraction_json_files(extraction_dir)] == [
        "KCNH2_PMID_26496715.json"
    ]


def test_extract_pmid_from_gene_pmid_json_filename():
    assert extract_pmid_from_filename("SCN5A_PMID_32533946.json") == "32533946"


def test_repair_promotes_top_level_metadata_and_normalizes_records():
    payload = {
        "pmid": "26705554",
        "title": "Clinical EP study",
        "variants": [
            {
                "gene_symbol": "RYR2",
                "protein_notation": "p.Pro2328Ser",
                "individual_records": [
                    {
                        "individual_id": "iPSC_donor",
                        "affected_status": "carrier",
                    },
                    "non-object row",
                ],
            }
        ],
    }

    repaired, repairs = repair_extraction_data(payload, "RYR2_PMID_26705554.json")

    assert repaired["paper_metadata"]["pmid"] == "26705554"
    assert repaired["paper_metadata"]["title"] == "Clinical EP study"
    assert repaired["variants"][0]["individual_records"] == [
        {"individual_id": "iPSC_donor", "affected_status": "uncertain"}
    ]
    assert any("dropped non-object" in repair for repair in repairs)


def test_normalize_affected_status_enum_values():
    assert normalize_affected_status("symptomatic") == "affected"
    assert normalize_affected_status("asymptomatic") == "unaffected"
    assert normalize_affected_status("carrier") == "uncertain"


def test_migration_normalizes_individual_record_status(tmp_path):
    extraction_dir = tmp_path / "extractions"
    extraction_dir.mkdir()
    payload = {
        "pmid": "26705554",
        "title": "Clinical EP study",
        "variants": [
            {
                "gene_symbol": "RYR2",
                "protein_notation": "p.Pro2328Ser",
                "individual_records": [
                    {
                        "individual_id": "iPSC_donor",
                        "affected_status": "carrier",
                    }
                ],
            }
        ],
    }
    (extraction_dir / "RYR2_PMID_26705554.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )

    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    stats = migrate_extraction_directory(conn, extraction_dir)

    assert stats["successful"] == 1
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    assert (
        conn.execute("SELECT affected_status FROM individual_records").fetchone()[0]
        == "uncertain"
    )
    conn.close()


def test_migration_repairs_pmcid_in_pmid_field_from_filename(tmp_path):
    extraction_dir = tmp_path / "extractions"
    extraction_dir.mkdir()
    payload = {
        "paper_metadata": {
            "pmid": "PMC9522753",
            "title": "Inherited channelopathies update",
        },
        "variants": [
            {
                "gene_symbol": "KCNH2",
                "protein_notation": "p.Lys897Thr",
                "penetrance_data": {
                    "total_carriers_observed": 8,
                    "affected_count": 5,
                    "unaffected_count": 3,
                },
            }
        ],
    }
    (extraction_dir / "KCNH2_PMID_34546463.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )

    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    stats = migrate_extraction_directory(conn, extraction_dir)

    assert stats["successful"] == 1
    assert conn.execute("SELECT pmid, pmc_id FROM papers").fetchone() == (
        "34546463",
        "PMC9522753",
    )
    assert conn.execute("SELECT DISTINCT pmid FROM penetrance_data").fetchone() == (
        "34546463",
    )
    conn.close()


def test_repair_keeps_different_valid_metadata_pmid():
    payload = {
        "paper_metadata": {"pmid": "12345678", "title": "Correct metadata"},
        "variants": [],
    }

    valid, errors, warnings = validate_extraction_data(
        payload, "KCNH2_PMID_34546463.json"
    )
    repaired, repairs = repair_extraction_data(payload, "KCNH2_PMID_34546463.json")

    assert valid
    assert not errors
    assert any("keeping metadata" in warning for warning in warnings)
    assert repaired["paper_metadata"]["pmid"] == "12345678"
    assert not any("Set pmid from filename" in repair for repair in repairs)


def test_validation_rejects_invalid_pmid_without_filename_fallback():
    payload = {
        "paper_metadata": {
            "pmid": "PMC9522753",
            "title": "Legacy PMCID metadata",
        },
        "variants": [],
    }

    valid, errors, _ = validate_extraction_data(payload, "legacy_extraction.json")

    assert not valid
    assert errors == ["paper_metadata.pmid is not a valid PMID: PMC9522753"]


def test_migration_normalizes_valid_pmid_whitespace(tmp_path):
    extraction_dir = tmp_path / "extractions"
    extraction_dir.mkdir()
    payload = {
        "paper_metadata": {
            "pmid": " 34546463 ",
            "title": "Whitespace-padded PMID",
        },
        "variants": [],
    }
    (extraction_dir / "KCNH2_PMID_34546463.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )

    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))
    stats = migrate_extraction_directory(conn, extraction_dir)

    assert stats["successful"] == 1
    assert conn.execute("SELECT pmid FROM papers").fetchone() == ("34546463",)
    conn.close()


def test_repair_still_recovers_missing_and_unknown_pmids_from_filename():
    for bad_pmid in (None, "", "UNKNOWN"):
        payload = {
            "paper_metadata": {"pmid": bad_pmid, "title": "Incomplete metadata"},
            "variants": [],
        }

        repaired, repairs = repair_extraction_data(payload, "KCNH2_PMID_34546463.json")

        assert repaired["paper_metadata"]["pmid"] == "34546463"
        assert any("Set pmid from filename" in repair for repair in repairs)


def test_repair_does_not_overwrite_conflicting_existing_pmc_id():
    payload = {
        "paper_metadata": {
            "pmid": "PMC9522753",
            "pmc_id": "PMC1111111",
            "title": "Conflicting legacy metadata",
        },
        "variants": [],
    }

    repaired, repairs = repair_extraction_data(payload, "KCNH2_PMID_34546463.json")

    assert repaired["paper_metadata"]["pmid"] == "34546463"
    assert repaired["paper_metadata"]["pmc_id"] == "PMC1111111"
    assert any("discarded conflicting" in repair for repair in repairs)


def test_targeted_migration_preserves_existing_pmid_evidence(tmp_path):
    db_path = tmp_path / "variants.db"
    conn = create_database_schema(str(db_path))

    initial_dir = tmp_path / "initial"
    initial_dir.mkdir()
    initial_payload = {
        "paper_metadata": {"pmid": "29925740", "title": "RYR2 baseline"},
        "variants": [
            {
                "gene_symbol": "RYR2",
                "protein_notation": "p.Ala100Val",
                "penetrance_data": {
                    "total_carriers_observed": 3,
                    "affected_count": 2,
                    "unaffected_count": 1,
                },
                "individual_records": [
                    {
                        "individual_id": "P1",
                        "affected_status": "affected",
                    }
                ],
            }
        ],
    }
    (initial_dir / "RYR2_PMID_29925740.json").write_text(
        json.dumps(initial_payload), encoding="utf-8"
    )

    assert migrate_extraction_directory(conn, initial_dir)["successful"] == 1

    targeted_dir = tmp_path / "targeted"
    targeted_dir.mkdir()
    targeted_payload = {
        "paper_metadata": {"pmid": "29925740", "title": "RYR2 targeted"},
        "variants": [
            {
                "gene_symbol": "RYR2",
                "protein_notation": "p.Ala100Val",
                "penetrance_data": {
                    "total_carriers_observed": 99,
                    "affected_count": 99,
                },
                "individual_records": [
                    {
                        "individual_id": "P1",
                        "affected_status": "affected",
                    }
                ],
            },
            {
                "gene_symbol": "RYR2",
                "protein_notation": "p.Gly200Arg",
            },
        ],
    }
    (targeted_dir / "RYR2_PMID_29925740.json").write_text(
        json.dumps(targeted_payload), encoding="utf-8"
    )

    stats = migrate_extraction_directory(
        conn,
        targeted_dir,
        replace_existing_paper=False,
    )

    assert stats["successful"] == 1
    assert conn.execute("SELECT COUNT(*) FROM papers").fetchone()[0] == 1
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 2
    assert conn.execute("SELECT COUNT(*) FROM variant_papers").fetchone()[0] == 2
    assert conn.execute("SELECT COUNT(*) FROM penetrance_data").fetchone()[0] == 1
    assert conn.execute("SELECT COUNT(*) FROM individual_records").fetchone()[0] == 1
    assert (
        conn.execute(
            """
            SELECT pd.total_carriers_observed
            FROM penetrance_data pd
            JOIN variants v ON pd.variant_id = v.variant_id
            WHERE v.protein_notation = 'p.Ala100Val'
            """
        ).fetchone()[0]
        == 3
    )
    conn.close()
