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
