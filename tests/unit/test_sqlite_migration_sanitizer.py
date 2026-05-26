"""Tests for SQLite migration notation guards."""

import json

from harvesting.migrate_to_sqlite import (
    create_database_schema,
    find_extraction_json_files,
    insert_paper_metadata,
    insert_variant_data,
    migrate_extraction_directory,
    normalize_affected_status,
    repair_extraction_data,
)
from utils.pmid_utils import extract_pmid_from_filename


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
