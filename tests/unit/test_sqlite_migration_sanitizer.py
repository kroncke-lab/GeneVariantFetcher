"""Tests for SQLite migration notation guards."""

import json

from harvesting.migrate_to_sqlite import (
    create_database_schema,
    find_extraction_json_files,
    insert_paper_metadata,
    insert_variant_data,
    migrate_extraction_directory,
)


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
