"""Tests for SQLite migration notation guards."""

from harvesting.migrate_to_sqlite import (
    create_database_schema,
    insert_paper_metadata,
    insert_variant_data,
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
