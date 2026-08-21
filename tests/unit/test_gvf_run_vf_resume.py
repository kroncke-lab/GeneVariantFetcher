from __future__ import annotations

import sqlite3

from cli.gvf_run import validate_prior_vf_enrichment


def _seed_vf_db(path, *, uncovered: bool) -> None:
    connection = sqlite3.connect(path)
    try:
        connection.executescript(
            """
            CREATE TABLE variants(variant_id INTEGER PRIMARY KEY);
            CREATE TABLE vf_enrichment(variant_id INTEGER PRIMARY KEY);
            CREATE TABLE quarantined_variants(variant_id INTEGER);
            INSERT INTO variants VALUES (1), (2);
            INSERT INTO vf_enrichment VALUES (1);
            INSERT INTO quarantined_variants VALUES (99);
            """
        )
        if not uncovered:
            connection.execute("INSERT INTO vf_enrichment VALUES (2)")
        connection.commit()
    finally:
        connection.close()


def test_prior_vf_enrichment_is_valid_only_with_full_live_variant_coverage(tmp_path):
    db = tmp_path / "covered.db"
    _seed_vf_db(db, uncovered=False)

    result = validate_prior_vf_enrichment(db)

    assert result == {
        "valid": True,
        "variants": 2,
        "enrichment_rows": 2,
        "uncovered": 0,
        "quarantined": 1,
        "reason": "all live variants retain prior VariantFeatures audit coverage",
    }


def test_prior_vf_enrichment_rejects_uncovered_live_variant(tmp_path):
    db = tmp_path / "uncovered.db"
    _seed_vf_db(db, uncovered=True)

    result = validate_prior_vf_enrichment(db)

    assert result["valid"] is False
    assert result["uncovered"] == 1
