"""Primary and all-observed source-layer contracts."""

import sqlite3

from harvesting.migrate_to_sqlite import create_database_schema
from scripts.backfill_source_layers import backfill
from utils.source_layers import (
    add_source_layer_witness,
    combine_source_layers,
    normalize_source_layer,
    source_layer_sql_case,
    source_layer_tokens,
)


def test_legacy_composite_keeps_origin_as_primary():
    assert normalize_source_layer("llm_table,figure") == "llm_table"
    assert normalize_source_layer("regex_text|figure") == "regex_text"
    assert normalize_source_layer("not-a-layer,figure") is None


def test_observed_layers_are_primary_first_and_deduplicated():
    assert (
        combine_source_layers("llm_table", "figure,llm_table", "clinvar", "figure")
        == "llm_table,figure,clinvar"
    )
    assert source_layer_tokens("llm_table", "figure,clinvar") == {
        "llm_table",
        "figure",
        "clinvar",
    }


def test_recovery_witness_does_not_claim_missing_legacy_primary():
    primary, observed = add_source_layer_witness(
        source_layer=None,
        observed_source_layers=None,
        witness_layer="clinvar",
        source_location="Table 2",
    )

    assert primary == "llm_table"
    assert observed == "llm_table,clinvar"


def test_sql_classifier_reads_legacy_composite_primary():
    con = sqlite3.connect(":memory:")
    expr = source_layer_sql_case("source_location", "source_layer")
    row = con.execute(
        f"""SELECT {expr}
            FROM (SELECT 'figure-reader' AS source_location,
                         'llm_table,figure' AS source_layer)"""
    ).fetchone()
    con.close()
    assert row == ("llm_table",)


def test_backfill_splits_legacy_primary_without_losing_witnesses(tmp_path):
    db = tmp_path / "legacy.db"
    con = create_database_schema(str(db))
    con.execute("INSERT INTO papers(pmid) VALUES ('111')")
    con.execute(
        """INSERT INTO variants(variant_id, gene_symbol, protein_notation)
           VALUES (1, 'KCNH2', 'p.Arg123His')"""
    )
    con.execute(
        """INSERT INTO variant_papers(
               variant_id, pmid, source_location, source_layer,
               observed_source_layers
           ) VALUES (1, '111', 'Table 2', 'llm_table,figure', NULL)"""
    )
    con.commit()
    con.close()

    backup, updated = backfill(db, backup_suffix="test")

    con = sqlite3.connect(db)
    row = con.execute(
        "SELECT source_layer, observed_source_layers FROM variant_papers"
    ).fetchone()
    con.close()
    assert updated == 1
    assert row == ("llm_table", "llm_table,figure")
    assert backup.is_file()
