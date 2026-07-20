"""Verbatim source notation is persisted alongside the normalized variant.

Criticism 7: the original paper notation must travel with the normalized form so
a curator can trace a variant back to the source and catch normalization errors.
"""

import sqlite3

from harvesting.migrate_to_sqlite import create_database_schema, insert_variant_data


def _conn(tmp_path):
    db = str(tmp_path / "t.db")
    create_database_schema(db)
    conn = sqlite3.connect(db)
    conn.execute("INSERT INTO papers (pmid, title) VALUES ('1', 't')")
    return conn


def _source_notation(conn, variant_id):
    return conn.execute(
        "SELECT source_notation FROM variant_papers WHERE variant_id = ? AND pmid = '1'",
        (variant_id,),
    ).fetchone()[0]


def test_explicit_source_notation_is_persisted(tmp_path):
    conn = _conn(tmp_path)
    try:
        vid = insert_variant_data(
            conn.cursor(),
            "1",
            {
                "gene_symbol": "BRCA1",
                "protein_notation": "p.Arg1443Ter",
                "source_notation": "R1443X",
            },
        )
        conn.commit()
        assert _source_notation(conn, vid) == "R1443X"
    finally:
        conn.close()


def test_source_notation_falls_back_to_presanitize_notation(tmp_path):
    conn = _conn(tmp_path)
    try:
        vid = insert_variant_data(
            conn.cursor(),
            "1",
            {"gene_symbol": "BRCA1", "protein_notation": "p.Cys61Gly"},
        )
        conn.commit()
        # No verbatim provided -> the pre-normalization notation is retained.
        assert _source_notation(conn, vid) == "p.Cys61Gly"
    finally:
        conn.close()
