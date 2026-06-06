"""Migration idempotency: exact-duplicate child rows must never accumulate.

Regression guard for the carriers-MAE blow-up where a variant represented across
several supplement-table cells (KCNQ1 PMID 32893267) produced N identical
penetrance rows that the scorer summed into an N-fold carrier count.
"""

import copy
import json
import sqlite3

from harvesting.migrate_to_sqlite import (
    create_database_schema,
    dedup_existing_rows,
    insert_variant_data,
    migrate_extraction_file,
)


def _carrier_rows(conn, pmid="32893267"):
    return conn.execute(
        "SELECT total_carriers_observed FROM penetrance_data WHERE pmid = ?",
        (pmid,),
    ).fetchall()


def test_insert_variant_data_collapses_exact_duplicate_penetrance(tmp_path):
    """One variant extracted from 4 identical table cells -> a single 25, not 100."""
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    variant = {
        "gene_symbol": "KCNQ1",
        "protein_notation": "p.Val254Met",
        "penetrance_data": {
            "total_carriers_observed": 25,
            "affected_count": 25,
            "unaffected_count": 0,
        },
    }
    for _ in range(4):
        insert_variant_data(cur, "32893267", copy.deepcopy(variant))
    conn.commit()

    rows = _carrier_rows(conn)
    assert len(rows) == 1
    assert rows[0][0] == 25
    conn.close()


def test_insert_variant_data_keeps_genuinely_different_penetrance(tmp_path):
    """Real sub-cohort splits (different counts) are preserved, not collapsed."""
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    base = {"gene_symbol": "KCNQ1", "protein_notation": "p.Val254Met"}
    for carriers in (18, 7):
        v = copy.deepcopy(base)
        v["penetrance_data"] = {"total_carriers_observed": carriers}
        insert_variant_data(cur, "32893267", v)
    conn.commit()

    rows = sorted(r[0] for r in _carrier_rows(conn))
    assert rows == [7, 18]
    conn.close()


def test_migrate_extraction_file_twice_is_idempotent(tmp_path):
    """Re-migrating the same file adds no duplicate child rows."""
    extraction = {
        "paper_metadata": {"pmid": "32893267", "title": "t", "journal": "j"},
        "variants": [
            {
                "gene_symbol": "KCNQ1",
                "protein_notation": "p.Val254Met",
                "penetrance_data": {"total_carriers_observed": 25},
                "functional_data": {"summary": "loss of function"},
                "patients": {"count": 25, "phenotype": "LQTS"},
            }
        ],
        "extraction_metadata": {"total_variants_found": 1},
    }
    jf = tmp_path / "KCNQ1_PMID_32893267.json"
    jf.write_text(json.dumps(extraction))

    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    for _ in range(2):
        ok, _msg = migrate_extraction_file(cur, jf, replace_existing_paper=True)
        assert ok
    conn.commit()

    assert len(_carrier_rows(conn)) == 1
    assert _carrier_rows(conn)[0][0] == 25
    for table in ("functional_data", "phenotypes"):
        n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        assert n == 1, f"{table} duplicated on re-migration"
    conn.close()


def test_dedup_existing_rows_collapses_then_is_idempotent(tmp_path):
    """Back-fill: exact dups already in a DB are collapsed; re-running removes 0."""
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute(
        "INSERT INTO variants (gene_symbol, protein_notation) VALUES ('KCNQ1','p.Val254Met')"
    )
    vid = cur.lastrowid
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    # 4 byte-for-byte identical penetrance rows (the pre-guard corruption)
    for _ in range(4):
        cur.execute(
            "INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed) "
            "VALUES (?, '32893267', 25)",
            (vid,),
        )
    # plus one genuinely different row that must survive
    cur.execute(
        "INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed) "
        "VALUES (?, '32893267', 7)",
        (vid,),
    )
    conn.commit()

    removed = dedup_existing_rows(conn)
    assert removed["penetrance_data"] == 3
    survivors = sorted(r[0] for r in _carrier_rows(conn))
    assert survivors == [7, 25]

    removed_again = dedup_existing_rows(conn)
    assert removed_again["penetrance_data"] == 0
    conn.close()
