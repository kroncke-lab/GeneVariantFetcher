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


def _age_rows(conn):
    return conn.execute(
        """
        SELECT pd.total_carriers_observed, adp.age_range, adp.carriers_in_range
        FROM age_dependent_penetrance adp
        JOIN penetrance_data pd ON pd.penetrance_id = adp.penetrance_id
        ORDER BY adp.age_range
        """
    ).fetchall()


def _age_row_count(conn):
    return conn.execute("SELECT COUNT(*) FROM age_dependent_penetrance").fetchone()[0]


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


def test_insert_variant_data_merges_age_bins_from_duplicate_penetrance(tmp_path):
    """Duplicate parent rows can still carry unique age-stratified child facts."""
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    base = {
        "gene_symbol": "KCNQ1",
        "protein_notation": "p.Val254Met",
        "penetrance_data": {"total_carriers_observed": 25},
    }
    first = copy.deepcopy(base)
    first["penetrance_data"]["age_dependent_penetrance"] = [
        {"age_range": "0-20", "carriers_in_range": 10}
    ]
    second = copy.deepcopy(base)
    second["penetrance_data"]["age_dependent_penetrance"] = [
        {"age_range": "20-40", "carriers_in_range": 15}
    ]

    insert_variant_data(cur, "32893267", first)
    insert_variant_data(cur, "32893267", second)
    insert_variant_data(cur, "32893267", second)
    conn.commit()

    assert len(_carrier_rows(conn)) == 1
    assert _age_rows(conn) == [(25, "0-20", 10), (25, "20-40", 15)]
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


def test_dedup_existing_rows_reparents_age_bins_before_parent_cascade(tmp_path):
    """Back-fill preserves unique age-bin rows before duplicate parents cascade."""
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute(
        "INSERT INTO variants (gene_symbol, protein_notation) VALUES ('KCNQ1','p.Val254Met')"
    )
    vid = cur.lastrowid
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    parent_ids = []
    for age_range, carriers in (("0-20", 10), ("20-40", 15)):
        cur.execute(
            "INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed) "
            "VALUES (?, '32893267', 25)",
            (vid,),
        )
        parent_ids.append(cur.lastrowid)
        cur.execute(
            """
            INSERT INTO age_dependent_penetrance (
                penetrance_id, age_range, carriers_in_range
            ) VALUES (?, ?, ?)
            """,
            (parent_ids[-1], age_range, carriers),
        )
    conn.commit()

    removed = dedup_existing_rows(conn)
    assert removed["penetrance_data"] == 1
    assert len(_carrier_rows(conn)) == 1
    assert _age_rows(conn) == [(25, "0-20", 10), (25, "20-40", 15)]

    removed_again = dedup_existing_rows(conn)
    assert removed_again["penetrance_data"] == 0
    assert _age_rows(conn) == [(25, "0-20", 10), (25, "20-40", 15)]
    conn.close()


def test_dedup_existing_rows_cleans_age_bin_orphans_when_fk_cascade_is_off(tmp_path):
    """Back-fill is robust even when PRAGMA foreign_keys is ignored in a txn."""
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute(
        "INSERT INTO variants (gene_symbol, protein_notation) VALUES ('KCNQ1','p.Val254Met')"
    )
    vid = cur.lastrowid
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    parent_ids = []
    for age_range, carriers in (("0-20", 10), ("20-40", 15)):
        cur.execute(
            "INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed) "
            "VALUES (?, '32893267', 25)",
            (vid,),
        )
        parent_ids.append(cur.lastrowid)
        cur.execute(
            """
            INSERT INTO age_dependent_penetrance (
                penetrance_id, age_range, carriers_in_range
            ) VALUES (?, ?, ?)
            """,
            (parent_ids[-1], age_range, carriers),
        )
    conn.commit()

    conn.execute("PRAGMA foreign_keys = OFF")
    conn.execute("BEGIN")
    removed = dedup_existing_rows(conn)

    assert removed["penetrance_data"] == 1
    assert removed["age_dependent_penetrance_orphans"] == 1
    assert len(_carrier_rows(conn)) == 1
    assert _age_rows(conn) == [(25, "0-20", 10), (25, "20-40", 15)]
    assert _age_row_count(conn) == 2
    conn.close()
