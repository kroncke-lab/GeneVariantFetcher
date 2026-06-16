"""Migration idempotency: exact-duplicate child rows must never accumulate.

Regression guard for the carriers-MAE blow-up where a variant represented across
several supplement-table cells (KCNQ1 PMID 32893267) produced N identical
penetrance rows that the scorer summed into an N-fold carrier count.
"""

import copy
import json
import sqlite3

from harvesting.migrate_to_sqlite import (
    OBSERVATION_PROVENANCE_KEYS,
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


def _fact_count(conn):
    return conn.execute("SELECT COUNT(*) FROM fact_provenance").fetchone()[0]


def _table_columns(conn, table):
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}


def test_observation_provenance_columns_exist_on_child_tables(tmp_path):
    conn = create_database_schema(str(tmp_path / "t.db"))

    for table in ("individual_records", "phenotypes"):
        columns = _table_columns(conn, table)
        for column in OBSERVATION_PROVENANCE_KEYS:
            assert column in columns
    conn.close()


def test_insert_variant_data_round_trips_observation_provenance(tmp_path):
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    insert_variant_data(
        cur,
        "32893267",
        {
            "gene_symbol": "KCNQ1",
            "protein_notation": "p.Val254Met",
            "source_location": "Supplementary Table 2, row 4",
            "patients": {
                "count": 1,
                "phenotype": "LQTS",
                "source_container": "supplement",
                "source_kind": "table",
                "source_ref": "Supplementary Table 2",
                "page_label": "e12",
                "pdf_page": 7,
                "row_label": "Patient 4",
                "row_ordinal": 4,
                "column_ref": "Phenotype",
                "locator_extra": {"cell_coords": {"row": 4, "col": 6}},
            },
            "individual_records": [
                {
                    "individual_id": "Patient 4",
                    "affected_status": "affected",
                    "phenotype_details": "QT prolongation",
                    "evidence_sentence": "Patient 4 had QT prolongation.",
                    "source_container": "supplement",
                    "source_kind": "table",
                    "source_ref": "Supplementary Table 2",
                    "page_label": "e12",
                    "pdf_page": 7,
                    "row_label": "Patient 4",
                    "row_ordinal": 4,
                    "column_ref": "Phenotype",
                    "locator_extra": {"bbox": [1, 2, 3, 4]},
                }
            ],
        },
    )
    conn.commit()

    individual = conn.execute(
        """
        SELECT source_container, source_kind, source_ref, page_label, pdf_page,
               row_label, row_ordinal, column_ref, figure_panel, source_record_id,
               locator_extra
        FROM individual_records
        """
    ).fetchone()
    assert individual[:9] == (
        "supplement",
        "table",
        "Supplementary Table 2",
        "e12",
        7,
        "Patient 4",
        4,
        "Phenotype",
        None,
    )
    assert len(individual[9]) == 64
    assert json.loads(individual[10]) == {"bbox": [1, 2, 3, 4]}

    phenotype = conn.execute(
        """
        SELECT source_container, source_kind, source_ref, page_label, pdf_page,
               row_label, row_ordinal, column_ref, source_record_id, locator_extra
        FROM phenotypes
        """
    ).fetchone()
    assert phenotype[:8] == (
        "supplement",
        "table",
        "Supplementary Table 2",
        "e12",
        7,
        "Patient 4",
        4,
        "Phenotype",
    )
    assert len(phenotype[8]) == 64
    assert json.loads(phenotype[9]) == {"cell_coords": {"col": 6, "row": 4}}
    conn.close()


def test_source_record_id_is_stable_across_reruns(tmp_path):
    variant = {
        "gene_symbol": "KCNQ1",
        "protein_notation": "p.Val254Met",
        "individual_records": [
            {
                "individual_id": "II-1",
                "affected_status": "affected",
                "source_ref": "Table 2",
                "row_label": "II-1",
            }
        ],
    }
    ids = []
    for db_name in ("first.db", "second.db"):
        conn = create_database_schema(str(tmp_path / db_name))
        cur = conn.cursor()
        cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
        insert_variant_data(cur, "32893267", copy.deepcopy(variant))
        conn.commit()
        ids.append(
            conn.execute("SELECT source_record_id FROM individual_records").fetchone()[
                0
            ]
        )
        conn.close()

    assert ids[0] == ids[1]


def test_insert_variant_data_allows_rows_without_provenance(tmp_path):
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    insert_variant_data(
        cur,
        "32893267",
        {
            "gene_symbol": "KCNQ1",
            "protein_notation": "p.Val254Met",
            "patients": {"count": 1, "phenotype": "LQTS"},
            "individual_records": [
                {"individual_id": "P1", "affected_status": "affected"}
            ],
        },
    )
    conn.commit()

    individual = conn.execute(
        """
        SELECT source_container, source_kind, source_ref, source_record_id
        FROM individual_records
        """
    ).fetchone()
    phenotype = conn.execute(
        "SELECT source_container, source_kind, source_ref, source_record_id FROM phenotypes"
    ).fetchone()
    assert individual[:3] == (None, None, None)
    assert phenotype[:3] == (None, None, None)
    assert len(individual[3]) == 64
    assert len(phenotype[3]) == 64
    conn.close()


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


def test_insert_variant_data_writes_source_layer(tmp_path):
    """variant_papers carries the stable provenance layer used by recall scoring."""
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    insert_variant_data(
        cur,
        "32893267",
        {
            "gene_symbol": "KCNQ1",
            "protein_notation": "p.Val254Met",
            "source_location": "Supplement Table S4",
        },
    )
    insert_variant_data(
        cur,
        "32893267",
        {
            "gene_symbol": "KCNQ1",
            "protein_notation": "p.Ala341Val",
            "source_layer": "pubtator",
        },
    )
    conn.commit()

    rows = conn.execute(
        """
        SELECT v.protein_notation, vp.source_layer
        FROM variant_papers vp
        JOIN variants v ON v.variant_id = vp.variant_id
        ORDER BY v.protein_notation
        """
    ).fetchall()
    assert rows == [("p.Ala341Val", "pubtator"), ("p.Val254Met", "llm_table")]
    conn.close()


def test_insert_variant_data_writes_standard_fact_provenance_idempotently(tmp_path):
    """Every standard fact gets source pointers, and replay does not duplicate them."""
    conn = create_database_schema(str(tmp_path / "t.db"))
    cur = conn.cursor()
    cur.execute("INSERT OR IGNORE INTO papers (pmid) VALUES ('32893267')")
    variant = {
        "gene_symbol": "KCNQ1",
        "protein_notation": "p.Val254Met",
        "source_location": "Table S1, row 20",
        "key_quotes": ["p.Val254Met had 3 carriers: 2 affected and 1 unaffected."],
        "patients": {"count": 3},
        "penetrance_data": {
            "total_carriers_observed": 3,
            "affected_count": 2,
            "unaffected_count": 1,
        },
        "count_provenance": {
            "carriers_column_label": "N carriers",
            "carriers_count_type": "per_variant_carrier",
            "affected_column_label": "Affected",
            "affected_count_type": "per_variant_carrier",
            "unaffected_column_label": "Unaffected",
            "unaffected_count_type": "per_variant_carrier",
        },
        "individual_records": [
            {
                "individual_id": "P1",
                "affected_status": "affected",
                "evidence_sentence": "P1 carried p.Val254Met and was affected.",
            }
        ],
        "fact_provenance": [
            {
                "fact_type": "affected_count",
                "fact_value": 2,
                "source_table": "Table S1",
                "source_row": "row 20",
                "source_column": "Affected",
                "evidence_quote": "2 affected",
            }
        ],
    }

    insert_variant_data(cur, "32893267", copy.deepcopy(variant))
    conn.commit()
    first_count = _fact_count(conn)
    insert_variant_data(cur, "32893267", copy.deepcopy(variant))
    conn.commit()

    assert _fact_count(conn) == first_count
    assert first_count >= 8
    assert (
        conn.execute(
            "SELECT COUNT(DISTINCT fact_hash) FROM fact_provenance"
        ).fetchone()[0]
        == first_count
    )

    total = conn.execute(
        """
        SELECT fact_value, source_table, source_row, source_column, count_type
        FROM fact_provenance
        WHERE fact_type = 'total_carriers_observed'
        """
    ).fetchone()
    assert total == ("3", "Table S1", "row 20", "N carriers", "per_variant_carrier")

    individual = conn.execute(
        """
        SELECT fact_value, individual_id, evidence_quote
        FROM fact_provenance
        WHERE fact_type = 'individual_affected_status'
        """
    ).fetchone()
    assert individual == ("affected", "P1", "P1 carried p.Val254Met and was affected.")

    explicit = conn.execute(
        """
        SELECT source_table, source_row, source_column, provenance_kind
        FROM fact_provenance
        WHERE fact_type = 'affected_count' AND evidence_quote = '2 affected'
        """
    ).fetchone()
    assert explicit == ("Table S1", "row 20", "Affected", "explicit")
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
    assert _fact_count(conn) == 3
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
