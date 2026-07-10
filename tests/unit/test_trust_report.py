"""Trust-tier report over a two-tier DB."""

from harvesting.migrate_to_sqlite import create_database_schema
from pipeline.trust_gate import apply_trust_gate
from scripts.trust_report import list_quarantined, summarize_trust


def _seed(conn, pid, pmid, carriers, affected=None, unaffected=None):
    conn.execute("INSERT OR IGNORE INTO papers (pmid) VALUES (?)", (pmid,))
    conn.execute(
        "INSERT INTO variants (variant_id, gene_symbol, protein_notation) "
        "VALUES (?, 'KCNH2', ?)",
        (pid, f"p.V{pid}M"),
    )
    conn.execute(
        "INSERT INTO variant_papers (variant_id, pmid) VALUES (?, ?)", (pid, pmid)
    )
    conn.execute(
        "INSERT INTO penetrance_data (penetrance_id, variant_id, pmid, "
        "total_carriers_observed, affected_count, unaffected_count) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        (pid, pid, pmid, carriers, affected, unaffected),
    )


def test_summarize_and_list(tmp_path):
    db = str(tmp_path / "t.db")
    conn = create_database_schema(db)
    try:
        _seed(conn, 1, "111", 10, affected=6, unaffected=3)  # clean → trusted
        _seed(conn, 2, "111", 250_000)  # population ceiling → quarantine
        _seed(conn, 3, "111", 12, affected=20, unaffected=5)  # arith → quarantine
        conn.commit()
    finally:
        conn.close()

    apply_trust_gate(db)

    summary = summarize_trust(db)
    assert summary["tiered"] is True
    assert summary["total"] == 3
    assert summary["trusted"] == 1
    assert summary["quarantine"] == 2
    assert abs(summary["quarantine_rate"] - 2 / 3) < 1e-6
    assert summary["by_reason"].get("population_count") == 1
    assert summary["by_reason"].get("arith_inconsistent") == 1

    listed = list_quarantined(db, limit=10)
    assert len(listed) == 2
    # highest carrier count (the population one) sorts first
    assert listed[0]["carriers"] == 250_000
    assert "population_count" in listed[0]["reasons"]


def test_untiered_db_reports_gracefully(tmp_path):
    """A DB without the trust_tier column reports tiered=False, not a crash."""
    import sqlite3

    db = str(tmp_path / "old.db")
    conn = sqlite3.connect(db)
    try:
        conn.execute(
            "CREATE TABLE penetrance_data (penetrance_id INTEGER PRIMARY KEY, "
            "total_carriers_observed INTEGER)"
        )
        conn.execute(
            "INSERT INTO penetrance_data (penetrance_id, total_carriers_observed) "
            "VALUES (1, 5)"
        )
        conn.commit()
    finally:
        conn.close()

    summary = summarize_trust(db)
    assert summary["tiered"] is False
    assert summary["total"] == 1
