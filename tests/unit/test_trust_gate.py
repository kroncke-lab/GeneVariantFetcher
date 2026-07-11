"""Gold-free per-fact trust gate: pure structural rules + DB soft-quarantine."""

import json
import sqlite3

from harvesting.migrate_to_sqlite import create_database_schema
from pipeline import trust_gate
from pipeline.trust_gate import apply_trust_gate, decide_tier, evaluate_fact


def test_clean_fact_is_trusted():
    assert evaluate_fact({"carriers": 10, "affected": 6, "unaffected": 3}) == []
    assert decide_tier([]) == "trusted"


def test_arith_inconsistent_parts_exceed_whole():
    # 8 affected + 5 unaffected cannot come from 10 carriers.
    reasons = evaluate_fact({"carriers": 10, "affected": 8, "unaffected": 5})
    assert "arith_inconsistent" in reasons
    assert decide_tier(reasons) == "quarantine"


def test_arith_ok_when_parts_within_total():
    assert evaluate_fact({"carriers": 20, "affected": 8, "unaffected": 5}) == []


def test_count_is_total_role():
    reasons = evaluate_fact(
        {"carriers": 400}, provenance={"carriers_count_type": "cohort_total"}
    )
    assert reasons == ["count_is_total"]


def test_population_count_by_label_and_by_ceiling():
    by_label = evaluate_fact(
        {"carriers": 250}, provenance={"carriers_column_label": "gnomAD allele count"}
    )
    assert "population_count" in by_label
    by_ceiling = evaluate_fact({"carriers": 250_000})
    assert "population_count" in by_ceiling


def test_paper_outlier():
    stats = {"carriers_median": 4, "carriers_n": 5}
    assert "paper_outlier" in evaluate_fact({"carriers": 92}, paper_stats=stats)
    # near the median -> not an outlier
    assert evaluate_fact({"carriers": 6}, paper_stats=stats) == []


def test_paper_outlier_needs_enough_rows():
    # too few carrier rows in the paper -> median is unreliable, no outlier call
    stats = {"carriers_median": 4, "carriers_n": 2}
    assert evaluate_fact({"carriers": 92}, paper_stats=stats) == []


def test_multiple_reasons_accumulate():
    reasons = evaluate_fact(
        {"carriers": 250_000, "affected": 300_000, "unaffected": 0},
        provenance={"carriers_count_type": "screened_n"},
    )
    assert set(reasons) >= {"arith_inconsistent", "count_is_total", "population_count"}


def test_rule_version_is_stable_and_tagged():
    version = trust_gate.rule_version()
    assert version.startswith("tg2-")
    assert version == trust_gate.rule_version()


def test_study_type_mismatch_quarantines_review_functional_gwas():
    reasons = evaluate_fact(
        {"carriers": 12, "affected": 8},
        study_context={"study_design": "review_meta"},
    )
    assert "study_type_mismatch" in reasons
    reasons_fn = evaluate_fact(
        {"carriers": 3},
        study_context={"study_design": "functional_invitro"},
    )
    assert "study_type_mismatch" in reasons_fn


def test_population_study_strengthens_population_count():
    reasons = evaluate_fact(
        {"carriers": 80},
        study_context={"study_design": "cohort_biobank"},
    )
    assert "population_count" in reasons
    # Small clinical counts on case series stay clean on this rule alone.
    assert (
        evaluate_fact(
            {"carriers": 5},
            study_context={"study_design": "case_series"},
        )
        == []
    )


def _seed(conn, pid, pmid, counts, count_provenance=None):
    conn.execute("INSERT OR IGNORE INTO papers (pmid) VALUES (?)", (pmid,))
    conn.execute(
        "INSERT INTO variants (variant_id, gene_symbol, protein_notation) "
        "VALUES (?, 'KCNH2', ?)",
        (pid, f"p.V{pid}M"),
    )
    conn.execute(
        "INSERT INTO variant_papers (variant_id, pmid, count_provenance) "
        "VALUES (?, ?, ?)",
        (pid, pmid, json.dumps(count_provenance) if count_provenance else None),
    )
    conn.execute(
        "INSERT INTO penetrance_data (penetrance_id, variant_id, pmid, "
        "total_carriers_observed, affected_count, unaffected_count) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        (
            pid,
            pid,
            pmid,
            counts.get("carriers"),
            counts.get("affected"),
            counts.get("unaffected"),
        ),
    )


def test_apply_trust_gate_soft_quarantines_and_preserves_counts(tmp_path):
    db = str(tmp_path / "t.db")
    conn = create_database_schema(db)
    try:
        _seed(conn, 1, "111", {"carriers": 10, "affected": 6, "unaffected": 3})
        _seed(conn, 2, "111", {"carriers": 200_000})  # population ceiling
        conn.commit()
    finally:
        conn.close()

    stats = apply_trust_gate(db)
    assert stats["trusted"] == 1
    assert stats["quarantine"] == 1
    assert stats["by_reason"].get("population_count") == 1
    assert stats["rule_version"].startswith("tg2-")

    conn = sqlite3.connect(db)
    try:
        rows = {
            r[0]: r
            for r in conn.execute(
                "SELECT penetrance_id, trust_tier, trust_reasons, "
                "total_carriers_observed FROM penetrance_data"
            )
        }
    finally:
        conn.close()

    assert rows[1][1] == "trusted"
    assert rows[2][1] == "quarantine"
    assert "population_count" in rows[2][2]
    # soft-quarantine: the raw count is preserved, never NULLed.
    assert rows[2][3] == 200_000


def test_apply_trust_gate_is_idempotent(tmp_path):
    db = str(tmp_path / "t.db")
    conn = create_database_schema(db)
    try:
        _seed(conn, 1, "111", {"carriers": 10, "affected": 6, "unaffected": 3})
        conn.commit()
    finally:
        conn.close()
    first = apply_trust_gate(db)
    second = apply_trust_gate(db)
    assert first["trusted"] == second["trusted"] == 1
    assert first["quarantine"] == second["quarantine"] == 0


def test_apply_trust_gate_dedupes_multiple_extraction_metadata_rows(tmp_path):
    """A pmid with >1 extraction_metadata row must not fan a fact out once per
    row (which double-counts stats and lets an arbitrary row win the tier). The
    latest metadata row (max extraction_id) decides the study context."""
    db = str(tmp_path / "t.db")
    conn = create_database_schema(db)
    try:
        _seed(conn, 1, "111", {"carriers": 10, "affected": 6, "unaffected": 3})
        # Two metadata rows for the same pmid. The OLDER (extraction_id=1) is a
        # review (study_type_mismatch -> quarantine); the NEWER (extraction_id=2)
        # is a case series (clean). The newer must win.
        conn.execute(
            "INSERT INTO extraction_metadata (extraction_id, pmid, study_design) "
            "VALUES (1, '111', 'review_meta')"
        )
        conn.execute(
            "INSERT INTO extraction_metadata (extraction_id, pmid, study_design) "
            "VALUES (2, '111', 'case_series')"
        )
        conn.commit()
    finally:
        conn.close()

    stats = apply_trust_gate(db)
    # Exactly one fact tiered — no fan-out (would be 2 before the fix).
    assert stats["trusted"] + stats["quarantine"] == 1
    # Latest metadata (case_series) wins -> trusted, not the older review_meta.
    assert stats["trusted"] == 1
    assert stats["quarantine"] == 0

    conn = sqlite3.connect(db)
    try:
        tier = conn.execute(
            "SELECT trust_tier FROM penetrance_data WHERE penetrance_id = 1"
        ).fetchone()[0]
    finally:
        conn.close()
    assert tier == "trusted"
