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
    assert version.startswith("tg3-")
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
    assert stats["rule_version"].startswith("tg3-")

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


def test_apply_trust_gate_preserves_composed_paper_field_mask(tmp_path):
    db = str(tmp_path / "t.db")
    conn = create_database_schema(db)
    try:
        _seed(conn, 1, "111", {"carriers": 10, "affected": 6, "unaffected": 3})
        conn.execute(
            "UPDATE penetrance_data SET trust_tier='quarantine', "
            "trust_reasons=?, trust_rule_version='tg2-old|pfcg1-seal', "
            "field_trust=?, trust_sources='[\"paper_final_check\"]' "
            "WHERE penetrance_id=1",
            (
                '["paper_final_check:column_role_mismatch:high:affected"]',
                '{"total_carriers":"trusted","affected":"quarantine",'
                '"unaffected":"trusted","uncertain":"trusted"}',
            ),
        )
        conn.commit()
    finally:
        conn.close()

    stats = apply_trust_gate(db)
    assert stats["quarantine"] == 1

    conn = sqlite3.connect(db)
    try:
        row = conn.execute(
            "SELECT affected_count, trust_reasons, trust_rule_version, "
            "field_trust, trust_sources FROM penetrance_data WHERE penetrance_id=1"
        ).fetchone()
    finally:
        conn.close()
    assert row[0] == 6
    assert "paper_final_check:column_role_mismatch:high:affected" in row[1]
    assert "pfcg1-seal" in row[2]
    assert json.loads(row[3])["affected"] == "quarantine"
    assert json.loads(row[3])["total_carriers"] == "trusted"
    assert json.loads(row[4]) == ["paper_final_check"]


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


def test_negative_count_is_quarantined():
    reasons = evaluate_fact({"carriers": 10, "affected": -1})
    assert "negative_count" in reasons
    assert "arith_inconsistent" in reasons
    assert decide_tier(reasons) == "quarantine"


def test_full_partition_must_equal_total():
    # affected + unaffected + uncertain all reported but sum != total.
    reasons = evaluate_fact(
        {"carriers": 10, "affected": 4, "unaffected": 3, "uncertain": 1}
    )
    assert "arith_inconsistent" in reasons


def test_full_partition_that_balances_is_clean():
    assert (
        evaluate_fact({"carriers": 10, "affected": 6, "unaffected": 3, "uncertain": 1})
        == []
    )


def test_partial_partition_under_total_is_left_alone():
    # uncertain is None -> not a full partition -> a residual of unphenotyped
    # carriers is fine, not an inconsistency.
    assert evaluate_fact({"carriers": 20, "affected": 8, "unaffected": 5}) == []


def _implied_zero_counts():
    return {"carriers": 12, "affected": 12, "unaffected": 0}


def test_implied_unaffected_zero_fires_for_cohort_design():
    reasons = evaluate_fact(
        _implied_zero_counts(),
        provenance={"unaffected_column_label": None, "unaffected_count_type": None},
        study_context={"study_design": "cohort_population"},
    )
    assert "implied_unaffected_zero" in reasons


def test_implied_unaffected_zero_fires_for_family_cascade_ascertainment():
    reasons = evaluate_fact(
        _implied_zero_counts(),
        study_context={"ascertainment": "family_cascade"},
    )
    assert "implied_unaffected_zero" in reasons


def test_implied_unaffected_zero_dormant_on_unknown_design():
    # Cardiac legacy default: no study_design -> the derived zero is trusted.
    assert "implied_unaffected_zero" not in evaluate_fact(
        _implied_zero_counts(), study_context={}
    )


def test_implied_unaffected_zero_dormant_for_case_report():
    assert "implied_unaffected_zero" not in evaluate_fact(
        _implied_zero_counts(), study_context={"study_design": "case_report"}
    )


def test_implied_unaffected_zero_dormant_when_unaffected_is_sourced():
    # A real unaffected column reporting 0 is trusted even in a cohort study.
    reasons = evaluate_fact(
        _implied_zero_counts(),
        provenance={
            "unaffected_column_label": "Unaffected carriers",
            "unaffected_count_type": "per_variant_carrier",
        },
        study_context={"study_design": "cohort_population"},
    )
    assert "implied_unaffected_zero" not in reasons


def test_implied_unaffected_zero_masks_only_unaffected_field(tmp_path):
    db = str(tmp_path / "t.db")
    create_database_schema(db)
    conn = sqlite3.connect(db)
    try:
        conn.execute(
            "INSERT INTO variants (variant_id, gene_symbol, protein_notation) "
            "VALUES (1, 'BRCA1', 'p.C61G')"
        )
        conn.execute("INSERT INTO papers (pmid, title) VALUES ('1', 't')")
        conn.execute(
            "INSERT INTO variant_papers (variant_id, pmid, count_provenance) "
            "VALUES (1, '1', ?)",
            (
                json.dumps(
                    {"unaffected_column_label": None, "unaffected_count_type": None}
                ),
            ),
        )
        conn.execute(
            "INSERT INTO penetrance_data (penetrance_id, variant_id, pmid, "
            "total_carriers_observed, affected_count, unaffected_count) "
            "VALUES (1, 1, '1', 12, 12, 0)"
        )
        conn.execute(
            "INSERT INTO extraction_metadata (pmid, study_design) VALUES ('1', 'cohort_population')"
        )
        conn.commit()
    finally:
        conn.close()

    apply_trust_gate(db)

    conn = sqlite3.connect(db)
    try:
        row = conn.execute(
            "SELECT trust_tier, field_trust, affected_count, unaffected_count "
            "FROM penetrance_data WHERE penetrance_id = 1"
        ).fetchone()
    finally:
        conn.close()
    tier, field_trust_raw, affected, unaffected = row
    field_trust = json.loads(field_trust_raw)
    assert tier == "quarantine"
    # Only the unaffected field is masked; affected/total stay trusted.
    assert field_trust["unaffected"] == "quarantine"
    assert field_trust["affected"] == "trusted"
    assert field_trust["total_carriers"] == "trusted"
    # Raw counts are never mutated.
    assert affected == 12 and unaffected == 0
