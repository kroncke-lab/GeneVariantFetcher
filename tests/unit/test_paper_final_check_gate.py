"""Actionable paper-final-check trust composition (no live model calls)."""

import json
import sqlite3

from pipeline.paper_final_check import SUMMARY_PROMPT_VERSION
from pipeline.paper_final_check_gate import apply_paper_final_check_gate


def _seed(db):
    conn = sqlite3.connect(db)
    conn.executescript(
        """
        CREATE TABLE papers (
            pmid TEXT PRIMARY KEY, title TEXT, gene_symbol TEXT
        );
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY, gene_symbol TEXT,
            cdna_notation TEXT, protein_notation TEXT, genomic_position TEXT
        );
        CREATE TABLE variant_papers (
            variant_id INTEGER, pmid TEXT, source_location TEXT,
            key_quotes TEXT, count_provenance TEXT
        );
        CREATE TABLE penetrance_data (
            penetrance_id INTEGER PRIMARY KEY,
            variant_id INTEGER,
            pmid TEXT,
            total_carriers_observed INTEGER,
            affected_count INTEGER,
            unaffected_count INTEGER,
            uncertain_count INTEGER,
            trust_tier TEXT,
            trust_reasons TEXT,
            trust_rule_version TEXT
        );
        CREATE TABLE paper_final_check (
            pmid TEXT PRIMARY KEY,
            verdict TEXT,
            flags_json TEXT,
            source_grounded INTEGER,
            n_missing INTEGER,
            check_version TEXT,
            prompt_version TEXT
        );
        """
    )
    conn.executemany(
        "INSERT INTO papers VALUES (?,?,?)",
        [("111", "Paper one", "BRCA2"), ("222", "Paper two", "BRCA2")],
    )
    conn.executemany(
        "INSERT INTO variants VALUES (?,?,?,?,?)",
        [
            (1, "BRCA2", "c.1A>G", "p.Arg1Gly", None),
            (2, "BRCA2", "c.2A>G", "p.Arg2Gly", None),
            (3, "BRCA2", "c.3A>G", "p.Arg3Gly", None),
        ],
    )
    conn.executemany(
        "INSERT INTO variant_papers VALUES (?,?,?,?,?)",
        [
            (1, "111", "Table 2", "[]", "{}"),
            (2, "111", "Table 2", "[]", "{}"),
            (3, "222", "Table 1", "[]", "{}"),
        ],
    )
    conn.executemany(
        "INSERT INTO penetrance_data VALUES (?,?,?,?,?,?,?,?,?,?)",
        [
            (
                10,
                1,
                "111",
                50000,
                0,
                0,
                0,
                "quarantine",
                json.dumps(["population_count"]),
                "tg2-old",
            ),
            (11, 2, "111", 3, 2, 1, 0, "trusted", "[]", "tg2-old"),
            (12, 3, "222", 5, 5, 0, 0, "trusted", "[]", "tg2-old"),
        ],
    )
    flags = [
        {
            "fact_ids": [11],
            "variant_ids": [2],
            "variant": "p.Val30Met",
            "fields": ["affected"],
            "reason_code": "wrong_column",
            "evidence_quote": "Age at diagnosis",
            "evidence_quote_verified": True,
            "issue": "Age at diagnosis was used as affected count",
            "severity": "high",
        },
        {
            "fact_ids": [12],
            "variant_ids": [3],
            "variant": "p.Arg100Ter",
            "fields": ["total_carriers"],
            "reason_code": "unsupported_count",
            "evidence_quote": "Weak provenance",
            "evidence_quote_verified": True,
            "issue": "Weak provenance",
            "severity": "high",
        },
        {
            "fact_ids": [999],
            "variant_ids": [999],
            "variant": "p.Missing",
            "fields": ["total_carriers"],
            "reason_code": "count_is_total",
            "evidence_quote": "Cohort denominator",
            "evidence_quote_verified": True,
            "issue": "Cohort denominator",
            "severity": "high",
        },
    ]
    from pipeline.paper_final_check import gather_paper_payloads, payload_content_hash

    payload_sigs = {
        payload["pmid"]: payload_content_hash(payload)[:12]
        for payload in gather_paper_payloads(conn)
    }
    conn.execute(
        "INSERT INTO paper_final_check VALUES (?,?,?,?,?,?,?)",
        (
            "111",
            "flag",
            json.dumps(flags),
            1,
            1,
            f"pfs9-test-source-{payload_sigs['111']}",
            SUMMARY_PROMPT_VERSION,
        ),
    )
    conn.execute(
        "INSERT INTO paper_final_check VALUES (?,?,?,?,?,?,?)",
        (
            "222",
            "flag",
            json.dumps([flags[1]]),
            1,
            0,
            f"pfs9-test-source-{payload_sigs['222']}",
            SUMMARY_PROMPT_VERSION,
        ),
    )
    conn.commit()
    conn.close()


def _trust_rows(db):
    conn = sqlite3.connect(db)
    rows = conn.execute(
        "SELECT penetrance_id, total_carriers_observed, affected_count, "
        "trust_tier, trust_reasons, field_trust, trust_sources, "
        "trust_rule_version FROM penetrance_data ORDER BY penetrance_id"
    ).fetchall()
    conn.close()
    return rows


def test_gate_quarantines_only_bound_fields_and_preserves_raw_counts(tmp_path):
    db = tmp_path / "gate.db"
    _seed(db)

    stats = apply_paper_final_check_gate(db)

    assert stats["applied_facts"] == 1
    assert stats["applied_bindings"] == 1
    assert stats["unresolved_actionable"] == 1
    assert stats["advisory_flags"] == 2
    assert stats["advisory_reason_flags"] == 2
    assert stats["actionable_flags"] == 2
    assert stats["missing_carriers"] == 1
    assert stats["fields_quarantined"]["affected"] == 1

    rows = _trust_rows(db)
    structural, flagged, clean = rows

    # Structural reasons compose and raw values are untouched.
    assert structural[1:4] == (50000, 0, "quarantine")
    assert json.loads(structural[4]) == ["population_count"]
    assert set(json.loads(structural[5]).values()) == {"quarantine"}

    assert flagged[1:4] == (3, 2, "quarantine")
    assert json.loads(flagged[4]) == ["paper_final_check:wrong_column:high:affected"]
    field_trust = json.loads(flagged[5])
    assert field_trust["affected"] == "quarantine"
    assert field_trust["total_carriers"] == "trusted"
    assert json.loads(flagged[6]) == ["paper_final_check"]

    assert clean[1:4] == (5, 5, "trusted")
    assert json.loads(clean[4]) == []

    conn = sqlite3.connect(db)
    actions = conn.execute(
        "SELECT match_status, action FROM paper_final_check_actions ORDER BY action_id"
    ).fetchall()
    conn.close()
    assert ("exact_fact_id", "quarantine_fields") in actions
    assert ("fact_variant_binding_mismatch", "not_applied") in actions


def test_gate_is_idempotent_and_db_only_flags_do_not_blanket(tmp_path):
    db = tmp_path / "gate.db"
    _seed(db)
    apply_paper_final_check_gate(db)
    first = _trust_rows(db)
    second_stats = apply_paper_final_check_gate(db)
    second = _trust_rows(db)

    assert first == second
    assert second_stats["applied_facts"] == 1

    conn = sqlite3.connect(db)
    conn.execute(
        "UPDATE paper_final_check SET source_grounded=0, n_missing=0 WHERE pmid='111'"
    )
    conn.commit()
    conn.close()

    stats = apply_paper_final_check_gate(db)
    assert stats["ungrounded_actionable"] == 2
    rows = _trust_rows(db)
    assert rows[1][3] == "trusted"
    assert json.loads(rows[1][4]) == []
    # The separate pre-existing structural quarantine remains intact.
    assert rows[0][3] == "quarantine"


def test_gate_refuses_stale_fact_ids_after_payload_change(tmp_path):
    db = tmp_path / "gate.db"
    _seed(db)
    apply_paper_final_check_gate(db)
    conn = sqlite3.connect(db)
    conn.execute(
        "UPDATE penetrance_data SET total_carriers_observed=99 WHERE penetrance_id=11"
    )
    conn.commit()
    conn.close()

    stats = apply_paper_final_check_gate(db)

    assert stats["stale_papers"] == 1
    assert stats["stale_actionable"] == 2
    assert stats["stale_missing_carriers"] == 1
    rows = _trust_rows(db)
    # The stale exact-ID action is not guessed onto the changed payload.
    assert rows[1][3] == "trusted"
    assert json.loads(rows[1][4]) == []


def test_gate_refuses_objective_flag_without_verified_source_quote(tmp_path):
    db = tmp_path / "gate.db"
    _seed(db)
    conn = sqlite3.connect(db)
    raw = conn.execute(
        "SELECT flags_json FROM paper_final_check WHERE pmid='111'"
    ).fetchone()[0]
    flags = json.loads(raw)
    flags[0]["evidence_quote_verified"] = False
    conn.execute(
        "UPDATE paper_final_check SET flags_json=? WHERE pmid='111'",
        (json.dumps(flags),),
    )
    conn.commit()
    conn.close()

    stats = apply_paper_final_check_gate(db)

    assert stats["unverified_actionable"] == 1
    rows = _trust_rows(db)
    assert rows[1][3] == "trusted"
    conn = sqlite3.connect(db)
    statuses = conn.execute(
        "SELECT match_status FROM paper_final_check_actions"
    ).fetchall()
    conn.close()
    assert ("source_quote_unverified",) in statuses


def test_gate_persists_below_threshold_flag_as_advisory(tmp_path):
    db = tmp_path / "gate.db"
    _seed(db)
    conn = sqlite3.connect(db)
    raw = conn.execute(
        "SELECT flags_json FROM paper_final_check WHERE pmid='111'"
    ).fetchone()[0]
    flags = json.loads(raw)
    flags[0]["severity"] = "medium"
    conn.execute(
        "UPDATE paper_final_check SET flags_json=? WHERE pmid='111'",
        (json.dumps(flags),),
    )
    conn.commit()
    conn.close()

    stats = apply_paper_final_check_gate(db, min_severity="high")

    assert stats["advisory_flags"] == 3
    rows = _trust_rows(db)
    assert rows[1][3] == "trusted"
    conn = sqlite3.connect(db)
    action = conn.execute(
        "SELECT severity, match_status, action FROM paper_final_check_actions "
        "WHERE reason_code='wrong_column'"
    ).fetchone()
    conn.close()
    assert action == ("medium", "not_evaluated", "advisory_below_threshold")


def test_gate_preserves_reasonless_legacy_quarantine_in_freshness_hash(tmp_path):
    db = tmp_path / "gate.db"
    _seed(db)
    conn = sqlite3.connect(db)
    conn.execute(
        "UPDATE penetrance_data SET trust_tier='quarantine', trust_reasons=NULL "
        "WHERE penetrance_id=11"
    )
    from pipeline.paper_final_check import gather_paper_payloads, payload_content_hash

    payload = next(p for p in gather_paper_payloads(conn) if p["pmid"] == "111")
    sig = payload_content_hash(payload)[:12]
    conn.execute(
        "UPDATE paper_final_check SET check_version=? WHERE pmid='111'",
        (f"{SUMMARY_PROMPT_VERSION}-test-source-{sig}",),
    )
    conn.commit()
    conn.close()

    stats = apply_paper_final_check_gate(db)

    assert stats["stale_papers"] == 0
    assert stats["applied_facts"] == 1
    rows = _trust_rows(db)
    assert rows[1][3] == "quarantine"


def test_gate_refuses_outdated_reviewer_generation(tmp_path):
    db = tmp_path / "gate.db"
    _seed(db)
    conn = sqlite3.connect(db)
    conn.execute(
        "UPDATE paper_final_check SET prompt_version='pfs-old' WHERE pmid='111'"
    )
    conn.commit()
    conn.close()

    stats = apply_paper_final_check_gate(db)

    assert stats["stale_papers"] == 1
    assert stats["stale_reviewer_version"] == 1
    assert stats["stale_actionable"] == 2
    rows = _trust_rows(db)
    assert rows[1][3] == "trusted"
    conn = sqlite3.connect(db)
    statuses = conn.execute(
        "SELECT match_status FROM paper_final_check_actions WHERE pmid='111'"
    ).fetchall()
    conn.close()
    assert ("stale_reviewer_version",) in statuses


def test_gate_refuses_missing_reviewer_generation(tmp_path):
    db = tmp_path / "gate.db"
    _seed(db)
    conn = sqlite3.connect(db)
    conn.execute("UPDATE paper_final_check SET prompt_version=NULL WHERE pmid='111'")
    conn.commit()
    conn.close()

    stats = apply_paper_final_check_gate(db)

    assert stats["stale_papers"] == 1
    assert stats["stale_reviewer_version"] == 1
    assert stats["stale_actionable"] == 2
    assert stats["applied_facts"] == 0
    rows = _trust_rows(db)
    assert rows[1][3] == "trusted"


def test_gate_enforces_source_verified_phenotype_contradiction(tmp_path):
    db = tmp_path / "gate.db"
    _seed(db)
    conn = sqlite3.connect(db)
    flag = {
        "fact_ids": [12],
        "variant_ids": [3],
        "variant": "p.Arg100Ter",
        "fields": ["affected", "unaffected"],
        "reason_code": "phenotype_misclassified",
        "evidence_quote": "No clinical evidence of the target disease",
        "evidence_quote_verified": True,
        "issue": "Captured affected status contradicts the reported phenotype",
        "severity": "high",
    }
    conn.execute(
        "UPDATE paper_final_check SET flags_json=? WHERE pmid='222'",
        (json.dumps([flag]),),
    )
    conn.commit()
    conn.close()

    stats = apply_paper_final_check_gate(db)

    assert stats["applied_facts"] == 2
    rows = _trust_rows(db)
    phenotype_fields = json.loads(rows[2][5])
    assert phenotype_fields["affected"] == "quarantine"
    assert phenotype_fields["unaffected"] == "quarantine"
    assert phenotype_fields["total_carriers"] == "trusted"
