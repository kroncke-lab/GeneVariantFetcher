"""Unit tests for the per-paper final check (sniff test).

Covers the pure helpers (version, normalize, gather) and the SQLite adapter with
an injected stub checker — no network, no live model.
"""

import json
import sqlite3

from pipeline.paper_final_check import (
    apply_paper_final_check,
    build_paper_check_prompt,
    check_version,
    gather_paper_payloads,
    normalize_result,
)


def _seed(path, *, with_trust=False, extra_facts_pmid111=0):
    trust_cols = ", trust_tier TEXT, trust_reasons TEXT" if with_trust else ""
    conn = sqlite3.connect(path)
    conn.executescript(
        f"""
        CREATE TABLE papers (pmid TEXT PRIMARY KEY, title TEXT, gene_symbol TEXT);
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY, gene_symbol TEXT,
            cdna_notation TEXT, protein_notation TEXT, genomic_position TEXT
        );
        CREATE TABLE variant_papers (
            variant_id INTEGER, pmid TEXT, source_location TEXT,
            additional_notes TEXT, key_quotes TEXT, count_provenance TEXT
        );
        CREATE TABLE penetrance_data (
            penetrance_id INTEGER PRIMARY KEY, variant_id INTEGER, pmid TEXT,
            total_carriers_observed INTEGER, affected_count INTEGER,
            unaffected_count INTEGER, uncertain_count INTEGER,
            penetrance_percentage REAL{trust_cols}
        );
        """
    )
    conn.executemany(
        "INSERT INTO papers VALUES (?,?,?)",
        [
            ("111", "A BRCA2 case series", "BRCA2"),
            ("222", "A cohort study", "BRCA2"),
            ("333", "No counts here", "BRCA2"),
        ],
    )
    conn.executemany(
        "INSERT INTO variants VALUES (?,?,?,?,?)",
        [
            (1, "BRCA2", None, "p.Val30Met", None),
            (2, "BRCA2", "c.68A>G", None, None),
            (3, "BRCA2", None, "p.Arg100Ter", None),
        ],
    )
    conn.executemany(
        "INSERT INTO variant_papers VALUES (?,?,?,?,?,?)",
        [
            (
                1,
                "111",
                "Table 2",
                "",
                "3 carriers, 2 affected",
                json.dumps(
                    {
                        "carriers_count_type": "per_variant_carrier",
                        "carriers_column_label": "Carriers",
                    }
                ),
            ),
            (
                2,
                "111",
                "Table 2",
                "",
                "",
                json.dumps(
                    {
                        "carriers_count_type": "cohort_total",
                        "carriers_column_label": "gnomAD AF",
                    }
                ),
            ),
            (3, "222", "Abstract", "", "", None),
        ],
    )

    def pen(pid, vid, pmid, tot, aff, unaff, unc):
        base = (pid, vid, pmid, tot, aff, unaff, unc, None)
        if with_trust:
            conn.execute(
                "INSERT INTO penetrance_data VALUES (?,?,?,?,?,?,?,?,?,?)",
                base + ("quarantine", json.dumps(["population_count"])),
            )
        else:
            conn.execute("INSERT INTO penetrance_data VALUES (?,?,?,?,?,?,?,?)", base)

    pen(10, 1, "111", 3, 2, 1, 0)
    pen(11, 2, "111", 50000, 0, 0, 0)
    for i in range(extra_facts_pmid111):
        pen(100 + i, 1, "111", 3, 2, 1, 0)
    pen(20, 3, "222", 5, 5, 0, 0)
    conn.commit()
    conn.close()


class _Stub:
    """Stub checker returning a normalized flag verdict for each paper."""

    def __init__(self):
        self.calls = 0

    def check(self, payload):
        self.calls += 1
        return normalize_result(
            {
                "paper_verdict": "flag",
                "confidence": 0.4,
                "flags": [
                    {
                        "variant": payload["facts"][0]["variant"],
                        "issue": "looks like a cohort total",
                        "severity": "low",
                    }
                ],
                "summary": "stub summary",
            },
            payload,
        )


def test_check_version_is_stable_and_effort_sensitive():
    a = check_version("azure_ai/gpt-5.6-sol", "xhigh")
    assert a == check_version("azure_ai/gpt-5.6-sol", "xhigh")
    assert a.startswith("pfc1-")
    assert a != check_version("azure_ai/gpt-5.6-sol", "high")
    assert a != check_version("anthropic/claude-sonnet-5", "xhigh")


def test_normalize_infers_verdict_from_flags_and_normalizes_severity():
    payload = {"pmid": "1", "gene": "BRCA2", "facts": [{"n": 1}]}
    r = normalize_result(
        {"flags": [{"variant": "p.V30M", "issue": "cohort total", "severity": "HIGH"}]},
        payload,
    )
    assert r["verdict"] == "flag"
    assert r["n_flagged"] == 1
    assert r["flags"][0]["severity"] == "high"
    assert r["n_facts"] == 1
    assert r["pmid"] == "1"


def test_normalize_ok_when_no_flags_and_bad_confidence():
    payload = {"pmid": "1", "gene": "BRCA2", "facts": []}
    r = normalize_result(
        {"paper_verdict": "weird", "confidence": "n/a", "flags": []}, payload
    )
    assert r["verdict"] == "ok"
    assert r["confidence"] is None


def test_normalize_handles_non_dict_and_clamps_confidence():
    payload = {"pmid": "1", "facts": []}
    assert normalize_result("not json", payload)["verdict"] == "ok"
    assert normalize_result({"confidence": 1.5}, payload)["confidence"] == 1.0
    assert normalize_result({"confidence": -3}, payload)["confidence"] == 0.0


def test_build_prompt_embeds_payload():
    payload = {
        "pmid": "111",
        "gene": "BRCA2",
        "facts": [{"n": 1, "variant": "p.V30M", "counts": {"total_carriers": 3}}],
    }
    prompt = build_paper_check_prompt(payload)
    assert "STRICT JSON" in prompt
    assert "p.V30M" in prompt
    assert "111" in prompt


def test_gather_builds_payloads_and_skips_countless_papers(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)
    conn = sqlite3.connect(db)
    payloads = gather_paper_payloads(conn)
    conn.close()
    by = {p["pmid"]: p for p in payloads}
    assert set(by) == {"111", "222"}  # 333 has no penetrance rows → skipped
    p111 = by["111"]
    assert p111["gene"] == "BRCA2"
    assert len(p111["facts"]) == 2
    f0, f1 = p111["facts"]
    assert f0["variant"] == "p.Val30Met"
    assert f0["counts"]["total_carriers"] == 3
    assert f0["quote"] == "3 carriers, 2 affected"
    assert f0["count_role"] == "per_variant_carrier"
    assert f1["variant"] == "c.68A>G"
    assert f1["count_label"] == "gnomAD AF"
    assert p111["truncated"] is False


def test_gather_truncates_facts(tmp_path):
    db = tmp_path / "t.db"
    _seed(db, extra_facts_pmid111=3)  # 2 + 3 = 5 facts on 111
    conn = sqlite3.connect(db)
    payloads = gather_paper_payloads(conn, max_facts=2)
    conn.close()
    p111 = next(p for p in payloads if p["pmid"] == "111")
    assert len(p111["facts"]) == 2
    assert p111["n_facts_total"] == 5
    assert p111["truncated"] is True


def test_gather_includes_trust_when_present(tmp_path):
    db = tmp_path / "t.db"
    _seed(db, with_trust=True)
    conn = sqlite3.connect(db)
    payloads = gather_paper_payloads(conn)
    conn.close()
    fact = next(p for p in payloads if p["pmid"] == "111")["facts"][0]
    assert fact["trust_tier"] == "quarantine"
    assert fact["trust_reasons"] == ["population_count"]


def test_apply_records_and_is_idempotent(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)
    stub = _Stub()
    s1 = apply_paper_final_check(db, model="m", reasoning_effort="xhigh", checker=stub)
    assert s1["papers"] == 2
    assert s1["checked"] == 2
    assert s1["flag"] == 2
    assert stub.calls == 2

    conn = sqlite3.connect(db)
    rows = conn.execute(
        "SELECT pmid, verdict, model, reasoning_effort, n_flagged, prompt_version "
        "FROM paper_final_check ORDER BY pmid"
    ).fetchall()
    conn.close()
    assert [r[0] for r in rows] == ["111", "222"]
    assert all(r[1] == "flag" and r[2] == "m" and r[3] == "xhigh" for r in rows)
    assert all(r[5] == "pfc1" for r in rows)

    # Re-run replaces (PRIMARY KEY on pmid), never duplicates.
    apply_paper_final_check(db, model="m", reasoning_effort="xhigh", checker=_Stub())
    conn = sqlite3.connect(db)
    n = conn.execute("SELECT COUNT(*) FROM paper_final_check").fetchone()[0]
    conn.close()
    assert n == 2


def test_apply_dry_run_makes_no_llm_call(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)

    class Boom:
        def check(self, payload):
            raise AssertionError("checker must not be called on dry_run")

    s = apply_paper_final_check(
        db, model="m", reasoning_effort="xhigh", checker=Boom(), dry_run=True
    )
    assert s["dry_run"] is True
    assert s["papers"] == 2
    assert "payloads" in s
    conn = sqlite3.connect(db)
    n = conn.execute("SELECT COUNT(*) FROM paper_final_check").fetchone()[0]
    conn.close()
    assert n == 0


def test_apply_records_error_and_continues(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)

    class Raiser:
        def check(self, payload):
            raise RuntimeError("llm down")

    s = apply_paper_final_check(
        db, model="m", reasoning_effort="xhigh", checker=Raiser()
    )
    assert s["error"] == 2
    assert s["checked"] == 2
    conn = sqlite3.connect(db)
    verdicts = [r[0] for r in conn.execute("SELECT verdict FROM paper_final_check")]
    conn.close()
    assert verdicts == ["error", "error"]
