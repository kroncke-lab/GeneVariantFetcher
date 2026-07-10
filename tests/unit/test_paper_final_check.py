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
    s1 = apply_paper_final_check(
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=False,
        max_source_chars=60000,
        checker=stub,
    )
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
    apply_paper_final_check(
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=False,
        max_source_chars=60000,
        checker=_Stub(),
    )
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
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=False,
        max_source_chars=60000,
        checker=Boom(),
        dry_run=True,
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
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=False,
        max_source_chars=60000,
        checker=Raiser(),
    )
    assert s["error"] == 2
    assert s["checked"] == 2
    conn = sqlite3.connect(db)
    verdicts = [r[0] for r in conn.execute("SELECT verdict FROM paper_final_check")]
    conn.close()
    assert verdicts == ["error", "error"]


# --------------------------------------------------------------------------
# Source-grounded summary + completeness
# --------------------------------------------------------------------------

from pipeline.paper_final_check import (  # noqa: E402
    _select_source_for_summary,
    build_paper_summary_prompt,
    normalize_summary_result,
    resolve_source_text,
    summary_check_version,
)


def _seed_source(path, pmid, source_file):
    conn = sqlite3.connect(path)
    conn.executescript(
        """
        CREATE TABLE papers (pmid TEXT PRIMARY KEY, title TEXT, gene_symbol TEXT);
        CREATE TABLE variants (variant_id INTEGER PRIMARY KEY, gene_symbol TEXT,
            cdna_notation TEXT, protein_notation TEXT, genomic_position TEXT);
        CREATE TABLE variant_papers (variant_id INTEGER, pmid TEXT,
            source_location TEXT, key_quotes TEXT, count_provenance TEXT);
        CREATE TABLE penetrance_data (penetrance_id INTEGER PRIMARY KEY,
            variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
            affected_count INTEGER, unaffected_count INTEGER, uncertain_count INTEGER);
        CREATE TABLE extraction_metadata (extraction_id INTEGER PRIMARY KEY,
            pmid TEXT, source_file TEXT, source_type TEXT);
        """
    )
    conn.execute("INSERT INTO papers VALUES (?, 'A BRCA2 study', 'BRCA2')", (pmid,))
    conn.execute("INSERT INTO variants VALUES (1,'BRCA2',NULL,'p.Val30Met',NULL)")
    conn.execute("INSERT INTO variant_papers (variant_id, pmid) VALUES (1, ?)", (pmid,))
    conn.execute("INSERT INTO penetrance_data VALUES (1,1,?,3,3,0,0)", (pmid,))
    conn.execute(
        "INSERT INTO extraction_metadata VALUES (1, ?, ?, 'fulltext')",
        (pmid, source_file),
    )
    conn.commit()
    conn.close()


_GROUPS = [
    {
        "variant": "p.Val30Met",
        "in_extraction": True,
        "evidence_quote": "Table 2 lists p.Val30Met",
        "total_carriers": 3,
        "disease": "breast cancer",
    },
    {
        "variant": "p.Arg100Ter",
        "in_extraction": False,
        "evidence_quote": "Suppl Table S3 lists p.Arg100Ter with 7 carriers",
        "total_carriers": 7,
        "affected": 7,
        "disease": "breast cancer",
        "cohort_source": "Suppl S3",
        "location_in_paper": "Suppl Table S3",
    },
]

_SOURCE = (
    "Results: Table 2 lists p.Val30Met (3 carriers). "
    "Suppl Table S3 lists p.Arg100Ter with 7 carriers, all affected with "
    "breast cancer. "
) * 3


class _StubSummary:
    def __init__(self, groups):
        self.groups = groups
        self.calls = 0

    def check(self, payload, source_text, truncated=False):
        self.calls += 1
        raw = {
            "paper_verdict": "flag",
            "confidence": 0.5,
            "study_summary": "A BRCA2 family study.",
            "cohort_sources": ["Family F1"],
            "carrier_groups": self.groups,
            "completeness": {"status": "gaps", "notes": "supplement gap"},
            "flags": [],
        }
        return normalize_summary_result(raw, payload, source_text)


def test_resolve_source_prefers_full_context(tmp_path):
    d = tmp_path / "pmc_fulltext"
    d.mkdir()
    (d / "111_FULL_CONTEXT.md").write_text("FULL " + "carrier data " * 50)
    (d / "111_CLEANED.md").write_text("CLEANED " + "x " * 50)
    text, kind, path = resolve_source_text("111", "BRCA2", run_dir=tmp_path)
    assert kind == "fulltext"
    assert text.startswith("FULL")
    assert path.endswith("111_FULL_CONTEXT.md")


def test_resolve_source_uses_source_file_hint(tmp_path):
    d = tmp_path / "run" / "pmc_fulltext"
    d.mkdir(parents=True)
    f = d / "222_CLEANED.md"
    f.write_text("HINTED SOURCE " + "carrier " * 60)
    text, kind, path = resolve_source_text("222", "BRCA2", source_file_hint=str(f))
    assert kind == "fulltext"
    assert text.startswith("HINTED")


def test_resolve_source_oversize_drops_to_cleaned(tmp_path):
    d = tmp_path / "pmc_fulltext"
    d.mkdir()
    (d / "333_FULL_CONTEXT.md").write_text("F" * 5000)
    (d / "333_CLEANED.md").write_text("CLEANED carriers " * 30)
    text, kind, path = resolve_source_text(
        "333", "G", run_dir=tmp_path, oversize_bytes=1000
    )
    assert path.endswith("333_CLEANED.md")
    assert text.startswith("CLEANED")


def test_resolve_source_abstract_json(tmp_path):
    aj = tmp_path / "abstract_json"
    aj.mkdir()
    (aj / "444.json").write_text(
        json.dumps({"abstract": "A study of BRCA2 carriers. " * 20})
    )
    text, kind, path = resolve_source_text("444", "BRCA2", run_dir=tmp_path)
    assert kind == "abstract_only"
    assert "BRCA2 carriers" in text


def test_resolve_source_none(tmp_path):
    assert resolve_source_text("999", "BRCA2", run_dir=tmp_path) == (None, "none", None)


def test_select_source_drops_dead_sections():
    src = "### Results\n5 carriers found.\n### References\n1. Foo et al.\n2. Bar et al."
    selected, truncated = _select_source_for_summary(src, max_chars=10000)
    assert "5 carriers found" in selected
    assert "Foo et al" not in selected
    assert truncated is False


def test_select_source_truncates_but_keeps_tables():
    prose = "narrative sentence. " * 2000
    table = "\n| Variant | Carriers |\n| p.Val30Met | 12 |\n"
    src = "### Results\n" + prose + table
    selected, truncated = _select_source_for_summary(src, max_chars=5000)
    assert truncated is True
    assert "| p.Val30Met | 12 |" in selected
    assert len(selected) <= 5000 + 300


def test_build_summary_prompt_embeds_source_and_extracted():
    payload = {
        "pmid": "1",
        "gene": "BRCA2",
        "facts": [{"n": 1, "variant": "p.Val30Met", "counts": {"total_carriers": 3}}],
    }
    prompt = build_paper_summary_prompt(payload, "SOURCE TEXT HERE with carriers")
    assert "SOURCE TEXT HERE" in prompt
    assert "ALREADY CAPTURED" in prompt
    assert "p.Val30Met" in prompt
    assert "carrier_groups" in prompt
    assert "in_extraction" in prompt


def test_normalize_summary_derives_missing_and_quote_verified():
    payload = {"pmid": "1", "gene": "BRCA2", "facts": []}
    r = normalize_summary_result(
        {
            "paper_verdict": "flag",
            "confidence": 0.6,
            "study_summary": "BRCA2 study.",
            "carrier_groups": _GROUPS,
            "completeness": {"status": "gaps"},
            "flags": [],
        },
        payload,
        _SOURCE,
    )
    assert r["source_grounded"] is True
    assert r["completeness"]["n_reported"] == 2
    assert r["completeness"]["n_missing"] == 1
    by = {g["variant"]: g for g in r["carrier_groups"]}
    assert by["p.Arg100Ter"]["status"] == "reported_missing"
    assert by["p.Arg100Ter"]["quote_verified"] == 1  # quote is in the source
    assert by["p.Val30Met"]["status"] == "reported_extracted"


def test_apply_source_grounded_writes_summary_and_missing(tmp_path):
    ft = tmp_path / "pmc_fulltext"
    ft.mkdir()
    src_file = ft / "111_FULL_CONTEXT.md"
    src_file.write_text(_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", str(src_file))
    stub = _StubSummary(_GROUPS)

    stats = apply_paper_final_check(
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        gene="BRCA2",
        summary_checker=stub,
    )
    assert stats["source_grounded"] == 1
    assert stats["missing_carriers"] == 1
    assert stub.calls == 1

    conn = sqlite3.connect(db)
    row = conn.execute(
        "SELECT source_grounded, completeness_status, n_missing, prompt_version "
        "FROM paper_final_check WHERE pmid='111'"
    ).fetchone()
    groups = dict(
        conn.execute(
            "SELECT variant, status FROM paper_carrier_groups WHERE pmid='111'"
        ).fetchall()
    )
    verified = conn.execute(
        "SELECT quote_verified FROM paper_carrier_groups WHERE variant='p.Arg100Ter'"
    ).fetchone()[0]
    conn.close()
    assert row == (1, "gaps", 1, "pfs1")
    assert groups["p.Arg100Ter"] == "reported_missing"
    assert groups["p.Val30Met"] == "reported_extracted"
    assert verified == 1


def test_apply_source_grounded_version_skip(tmp_path):
    ft = tmp_path / "pmc_fulltext"
    ft.mkdir()
    (ft / "111_FULL_CONTEXT.md").write_text(_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", str(ft / "111_FULL_CONTEXT.md"))
    common = dict(
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        gene="BRCA2",
    )
    apply_paper_final_check(db, summary_checker=_StubSummary(_GROUPS), **common)
    stub2 = _StubSummary(_GROUPS)
    stats2 = apply_paper_final_check(db, summary_checker=stub2, **common)
    assert stats2["skipped"] == 1
    assert stub2.calls == 0
    conn = sqlite3.connect(db)
    n = conn.execute("SELECT COUNT(*) FROM paper_carrier_groups").fetchone()[0]
    conn.close()
    assert n == 2  # no duplication on re-run


def test_apply_source_grounded_replaces_stale_groups(tmp_path):
    ft = tmp_path / "pmc_fulltext"
    ft.mkdir()
    (ft / "111_FULL_CONTEXT.md").write_text(_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", str(ft / "111_FULL_CONTEXT.md"))
    conn = sqlite3.connect(db)
    from pipeline.paper_final_check import ensure_summary_schema

    ensure_summary_schema(conn)
    conn.execute(
        "INSERT INTO paper_carrier_groups (pmid, variant, status) "
        "VALUES ('111','STALE','reported_missing')"
    )
    conn.commit()
    conn.close()

    apply_paper_final_check(
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        gene="BRCA2",
        summary_checker=_StubSummary(_GROUPS),
    )
    conn = sqlite3.connect(db)
    variants = {
        r[0]
        for r in conn.execute(
            "SELECT variant FROM paper_carrier_groups WHERE pmid='111'"
        )
    }
    conn.close()
    assert "STALE" not in variants  # DELETE-before-insert cleared it
    assert variants == {"p.Val30Met", "p.Arg100Ter"}


def test_apply_no_source_degrades_to_db_only(tmp_path):
    db = tmp_path / "t.db"
    _seed_source(db, "111", str(tmp_path / "missing.md"))  # pointer to nonexistent file
    db_stub = _Stub()
    stats = apply_paper_final_check(
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        gene="BRCA2",
        checker=db_stub,
    )
    assert stats["source_grounded"] == 0
    assert stats["source_absent"] == 1
    assert db_stub.calls == 1
    conn = sqlite3.connect(db)
    row = conn.execute(
        "SELECT source_grounded, prompt_version FROM paper_final_check WHERE pmid='111'"
    ).fetchone()
    conn.close()
    assert row == (0, "pfc1")


def test_summary_check_version_sensitive_to_budget():
    a = summary_check_version("m", "xhigh", 60000)
    assert a.startswith("pfs1-")
    assert a != summary_check_version("m", "xhigh", 40000)
    assert a != summary_check_version("m", "high", 60000)


# --- review fixes ---------------------------------------------------------


class _FlakySummary:
    """Raises on the first call, succeeds after (for retry testing)."""

    def __init__(self, groups):
        self.groups = groups
        self.calls = 0

    def check(self, payload, source_text, truncated=False):
        self.calls += 1
        if self.calls == 1:
            raise RuntimeError("transient llm error")
        return normalize_summary_result(
            {
                "paper_verdict": "flag",
                "confidence": 0.5,
                "study_summary": "ok",
                "carrier_groups": self.groups,
                "completeness": {"status": "gaps"},
                "flags": [],
            },
            payload,
            source_text,
            truncated,
        )


def test_as_bool_string_false_is_a_miss():
    from pipeline.paper_final_check import _as_bool

    assert _as_bool("false") is False
    assert _as_bool("true") is True
    assert _as_bool(0) is False
    assert _as_bool(1) is True
    r = normalize_summary_result(
        {
            "carrier_groups": [
                {"variant": "p.X", "in_extraction": "false", "evidence_quote": "q"}
            ]
        },
        {"pmid": "1", "gene": "G", "facts": []},
        "src",
    )
    assert r["carrier_groups"][0]["status"] == "reported_missing"
    assert r["completeness"]["n_missing"] == 1


def test_normalize_status_reconciled_with_derived_misses_and_truncation():
    payload = {"pmid": "1", "gene": "G", "facts": []}
    # Model claims "complete" but a miss exists → forced to "gaps".
    r = normalize_summary_result(
        {
            "carrier_groups": [
                {"variant": "p.X", "in_extraction": False, "evidence_quote": "q"}
            ],
            "completeness": {"status": "complete"},
        },
        payload,
        "src",
    )
    assert r["completeness"]["status"] == "gaps"
    # Truncated source with no miss cannot claim "complete" → "unsure".
    r2 = normalize_summary_result(
        {
            "carrier_groups": [
                {"variant": "p.Y", "in_extraction": True, "evidence_quote": "q"}
            ],
            "completeness": {"status": "complete"},
        },
        payload,
        "src",
        truncated=True,
    )
    assert r2["completeness"]["status"] == "unsure"


def test_resolve_prefers_fulltext_over_abstract_hint(tmp_path):
    aj = tmp_path / "abstract_json"
    aj.mkdir()
    hint = aj / "555.json"
    hint.write_text(json.dumps({"abstract": "abstract text " * 30}))
    ft = tmp_path / "pmc_fulltext"
    ft.mkdir()
    (ft / "555_FULL_CONTEXT.md").write_text("FULL BODY " + "carrier " * 60)
    text, kind, path = resolve_source_text(
        "555", "G", run_dir=tmp_path, source_file_hint=str(hint)
    )
    assert kind == "fulltext"
    assert path.endswith("555_FULL_CONTEXT.md")


def test_resolve_empty_full_falls_to_cleaned(tmp_path):
    d = tmp_path / "pmc_fulltext"
    d.mkdir()
    (d / "666_FULL_CONTEXT.md").write_text("stub")  # <200 chars → unusable
    (d / "666_CLEANED.md").write_text("CLEANED carriers " * 40)
    text, kind, path = resolve_source_text("666", "G", run_dir=tmp_path)
    assert kind == "fulltext"
    assert path.endswith("666_CLEANED.md")
    assert text.startswith("CLEANED")


def test_apply_errored_paper_is_retryable(tmp_path):
    ft = tmp_path / "pmc_fulltext"
    ft.mkdir()
    (ft / "111_FULL_CONTEXT.md").write_text(_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", str(ft / "111_FULL_CONTEXT.md"))
    flaky = _FlakySummary(_GROUPS)
    common = dict(
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        gene="BRCA2",
    )
    s1 = apply_paper_final_check(db, summary_checker=flaky, **common)
    assert s1["error"] == 1
    s2 = apply_paper_final_check(db, summary_checker=flaky, **common)
    assert s2["skipped"] == 0  # error rows are NOT version-skipped
    assert s2["source_grounded"] == 1  # retry succeeded
    assert flaky.calls == 2


def test_apply_stale_groups_cleared_when_degraded_to_db_only(tmp_path):
    ft = tmp_path / "pmc_fulltext"
    ft.mkdir()
    (ft / "111_FULL_CONTEXT.md").write_text(_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", str(ft / "111_FULL_CONTEXT.md"))
    apply_paper_final_check(
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        gene="BRCA2",
        summary_checker=_StubSummary(_GROUPS),
    )
    conn = sqlite3.connect(db)
    assert conn.execute("SELECT COUNT(*) FROM paper_carrier_groups").fetchone()[0] == 2
    conn.close()

    # Re-check with source_grounded off → DB-only → stale groups must be cleared.
    apply_paper_final_check(
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=False,
        max_source_chars=60000,
        gene="BRCA2",
        checker=_Stub(),
    )
    conn = sqlite3.connect(db)
    n = conn.execute("SELECT COUNT(*) FROM paper_carrier_groups").fetchone()[0]
    conn.close()
    assert n == 0
