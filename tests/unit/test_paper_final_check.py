"""Unit tests for the per-paper final check (sniff test).

Covers the pure helpers (version, normalize, gather) and the SQLite adapter with
an injected stub checker — no network, no live model.
"""

import json
import sqlite3

import pytest

from pipeline.paper_final_check import (
    apply_paper_final_check,
    build_paper_check_prompt,
    check_version,
    gather_paper_payloads,
    normalize_result,
    payload_content_hash,
)


def _seed(path, *, with_trust=False, extra_facts_pmid111=0):
    trust_cols = ", trust_tier TEXT, trust_reasons TEXT" if with_trust else ""
    conn = sqlite3.connect(path)
    conn.executescript(
        f"""
        CREATE TABLE papers (pmid TEXT PRIMARY KEY, title TEXT, gene_symbol TEXT);
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY, gene_symbol TEXT,
            cdna_notation TEXT, protein_notation TEXT, genomic_position TEXT,
            structural_description TEXT, variant_class TEXT
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
        "INSERT INTO variants VALUES (?,?,?,?,?,?,?)",
        [
            (1, "BRCA2", None, "p.Val30Met", None, None, "missense"),
            (2, "BRCA2", "c.68A>G", None, None, None, "missense"),
            (3, "BRCA2", None, "p.Arg100Ter", None, None, "nonsense"),
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
    assert a.startswith("pfc7-")
    assert a != check_version("azure_ai/gpt-5.6-sol", "high")
    assert a != check_version("anthropic/claude-sonnet-5", "xhigh")


def test_normalize_flag_response_normalizes_severity():
    payload = {
        "pmid": "1",
        "gene": "BRCA2",
        "facts": [{"fact_id": 11, "variant_id": 101}],
        "captured_fact_index": [
            {"fact_id": 11, "variant_id": 101},
            {"fact_id": 12, "variant_id": 102},
        ],
    }
    r = normalize_result(
        {
            "paper_verdict": "flag",
            "confidence": 0.6,
            "flags": [
                {
                    "fact_ids": [11, "11", 12],
                    "variant": "p.V30M",
                    "fields": [
                        "affected",
                        "affected",
                        "total_carriers_observed",
                        "bogus",
                    ],
                    "reason_code": "WRONG_COLUMN",
                    "issue": "cohort total",
                    "severity": "HIGH",
                }
            ],
            "summary": "The extracted count looks like a cohort total.",
        },
        payload,
    )
    assert r["verdict"] == "flag"
    assert r["n_flagged"] == 1
    assert r["flags"][0]["severity"] == "high"
    assert r["flags"][0]["fact_ids"] == [11, 12]
    assert r["flags"][0]["variant_ids"] == [101, 102]
    assert r["flags"][0]["fields"] == ["affected", "total_carriers"]
    assert r["flags"][0]["reason_code"] == "wrong_column"
    assert r["n_facts"] == 1
    assert r["pmid"] == "1"


def test_normalize_accepts_explicit_valid_ok_response():
    payload = {"pmid": "1", "gene": "BRCA2", "facts": []}
    r = normalize_result(
        {
            "paper_verdict": "ok",
            "confidence": 0.9,
            "flags": [],
            "summary": "No suspect extracted counts were found.",
        },
        payload,
    )
    assert r["verdict"] == "ok"
    assert r["confidence"] == 0.9


@pytest.mark.parametrize(
    "raw",
    [
        {},
        {"paper_verdict": "ok"},
        {"paper_verdict": "ok", "confidence": 0.8, "flags": []},
    ],
)
def test_normalize_rejects_empty_or_schema_incomplete_response(raw):
    payload = {"pmid": "1", "facts": []}
    with pytest.raises(ValueError):
        normalize_result(raw, payload)


@pytest.mark.parametrize(
    "raw",
    [
        "not json",
        {
            "paper_verdict": "ok",
            "confidence": 1.5,
            "flags": [],
            "summary": "Explicit but out-of-range confidence.",
        },
        {
            "paper_verdict": "ok",
            "confidence": -3,
            "flags": [],
            "summary": "Explicit but out-of-range confidence.",
        },
    ],
)
def test_normalize_rejects_non_object_or_invalid_confidence(raw):
    with pytest.raises(ValueError):
        normalize_result(raw, {"pmid": "1", "facts": []})


def test_build_prompt_embeds_payload():
    payload = {
        "pmid": "111",
        "gene": "BRCA2",
        "facts": [{"n": 1, "variant": "p.V30M", "counts": {"total_carriers": 3}}],
        "captured_fact_index": [
            {"n": 1, "variant": "p.V30M", "counts": {"total_carriers": 3}}
        ],
    }
    prompt = build_paper_check_prompt(payload)
    assert "STRICT JSON" in prompt
    assert "p.V30M" in prompt
    assert "111" in prompt
    assert "COMPLETE compact" in prompt
    assert "captured_fact_index" in prompt


def test_gather_builds_payloads_and_includes_countless_papers(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)
    conn = sqlite3.connect(db)
    payloads = gather_paper_payloads(conn)
    conn.close()
    by = {p["pmid"]: p for p in payloads}
    assert set(by) == {"111", "222", "333"}
    assert by["333"]["facts"] == []
    assert by["333"]["captured_fact_index"] == []
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
    assert len(p111["captured_fact_index"]) == 2


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
    assert len(p111["captured_fact_index"]) == 5
    assert all(
        "variant" in fact and "counts" in fact for fact in p111["captured_fact_index"]
    )


def test_gather_keeps_complete_index_beyond_sixty_rows(tmp_path):
    db = tmp_path / "t.db"
    _seed(db, extra_facts_pmid111=63)  # 65 facts total on PMID 111
    conn = sqlite3.connect(db)
    payloads = gather_paper_payloads(conn)
    conn.close()

    p111 = next(p for p in payloads if p["pmid"] == "111")
    assert len(p111["facts"]) == 60
    assert p111["truncated"] is True
    assert p111["n_facts_total"] == 65
    assert len(p111["captured_fact_index"]) == 65
    prompt = build_paper_check_prompt(p111)
    assert '"fact_id": 162' in prompt


def test_gather_prioritizes_targeted_table_evidence_beyond_fact_cap(tmp_path):
    """A suspicious late row keeps its header/row evidence in the LLM payload."""
    db = tmp_path / "t.db"
    _seed(db, extra_facts_pmid111=63)
    conn = sqlite3.connect(db)
    conn.executescript(
        """
        CREATE TABLE fact_provenance (
            provenance_id INTEGER PRIMARY KEY,
            variant_id INTEGER,
            pmid TEXT,
            fact_type TEXT,
            fact_value TEXT,
            source_table TEXT,
            source_row TEXT,
            source_column TEXT,
            evidence_quote TEXT,
            count_type TEXT,
            source_layer TEXT
        );
        INSERT INTO penetrance_data VALUES (999,1,'111',43,43,0,0,NULL);
        """
    )
    evidence = (
        "Header: | Gene | Variant | Mean age diagnosis | No. of cases |\n"
        "Target row: | BRCA1 | 180 delA | 43 | 1 |\n"
        "Next row: |  | 185 delAG | 33 | 1 |"
    )
    conn.execute(
        "INSERT INTO fact_provenance VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        (
            1,
            1,
            "111",
            "affected_count",
            "43",
            "Table 5",
            "1",
            "No. of cases",
            evidence,
            "per_variant_carrier",
            "regex_table",
        ),
    )
    conn.commit()

    payload = next(
        value
        for value in gather_paper_payloads(conn, max_facts=60)
        if value["pmid"] == "111"
    )
    conn.close()

    late = next(
        fact for fact in payload["facts"] if fact["counts"].get("affected") == 43
    )
    assert late["fact_id"] == 999
    assert late["table_evidence"][0]["source_table"] == "Table 5"
    assert "Mean age diagnosis" in late["table_evidence"][0]["evidence_quote"]
    assert len(payload["facts"]) == 60
    assert len(payload["captured_fact_index"]) == 66


def test_gather_names_structural_only_variant(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)
    conn = sqlite3.connect(db)
    conn.execute(
        "INSERT INTO variants VALUES (4,'BRCA2',NULL,NULL,NULL,"
        "'deletion of exons 3-5','exon_deletion')"
    )
    conn.execute("INSERT INTO variant_papers VALUES (4,'222','Table 4','','','')")
    conn.execute("INSERT INTO penetrance_data VALUES (30,4,'222',2,2,0,0,NULL)")
    conn.commit()
    payloads = gather_paper_payloads(conn)
    conn.close()

    p222 = next(p for p in payloads if p["pmid"] == "222")
    structural = next(
        fact
        for fact in p222["captured_fact_index"]
        if fact["variant"] == "deletion of exons 3-5"
    )
    assert structural["identifiers"]["variant_class"] == "exon_deletion"
    assert (
        structural["identifiers"]["structural_description"] == "deletion of exons 3-5"
    )


def test_payload_content_hash_is_order_independent_and_content_sensitive():
    first = {
        "pmid": "1",
        "gene": "G",
        "title": "T",
        "facts": [{"n": 1, "variant": "p.A1V", "counts": {"total_carriers": 2}}],
        "captured_fact_index": [
            {"n": 1, "variant": "p.A1V", "counts": {"total_carriers": 2}},
            {"n": 2, "variant": "p.B2C", "counts": {"total_carriers": 3}},
        ],
        "n_facts_total": 2,
    }
    reordered = dict(first)
    reordered["captured_fact_index"] = list(reversed(first["captured_fact_index"]))
    changed = json.loads(json.dumps(first))
    changed["captured_fact_index"][0]["counts"]["total_carriers"] = 9

    assert payload_content_hash(first) == payload_content_hash(reordered)
    assert payload_content_hash(first) != payload_content_hash(changed)


def test_gather_includes_trust_when_present(tmp_path):
    db = tmp_path / "t.db"
    _seed(db, with_trust=True)
    conn = sqlite3.connect(db)
    payloads = gather_paper_payloads(conn)
    conn.close()
    fact = next(p for p in payloads if p["pmid"] == "111")["facts"][0]
    assert fact["fact_id"] == 10
    assert fact["variant_id"] == 1
    assert fact["trust_tier"] == "quarantine"
    assert fact["trust_reasons"] == ["population_count"]


def test_gather_excludes_prior_paper_check_from_reviewer_payload(tmp_path):
    db = tmp_path / "t.db"
    _seed(db, with_trust=True)
    conn = sqlite3.connect(db)
    conn.execute(
        "UPDATE penetrance_data SET trust_tier='quarantine', trust_reasons=? "
        "WHERE penetrance_id=10",
        ('["paper_final_check:wrong_column:high:affected"]',),
    )
    conn.commit()
    payloads = gather_paper_payloads(conn)
    conn.close()
    fact = next(p for p in payloads if p["pmid"] == "111")["facts"][0]
    assert fact["trust_tier"] == "trusted"
    assert "trust_reasons" not in fact


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
    assert s1["papers"] == 3
    assert s1["checked"] == 2
    assert s1["skipped_empty_no_source"] == 1
    assert s1["flag"] == 2
    assert stub.calls == 2

    conn = sqlite3.connect(db)
    rows = conn.execute(
        "SELECT pmid, verdict, model, reasoning_effort, n_flagged, prompt_version "
        "FROM paper_final_check ORDER BY pmid"
    ).fetchall()
    conn.close()
    assert [r[0] for r in rows] == ["111", "222", "333"]
    assert all(r[1] == "flag" and r[2] == "m" and r[3] == "xhigh" for r in rows[:2])
    assert rows[2][1] == "skipped"
    assert all(r[5] == "pfc7" for r in rows)

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
    assert n == 3


def test_db_only_version_skip_invalidates_when_payload_changes(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)
    common = dict(
        model="m",
        reasoning_effort="xhigh",
        source_grounded=False,
        max_source_chars=60000,
        pmids=["111"],
    )
    first = _Stub()
    apply_paper_final_check(db, checker=first, **common)
    conn = sqlite3.connect(db)
    version_before = conn.execute(
        "SELECT check_version FROM paper_final_check WHERE pmid='111'"
    ).fetchone()[0]
    conn.execute(
        "UPDATE penetrance_data SET total_carriers_observed=4 WHERE penetrance_id=10"
    )
    conn.commit()
    conn.close()

    second = _Stub()
    stats = apply_paper_final_check(db, checker=second, **common)
    conn = sqlite3.connect(db)
    version_after = conn.execute(
        "SELECT check_version FROM paper_final_check WHERE pmid='111'"
    ).fetchone()[0]
    conn.close()

    assert second.calls == 1
    assert stats["checked"] == 1
    assert stats["skipped"] == 0
    assert version_after != version_before


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
    assert s["papers"] == 3
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
    assert verdicts == ["error", "error", "skipped"]


def test_empty_db_checker_result_is_error_and_retryable(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)

    class Empty:
        def check(self, payload):
            return {}

    common = dict(
        model="m",
        reasoning_effort="xhigh",
        source_grounded=False,
        max_source_chars=60000,
        pmids=["111"],
    )
    first = apply_paper_final_check(db, checker=Empty(), **common)
    assert first["error"] == 1

    retry = _Stub()
    second = apply_paper_final_check(db, checker=retry, **common)
    assert second["skipped"] == 0
    assert second["flag"] == 1
    assert retry.calls == 1


def test_empty_fact_paper_without_source_is_durably_skipped(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)

    class Boom:
        def check(self, payload):
            raise AssertionError("empty/no-source paper must not call the DB checker")

    stats = apply_paper_final_check(
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        pmids=["333"],
        checker=Boom(),
    )
    assert stats["papers"] == 1
    assert stats["checked"] == 0
    assert stats["skipped"] == 1
    assert stats["skipped_empty_no_source"] == 1
    conn = sqlite3.connect(db)
    row = conn.execute(
        "SELECT verdict, summary FROM paper_final_check WHERE pmid='333'"
    ).fetchone()
    conn.close()
    assert row[0] == "skipped"
    assert "no extracted count facts" in row[1]


# --------------------------------------------------------------------------
# Source-grounded summary + completeness
# --------------------------------------------------------------------------

from pipeline.paper_final_check import (  # noqa: E402
    _select_source_for_summary,
    build_paper_summary_prompt,
    normalize_summary_result,
    resolve_summary_source,
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


def _valid_summary_response(**overrides):
    raw = {
        "paper_verdict": "ok",
        "confidence": 0.9,
        "study_summary": "A source-grounded BRCA2 carrier study.",
        "cohort_sources": [],
        "carrier_groups": [],
        "completeness": {"status": "complete", "notes": "No gaps found."},
        "flags": [],
    }
    raw.update(overrides)
    return raw


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


def _stage_scout(base, pmid, zones=None, abstract=None):
    """Stage the abstract + Data-Scout DATA_ZONES that resolve_summary_source reads."""
    if zones is not None:
        sc = base / "scout_output"
        sc.mkdir(parents=True, exist_ok=True)
        (sc / f"{pmid}_DATA_ZONES.md").write_text(zones)
    if abstract is not None:
        aj = base / "abstract_json"
        aj.mkdir(parents=True, exist_ok=True)
        (aj / f"{pmid}.json").write_text(json.dumps({"abstract": abstract}))


def test_resolve_summary_combines_abstract_and_scout(tmp_path):
    _stage_scout(
        tmp_path,
        "111",
        zones="## Data Zones\ncarrier data " * 20,
        abstract="A study of BRCA2 carriers. " * 10,
    )
    text, kind, path = resolve_summary_source("111", "BRCA2", run_dir=tmp_path)
    assert kind == "abstract_scout"
    assert "BRCA2 carriers" in text and "carrier data" in text
    assert "## Abstract" in text and "## Data-scouted zones" in text
    assert path.endswith("111_DATA_ZONES.md")


def test_resolve_summary_scout_only(tmp_path):
    _stage_scout(tmp_path, "222", zones="scouted carrier rows " * 20)
    text, kind, path = resolve_summary_source("222", "BRCA2", run_dir=tmp_path)
    assert kind == "scout_only"
    assert "scouted carrier rows" in text
    assert path.endswith("222_DATA_ZONES.md")


def test_resolve_summary_abstract_only(tmp_path):
    _stage_scout(tmp_path, "444", abstract="A study of BRCA2 carriers. " * 20)
    text, kind, path = resolve_summary_source("444", "BRCA2", run_dir=tmp_path)
    assert kind == "abstract_only"
    assert "BRCA2 carriers" in text


def test_resolve_summary_abstract_from_json_hint(tmp_path):
    f = tmp_path / "555.json"
    f.write_text(json.dumps({"abstract": "Hinted BRCA2 abstract. " * 20}))
    text, kind, path = resolve_summary_source("555", "BRCA2", source_file_hint=str(f))
    assert kind == "abstract_only"
    assert "Hinted BRCA2 abstract" in text


def test_resolve_summary_never_loads_full_text(tmp_path):
    # A FULL_CONTEXT full-text file must be IGNORED: only abstract + scout zones load.
    ft = tmp_path / "pmc_fulltext"
    ft.mkdir()
    (ft / "666_FULL_CONTEXT.md").write_text("FULL TEXT carrier data " * 50)
    assert resolve_summary_source("666", "BRCA2", run_dir=tmp_path) == (
        None,
        "none",
        None,
    )


def test_resolve_summary_none(tmp_path):
    assert resolve_summary_source("999", "BRCA2", run_dir=tmp_path) == (
        None,
        "none",
        None,
    )


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


def test_select_source_caps_long_repetitive_tables():
    header = "| Variant | Carriers |\n| --- | --- |\n"
    rows = "".join(f"| p.V{i} | {i} |\n" for i in range(200))
    src = "### Results\n" + header + rows
    # Well under the char budget, but the 200-row table is row-sampled, so the
    # model sees only a subset of the rows → the source is TRUNCATED. This must
    # propagate so a "complete" verdict cannot be claimed over a partial view.
    selected, truncated = _select_source_for_summary(src, max_chars=100000)
    assert truncated is True
    assert "omitted for the sniff test" in selected
    assert selected.count("| p.V") <= 25  # head + tail sample, not all 200
    assert "| p.V0 |" in selected and "| p.V199 |" in selected  # head+tail preserved


def test_select_source_short_table_not_truncated():
    # A short table passes through untouched → NOT truncated.
    header = "| Variant | Carriers |\n| --- | --- |\n"
    rows = "".join(f"| p.V{i} | {i} |\n" for i in range(5))
    src = "### Results\n" + header + rows
    selected, truncated = _select_source_for_summary(src, max_chars=100000)
    assert truncated is False
    assert "omitted for the sniff test" not in selected
    assert selected.count("| p.V") == 5  # every row preserved


def test_sampled_long_table_cannot_record_complete():
    # A paper whose only long table was row-sampled must yield truncated=True AND
    # must never persist completeness="complete": the model saw a subset of the
    # table, so full-paper completeness is unprovable. Mirrors the abstract-only
    # truncation guard.
    header = "| Variant | Carriers |\n| --- | --- |\n"
    rows = "".join(f"| p.V{i} | {i} |\n" for i in range(200))
    src = "### Results\n" + header + rows
    selected, truncated = _select_source_for_summary(src, max_chars=100000)
    assert truncated is True

    # The model returns a clean "complete" with no missed carriers; the sampled
    # -table truncation must downgrade it to "unsure".
    result = normalize_summary_result(
        _valid_summary_response(),
        {"pmid": "1", "gene": "G", "facts": []},
        selected,
        truncated=truncated,
    )
    assert result["completeness"]["status"] != "complete"
    assert result["completeness"]["status"] == "unsure"


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
    assert "COMPLETE ``captured_fact_index``" in prompt


def test_zero_fact_paper_with_source_gets_source_grounded_audit(tmp_path):
    db = tmp_path / "t.db"
    _seed(db)
    quote = "Table S1 reports p.Arg77Ter in four affected carriers"
    _stage_scout(tmp_path, "333", zones=(quote + ". ") * 8)
    summary = _StubSummary(
        [
            {
                "variant": "p.Arg77Ter",
                "in_extraction": False,
                "evidence_quote": quote,
                "total_carriers": 4,
                "affected": 4,
            }
        ]
    )

    stats = apply_paper_final_check(
        db,
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        pmids=["333"],
        summary_checker=summary,
    )

    assert summary.calls == 1
    assert stats["papers"] == 1
    assert stats["checked"] == 1
    assert stats["source_grounded"] == 1
    assert stats["missing_carriers"] == 1
    conn = sqlite3.connect(db)
    row = conn.execute(
        "SELECT n_facts, source_grounded, n_missing FROM paper_final_check "
        "WHERE pmid='333'"
    ).fetchone()
    conn.close()
    assert row == (0, 1, 1)


def test_normalize_summary_derives_missing_and_quote_verified():
    payload = {"pmid": "1", "gene": "BRCA2", "facts": []}
    r = normalize_summary_result(
        _valid_summary_response(
            paper_verdict="flag",
            confidence=0.6,
            study_summary="BRCA2 study.",
            carrier_groups=_GROUPS,
            completeness={"status": "gaps"},
            flags=[],
        ),
        payload,
        _SOURCE,
    )
    assert r["source_grounded"] is True
    assert r["completeness"]["n_reported"] == 2
    assert r["completeness"]["n_missing"] == 1
    # Trust is orthogonal to completeness: a grounded miss with NO per-variant
    # trust flag yields verdict 'ok' + completeness 'gaps' (model 'flag' dropped).
    assert r["verdict"] == "ok"
    assert r["completeness"]["status"] == "gaps"
    by = {g["variant"]: g for g in r["carrier_groups"]}
    assert by["p.Arg100Ter"]["status"] == "reported_missing"
    assert by["p.Arg100Ter"]["quote_verified"] == 1  # quote is in the source
    assert by["p.Val30Met"]["status"] == "reported_extracted"


@pytest.mark.parametrize(
    "raw",
    [
        {},
        {"paper_verdict": "ok", "confidence": 0.8, "flags": []},
        {
            "paper_verdict": "ok",
            "confidence": 0.8,
            "flags": [],
            "study_summary": "Study summary.",
            "cohort_sources": [],
            "carrier_groups": [],
        },
    ],
)
def test_normalize_summary_rejects_empty_or_incomplete_response(raw):
    with pytest.raises(ValueError):
        normalize_summary_result(raw, {"pmid": "1", "gene": "G"}, "source")


def test_normalize_summary_accepts_explicit_valid_ok_response():
    result = normalize_summary_result(
        _valid_summary_response(),
        {"pmid": "1", "gene": "G", "facts": []},
        "source",
    )
    assert result["verdict"] == "ok"
    assert result["completeness"]["status"] == "complete"
    assert result["carrier_groups"] == []


def test_apply_source_grounded_writes_summary_and_missing(tmp_path):
    _stage_scout(tmp_path, "111", zones=_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", "")
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
    assert row == (1, "gaps", 1, "pfs12")
    assert groups["p.Arg100Ter"] == "reported_missing"
    assert groups["p.Val30Met"] == "reported_extracted"
    assert verified == 1


def test_apply_source_grounded_version_skip(tmp_path):
    _stage_scout(tmp_path, "111", zones=_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", "")
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


def test_source_grounded_version_skip_invalidates_when_payload_changes(tmp_path):
    _stage_scout(tmp_path, "111", zones=_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", "")
    common = dict(
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        gene="BRCA2",
    )
    apply_paper_final_check(db, summary_checker=_StubSummary(_GROUPS), **common)
    conn = sqlite3.connect(db)
    version_before = conn.execute(
        "SELECT check_version FROM paper_final_check WHERE pmid='111'"
    ).fetchone()[0]
    conn.execute(
        "UPDATE penetrance_data SET total_carriers_observed=4 WHERE penetrance_id=1"
    )
    conn.commit()
    conn.close()

    second = _StubSummary(_GROUPS)
    stats = apply_paper_final_check(db, summary_checker=second, **common)
    conn = sqlite3.connect(db)
    version_after = conn.execute(
        "SELECT check_version FROM paper_final_check WHERE pmid='111'"
    ).fetchone()[0]
    conn.close()

    assert second.calls == 1
    assert stats["checked"] == 1
    assert stats["skipped"] == 0
    assert version_after != version_before


def test_apply_source_grounded_replaces_stale_groups(tmp_path):
    _stage_scout(tmp_path, "111", zones=_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", "")
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
    assert row == (0, "pfc7")


def test_summary_check_version_sensitive_to_budget():
    a = summary_check_version("m", "xhigh", 60000)
    assert a.startswith("pfs12-")
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
            _valid_summary_response(
                paper_verdict="flag",
                confidence=0.5,
                study_summary="ok",
                cohort_sources=[],
                carrier_groups=self.groups,
                completeness={"status": "gaps"},
                flags=[],
            ),
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
        _valid_summary_response(
            carrier_groups=[
                {
                    "variant": "p.X",
                    "in_extraction": "false",
                    "evidence_quote": "p.X seen in five carriers",
                }
            ]
        ),
        {"pmid": "1", "gene": "G", "facts": []},
        "Table 2 reports p.X seen in five carriers.",
    )
    assert r["carrier_groups"][0]["status"] == "reported_missing"
    assert r["completeness"]["n_missing"] == 1  # grounded: quote is in the source


def test_normalize_status_reconciled_with_derived_misses_and_truncation():
    payload = {"pmid": "1", "gene": "G", "facts": []}
    # Model claims "complete" but a GROUNDED miss exists → forced to "gaps".
    r = normalize_summary_result(
        _valid_summary_response(
            carrier_groups=[
                {
                    "variant": "p.X",
                    "in_extraction": False,
                    "evidence_quote": "p.X seen in five carriers",
                }
            ],
            completeness={"status": "complete"},
        ),
        payload,
        "Table 2 reports p.X seen in five carriers.",
    )
    assert r["completeness"]["status"] == "gaps"
    # Truncated source with no miss cannot claim "complete" → "unsure".
    r2 = normalize_summary_result(
        _valid_summary_response(
            carrier_groups=[
                {"variant": "p.Y", "in_extraction": True, "evidence_quote": "q"}
            ],
            completeness={"status": "complete"},
        ),
        payload,
        "src",
        truncated=True,
    )
    assert r2["completeness"]["status"] == "unsure"


def test_ungrounded_missing_is_unsure_not_gaps_and_ok_verdict():
    # A reported-missing group whose quote is NOT in the source is the model's
    # over-report: it must NOT create a gap or a flag (the over-sensitivity fix).
    r = normalize_summary_result(
        _valid_summary_response(
            paper_verdict="flag",
            carrier_groups=[
                {
                    "variant": "p.Z",
                    "in_extraction": False,
                    "evidence_quote": "a carrier claim absent from the source",
                }
            ],
            completeness={"status": "gaps"},
        ),
        {"pmid": "1", "gene": "G", "facts": []},
        "Unrelated source text with no such carrier claim.",
    )
    assert r["completeness"]["n_missing"] == 0
    assert r["completeness"]["n_missing_unverified"] == 1
    assert r["completeness"]["status"] == "unsure"
    assert r["verdict"] == "ok"  # over-eager model 'flag' dropped: nothing grounded


def test_ok_verdict_with_real_flag_is_reconciled_to_flag():
    # A model 'ok' that ships a per-variant flag is reconciled to 'flag'.
    r = normalize_summary_result(
        _valid_summary_response(
            paper_verdict="ok",
            flags=[
                {
                    "variant": "p.A",
                    "issue": "count is a cohort total",
                    "severity": "high",
                }
            ],
        ),
        {"pmid": "1", "gene": "G", "facts": []},
        "source",
    )
    assert r["n_flagged"] == 1
    assert r["verdict"] == "flag"


def test_source_flag_verifies_noncontiguous_header_and_target_rows():
    source = """| Gene | Variant | gnomAD allele count |
|---|---|---|
| OTHER | p.X | 9 |
| KCNH2 | p.Arg744* | 0 |
"""
    evidence = """| Gene | Variant | gnomAD allele count |
|---|---|---|
| KCNH2 | p.Arg744* | 0 |"""
    r = normalize_summary_result(
        _valid_summary_response(
            flags=[
                {
                    "fact_ids": [1],
                    "variant": "p.Arg744*",
                    "fields": ["total_carriers"],
                    "reason_code": "population_count",
                    "evidence_quote": evidence,
                    "issue": "population column",
                    "severity": "high",
                }
            ]
        ),
        {
            "pmid": "1",
            "gene": "KCNH2",
            "facts": [],
            "captured_fact_index": [{"fact_id": 1, "variant_id": 2}],
        },
        source,
    )
    assert r["flags"][0]["evidence_quote_verified"] is True


def test_source_flag_rejects_fragments_stitched_across_tables():
    source = """| Gene | Variant | gnomAD allele count |
|---|---|---|
| OTHER | p.X | 9 |

| Age at diagnosis | Cases |
|---|---|
| 45 | 12 |
"""
    evidence = """| Gene | Variant | gnomAD allele count |
| 45 | 12 |"""
    r = normalize_summary_result(
        _valid_summary_response(
            flags=[
                {
                    "fact_ids": [1],
                    "variant": "p.Arg744*",
                    "fields": ["total_carriers"],
                    "reason_code": "population_count",
                    "evidence_quote": evidence,
                    "issue": "stitched evidence",
                    "severity": "high",
                }
            ]
        ),
        {
            "pmid": "1",
            "gene": "KCNH2",
            "facts": [],
            "captured_fact_index": [{"fact_id": 1, "variant_id": 2}],
        },
        source,
    )
    assert r["flags"][0]["evidence_quote_verified"] is False


def test_source_flag_accepts_title_and_rows_from_same_sampled_table():
    source = """### Table 2

Examples of BRCA1 variants in affected families

| Class | Variant | Carrier |
|---|---|---|
| Pathogenic | c.1del | 18 |
| … 26 more rows (same columns) omitted for the sniff test … |
| Novel | c.99del | 2 |
"""
    evidence = """### Table 2
Examples of BRCA1 variants in affected families
| Class | Variant | Carrier |
| Pathogenic | c.1del | 18 |
| Novel | c.99del | 2 |"""
    r = normalize_summary_result(
        _valid_summary_response(
            flags=[
                {
                    "fact_ids": [1],
                    "variant": "wrong-gene table",
                    "fields": ["total_carriers"],
                    "reason_code": "wrong_gene",
                    "evidence_quote": evidence,
                    "issue": "BRCA1 table captured for BRCA2",
                    "severity": "high",
                }
            ]
        ),
        {
            "pmid": "1",
            "gene": "BRCA2",
            "facts": [],
            "captured_fact_index": [{"fact_id": 1, "variant_id": 2}],
        },
        source,
    )
    assert r["flags"][0]["evidence_quote_verified"] is True


def test_source_flag_rejects_distant_title_attached_to_other_table():
    source = """Examples of BRCA1 variants in affected families

| Class | Variant | Carrier |
|---|---|---|
| Pathogenic | c.1del | 18 |

---

Examples of BRCA2 variants in affected families

| Class | Variant | Carrier |
|---|---|---|
| Novel | c.99del | 2 |
"""
    evidence = """Examples of BRCA1 variants in affected families
| Class | Variant | Carrier |
| Novel | c.99del | 2 |"""
    r = normalize_summary_result(
        _valid_summary_response(
            flags=[
                {
                    "fact_ids": [1],
                    "variant": "stitched title",
                    "fields": ["total_carriers"],
                    "reason_code": "wrong_gene",
                    "evidence_quote": evidence,
                    "issue": "title and row come from different tables",
                    "severity": "high",
                }
            ]
        ),
        {
            "pmid": "1",
            "gene": "BRCA2",
            "facts": [],
            "captured_fact_index": [{"fact_id": 1, "variant_id": 2}],
        },
        source,
    )
    assert r["flags"][0]["evidence_quote_verified"] is False


def test_apply_errored_paper_is_retryable(tmp_path):
    _stage_scout(tmp_path, "111", zones=_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", "")
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


def test_transient_error_preserves_prior_grounded_findings_and_records_attempt(
    tmp_path,
):
    _stage_scout(tmp_path, "111", zones=_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", "")
    common = dict(
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        gene="BRCA2",
    )
    apply_paper_final_check(db, summary_checker=_StubSummary(_GROUPS), **common)
    conn = sqlite3.connect(db)
    prior = conn.execute(
        "SELECT check_version, verdict, n_missing FROM paper_final_check "
        "WHERE pmid='111'"
    ).fetchone()
    conn.execute(
        "UPDATE penetrance_data SET total_carriers_observed=4 WHERE penetrance_id=1"
    )
    conn.commit()
    conn.close()

    class Raiser:
        def check(self, payload, source_text, truncated=False):
            raise RuntimeError("transient reviewer outage")

    stats = apply_paper_final_check(db, summary_checker=Raiser(), **common)
    assert stats["error"] == 1

    conn = sqlite3.connect(db)
    current = conn.execute(
        "SELECT check_version, verdict, n_missing FROM paper_final_check "
        "WHERE pmid='111'"
    ).fetchone()
    groups = conn.execute(
        "SELECT COUNT(*) FROM paper_carrier_groups WHERE pmid='111'"
    ).fetchone()[0]
    attempt = conn.execute(
        "SELECT verdict, summary FROM paper_final_check_attempts "
        "WHERE pmid='111' ORDER BY attempt_id DESC LIMIT 1"
    ).fetchone()
    conn.close()

    assert current == prior
    assert current[1:] == ("ok", 1)
    assert groups == 2
    assert attempt[0] == "error"
    assert "transient reviewer outage" in attempt[1]


def test_empty_summary_checker_result_is_error_and_retryable(tmp_path):
    _stage_scout(tmp_path, "111", zones=_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", "")

    class EmptySummary:
        def check(self, payload, source_text, truncated=False):
            return {}

    common = dict(
        model="m",
        reasoning_effort="xhigh",
        source_grounded=True,
        max_source_chars=60000,
        run_dir=tmp_path,
        gene="BRCA2",
    )
    first = apply_paper_final_check(db, summary_checker=EmptySummary(), **common)
    assert first["error"] == 1

    retry = _StubSummary(_GROUPS)
    second = apply_paper_final_check(db, summary_checker=retry, **common)
    assert second["skipped"] == 0
    assert second["source_grounded"] == 1
    assert retry.calls == 1


def test_apply_stale_groups_cleared_when_degraded_to_db_only(tmp_path):
    _stage_scout(tmp_path, "111", zones=_SOURCE)
    db = tmp_path / "t.db"
    _seed_source(db, "111", "")
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
