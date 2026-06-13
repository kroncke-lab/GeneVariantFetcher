import json

from scripts.recall_audit import summarize_acquisition_outcome as summary_mod


def _fulltext(body: str) -> str:
    return body + "\n" + ("full text methods results cohort table. " * 40)


def test_acquisition_outcome_reports_selected_and_successful_pmid_recall(tmp_path):
    usable = tmp_path / "123_FULL_CONTEXT.md"
    usable.write_text(_fulltext("# MAIN TEXT\n\nKCNH2 p.Arg1His\n"), encoding="utf-8")
    stub = tmp_path / "456_FULL_CONTEXT.md"
    stub.write_text(
        "# ABSTRACT-ONLY FALLBACK\n\n"
        "> **WARNING:** Full text could not be retrieved for PMID 456.\n"
        "> This document contains only the PubMed abstract and metadata.\n",
        encoding="utf-8",
    )
    worklist_rows = [
        {
            "gene": "KCNH2",
            "pmid": "123",
            "action": "fetch",
            "missing_rows": "3",
            "missing_distinct_variants": "2",
        },
        {
            "gene": "KCNH2",
            "pmid": "456",
            "action": "fetch",
            "missing_rows": "5",
            "missing_distinct_variants": "4",
        },
        {
            "gene": "KCNH2",
            "pmid": "789",
            "action": "refresh_replay",
            "missing_rows": "7",
            "missing_distinct_variants": "6",
        },
    ]
    fetch_rows = [
        {
            "pmid": "123",
            "outcome": "success_via_elsevier_api",
            "path": str(usable),
            "reason": "ok",
        },
        {
            "pmid": "456",
            "outcome": "success",
            "path": str(stub),
            "reason": "stub slipped through",
        },
    ]

    summary, source_override_rows = summary_mod.build_summary(
        gene="KCNH2",
        worklist_rows=worklist_rows,
        fetch_rows=fetch_rows,
        gold_pmids={"123", "456", "789", "999"},
    )

    pmid_recall = summary["pmid_recall"]
    assert pmid_recall["selected_for_fetch_download"] == {
        "pmids": 2,
        "gold_pmids": 4,
        "recall": 0.5,
        "missing_rows": 8,
        "missing_distinct_variants": 6,
    }
    assert pmid_recall["selected_for_source_acquisition_or_binding"]["pmids"] == 3
    assert pmid_recall["fetch_attempted"]["pmids"] == 2
    assert pmid_recall["usable_fulltext_downloaded"] == {
        "pmids": 1,
        "gold_pmids": 4,
        "recall": 0.25,
        "missing_rows": 3,
        "missing_distinct_variants": 2,
    }
    assert source_override_rows == [
        {
            "gene": "KCNH2",
            "pmid": "123",
            "action": "refresh_replay",
            "route": "fetched_fulltext",
            "available_context_path": str(usable),
            "available_context_bytes": usable.stat().st_size,
            "fetch_outcome": "success_via_elsevier_api",
            "fetch_reason": "ok",
            "missing_rows": "3",
            "missing_distinct_variants": "2",
            "notes": "usable full text downloaded/extracted by fetch_paywalled",
        }
    ]


def test_acquisition_outcome_can_report_no_gold_worklist_coverage(tmp_path):
    usable = tmp_path / "123_FULL_CONTEXT.md"
    usable.write_text(_fulltext("# MAIN TEXT\n\nKCNH2 p.Arg1His\n"), encoding="utf-8")
    worklist_rows = [
        {"gene": "KCNH2", "pmid": "123", "action": "fetch"},
        {"gene": "KCNH2", "pmid": "456", "action": "fetch"},
        {"gene": "KCNH2", "pmid": "789", "action": "refresh_replay"},
    ]
    fetch_rows = [
        {"pmid": "123", "outcome": "success", "path": str(usable)},
        {"pmid": "456", "outcome": "failed", "reason": "paywall"},
    ]

    summary, source_override_rows = summary_mod.build_summary(
        gene="KCNH2",
        worklist_rows=worklist_rows,
        fetch_rows=fetch_rows,
        gold_pmids=None,
    )

    assert summary["gold_pmids"] is None
    assert summary["denominator"] == "worklist_pmids"
    assert summary["pmid_recall"]["selected_for_fetch_download"]["recall"] is None
    assert summary["pmid_coverage"]["selected_for_fetch_download"] == {
        "pmids": 2,
        "total_pmids": 3,
        "coverage": 0.666667,
        "missing_rows": 0,
        "missing_distinct_variants": 0,
    }
    assert summary["pmid_coverage"]["usable_fulltext_downloaded"]["pmids"] == 1
    assert source_override_rows[0]["pmid"] == "123"


def test_acquisition_outcome_counts_publisher_api_successes(tmp_path):
    wiley = tmp_path / "111_FULL_CONTEXT.md"
    springer = tmp_path / "222_FULL_CONTEXT.md"
    generic = tmp_path / "333_FULL_CONTEXT.md"
    wiley.write_text(_fulltext("# MAIN TEXT\n\nSCN5A p.Arg18Gln\n"), encoding="utf-8")
    springer.write_text(
        _fulltext("# MAIN TEXT\n\nSCN5A p.Glu1784Lys\n"), encoding="utf-8"
    )
    generic.write_text(
        _fulltext("# MAIN TEXT\n\nSCN5A p.Asn406Ser\n"), encoding="utf-8"
    )
    worklist_rows = [
        {"gene": "SCN5A", "pmid": "111", "action": "fetch"},
        {"gene": "SCN5A", "pmid": "222", "action": "fetch"},
        {"gene": "SCN5A", "pmid": "333", "action": "fetch"},
        {"gene": "SCN5A", "pmid": "555", "action": "fetch"},
    ]
    fetch_rows = [
        {
            "pmid": "111",
            "outcome": "success_via_wiley_api",
            "path": str(wiley),
            "reason": "tdm api",
        },
        {
            "pmid": "222",
            "outcome": "success_via_springer_api",
            "path": str(springer),
            "reason": "openaccess api",
        },
        {
            "pmid": "333",
            "outcome": "success_via_publisher_api",
            "path": str(generic),
            "reason": "generic publisher api",
        },
        {
            "pmid": "555",
            "outcome": "failed",
            "reason": "cloudflare",
        },
    ]

    summary, source_override_rows = summary_mod.build_summary(
        gene="SCN5A",
        worklist_rows=worklist_rows,
        fetch_rows=fetch_rows,
        gold_pmids={"111", "222", "333", "444", "555", "666"},
    )

    assert summary["pmid_recall"]["usable_fulltext_downloaded"]["pmids"] == 3
    assert summary["pmid_recall"]["usable_fulltext_downloaded"]["recall"] == 0.5
    assert summary["usable_fulltext_outcomes"] == {
        "success_via_publisher_api": 1,
        "success_via_springer_api": 1,
        "success_via_wiley_api": 1,
    }
    assert {row["pmid"] for row in source_override_rows} == {"111", "222", "333"}


def test_acquisition_outcome_can_read_interrupted_fetch_output_dir(tmp_path):
    output_dir = tmp_path / "fetch"
    output_dir.mkdir()
    usable = output_dir / "123_FULL_CONTEXT.md"
    usable.write_text(_fulltext("# MAIN TEXT\n\nKCNH2 p.Arg1His\n"), encoding="utf-8")
    result_dir = output_dir / "123"
    result_dir.mkdir()
    (result_dir / "result.json").write_text(
        json.dumps(
            {
                "pmid": "123",
                "publisher": "aha",
                "canonical_full_context_path": str(usable),
                "error": None,
                "notes": [],
            }
        ),
        encoding="utf-8",
    )
    failed_dir = output_dir / "456"
    failed_dir.mkdir()
    (failed_dir / "result.json").write_text(
        json.dumps(
            {
                "pmid": "456",
                "publisher": "oxford",
                "canonical_full_context_path": str(output_dir / "456_FULL_CONTEXT.md"),
                "error": "navigation failed",
                "notes": [],
            }
        ),
        encoding="utf-8",
    )

    rows = summary_mod._load_fetch_output_dir(output_dir)
    summary, source_override_rows = summary_mod.build_summary(
        gene="KCNH2",
        worklist_rows=[
            {"gene": "KCNH2", "pmid": "123", "action": "fetch"},
            {"gene": "KCNH2", "pmid": "456", "action": "fetch"},
        ],
        fetch_rows=rows,
        gold_pmids={"123", "456", "789"},
    )

    assert summary["fetch_outcomes"] == {"error": 1, "success_from_output_dir": 1}
    assert summary["pmid_recall"]["fetch_attempted"]["pmids"] == 2
    assert summary["pmid_recall"]["usable_fulltext_downloaded"]["pmids"] == 1
    assert source_override_rows[0]["pmid"] == "123"


def test_acquisition_outcome_rejects_short_flat_context_from_output_dir(tmp_path):
    output_dir = tmp_path / "fetch"
    output_dir.mkdir()
    short = output_dir / "123_FULL_CONTEXT.md"
    short.write_text("placeholder", encoding="utf-8")
    result_dir = output_dir / "123"
    result_dir.mkdir()
    (result_dir / "result.json").write_text(
        json.dumps(
            {
                "pmid": "123",
                "canonical_full_context_path": str(short),
                "error": None,
                "notes": [],
            }
        ),
        encoding="utf-8",
    )

    rows = summary_mod._load_fetch_output_dir(output_dir)
    summary, source_override_rows = summary_mod.build_summary(
        gene="KCNH2",
        worklist_rows=[{"gene": "KCNH2", "pmid": "123", "action": "fetch"}],
        fetch_rows=rows,
        gold_pmids={"123"},
    )

    assert summary["fetch_outcomes"] == {"empty": 1}
    assert summary["pmid_recall"]["usable_fulltext_downloaded"]["pmids"] == 0
    assert source_override_rows == []


def test_acquisition_outcome_reports_refresh_successful_pmid_recall(tmp_path):
    candidates = tmp_path / "replay_candidates.csv"
    candidates.write_text(
        "\n".join(
            [
                "pmid,source_file,output_file,current_variants,deterministic_variants,reasons",
                "123,a,b,0,0,forced_pmid",
                "456,a,b,0,0,forced_pmid",
                "789,a,b,0,0,forced_pmid",
            ]
        ),
        encoding="utf-8",
    )
    refresh_summary = tmp_path / "refresh_summary.json"
    refresh_summary.write_text(
        json.dumps(
            {
                "candidates_csv": str(candidates),
                "replay": {
                    "attempted": 3,
                    "successful": 1,
                    "failed": 1,
                    "gated": 1,
                    "errors": [{"pmid": "456", "error": "too short"}],
                    "gated_regressions": [{"pmid": "789"}],
                    "gated_explosions": [],
                },
            }
        ),
        encoding="utf-8",
    )
    refresh_attempted, refresh_success = summary_mod._pmids_from_refresh_summary(
        refresh_summary
    )

    summary, _source_override_rows = summary_mod.build_summary(
        gene="KCNH2",
        worklist_rows=[
            {
                "gene": "KCNH2",
                "pmid": "123",
                "action": "fetch",
                "missing_rows": "3",
                "missing_distinct_variants": "2",
            },
            {
                "gene": "KCNH2",
                "pmid": "456",
                "action": "fetch",
                "missing_rows": "5",
                "missing_distinct_variants": "4",
            },
            {
                "gene": "KCNH2",
                "pmid": "789",
                "action": "fetch",
                "missing_rows": "7",
                "missing_distinct_variants": "6",
            },
        ],
        fetch_rows=[],
        refresh_attempted_pmids=refresh_attempted,
        refresh_success_pmids=refresh_success,
        gold_pmids={"123", "456", "789", "999"},
    )

    assert summary["pmid_recall"]["source_refresh_attempted"] == {
        "pmids": 3,
        "gold_pmids": 4,
        "recall": 0.75,
        "missing_rows": 15,
        "missing_distinct_variants": 12,
    }
    assert summary["pmid_recall"]["source_refresh_successful"] == {
        "pmids": 1,
        "gold_pmids": 4,
        "recall": 0.25,
        "missing_rows": 3,
        "missing_distinct_variants": 2,
    }


def test_acquisition_outcome_cli_writes_summary_and_source_override(
    tmp_path, monkeypatch
):
    report = tmp_path / "paper_disagreement_report.csv"
    report.write_text(
        "\n".join(
            [
                "gene,pmid",
                "KCNH2,123",
                "KCNH2,456",
                "SCN5A,999",
            ]
        ),
        encoding="utf-8",
    )
    worklist = tmp_path / "worklist.csv"
    worklist.write_text(
        "\n".join(
            [
                "gene,pmid,action,missing_rows,missing_distinct_variants",
                "KCNH2,123,fetch,3,2",
                "KCNH2,456,manual_or_blocked,5,4",
            ]
        ),
        encoding="utf-8",
    )
    usable = tmp_path / "123_FULL_CONTEXT.md"
    usable.write_text(_fulltext("# MAIN TEXT\n\nKCNH2 p.Arg1His\n"), encoding="utf-8")
    fetch_summary = tmp_path / "summary.json"
    fetch_summary.write_text(
        json.dumps(
            [
                {
                    "pmid": "123",
                    "outcome": "success",
                    "canonical_path": str(usable),
                    "reason": "ok",
                }
            ]
        ),
        encoding="utf-8",
    )
    out = tmp_path / "out.json"
    override = tmp_path / "source_override.csv"

    monkeypatch.setattr(
        "sys.argv",
        [
            "summarize_acquisition_outcome.py",
            "--gene",
            "KCNH2",
            "--worklist",
            str(worklist),
            "--fetch-summary",
            str(fetch_summary),
            "--report",
            str(report),
            "--out",
            str(out),
            "--source-override-out",
            str(override),
        ],
    )
    assert summary_mod.main() == 0

    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["gold_pmids"] == 2
    assert payload["pmid_recall"]["selected_for_fetch_download"]["pmids"] == 1
    assert payload["pmid_recall"]["usable_fulltext_downloaded"]["pmids"] == 1
    assert "123" in override.read_text(encoding="utf-8")
