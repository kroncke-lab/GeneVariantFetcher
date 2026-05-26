import json

import pytest

import scripts.refresh_run_db as refresh_run_db
from scripts.refresh_run_db import select_replay_candidates


def test_selects_stale_abstract_only_when_fulltext_exists(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "11111111"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# MAIN TEXT\n\nKCNH2 c.1601G>A R534Q\n" * 40,
        encoding="utf-8",
    )
    (extraction_dir / f"KCNH2_PMID_{pmid}.json").write_text(
        json.dumps(
            {
                "variants": [],
                "extraction_metadata": {
                    "abstract_only": True,
                    "notes": "Abstract-only extraction; full text not available",
                },
            }
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="KCNH2",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=5,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
        replay_unbound_source=False,
    )

    assert [c.pmid for c in candidates] == [pmid]
    assert "stale_abstract_only" in candidates[0].reasons


def test_does_not_select_stale_abstract_only_when_only_fallback_exists(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "33333333"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# ABSTRACT-ONLY FALLBACK\n\n"
        "> **WARNING:** Full text could not be retrieved for PMID 33333333.\n"
        "> This document contains only the PubMed abstract and metadata.\n",
        encoding="utf-8",
    )
    (extraction_dir / f"KCNH2_PMID_{pmid}.json").write_text(
        json.dumps(
            {
                "variants": [],
                "extraction_metadata": {
                    "abstract_only": True,
                    "notes": "Abstract-only extraction; full text not available",
                },
            }
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="KCNH2",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=5,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
        replay_unbound_source=False,
    )

    assert candidates == []


def test_selects_deterministic_vertical_table_lift(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "22222222"
    rows = "\n".join(f"RYR2\n{1258 + idx}c>t\nR{420 + idx}W\nNT" for idx in range(6))
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "Table S1. Subject Clinical and Genetic Characteristics\n" + rows,
        encoding="utf-8",
    )
    (extraction_dir / f"RYR2_PMID_{pmid}.json").write_text(
        json.dumps(
            {"variants": [], "extraction_metadata": {"total_variants_found": 0}}
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="RYR2",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=5,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
        replay_unbound_source=False,
    )

    assert [c.pmid for c in candidates] == [pmid]
    assert candidates[0].deterministic_variants == 6
    assert any(r.startswith("deterministic_parser_lift") for r in candidates[0].reasons)


def test_selects_large_absolute_deterministic_lift_even_below_ratio(
    tmp_path, monkeypatch
):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "32893267"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# MAIN TEXT\n\nKCNQ1 variant table\n",
        encoding="utf-8",
    )
    (extraction_dir / f"KCNQ1_PMID_{pmid}.json").write_text(
        json.dumps(
            {
                "variants": [{} for _ in range(356)],
                "extraction_metadata": {"total_variants_found": 356},
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        refresh_run_db, "deterministic_variant_count", lambda *_args: 412
    )

    candidates = select_replay_candidates(
        gene="KCNQ1",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=20,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
        replay_unbound_source=False,
    )

    assert [c.pmid for c in candidates] == [pmid]
    assert candidates[0].current_variants == 356
    assert candidates[0].deterministic_variants == 412
    assert any(
        r.startswith("deterministic_parser_absolute_lift")
        for r in candidates[0].reasons
    )


def test_selects_unbound_source_metadata_when_fulltext_exists(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "20129283"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# MAIN TEXT\n\nSCN5A c.123A>G p.Lys41Arg\n" * 40,
        encoding="utf-8",
    )
    output_file = extraction_dir / f"SCN5A_PMID_{pmid}.json"
    output_file.write_text(
        json.dumps(
            {
                "variants": [{"protein_notation": "p.Lys41Arg"}],
                "extraction_metadata": {
                    "total_variants_found": 1,
                    "source_file": str(output_file),
                },
            }
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="SCN5A",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=20,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
        replay_unbound_source=True,
    )

    assert [c.pmid for c in candidates] == [pmid]
    assert "unbound_source_metadata" in candidates[0].reasons


def test_force_pmids_selects_available_source_without_other_reasons(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "28404607"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# MAIN TEXT\n\nRYR2 p.Arg420Trp\n" * 40,
        encoding="utf-8",
    )
    (extraction_dir / f"RYR2_PMID_{pmid}.json").write_text(
        json.dumps(
            {
                "variants": [{"protein_notation": "p.Arg420Trp"}],
                "extraction_metadata": {
                    "total_variants_found": 1,
                    "model_used": "azure_ai/grok-4-20-reasoning",
                },
            }
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="RYR2",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=20,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
        replay_unbound_source=False,
        force_pmids={pmid},
    )

    assert [c.pmid for c in candidates] == [pmid]
    assert "forced_pmid" in candidates[0].reasons


def test_source_overrides_prefer_report_available_context(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    recovered_dir = tmp_path / "recovered"
    harvest_dir.mkdir()
    extraction_dir.mkdir()
    recovered_dir.mkdir()

    pmid = "32893267"
    stale_source = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
    override_source = recovered_dir / f"{pmid}_FULL_CONTEXT.md"
    stale_source.write_text("# MAIN TEXT\n\nKCNQ1 p.Ala1Val\n" * 40, encoding="utf-8")
    override_source.write_text(
        "# MAIN TEXT\n\nKCNQ1 p.Cys2Trp\nKCNQ1 p.Arg3His\n" * 40,
        encoding="utf-8",
    )
    (extraction_dir / f"KCNQ1_PMID_{pmid}.json").write_text(
        json.dumps(
            {
                "variants": [{"protein_notation": "p.Ala1Val"}],
                "extraction_metadata": {"total_variants_found": 1},
            }
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="KCNQ1",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=20,
        min_deterministic_lift=5,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
        replay_unbound_source=False,
        force_pmids={pmid},
        source_overrides={pmid: override_source},
    )

    assert [c.pmid for c in candidates] == [pmid]
    assert candidates[0].source_file == override_source


def test_load_report_pmids_filters_gene_class_and_missing_rows(tmp_path):
    report = tmp_path / "paper_disagreement_report.csv"
    report.write_text(
        "\n".join(
            [
                "gene,pmid,failure_class,row_recall,missing_rows",
                "KCNQ1,32893267,source_unbound_available,0.041,235",
                "KCNQ1,19716085,source_missing_or_stub,0.039,173",
                "SCN5A,20129283,source_unbound_available,0.066,188",
                "KCNQ1,30758498,source_unbound_available,0.700,93",
            ]
        ),
        encoding="utf-8",
    )

    pmids = refresh_run_db.load_report_pmids(
        report=report,
        gene="KCNQ1",
        failure_classes={"source_unbound_available"},
        min_missing_rows=100,
        max_row_recall=0.5,
    )

    assert pmids == {"32893267"}


def test_load_report_pmids_can_filter_high_existing_row_recall(tmp_path):
    report = tmp_path / "paper_disagreement_report.csv"
    report.write_text(
        "\n".join(
            [
                "gene,pmid,failure_class,row_recall,missing_rows",
                "KCNQ1,32893267,source_unbound_available,0.041,235",
                "KCNQ1,17470695,source_unbound_available,0.839,9",
                "KCNQ1,32145446,source_unbound_available,0.000,1",
            ]
        ),
        encoding="utf-8",
    )

    pmids = refresh_run_db.load_report_pmids(
        report=report,
        gene="KCNQ1",
        failure_classes={"source_unbound_available"},
        min_missing_rows=1,
        max_row_recall=0.5,
    )

    assert pmids == {"32893267", "32145446"}


def test_load_report_available_contexts_filters_and_uses_largest(tmp_path):
    report = tmp_path / "paper_disagreement_report.csv"
    small_context = tmp_path / "small_FULL_CONTEXT.md"
    large_context = tmp_path / "large_FULL_CONTEXT.md"
    other_context = tmp_path / "other_FULL_CONTEXT.md"
    small_context.write_text("# MAIN TEXT\n\nKCNQ1 p.Ala1Val\n" * 40, encoding="utf-8")
    large_context.write_text("# MAIN TEXT\n\nKCNQ1 p.Cys2Trp\n" * 80, encoding="utf-8")
    other_context.write_text("# MAIN TEXT\n\nSCN5A p.Arg3His\n" * 80, encoding="utf-8")
    report.write_text(
        "\n".join(
            [
                "gene,pmid,failure_class,missing_rows,available_context_path",
                f"KCNQ1,32893267,source_unbound_available,235,{small_context}",
                f"KCNQ1,32893267,source_unbound_available,235,{large_context}",
                f"KCNQ1,19716085,source_missing_or_stub,173,{other_context}",
                f"SCN5A,20129283,source_unbound_available,188,{other_context}",
            ]
        ),
        encoding="utf-8",
    )

    contexts = refresh_run_db.load_report_available_contexts(
        report=report,
        gene="KCNQ1",
        failure_classes={"source_unbound_available"},
        min_missing_rows=100,
    )

    assert contexts == {"32893267": large_context}


def test_load_report_available_contexts_can_search_larger_recovered_context(tmp_path):
    report = tmp_path / "paper_disagreement_report.csv"
    run_context = tmp_path / "run" / "20129283_FULL_CONTEXT.md"
    recovered_context = (
        tmp_path / "targeted_table_recovery" / "20129283_FULL_CONTEXT.md"
    )
    run_context.parent.mkdir()
    recovered_context.parent.mkdir()
    run_context.write_text("# MAIN TEXT\n\nSCN5A p.Arg1His\n" * 40, encoding="utf-8")
    recovered_context.write_text(
        "# MAIN TEXT\n\nSCN5A p.Arg1His\nSCN5A p.Arg2His\n" * 80,
        encoding="utf-8",
    )
    report.write_text(
        "\n".join(
            [
                "gene,pmid,failure_class,row_recall,missing_rows,available_context_path",
                f"SCN5A,20129283,source_unbound_available,0.45,188,{run_context}",
            ]
        ),
        encoding="utf-8",
    )

    contexts = refresh_run_db.load_report_available_contexts(
        report=report,
        gene="SCN5A",
        failure_classes={"source_unbound_available"},
        min_missing_rows=1,
        max_row_recall=0.5,
        context_search_roots=[tmp_path],
    )

    assert contexts == {"20129283": recovered_context}


def test_only_forced_pmids_requires_forced_inputs(tmp_path, monkeypatch):
    run_dir = tmp_path / "run"
    (run_dir / "pmc_fulltext").mkdir(parents=True)
    (run_dir / "extractions").mkdir()

    monkeypatch.setattr(
        "sys.argv",
        [
            "refresh_run_db.py",
            "--gene",
            "KCNQ1",
            "--run-dir",
            str(run_dir),
            "--only-forced-pmids",
            "--dry-run",
        ],
    )

    with pytest.raises(SystemExit, match="--only-forced-pmids requires"):
        refresh_run_db.main()
