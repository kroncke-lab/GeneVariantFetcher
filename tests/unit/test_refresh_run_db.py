import json

import pytest

import scripts.refresh_run_db as refresh_run_db
from scripts.refresh_run_db import (
    ReplayCandidate,
    _is_variant_explosion,
    replay_candidates,
    select_replay_candidates,
)
from utils.models import ExtractionResult


def _long_fulltext(body: str) -> str:
    return body + "\n" + ("methods results cohort variant table. " * 40)


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
        _long_fulltext(
            "Table S1. Subject Clinical and Genetic Characteristics\n" + rows
        ),
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


def test_selects_pdf_linearized_table_lift(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "30758498"
    linearized_table = """
eTable 1. LQT1 Mutations or Rare Variants
Mutation
site
Site
N
Female
(n)
Proband
(n)
Mean QTc
(proband)
Syncope
(n)
CA/VF
(n)
c.521G>A  p.R147H
MS
non-pore MS
2
1
1
444
1
0
c.502G>A  p.G168R
MS
non-pore MS
8
7
2
500
0
0
c.520C>T  p.R174C
C-loop
non-pore MS
7
4
3
466
1
0
"""
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        _long_fulltext(linearized_table),
        encoding="utf-8",
    )
    (extraction_dir / f"KCNQ1_PMID_{pmid}.json").write_text(
        json.dumps(
            {"variants": [], "extraction_metadata": {"total_variants_found": 0}}
        ),
        encoding="utf-8",
    )

    candidates = select_replay_candidates(
        gene="KCNQ1",
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=3,
        min_deterministic_lift=3,
        deterministic_lift_ratio=1.2,
        include_source_newer=False,
        replay_missing_fingerprint=False,
        replay_unbound_source=False,
    )

    assert [c.pmid for c in candidates] == [pmid]
    assert candidates[0].deterministic_variants == 3
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
        _long_fulltext("# MAIN TEXT\n\nKCNQ1 variant table\n"),
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
        refresh_run_db,
        "deterministic_variant_list",
        lambda *_args: [{"cdna_notation": f"c.{i}A>G"} for i in range(412)],
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


class _StubExtractor:
    """Test double for ExpertExtractor that returns a preconfigured result."""

    def __init__(self, *, models=None, tier_threshold, fulltext_dir):
        self.models = models
        del tier_threshold, fulltext_dir

    def _set_result(self, result):
        self._result = result

    def extract(self, paper):
        return self._result


def _make_replay_setup(tmp_path, *, prior_variant_count: int, new_variant_count: int):
    """Wire a minimal replay scenario; returns (candidate, harvest_dir, backup_dir, prior_payload)."""
    gene = "KCNQ1"
    pmid = "12345678"
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    backup_dir = tmp_path / "refresh_X" / "extraction_json_backup"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    source_file = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
    source_file.write_text(
        "# MAIN TEXT\n\nKCNQ1 c.123A>G p.Lys41Arg\n" * 40,
        encoding="utf-8",
    )

    prior_payload = {
        "variants": [
            {"protein_notation": f"p.X{i}Y"} for i in range(prior_variant_count)
        ],
        "extraction_metadata": {"total_variants_found": prior_variant_count},
    }
    output_file = extraction_dir / f"{gene}_PMID_{pmid}.json"
    output_file.write_text(json.dumps(prior_payload), encoding="utf-8")

    candidate = ReplayCandidate(
        pmid=pmid,
        source_file=source_file,
        output_file=output_file,
        current_variants=prior_variant_count,
        deterministic_variants=new_variant_count,
        reasons=["forced_pmid"],
    )

    new_payload = {
        "variants": [
            {"protein_notation": f"p.Z{i}W"} for i in range(new_variant_count)
        ],
        "extraction_metadata": {"total_variants_found": new_variant_count},
    }
    extraction_result = ExtractionResult(
        pmid=pmid,
        success=True,
        extracted_data=new_payload,
        model_used="stub-extractor",
    )
    return {
        "gene": gene,
        "pmid": pmid,
        "candidate": candidate,
        "harvest_dir": harvest_dir,
        "backup_dir": backup_dir,
        "output_file": output_file,
        "prior_payload": prior_payload,
        "extraction_result": extraction_result,
    }


def _install_stub_extractor(monkeypatch, extraction_result):
    """Replace ExpertExtractor + is_usable_fulltext_source for one test."""
    stub = _StubExtractor(tier_threshold=0, fulltext_dir="")
    stub._set_result(extraction_result)
    captured = {}

    def _factory(*, models=None, tier_threshold, fulltext_dir):
        captured["models"] = models
        stub.models = models
        return stub

    monkeypatch.setattr(refresh_run_db, "ExpertExtractor", _factory)
    monkeypatch.setattr(
        refresh_run_db,
        "is_usable_fulltext_source",
        lambda _path: True,
    )
    monkeypatch.setattr(
        refresh_run_db,
        "_source_metadata",
        lambda _path: {"source_file": "stub"},
    )
    return captured


def test_replay_gates_per_pmid_regression(tmp_path, monkeypatch):
    """When re-extraction yields fewer variants, the backup must be restored."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=10, new_variant_count=3)
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
    )

    assert stats["attempted"] == 1
    assert stats["successful"] == 0
    assert stats["gated"] == 1
    assert stats["gated_regressions"] == [
        {
            "pmid": setup["pmid"],
            "prior_variant_count": 10,
            "new_variant_count": 3,
            "delta": -7,
        }
    ]
    # File on disk must match the prior (10-variant) payload, not the new (3) payload.
    payload = json.loads(setup["output_file"].read_text())
    assert len(payload["variants"]) == 10
    # Backup must still be present.
    assert (setup["backup_dir"] / setup["output_file"].name).exists()


def test_replay_accepts_fewer_but_paired_over_cdna_only_prior(tmp_path, monkeypatch):
    """A cleaner paired (cDNA+protein) re-extraction must not be regression-gated
    by a stale, over-counted cDNA-only prior, even with fewer total rows."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=10, new_variant_count=6)
    # Prior: 10 rows over 6 distinct codons (4 are exact duplicates — the stale
    # over-count) -> count 10, quality 10, 6 distinct positions.
    base = [100, 103, 106, 109, 112, 115]  # codons 34..39
    prior_payload = {
        "variants": [{"cdna_notation": f"c.{n}A>G"} for n in base + base[:4]],
        "extraction_metadata": {"total_variants_found": 10},
    }
    setup["output_file"].write_text(json.dumps(prior_payload), encoding="utf-8")
    # New: the same 6 distinct variants, now PAIRED -> count 6, quality 12, and it
    # covers all 6 prior positions (coverage 1.0), so it is a clean quality lift.
    new_payload = {
        "variants": [
            {"cdna_notation": f"c.{n}A>G", "protein_notation": f"p.X{(n + 2) // 3}Y"}
            for n in base
        ],
        "extraction_metadata": {"total_variants_found": 6},
    }
    setup["extraction_result"].extracted_data = new_payload
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
    )

    assert stats["successful"] == 1
    assert stats["gated"] == 0
    # The cleaner 6-variant paired payload must now be on disk.
    payload = json.loads(setup["output_file"].read_text())
    assert len(payload["variants"]) == 6
    assert all(v.get("protein_notation") for v in payload["variants"])


def test_replay_gates_pairing_that_drops_positions(tmp_path, monkeypatch):
    """No-gold safety net: a re-extraction that pairs rows but DROPS distinct
    positions is gated even when its quality score held — pairing must not mask a
    recall loss on a gene with no gold rescore to catch it."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=10, new_variant_count=6)
    # Prior: 10 cDNA-only at 10 distinct codons.
    prior_payload = {
        "variants": [{"cdna_notation": f"c.{100 + 3 * i}A>G"} for i in range(10)],
        "extraction_metadata": {"total_variants_found": 10},
    }
    setup["output_file"].write_text(json.dumps(prior_payload), encoding="utf-8")
    # New: 6 PAIRED rows (quality 12 >= prior 10) but only 4 distinct codons
    # (2 duplicated) -> covers 4/10 = 40% of prior positions, below the 85% gate.
    new_payload = {
        "variants": [
            {"cdna_notation": f"c.{n}A>G", "protein_notation": f"p.X{(n + 2) // 3}Y"}
            for n in (100, 103, 106, 109, 100, 103)
        ],
        "extraction_metadata": {"total_variants_found": 6},
    }
    setup["extraction_result"].extracted_data = new_payload
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
    )

    assert stats["gated"] == 1
    assert stats["successful"] == 0
    # The lossy paired re-extraction is rejected; the prior (10) is preserved.
    assert len(json.loads(setup["output_file"].read_text())["variants"]) == 10


def test_replay_accepts_faithful_re_extraction_over_overcounted_prior(
    tmp_path, monkeypatch
):
    """The 30758498 class: the prior is OVER-counted vs the trusted deterministic
    TABLE parse. A re-extraction faithful to that table parse is accepted even with
    far fewer rows than the prior — the table is the structural ground truth, so
    the prior's excess is droppable. This is the gold-free recovery the safety
    guard must still allow."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=18, new_variant_count=6)
    # Prior: 18 cDNA-only rows (over-counted, unpaired).
    prior_payload = {
        "variants": [{"cdna_notation": f"c.{100 + 3 * i}A>G"} for i in range(18)],
        "extraction_metadata": {"total_variants_found": 18},
    }
    setup["output_file"].write_text(json.dumps(prior_payload), encoding="utf-8")
    # The deterministic TABLE parse found 6 variants at 6 codons (the real table).
    det_codons = [34, 35, 36, 37, 38, 39]
    setup["candidate"].deterministic_positions = frozenset(det_codons)
    # Re-extraction returns those 6 PAIRED, covering 100% of the table -> accept.
    new_payload = {
        "variants": [
            {"cdna_notation": f"c.{3 * c - 2}A>G", "protein_notation": f"p.X{c}Y"}
            for c in det_codons
        ],
        "extraction_metadata": {"total_variants_found": 6},
    }
    setup["extraction_result"].extracted_data = new_payload
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
    )

    assert stats["successful"] == 1
    assert stats["gated"] == 0
    payload = json.loads(setup["output_file"].read_text())
    assert len(payload["variants"]) == 6
    assert all(v.get("protein_notation") for v in payload["variants"])


def test_variant_positions_unifies_cdna_and_protein():
    # cDNA implied codon ((nt+2)//3) and protein residue collapse to one position.
    assert refresh_run_db._variant_positions([{"cdna_notation": "c.153C>A"}]) == {51}
    assert refresh_run_db._variant_positions([{"protein_notation": "p.Y51X"}]) == {51}
    assert refresh_run_db._variant_positions(
        [{"cdna_notation": "c.153C>A", "protein_notation": "p.Tyr51Ter"}]
    ) == {51}


def test_position_coverage_detects_dropped_variants():
    prior = [{"protein_notation": f"p.X{r}Y"} for r in (10, 20, 30, 40)]
    full = [{"protein_notation": f"p.X{r}Z"} for r in (10, 20, 30, 40)]
    half = [{"protein_notation": f"p.X{r}Z"} for r in (10, 20)]
    assert refresh_run_db._position_coverage(full, prior) == 1.0
    assert refresh_run_db._position_coverage(half, prior) == 0.5
    assert refresh_run_db._position_coverage([], []) == 1.0  # nothing to lose


def test_selects_quality_lift_without_count_lift(tmp_path, monkeypatch):
    """The stale-over-count case: the deterministic parse has FEWER rows than the
    current extraction but better pairing while covering its positions -> selected
    by the quality-lift trigger that the count-only selector would miss."""
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()
    pmid = "30758498"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        _long_fulltext("# MAIN TEXT\n\nKCNQ1 variant table\n"), encoding="utf-8"
    )
    base = [100, 103, 106, 109, 112, 115]  # 6 distinct codons
    # Current DB: 12 cDNA-only rows over those 6 codons (over-counted, quality 6).
    (extraction_dir / f"KCNQ1_PMID_{pmid}.json").write_text(
        json.dumps(
            {
                "variants": [{"cdna_notation": f"c.{n}A>G"} for n in base + base],
                "extraction_metadata": {"total_variants_found": 12},
            }
        ),
        encoding="utf-8",
    )
    # Deterministic parse: 6 PAIRED rows over the same 6 codons (quality 12, fewer
    # rows) -> no count lift, but a clean quality lift covering all positions.
    monkeypatch.setattr(
        refresh_run_db,
        "deterministic_variant_list",
        lambda *_args: [
            {"cdna_notation": f"c.{n}A>G", "protein_notation": f"p.X{(n + 2) // 3}Y"}
            for n in base
        ],
    )

    candidates = select_replay_candidates(
        gene="KCNQ1",
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
    reasons = candidates[0].reasons
    assert any(r.startswith("deterministic_quality_lift") for r in reasons)
    # Raw-count lift must NOT fire (6 deterministic < 12 current).
    assert not any(r.startswith("deterministic_parser_lift") for r in reasons)


def test_replay_gates_regression_when_prior_count_is_metadata_only(
    tmp_path, monkeypatch
):
    """Older JSONs may only have total_variants_found; those still need gating."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=10, new_variant_count=3)
    prior_payload = {
        "extraction_metadata": {
            "total_variants_found": 10,
            "notes": "metadata-only prior",
        }
    }
    setup["output_file"].write_text(json.dumps(prior_payload), encoding="utf-8")
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
    )

    assert stats["successful"] == 0
    assert stats["gated"] == 1
    assert stats["gated_regressions"] == [
        {
            "pmid": setup["pmid"],
            "prior_variant_count": 10,
            "new_variant_count": 3,
            "delta": -7,
        }
    ]
    payload = json.loads(setup["output_file"].read_text())
    assert "variants" not in payload
    assert payload["extraction_metadata"]["total_variants_found"] == 10


def test_replay_accepts_per_pmid_improvement(tmp_path, monkeypatch):
    """When re-extraction yields more variants, the new JSON must be written."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=3, new_variant_count=10)
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
    )

    assert stats["successful"] == 1
    assert stats["gated"] == 0
    assert stats["gated_regressions"] == []
    payload = json.loads(setup["output_file"].read_text())
    assert len(payload["variants"]) == 10


def test_no_gate_regressions_disables_per_pmid_acceptance(tmp_path, monkeypatch):
    """--no-gate-regressions must overwrite even when the new JSON has fewer variants."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=10, new_variant_count=3)
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
        gate_regressions=False,
    )

    assert stats["successful"] == 1
    assert stats["gated"] == 0
    payload = json.loads(setup["output_file"].read_text())
    assert len(payload["variants"]) == 3  # overwritten with fewer variants


def test_replay_gates_variant_explosion(tmp_path, monkeypatch):
    """A suspicious count blow-up must restore the backup and be recorded."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=10, new_variant_count=600)
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
    )

    assert stats["successful"] == 0
    assert stats["gated"] == 1
    assert stats["gated_explosions"] == [
        {
            "pmid": setup["pmid"],
            "prior_variant_count": 10,
            "new_variant_count": 600,
            "delta": 590,
        }
    ]
    # Backup (10-variant) payload must be restored, not the 600-variant explosion.
    payload = json.loads(setup["output_file"].read_text())
    assert len(payload["variants"]) == 10
    assert (setup["backup_dir"] / setup["output_file"].name).exists()


def test_replay_accepts_large_legitimate_recovery(tmp_path, monkeypatch):
    """A real supplement recovery (5 -> 28) is well under the floor and accepted."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=5, new_variant_count=28)
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
    )

    assert stats["successful"] == 1
    assert stats["gated"] == 0
    assert stats["gated_explosions"] == []
    payload = json.loads(setup["output_file"].read_text())
    assert len(payload["variants"]) == 28


def test_no_gate_explosions_disables_explosion_gate(tmp_path, monkeypatch):
    """--no-gate-explosions must accept the larger extraction."""
    setup = _make_replay_setup(tmp_path, prior_variant_count=10, new_variant_count=600)
    _install_stub_extractor(monkeypatch, setup["extraction_result"])

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
        gate_explosions=False,
    )

    assert stats["successful"] == 1
    assert stats["gated"] == 0
    payload = json.loads(setup["output_file"].read_text())
    assert len(payload["variants"]) == 600


def test_replay_model_override_is_passed_to_extractor(tmp_path, monkeypatch):
    setup = _make_replay_setup(tmp_path, prior_variant_count=1, new_variant_count=2)
    captured = _install_stub_extractor(monkeypatch, setup["extraction_result"])
    replay_models = ["anthropic/claude-sonnet-4-6", "anthropic/claude-opus-4-7"]

    stats = replay_candidates(
        candidates=[setup["candidate"]],
        gene=setup["gene"],
        harvest_dir=setup["harvest_dir"],
        backup_dir=setup["backup_dir"],
        tier_threshold=0,
        dry_run=False,
        replay_models=replay_models,
    )

    assert captured["models"] == replay_models
    assert stats["replay_models"] == replay_models


def test_split_model_args_accepts_repeatable_and_comma_forms():
    assert refresh_run_db._split_model_args(
        [
            "anthropic/claude-sonnet-4-6, anthropic/claude-opus-4-7",
            "azure_ai/gpt-5.3-codex-1",
        ]
    ) == [
        "anthropic/claude-sonnet-4-6",
        "anthropic/claude-opus-4-7",
        "azure_ai/gpt-5.3-codex-1",
    ]


def test_is_variant_explosion_thresholds():
    """Pure-helper boundary checks for the explosion gate."""
    # Recovery from a 0/None prior is never an explosion.
    assert (
        _is_variant_explosion(None, 999, ratio=10.0, min_new=400, min_delta=300)
        is False
    )
    assert (
        _is_variant_explosion(0, 999, ratio=10.0, min_new=400, min_delta=300) is False
    )
    # Legitimate modest recovery: large multiple but small absolute -> not gated.
    assert _is_variant_explosion(5, 28, ratio=10.0, min_new=400, min_delta=300) is False
    # Large absolute but small multiple -> not gated.
    assert (
        _is_variant_explosion(100, 450, ratio=10.0, min_new=400, min_delta=300) is False
    )
    # All three conditions met -> gated.
    assert (
        _is_variant_explosion(10, 600, ratio=10.0, min_new=400, min_delta=300) is True
    )


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


def test_load_source_override_csv_filters_to_usable_refresh_rows(tmp_path):
    usable = tmp_path / "123_FULL_CONTEXT.md"
    usable.write_text(
        _long_fulltext("# MAIN TEXT\n\nKCNH2 p.Arg1His\n"), encoding="utf-8"
    )
    abstract = tmp_path / "456_FULL_CONTEXT.md"
    abstract.write_text(
        "# ABSTRACT-ONLY FALLBACK\n\n"
        "> **WARNING:** Full text could not be retrieved for PMID 456.\n"
        "> This document contains only the PubMed abstract and metadata.\n",
        encoding="utf-8",
    )
    csv_path = tmp_path / "source_override.csv"
    csv_path.write_text(
        "\n".join(
            [
                "pmid,action,available_context_path",
                f"123,refresh_replay,{usable}",
                f"456,refresh_replay,{abstract}",
                f"789,fetch,{usable}",
            ]
        ),
        encoding="utf-8",
    )

    selected = refresh_run_db.load_source_override_csv(csv_path)

    assert selected == {"123": usable}


def test_stage_extractions_uses_copy_and_leaves_active_json(tmp_path, monkeypatch):
    run_dir = tmp_path / "run"
    harvest_dir = run_dir / "pmc_fulltext"
    extraction_dir = run_dir / "extractions"
    harvest_dir.mkdir(parents=True)
    extraction_dir.mkdir()

    pmid = "12345678"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# MAIN TEXT\n\nKCNH2 p.Arg1His\n",
        encoding="utf-8",
    )
    active_json = extraction_dir / f"KCNH2_PMID_{pmid}.json"
    active_json.write_text(
        json.dumps(
            {
                "variants": [{"protein_notation": "p.Arg1His"}],
                "extraction_metadata": {"total_variants_found": 1},
            }
        ),
        encoding="utf-8",
    )
    pmids_file = tmp_path / "pmids.txt"
    pmids_file.write_text(f"{pmid}\n", encoding="utf-8")

    monkeypatch.setattr(
        "sys.argv",
        [
            "refresh_run_db.py",
            "--gene",
            "KCNH2",
            "--run-dir",
            str(run_dir),
            "--pmids-file",
            str(pmids_file),
            "--only-forced-pmids",
            "--stage-extractions",
            "--skip-recovery",
            "--dry-run",
        ],
    )

    assert refresh_run_db.main() == 0
    refresh_dirs = list(run_dir.glob("refresh_*"))
    assert len(refresh_dirs) == 1
    summary = json.loads((refresh_dirs[0] / "refresh_summary.json").read_text())
    staged_dir = refresh_dirs[0] / "staged_extractions"

    assert summary["staged_extractions"] is True
    assert summary["original_extraction_dir"] == str(extraction_dir)
    assert summary["extraction_dir"] == str(staged_dir)
    assert (staged_dir / active_json.name).read_text(encoding="utf-8") == (
        active_json.read_text(encoding="utf-8")
    )


def test_fold_on_disk_supplements_grows_full_context_and_marks_pmid(tmp_path):
    from harvesting.supplement_fold import FOLD_BEGIN

    pmid = "30758498"
    (tmp_path / f"{pmid}_FULL_CONTEXT.md").write_text(
        _long_fulltext("# MAIN\n\nKCNQ1 LQTS cohort body"), encoding="utf-8"
    )
    supp = tmp_path / f"{pmid}_supplements"
    supp.mkdir()
    (supp / "tableS1.csv").write_text("variant,carriers\nc.1A>G,3\n", encoding="utf-8")

    folded = refresh_run_db.fold_on_disk_supplements(tmp_path)

    assert pmid in folded
    fc = (tmp_path / f"{pmid}_FULL_CONTEXT.md").read_text(encoding="utf-8")
    assert FOLD_BEGIN in fc and "c.1A>G,3" in fc


def test_discover_prefers_full_context_when_condensed_is_stale(tmp_path):
    from harvesting.supplement_fold import FOLD_BEGIN, FOLD_END

    pmid = "12345678"
    folded_block = f"\n\n{FOLD_BEGIN}\nsupp table c.1A>G,3\n{FOLD_END}\n"
    (tmp_path / f"{pmid}_FULL_CONTEXT.md").write_text(
        _long_fulltext("# MAIN\n\nKCNQ1 body") + folded_block, encoding="utf-8"
    )
    # Condensed source predates the fold (no sentinel) -> stale -> must be skipped
    (tmp_path / f"{pmid}_DATA_ZONES.md").write_text(
        _long_fulltext("# zones\n\nKCNQ1 condensed"), encoding="utf-8"
    )

    src = refresh_run_db.discover_source_files(tmp_path)
    assert src[pmid].name.endswith("_FULL_CONTEXT.md")

    # Once the condensed form carries the sentinel (regenerated), it wins again.
    (tmp_path / f"{pmid}_DATA_ZONES.md").write_text(
        _long_fulltext("# zones\n\nKCNQ1 condensed") + folded_block, encoding="utf-8"
    )
    src2 = refresh_run_db.discover_source_files(tmp_path)
    assert src2[pmid].name.endswith("_DATA_ZONES.md")
