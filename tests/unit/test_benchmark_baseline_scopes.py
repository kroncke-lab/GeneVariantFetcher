"""Exact-scope baseline selection and benchmark gate behavior."""

from pathlib import Path

import pytest

from benchmarks.curated_extraction_eval.run_benchmark import (
    baseline_scope_key,
    gate_result,
    require_clean_baseline_write,
    require_complete_gene_set,
    select_baseline_profile,
    update_baseline_document,
)


ALL = ["APOE", "KCNH2", "KCNQ1", "RYR2", "SCN5A"]
CARDIAC = ["KCNH2", "KCNQ1", "RYR2", "SCN5A"]
ALL_METRICS = {"recall": {"unique_variants": {"recall": 0.91}}}
CARDIAC_METRICS = {"recall": {"unique_variants": {"recall": 0.92}}}


def _profile(genes, metrics):
    return {"genes": sorted(genes), **metrics}


def _v2():
    return {
        "schema_version": 2,
        "profiles": {
            baseline_scope_key(ALL): _profile(ALL, ALL_METRICS),
            baseline_scope_key(CARDIAC): _profile(CARDIAC, CARDIAC_METRICS),
        },
    }


def test_v2_baseline_selects_only_the_exact_gene_set():
    baseline, note = select_baseline_profile(_v2(), CARDIAC, ALL)

    assert note is None
    assert baseline["recall"]["unique_variants"]["recall"] == 0.92


def test_v2_baseline_rejects_changed_fixture_with_same_gene_scope():
    document = _v2()
    document["profiles"][baseline_scope_key(CARDIAC)]["fixture_sha256"] = "old"

    baseline, note = select_baseline_profile(
        document, CARDIAC, ALL, fixture_hash="current"
    )

    assert baseline is None
    assert "fixture fingerprint differs" in note


def test_v2_all_gene_profile_does_not_match_an_unregistered_subset():
    document = {
        "schema_version": 2,
        "profiles": {baseline_scope_key(ALL): _profile(ALL, ALL_METRICS)},
    }

    baseline, note = select_baseline_profile(document, CARDIAC, ALL)

    assert baseline is None
    assert "no exact baseline profile" in note


def test_legacy_baseline_never_silently_matches_subset():
    baseline, note = select_baseline_profile(ALL_METRICS, CARDIAC, ALL)

    assert baseline is None
    assert "legacy baseline has no gene scope" in note


def test_legacy_baseline_remains_usable_for_exact_full_manifest():
    baseline, note = select_baseline_profile(ALL_METRICS, ALL, ALL)

    assert baseline == ALL_METRICS
    assert "legacy all-gene" in note


def test_write_baseline_replaces_only_selected_scope():
    document = _v2()
    original_all = document["profiles"][baseline_scope_key(ALL)]
    replacement = {"recall": {"unique_variants": {"recall": 0.95}}}

    updated = update_baseline_document(
        document, CARDIAC, replacement, ALL, fixture_hash="new-fixture"
    )

    assert updated["profiles"][baseline_scope_key(ALL)] == original_all
    assert (
        updated["profiles"][baseline_scope_key(CARDIAC)]["recall"]
        == replacement["recall"]
    )
    assert (
        updated["profiles"][baseline_scope_key(CARDIAC)]["fixture_sha256"]
        == "new-fixture"
    )


def test_missing_requested_gene_refuses_partial_score(tmp_path: Path):
    db = tmp_path / "KCNH2.db"
    db.write_bytes(b"sqlite")

    with pytest.raises(SystemExit, match="missing requested gene DB.*KCNQ1"):
        require_complete_gene_set({"KCNH2": db}, ["KCNH2", "KCNQ1"])


def test_regression_gate_fails_closed_on_degraded_extract():
    summary = {
        "aggregate_recall": {
            "unique_variants": {"matched": 92, "gold": 100, "recall": 0.92}
        }
    }
    baseline = {"recall": {"unique_variants": {"recall": 0.92}}}

    problems = gate_result(
        summary,
        baseline,
        0.005,
        0.05,
        {"KCNQ1": "source recovery failed"},
    )

    assert any("extract.KCNQ1" in problem for problem in problems)


def test_degraded_extract_cannot_write_baseline():
    with pytest.raises(SystemExit, match="refusing to write baseline.*KCNQ1"):
        require_clean_baseline_write({"KCNQ1": "source recovery failed"})
