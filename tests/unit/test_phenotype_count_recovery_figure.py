"""Regression tests for the automatic phenotype-count recovery figure."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmarks.codex_paper_eval import run_eval
from scripts import build_phenotype_count_recovery as figure
from scripts import build_stratified_phenotype_count_recovery as stratified


def test_all_identity_misses_render_without_fabricating_count_agreement(tmp_path):
    missing = {
        field: {"count_extracted": 0, "zero_filled": 0, "identity_miss": 20}
        for field in figure.FIELDS
    }
    report = {
        "overall": {
            "count": {
                field: {"gold_asserted": 20, "predicted": 0} for field in figure.FIELDS
            }
        }
    }
    summary = figure.summarize([], missing, report)
    for record in summary.values():
        assert record["identity_coverage"] == 0
        assert record["identity_miss"] == 20
        assert record["mae"] is None
        assert record["exact_fraction"] is None
        assert record["nonzero_exact_fraction"] is None
    svg = figure.render_zero_filled_svg([], summary, 0, "empty", "2026-09-05")
    assert svg.count("No identity-matched count rows") == 2
    path = tmp_path / "rows.csv"
    figure.write_csv([], path)
    assert len(path.read_text().splitlines()) == 1
    assert "manual_count,ai_count_raw" in path.read_text()


def test_marker_size_key_uses_the_same_radius_function_as_plot_bubbles():
    radii = [figure.bubble_radius(count) for count in figure.MARKER_SIZE_KEY_COUNTS]

    assert radii == sorted(radii)
    assert len(set(radii)) == len(radii)
    key = "\n".join(figure.draw_marker_size_key(100.0, 190.0))
    assert "Marker size · variant–paper pairs" in key
    for count, radius in zip(figure.MARKER_SIZE_KEY_COUNTS, radii, strict=True):
        assert f'r="{radius:.1f}"' in key
        assert f">{count}</text>" in key


def test_unavailable_preferred_example_is_omitted_without_breaking_later_runs():
    assert figure.draw_example_annotations([], "affected", 100.0, 190.0, 900.0) == []


def test_empty_custom_annotation_set_disables_default_examples():
    row = {
        "measure": "affected",
        "gene": "RYR2",
        "pmid": "15466642",
        "manual_variant": "I4848V",
        "manual_count": 4,
        "ai_count_evaluation": 4,
    }

    assert (
        figure.draw_example_annotations(
            [row],
            "affected",
            100.0,
            190.0,
            900.0,
            annotations={},
            docks={},
        )
        == []
    )


def test_annotation_connectors_remain_diagonal_and_dock_on_title_block():
    kcnq1 = figure.annotation_connector_anchor(
        point_x=790.6,
        point_y=1044.0,
        label_x=590.6,
        label_y=674.0,
        gene="KCNQ1",
        variant="V254M",
        pmid="14678125",
    )
    assert min(abs(790.6 - kcnq1[0]), abs(1044.0 - kcnq1[1])) >= 24.0
    assert kcnq1[0] < 590.6

    ryr2 = figure.annotation_connector_anchor(
        point_x=804.2,
        point_y=285.2,
        label_x=609.2,
        label_y=360.2,
        gene="RYR2",
        variant="G357S",
        pmid="25814417",
    )
    assert min(abs(804.2 - ryr2[0]), abs(285.2 - ryr2[1])) >= 24.0
    assert ryr2[0] < 609.2

    i4848v = figure.annotation_connector_anchor(
        point_x=300.1,
        point_y=855.3,
        label_x=450.1,
        label_y=890.3,
        gene="RYR2",
        variant="I4848V",
        pmid="15466642",
    )
    assert i4848v == (438.1, 907.3)


def test_scoring_artifact_builder_writes_run_scoped_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    run_dir = tmp_path / "run"
    run_dir.mkdir()

    def fake_run(command, **_kwargs):
        for flag in ("--svg-out", "--png-out", "--pdf-out", "--csv-out", "--json-out"):
            path = Path(command[command.index(flag) + 1])
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("fixture")
        return SimpleNamespace(returncode=0, stdout="generated", stderr="")

    monkeypatch.setattr(run_eval.subprocess, "run", fake_run)
    report = {
        "overall": {
            "count": {
                "affected": {"gold_asserted": 1},
                "unaffected": {"gold_asserted": 1},
            }
        }
    }

    result = run_eval.build_phenotype_count_recovery_artifacts(
        run_dir, tmp_path / "gold", report
    )

    assert result["status"] == "generated"
    assert result["files"] == {
        "svg": "figures/phenotype_count_recovery.svg",
        "png": "figures/phenotype_count_recovery.png",
        "pdf": "figures/phenotype_count_recovery.pdf",
        "csv": "figures/data/phenotype_count_recovery.csv",
        "json": "figures/data/phenotype_count_recovery.json",
    }


def test_candidate_score_refreshes_both_cumulative_figures(
    monkeypatch: pytest.MonkeyPatch,
):
    commands: list[list[str]] = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout="generated", stderr="")

    monkeypatch.setattr(run_eval.subprocess, "run", fake_run)

    result = run_eval.build_cumulative_phenotype_count_artifacts(
        {"tranche_consumption": {"comparison_arm": "candidate"}}
    )

    assert result["status"] == "generated"
    assert [Path(command[1]).name for command in commands] == [
        "build_combined_phenotype_count_recovery.py",
        "build_stratified_phenotype_count_recovery.py",
    ]


def test_baseline_score_does_not_refresh_cumulative_figures():
    result = run_eval.build_cumulative_phenotype_count_artifacts(
        {"tranche_consumption": {"comparison_arm": "baseline"}}
    )

    assert result["status"] == "not_applicable"


def test_stratified_cohort_attempt_labels_are_derived_not_hard_coded():
    svg = stratified.render_svg(
        [],
        [],
        historical_attempt_count=118,
        mixed_attempt_count=218,
    )

    assert "Historical cardiac evaluation · 118 gene–paper attempts" in svg
    assert "Opened mixed-gold candidate arms · 218 gene–paper attempts" in svg
    assert "Opened mixed-gold calibration · 98 gene–paper attempts" not in svg
