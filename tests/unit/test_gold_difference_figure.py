"""Regression tests for the gold difference figure (automated minus reference)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmarks.codex_paper_eval import run_eval
from scripts import build_gold_difference_figure as figure

GOLD_CSV = (
    "variant,pmid,carriers,affected,unaffected,gold_v2_carriers,gold_v2_affected,"
    "gold_v2_unaffected,gold_v2_status,gold_v2_note,gold_v2_source\n"
    "A1V,1,4,2,2,,,,,,\n"
    "R2W,1,3,3,0,,,,,,\n"
    "G3D,1,5,,,,,,,,\n"
    "DUP,1,9,9,0,,,,excluded_duplicate_current_cohort,,\n"
    "V4A,2,10,6,4,,,,,,\n"
)
LOCKED_AT = "2026-09-04T00:00:00+00:00"

# Candidate: A1V exact carriers/affected but abstains on unaffected; R2W is never
# found; G3D supplies a wrong carrier count; X9Y is an extra; V4A is matched but
# supplies nothing.
CANDIDATE = {
    "1": [
        {"variant": "A1V", "carriers": 4, "affected": 2, "unaffected": None},
        {"variant": "G3D", "carriers": 7, "affected": None, "unaffected": None},
        {"variant": "X9Y", "carriers": 1, "affected": None, "unaffected": None},
    ],
    "2": [{"variant": "V4A", "carriers": None, "affected": None, "unaffected": None}],
}
# Baseline: A1V exact; R2W matched without a count; G3D never found; V4A exact.
BASELINE = {
    "1": [
        {"variant": "A1V", "carriers": 4, "affected": 2, "unaffected": 2},
        {"variant": "R2W", "carriers": None, "affected": None, "unaffected": None},
    ],
    "2": [{"variant": "V4A", "carriers": 10, "affected": 6, "unaffected": 4}],
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _gold_root(tmp_path: Path) -> Path:
    gold = tmp_path / "gold"
    gold.mkdir(exist_ok=True)
    (gold / "SCN5A_recall_input.csv").write_text(GOLD_CSV)
    return gold


def _write_locked_run(
    run_dir: Path,
    gold_root: Path,
    run_id: str,
    variants_by_paper: dict[str, list[dict]],
    *,
    consumption: dict | None = None,
) -> Path:
    """Score a synthetic run with the real scorer, then lock it."""
    run_dir.mkdir(parents=True)
    selection = {
        "run_id": run_id,
        "seed": 1,
        "papers": [
            {"gene": "SCN5A", "pmid": "1", "gold_provenance": "fixture"},
            {"gene": "SCN5A", "pmid": "2", "gold_provenance": "fixture"},
        ],
    }
    predictions = {
        "run_id": run_id,
        "papers": [
            {"gene": "SCN5A", "pmid": pmid, "tool": "text", "variants": rows}
            for pmid, rows in variants_by_paper.items()
        ],
    }
    (run_dir / "selection.json").write_text(json.dumps(selection))
    (run_dir / "predictions.json").write_text(json.dumps(predictions))
    scores, _by_gene, _by_provenance = run_eval.score_prediction_lane(
        selection, predictions, gold_root, "primary"
    )
    report = {
        "run_id": run_id,
        "locked_at": LOCKED_AT,
        "overall": run_eval.aggregate(scores),
        "papers": scores,
        "gold_sources": {"SCN5A": str(run_eval.gold_csv_path(gold_root, "SCN5A"))},
    }
    if consumption:
        report["tranche_consumption"] = consumption
    (run_dir / "report.json").write_text(json.dumps(report))
    (run_dir / "LOCK.json").write_text(
        json.dumps(
            {
                "locked_at": LOCKED_AT,
                "selection_sha256": _sha(run_dir / "selection.json"),
                "predictions_sha256": _sha(run_dir / "predictions.json"),
            }
        )
    )
    return run_dir


def _by_variant(rows: list[dict], measure: str) -> dict[str, dict]:
    return {row["gold_variant"]: row for row in rows if row["measure"] == measure}


def test_rows_classify_every_asserted_gold_row_and_reproduce_the_report(tmp_path):
    gold = _gold_root(tmp_path)
    run_dir = _write_locked_run(tmp_path / "candidate", gold, "cand", CANDIDATE)

    arm = figure.load_arm("run", run_dir, None)

    carriers = _by_variant(arm["rows"], "carriers")
    assert set(carriers) == {"A1V", "R2W", "G3D", "V4A"}  # DUP is excluded
    assert (carriers["A1V"]["status"], carriers["A1V"]["difference"]) == ("supplied", 0)
    assert (carriers["R2W"]["status"], carriers["R2W"]["difference"]) == (
        "identity_miss",
        -3,
    )
    assert (carriers["G3D"]["status"], carriers["G3D"]["difference"]) == ("supplied", 2)
    assert (carriers["V4A"]["status"], carriers["V4A"]["difference"]) == (
        "abstained",
        -10,
    )
    assert carriers["V4A"]["automated_count_raw"] is None
    assert carriers["V4A"]["automated_count_evaluated"] == 0
    # G3D asserts only carriers, so it never enters the affected/unaffected panels.
    assert set(_by_variant(arm["rows"], "affected")) == {"A1V", "R2W", "V4A"}

    metrics = arm["metrics"]["carriers"]
    assert metrics["gold_rows"] == 4
    assert metrics["identity_matched"] == 3
    assert metrics["identity_miss"] == 1
    assert (metrics["supplied"], metrics["abstained"]) == (2, 1)
    assert metrics["exact"] == 1
    assert metrics["end_to_end_mae"] == pytest.approx((0 + 3 + 2 + 10) / 4)
    assert metrics["conditional_mae"] == pytest.approx(1.0)
    assert metrics["bias"] == pytest.approx((0 - 3 + 2 - 10) / 4)
    assert metrics["coverage_on_matched"] == pytest.approx(2 / 3)
    # A gold zero that nothing supplied is exact for free and is reported as such.
    unaffected = arm["metrics"]["unaffected"]
    assert unaffected["zero_baseline_exact"] == 1
    assert unaffected["exact"] == 1
    assert arm["identity"] == {
        "tp": 3,
        "fp": 1,
        "fn": 1,
        "recall": pytest.approx(0.75),
        "precision": pytest.approx(0.75),
        "primary_score_lane": "primary",
    }


def test_paired_arms_report_row_level_changes_and_status_transitions(tmp_path):
    gold = _gold_root(tmp_path)
    candidate_dir = _write_locked_run(tmp_path / "candidate", gold, "cand", CANDIDATE)
    baseline_dir = _write_locked_run(tmp_path / "baseline", gold, "base", BASELINE)

    result = figure.build_figure(
        candidate_dir, baseline_run_dir=baseline_dir, auto_baseline=False
    )

    paired = result["summary"]["paired"]
    carriers = paired["by_measure"]["carriers"]
    # G3D: miss (-5) -> supplied (+2) improved; V4A: exact -> abstained worsened;
    # A1V unchanged; R2W abstained (-3) -> miss (-3) same size, new status.
    assert (
        carriers["rows_improved"],
        carriers["rows_worsened"],
        carriers["rows_unchanged"],
    ) == (1, 1, 2)
    assert carriers["status_transitions"]["identity_miss->supplied"] == 1
    assert carriers["status_transitions"]["abstained->identity_miss"] == 1
    assert carriers["delta_candidate_minus_baseline"][
        "end_to_end_mae"
    ] == pytest.approx((0 + 3 + 2 + 10) / 4 - (0 + 3 + 5 + 0) / 4)
    assert paired["identity"]["recall"][
        "delta_candidate_minus_baseline"
    ] == pytest.approx(0.0)
    assert set(result["summary"]["arms"]) == {"baseline", "candidate"}
    assert {row["arm"] for row in result["rows"]} == {"baseline", "candidate"}
    assert len(result["rows"]) == 2 * (4 + 3 + 3)

    svg = result["svg"]
    for text in (
        "Baseline arm",
        "Candidate arm",
        "Candidate − baseline on the same",
        "Carriers",
        "Affected individuals",
        "Unaffected individuals",
        "no registered comparison found",
        'class="bias-line"',
        'class="floor-line"',
    ):
        assert text in svg
    # The axes come from the data: the largest reference count is 10.
    assert result["contract"]["x_max"] == 20
    assert result["contract"]["y_max"] == 20


def test_single_run_renders_without_a_delta_strip(tmp_path):
    gold = _gold_root(tmp_path)
    run_dir = _write_locked_run(tmp_path / "solo", gold, "solo", CANDIDATE)

    result = figure.build_figure(run_dir, auto_baseline=False)

    assert list(result["summary"]["arms"]) == ["run"]
    assert result["summary"]["paired"] is None
    assert "Candidate − baseline" not in result["svg"]
    assert "Run · solo" in result["svg"]


def test_baseline_is_discovered_from_the_tranche_consumption_log(tmp_path):
    gold = _gold_root(tmp_path)
    registry_dir = tmp_path / "suite"
    registry_dir.mkdir()
    (registry_dir / "registry.json").write_text("{}")
    baseline_dir = _write_locked_run(tmp_path / "runs" / "base", gold, "base", BASELINE)
    consumption = {"tier_id": "t01", "comparison_arm": "candidate", "run_id": "cand"}
    candidate_dir = _write_locked_run(
        tmp_path / "runs" / "cand", gold, "cand", CANDIDATE, consumption=consumption
    )
    (candidate_dir / "setup.json").write_text(
        json.dumps(
            {
                "cohort": {
                    "registry": str(registry_dir / "registry.json"),
                    "consumption_log": {"path": "consumption_log.jsonl"},
                }
            }
        )
    )
    (registry_dir / "consumption_log.jsonl").write_text(
        json.dumps(
            {
                "tier_id": "t01",
                "comparison_arm": "baseline",
                "run_id": "base",
                "scored_at": "2026-09-03T00:00:00+00:00",
            }
        )
        + "\n"
        + json.dumps({"tier_id": "t02", "comparison_arm": "baseline", "run_id": "x"})
        + "\n"
    )
    report = json.loads((candidate_dir / "report.json").read_text())

    assert figure.discover_baseline(candidate_dir, report) == baseline_dir
    baseline_report = json.loads((baseline_dir / "report.json").read_text())
    assert figure.discover_baseline(baseline_dir, baseline_report) is None


def test_moved_predictions_are_refused(tmp_path):
    gold = _gold_root(tmp_path)
    run_dir = _write_locked_run(tmp_path / "moved", gold, "moved", CANDIDATE)
    (run_dir / "predictions.json").write_text("{}")

    with pytest.raises(SystemExit, match="lock mismatch"):
        figure.load_arm("run", run_dir, None)


def test_square_page_pads_a_portrait_page_to_a_square():
    svg = (
        '<?xml version="1.0"?>\n<svg xmlns="http://www.w3.org/2000/svg" '
        'width="8.00in" height="8.63in" viewBox="0 0 2400 2590"><g/></svg>'
    )
    out = figure.square_page(svg, 2400, 2590, 2590)
    assert 'viewBox="0 0 2590 2590"' in out
    assert 'width="8.63in"' in out and 'height="8.63in"' in out
    assert out.startswith('<?xml version="1.0"?>')


def test_score_artifact_builder_records_outputs_and_never_blocks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    commands: list[list[str]] = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        out_dir = Path(command[command.index("--out-dir") + 1])
        (out_dir / "data").mkdir(parents=True, exist_ok=True)
        (out_dir / "gold_difference.svg").write_text("<svg/>")
        (out_dir / "data" / "gold_difference.json").write_text("{}")
        return SimpleNamespace(returncode=0, stdout="wrote", stderr="")

    monkeypatch.setattr(run_eval.subprocess, "run", fake_run)
    report = {"overall": {"count": {"carriers": {"gold_asserted": 3}}}}

    result = run_eval.build_gold_difference_artifacts(
        run_dir, tmp_path / "gold", report
    )

    assert result["status"] == "generated"
    assert result["files"] == {
        "svg": "figures/gold_difference.svg",
        "json": "figures/data/gold_difference.json",
    }
    assert Path(commands[0][1]).name == "build_gold_difference_figure.py"
    assert commands[0][commands[0].index("--gold-root") + 1] == str(tmp_path / "gold")

    def failing_run(command, **_kwargs):
        return SimpleNamespace(returncode=2, stdout="", stderr="renderer exploded")

    monkeypatch.setattr(run_eval.subprocess, "run", failing_run)
    failed = run_eval.build_gold_difference_artifacts(
        run_dir, tmp_path / "gold", report
    )
    assert failed["status"] == "failed"
    assert "renderer exploded" in failed["detail"]

    empty = run_eval.build_gold_difference_artifacts(
        run_dir, tmp_path / "gold", {"overall": {"count": {"carriers": {}}}}
    )
    assert empty["status"] == "not_applicable"


def test_compare_refresh_passes_the_baseline_and_the_comparison(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    baseline_dir = tmp_path / "base"
    candidate_dir = tmp_path / "cand"
    for directory in (baseline_dir, candidate_dir):
        directory.mkdir()
        (directory / "predictions.json").write_text("{}")
        (directory / "LOCK.json").write_text("{}")
        (directory / "report.json").write_text("{}")
    comparison = tmp_path / "cand" / "compare.json"
    comparison.write_text("{}")
    commands: list[list[str]] = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(run_eval.subprocess, "run", fake_run)

    result = run_eval.refresh_gold_difference_after_compare(
        baseline_dir / "report.json",
        candidate_dir / "report.json",
        comparison,
        tmp_path / "gold",
    )

    assert result["status"] == "generated"
    command = commands[0]
    assert command[command.index("--run-dir") + 1] == str(candidate_dir)
    assert command[command.index("--baseline-run-dir") + 1] == str(baseline_dir)
    assert command[command.index("--comparison-json") + 1] == str(comparison.resolve())
    assert command[command.index("--out-dir") + 1] == str(candidate_dir / "figures")

    (baseline_dir / "LOCK.json").unlink()
    skipped = run_eval.refresh_gold_difference_after_compare(
        baseline_dir / "report.json", candidate_dir / "report.json", comparison, None
    )
    assert skipped["status"] == "not_applicable"
