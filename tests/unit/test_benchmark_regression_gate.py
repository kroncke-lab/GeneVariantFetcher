"""run_benchmark.py --fail-on-regression must catch recall drops / MAE rises."""

from benchmarks.curated_extraction_eval.run_benchmark import (
    check_regression,
    gate_result,
)


def _summary(uniq_recall, carriers_mae):
    return {
        "aggregate_recall": {
            "unique_variants": {"matched": 90, "gold": 100, "recall": uniq_recall}
        },
        "aggregate_mae": {"carriers": {"mae": carriers_mae}},
    }


BASELINE = {
    "recall": {"unique_variants": {"recall": 0.90}},
    "mae": {"carriers": 0.30},
}


def test_no_regression_within_tolerance():
    assert check_regression(_summary(0.899, 0.34), BASELINE, 0.005, 0.05) == []


def test_recall_regression_flagged():
    problems = check_regression(_summary(0.88, 0.30), BASELINE, 0.005, 0.05)
    assert any("recall.unique_variants" in p for p in problems)


def test_mae_regression_flagged():
    problems = check_regression(_summary(0.90, 0.40), BASELINE, 0.005, 0.05)
    assert any("mae.carriers" in p for p in problems)


def test_no_baseline_is_noop():
    assert check_regression(_summary(0.10, 9.9), None, 0.005, 0.05) == []


def test_missing_baseline_fails_closed():
    """--fail-on-regression with no baseline must FAIL, not silently pass."""
    problems = gate_result(_summary(0.10, 9.9), None, 0.005, 0.05)
    assert problems, "a gate with no baseline must not report success"


def test_missing_current_metric_is_regression():
    """A metric present in the baseline but absent from the run is fail-closed."""
    summary = {
        "aggregate_recall": {"unique_variants": {"recall": None}},
        "aggregate_mae": {"carriers": {"mae": 0.30}},
    }
    problems = check_regression(summary, BASELINE, 0.005, 0.05)
    assert any("recall.unique_variants" in p and "missing" in p for p in problems)


def test_e2e_count_error_regression_flagged():
    """End-to-end count error rising above baseline trips the gate."""
    summary = {
        "aggregate_recall": {"unique_variants": {"recall": 0.90}},
        "aggregate_mae": {"carriers": {"mae": 0.30}},
        "aggregate_count_error_end_to_end": {"carriers": {"mae": 5.0}},
    }
    baseline = {**BASELINE, "count_error_end_to_end": {"carriers": 1.0}}
    problems = check_regression(summary, baseline, 0.005, 0.05)
    assert any("e2e.carriers" in p for p in problems)
