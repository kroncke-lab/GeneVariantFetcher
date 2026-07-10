"""End-to-end count error: misses count as zero, tail is reported.

Guards the honest count metric that complements matched-only MAE. A variant that
is missed entirely, or matched but never assigned a count, must contribute its
full gold count to the error -- so the metric cannot be flattered by skipping
hard cases.
"""

import pytest

from cli.compare_variants import (
    ComparisonRow,
    compute_end_to_end_count_error,
    compute_rows_mae,
)


def _row(**kw) -> ComparisonRow:
    base = dict(
        pmid="1",
        excel_variant_raw="p.V1M",
        excel_variant_norm="p.V1M",
        sqlite_variant_raw=None,
        sqlite_variant_norm=None,
        match_type="none",
        match_score=None,
        excel_carriers_total=None,
        sqlite_carriers_total=None,
        carriers_diff=None,
        excel_affected=None,
        sqlite_affected=None,
        affected_diff=None,
        excel_unaffected=None,
        sqlite_unaffected=None,
        unaffected_diff=None,
    )
    base.update(kw)
    return ComparisonRow(**base)


def test_misses_count_as_zero_and_tail_is_reported():
    results = [
        # matched, both counts present -> error 2
        _row(match_type="exact", excel_carriers_total=10, sqlite_carriers_total=8),
        # missed entirely -> extracted treated as 0 -> error 50
        _row(match_type="none", missing_in_sqlite=True, excel_carriers_total=50),
        # matched but never assigned a count -> extracted 0 -> error 5
        _row(match_type="exact", excel_carriers_total=5, sqlite_carriers_total=None),
        # DB-only extra: not a gold row -> excluded
        _row(match_type="none", missing_in_excel=True, sqlite_carriers_total=100),
        # gold row with no gold carrier count -> skipped for carriers
        _row(match_type="exact", excel_affected=3, sqlite_affected=3),
    ]

    e2e = compute_end_to_end_count_error(results)["carriers"]
    assert e2e["n"] == 3
    assert e2e["n_missed"] == 2  # the miss + the matched-but-uncounted row
    assert e2e["sum_abs_error"] == 57
    assert e2e["mae"] == 19.0
    assert e2e["median"] == 5.0
    assert e2e["p95"] == pytest.approx(45.5)
    assert e2e["max"] == 50

    # Matched-only MAE sees only the one both-present row (error 2), so it reads
    # far lower than the end-to-end error -- the exact gap this metric exposes.
    matched_only = compute_rows_mae(results)["carriers"]
    assert matched_only["mae"] == 2.0
    assert matched_only["n_matched"] == 1


def test_empty_fields_are_none():
    e2e = compute_end_to_end_count_error([])["unaffected"]
    assert e2e == {
        "n": 0,
        "n_missed": 0,
        "sum_abs_error": 0,
        "mae": None,
        "median": None,
        "p95": None,
        "max": None,
    }
