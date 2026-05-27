"""Unit tests for pipeline.count_outlier_guard."""

import pytest

from pipeline.count_outlier_guard import (
    apply_outlier_policy,
    detect_count_outliers,
)


def _v(carriers=None, affected=None, unaffected=None):
    variant = {}
    if carriers is not None:
        variant["patients"] = {"count": carriers}
    if affected is not None or unaffected is not None:
        variant["penetrance_data"] = {}
        if affected is not None:
            variant["penetrance_data"]["affected_count"] = affected
        if unaffected is not None:
            variant["penetrance_data"]["unaffected_count"] = unaffected
    return variant


def test_detects_kcnq1_g589d_style_outlier():
    """KCNQ1 29622001 G589D pattern: one row at 453, others at ~5."""
    variants = [
        _v(carriers=453),
        _v(carriers=4),
        _v(carriers=5),
        _v(carriers=6),
        _v(carriers=3),
        _v(carriers=8),
    ]
    annotations = detect_count_outliers(variants)
    assert len(annotations) == 1
    ann = annotations[0]
    assert ann.variant_index == 0
    assert ann.field == "carriers"
    assert ann.value == 453


def test_no_outliers_when_paper_too_small():
    """A 3-variant paper is too small to compute a robust median; skip."""
    variants = [_v(carriers=453), _v(carriers=4), _v(carriers=5)]
    assert detect_count_outliers(variants) == []


def test_no_outliers_when_value_below_absolute_threshold():
    """A value above the multiplier ratio but <= 50 absolute is not flagged."""
    variants = [
        _v(carriers=49),  # 49 > 10*4=40 ratio, but absolute_threshold=50 blocks it
        _v(carriers=4),
        _v(carriers=4),
        _v(carriers=5),
        _v(carriers=3),
    ]
    assert detect_count_outliers(variants) == []


def test_no_outliers_when_legitimate_consortium_table():
    """When several rows have large counts, the median is high and nothing fires."""
    variants = [
        _v(carriers=85),
        _v(carriers=120),
        _v(carriers=137),
        _v(carriers=92),
        _v(carriers=78),
    ]
    assert detect_count_outliers(variants) == []


def test_independent_fields_flagged_independently():
    """Carriers outlier on row 0; affected outlier on row 1."""
    variants = [
        _v(carriers=200, affected=2),
        _v(carriers=3, affected=180),
        _v(carriers=4, affected=2),
        _v(carriers=5, affected=3),
        _v(carriers=6, affected=2),
    ]
    annotations = detect_count_outliers(variants)
    fields_by_idx = {(a.variant_index, a.field) for a in annotations}
    assert (0, "carriers") in fields_by_idx
    assert (1, "affected") in fields_by_idx


def test_apply_policy_flag_preserves_value_and_adds_metadata():
    variants = [
        _v(carriers=453),
        _v(carriers=4),
        _v(carriers=5),
        _v(carriers=6),
        _v(carriers=3),
    ]
    anns = detect_count_outliers(variants)
    result = apply_outlier_policy(variants, anns, policy="flag")
    assert result.flagged == 1
    assert result.cleared == 0
    # Raw count preserved
    assert variants[0]["patients"]["count"] == 453
    # Flag metadata recorded
    flags = variants[0]["count_outlier_flags"]
    assert flags["carriers"]["raw"] == 453
    assert flags["carriers"]["policy"] == "flag"


def test_apply_policy_clear_zeros_value_and_records_raw():
    variants = [
        _v(carriers=453),
        _v(carriers=4),
        _v(carriers=5),
        _v(carriers=6),
        _v(carriers=3),
    ]
    anns = detect_count_outliers(variants)
    result = apply_outlier_policy(variants, anns, policy="clear")
    assert result.flagged == 1
    assert result.cleared == 1
    # Count is now None
    assert variants[0]["patients"]["count"] is None
    # Raw is preserved under flags
    assert variants[0]["count_outlier_flags"]["carriers"]["raw"] == 453


def test_apply_policy_off_is_noop():
    variants = [_v(carriers=453)] + [_v(carriers=4) for _ in range(5)]
    anns = detect_count_outliers(variants)
    result = apply_outlier_policy(variants, anns, policy="off")
    assert result.flagged == 0
    assert variants[0]["patients"]["count"] == 453
    assert "count_outlier_flags" not in variants[0]


def test_apply_policy_unknown_raises():
    with pytest.raises(ValueError):
        apply_outlier_policy([], [], policy="invalid")


def test_detector_tolerates_missing_count_fields():
    """Variants without a patients/penetrance_data block are skipped, not errored."""
    variants = [
        _v(carriers=453),
        {},  # no fields
        _v(),  # patients absent
        _v(carriers=4),
        _v(carriers=5),
        _v(carriers=6),
    ]
    anns = detect_count_outliers(variants)
    # Still triggers on the 453 row; the empty variants are silently skipped
    assert any(a.field == "carriers" and a.variant_index == 0 for a in anns)


def test_zero_and_none_counts_excluded_from_median():
    """Zeros and Nones should not depress the median artificially."""
    variants = [
        _v(carriers=137),
        _v(carriers=0),
        _v(carriers=None),
        _v(carriers=2),
        _v(carriers=3),
        _v(carriers=2),
    ]
    # Non-zero/non-null values: [137, 2, 3, 2] → median 2.5; 137 > 25 (10x) and > 50.
    anns = detect_count_outliers(variants)
    assert any(a.value == 137 for a in anns)
