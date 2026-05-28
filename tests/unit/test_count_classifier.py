"""Tests for pipeline.count_classifier (#5)."""

import pytest

from pipeline.count_classifier import (
    ACCEPTED_COUNT_TYPE,
    detect_misclassified_counts,
    enforce_per_variant_policy,
)


def _variant(carriers=None, affected=None, unaffected=None, provenance=None):
    """Build a variant dict matching the production schema, including the
    new count_provenance block when provided.
    """
    v = {}
    if carriers is not None:
        v["patients"] = {"count": carriers}
        v.setdefault("penetrance_data", {})["total_carriers_observed"] = carriers
    if affected is not None or unaffected is not None:
        v.setdefault("penetrance_data", {})
        if affected is not None:
            v["penetrance_data"]["affected_count"] = affected
        if unaffected is not None:
            v["penetrance_data"]["unaffected_count"] = unaffected
    if provenance is not None:
        v["count_provenance"] = provenance
    return v


# ---- detection -------------------------------------------------------------


def test_detects_cohort_total_carriers_assignment():
    """KCNQ1 29622001 G589D pattern: LLM emits count_type=cohort_total
    but still wrote 453 into patients.count — refuse it."""
    variant = _variant(
        carriers=453,
        provenance={
            "carriers_column_label": "Total N",
            "carriers_count_type": "cohort_total",
        },
    )
    anns = detect_misclassified_counts([variant])
    assert len(anns) == 1
    ann = anns[0]
    assert ann.variant_index == 0
    assert ann.field == "carriers"
    assert ann.value == 453
    assert ann.declared_count_type == "cohort_total"
    assert ann.column_label == "Total N"


def test_passes_through_per_variant_carrier():
    """Properly declared per-variant counts must not be flagged."""
    variant = _variant(
        carriers=5,
        provenance={
            "carriers_column_label": "No. of patients",
            "carriers_count_type": ACCEPTED_COUNT_TYPE,
        },
    )
    assert detect_misclassified_counts([variant]) == []


def test_silently_skips_variants_without_provenance():
    """Backward compat: old extraction JSONs have no count_provenance block.
    The classifier must not touch them — only the outlier guard handles that case."""
    variant = _variant(carriers=453)  # no provenance at all
    assert detect_misclassified_counts([variant]) == []


def test_skips_when_count_already_null():
    """No value to refuse → no annotation."""
    variant = _variant(
        provenance={
            "carriers_column_label": "Total N",
            "carriers_count_type": "cohort_total",
        }
    )
    assert detect_misclassified_counts([variant]) == []


def test_each_field_classified_independently():
    """Carriers may be per_variant while affected is cohort_total, etc."""
    variant = _variant(
        carriers=5,
        affected=120,
        unaffected=80,
        provenance={
            "carriers_column_label": "No. of patients",
            "carriers_count_type": "per_variant_carrier",
            "affected_column_label": "All affected",
            "affected_count_type": "cohort_total",
            "unaffected_column_label": "Controls",
            "unaffected_count_type": "control",
        },
    )
    anns = detect_misclassified_counts([variant])
    fields_flagged = {a.field for a in anns}
    assert fields_flagged == {"affected", "unaffected"}


def test_unknown_count_type_is_refused():
    """The closed vocabulary includes 'unknown' as a non-per-variant type."""
    variant = _variant(
        carriers=100,
        provenance={"carriers_count_type": "unknown"},
    )
    anns = detect_misclassified_counts([variant])
    assert len(anns) == 1
    assert anns[0].declared_count_type == "unknown"


def test_count_type_canonicalization_handles_capitalization():
    """LLMs sometimes return 'Cohort_Total' or 'COHORT_TOTAL'; still refused."""
    variant = _variant(
        carriers=200,
        provenance={"carriers_count_type": "Cohort_Total"},
    )
    anns = detect_misclassified_counts([variant])
    assert len(anns) == 1
    assert anns[0].declared_count_type == "cohort_total"


# ---- policy enforcement ----------------------------------------------------


def test_flag_policy_preserves_value_and_records_flag():
    variant = _variant(
        carriers=453,
        provenance={"carriers_count_type": "cohort_total"},
    )
    anns = detect_misclassified_counts([variant])
    result = enforce_per_variant_policy([variant], anns, policy="flag")
    assert result.flagged == 1
    assert result.cleared == 0
    # Both mirrors preserved
    assert variant["patients"]["count"] == 453
    assert variant["penetrance_data"]["total_carriers_observed"] == 453
    # Flag metadata present
    flags = variant["count_classifier_flags"]["carriers"]
    assert flags["raw"] == 453
    assert flags["declared_count_type"] == "cohort_total"
    assert flags["policy"] == "flag"


def test_clear_policy_nulls_both_mirrors():
    """Same mirror invariant as the outlier guard."""
    variant = _variant(
        carriers=453,
        provenance={"carriers_count_type": "cohort_total"},
    )
    anns = detect_misclassified_counts([variant])
    result = enforce_per_variant_policy([variant], anns, policy="clear")
    assert result.cleared == 1
    assert variant["patients"]["count"] is None
    assert variant["penetrance_data"]["total_carriers_observed"] is None
    assert variant["count_classifier_flags"]["carriers"]["raw"] == 453


def test_off_policy_is_noop():
    variant = _variant(
        carriers=453,
        provenance={"carriers_count_type": "cohort_total"},
    )
    anns = detect_misclassified_counts([variant])
    result = enforce_per_variant_policy([variant], anns, policy="off")
    assert result.flagged == 0
    assert variant["patients"]["count"] == 453
    assert "count_classifier_flags" not in variant


def test_unknown_policy_raises():
    with pytest.raises(ValueError):
        enforce_per_variant_policy([], [], policy="invalid")


def test_composes_with_outlier_guard():
    """A value that is BOTH provenance-misclassified AND a statistical outlier
    should be safely cleared by both guards in sequence (idempotent on a
    null count)."""
    from pipeline.count_outlier_guard import (
        apply_outlier_policy,
        detect_count_outliers,
    )

    variants = [
        _variant(
            carriers=453,
            provenance={"carriers_count_type": "cohort_total"},
        ),
        _variant(carriers=4),
        _variant(carriers=5),
        _variant(carriers=6),
        _variant(carriers=3),
    ]
    # Classifier runs first
    ca = detect_misclassified_counts(variants)
    enforce_per_variant_policy(variants, ca, policy="clear")
    # Outlier guard now runs on the post-classifier state; the 453 row is
    # now null and should not be re-flagged.
    oa = detect_count_outliers(variants)
    assert not any(a.variant_index == 0 for a in oa)
    apply_outlier_policy(variants, oa, policy="clear")
    # Final state: both mirrors null, classifier flag present.
    assert variants[0]["patients"]["count"] is None
    assert variants[0]["penetrance_data"]["total_carriers_observed"] is None
    assert variants[0]["count_classifier_flags"]["carriers"]["raw"] == 453
