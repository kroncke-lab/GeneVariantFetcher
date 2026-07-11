"""Study-context rules for the count classifier (Stage 5 A3)."""

from pipeline.count_classifier import detect_misclassified_counts


def _variant(carriers: int, declared: str | None = "per_variant_carrier"):
    return {
        "patients": {"count": carriers},
        "penetrance_data": {"total_carriers_observed": carriers},
        "count_provenance": {
            "carriers_column_label": "N",
            "carriers_count_type": declared,
        },
    }


def test_biobank_large_carrier_flagged_by_study_context():
    anns = detect_misclassified_counts(
        [_variant(200)],
        study_design="cohort_biobank",
    )
    assert anns
    assert "population" in anns[0].reason or "cohort" in anns[0].reason


def test_review_meta_carrier_is_trap():
    anns = detect_misclassified_counts(
        [_variant(5)],
        study_design="review_meta",
    )
    assert anns
    assert "trap" in anns[0].reason or "review" in anns[0].reason


def test_case_series_small_count_clean():
    anns = detect_misclassified_counts(
        [_variant(3)],
        study_design="case_series",
    )
    assert anns == []
