"""Tests for recall-audit claim-verification pilot helpers."""

from scripts.recall_audit.run_claim_verification_pilot import (
    best_gold_row,
    gold_value,
    gold_value_source,
)


def test_best_gold_row_preserves_duplicate_variant_cohorts():
    gold_rows = [
        {"variant": "R176W", "carriers": "112", "affected": "112", "unaffected": "0"},
        {"variant": "R176W", "carriers": "16", "affected": "0", "unaffected": "16"},
    ]

    selected = best_gold_row(
        gold_rows,
        {"total_carriers": 112, "affected": 18, "unaffected": 94},
    )

    assert selected["carriers"] == "112"


def test_gold_v2_values_preserve_original_columns():
    row = {
        "variant": "S1103Y",
        "carriers": "85",
        "affected": "39",
        "unaffected": "46",
        "gold_v2_carriers": "26",
        "gold_v2_affected": "17",
        "gold_v2_unaffected": "9",
        "gold_v2_status": "adjudicated_variant_carrier_count",
    }

    assert gold_value(row, "total_carriers", gold_value_set="v2") == 26
    assert gold_value(row, "affected", gold_value_set="v2") == 17
    assert gold_value(row, "unaffected", gold_value_set="v2") == 9
    assert gold_value(row, "total_carriers", gold_value_set="original") == 85
    assert gold_value_source(row, gold_value_set="v2") == "gold_v2"


def test_gold_v2_blank_count_is_explicit_null_when_status_populated():
    row = {
        "variant": "c.3599-9delT",
        "carriers": "2",
        "affected": "1",
        "unaffected": "1",
        "gold_v2_carriers": "1",
        "gold_v2_affected": "1",
        "gold_v2_unaffected": "",
        "gold_v2_status": "adjudicated_null_unaffected",
    }

    assert gold_value(row, "unaffected", gold_value_set="v2") is None
    assert gold_value(row, "unaffected", gold_value_set="original") == 1


def test_best_gold_row_scores_against_adjudicated_values_when_present():
    gold_rows = [
        {
            "variant": "S1103Y",
            "carriers": "85",
            "affected": "39",
            "unaffected": "46",
            "gold_v2_carriers": "26",
            "gold_v2_affected": "17",
            "gold_v2_unaffected": "9",
            "gold_v2_status": "adjudicated_variant_carrier_count",
        },
        {"variant": "S1103Y", "carriers": "1", "affected": "1", "unaffected": "0"},
    ]

    selected = best_gold_row(
        gold_rows,
        {"total_carriers": 26, "affected": 17, "unaffected": 9},
        gold_value_set="v2",
    )

    assert selected["gold_v2_carriers"] == "26"
