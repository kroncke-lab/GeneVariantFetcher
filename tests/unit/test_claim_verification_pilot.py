"""Tests for recall-audit claim-verification pilot helpers."""

from scripts.recall_audit.run_claim_verification_pilot import best_gold_row


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
