"""Tests for pipeline.evidence_card (#6)."""

import pytest

from pipeline.evidence_card import (
    GOLD_ERROR_THRESHOLD_DEFAULT,
    LARGE_VALUE_THRESHOLD_DEFAULT,
    VERDICT_REJECT,
    VERDICT_WITHHOLD,
    build_evidence_cards,
)


def _variant(
    *,
    protein="p.Gly1Asp",
    carriers=None,
    affected=None,
    provenance=None,
    outlier_flags=None,
    classifier_flags=None,
):
    v = {"protein_notation": protein}
    if carriers is not None:
        v["patients"] = {"count": carriers}
        v.setdefault("penetrance_data", {})["total_carriers_observed"] = carriers
    if affected is not None:
        v.setdefault("penetrance_data", {})["affected_count"] = affected
    if provenance is not None:
        v["count_provenance"] = provenance
    if outlier_flags is not None:
        v["count_outlier_flags"] = outlier_flags
    if classifier_flags is not None:
        v["count_classifier_flags"] = classifier_flags
    return v


# ---- triggers --------------------------------------------------------------


def test_gold_error_trigger_fires_when_diff_at_threshold():
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[_variant(carriers=453)],
        gold_counts_by_variant={"p.Gly1Asp": {"carriers": 7}},
    )
    assert len(cards) == 1
    c = cards[0]
    assert c.error == 446
    assert any(t.startswith("gold_error_ge_") for t in c.triggers)


def test_gold_error_trigger_does_not_fire_below_threshold():
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[_variant(carriers=8)],
        gold_counts_by_variant={"p.Gly1Asp": {"carriers": 7}},
    )
    # |8-7|=1 < 10 default; nothing else triggers; no card
    assert cards == []


def test_outlier_guard_flag_produces_card_in_no_gold_mode():
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[_variant(carriers=453, outlier_flags={"carriers": {"raw": 453}})],
        # no gold_counts_by_variant -> no-gold mode
    )
    assert len(cards) == 1
    assert "outlier_guard" in cards[0].triggers


def test_classifier_flag_produces_card_in_no_gold_mode():
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[
            _variant(
                carriers=453,
                classifier_flags={"carriers": {"raw": 453}},
            )
        ],
    )
    assert len(cards) == 1
    assert "classifier" in cards[0].triggers


def test_unknown_provenance_large_value_triggers():
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[
            _variant(
                carriers=100,
                provenance={"carriers_count_type": "unknown"},
            )
        ],
    )
    assert len(cards) == 1
    assert "unknown_provenance_large" in cards[0].triggers


def test_unknown_provenance_below_threshold_does_not_trigger():
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[
            _variant(
                carriers=10,  # < 50
                provenance={"carriers_count_type": "unknown"},
            )
        ],
    )
    assert cards == []


def test_per_variant_provenance_blocks_unknown_trigger():
    """Properly declared per-variant counts must not pull cards."""
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[
            _variant(
                carriers=200,
                provenance={"carriers_count_type": "per_variant_carrier"},
            )
        ],
    )
    assert cards == []


def test_multiple_triggers_combine_on_one_card():
    """When a row hits both outlier guard AND gold error, ONE card lists both."""
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[
            _variant(
                carriers=453,
                outlier_flags={"carriers": {"raw": 453}},
            )
        ],
        gold_counts_by_variant={"p.Gly1Asp": {"carriers": 7}},
    )
    assert len(cards) == 1
    assert {"outlier_guard", f"gold_error_ge_{GOLD_ERROR_THRESHOLD_DEFAULT}"} <= set(
        cards[0].triggers
    )


# ---- excerpts --------------------------------------------------------------


def test_source_excerpt_finds_value_in_text():
    source = "\n".join(
        [
            "Methods section.",
            "We screened 200 patients and identified the variant.",
            "Table 1 shows the breakdown.",
            "Total carriers: 453 across the cohort.",
            "Discussion follows.",
        ]
    )
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[_variant(carriers=453, outlier_flags={"carriers": {"raw": 453}})],
        source_text=source,
    )
    assert "453" in cards[0].source_excerpt
    assert "Total carriers" in cards[0].source_excerpt


# ---- heuristic verdict -----------------------------------------------------


def test_heuristic_rejects_when_count_type_is_cohort_total():
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[
            _variant(
                carriers=453,
                provenance={"carriers_count_type": "cohort_total"},
                classifier_flags={"carriers": {"raw": 453}},
            )
        ],
    )
    assert cards
    cards[0].run_heuristic_verdict()
    assert cards[0].verdict == VERDICT_REJECT
    assert "cohort_total" in (cards[0].verdict_reason or "")


def test_heuristic_rejects_when_source_text_matches_cohort_pattern():
    source = (
        "Methods: Total cohort comprised 453 patients across centers. "
        "Variant frequencies are reported below.\nVariant counts: 5, 6, 8."
    )
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[_variant(carriers=453, outlier_flags={"carriers": {"raw": 453}})],
        source_text=source,
    )
    cards[0].run_heuristic_verdict()
    assert cards[0].verdict == VERDICT_REJECT


def test_heuristic_rejects_when_extracted_value_matches_explicit_N():
    """The source 'N=137 patients' matches the N=<digits> reject pattern;
    rule 2 fires before rule 3, but the outcome (REJECT) is what matters."""
    source = "Discussion: We studied N=137 patients overall."
    cards = build_evidence_cards(
        pmid="16453024",
        variants=[_variant(carriers=137, outlier_flags={"carriers": {"raw": 137}})],
        source_text=source,
    )
    cards[0].run_heuristic_verdict()
    assert cards[0].verdict == VERDICT_REJECT
    # Either rule 2 (N=N regex) or rule 3 (explicit N=value match) is fine
    assert cards[0].verdict_reason  # non-empty


def test_heuristic_withholds_without_positive_evidence():
    source = "Family A: variant p.Gly1Asp present in proband and 4 relatives."
    cards = build_evidence_cards(
        pmid="11111111",
        variants=[_variant(carriers=200, outlier_flags={"carriers": {"raw": 200}})],
        source_text=source,
    )
    cards[0].run_heuristic_verdict()
    # Heuristic never auto-confirms; defer to downstream verifier
    assert cards[0].verdict == VERDICT_WITHHOLD


def test_as_dict_roundtrip_shape():
    """as_dict() should return a JSON-serializable mapping for downstream."""
    import json

    cards = build_evidence_cards(
        pmid="29622001",
        variants=[
            _variant(
                carriers=453,
                provenance={"carriers_count_type": "cohort_total"},
                classifier_flags={"carriers": {"raw": 453}},
            )
        ],
    )
    cards[0].run_heuristic_verdict()
    d = cards[0].as_dict()
    # Round-trip through json — catches non-serializable values
    json.dumps(d)
    assert d["pmid"] == "29622001"
    assert d["field"] == "carriers"
    assert d["verdict"] == VERDICT_REJECT


# ---- thresholds ------------------------------------------------------------


def test_large_value_threshold_is_configurable():
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[_variant(carriers=60, provenance={"carriers_count_type": "unknown"})],
        large_value_threshold=100,  # raise the bar
    )
    # 60 < 100 → no unknown_provenance_large trigger
    assert cards == []


def test_gold_error_threshold_is_configurable():
    cards = build_evidence_cards(
        pmid="29622001",
        variants=[_variant(carriers=12)],
        gold_counts_by_variant={"p.Gly1Asp": {"carriers": 5}},
        gold_error_threshold=5,  # lower the bar
    )
    # |12-5|=7 >= 5
    assert cards and "gold_error_ge_5" in cards[0].triggers
