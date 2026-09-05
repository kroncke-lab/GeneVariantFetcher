"""Protect the diagnostic's denominators, missing-value semantics, and lanes."""

from scripts.recall_audit.phenotype_failure_panel import summarize_paper


def count_row(status, gold, raw, measure="affected"):
    delta = (raw if raw is not None else 0) - gold
    return dict(
        status=status,
        gold_count=gold,
        automated_count_raw=raw,
        measure=measure,
        difference=delta,
        absolute_difference=abs(delta),
    )


def test_missing_gold_zero_is_not_a_supplied_exact_or_count_false_positive():
    rows = [count_row("abstained", 0, None), count_row("identity_miss", 7, None)]
    result = summarize_paper(
        rows, dict(gene="TEST", pmid="1", tp=1, fp=8, fn=1), "run", "paper"
    )
    assert result["affected_supplied_exact"] == 0
    assert result["affected_missing_gold_zero_fields"] == 1
    assert result["affected_identity_miss_error"] == 7
    assert result["affected_overcount_units"] == 0
    assert result["fp"] == 8
    assert rows[0]["automated_count_raw"] is None


def test_wrong_supplied_and_missing_counts_have_separate_error_terms():
    rows = [
        count_row("supplied", 2, 5),
        count_row("supplied", 4, 1),
        count_row("abstained", 6, None),
        count_row("supplied", 0, 0, "unaffected"),
    ]
    result = summarize_paper(
        rows, dict(gene="TEST", pmid="1", tp=3, fp=0, fn=0), "run", "paper"
    )
    assert result["phenotype_capture"] == "partial"
    assert result["phenotype_abs_error"] == 12
    assert result["affected_overcount_units"] == 3
    assert result["affected_supplied_undercount_units"] == 3
    assert result["affected_abstained_error"] == 6
    assert result["unaffected_supplied_exact"] == 1
    assert result["phenotype_positive_supplied_fields"] == 2


def test_same_pmid_in_two_lanes_remains_two_distinct_records():
    paper = dict(gene="TEST", pmid="1", tp=1, fp=3, fn=0)
    a = summarize_paper([count_row("supplied", 4, 4)], paper, "legacy", "linkage")
    b = summarize_paper(
        [count_row("identity_miss", 4, None)],
        {**paper, "tp": 0, "fp": 0, "fn": 1},
        "new",
        "paper",
    )
    assert a["run_id"] != b["run_id"]
    assert a["score_lane"] != b["score_lane"]
    assert a["phenotype_abs_error"] == 0
    assert b["phenotype_abs_error"] == 4


def test_unasserted_phenotype_is_not_all_na_and_carriers_do_not_supply_it():
    result = summarize_paper(
        [count_row("supplied", 4, 4, "carriers")],
        dict(gene="TEST", pmid="1", tp=1, fp=0, fn=0),
        "run",
        "paper",
    )
    assert result["phenotype_capture"] == "none_asserted"
    assert result["phenotype_supplied_fields"] == 0
