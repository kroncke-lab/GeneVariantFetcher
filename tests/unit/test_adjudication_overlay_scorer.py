"""The scorer must consume the human adjudication overlay when supplied.

Regression guard for the closed feedback loop: corrections produced by
scripts/ingest_review_adjudications.py were previously never read back by the
scorer. apply_adjudication_overlay folds verdicts into the comparison rows
before recall/MAE/precision are computed.
"""

import csv

from cli.compare_variants import (
    ComparisonRow,
    apply_adjudication_overlay,
    compute_precision_summary,
    compute_recall_summary,
    compute_rows_mae,
    load_adjudication_overlay,
)


def _row(**kw) -> ComparisonRow:
    base = dict(
        pmid="1",
        excel_variant_raw="p.V1M",
        excel_variant_norm="p.V1M",
        sqlite_variant_raw="p.V1M",
        sqlite_variant_norm="p.V1M",
        match_type="exact",
        match_score=1.0,
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


def _overlay(**over):
    base = {
        "gene": "KCNQ1",
        "pmid": "1",
        "source_notation": "p.V1M",
        "verdict": "",
        "action": "",
        "corrected_affected": "",
        "corrected_unaffected": "",
        "corrected_total": "",
    }
    base.update(over)
    return {("1", "p.V1M"): base}


def test_empty_overlay_is_noop():
    rows = [_row(excel_affected=3, sqlite_affected=5, affected_diff=-2)]
    out = apply_adjudication_overlay(rows, {})
    assert out is rows


def test_count_override_corrects_extracted_counts_and_mae():
    # Extracted 92 but gold says 2; before override MAE = 90.
    rows = [
        _row(
            excel_affected=2,
            sqlite_affected=92,
            affected_diff=-90,
            count_mismatch=True,
        )
    ]
    assert compute_rows_mae(rows)["affected"]["mae"] == 90.0

    overlay = _overlay(action="count_override", corrected_affected="2")
    out = apply_adjudication_overlay(rows, overlay)
    assert out[0].sqlite_affected == 2
    assert out[0].affected_diff == 0
    assert out[0].count_mismatch is False
    assert compute_rows_mae(out)["affected"]["mae"] == 0.0


def test_wrong_paper_excludes_row():
    rows = [_row(excel_affected=3, sqlite_affected=3, affected_diff=0)]
    out = apply_adjudication_overlay(rows, _overlay(verdict="wrong_paper"))
    assert out == []


def test_wrong_paper_excludes_missed_gold_rows_from_recall():
    # A paper the curator judged out of scope (wrong_paper) that carries BOTH an
    # extracted+matched variant AND a gold variant GVF missed. The missed-gold
    # row has no DB source_notation, so it never keys into the overlay -- only
    # the matched variant's notation appears there. Excluding the paper must drop
    # BOTH rows; dropping only the keyed one leaves the missed-gold variant in
    # the recall denominator for a paper that is out of scope for recall.
    matched = _row(
        pmid="1",
        excel_variant_raw="p.V1M",
        excel_variant_norm="p.V1M",
        sqlite_variant_raw="p.V1M",
        sqlite_variant_norm="p.V1M",
        match_type="exact",
        excel_affected=3,
        sqlite_affected=3,
        affected_diff=0,
    )
    missed_gold = _row(
        pmid="1",
        excel_variant_raw="p.V2M",
        excel_variant_norm="p.V2M",
        sqlite_variant_raw=None,
        sqlite_variant_norm=None,
        match_type="none",
        match_score=None,
        missing_in_sqlite=True,
        excel_affected=2,
    )
    rows = [matched, missed_gold]

    # Both variants sit in the gold recall denominator before adjudication.
    before = compute_recall_summary(rows)
    assert before["variant_rows"]["gold"] == 2
    assert before["unique_variants"]["gold"] == 2

    # The overlay only keys the matched variant (the one with a DB notation).
    overlay = _overlay(verdict="wrong_paper")  # keyed ("1", "p.V1M")
    out = apply_adjudication_overlay(rows, overlay)

    # The whole paper drops: no rows survive and the recall denominator empties.
    assert out == []
    after = compute_recall_summary(out)
    assert after["variant_rows"]["gold"] == 0
    assert after["unique_variants"]["gold"] == 0
    assert after["pmids"]["gold"] == 0


def test_confirmed_true_positive_extra_leaves_precision_denominator():
    # A DB-only extra on a gold PMID: precision proxy counts it as a possible FP.
    extra = _row(
        excel_variant_raw="",
        excel_variant_norm="",
        match_type="none",
        missing_in_excel=True,
        sqlite_affected=4,
    )
    gold_match = _row(pmid="1", excel_affected=3, sqlite_affected=3, affected_diff=0)
    rows = [gold_match, extra]
    assert compute_precision_summary(rows)["counted_extra_on_gold_pmids"] == 1

    overlay = {("1", "p.V1M"): {"source_notation": "p.V1M", "verdict": "confirm"}}
    out = apply_adjudication_overlay(rows, overlay)
    # The confirmed-real extra is removed from the false-positive set.
    assert compute_precision_summary(out)["counted_extra_on_gold_pmids"] == 0


def test_false_positive_extra_is_kept_as_confirmed_fp():
    extra = _row(
        excel_variant_norm="",
        match_type="none",
        missing_in_excel=True,
        sqlite_affected=4,
    )
    overlay = {("1", "p.V1M"): {"source_notation": "p.V1M", "verdict": "wrong_variant"}}
    out = apply_adjudication_overlay([extra], overlay)
    assert len(out) == 1


def test_load_overlay_round_trips_csv(tmp_path):
    path = tmp_path / "KCNQ1_review_adjudications.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene", "pmid", "source_notation", "verdict", "action"])
        w.writerow(["KCNQ1", "1", "p.V1M", "correct_counts", "count_override"])
    overlay = load_adjudication_overlay(path)
    # Keyed on (pmid, canonical variant) -- same canonicalization the scorer
    # applies to DB rows, so a loaded overlay lines up with real ComparisonRows.
    assert len(overlay) == 1
    (key,) = overlay
    assert key[0] == "1"
    assert overlay[key]["action"] == "count_override"
