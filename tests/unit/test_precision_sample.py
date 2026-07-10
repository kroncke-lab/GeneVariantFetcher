"""Precision sampler: population selection, Wilson CI, and verdict tally."""

import csv

from cli.compare_variants import ComparisonRow
from scripts.precision_sample import (
    select_counted_extras,
    tally_verdicts,
    wilson_interval,
)


def _row(**kw) -> ComparisonRow:
    base = dict(
        pmid="1",
        excel_variant_raw="",
        excel_variant_norm="",
        sqlite_variant_raw="p.V1M",
        sqlite_variant_norm="V1M",
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


def test_select_counted_extras_matches_precision_population():
    rows = [
        # gold-matched row establishes PMID 1 as a gold PMID
        _row(pmid="1", excel_variant_norm="V2L", missing_in_excel=False),
        # counted extra on a gold PMID -> selected
        _row(pmid="1", missing_in_excel=True, sqlite_affected=4),
        # extra on a gold PMID but no count -> excluded
        _row(pmid="1", missing_in_excel=True),
        # counted extra on a NON-gold PMID -> excluded (not judgeable)
        _row(pmid="999", missing_in_excel=True, sqlite_carriers_total=7),
    ]
    extras = select_counted_extras(rows)
    assert len(extras) == 1
    assert extras[0].pmid == "1"
    assert extras[0].sqlite_affected == 4


def test_wilson_interval_basic():
    p, lo, hi = wilson_interval(8, 10)
    assert p == 0.8
    assert 0.0 <= lo < 0.8 < hi <= 1.0
    assert wilson_interval(0, 0) == (None, None, None)


def test_tally_verdicts_excludes_unsure(tmp_path):
    path = tmp_path / "filled.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pmid", "variant", "verdict"])
        for verdict in [
            "real_tp",
            "real_tp",
            "real_tp",
            "false_positive",
            "unsure",
            "",
        ]:
            w.writerow(["1", "p.V1M", verdict])
    result = tally_verdicts(path)
    assert result["true_positives"] == 3
    assert result["false_positives"] == 1
    assert result["unsure"] == 1
    assert result["blank"] == 1
    assert result["adjudicated"] == 4
    assert result["true_precision"] == 0.75
