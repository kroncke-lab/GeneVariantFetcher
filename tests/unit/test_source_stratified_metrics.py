"""Tests for scripts/recall_audit/source_stratified_metrics.py.

The point of the script is that a single headline recall number conflates
acquisition failures (the paper was never downloaded) with extraction failures
(the text was there and we missed it). These pin the split.
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "source_stratified_metrics",
    REPO / "scripts" / "recall_audit" / "source_stratified_metrics.py",
)
ssm = importlib.util.module_from_spec(_spec)
sys.modules["source_stratified_metrics"] = ssm
_spec.loader.exec_module(ssm)

FIELDS = [
    "pmid",
    "excel_variant_raw",
    "sqlite_variant_raw",
    "match_type",
    "excel_carriers_total",
    "sqlite_carriers_total",
    "excel_affected",
    "sqlite_affected",
    "excel_unaffected",
    "sqlite_unaffected",
    "missing_in_sqlite",
    "missing_in_excel",
]


def _rows(path: Path, rows: list[dict]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})
    return path


def _paper(corpus: Path, gene: str, pmid: str, *, body="x" * 20000, links=""):
    d = corpus / gene / pmid
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{pmid}_FULL_CONTEXT.md").write_text(body + links)
    return d


# ---------------------------------------------------------------------------
# Source classification
# ---------------------------------------------------------------------------


def test_a_full_text_paper_with_no_links_is_source_complete(tmp_path):
    _paper(tmp_path, "KCNQ1", "1")
    src = ssm.classify_paper(tmp_path, "KCNQ1", "1")
    assert src.has_full_text and src.is_source_complete
    assert src.reason() == "complete"


def test_the_abstract_only_fallback_is_not_full_text(tmp_path):
    _paper(
        tmp_path,
        "KCNQ1",
        "1",
        body=ssm.ABSTRACT_ONLY_MARKER + "\n\n> WARNING\n" + "x" * 20000,
    )
    src = ssm.classify_paper(tmp_path, "KCNQ1", "1")
    assert not src.has_full_text
    assert src.reason() == "abstract-only fallback"


def test_a_stub_shorter_than_a_body_is_not_full_text(tmp_path):
    _paper(tmp_path, "KCNQ1", "1", body="tiny")
    assert ssm.classify_paper(tmp_path, "KCNQ1", "1").reason() == (
        "source too short to be a body"
    )


def test_a_missing_paper_is_reported_not_crashed(tmp_path):
    src = ssm.classify_paper(tmp_path, "KCNQ1", "404")
    assert not src.exists and src.reason() == "no source on disk"


def test_an_advertised_but_unfetched_supplement_breaks_completeness(tmp_path):
    _paper(tmp_path, "KCNQ1", "1", links="\n\n### T1\n\n_link_: mmc1.xls\n")
    src = ssm.classify_paper(tmp_path, "KCNQ1", "1")
    assert src.has_full_text
    assert not src.is_source_complete
    assert src.missing_supplements == 1
    assert "advertised supplement" in src.reason()


def test_a_fetched_supplement_restores_completeness(tmp_path):
    d = _paper(tmp_path, "KCNQ1", "1", links="\n\n### T1\n\n_link_: mmc1.xls\n")
    supp = d / "1_supplements"
    supp.mkdir()
    (supp / "mmc1.xls").write_bytes(b"data")
    assert ssm.classify_paper(tmp_path, "KCNQ1", "1").is_source_complete


def test_a_supplement_extracted_into_a_nested_archive_dir_still_counts(tmp_path):
    """Europe PMC archives unpack into <stem>/, so the scan must recurse."""
    d = _paper(tmp_path, "KCNQ1", "1", links="\n\n### T1\n\n_link_: mmc1.xls\n")
    nested = d / "1_supplements" / "PMC9_supplements"
    nested.mkdir(parents=True)
    (nested / "mmc1.xls").write_bytes(b"data")
    assert ssm.classify_paper(tmp_path, "KCNQ1", "1").is_source_complete


# ---------------------------------------------------------------------------
# Stratification
# ---------------------------------------------------------------------------


def test_the_two_subsets_partition_all_gold(tmp_path):
    corpus = tmp_path / "corpus"
    _paper(corpus, "KCNQ1", "1")  # complete
    _paper(corpus, "KCNQ1", "2", body="tiny")  # incomplete
    rows = _rows(
        tmp_path / "discrepancies.csv",
        [
            {
                "pmid": "1",
                "match_type": "exact",
                "excel_carriers_total": "5",
                "sqlite_carriers_total": "5",
            },
            {"pmid": "1", "match_type": "none", "missing_in_sqlite": "True"},
            {
                "pmid": "2",
                "match_type": "exact",
                "excel_carriers_total": "3",
                "sqlite_carriers_total": "4",
            },
        ],
    )
    out = ssm.score_gene("KCNQ1", rows, corpus)
    s = out["strata"]
    assert s["all_gold"]["gold_rows"] == 3
    assert s["source_complete"]["gold_rows"] == 2
    assert s["source_incomplete"]["gold_rows"] == 1
    # Disjoint and exhaustive.
    assert (
        s["source_complete"]["gold_rows"] + s["source_incomplete"]["gold_rows"]
        == s["all_gold"]["gold_rows"]
    )
    assert (
        s["source_complete"]["matched_rows"] + s["source_incomplete"]["matched_rows"]
        == s["all_gold"]["matched_rows"]
    )


def test_mae_and_rmse_are_computed_over_matched_rows(tmp_path):
    corpus = tmp_path / "corpus"
    _paper(corpus, "KCNQ1", "1")
    rows = _rows(
        tmp_path / "d.csv",
        [
            {
                "pmid": "1",
                "match_type": "exact",
                "excel_carriers_total": "10",
                "sqlite_carriers_total": "13",
            },  # err +3
            {
                "pmid": "1",
                "match_type": "exact",
                "excel_carriers_total": "10",
                "sqlite_carriers_total": "9",
            },  # err -1
        ],
    )
    c = ssm.score_gene("KCNQ1", rows, corpus)["strata"]["all_gold"]["counts"]
    assert c["carriers"]["n_matched"] == 2
    assert c["carriers"]["mae"] == 2.0  # (3 + 1) / 2
    assert abs(c["carriers"]["rmse"] - (10 / 2) ** 0.5) < 1e-9  # sqrt((9+1)/2)


def test_extras_only_count_on_gold_curated_pmids(tmp_path):
    """A DB row on a PMID gold never touched is not judgeable."""
    corpus = tmp_path / "corpus"
    _paper(corpus, "KCNQ1", "1")
    _paper(corpus, "KCNQ1", "99")
    rows = _rows(
        tmp_path / "d.csv",
        [
            {
                "pmid": "1",
                "match_type": "exact",
                "excel_carriers_total": "1",
                "sqlite_carriers_total": "1",
            },
            {"pmid": "1", "missing_in_excel": "True", "sqlite_carriers_total": "2"},
            {"pmid": "99", "missing_in_excel": "True", "sqlite_carriers_total": "7"},
        ],
    )
    s = ssm.score_gene("KCNQ1", rows, corpus)["strata"]["all_gold"]
    assert s["extra_rows_on_gold_pmids"] == 1  # PMID 99 excluded
    assert s["counted_extra_rows_on_gold_pmids"] == 1
    assert s["precision_vs_counted_gold_pmids"] == 0.5  # 1 matched / (1 + 1)


def test_an_extra_without_any_count_is_excluded_from_headline_precision(tmp_path):
    """Zero-count linkage rows dominate the extras and are not real FPs."""
    corpus = tmp_path / "corpus"
    _paper(corpus, "KCNQ1", "1")
    rows = _rows(
        tmp_path / "d.csv",
        [
            {
                "pmid": "1",
                "match_type": "exact",
                "excel_carriers_total": "1",
                "sqlite_carriers_total": "1",
            },
            {"pmid": "1", "missing_in_excel": "True"},  # no counts at all
        ],
    )
    s = ssm.score_gene("KCNQ1", rows, corpus)["strata"]["all_gold"]
    assert s["extra_rows_on_gold_pmids"] == 1
    assert s["counted_extra_rows_on_gold_pmids"] == 0
    assert s["precision_vs_counted_gold_pmids"] == 1.0
    assert s["loose_precision_vs_gold_pmids"] == 0.5
