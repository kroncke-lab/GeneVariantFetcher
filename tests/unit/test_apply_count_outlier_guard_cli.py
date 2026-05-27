"""Smoke tests for scripts.apply_count_outlier_guard CLI."""

import json
import sys

import scripts.apply_count_outlier_guard as cli


def _make_extraction_json(path, variants):
    payload = {"variants": variants, "extraction_metadata": {}}
    path.write_text(json.dumps(payload), encoding="utf-8")


def _seed_extractions(extractions_dir):
    """Create one PMID with a clear outlier, one with no outliers."""
    extractions_dir.mkdir(parents=True, exist_ok=True)

    # PMID with study-wide-N reuse on one row
    outlier_variants = [
        {"patients": {"count": 453}},
        {"patients": {"count": 4}},
        {"patients": {"count": 5}},
        {"patients": {"count": 6}},
        {"patients": {"count": 3}},
        {"patients": {"count": 8}},
    ]
    _make_extraction_json(
        extractions_dir / "KCNQ1_PMID_29622001.json", outlier_variants
    )

    # PMID with consistent counts; should not flag
    normal_variants = [
        {"patients": {"count": 7}},
        {"patients": {"count": 8}},
        {"patients": {"count": 6}},
        {"patients": {"count": 10}},
        {"patients": {"count": 9}},
    ]
    _make_extraction_json(extractions_dir / "KCNQ1_PMID_99999999.json", normal_variants)


def test_dry_run_detects_without_modifying(tmp_path, monkeypatch, capsys):
    extractions = tmp_path / "extractions"
    _seed_extractions(extractions)
    report = tmp_path / "report.json"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "apply_count_outlier_guard.py",
            "--extractions",
            str(extractions),
            "--policy",
            "off",
            "--report-out",
            str(report),
            "--dry-run",
        ],
    )
    rc = cli.main()
    assert rc == 0

    summary = json.loads(report.read_text())
    assert summary["totals"]["pmids_processed"] == 2
    assert summary["totals"]["pmids_with_outliers"] == 1
    # Outlier PMID present in rows; normal PMID absent (no annotations)
    files = {r["file"] for r in summary["rows"]}
    assert "KCNQ1_PMID_29622001.json" in files

    # JSON on disk untouched
    outlier_json = json.loads((extractions / "KCNQ1_PMID_29622001.json").read_text())
    assert outlier_json["variants"][0]["patients"]["count"] == 453
    assert "count_outlier_guard" not in outlier_json


def test_clear_policy_zeros_and_backs_up(tmp_path, monkeypatch):
    extractions = tmp_path / "extractions"
    _seed_extractions(extractions)
    report = tmp_path / "report.json"
    backup = tmp_path / "backup"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "apply_count_outlier_guard.py",
            "--extractions",
            str(extractions),
            "--policy",
            "clear",
            "--backup-dir",
            str(backup),
            "--report-out",
            str(report),
        ],
    )
    rc = cli.main()
    assert rc == 0

    summary = json.loads(report.read_text())
    assert summary["totals"]["total_cleared_values"] == 1

    # Outlier file mutated; raw preserved in flags
    outlier_json = json.loads((extractions / "KCNQ1_PMID_29622001.json").read_text())
    assert outlier_json["variants"][0]["patients"]["count"] is None
    assert outlier_json["variants"][0]["count_outlier_flags"]["carriers"]["raw"] == 453
    assert outlier_json["count_outlier_guard"]["policy"] == "clear"

    # Backup contains the original
    backup_payload = json.loads((backup / "KCNQ1_PMID_29622001.json").read_text())
    assert backup_payload["variants"][0]["patients"]["count"] == 453


def test_flag_policy_preserves_value_and_records_flag(tmp_path, monkeypatch):
    extractions = tmp_path / "extractions"
    _seed_extractions(extractions)
    report = tmp_path / "report.json"
    backup = tmp_path / "backup"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "apply_count_outlier_guard.py",
            "--extractions",
            str(extractions),
            "--policy",
            "flag",
            "--backup-dir",
            str(backup),
            "--report-out",
            str(report),
        ],
    )
    rc = cli.main()
    assert rc == 0

    outlier_json = json.loads((extractions / "KCNQ1_PMID_29622001.json").read_text())
    # Value preserved (not cleared)
    assert outlier_json["variants"][0]["patients"]["count"] == 453
    # Flag metadata present
    assert outlier_json["variants"][0]["count_outlier_flags"]["carriers"]["raw"] == 453
    assert outlier_json["count_outlier_guard"]["policy"] == "flag"
