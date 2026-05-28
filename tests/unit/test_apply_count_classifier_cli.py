"""Smoke tests for scripts.apply_count_classifier CLI."""

import json
import sys

import scripts.apply_count_classifier as cli


def _write(path, variants):
    path.write_text(
        json.dumps({"variants": variants, "extraction_metadata": {}}),
        encoding="utf-8",
    )


def _seed(extractions_dir):
    extractions_dir.mkdir(parents=True, exist_ok=True)

    # Misclassified: count_type=cohort_total but value populated
    misclassified = [
        {
            "patients": {"count": 453},
            "penetrance_data": {"total_carriers_observed": 453},
            "count_provenance": {
                "carriers_column_label": "Total N",
                "carriers_count_type": "cohort_total",
            },
        }
    ]
    _write(extractions_dir / "KCNQ1_PMID_29622001.json", misclassified)

    # Properly classified: count_type=per_variant_carrier; do not touch
    proper = [
        {
            "patients": {"count": 5},
            "penetrance_data": {"total_carriers_observed": 5},
            "count_provenance": {
                "carriers_column_label": "No. of patients",
                "carriers_count_type": "per_variant_carrier",
            },
        }
    ]
    _write(extractions_dir / "KCNQ1_PMID_99999999.json", proper)

    # Legacy: no provenance at all; classifier must skip silently
    legacy = [{"patients": {"count": 200}}]
    _write(extractions_dir / "KCNQ1_PMID_88888888.json", legacy)


def test_dry_run_detects_only_misclassified_pmids(tmp_path, monkeypatch):
    extractions = tmp_path / "extractions"
    _seed(extractions)
    report = tmp_path / "report.json"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "apply_count_classifier.py",
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
    assert summary["totals"]["pmids_processed"] == 3
    assert summary["totals"]["pmids_with_misclassified"] == 1
    assert summary["totals"]["pmids_with_provenance"] == 2  # legacy excluded


def test_clear_policy_nulls_both_mirrors_and_backs_up(tmp_path, monkeypatch):
    extractions = tmp_path / "extractions"
    _seed(extractions)
    report = tmp_path / "report.json"
    backup = tmp_path / "backup"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "apply_count_classifier.py",
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

    # Misclassified file mutated
    cleared = json.loads((extractions / "KCNQ1_PMID_29622001.json").read_text())
    assert cleared["variants"][0]["patients"]["count"] is None
    assert cleared["variants"][0]["penetrance_data"]["total_carriers_observed"] is None
    assert cleared["variants"][0]["count_classifier_flags"]["carriers"]["raw"] == 453
    assert cleared["count_classifier"]["policy"] == "clear"

    # Proper file untouched
    proper = json.loads((extractions / "KCNQ1_PMID_99999999.json").read_text())
    assert proper["variants"][0]["patients"]["count"] == 5
    assert "count_classifier" not in proper

    # Legacy file untouched (no provenance to act on)
    legacy = json.loads((extractions / "KCNQ1_PMID_88888888.json").read_text())
    assert legacy["variants"][0]["patients"]["count"] == 200
    assert "count_classifier" not in legacy

    # Backup contains the original
    backup_payload = json.loads((backup / "KCNQ1_PMID_29622001.json").read_text())
    assert backup_payload["variants"][0]["patients"]["count"] == 453
