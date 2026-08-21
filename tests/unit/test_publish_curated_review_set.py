import json

import pytest

from scripts import publish_curated_review_set as publish


def test_active_review_scope_keeps_only_lead_approved_human_brca2_papers():
    manifest = publish.load_manifest()

    assert manifest["BRCA2"] == ["26833046", "26848529"]
    assert publish.load_review_scope_exclusions() == {
        ("BRCA2", "10398279"),
        ("BRCA2", "15365993"),
        ("BRCA2", "18489799"),
        ("BRCA2", "19944633"),
        ("BRCA2", "21356067"),
        ("BRCA2", "22655046"),
        ("BRCA2", "25802882"),
    }


def test_brca2_review_pmid_file_matches_active_scope():
    expected = publish.load_manifest()["BRCA2"]
    actual = [
        line.strip()
        for line in (publish.REVIEW_PMID_DIR / "BRCA2.txt").read_text().splitlines()
        if line.strip()
    ]

    assert actual == expected


def test_relative_db_override_is_resolved_before_cross_repo_publish(tmp_path):
    relative = tmp_path / "nested" / "BRCA2.db"
    relative.parent.mkdir()
    relative.touch()

    resolved = publish.db_path_for(
        "BRCA2", db_overrides={"BRCA2": relative}, extract_root=None
    ).resolve()

    assert resolved.is_absolute()
    assert resolved == relative.resolve()


def test_extract_root_uses_run_status_active_db_not_newest_backup(tmp_path):
    run = tmp_path / "BRCA2" / "run"
    run.mkdir(parents=True)
    active = run / "BRCA2.db"
    active.touch()
    backup = run / "BRCA2.before_fp_quarantine.db"
    backup.touch()
    (run / "RUN_STATUS.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "exit_code": 0,
                "stage_failures": [],
                "active_db": "BRCA2.db",
            }
        )
    )
    backup.touch()

    assert publish.latest_gene_db(tmp_path, "BRCA2") == active.resolve()


def test_extract_root_refuses_ambiguous_completed_runs(tmp_path):
    for name in ("first", "second"):
        run = tmp_path / "BRCA2" / name
        run.mkdir(parents=True)
        (run / "BRCA2.db").touch()
        (run / "RUN_STATUS.json").write_text(
            json.dumps(
                {
                    "status": "completed",
                    "exit_code": 0,
                    "stage_failures": [],
                    "active_db": "BRCA2.db",
                }
            )
        )

    with pytest.raises(SystemExit, match="multiple completed active DBs"):
        publish.latest_gene_db(tmp_path, "BRCA2")
