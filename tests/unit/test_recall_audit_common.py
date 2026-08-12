from pathlib import Path

import pytest

from scripts.recall_audit import common


def _status_doc(tmp_path: Path) -> tuple[Path, Path, Path]:
    score = tmp_path / "recall_metrics" / "current"
    score.mkdir(parents=True)
    db = tmp_path / "validation_runs" / "canonical_baseline" / "KCNH2.db"
    db.parent.mkdir(parents=True)
    db.touch()
    status = tmp_path / "RECALL_STATUS.md"
    status.write_text(
        "\n".join(
            [
                "Historical baseline: `/obsolete/KCNH2.db`.",
                "",
                f"Current local scored artifact: `{score}/`.",
                "",
                "## Current Canonical Baseline",
                f"- `{db}`",
            ]
        ),
        encoding="utf-8",
    )
    return status, score, db


def test_status_resolves_current_score_and_canonical_db(tmp_path, monkeypatch):
    status, score, db = _status_doc(tmp_path)
    monkeypatch.setattr(common, "STATUS_DOC", status)
    monkeypatch.delenv(common.RUN_DIR_ENV, raising=False)

    assert common.resolve_recall_score(None) == score
    assert common.resolve_gene_db("kcnh2") == db


def test_explicit_inputs_override_status(tmp_path, monkeypatch):
    status, _, _ = _status_doc(tmp_path)
    monkeypatch.setattr(common, "STATUS_DOC", status)
    explicit_score = tmp_path / "other-score"
    explicit_db = tmp_path / "other.db"

    assert common.resolve_recall_score(explicit_score) == explicit_score
    assert common.resolve_gene_db("KCNH2", explicit_db) == explicit_db


def test_run_directory_override_keeps_the_self_contained_layout(tmp_path, monkeypatch):
    status, score, db = _status_doc(tmp_path)
    run = tmp_path / "other-run"
    run.mkdir()
    monkeypatch.setattr(common, "STATUS_DOC", status)
    monkeypatch.setenv(common.RUN_DIR_ENV, str(run))

    assert common.resolve_recall_score(None) == run / "recall_score"
    assert common.resolve_gene_db("KCNH2") == run / "dbs" / "KCNH2.db"
    assert common.resolve_status_gene_db("KCNH2") == db
    assert common.resolve_recall_score(score) == score
    assert common.resolve_gene_db("KCNH2", db) == db


def test_missing_status_target_fails_with_actionable_message(tmp_path, monkeypatch):
    status = tmp_path / "RECALL_STATUS.md"
    status.write_text(
        "Current local scored artifact: `missing/current/`.\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(common, "STATUS_DOC", status)
    monkeypatch.setattr(common, "REPO_ROOT", tmp_path)
    monkeypatch.delenv(common.RUN_DIR_ENV, raising=False)

    with pytest.raises(SystemExit) as excinfo:
        common.resolve_recall_score(None)

    assert "does not exist" in str(excinfo.value)
    assert "--recall-score" in str(excinfo.value)
