"""Tests for the opt-in Variant_Browser publish step (cli/gvf_run.py).

Covers review-repo resolution and the best-effort ``step_publish_review``
contract: it builds the right ``gvf_publish.sh`` command, and it never raises —
a missing repo, a non-zero exit, or a timeout all warn and return False so the
GVF run is unaffected.
"""

import subprocess
import types

import pytest

import cli.gvf_run as gvf_run


def _fake_repo(tmp_path):
    """A directory shaped like a Variant_Browser checkout with the publish script."""
    repo = tmp_path / "Variant_Browser"
    (repo / "scripts").mkdir(parents=True)
    (repo / "scripts" / "gvf_publish.sh").write_text("#!/usr/bin/env bash\nexit 0\n")
    return repo


def test_find_review_repo_explicit(tmp_path):
    repo = _fake_repo(tmp_path)
    assert gvf_run._find_review_repo(repo) == repo


def test_find_review_repo_env(tmp_path, monkeypatch):
    repo = _fake_repo(tmp_path)
    monkeypatch.setenv("GVF_REVIEW_REPO", str(repo))
    assert gvf_run._find_review_repo(None) == repo


def test_find_review_repo_missing(tmp_path, monkeypatch):
    monkeypatch.delenv("GVF_REVIEW_REPO", raising=False)
    monkeypatch.delenv("VARIANT_BROWSER_DIR", raising=False)
    # Point the sibling-default fallback at a parent with no Variant_Browser so
    # every candidate misses; then None means "warn and skip", not a real repo.
    monkeypatch.setattr(gvf_run, "REPO_ROOT", tmp_path / "repo")
    assert gvf_run._find_review_repo(tmp_path / "nope") is None


def test_publish_builds_command_with_disease(tmp_path, monkeypatch):
    repo = _fake_repo(tmp_path)
    db = tmp_path / "KCNH2.db"
    db.write_text("")
    captured = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        captured["env"] = kwargs["env"]
        return types.SimpleNamespace(returncode=0, stdout="Publishing KCNH2", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    manifest = tmp_path / "KCNH2.txt"
    manifest.write_text("12345\n")
    ok = gvf_run.step_publish_review(
        gene="KCNH2",
        db=db,
        disease="Long QT type 2",
        review_repo=repo,
        pmid_file=manifest,
    )
    assert ok is True
    assert captured["cmd"] == [
        "bash",
        str(repo / "scripts" / "gvf_publish.sh"),
        "KCNH2",
        str(db),
        "Long QT type 2",
    ]
    assert captured["env"]["GVF_PMID_FILE"] == str(manifest.resolve())


def test_publish_omits_disease_when_absent(tmp_path, monkeypatch):
    repo = _fake_repo(tmp_path)
    db = tmp_path / "SCN5A.db"
    db.write_text("")
    captured = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setenv("GVF_FULL_DB_REPLACE", "1")
    gvf_run.step_publish_review(gene="SCN5A", db=db, disease=None, review_repo=repo)
    assert captured["cmd"][-1] == str(db)  # no trailing disease arg


def test_publish_refuses_without_manifest_or_explicit_full_replace(
    tmp_path, monkeypatch
):
    repo = _fake_repo(tmp_path)
    monkeypatch.delenv("GVF_FULL_DB_REPLACE", raising=False)

    def boom(*args, **kwargs):
        raise AssertionError("subprocess.run should not be called")

    monkeypatch.setattr(subprocess, "run", boom)
    assert (
        gvf_run.step_publish_review(
            gene="SCN5A",
            db=tmp_path / "SCN5A.db",
            disease=None,
            review_repo=repo,
        )
        is False
    )


def test_publish_missing_repo_is_noop(tmp_path, monkeypatch):
    monkeypatch.delenv("GVF_REVIEW_REPO", raising=False)
    monkeypatch.delenv("VARIANT_BROWSER_DIR", raising=False)

    def boom(*a, **k):  # subprocess must never be called when no repo is found
        raise AssertionError("subprocess.run should not be called")

    monkeypatch.setattr(subprocess, "run", boom)
    ok = gvf_run.step_publish_review(
        gene="KCNH2", db=tmp_path / "x.db", disease=None, review_repo=tmp_path / "nope"
    )
    assert ok is False


def test_publish_nonzero_exit_warns_not_raises(tmp_path, monkeypatch):
    repo = _fake_repo(tmp_path)
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **k: types.SimpleNamespace(
            returncode=1, stdout="", stderr="pair not found"
        ),
    )
    ok = gvf_run.step_publish_review(
        gene="KCNH2", db=tmp_path / "x.db", disease=None, review_repo=repo
    )
    assert ok is False  # did not raise


def test_publish_timeout_is_swallowed(tmp_path, monkeypatch):
    repo = _fake_repo(tmp_path)

    def timeout(*a, **k):
        raise subprocess.TimeoutExpired(cmd="gvf_publish.sh", timeout=1)

    monkeypatch.setattr(subprocess, "run", timeout)
    ok = gvf_run.step_publish_review(
        gene="KCNH2", db=tmp_path / "x.db", disease=None, review_repo=repo
    )
    assert ok is False
