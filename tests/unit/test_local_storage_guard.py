"""The `corpus/` external-storage guard must fail closed, not create local storage.

`corpus/` is an absolute symlink to an external volume on Brett's workstation.
The failure this file pins is the *missing* link (a fresh checkout, or a rename):
`build_source_corpus.py` called `out.mkdir(parents=True, exist_ok=True)`, which
would silently build a second corpus on the internal disk that no later run
reads, re-fetching paywalled source already cached on the volume.

This is the write-side counterpart to test_local_data_guard.py, which covers the
read side (the offline suite must not discover on-disk local data).

The guard is anchored on the real repo root, so these tests relocate that anchor
to `tmp_path` rather than touching the developer's actual `corpus/` link.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

import cli.gvf_run as gvf_run
import pipeline.steps as steps
import utils.gene_metadata as gene_metadata
import utils.local_storage as local_storage
from utils.local_storage import (
    ALLOW_LOCAL_ENV,
    LocalStorageError,
    external_path_state,
    require_external_storage,
)


@pytest.fixture
def fake_repo(monkeypatch, tmp_path):
    """Repoint the guard's repo-root anchor at an empty tmp_path checkout."""
    monkeypatch.setattr(local_storage, "_ANCHORS", (tmp_path,))
    monkeypatch.delenv(ALLOW_LOCAL_ENV, raising=False)
    return tmp_path


# ---------------------------------------------------------------------------
# external_path_state
# ---------------------------------------------------------------------------


def test_state_is_linked_for_a_mounted_symlink(fake_repo):
    target = fake_repo / "external"
    target.mkdir()
    (fake_repo / "corpus").symlink_to(target)

    assert external_path_state("corpus") == "linked"


def test_state_distinguishes_a_real_local_directory(fake_repo):
    (fake_repo / "corpus").mkdir()

    assert external_path_state("corpus") == "local"


def test_state_is_dangling_when_the_volume_is_unmounted(fake_repo):
    (fake_repo / "corpus").symlink_to(fake_repo / "never_mounted")

    assert external_path_state("corpus") == "dangling"


def test_state_is_absent_on_a_fresh_checkout(fake_repo):
    assert external_path_state("corpus") == "absent"


# ---------------------------------------------------------------------------
# require_external_storage
# ---------------------------------------------------------------------------


def test_absent_link_is_refused(fake_repo):
    """The regression: a fresh checkout must not grow a local corpus/ tree."""
    with pytest.raises(LocalStorageError, match="there is no corpus/ in the repo root"):
        require_external_storage(fake_repo / "corpus")

    assert not (fake_repo / "corpus").exists(), "the guard must not create anything"


def test_dangling_link_is_refused_with_the_mount_reason(fake_repo):
    (fake_repo / "corpus").symlink_to(fake_repo / "never_mounted")

    with pytest.raises(LocalStorageError, match="target is unreachable"):
        require_external_storage(fake_repo / "corpus" / "KCNH2")


def test_real_local_directory_requires_the_explicit_escape_hatch(fake_repo):
    (fake_repo / "corpus").mkdir()

    with pytest.raises(LocalStorageError, match="real local directory"):
        require_external_storage(fake_repo / "corpus" / "KCNH2")


def test_detached_drive_names_the_volume_to_attach(fake_repo):
    """A dangling link names its own volume, read from the link, not hardcoded."""
    (fake_repo / "corpus").symlink_to("/Volumes/SomeOtherDisk/research/corpus")

    with pytest.raises(LocalStorageError) as excinfo:
        require_external_storage(fake_repo / "corpus")

    message = str(excinfo.value)
    assert "SomeOtherDisk" in message
    assert "/Volumes/SomeOtherDisk" in message


def test_fresh_clone_is_not_told_to_mount_someone_elses_drive(fake_repo):
    """A collaborator has no Ezekers; naming it reads as a broken setup."""
    with pytest.raises(LocalStorageError) as excinfo:
        require_external_storage(fake_repo / "corpus")

    message = str(excinfo.value)
    assert "Volumes" not in message
    assert "mount" not in message.lower()
    assert "fresh clone" in message
    assert "ln -s" in message
    assert ALLOW_LOCAL_ENV in message


def test_nested_paths_under_a_mounted_link_pass(fake_repo):
    target = fake_repo / "external"
    target.mkdir()
    (fake_repo / "corpus").symlink_to(target)

    require_external_storage(fake_repo / "corpus" / "KCNH2" / "12345678")


def test_paths_outside_the_guarded_tree_pass(fake_repo, tmp_path):
    """An explicit GVF_CORPUS_DIR elsewhere is not the guard's business."""
    require_external_storage(tmp_path / "somewhere_else")
    require_external_storage(fake_repo / "results")
    require_external_storage(
        Path("/Volumes/Ezekers/ResearchData/GeneVariantFetcher/corpus")
    )


def test_corpus_prefixed_sibling_is_not_treated_as_corpus(fake_repo):
    """`is_relative_to` must not match a sibling that merely starts with 'corpus'."""
    require_external_storage(fake_repo / "corpus_backups")


def test_escape_hatch_allows_local_storage(fake_repo, monkeypatch):
    monkeypatch.setenv(ALLOW_LOCAL_ENV, "1")

    require_external_storage(fake_repo / "corpus")


@pytest.mark.parametrize("raw", ["0", "", "no", "maybe"])
def test_escape_hatch_stays_closed_for_non_truthy_values(fake_repo, monkeypatch, raw):
    monkeypatch.setenv(ALLOW_LOCAL_ENV, raw)

    with pytest.raises(LocalStorageError):
        require_external_storage(fake_repo / "corpus")


# ---------------------------------------------------------------------------
# build_source_corpus wiring
# ---------------------------------------------------------------------------


def test_corpus_builder_exits_nonzero_instead_of_creating_a_local_corpus(
    fake_repo, monkeypatch, capsys
):
    """`--apply` against a missing link must refuse before it writes anything."""
    import scripts.build_source_corpus as build_source_corpus

    out = fake_repo / "corpus"
    monkeypatch.setattr(
        "sys.argv", ["build_source_corpus.py", "--apply", "--out", str(out)]
    )

    assert build_source_corpus.main() == 2
    assert not out.exists(), "the builder must not create a local corpus"
    assert "cannot use corpus/" in capsys.readouterr().err


def test_corpus_builder_guards_before_rebuild_deletes_anything(
    fake_repo, monkeypatch, capsys
):
    """`--rebuild` rmtree's the corpus; the guard must run before that, not after."""
    import scripts.build_source_corpus as build_source_corpus

    out = fake_repo / "corpus"
    out.symlink_to(fake_repo / "never_mounted")
    monkeypatch.setattr(
        "sys.argv",
        ["build_source_corpus.py", "--apply", "--rebuild", "--out", str(out)],
    )

    assert build_source_corpus.main() == 2
    assert out.is_symlink(), "the dangling link must be left in place for the operator"


# ---------------------------------------------------------------------------
# pipeline.steps._resolve_corpus_dir — silent corpus loss
# ---------------------------------------------------------------------------


def test_unreachable_corpus_link_warns_instead_of_silently_running_cold(
    monkeypatch, tmp_path, caplog, allow_local_data_discovery
):
    """An unmounted volume disables corpus reuse; that must be said out loud.

    Needs `allow_local_data_discovery`: under the offline suite's default guard
    `_resolve_corpus_dir` returns None before it ever looks at `<repo>/corpus`,
    so without lifting it this test would pass without reaching the warning.
    """
    monkeypatch.setattr(steps, "__file__", str(tmp_path / "pipeline" / "steps.py"))
    monkeypatch.delenv("GVF_CORPUS_DIR", raising=False)
    (tmp_path / "corpus").symlink_to("/Volumes/SomeOtherDisk/research/corpus")

    with caplog.at_level(logging.WARNING, logger=steps.logger.name):
        assert steps._resolve_corpus_dir() is None

    assert "corpus reuse DISABLED" in caplog.text
    # Named from the link, not hardcoded, so it is right on any machine.
    assert "SomeOtherDisk" in caplog.text


def test_misconfigured_corpus_env_var_warns(monkeypatch, tmp_path, caplog):
    monkeypatch.setenv("GVF_CORPUS_DIR", str(tmp_path / "does_not_exist"))

    with caplog.at_level(logging.WARNING, logger=steps.logger.name):
        assert steps._resolve_corpus_dir() is None

    assert "corpus reuse DISABLED" in caplog.text


def test_absent_corpus_on_a_fresh_checkout_stays_quiet(
    monkeypatch, tmp_path, caplog, allow_local_data_discovery
):
    """Nothing is wrong with having no corpus yet — do not cry wolf."""
    monkeypatch.setattr(steps, "__file__", str(tmp_path / "pipeline" / "steps.py"))
    monkeypatch.delenv("GVF_CORPUS_DIR", raising=False)

    with caplog.at_level(logging.WARNING, logger=steps.logger.name):
        assert steps._resolve_corpus_dir() is None

    assert caplog.text == ""


# ---------------------------------------------------------------------------
# doctor — an unmounted volume must stop the run at Step 1
# ---------------------------------------------------------------------------


def _fresh_status() -> dict:
    return {"required": {}, "ok": True}


def test_doctor_blocks_when_the_volume_is_unmounted(fake_repo):
    """The user's scenario: the drive was there, now it is not. Stop at Step 1."""
    (fake_repo / "corpus").symlink_to("/Volumes/SomeOtherDisk/research/corpus")
    status = _fresh_status()

    gvf_run._apply_local_storage_checks(status)

    assert status["ok"] is False
    assert status["local_storage"] == {"corpus": "dangling"}
    blocked = [k for k, v in status["required"].items() if not v]
    assert len(blocked) == 1
    assert "SomeOtherDisk" in blocked[0]


def test_doctor_does_not_block_a_mounted_volume(fake_repo):
    target = fake_repo / "external"
    target.mkdir()
    (fake_repo / "corpus").symlink_to(target)
    status = _fresh_status()

    gvf_run._apply_local_storage_checks(status)

    assert status["ok"] is True
    assert status["required"] == {}


def test_doctor_requires_opt_in_for_a_real_local_corpus(fake_repo):
    (fake_repo / "corpus").mkdir()
    status = _fresh_status()

    gvf_run._apply_local_storage_checks(status)

    assert status["ok"] is False
    assert status["local_storage"] == {"corpus": "local"}
    assert any("GVF_ALLOW_LOCAL_CORPUS" in key for key in status["required"])


def test_doctor_honors_a_deliberate_corpus_override(fake_repo, monkeypatch, tmp_path):
    """A configured output bypasses the unused repo-local corpus path."""
    (fake_repo / "corpus").mkdir()
    monkeypatch.setenv("GVF_CORPUS_DIR", str(tmp_path / "configured-corpus"))
    status = _fresh_status()

    gvf_run._apply_local_storage_checks(status)

    assert status["ok"] is True
    assert status["required"] == {}
    assert status["local_storage"] == {"corpus": "local"}


def test_doctor_does_not_block_a_fresh_checkout_with_no_corpus(fake_repo):
    """A collaborator with no corpus yet must still be able to run."""
    status = _fresh_status()

    gvf_run._apply_local_storage_checks(status)

    assert status["ok"] is True
    assert status["local_storage"] == {"corpus": "absent"}


def test_doctor_check_preserves_an_already_failing_status(fake_repo):
    """Folding storage in must not resurrect an `ok` that something else cleared."""
    target = fake_repo / "external"
    target.mkdir()
    (fake_repo / "corpus").symlink_to(target)
    status = {"required": {"NCBI_EMAIL": False}, "ok": False}

    gvf_run._apply_local_storage_checks(status)

    assert status["ok"] is False


# ---------------------------------------------------------------------------
# gvf-run corpus sync — configured output must reach the builder
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("configured", [False, True])
def test_corpus_sync_passes_the_configured_output(monkeypatch, tmp_path, configured):
    commands = []

    class Completed:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(command, **_kwargs):
        commands.append(command)
        return Completed()

    monkeypatch.setattr(gvf_run.subprocess, "run", fake_run)
    if configured:
        output = tmp_path / "configured corpus"
        monkeypatch.setenv("GVF_CORPUS_DIR", str(output))
    else:
        output = None
        monkeypatch.delenv("GVF_CORPUS_DIR", raising=False)

    gvf_run.step_corpus_sync(tmp_path / "run", gene="KCNH2")

    assert len(commands) == 1
    command = commands[0]
    assert command[command.index("--assume-gene") + 1] == "KCNH2"
    if output is None:
        assert "--out" not in command
    else:
        assert command[command.index("--out") + 1] == str(output)


# ---------------------------------------------------------------------------
# utils.gene_metadata — a configured-but-unreachable VariantFeatures DB
# ---------------------------------------------------------------------------


def test_unreachable_configured_variantfeatures_db_warns(monkeypatch, tmp_path, caplog):
    """An unmounted volume must not silently downgrade to built-in metadata."""
    monkeypatch.setenv("VARIANTFEATURES_DB", "/Volumes/SomeOtherDisk/vf/variants.db")

    with caplog.at_level(logging.WARNING, logger=gene_metadata.logger.name):
        assert gene_metadata.default_variantfeatures_db_path() is None

    assert "VARIANTFEATURES_DB" in caplog.text
    assert "SomeOtherDisk" in caplog.text


def test_reachable_configured_variantfeatures_db_stays_quiet(
    monkeypatch, tmp_path, caplog
):
    db = tmp_path / "variants.db"
    db.write_bytes(b"SQLite format 3\x00")
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db))

    with caplog.at_level(logging.WARNING, logger=gene_metadata.logger.name):
        assert gene_metadata.default_variantfeatures_db_path() == db

    assert caplog.text == ""
