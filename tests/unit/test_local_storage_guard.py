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

import pipeline.steps as steps
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


def test_state_is_linked_for_a_real_local_directory(fake_repo):
    (fake_repo / "corpus").mkdir()

    assert external_path_state("corpus") == "linked"


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
    with pytest.raises(LocalStorageError, match="no 'corpus' entry"):
        require_external_storage(fake_repo / "corpus")

    assert not (fake_repo / "corpus").exists(), "the guard must not create anything"


def test_dangling_link_is_refused_with_the_mount_reason(fake_repo):
    (fake_repo / "corpus").symlink_to(fake_repo / "never_mounted")

    with pytest.raises(LocalStorageError, match="target is unreachable"):
        require_external_storage(fake_repo / "corpus" / "KCNH2")


def test_error_names_the_volume_and_the_escape_hatch(fake_repo):
    with pytest.raises(LocalStorageError) as excinfo:
        require_external_storage(fake_repo / "corpus")

    message = str(excinfo.value)
    assert "Ezekers" in message
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
    assert "refusing to create" in capsys.readouterr().err


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
    (tmp_path / "corpus").symlink_to(tmp_path / "never_mounted")

    with caplog.at_level(logging.WARNING, logger=steps.logger.name):
        assert steps._resolve_corpus_dir() is None

    assert "corpus reuse DISABLED" in caplog.text
    assert "Ezekers" in caplog.text


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
