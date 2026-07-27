"""The offline suite must not discover a developer's on-disk local data.

Each helper covered here guesses a path when nothing is configured. The guess
resolves in the main checkout and not in a side worktree or CI, so without the
`GVF_DISABLE_LOCAL_DATA` guard the same suite exercises different code paths in
different checkouts. These tests pin all three branches per helper: guard on
(guess skipped), guard off (guess used), and explicitly configured (wins over
the guard).

The "guard off" cases repoint each module's own repo-root anchor at `tmp_path`
so they never depend on what is actually checked out next to this repo.
"""

from __future__ import annotations

import pytest

import benchmarks.cold_start_eval.run_cold_start as run_cold_start
import cli.gvf_run as gvf_run
import pipeline.paper_final_check as paper_final_check
import pipeline.steps as steps
import utils.gene_metadata as gene_metadata
from utils.env_utils import (
    DISABLE_LOCAL_DATA_ENV,
    get_env_bool,
    local_data_discovery_disabled,
)

# `GVF_DISABLE_LOCAL_DATA=0` is the documented escape hatch for reproducing the
# old discover-whatever-is-on-disk behaviour (see tests/unit/conftest.py). Under
# it these assertions are inapplicable rather than broken, so skip them instead
# of handing the developer seven red tests. Evaluated at collection time, after
# the conftest's `pytest_configure` has installed the default.
requires_guard = pytest.mark.skipif(
    not get_env_bool(DISABLE_LOCAL_DATA_ENV, True),
    reason=f"{DISABLE_LOCAL_DATA_ENV} explicitly off — the offline guard is lifted",
)


# ---------------------------------------------------------------------------
# get_env_bool / local_data_discovery_disabled
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("1", True),
        ("true", True),
        ("TRUE", True),
        ("  Yes  ", True),
        ("on", True),
        ("0", False),
        ("false", False),
        ("no", False),
        ("off", False),
    ],
)
def test_get_env_bool_reads_both_polarities(monkeypatch, raw, expected):
    monkeypatch.setenv("GVF_SOME_FLAG", raw)
    assert get_env_bool("GVF_SOME_FLAG", default=not expected) is expected


@pytest.mark.parametrize("raw", ["", "   ", "maybe", "2"])
def test_get_env_bool_falls_back_on_unusable_values(monkeypatch, raw):
    """A typo must not crash module import, and must not silently flip a flag."""
    monkeypatch.setenv("GVF_SOME_FLAG", raw)
    assert get_env_bool("GVF_SOME_FLAG", default=True) is True
    assert get_env_bool("GVF_SOME_FLAG", default=False) is False


def test_get_env_bool_unset_uses_default(monkeypatch):
    monkeypatch.delenv("GVF_SOME_FLAG", raising=False)
    assert get_env_bool("GVF_SOME_FLAG", default=True) is True
    assert get_env_bool("GVF_SOME_FLAG", default=False) is False


@requires_guard
def test_offline_suite_runs_with_local_data_discovery_disabled():
    """The autouse conftest guard is in force for every test in this suite."""
    assert local_data_discovery_disabled() is True


def test_allow_local_data_discovery_fixture_lifts_the_guard(
    allow_local_data_discovery,
):
    assert local_data_discovery_disabled() is False


# ---------------------------------------------------------------------------
# utils.gene_metadata.default_variantfeatures_db_path
# ---------------------------------------------------------------------------


def _fake_sibling_variantfeatures_db(monkeypatch, tmp_path):
    """Anchor gene_metadata's repo-root guess at tmp_path and plant a DB there.

    ``default_variantfeatures_db_path`` derives its candidates from
    ``Path(__file__).resolve().parents[1]``, read from the module global at call
    time, so repointing ``__file__`` relocates both guesses under ``tmp_path``.
    """
    monkeypatch.setattr(
        gene_metadata, "__file__", str(tmp_path / "utils" / "gene_metadata.py")
    )
    db_path = tmp_path / "data" / "variants.db"
    db_path.parent.mkdir(parents=True, exist_ok=True)
    db_path.write_bytes(b"SQLite format 3\x00")
    return db_path


@requires_guard
def test_variantfeatures_db_guess_is_skipped_under_the_guard(monkeypatch, tmp_path):
    """A sibling database on disk must not change offline-suite behaviour."""
    _fake_sibling_variantfeatures_db(monkeypatch, tmp_path)

    assert gene_metadata.default_variantfeatures_db_path() is None


def test_variantfeatures_db_guess_is_used_when_discovery_is_allowed(
    monkeypatch, tmp_path, allow_local_data_discovery
):
    """Production behaviour is unchanged: the guess still resolves."""
    db_path = _fake_sibling_variantfeatures_db(monkeypatch, tmp_path)

    assert gene_metadata.default_variantfeatures_db_path() == db_path


@requires_guard
def test_explicit_variantfeatures_db_wins_over_the_guard(monkeypatch, tmp_path):
    """Tests that build their own fixture database keep working."""
    _fake_sibling_variantfeatures_db(monkeypatch, tmp_path)
    explicit = tmp_path / "explicit.db"
    explicit.write_bytes(b"SQLite format 3\x00")
    monkeypatch.setenv("VARIANTFEATURES_DB", str(explicit))

    assert local_data_discovery_disabled() is True
    assert gene_metadata.default_variantfeatures_db_path() == explicit


def test_legacy_variant_features_db_alias_also_wins(monkeypatch, tmp_path):
    explicit = tmp_path / "explicit.db"
    explicit.write_bytes(b"SQLite format 3\x00")
    monkeypatch.setenv("VARIANT_FEATURES_DB", str(explicit))

    assert gene_metadata.default_variantfeatures_db_path() == explicit


@requires_guard
def test_gene_metadata_falls_back_to_builtins_under_the_guard():
    """No database reachable means built-in metadata only — never a stale slice."""
    metadata = gene_metadata.get_gene_metadata("KCNH2")

    assert metadata.protein_length == 1159
    assert metadata.sources == ("builtin",)
    assert metadata.variantfeatures_db is None


# ---------------------------------------------------------------------------
# pipeline.steps._resolve_corpus_dir
# ---------------------------------------------------------------------------


@requires_guard
def test_corpus_dir_guess_is_skipped_under_the_guard(monkeypatch, tmp_path):
    monkeypatch.setattr(steps, "__file__", str(tmp_path / "pipeline" / "steps.py"))
    (tmp_path / "corpus").mkdir()

    assert steps._resolve_corpus_dir() is None


def test_corpus_dir_guess_is_used_when_discovery_is_allowed(
    monkeypatch, tmp_path, allow_local_data_discovery
):
    monkeypatch.setattr(steps, "__file__", str(tmp_path / "pipeline" / "steps.py"))
    corpus = tmp_path / "corpus"
    corpus.mkdir()

    assert steps._resolve_corpus_dir() == corpus


def test_explicit_corpus_dir_wins_over_the_guard(monkeypatch, tmp_path):
    corpus = tmp_path / "explicit_corpus"
    corpus.mkdir()
    monkeypatch.setenv("GVF_CORPUS_DIR", str(corpus))

    assert steps._resolve_corpus_dir() == corpus


# ---------------------------------------------------------------------------
# pipeline.paper_final_check._load_scout_zones
# ---------------------------------------------------------------------------


def _plant_corpus_scout_zones(monkeypatch, tmp_path, pmid="12345678", gene="KCNH2"):
    monkeypatch.setattr(paper_final_check, "REPO_ROOT", tmp_path)
    paper_dir = tmp_path / "corpus" / gene / pmid
    paper_dir.mkdir(parents=True)
    zones = paper_dir / f"{pmid}_DATA_ZONES.md"
    zones.write_text("zone text " * 100, encoding="utf-8")
    return zones


@requires_guard
def test_scout_zones_corpus_lookup_is_skipped_under_the_guard(monkeypatch, tmp_path):
    _plant_corpus_scout_zones(monkeypatch, tmp_path)

    assert paper_final_check._load_scout_zones("12345678", "KCNH2", None) == (
        None,
        None,
    )


def test_scout_zones_corpus_lookup_is_used_when_discovery_is_allowed(
    monkeypatch, tmp_path, allow_local_data_discovery
):
    zones = _plant_corpus_scout_zones(monkeypatch, tmp_path)

    text, path = paper_final_check._load_scout_zones("12345678", "KCNH2", None)

    assert path == str(zones)
    assert text is not None and text.startswith("zone text")


def test_scout_zones_still_read_from_an_explicit_run_dir_under_the_guard(
    monkeypatch, tmp_path
):
    """The guard only suppresses guessing; an explicit run dir is still read."""
    _plant_corpus_scout_zones(monkeypatch, tmp_path)
    run_path = tmp_path / "run"
    scout_dir = run_path / "scout_output"
    scout_dir.mkdir(parents=True)
    (scout_dir / "12345678_DATA_ZONES.md").write_text(
        "run zones " * 100, encoding="utf-8"
    )

    text, path = paper_final_check._load_scout_zones("12345678", "KCNH2", run_path)

    assert path == str(scout_dir / "12345678_DATA_ZONES.md")
    assert text is not None and text.startswith("run zones")


# ---------------------------------------------------------------------------
# cli.gvf_run._find_review_repo
# ---------------------------------------------------------------------------


def _plant_sibling_review_repo(monkeypatch, tmp_path):
    """A real-shaped sibling Variant_Browser checkout next to a fake repo root."""
    monkeypatch.setattr(gvf_run, "REPO_ROOT", tmp_path / "GeneVariantFetcher")
    repo = tmp_path / "Variant_Browser"
    (repo / "scripts").mkdir(parents=True)
    (repo / "scripts" / "gvf_publish.sh").write_text("#!/usr/bin/env bash\nexit 0\n")
    return repo


@requires_guard
def test_review_repo_sibling_guess_is_skipped_under_the_guard(monkeypatch, tmp_path):
    """Offline tests must not find the developer's real Variant_Browser checkout.

    `step_publish_review` swallows every exception from `subprocess.run`, so a
    test asserting "subprocess is never called" passed for the wrong reason
    whenever this guess resolved.
    """
    _plant_sibling_review_repo(monkeypatch, tmp_path)

    assert gvf_run._find_review_repo(None) is None


def test_review_repo_sibling_guess_is_used_when_discovery_is_allowed(
    monkeypatch, tmp_path, allow_local_data_discovery
):
    repo = _plant_sibling_review_repo(monkeypatch, tmp_path)

    assert gvf_run._find_review_repo(None) == repo


def test_explicit_review_repo_wins_over_the_guard(monkeypatch, tmp_path):
    _plant_sibling_review_repo(monkeypatch, tmp_path)
    explicit = tmp_path / "explicit_browser"
    (explicit / "scripts").mkdir(parents=True)
    (explicit / "scripts" / "gvf_publish.sh").write_text("#!/usr/bin/env bash\n")

    assert gvf_run._find_review_repo(explicit) == explicit
    monkeypatch.setenv("GVF_REVIEW_REPO", str(explicit))
    assert gvf_run._find_review_repo(None) == explicit


# ---------------------------------------------------------------------------
# benchmarks.cold_start_eval.run_cold_start._covered_genes
# ---------------------------------------------------------------------------


@requires_guard
def test_cold_start_coverage_ignores_the_local_corpus(monkeypatch, tmp_path):
    """Whether a gene counts as a cold start must not depend on local runs.

    `_covered_genes` unions in the top-level `corpus/` gene dirs. corpus/ is
    gitignored, so locally that set grows with every gene the developer has run
    — the day an LDLR run lands, `is_cold_start_gene("LDLR")` flips to False
    here and stays True in CI.
    """
    monkeypatch.setattr(run_cold_start, "BASE_DIR", tmp_path)
    (tmp_path / "corpus" / "LDLR").mkdir(parents=True)

    assert run_cold_start.is_cold_start_gene("LDLR")[0] is True
    assert "LDLR" not in run_cold_start._covered_genes()


def test_cold_start_coverage_reads_the_corpus_when_discovery_is_allowed(
    monkeypatch, tmp_path, allow_local_data_discovery
):
    """A real cold-start run still refuses a gene already cached in corpus/."""
    monkeypatch.setattr(run_cold_start, "BASE_DIR", tmp_path)
    (tmp_path / "corpus" / "LDLR").mkdir(parents=True)

    assert "LDLR" in run_cold_start._covered_genes()
    assert run_cold_start.is_cold_start_gene("LDLR")[0] is False


def test_cold_start_coverage_keeps_its_tracked_sources_under_the_guard():
    """The guard must not blind the harness to its committed sources."""
    covered = run_cold_start._covered_genes()

    assert {"KCNH2", "BRCA1", "MYBPC3"} <= covered


# ---------------------------------------------------------------------------
# The regression this whole file exists for
# ---------------------------------------------------------------------------


def test_default_db_path_is_checkout_independent(monkeypatch):
    """Same answer from the main checkout, a side worktree, and CI.

    The three anchors below are the real shapes: the main checkout (a sibling
    ``variantFeatures`` really is next to it), a `.claude/worktrees/<name>`
    worktree, and a CI runner.
    """
    for anchor in (
        "/Users/dev/GitRepos/GeneVariantFetcher/utils/gene_metadata.py",
        "/Users/dev/GitRepos/GeneVariantFetcher/.claude/worktrees/wt/utils/gene_metadata.py",
        "/home/runner/work/GeneVariantFetcher/GeneVariantFetcher/utils/gene_metadata.py",
    ):
        monkeypatch.setattr(gene_metadata, "__file__", anchor)
        assert gene_metadata.default_variantfeatures_db_path() is None
