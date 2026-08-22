"""Offline tests for the corpus read-cache (idempotent cross-run reuse).

`_consolidate_from_corpus` is the gate that makes new runs reuse already-fetched
source from corpus/<GENE>/<PMID>/ instead of re-downloading — but only when the
cached full text is usable, so a stub/compromised cached copy falls through to a
fresh fetch (e.g. after a new publisher key is added).
"""

from __future__ import annotations

from pathlib import Path

import pipeline.steps as steps
from harvesting.supplement_fold import FOLD_BEGIN
from pipeline.source_quality import SUPPLEMENT_SURFACE_STATUS_MARKER


def _make_corpus_paper(
    corpus: Path,
    gene: str,
    pmid: str,
    *,
    body: str,
    figures: bool = False,
    supplements: bool = False,
) -> None:
    d = corpus / gene / pmid
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{pmid}_FULL_CONTEXT.md").write_text(body, encoding="utf-8")
    if figures:
        (d / f"{pmid}_figures").mkdir()
        (d / f"{pmid}_figures" / "fig_pmc_1.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    if supplements:
        (d / f"{pmid}_supplements").mkdir()
        (d / f"{pmid}_supplements" / "table_s1.csv").write_text("a,b\n1,2\n")


def test_corpus_cache_reuses_usable_skips_stub_and_missing(tmp_path: Path):
    corpus = tmp_path / "corpus"
    harvest = tmp_path / "run" / "pmc_fulltext"
    harvest.mkdir(parents=True)

    # Usable full text (well over MIN_EXTRACTION_INPUT_SIZE, no fallback marker)
    _make_corpus_paper(
        corpus,
        "KCNH2",
        "111",
        body="# Real paper\n" + "body text. " * 200,
        figures=True,
        supplements=True,
    )
    # Stub: tiny -> is_usable_fulltext_source == False
    _make_corpus_paper(corpus, "KCNH2", "222", body="too short")
    # Stub: abstract-only fallback marker -> not usable even though long
    _make_corpus_paper(
        corpus, "KCNH2", "333", body="# Abstract-only fallback\n" + "x " * 500
    )
    # A recovered main body with an uninspected supplement surface is useful for
    # the current extraction but must be re-fetched on a future run.
    _make_corpus_paper(
        corpus,
        "KCNH2",
        "555",
        body=(
            f"<!-- {SUPPLEMENT_SURFACE_STATUS_MARKER} unavailable -->\n\n"
            + "body text. " * 200
        ),
    )

    recovered = steps._consolidate_from_corpus(
        ["111", "222", "333", "444", "555"], harvest, "KCNH2", corpus
    )

    assert recovered == {"111"}, "only the usable cached paper should be reused"
    assert (harvest / "111_FULL_CONTEXT.md").exists()
    assert (harvest / "111_figures" / "fig_pmc_1.png").exists()
    assert (harvest / "111_supplements" / "table_s1.csv").exists()
    reused_text = (harvest / "111_FULL_CONTEXT.md").read_text(encoding="utf-8")
    assert FOLD_BEGIN in reused_text
    assert "a,b\n1,2" in reused_text
    # stub / abstract-only / missing must NOT be copied (so they get re-fetched)
    assert not (harvest / "222_FULL_CONTEXT.md").exists()
    assert not (harvest / "333_FULL_CONTEXT.md").exists()
    assert not (harvest / "444_FULL_CONTEXT.md").exists()
    assert not (harvest / "555_FULL_CONTEXT.md").exists()


def test_corpus_cache_stages_cleaned_sibling_of_empty_full_context(tmp_path: Path):
    """KCNQ1 27114410: 0-byte cached FULL_CONTEXT beside a 17.7 KB CLEANED.

    The stub is still not reused (the harvester must re-attempt the fetch),
    but the populated CLEANED sibling is staged so extraction keeps the paper
    when that re-fetch fails again.
    """
    corpus = tmp_path / "corpus"
    harvest = tmp_path / "run" / "pmc_fulltext"
    harvest.mkdir(parents=True)

    _make_corpus_paper(corpus, "KCNQ1", "27114410", body="")
    cleaned_body = "# MAIN TEXT\n" + "body text. " * 200
    (corpus / "KCNQ1" / "27114410" / "27114410_CLEANED.md").write_text(
        cleaned_body, encoding="utf-8"
    )

    recovered = steps._consolidate_from_corpus(["27114410"], harvest, "KCNQ1", corpus)

    assert recovered == set(), "the empty stub itself must not count as reused"
    assert not (harvest / "27114410_FULL_CONTEXT.md").exists()
    staged = harvest / "27114410_CLEANED.md"
    assert staged.read_text(encoding="utf-8") == cleaned_body
    # The corpus copy is untouched.
    corpus_fc = corpus / "KCNQ1" / "27114410" / "27114410_FULL_CONTEXT.md"
    assert corpus_fc.stat().st_size == 0


def test_prior_run_cache_retries_body_only_source(tmp_path: Path):
    output_base = tmp_path / "results"
    prior = output_base / "KCNH2" / "old" / "pmc_fulltext"
    prior.mkdir(parents=True)
    (prior / "555_FULL_CONTEXT.md").write_text(
        f"<!-- {SUPPLEMENT_SURFACE_STATUS_MARKER} scrape_error -->\n\n"
        + "body text. " * 1000,
        encoding="utf-8",
    )
    harvest = output_base / "KCNH2" / "new" / "pmc_fulltext"
    harvest.mkdir(parents=True)

    recovered = steps._consolidate_prior_downloads(["555"], harvest, output_base)

    assert recovered == set()
    assert not (harvest / "555_FULL_CONTEXT.md").exists()


def test_corpus_cache_no_corpus_dir_is_noop(tmp_path: Path):
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    assert steps._consolidate_from_corpus(["111"], harvest, "KCNH2", None) == set()
    # missing gene subdir -> empty
    (tmp_path / "corpus").mkdir()
    assert (
        steps._consolidate_from_corpus(["111"], harvest, "KCNH2", tmp_path / "corpus")
        == set()
    )


def test_resolve_corpus_dir_honors_env(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("GVF_CORPUS_DIR", str(tmp_path))
    assert steps._resolve_corpus_dir() == tmp_path
    monkeypatch.setenv("GVF_CORPUS_DIR", str(tmp_path / "does_not_exist"))
    assert steps._resolve_corpus_dir() is None
