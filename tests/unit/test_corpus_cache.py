"""Offline tests for the corpus read-cache (idempotent cross-run reuse).

`_consolidate_from_corpus` is the gate that makes new runs reuse already-fetched
source from corpus/<GENE>/<PMID>/ instead of re-downloading — but only when the
cached full text is usable, so a stub/compromised cached copy falls through to a
fresh fetch (e.g. after a new publisher key is added).
"""

from __future__ import annotations

from pathlib import Path

import pipeline.steps as steps


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

    recovered = steps._consolidate_from_corpus(
        ["111", "222", "333", "444"], harvest, "KCNH2", corpus
    )

    assert recovered == {"111"}, "only the usable cached paper should be reused"
    assert (harvest / "111_FULL_CONTEXT.md").exists()
    assert (harvest / "111_figures" / "fig_pmc_1.png").exists()
    assert (harvest / "111_supplements" / "table_s1.csv").exists()
    # stub / abstract-only / missing must NOT be copied (so they get re-fetched)
    assert not (harvest / "222_FULL_CONTEXT.md").exists()
    assert not (harvest / "333_FULL_CONTEXT.md").exists()
    assert not (harvest / "444_FULL_CONTEXT.md").exists()


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
