"""Unit test for the corpus->flat-harvest bridge."""

from __future__ import annotations

from pathlib import Path

from scripts.corpus_to_harvest import build


def test_build_copies_fulltext_and_symlinks_supplements(tmp_path):
    corpus = tmp_path / "corpus"
    pdir = corpus / "SCN5A" / "29325976"
    (pdir / "29325976_supplements").mkdir(parents=True)
    (pdir / "29325976_FULL_CONTEXT.md").write_text("# body\n", encoding="utf-8")
    (pdir / "29325976_supplements" / "mmc1.docx").write_text("x", encoding="utf-8")
    # condensed forms that must NOT be copied (so folded FULL_CONTEXT is the source)
    (pdir / "29325976_DATA_ZONES.md").write_text("zones", encoding="utf-8")
    # a second paper to confirm --pmids filtering
    p2 = corpus / "SCN5A" / "11111111"
    p2.mkdir(parents=True)
    (p2 / "11111111_FULL_CONTEXT.md").write_text("# other\n", encoding="utf-8")

    out = tmp_path / "flat"
    n_papers, n_supp = build("SCN5A", corpus, out, {"29325976"})

    assert n_papers == 1 and n_supp == 1  # only the requested PMID
    fc = out / "29325976_FULL_CONTEXT.md"
    assert fc.is_file() and not fc.is_symlink()  # copied (fold won't mutate corpus)
    supp = out / "29325976_supplements"
    assert supp.is_symlink()  # symlinked (no duplication)
    assert (supp / "mmc1.docx").read_text() == "x"
    assert not (out / "29325976_DATA_ZONES.md").exists()  # condensed not bridged
    assert not (out / "11111111_FULL_CONTEXT.md").exists()  # filtered out
