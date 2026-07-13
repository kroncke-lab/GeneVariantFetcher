"""Unit tests for the on-disk supplement fold (C2)."""

import pytest

from harvesting.supplement_fold import (
    FOLD_BEGIN,
    FOLD_END,
    _strip_folded_block,
    build_supplement_markdown,
    fold_supplements_into_full_context,
)
from scripts.fold_supplements import discover_corpus_papers, discover_pmids, main

# A dummy converter is enough: .csv/.txt are read directly by _convert_supplement
# and never touch the converter object.
_DUMMY = object()


def _write_supp(supp_dir, name, text):
    supp_dir.mkdir(parents=True, exist_ok=True)
    (supp_dir / name).write_text(text, encoding="utf-8")


def test_strip_folded_block_roundtrip():
    base = "# MAIN\n\nbody\n"
    assert _strip_folded_block(base) == base  # nothing to strip
    folded = base.rstrip() + f"\n\n{FOLD_BEGIN}\n\nstuff\n\n{FOLD_END}\n"
    assert _strip_folded_block(folded).rstrip() == "# MAIN\n\nbody".rstrip()
    # Truncated end marker: drop from the begin marker onward.
    truncated = base.rstrip() + f"\n\n{FOLD_BEGIN}\n\norphan"
    assert FOLD_BEGIN not in _strip_folded_block(truncated)


def test_build_supplement_markdown_converts_csv_and_txt(tmp_path):
    supp = tmp_path / "12345678_supplements"
    _write_supp(supp, "tableS1.csv", "variant,carriers\nc.1A>G,3\n")
    _write_supp(supp, "notes.txt", "extra cohort detail\n")

    md, converted = build_supplement_markdown(supp, converter=_DUMMY)

    assert converted == 2
    assert "c.1A>G,3" in md
    assert "extra cohort detail" in md
    assert "SUPPLEMENTAL FILE" in md


def test_build_supplement_markdown_folds_nested_zip_extracted_files(tmp_path):
    # Files extracted from a .zip land in a subdirectory; the recursive walk must
    # fold them (the .zip itself stays excluded), and label them by relative path.
    supp = tmp_path / "12345678_supplements"
    _write_supp(supp, "top.csv", "variant,carriers\nc.1A>G,3\n")
    _write_supp(supp / "mmc1", "nested_table.csv", "variant,carriers\nc.9G>T,5\n")
    # cruft that must be ignored
    _write_supp(supp / "__MACOSX", "._junk.csv", "garbage\n")
    _write_supp(supp, ".DS_Store.csv", "garbage\n")
    _write_supp(supp / ".cache", "visible_name.csv", "garbage\n")

    md, converted = build_supplement_markdown(supp, converter=_DUMMY)

    assert converted == 2  # top + nested, NOT the cruft
    assert "c.1A>G,3" in md and "c.9G>T,5" in md
    assert "mmc1/nested_table.csv" in md  # nested provenance label
    assert "garbage" not in md


def test_build_supplement_markdown_empty_or_missing(tmp_path):
    assert build_supplement_markdown(tmp_path / "nope", converter=_DUMMY) == ("", 0)
    empty = tmp_path / "999_supplements"
    empty.mkdir()
    assert build_supplement_markdown(empty, converter=_DUMMY) == ("", 0)


def test_fold_is_idempotent_and_nondestructive(tmp_path):
    pmid = "12345678"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    original = "# MAIN TEXT\n\nKCNQ1 body text\n"
    fc.write_text(original, encoding="utf-8")
    _write_supp(
        harvest / f"{pmid}_supplements", "tableS1.csv", "variant,carriers\nc.1A>G,3\n"
    )

    # First fold.
    out = fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY)
    assert out == fc
    folded_once = fc.read_text(encoding="utf-8")
    assert "# MAIN TEXT" in folded_once
    assert "c.1A>G,3" in folded_once
    assert folded_once.count(FOLD_BEGIN) == 1
    assert folded_once.count(FOLD_END) == 1

    # Backup holds the true pre-fold original.
    backup = harvest / f"{pmid}_FULL_CONTEXT.md.pre_fold_bak"
    assert backup.read_text(encoding="utf-8") == original

    # Second fold is a no-op on content (no double-append) and keeps the backup.
    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) is None
    folded_twice = fc.read_text(encoding="utf-8")
    assert folded_twice == folded_once
    assert folded_twice.count("# MAIN TEXT") == 1
    assert backup.read_text(encoding="utf-8") == original


def test_fold_migrates_covered_legacy_inline_block_without_duplication(tmp_path):
    pmid = "44444444"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text(
        "# MAIN\n\nbody\n\n# SUPPLEMENTAL FILE 1: current.csv\n\n"
        "variant,carriers\nc.2A>G,4\n",
        encoding="utf-8",
    )
    _write_supp(
        harvest / f"{pmid}_supplements",
        "current.csv",
        "variant,carriers\nc.2A>G,4\n",
    )

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) == fc
    folded = fc.read_text(encoding="utf-8")
    assert folded.count("# MAIN") == 1
    assert folded.count("c.2A>G,4") == 1
    assert folded.count(FOLD_BEGIN) == 1

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) is None
    assert fc.read_text(encoding="utf-8") == folded


def test_fold_retains_legacy_content_missing_from_disk(tmp_path):
    pmid = "55555555"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text(
        "# MAIN\n\nbody\n\n# SUPPLEMENTAL FILE 1: missing.csv\n\nlegacy,only\n",
        encoding="utf-8",
    )
    _write_supp(
        harvest / f"{pmid}_supplements",
        "current.csv",
        "variant,carriers\nc.2A>G,4\n",
    )

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) == fc
    folded = fc.read_text(encoding="utf-8")
    assert "legacy,only" in folded
    assert "c.2A>G,4" in folded
    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) is None
    assert fc.read_text(encoding="utf-8") == folded


def test_first_fold_uses_good_files_when_one_conversion_fails(tmp_path, monkeypatch):
    pmid = "66666666"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    original = "# MAIN\n\nbody\n"
    fc.write_text(original, encoding="utf-8")
    supp_dir = harvest / f"{pmid}_supplements"
    _write_supp(supp_dir, "good.txt", "good supplement\n")
    _write_supp(supp_dir, "bad.txt", "bad supplement\n")

    from harvesting import supplement_fold

    real_convert = supplement_fold._convert_supplement

    def fail_one(**kwargs):
        if kwargs["file_path"].name == "bad.txt":
            raise ValueError("conversion failed")
        return real_convert(**kwargs)

    monkeypatch.setattr(supplement_fold, "_convert_supplement", fail_one)

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) == fc
    folded = fc.read_text(encoding="utf-8")
    assert "good supplement" in folded
    assert "bad supplement" not in folded
    assert (harvest / f"{pmid}_FULL_CONTEXT.md.pre_fold_bak").read_text() == original


def test_refold_conversion_failure_preserves_previous_block(tmp_path, monkeypatch):
    pmid = "67676767"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text("# MAIN\n\nbody\n", encoding="utf-8")
    supp_dir = harvest / f"{pmid}_supplements"
    _write_supp(supp_dir, "good.txt", "good supplement\n")
    _write_supp(supp_dir, "bad.txt", "bad supplement\n")

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) == fc
    folded = fc.read_text(encoding="utf-8")

    from harvesting import supplement_fold

    real_convert = supplement_fold._convert_supplement

    def fail_previously_folded_file(**kwargs):
        if kwargs["file_path"].name == "bad.txt":
            raise ValueError("conversion failed")
        return real_convert(**kwargs)

    monkeypatch.setattr(
        supplement_fold, "_convert_supplement", fail_previously_folded_file
    )

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) is None
    assert fc.read_text(encoding="utf-8") == folded


def test_fold_keeps_good_tables_when_another_file_converts_empty(tmp_path, monkeypatch):
    pmid = "77777777"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text("# MAIN\n\nbody\n", encoding="utf-8")
    supp_dir = harvest / f"{pmid}_supplements"
    _write_supp(supp_dir, "table.txt", "variant,carriers\nc.2A>G,4\n")
    _write_supp(supp_dir, "image_only.pdf", "not real PDF")

    from harvesting import supplement_fold

    real_convert = supplement_fold._convert_supplement

    def empty_one(**kwargs):
        if kwargs["file_path"].name == "image_only.pdf":
            return "", [], []
        return real_convert(**kwargs)

    monkeypatch.setattr(supplement_fold, "_convert_supplement", empty_one)

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) == fc
    assert "c.2A>G,4" in fc.read_text(encoding="utf-8")


def test_fold_returns_none_without_supplements(tmp_path):
    pmid = "22222222"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    original = "# MAIN\n\nbody\n"
    fc.write_text(original, encoding="utf-8")
    (harvest / f"{pmid}_supplements").mkdir()  # present but empty

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) is None
    assert fc.read_text(encoding="utf-8") == original  # untouched
    assert not (harvest / f"{pmid}_FULL_CONTEXT.md.pre_fold_bak").exists()  # no backup


def test_fold_returns_none_without_full_context(tmp_path):
    harvest = tmp_path / "pmc_fulltext"
    _write_supp(harvest / "33333333_supplements", "s.csv", "v,c\nc.1A>G,1\n")
    assert (
        fold_supplements_into_full_context("33333333", harvest, converter=_DUMMY)
        is None
    )


def test_discover_pmids_requires_both_artifacts(tmp_path):
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    # Has both -> discovered.
    (harvest / "111_FULL_CONTEXT.md").write_text("x", encoding="utf-8")
    _write_supp(harvest / "111_supplements", "a.csv", "v,c\nc.1A>G,1\n")
    # Supplements dir but no FULL_CONTEXT -> skipped.
    _write_supp(harvest / "222_supplements", "a.csv", "v,c\nc.2A>G,1\n")
    # FULL_CONTEXT but no supplements dir -> skipped.
    (harvest / "333_FULL_CONTEXT.md").write_text("x", encoding="utf-8")

    assert discover_pmids(harvest) == ["111"]


def test_discover_corpus_papers_scopes_genes(tmp_path):
    for gene, pmid in (("SCN5A", "111"), ("KCNH2", "222")):
        paper = tmp_path / "corpus" / gene / pmid
        paper.mkdir(parents=True)
        (paper / f"{pmid}_FULL_CONTEXT.md").write_text("body")
        (paper / f"{pmid}_supplements").mkdir()

    papers = discover_corpus_papers(tmp_path / "corpus", ["SCN5A"])

    assert papers == [("111", tmp_path / "corpus" / "SCN5A" / "111")]


def test_corpus_mode_rejects_missing_directory(tmp_path, monkeypatch, capsys):
    missing = tmp_path / "missing-corpus"
    monkeypatch.setattr("sys.argv", ["fold_supplements.py", "--corpus", str(missing)])

    with pytest.raises(SystemExit, match="2"):
        main()
    assert f"--corpus directory does not exist: {missing}" in capsys.readouterr().err
