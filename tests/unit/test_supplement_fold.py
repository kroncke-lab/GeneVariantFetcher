"""Unit tests for the on-disk supplement fold (C2)."""

from harvesting.supplement_fold import (
    FOLD_BEGIN,
    FOLD_END,
    _strip_folded_block,
    build_supplement_markdown,
    fold_supplements_into_full_context,
)
from scripts.fold_supplements import discover_pmids

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
    fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY)
    folded_twice = fc.read_text(encoding="utf-8")
    assert folded_twice == folded_once
    assert folded_twice.count("# MAIN TEXT") == 1
    assert backup.read_text(encoding="utf-8") == original


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
