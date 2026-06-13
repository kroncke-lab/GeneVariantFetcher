"""Unit test for the supplement-fold-gap QC (flat harvest-dir mode)."""

from __future__ import annotations

from pathlib import Path

from scripts.recall_audit.supplement_fold_gap import (
    _full_context_reflects_supplements,
    _has_convertible_supplements,
)


def _supp(harvest: Path, pmid: str, name: str, text: str = "v,c\nc.1A>G,1\n") -> None:
    d = harvest / f"{pmid}_supplements"
    d.mkdir(parents=True, exist_ok=True)
    (d / name).write_text(text, encoding="utf-8")


def test_convertible_count_skips_zip_and_cruft(tmp_path):
    _supp(tmp_path, "111", "tableS1.csv")
    _supp(tmp_path, "111", "raw.zip", "binary")
    (tmp_path / "111_supplements" / "__MACOSX").mkdir()
    (tmp_path / "111_supplements" / "__MACOSX" / "._x.csv").write_text("junk")
    # nested (extracted-zip) convertible file counts
    (tmp_path / "111_supplements" / "sub").mkdir()
    (tmp_path / "111_supplements" / "sub" / "nested.xlsx").write_text("x")

    assert (
        _has_convertible_supplements(tmp_path / "111_supplements") == 2
    )  # csv + nested xlsx
    assert _has_convertible_supplements(tmp_path / "999_supplements") == 0


def test_full_context_reflects_supplements_marker(tmp_path):
    folded = tmp_path / "folded_FULL_CONTEXT.md"
    folded.write_text(
        "# MAIN\n\n# SUPPLEMENTAL FILE 1: tableS1.csv\n\nv,c\n", encoding="utf-8"
    )
    # body merely *mentions* "supplementary" but never folded the file content
    unparsed = tmp_path / "unparsed_FULL_CONTEXT.md"
    unparsed.write_text(
        "# MAIN\n\nSee Supplementary Table 1 (not included).\n", encoding="utf-8"
    )

    assert _full_context_reflects_supplements(folded) is True
    assert _full_context_reflects_supplements(unparsed) is False
