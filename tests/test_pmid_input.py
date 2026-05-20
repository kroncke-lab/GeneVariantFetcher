"""Tests for explicit PMID input parsing in the extract CLI.

Covers:
- `_parse_pmid_arg` (comma-separated and `@file` syntax)
- `_read_pmid_file` (one-per-line files with `#` comments)
- The shipped KCNH2 gold PMID fixture
"""

from __future__ import annotations

from pathlib import Path

import pytest
import typer

from cli import _parse_pmid_arg, _read_pmid_file


class TestParsePmidArg:
    def test_comma_separated_basic(self):
        assert _parse_pmid_arg("123,456,789") == ["123", "456", "789"]

    def test_comma_separated_dedupes_preserving_order(self):
        assert _parse_pmid_arg("789,123,456,123,789") == ["789", "123", "456"]

    def test_drops_non_digit_tokens_silently(self):
        # Non-digit tokens (whitespace, prose, accidental headers) drop without error.
        assert _parse_pmid_arg("123,abc,456,,789") == ["123", "456", "789"]

    def test_strips_whitespace_around_tokens(self):
        assert _parse_pmid_arg(" 123 , 456 ,  789 ") == ["123", "456", "789"]

    def test_at_file_syntax(self, tmp_path: Path):
        f = tmp_path / "list.txt"
        f.write_text("# header\n\n123\n456\n# mid comment\n789\n", encoding="utf-8")
        assert _parse_pmid_arg(f"@{f}") == ["123", "456", "789"]

    def test_at_file_missing_raises_bad_parameter(self, tmp_path: Path):
        with pytest.raises(typer.BadParameter):
            _parse_pmid_arg(f"@{tmp_path / 'nope.txt'}")


class TestReadPmidFile:
    def test_skips_blank_lines_and_comments(self, tmp_path: Path):
        f = tmp_path / "list.txt"
        f.write_text(
            "# top comment\n\n12345\n   \n# another comment\n67890\n",
            encoding="utf-8",
        )
        assert _read_pmid_file(f) == ["12345", "67890"]

    def test_strips_trailing_inline_comment(self, tmp_path: Path):
        # Allow `12345 # paper title` style annotations on PMID lines.
        f = tmp_path / "list.txt"
        f.write_text("12345 # KCNH2 cohort paper\n67890#tagged\n", encoding="utf-8")
        assert _read_pmid_file(f) == ["12345", "67890"]

    def test_missing_file_raises_bad_parameter(self, tmp_path: Path):
        with pytest.raises(typer.BadParameter):
            _read_pmid_file(tmp_path / "missing.txt")

    def test_passes_through_non_digit_tokens(self, tmp_path: Path):
        # _read_pmid_file does not validate digit-only — that's `_parse_pmid_arg`'s job.
        # Keep the contract: this function returns raw tokens, validation happens upstream.
        f = tmp_path / "list.txt"
        f.write_text("12345\nbadtoken\n67890\n", encoding="utf-8")
        assert _read_pmid_file(f) == ["12345", "badtoken", "67890"]


class TestGoldStandardFixture:
    """The shipped KCNH2 PMID fixture stays valid for explicit-PMID smoke tests."""

    GOLD_PATH = Path("tests/fixtures/pmids/kcnh2_gold_pmids.txt")

    def test_file_exists(self):
        assert self.GOLD_PATH.exists(), (
            "tests/fixtures/pmids/kcnh2_gold_pmids.txt is shipped in-repo "
            "as a reusable explicit-PMID input fixture."
        )

    def test_parses_to_nonempty_pmid_list(self):
        pmids = _parse_pmid_arg(f"@{self.GOLD_PATH}")
        assert len(pmids) > 0
        assert all(p.isdigit() for p in pmids)
        assert len(set(pmids)) == len(pmids), "no duplicates allowed"

    def test_count_matches_expected_baseline(self):
        # 262 distinct PMIDs derived from the curated KCNH2 Excel database
        # (rows where excel_variant_raw is non-empty in discrepancies.csv).
        # If this number changes, regenerate the fixture from the normalized
        # KCNH2 gold input or scored discrepancy artifacts and update this test.
        pmids = _parse_pmid_arg(f"@{self.GOLD_PATH}")
        assert len(pmids) == 262
