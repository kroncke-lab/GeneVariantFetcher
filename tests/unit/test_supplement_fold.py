"""Unit tests for the on-disk supplement fold (C2)."""

import json
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


def test_converter_error_placeholder_is_not_counted_as_source(tmp_path):
    class InvalidPdfConverter:
        @staticmethod
        def pdf_to_markdown(_path):
            return "[Invalid PDF file: supplement.pdf]\n\n"

    supp = tmp_path / "12345678_supplements"
    supp.mkdir()
    (supp / "supplement.pdf").write_bytes(b"<html>sign-in page</html>")

    markdown, converted = build_supplement_markdown(
        supp, converter=InvalidPdfConverter()
    )

    assert markdown == ""
    assert converted == 0


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


def test_fold_does_not_duplicate_members_already_expanded_from_archive(tmp_path):
    pmid = "45454545"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    original = (
        "# MAIN\n\nbody\n\n"
        "# SUPPLEMENTAL FILE 1: archive.zip\n\n"
        "##### Nested file: table_s1.csv\n\n"
        "variant,carriers\nc.2A>G,4\n"
    )
    fc.write_text(original, encoding="utf-8")
    _write_supp(
        harvest / f"{pmid}_supplements" / "archive",
        "table_s1.csv",
        "variant,carriers\nc.2A>G,4\n",
    )

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) is None
    assert fc.read_text(encoding="utf-8") == original
    assert fc.read_text(encoding="utf-8").count("c.2A>G,4") == 1


def test_refold_removes_redundant_member_block_after_archive_expansion(tmp_path):
    pmid = "46464646"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    base = (
        "# MAIN\n\nbody\n\n"
        "# SUPPLEMENTAL FILE 1: archive.zip\n\n"
        "##### Nested file: table_s1.csv\n\n"
        "variant,carriers\nc.2A>G,4\n"
    )
    fc.write_text(
        base.rstrip()
        + f"\n\n{FOLD_BEGIN}\n\n# FOLDED SUPPLEMENTS (re-extraction aid)\n"
        "# SUPPLEMENTAL FILE 1: archive/table_s1.csv\n\n"
        "variant,carriers\nc.2A>G,4\n\n"
        f"{FOLD_END}\n",
        encoding="utf-8",
    )
    _write_supp(
        harvest / f"{pmid}_supplements" / "archive",
        "table_s1.csv",
        "variant,carriers\nc.2A>G,4\n",
    )

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) == fc
    cleaned = fc.read_text(encoding="utf-8")
    assert cleaned == base
    assert cleaned.count("c.2A>G,4") == 1
    assert FOLD_BEGIN not in cleaned


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


def test_refold_drops_old_placeholder_and_adds_new_valid_supplement(tmp_path):
    class InvalidPdfConverter:
        @staticmethod
        def pdf_to_markdown(_path):
            return "[Invalid PDF file: bad.pdf]\n\n"

    pmid = "68686868"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text(
        "# MAIN\n\nbody\n\n"
        f"{FOLD_BEGIN}\n\n# FOLDED SUPPLEMENTS (re-extraction aid)\n"
        "# SUPPLEMENTAL FILE 1: bad.pdf\n\n"
        "[Invalid PDF file: bad.pdf]\n\n"
        f"{FOLD_END}\n",
        encoding="utf-8",
    )
    supp_dir = harvest / f"{pmid}_supplements"
    supp_dir.mkdir()
    (supp_dir / "bad.pdf").write_bytes(b"not a real PDF")
    _write_supp(supp_dir, "table_s1.csv", "variant,carriers\nc.2A>G,4\n")

    assert (
        fold_supplements_into_full_context(
            pmid, harvest, converter=InvalidPdfConverter()
        )
        == fc
    )
    refolded = fc.read_text(encoding="utf-8")
    assert "c.2A>G,4" in refolded
    assert "[Invalid PDF file" not in refolded
    assert refolded.count(FOLD_BEGIN) == 1


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


def test_discover_corpus_papers_ignores_hidden_genes_and_papers(tmp_path):
    corpus = tmp_path / "corpus"
    visible = corpus / "SCN5A" / "111"
    visible.mkdir(parents=True)
    (visible / "111_FULL_CONTEXT.md").write_text("body")
    (visible / "111_supplements").mkdir()
    hidden_gene_paper = corpus / ".cache" / "222"
    hidden_gene_paper.mkdir(parents=True)
    (hidden_gene_paper / "222_FULL_CONTEXT.md").write_text("body")
    (hidden_gene_paper / "222_supplements").mkdir()
    hidden_paper = corpus / "SCN5A" / ".cache"
    hidden_paper.mkdir()
    (hidden_paper / ".cache_FULL_CONTEXT.md").write_text("body")
    (hidden_paper / ".cache_supplements").mkdir()

    assert discover_corpus_papers(corpus, []) == [("111", visible)]


def test_corpus_mode_rejects_missing_directory(tmp_path, monkeypatch, capsys):
    missing = tmp_path / "missing-corpus"
    monkeypatch.setattr("sys.argv", ["fold_supplements.py", "--corpus", str(missing)])

    with pytest.raises(SystemExit, match="2"):
        main()
    assert f"--corpus directory does not exist: {missing}" in capsys.readouterr().err


def test_cli_rejects_missing_pmids_file(tmp_path, monkeypatch, capsys):
    corpus = tmp_path / "corpus"
    corpus.mkdir()
    missing = tmp_path / "missing-pmids.txt"
    monkeypatch.setattr(
        "sys.argv",
        [
            "fold_supplements.py",
            "--corpus",
            str(corpus),
            "--pmids-file",
            str(missing),
        ],
    )

    with pytest.raises(SystemExit, match="2"):
        main()
    assert f"--pmids-file does not exist: {missing}" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# PDF identity gate at fold time (KCNQ1 31520628: cited CDC vital-statistics
# PDFs were captured as supplements and folded into FULL_CONTEXT)
# ---------------------------------------------------------------------------

#: A cited foreign report prints its OWN doi, which contradicts the parent
#: article's. That contradiction is what condemns it — a numbered-serial
#: masthead alone only annotates, because a paper's own article PDF opens the
#: same way (LMNA 17334235).
_UNRELATED_PDF_TEXT = (
    "National Vital Statistics Reports Volume 60, Number 7. "
    "doi:10.15620/cdc.nvsr6007. "
    "Deaths: preliminary statistics for 2011. Division of Vital Statistics."
)


def _convert_pdfs_as(text):
    """A _convert_supplement stand-in that renders every .pdf as ``text``."""
    from harvesting import supplement_fold

    real_convert = supplement_fold._convert_supplement

    def fake(**kwargs):
        if kwargs["file_path"].suffix == ".pdf":
            return text, 0, []
        return real_convert(**kwargs)

    return fake


def test_fold_excludes_identity_unverified_pdf(tmp_path, monkeypatch, caplog):
    import logging

    from harvesting import supplement_fold

    pmid = "31520628"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text("# MAIN\n\nbody\n", encoding="utf-8")
    supp_dir = harvest / f"{pmid}_supplements"
    _write_supp(supp_dir, "nvsr60_07.pdf", "raw pdf bytes stand-in")
    _write_supp(supp_dir, "tableS2.csv", "variant,carriers\nc.1A>G,3\n")
    # The gate needs the parent DOI to recognize the report's own DOI as a
    # contradiction; without one, a foreign document is merely unidentified.
    (harvest / f"{pmid}_artifacts.json").write_text(
        json.dumps({"doi": "10.1016/j.ajog.2019.09.004"}), encoding="utf-8"
    )

    monkeypatch.setattr(
        supplement_fold, "_convert_supplement", _convert_pdfs_as(_UNRELATED_PDF_TEXT)
    )

    with caplog.at_level(logging.WARNING, logger="harvesting.supplement_fold"):
        out = fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY)

    assert out == fc
    folded = fc.read_text(encoding="utf-8")
    assert "c.1A>G" in folded  # the real supplement still folds
    assert "National Vital Statistics" not in folded
    assert "nvsr60_07.pdf" not in folded
    assert any(
        "supplement_identity_unverified" in record.getMessage()
        for record in caplog.records
    )
    # Never delete downloaded source.
    assert (supp_dir / "nvsr60_07.pdf").exists()


def test_fold_keeps_pdf_verified_by_manifest_doi(tmp_path, monkeypatch):
    import json

    from harvesting import supplement_fold

    pmid = "31520628"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text("# MAIN\n\nbody\n", encoding="utf-8")
    supp_dir = harvest / f"{pmid}_supplements"
    _write_supp(supp_dir, "cohort_tables.pdf", "raw pdf bytes stand-in")
    (harvest / f"{pmid}_artifacts.json").write_text(
        json.dumps({"pmid": pmid, "doi": "10.1016/j.ajog.2019.09.004"}),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        supplement_fold,
        "_convert_supplement",
        _convert_pdfs_as("Cohort tables. doi:10.1016/j.ajog.2019.09.004"),
    )

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) == fc
    folded = fc.read_text(encoding="utf-8")
    assert "cohort_tables.pdf" in folded
    assert "10.1016/j.ajog.2019.09.004" in folded


def test_refold_strips_previously_folded_misbound_pdf(tmp_path, monkeypatch):
    """A stale fold whose only content is now-quarantined must be rebuilt.

    31520628's supplements dir holds nothing but the two mis-bound CDC PDFs,
    so the refold produces no markdown — the old block must still go.
    """
    from harvesting import supplement_fold

    pmid = "31520628"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    stale = (
        "# MAIN\n\nbody\n"
        f"\n{FOLD_BEGIN}\n\n"
        "# FOLDED SUPPLEMENTS (re-extraction aid)\n"
        "\n\n# SUPPLEMENTAL FILE 1: nvsr60_07.pdf\n\n"
        f"{_UNRELATED_PDF_TEXT}\n"
        f"\n{FOLD_END}\n"
    )
    fc.write_text(stale, encoding="utf-8")
    supp_dir = harvest / f"{pmid}_supplements"
    _write_supp(supp_dir, "nvsr60_07.pdf", "raw pdf bytes stand-in")
    (harvest / f"{pmid}_artifacts.json").write_text(
        json.dumps({"doi": "10.1016/j.ajog.2019.09.004"}), encoding="utf-8"
    )

    monkeypatch.setattr(
        supplement_fold, "_convert_supplement", _convert_pdfs_as(_UNRELATED_PDF_TEXT)
    )

    assert fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY) == fc
    rebuilt = fc.read_text(encoding="utf-8")
    assert FOLD_BEGIN not in rebuilt
    assert "National Vital Statistics" not in rebuilt
    assert "# MAIN" in rebuilt
    # Non-destructive: the pre-fold backup preserves what was overwritten.
    assert (harvest / f"{pmid}_FULL_CONTEXT.md.pre_fold_bak").exists()


def test_fold_includes_unverified_pdf_and_says_so(tmp_path, monkeypatch, caplog):
    """Absence of a marker is not evidence of a foreign document.

    A scanned supplement named ``download.pdf`` has no filename marker, no
    DOI, and no 'supplement' wording — the fail-closed gate dropped exactly
    this shape. It must fold, with the flag recorded in the log.
    """
    import logging

    from harvesting import supplement_fold

    pmid = "27114410"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text("# MAIN\n\nbody\n", encoding="utf-8")
    supp_dir = harvest / f"{pmid}_supplements"
    _write_supp(supp_dir, "download.pdf", "raw pdf bytes stand-in")

    monkeypatch.setattr(
        supplement_fold,
        "_convert_supplement",
        _convert_pdfs_as("Table 1\n\nKCNQ1 c.1A>G 3 carriers\n"),
    )

    with caplog.at_level(logging.WARNING, logger="harvesting.supplement_fold"):
        out = fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY)

    assert out == fc
    folded = fc.read_text(encoding="utf-8")
    assert "KCNQ1 c.1A>G 3 carriers" in folded
    assert "download.pdf" in folded
    assert any(
        "folding download.pdf unverified" in record.getMessage()
        for record in caplog.records
    )


def test_fold_excludes_pdf_harvested_from_reference_markup(tmp_path, monkeypatch):
    """Provenance recorded at scrape time is read back from the artifacts file."""
    import json

    from harvesting import supplement_fold

    pmid = "31520628"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text("# MAIN\n\nbody\n", encoding="utf-8")
    supp_dir = harvest / f"{pmid}_supplements"
    _write_supp(supp_dir, "report_2011.pdf", "raw pdf bytes stand-in")
    (harvest / f"{pmid}_artifacts.json").write_text(
        json.dumps(
            {
                "pmid": pmid,
                "doi": "10.1016/j.ajog.2019.09.004",
                "supplements": [
                    {
                        "filename": "report_2011.pdf",
                        "url": "https://stats.example.gov/report_2011.pdf",
                        "source": "scraper_reference_section",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        supplement_fold,
        "_convert_supplement",
        _convert_pdfs_as("Deaths: preliminary statistics. " + "text " * 200),
    )

    out = fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY)

    assert out is None  # nothing foldable, nothing to rewrite
    assert "report_2011.pdf" not in fc.read_text(encoding="utf-8")
    assert (supp_dir / "report_2011.pdf").exists()


# ---------------------------------------------------------------------------
# Fold size caps and the binary-garbage backstop (RYR2 PMID 32508047: a
# 13.6 MB .doc whose conversion leaked 14.5 MB of raw OLE bytes into
# FULL_CONTEXT and wedged the data scout for 2h45m+).
# ---------------------------------------------------------------------------


def test_fold_skips_binary_garbage_conversion(tmp_path, monkeypatch, caplog):
    from harvesting import supplement_fold

    pmid = "32508047"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text("# MAIN TEXT\n\nRYR2 body text\n", encoding="utf-8")
    supp = harvest / f"{pmid}_supplements"
    _write_supp(supp, "good.csv", "variant,carriers\nc.1A>G,3\n")
    _write_supp(supp, "leaky.txt", "placeholder")

    soup = ("\x00\x01g\x00\x00\x0b|ˇˇ\x02\x03\x04 " * 200)[:4000]
    real_convert = supplement_fold._convert_supplement

    def leak_one(*, file_path, **kwargs):
        if file_path.name == "leaky.txt":
            return soup, 0, []
        return real_convert(file_path=file_path, **kwargs)

    monkeypatch.setattr(supplement_fold, "_convert_supplement", leak_one)

    with caplog.at_level("WARNING"):
        out = fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY)

    assert out == fc
    folded = fc.read_text(encoding="utf-8")
    assert "c.1A>G,3" in folded
    assert "leaky.txt" not in folded
    assert "\x00" not in folded
    assert any("raw binary" in r.getMessage() for r in caplog.records)
    assert (supp / "leaky.txt").exists()  # source file untouched


def test_fold_skips_file_over_per_file_cap_and_can_shrink_old_fold(
    tmp_path, monkeypatch, caplog
):
    pmid = "11112222"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    # Pre-existing fold block already carries the oversized file's label, the
    # way 32508047's corpus FULL_CONTEXT carried the 14.5 MB soup block.
    fc.write_text(
        "# MAIN TEXT\n\nbody\n\n"
        f"{FOLD_BEGIN}\n\n"
        "# FOLDED SUPPLEMENTS (re-extraction aid)\n"
        "\n\n# SUPPLEMENTAL FILE 1: huge.txt\n\nOLD GIANT CONTENT\n"
        "\n\n# SUPPLEMENTAL FILE 2: small.txt\n\nsmall detail\n"
        f"\n\n{FOLD_END}\n",
        encoding="utf-8",
    )
    supp = harvest / f"{pmid}_supplements"
    _write_supp(supp, "huge.txt", "x" * 500)
    _write_supp(supp, "small.txt", "small detail\n")
    monkeypatch.setenv("GVF_FOLD_MAX_FILE_CHARS", "100")

    with caplog.at_level("WARNING"):
        out = fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY)

    assert out == fc
    folded = fc.read_text(encoding="utf-8")
    # The oversized file is out of the fold (not a refold refusal), the
    # small one is still there, and the on-disk file survives.
    assert "OLD GIANT CONTENT" not in folded
    assert "huge.txt" not in folded
    assert "small detail" in folded
    assert any(
        "oversized converted supplement" in r.getMessage() for r in caplog.records
    )
    assert (supp / "huge.txt").exists()


def test_fold_total_cap_skips_file_that_would_overflow(tmp_path, monkeypatch, caplog):
    pmid = "33334444"
    harvest = tmp_path / "pmc_fulltext"
    harvest.mkdir()
    fc = harvest / f"{pmid}_FULL_CONTEXT.md"
    fc.write_text("# MAIN TEXT\n\nbody\n", encoding="utf-8")
    supp = harvest / f"{pmid}_supplements"
    # Alphabetical fold order: a.txt (60 chars) fits, b.txt (60 chars) would
    # overflow the 100-char total budget and is skipped, c.txt (20) still fits.
    _write_supp(supp, "a.txt", "A" * 60)
    _write_supp(supp, "b.txt", "B" * 60)
    _write_supp(supp, "c.txt", "C" * 20)
    monkeypatch.setenv("GVF_FOLD_MAX_TOTAL_CHARS", "100")

    with caplog.at_level("WARNING"):
        out = fold_supplements_into_full_context(pmid, harvest, converter=_DUMMY)

    assert out == fc
    folded = fc.read_text(encoding="utf-8")
    assert "A" * 60 in folded
    assert "B" * 60 not in folded
    assert "C" * 20 in folded
    assert any("fold size budget" in r.getMessage() for r in caplog.records)
