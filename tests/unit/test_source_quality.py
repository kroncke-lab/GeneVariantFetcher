import logging

from pipeline.source_quality import (
    SUPPLEMENT_SURFACE_STATUS_MARKER,
    demote_empty_full_context,
    is_reusable_fulltext_source,
    is_usable_fulltext_source,
)

_REAL_BODY = "# MAIN TEXT\n\n" + ("methods results discussion variants. " * 80)


def test_source_quality_rejects_short_placeholder(tmp_path):
    source = tmp_path / "123_FULL_CONTEXT.md"
    source.write_text("not enough text", encoding="utf-8")

    assert is_usable_fulltext_source(source) is False


def test_source_quality_rejects_abstract_only_fallback(tmp_path):
    source = tmp_path / "123_FULL_CONTEXT.md"
    source.write_text(
        "# ABSTRACT-ONLY FALLBACK\n\n"
        "> **WARNING:** Full text could not be retrieved for PMID 123.\n"
        "> This document contains only the PubMed abstract and metadata.\n"
        + ("abstract text " * 80),
        encoding="utf-8",
    )

    assert is_usable_fulltext_source(source) is False


def test_source_quality_accepts_long_non_fallback_markdown(tmp_path):
    source = tmp_path / "123_FULL_CONTEXT.md"
    source.write_text(
        "# MAIN TEXT\n\nKCNH2 p.Arg1His cohort table.\n"
        + ("methods results discussion variants. " * 80),
        encoding="utf-8",
    )

    assert is_usable_fulltext_source(source) is True
    assert is_reusable_fulltext_source(source) is True


def test_source_quality_rejects_large_abstract_reference_shell(tmp_path):
    source = tmp_path / "123_FULL_CONTEXT.md"
    source.write_text(
        "## Abstract\n"
        + ("cohort abstract text. " * 200)
        + "\n## References\n"
        + ("citation. " * 300),
        encoding="utf-8",
    )

    assert is_usable_fulltext_source(source) is False


def test_body_only_source_is_extractable_but_not_reusable(tmp_path):
    source = tmp_path / "123_FULL_CONTEXT.md"
    source.write_text(
        f"<!-- {SUPPLEMENT_SURFACE_STATUS_MARKER} unavailable -->\n\n"
        "# MAIN TEXT\n\n" + ("methods results discussion variants. " * 80),
        encoding="utf-8",
    )

    assert is_usable_fulltext_source(source) is True
    assert is_reusable_fulltext_source(source) is False


def test_demote_empty_full_context_returns_populated_cleaned(tmp_path, caplog):
    """The KCNQ1 27114410 shape: 0-byte FULL_CONTEXT beside a real CLEANED."""
    fc = tmp_path / "27114410_FULL_CONTEXT.md"
    fc.write_text("", encoding="utf-8")
    cleaned = tmp_path / "27114410_CLEANED.md"
    cleaned.write_text(_REAL_BODY, encoding="utf-8")

    with caplog.at_level(logging.WARNING, logger="pipeline.source_quality"):
        selected = demote_empty_full_context(fc)

    assert selected == cleaned
    # The warning names both files.
    joined = " ".join(record.getMessage() for record in caplog.records)
    assert "27114410_FULL_CONTEXT.md" in joined
    assert "27114410_CLEANED.md" in joined
    # Neither file was modified.
    assert fc.stat().st_size == 0
    assert cleaned.read_text(encoding="utf-8") == _REAL_BODY


def test_demote_empty_full_context_keeps_populated_full_context(tmp_path):
    fc = tmp_path / "111_FULL_CONTEXT.md"
    fc.write_text(_REAL_BODY, encoding="utf-8")
    (tmp_path / "111_CLEANED.md").write_text(_REAL_BODY, encoding="utf-8")

    assert demote_empty_full_context(fc) == fc


def test_demote_empty_full_context_without_usable_sibling_is_noop(tmp_path):
    fc = tmp_path / "222_FULL_CONTEXT.md"
    fc.write_text("", encoding="utf-8")
    # A sibling below the usability bar must not win either.
    (tmp_path / "222_CLEANED.md").write_text("tiny", encoding="utf-8")

    assert demote_empty_full_context(fc) == fc


def test_demote_empty_full_context_falls_back_to_data_zones(tmp_path):
    fc = tmp_path / "333_FULL_CONTEXT.md"
    fc.write_text("", encoding="utf-8")
    zones = tmp_path / "333_DATA_ZONES.md"
    zones.write_text(_REAL_BODY, encoding="utf-8")

    assert demote_empty_full_context(fc) == zones


_ACCESS_DENIED_BODY = (
    "# 403 Forbidden\n\nAccess denied. Your institution does not have a "
    "subscription to this journal. " * 12
)


def test_demote_empty_full_context_skips_junk_sibling(tmp_path, caplog):
    """Above the floor is not the same as usable: an access-denied page is junk.

    The floor plus the abstract-only marker both pass here; the harvest-time
    content-quality validator is what refuses to stage it.
    """
    fc = tmp_path / "444_FULL_CONTEXT.md"
    fc.write_text("", encoding="utf-8")
    cleaned = tmp_path / "444_CLEANED.md"
    cleaned.write_text(_ACCESS_DENIED_BODY, encoding="utf-8")

    with caplog.at_level(logging.WARNING, logger="pipeline.source_quality"):
        selected = demote_empty_full_context(fc)

    assert selected == fc
    joined = " ".join(record.getMessage() for record in caplog.records)
    assert "not staged" in joined
    assert "444_CLEANED.md" in joined


def test_demote_empty_full_context_stages_article_shaped_sibling(tmp_path):
    fc = tmp_path / "555_FULL_CONTEXT.md"
    fc.write_text("", encoding="utf-8")
    cleaned = tmp_path / "555_CLEANED.md"
    cleaned.write_text(
        "# Abstract\n\nWe screened 212 probands.\n\n## Methods\n\n"
        + ("Sanger sequencing of KCNQ1 exons identified p.Ala561Val. " * 40)
        + "\n\n## Results\n\nCarriers were genotyped.\n\n## Discussion\n\ntext\n",
        encoding="utf-8",
    )

    assert demote_empty_full_context(fc) == cleaned


def test_demote_empty_full_context_skips_junk_then_takes_data_zones(tmp_path):
    """A failing CLEANED must not shadow a usable DATA_ZONES rendering."""
    fc = tmp_path / "666_FULL_CONTEXT.md"
    fc.write_text("", encoding="utf-8")
    (tmp_path / "666_CLEANED.md").write_text(_ACCESS_DENIED_BODY, encoding="utf-8")
    zones = tmp_path / "666_DATA_ZONES.md"
    zones.write_text(_REAL_BODY, encoding="utf-8")

    assert demote_empty_full_context(fc) == zones
