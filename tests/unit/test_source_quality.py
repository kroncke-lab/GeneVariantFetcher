from pipeline.source_quality import (
    SUPPLEMENT_SURFACE_STATUS_MARKER,
    is_reusable_fulltext_source,
    is_usable_fulltext_source,
)


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


def test_body_only_source_is_extractable_but_not_reusable(tmp_path):
    source = tmp_path / "123_FULL_CONTEXT.md"
    source.write_text(
        f"<!-- {SUPPLEMENT_SURFACE_STATUS_MARKER} unavailable -->\n\n"
        "# MAIN TEXT\n\n" + ("methods results discussion variants. " * 80),
        encoding="utf-8",
    )

    assert is_usable_fulltext_source(source) is True
    assert is_reusable_fulltext_source(source) is False
