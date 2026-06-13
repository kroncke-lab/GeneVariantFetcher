from pipeline.source_quality import is_usable_fulltext_source


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
