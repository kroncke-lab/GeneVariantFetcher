"""Tests for preprocessing behavior in pipeline steps."""

from pipeline.steps import preprocess_papers


def test_preprocess_writes_cleaned_file_without_overwriting_source(tmp_path):
    """Preprocessing should write *_CLEANED.md and leave *_FULL_CONTEXT.md intact."""
    harvest_dir = tmp_path / "pmc_fulltext"
    harvest_dir.mkdir()

    source_file = harvest_dir / "12345678_FULL_CONTEXT.md"
    original_text = "## Introduction\n\nVariant p.Arg534Cys was reported.\n"
    source_file.write_text(original_text, encoding="utf-8")

    result = preprocess_papers(harvest_dir=harvest_dir, gene_symbol="KCNH2")

    cleaned_file = harvest_dir / "12345678_CLEANED.md"
    assert result.success is True
    assert result.stats.get("processed") == 1
    assert result.stats.get("cleaned_files_written") == 1
    assert cleaned_file.exists()

    # Source remains untouched (non-destructive mode).
    assert source_file.read_text(encoding="utf-8") == original_text
    assert cleaned_file.read_text(encoding="utf-8")
