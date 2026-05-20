"""Tests for the Data Scout CLI wrapper."""

from cli.scout import run_scout
from pipeline.data_scout import select_scout_source_path
from utils.manifest import Status


def test_run_scout_accepts_string_paths(tmp_path):
    """Workflow callers pass strings; run_scout should normalize them to Paths."""

    input_dir = tmp_path / "fulltext"
    output_dir = tmp_path / "scout"
    input_dir.mkdir()

    (input_dir / "12345_FULL_CONTEXT.md").write_text(
        "\n".join(
            [
                "# Results",
                "Table 1. KCNH2 variants in patients",
                "| Variant | Carriers | Affected |",
                "|---|---:|---:|",
                "| A561V | 2 | 2 |",
                "The KCNH2 proband carried A561V.",
            ]
        ),
        encoding="utf-8",
    )

    manifest = run_scout(
        input_path=str(input_dir),
        output_dir=str(output_dir),
        gene="KCNH2",
    )

    assert manifest.entries
    assert manifest.entries[0].status == Status.SUCCESS
    assert (output_dir / "12345_DATA_ZONES.md").exists()
    assert (output_dir / "scout_manifest.json").exists()


def test_select_scout_source_prefers_cleaned_for_oversized_context(tmp_path):
    """Large raw FULL_CONTEXT files should not bypass deterministic cleanup."""

    full_context = tmp_path / "12345_FULL_CONTEXT.md"
    cleaned = tmp_path / "12345_CLEANED.md"
    full_context.write_text("raw raw raw", encoding="utf-8")
    cleaned.write_text("clean KCNH2 variant table", encoding="utf-8")

    selected = select_scout_source_path(full_context, prefer_cleaned_above_chars=5)

    assert selected == cleaned
