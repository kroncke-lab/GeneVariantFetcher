import json

from harvesting.orchestrator import full_context_needs_retry


def test_abstract_only_full_context_is_retried(tmp_path):
    full_context = tmp_path / "15840476_FULL_CONTEXT.md"
    full_context.write_text(
        "# ABSTRACT-ONLY FALLBACK\n\n## Abstract\nVariants listed only in abstract.",
        encoding="utf-8",
    )

    assert full_context_needs_retry(full_context, tmp_path) is True


def test_abstract_only_status_is_retried_even_without_marker(tmp_path):
    status_dir = tmp_path / "pmid_status"
    status_dir.mkdir()
    (status_dir / "15840476.json").write_text(
        json.dumps({"status": "abstract_only", "source": "pubmed_abstract"}),
        encoding="utf-8",
    )
    full_context = tmp_path / "15840476_FULL_CONTEXT.md"
    full_context.write_text(
        "# MAIN TEXT\n\n" + ("abstract text\n" * 100), encoding="utf-8"
    )

    assert full_context_needs_retry(full_context, tmp_path) is True


def test_thin_publisher_shell_full_context_is_retried(tmp_path):
    full_context = tmp_path / "10973849_FULL_CONTEXT.md"
    full_context.write_text(
        "\n".join(
            [
                "# MAIN TEXT",
                "",
                "## Spectrum of Mutations in Long",
                "",
                "### Abstract",
                "",
                "### Methods",
                "",
                "eLetters should relate to an article recently published in the journal.",
            ]
        ),
        encoding="utf-8",
    )

    assert full_context_needs_retry(full_context, tmp_path) is True


def test_complete_full_context_is_not_retried(tmp_path):
    full_context = tmp_path / "23098067_FULL_CONTEXT.md"
    full_context.write_text(
        "# MAIN TEXT\n\n## Results\n\n"
        + ("full text with variant table content\n" * 400),
        encoding="utf-8",
    )

    assert full_context_needs_retry(full_context, tmp_path) is False


def test_complete_full_context_wins_over_stale_abstract_status(tmp_path):
    status_dir = tmp_path / "pmid_status"
    status_dir.mkdir()
    (status_dir / "19038855.json").write_text(
        json.dumps({"status": "abstract_only", "source": "pubmed_abstract"}),
        encoding="utf-8",
    )
    full_context = tmp_path / "19038855_FULL_CONTEXT.md"
    full_context.write_text(
        "# REHARVESTED FULL TEXT\n\n# MAIN TEXT\n\n"
        + ("full text with methods results and tables\n" * 300),
        encoding="utf-8",
    )

    assert full_context_needs_retry(full_context, tmp_path) is False
