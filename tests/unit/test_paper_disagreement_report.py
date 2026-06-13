import pytest

from scripts.recall_audit.common import find_full_contexts
from scripts.recall_audit.paper_disagreement_report import (
    _build_context_index,
    _failure_class,
    _source_info,
)


def _write_context(path, marker, repeats=220):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        (f"# MAIN TEXT\n\nAbstract\nIntroduction\nMethods\nResults\n{marker}\n")
        * repeats,
        encoding="utf-8",
    )


def test_source_info_uses_consumed_cleaned_source_before_latest_context(tmp_path):
    results_dir = tmp_path / "results"
    cleaned = tmp_path / "run" / "pmc_fulltext" / "123_CLEANED.md"
    latest = results_dir / "KCNH2" / "newer" / "pmc_fulltext" / "123_FULL_CONTEXT.md"
    _write_context(cleaned, "older recovered browser text", repeats=120)
    _write_context(latest, "newer recovered full context with extra table rows")
    context_index = _build_context_index(results_dir)

    info = _source_info(
        gene="KCNH2",
        pmid="123",
        metadata={"source_file": str(cleaned)},
        results_dir=results_dir,
        context_index=context_index,
    )

    assert info["context_path"] == str(cleaned)
    assert info["available_context_path"] == str(latest)
    assert info["context_bytes"] == cleaned.stat().st_size
    assert info["available_context_bytes"] == latest.stat().st_size
    assert info["data_available"] is True
    assert info["source_desync"] is True


def test_find_full_contexts_skips_broken_context_symlink(tmp_path):
    results_dir = tmp_path / "results"
    context_dir = results_dir / "KCNH2" / "run" / "pmc_fulltext"
    context_dir.mkdir(parents=True)
    broken = context_dir / "123_FULL_CONTEXT.md"
    try:
        broken.symlink_to(context_dir / "missing_FULL_CONTEXT.md")
    except OSError:
        pytest.skip("symlinks unavailable on this filesystem")
    good = results_dir / "KCNH2" / "newer" / "pmc_fulltext" / "123_FULL_CONTEXT.md"
    _write_context(good, "usable context")

    assert find_full_contexts(results_dir, "KCNH2", "123") == [good]


def test_source_info_marks_abstract_json_with_later_context_as_desync(tmp_path):
    results_dir = tmp_path / "results"
    extraction_json = tmp_path / "run" / "extractions" / "KCNH2_PMID_456.json"
    extraction_json.parent.mkdir(parents=True)
    extraction_json.write_text("{}", encoding="utf-8")
    latest = results_dir / "KCNH2" / "newer" / "pmc_fulltext" / "456_FULL_CONTEXT.md"
    _write_context(latest, "later full text source")
    context_index = _build_context_index(results_dir)

    info = _source_info(
        gene="KCNH2",
        pmid="456",
        metadata={
            "source_file": str(extraction_json),
            "abstract_only": True,
            "notes": "Abstract-only extraction; full text unavailable",
        },
        results_dir=results_dir,
        context_index=context_index,
    )

    assert info["source_status"] == "abstract_only"
    assert info["data_available"] is False
    assert info["context_path"] == ""
    assert info["available_context_path"] == str(latest)
    assert info["source_desync"] is True
    assert info["source_unbound"] is True
    assert (
        _failure_class(
            {
                **info,
                "missing_rows": 10,
                "gold_rows": 12,
                "count_mismatch_rows": 0,
                "extra_sqlite_rows": 0,
                "abs_affected_diff": 0,
                "abs_unaffected_diff": 0,
            }
        )
        == "stale_source_desync"
    )


def test_source_info_marks_unbound_nonabstract_source_without_calling_it_stale(
    tmp_path,
):
    results_dir = tmp_path / "results"
    extraction_json = tmp_path / "run" / "extractions" / "KCNH2_PMID_654.json"
    extraction_json.parent.mkdir(parents=True)
    extraction_json.write_text("{}", encoding="utf-8")
    latest = results_dir / "KCNH2" / "newer" / "pmc_fulltext" / "654_FULL_CONTEXT.md"
    _write_context(latest, "later full text source")
    context_index = _build_context_index(results_dir)

    info = _source_info(
        gene="KCNH2",
        pmid="654",
        metadata={
            "source_file": str(extraction_json),
            "notes": "Extractor reported a variant-rich clinical study",
        },
        results_dir=results_dir,
        context_index=context_index,
    )

    assert info["source_status"] == "not_attempted"
    assert info["source_unbound"] is True
    assert info["source_desync"] is False
    assert (
        _failure_class(
            {
                **info,
                "missing_rows": 10,
                "gold_rows": 12,
                "count_mismatch_rows": 0,
                "extra_sqlite_rows": 0,
                "abs_affected_diff": 0,
                "abs_unaffected_diff": 0,
            }
        )
        == "source_unbound_available"
    )


def test_source_info_finds_sibling_context_for_unbound_extraction_json(tmp_path):
    results_dir = tmp_path / "unrelated_results"
    extraction_json = tmp_path / "run" / "extractions" / "KCNQ1_PMID_32893267.json"
    extraction_json.parent.mkdir(parents=True)
    extraction_json.write_text("{}", encoding="utf-8")
    sibling = tmp_path / "run" / "pmc_fulltext" / "32893267_FULL_CONTEXT.md"
    _write_context(sibling, "sibling run source")

    info = _source_info(
        gene="KCNQ1",
        pmid="32893267",
        metadata={"source_file": str(extraction_json)},
        results_dir=results_dir,
        context_index={},
    )

    assert info["context_path"] == ""
    assert info["available_context_path"] == str(sibling)
    assert info["available_source_status"] == "recovered_pmc"
    assert info["source_unbound"] is True
    assert info["source_desync"] is False


def test_source_info_does_not_prefer_larger_abstract_fallback(tmp_path):
    results_dir = tmp_path / "results"
    cleaned = tmp_path / "run" / "pmc_fulltext" / "789_CLEANED.md"
    _write_context(cleaned, "usable cleaned text", repeats=120)
    fallback = results_dir / "KCNH2" / "newer" / "pmc_fulltext" / "789_FULL_CONTEXT.md"
    fallback.parent.mkdir(parents=True, exist_ok=True)
    fallback.write_text(
        (
            "# ABSTRACT-ONLY FALLBACK\n\n"
            "> **WARNING:** Full text could not be retrieved for PMID 789.\n"
            "> This document contains only the PubMed abstract and metadata.\n"
        )
        * 260,
        encoding="utf-8",
    )
    context_index = _build_context_index(results_dir)

    info = _source_info(
        gene="KCNH2",
        pmid="789",
        metadata={"source_file": str(cleaned)},
        results_dir=results_dir,
        context_index=context_index,
    )

    assert info["context_path"] == str(cleaned)
    assert info["available_context_path"] == str(fallback)
    assert info["available_source_status"] == "abstract_only"
    assert info["source_desync"] is False
