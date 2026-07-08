from pathlib import Path

from pipeline.extraction_priority import ExtractionCandidate, PriorityResult
from pipeline.extraction_triage import (
    apply_triage_filter,
    deterministic_triage,
    triage_priority_result,
)


def _candidate(pmid, *, signals, score=100, title="Target gene cohort"):
    return ExtractionCandidate(
        pmid=pmid,
        source_kind="fulltext",
        source_file=f"/tmp/{pmid}_FULL_CONTEXT.md",
        score=score,
        title=title,
        rank=int(pmid),
        selected=True,
        signals=signals,
    )


def test_deterministic_triage_extracts_gene_variant_count_tables():
    decision = deterministic_triage(
        _candidate(
            "1",
            signals={
                "gene_title_mentions": 1,
                "gene_abstract_mentions": 2,
                "gene_body_mentions": 20,
                "gene_variant_lines": 4,
                "gene_carrier_lines": 3,
                "table_mentions": 2,
                "original_data_mentions": 5,
            },
        )
    )

    assert decision.decision == "extract_now"
    assert "target_gene_variants" in decision.useful_info_types


def test_deterministic_triage_skips_review_without_variant_evidence():
    decision = deterministic_triage(
        _candidate(
            "2",
            title="A review of disease mechanisms",
            signals={
                "gene_title_mentions": 1,
                "gene_body_mentions": 10,
                "gene_variant_lines": 0,
                "review_title_penalty": True,
            },
        )
    )

    assert decision.decision == "skip"
    assert "review_or_protocol" in decision.false_positive_risks


def test_triage_filter_keeps_only_extract_now(tmp_path):
    keep = tmp_path / "1_FULL_CONTEXT.md"
    drop = tmp_path / "2_FULL_CONTEXT.md"
    keep.write_text("BRCA1 c.68_69delAG carriers n=12", encoding="utf-8")
    drop.write_text("review", encoding="utf-8")
    candidates = [
        _candidate(
            "1",
            signals={
                "gene_title_mentions": 1,
                "gene_body_mentions": 5,
                "gene_variant_lines": 2,
                "gene_carrier_lines": 2,
                "table_mentions": 1,
                "original_data_mentions": 1,
            },
        ),
        _candidate(
            "2",
            signals={"gene_title_mentions": 1, "review_title_penalty": True},
            title="BRCA1 review",
        ),
    ]
    priority = PriorityResult(
        selected_markdown_files=[keep, drop],
        selected_abstract_papers=[],
        candidates=candidates,
    )

    triage = triage_priority_result(
        priority,
        gene_symbol="BRCA1",
        mode="deterministic",
        report_dir=tmp_path / "triage",
    )
    filtered = apply_triage_filter(priority, triage)

    assert [p.name for p in filtered.selected_markdown_files] == ["1_FULL_CONTEXT.md"]
    assert (tmp_path / "triage" / "triage_candidates.tsv").exists()
    assert (tmp_path / "triage" / "triage_pmids.txt").read_text(
        encoding="utf-8"
    ) == "1\n"
