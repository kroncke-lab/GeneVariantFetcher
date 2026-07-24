from __future__ import annotations

import hashlib
import json
from types import SimpleNamespace

import pytest

from pipeline.sol_extractor import (
    EvidenceIndex,
    FULL_CONTEXT_MAX_CHARS,
    SolExtractor,
    build_compact_context,
    extract_paper,
    strip_irrelevant_context,
)


def _record(
    *,
    gene: str,
    protein: str | None,
    cdna: str | None,
    total: int | None,
    locator: str,
    quote: str,
) -> dict:
    return {
        "gene_symbol": gene,
        "cdna_notation": cdna,
        "protein_notation": protein,
        "genomic_position": None,
        "clinical_significance": "pathogenic",
        "phenotype": "hereditary cancer",
        "total_carriers": total,
        "affected": None,
        "unaffected": None,
        "uncertain": None,
        "count_semantics": {
            "total_carriers": ("human_variant_carriers" if total is not None else None),
            "affected": None,
            "unaffected": None,
            "uncertain": None,
        },
        "locator": locator,
        "evidence_quote": quote,
    }


class StaticRunner:
    def __init__(self, value: dict) -> None:
        self.value = value
        self.calls = 0
        self.kwargs = None

    def run_tool_loop(self, **kwargs):
        self.calls += 1
        self.kwargs = kwargs
        return SimpleNamespace(
            ok=True,
            status="completed",
            value=self.value,
            output_text=None,
            error=None,
            response_id="resp_test",
            usage=None,
            telemetry=None,
        )


class FailingIfCalledRunner:
    def __init__(self) -> None:
        self.calls = 0

    def run_tool_loop(self, **kwargs):
        self.calls += 1
        raise AssertionError("provider must not be called for an animal-only paper")


class ParseAndAcceptRunner:
    def __init__(self, *, row_end: int) -> None:
        self.row_end = row_end
        self.calls = 0

    def run_tool_loop(self, **kwargs):
        self.calls += 1
        parsed = kwargs["tool_handlers"]["parse_variant_table"](
            {
                "table_id": "T1",
                "row_start": 1,
                "row_end": self.row_end,
                "column_mapping": {
                    "affected": None,
                    "cdna": None,
                    "clinical_significance": None,
                    "gene": 0,
                    "patient_count": 2,
                    "phenotype": None,
                    "protein": 1,
                    "unaffected": None,
                    "uncertain": None,
                },
                "human_subjects": True,
                "per_variant": True,
                "count_semantics": {
                    "affected": None,
                    "patient_count": "human_variant_carriers",
                    "unaffected": None,
                    "uncertain": None,
                },
            }
        )
        return SimpleNamespace(
            ok=True,
            status="completed",
            value={
                "eligibility": {
                    "human_clinical": True,
                    "reason": "human carrier table",
                    "species": "human",
                },
                "accepted_table_batches": [{"batch_id": parsed["batch_id"]}],
                "evidence_records": [],
            },
            output_text=None,
            error=None,
            response_id="resp_parse_test",
            usage=None,
            telemetry=None,
        )


def test_strip_references_preserves_later_supplement_and_table() -> None:
    source = """# Paper title

## Results
BRCA2 c.68_69delAG was observed in two patients.

## References
1. Irrelevant BRCA2 review.
2. Another citation.

### Reference notes
This nested reference material should also disappear.

## Supplementary Information
Supplement Table S1 remains evidence.

| Gene | Variant | Carriers |
| --- | --- | --- |
| BRCA2 | c.68_69delAG | 2 |
"""

    cleaned = strip_irrelevant_context(source)

    assert "Irrelevant BRCA2 review" not in cleaned
    assert "nested reference material" not in cleaned
    assert "Supplementary Information" in cleaned
    assert "Supplement Table S1 remains evidence" in cleaned
    assert "| BRCA2 | c.68_69delAG | 2 |" in cleaned


def test_evidence_index_table_chunks_are_bounded_and_stable() -> None:
    rows = "\n".join(f"| BRCA2 | p.Arg{row}Trp | {row} |" for row in range(1, 106))
    source = f"""# Human BRCA2 cohort

Table 1. Human variant carriers
| Gene | Protein | Carriers |
| --- | --- | --- |
{rows}
"""
    index = EvidenceIndex(source)

    assert [table.table_id for table in index.tables] == ["T1"]
    first = index.read_table(table_id="T1", row_start=2, row_end=105)
    second = index.read_table(table_id="T1", row_start=2, row_end=105)

    assert first == second
    assert first["row_start"] == 2
    assert first["row_end"] == 101
    assert first["truncated"] is True
    assert len(first["rows"]) == EvidenceIndex.MAX_TABLE_ROWS
    assert first["rows"][0]["locator"] == "T1:R2"
    assert index.resolve_locator("T1:R2") == "| BRCA2 | p.Arg2Trp | 2 |"
    source_locator = first["rows"][0]["source_locator"]
    assert index.resolve_locator(source_locator) == "| BRCA2 | p.Arg2Trp | 2 |"


def test_normal_source_is_embedded_in_full_with_resolvable_span_labels() -> None:
    source = "Human BRCA2 cohort\n" + ("A" * (EvidenceIndex.MAX_SOURCE_SPAN + 7))
    index = EvidenceIndex(source)

    payload = json.loads(
        build_compact_context(
            index=index,
            gene_symbol="BRCA2",
            pmid="12345678",
        )
    )

    assert payload["context_delivery"] == {
        "mode": "full_cleaned_source",
        "complete": True,
        "source_span_count": 2,
    }
    assert "".join(span["text"] for span in payload["source_spans"]) == source
    for span in payload["source_spans"]:
        assert len(span["text"]) <= EvidenceIndex.MAX_SOURCE_SPAN
        assert index.resolve_locator(span["locator"]) == span["text"]


def test_oversized_source_stays_behind_bounded_tools() -> None:
    source = "X" * (FULL_CONTEXT_MAX_CHARS + 1)
    payload = json.loads(
        build_compact_context(
            index=EvidenceIndex(source),
            gene_symbol="BRCA2",
            pmid="12345678",
        )
    )

    assert payload["context_delivery"] == {
        "mode": "indexed_tools",
        "complete": False,
        "source_span_count": 0,
    }
    assert payload["source_spans"] == []
    assert "exceeds the direct-context threshold" in payload["note"]


def test_canine_brca2_is_short_circuited_before_provider_call() -> None:
    source = """# Germline BRCA2 variation in canine mammary tumors

## Abstract
This veterinary study sequenced female dogs with canine mammary tumors.
Canine blood and dog tumor samples were analyzed. No human subjects enrolled.

| Dog | BRCA2 variant | Tumor count |
| --- | --- | --- |
| D1 | p.Lys805Arg | 6 |
"""
    runner = FailingIfCalledRunner()

    output = SolExtractor(runner).extract(
        source_text=source,
        gene_symbol="BRCA2",
        pmid="19944633",
    )

    assert runner.calls == 0
    assert output["variants"] == []
    metadata = output["extraction_metadata"]
    assert metadata["eligibility"]["human_clinical"] is False
    assert metadata["eligibility"]["deterministic_short_circuit"] is True
    assert metadata["total_variants_found"] == 0


def test_hallucinated_count_is_cleared_but_grounded_identity_is_retained() -> None:
    quote = "BRCA2 c.68_69delAG was found in 2 human carriers."
    source = f"""# Human BRCA2 cohort

## Abstract
{quote}
"""
    start = source.index(quote)
    runner = StaticRunner(
        {
            "eligibility": {
                "human_clinical": True,
                "reason": "human cohort",
                "species": "human",
            },
            "accepted_table_batches": [],
            "evidence_records": [
                _record(
                    gene="BRCA2",
                    protein=None,
                    cdna="c.68_69delAG",
                    total=9,
                    locator=f"S:{start}:{start + len(quote)}",
                    quote=quote,
                )
            ],
        }
    )

    output = SolExtractor(runner).extract(
        source_text=source,
        gene_symbol="BRCA2",
        pmid="12345678",
    )

    assert runner.calls == 1
    assert len(output["variants"]) == 1
    assert output["variants"][0]["cdna_notation"] == "c.68_69delAG"
    assert output["variants"][0]["penetrance_data"]["total_carriers_observed"] is None
    rejected = output["extraction_metadata"]["sol_audit"]["rejected_raw_facts"]
    assert len(rejected) == 1
    assert rejected[0]["reason"] == "unsupported_count_fields_cleared"
    assert "total_carriers=9 cleared" in rejected[0]["details"][0]


def test_hallucinated_quote_is_rejected() -> None:
    source = """# Human BRCA2 cohort

## Abstract
BRCA2 p.Arg176Trp was found in 2 human carriers.
"""
    runner = StaticRunner(
        {
            "eligibility": {
                "human_clinical": True,
                "reason": "human cohort",
                "species": "human",
            },
            "accepted_table_batches": [],
            "evidence_records": [
                _record(
                    gene="BRCA2",
                    protein="p.Arg176Trp",
                    cdna=None,
                    total=200,
                    locator="S:0:50",
                    quote="p.Arg176Trp appeared in 200 carriers.",
                )
            ],
        }
    )

    output = SolExtractor(runner).extract(
        source_text=source,
        gene_symbol="BRCA2",
        pmid="12345678",
    )

    assert output["variants"] == []
    rejected = output["extraction_metadata"]["sol_audit"]["rejected_raw_facts"]
    assert "not an exact substring" in rejected[0]["reason"]


def test_parsed_duplicate_restatements_are_not_blindly_summed() -> None:
    source = """# Human BRCA2 cohort

Table 1. Human BRCA2 variant carriers
| Gene | Protein | Human carriers |
| --- | --- | --- |
| BRCA2 | p.Arg176Trp | 2 |
| BRCA2 | p.Arg176Trp | 3 |
"""
    runner = ParseAndAcceptRunner(row_end=2)

    output = SolExtractor(runner).extract(
        source_text=source,
        gene_symbol="BRCA2",
        pmid="12345678",
    )

    assert runner.calls == 1
    assert len(output["variants"]) == 1
    variant = output["variants"][0]
    assert variant["penetrance_data"]["total_carriers_observed"] is None
    assert variant["patients"]["count"] is None
    rejected = output["extraction_metadata"]["sol_audit"]["rejected_raw_facts"]
    assert any(
        item["reason"] == "conflicting_duplicate_counts_not_summed" for item in rejected
    )


def test_wrong_direct_count_semantics_clears_count_not_variant_identity() -> None:
    quote = "BRCA2 p.Arg176Trp was observed in 2 human carriers."
    source = f"# Human cohort\n\n## Abstract\n{quote}\n"
    start = source.index(quote)
    raw = _record(
        gene="BRCA2",
        protein="p.Arg176Trp",
        cdna=None,
        total=2,
        locator=f"S:{start}:{start + len(quote)}",
        quote=quote,
    )
    raw["count_semantics"]["total_carriers"] = "human_variant_affected"
    runner = StaticRunner(
        {
            "eligibility": {
                "human_clinical": True,
                "reason": "human cohort",
                "species": "human",
            },
            "accepted_table_batches": [],
            "evidence_records": [raw],
        }
    )

    output = SolExtractor(runner).extract(
        source_text=source,
        gene_symbol="BRCA2",
        pmid="12345678",
    )

    assert len(output["variants"]) == 1
    assert output["variants"][0]["protein_notation"] == "p.Arg176Trp"
    assert output["variants"][0]["penetrance_data"]["total_carriers_observed"] is None
    audit = output["extraction_metadata"]["sol_audit"]["rejected_raw_facts"]
    assert audit[0]["reason"] == "unsupported_count_fields_cleared"
    assert "count_semantics must be 'human_variant_carriers'" in audit[0]["details"][0]


def test_extract_paper_returns_driver_telemetry_and_checks_source_hash() -> None:
    source = "# Human BRCA2 cohort\n\n## Abstract\nNo target variants were reported.\n"
    source_hash = hashlib.sha256(source.encode("utf-8")).hexdigest()
    runner = StaticRunner(
        {
            "eligibility": {
                "human_clinical": True,
                "reason": "human cohort without a reported target variant",
                "species": "human",
            },
            "accepted_table_batches": [],
            "evidence_records": [],
        }
    )

    extraction, telemetry = extract_paper(
        gene="BRCA2",
        pmid="12345678",
        source_text=source,
        source_sha256=source_hash,
        reasoning_effort="high",
        runner=runner,
    )

    assert extraction["tables_processed"] == []
    assert telemetry["model"] == "azure_ai/gpt-5.6-sol"
    assert telemetry["requested_reasoning_effort"] == "high"
    assert telemetry["effective_reasoning_effort"] == "high"
    assert telemetry["source_sha256"] == source_hash
    assert telemetry["cleaned_source_sha256"]
    assert telemetry["status"] == "success"
    assert runner.kwargs["max_output_tokens"] == 100_000
    limits = runner.kwargs["limits"]
    assert limits.max_rounds == 32
    assert limits.max_tool_calls == 128
    assert limits.max_tool_output_chars == 30_000
    assert limits.max_total_tool_output_chars == 640_000
    initial = json.loads(runner.kwargs["initial_input"])
    assert initial["context_delivery"]["mode"] == "full_cleaned_source"
    assert "No target variants were reported." in initial["source_spans"][0]["text"]

    with pytest.raises(ValueError, match="source_sha256 mismatch"):
        extract_paper(
            gene="BRCA2",
            pmid="12345678",
            source_text=source,
            source_sha256="0" * 64,
            reasoning_effort="high",
            runner=runner,
        )
