"""Tests for pipeline.table_router (table enumeration + router + parser)."""

import json
from pathlib import Path

import pytest

from pipeline.table_router import (
    MarkdownTable,
    enumerate_markdown_tables,
    extract_via_router,
    parse_routed_table,
    parse_router_response,
)


SAMPLE_PAPER = """# Title

Some prose. There is no table here.

Table 1. Variant carriers in the cohort

| cdna | protein | n carriers | affected | unaffected |
|------|---------|------------|----------|------------|
| c.2690A>C | p.K897T | 7 | 3 | 4 |
| c.3140G>A | p.R1047L | 14 | 4 | 10 |

Some more discussion text.

Table 2. Functional assays (not a variant table)

| construct | tail current | activation V1/2 |
|-----------|--------------|-----------------|
| WT | 100% | -32 |
| K897T | 95% | -30 |
"""


def _stub_response(content: str):
    """Build a minimal object that quacks like litellm's response shape."""

    class _Choice:
        def __init__(self, content: str) -> None:
            self.message = type("M", (), {"content": content})

    class _Resp:
        def __init__(self, content: str) -> None:
            self.choices = [_Choice(content)]

    return _Resp(content)


def test_enumerate_finds_captioned_tables():
    # Both Table 1 and Table 2 have "Table N." captions, so both survive the
    # cheap Python-side filter. The downstream LLM router is the actual judge
    # of which table contains variant data.
    tables = enumerate_markdown_tables(SAMPLE_PAPER)
    assert len(tables) == 2
    t = tables[0]
    assert t.caption is not None and t.caption.lower().startswith("table 1.")
    assert len(t.header_cells) == 5
    assert len(t.data_lines) == 2


def test_enumerate_drops_caption_less_pseudo_tables():
    # Pseudo-table (no caption, no variant-ish keyword in header) gets filtered.
    pseudo = """Some prose that contains | pipes | for | formatting

| word1 | word2 | word3 |
|-------|-------|-------|
| foo   | bar   | baz   |
"""
    assert enumerate_markdown_tables(pseudo) == []


def test_enumerate_disable_filter_keeps_all():
    tables = enumerate_markdown_tables(SAMPLE_PAPER, only_variant_like=False)
    assert len(tables) == 2


def test_parse_router_response_handles_fenced_json():
    raw = """Some prose preface.\n```json\n{"variant_tables":[\n  {"table_id":"T1",\n   "column_mapping":{"cdna":0,"protein":1,"patient_count":2,"affected":3,"unaffected":4},\n   "confidence":0.95}\n]}\n```\nTrailing chatter."""
    routed = parse_router_response(raw)
    assert len(routed) == 1
    assert routed[0].table_id == "T1"
    assert routed[0].column_mapping == {
        "cdna": 0,
        "protein": 1,
        "patient_count": 2,
        "affected": 3,
        "unaffected": 4,
    }


def test_parse_router_response_skips_non_json_braced_preface():
    raw = """I considered table {ID: T1}, then chose the strict schema.

{"variant_tables":[
  {"table_id":"T1",
   "column_mapping":{"cdna":0,"protein":1,"patient_count":2},
   "confidence":0.95}
]}"""
    routed = parse_router_response(raw)
    assert len(routed) == 1
    assert routed[0].table_id == "T1"


def test_parse_router_response_drops_invalid_mappings():
    raw = json.dumps(
        {
            "variant_tables": [
                {"table_id": "T1", "column_mapping": {"junk_field": 0}},
                {
                    "table_id": "T2",
                    "column_mapping": {"cdna": "not_an_int"},
                },
                {"table_id": "T3", "column_mapping": {"cdna": 0}},
            ]
        }
    )
    routed = parse_router_response(raw)
    # T1: no allowed fields → drop
    # T2: cdna is not int → drop after coercion fails
    # T3: cdna only, no count → drop (must have notation AND count)
    assert routed == []


def test_parse_routed_table_extracts_variants():
    tables = enumerate_markdown_tables(SAMPLE_PAPER)
    assert tables, "expected at least one variant-like table"
    table = tables[0]
    mapping = {
        "cdna": 0,
        "protein": 1,
        "patient_count": 2,
        "affected": 3,
        "unaffected": 4,
    }
    variants = parse_routed_table(table, mapping, "KCNH2")
    assert len(variants) == 2
    by_protein = {v["protein_notation"]: v for v in variants}
    assert "p.K897T" in by_protein
    assert "p.R1047L" in by_protein
    pen = by_protein["p.R1047L"]["penetrance_data"]
    assert pen["total_carriers_observed"] == 14
    assert pen["affected_count"] == 4
    assert pen["unaffected_count"] == 10


def test_parse_routed_table_filters_multigene_and_derives_total():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Arrhythmia panel carriers",
        header_line="| Gene | Protein | Affected | Unaffected |",
        header_cells=["Gene", "Protein", "Affected", "Unaffected"],
        data_lines=[
            "| KCNQ1 | p.Met1? | 2 | 3 |",
            "| SCN5A | p.Ala226Asp | 9 | 0 |",
        ],
        char_start=0,
        char_end=160,
    )
    mapping = {"gene": 0, "protein": 1, "affected": 2, "unaffected": 3}

    variants = parse_routed_table(table, mapping, "KCNQ1")

    assert len(variants) == 1
    assert variants[0]["protein_notation"] == "p.Met1?"
    pen = variants[0]["penetrance_data"]
    assert pen["total_carriers_observed"] == 5
    assert pen["affected_count"] == 2
    assert pen["unaffected_count"] == 3


def test_extract_via_router_infers_one_carrier_per_clinical_row_without_count():
    def stub(**_):
        raise AssertionError("LLM router should not run for row-level clinical tables")

    paper = """Table 1. KCNH2 mutation-positive probands

| Patient | Mutation | QTc | Phenotype |
|---------|----------|-----|-----------|
| P1 | p.Arg176Trp | 510 | syncope |
| P2 | p.Asn629Ser | 430 | unaffected |
"""

    result = extract_via_router(
        paper, "KCNH2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )

    assert result["used_fallback"] is False
    assert len(result["variants"]) == 2
    by_protein = {v["protein_notation"]: v for v in result["variants"]}
    affected = by_protein["p.Arg176Trp"]["penetrance_data"]
    unaffected = by_protein["p.Asn629Ser"]["penetrance_data"]
    assert affected["total_carriers_observed"] == 1
    assert affected["affected_count"] == 1
    assert affected["unaffected_count"] == 0
    assert unaffected["total_carriers_observed"] == 1
    assert unaffected["affected_count"] == 0
    assert unaffected["unaffected_count"] == 1


def test_row_level_inference_requires_clinical_context():
    def stub(**_):
        return _stub_response(json.dumps({"variant_tables": []}))

    paper = """Table 1. Protein nomenclature reference

| Variant | Domain |
|---------|--------|
| p.Arg176Trp | N terminus |
| p.Asn629Ser | pore |
"""

    result = extract_via_router(
        paper, "KCNH2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )

    assert result["variants"] == []


def test_parse_routed_table_rejects_misrouted_gwas_allele_rows():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. GWAS lead SNPs",
        header_line="| Locus | SNV | AA | n |",
        header_cells=["Locus", "SNV", "AA", "n"],
        data_lines=[
            "| KCNH2 | rs113843864 | A | 29762 |",
            "| KCNH2 | rs2072412 | C | 31607149 |",
        ],
        char_start=0,
        char_end=100,
    )
    mapping = {"gene": 0, "protein": 2, "patient_count": 3}

    variants = parse_routed_table(table, mapping, "KCNH2")

    assert variants == []


def test_extract_via_router_end_to_end_with_stub():
    captured = {}

    def stub(*, model, messages, temperature, max_tokens):
        captured["model"] = model
        # Verify the prompt contains both detected tables only after filter
        captured["prompt"] = messages[-1]["content"]
        return _stub_response(
            json.dumps(
                {
                    "variant_tables": [
                        {
                            "table_id": "T1",
                            "column_mapping": {
                                "cdna": 0,
                                "protein": 1,
                                "patient_count": 2,
                                "affected": 3,
                                "unaffected": 4,
                            },
                            "confidence": 0.99,
                        }
                    ]
                }
            )
        )

    result = extract_via_router(
        SAMPLE_PAPER, "KCNH2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )
    assert captured["model"] == "azure_ai/Kimi-K2.6-1"
    # T1 is deterministically mapped and therefore not sent to the LLM router.
    assert "Variant carriers in the cohort" not in captured["prompt"]
    # The assay table is still offered to the router because it has a caption.
    assert "Functional assays" in captured["prompt"]
    assert len(result["variants"]) == 2
    assert result["used_fallback"] is False
    assert result["error"] is None


def test_extract_via_router_skips_llm_for_unambiguous_clinical_table():
    def stub(**_):
        raise AssertionError("LLM router should not run for unambiguous headers")

    paper = """Table 1. KCNH2 clinical carriers

| Gene | Protein | Affected | Unaffected |
|------|---------|----------|------------|
| KCNH2 | p.Gly572Ser | 2 | 1 |
"""

    result = extract_via_router(
        paper, "KCNH2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )

    assert result["used_fallback"] is False
    assert len(result["variants"]) == 1
    assert result["variants"][0]["protein_notation"] == "p.Gly572Ser"


def test_extract_via_router_falls_back_when_no_tables():
    def stub(**_):
        raise AssertionError("LLM should not be called when no tables detected")

    result = extract_via_router(
        "Just prose. No tables at all.",
        "KCNH2",
        model="azure_ai/Kimi-K2.6-1",
        llm_caller=stub,
    )
    assert result["used_fallback"] is True
    assert result["variants"] == []


def test_extract_via_router_falls_back_when_router_returns_empty():
    def stub(**_):
        return _stub_response('{"variant_tables": []}')

    paper = """Table 1. Functional assays

| construct | tail current | activation V1/2 |
|-----------|--------------|-----------------|
| WT | 100% | -32 |
| K897T | 95% | -30 |
"""

    result = extract_via_router(
        paper, "KCNH2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )
    assert result["used_fallback"] is True
    assert result["variants"] == []
    assert result["error"] is None
