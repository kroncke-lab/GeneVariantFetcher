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
    assert "T1" in captured["prompt"]
    # The assay table is offered to the router (since it has a "Table N."
    # caption); the stub picks only T1, so only those variants come through.
    assert "Functional assays" in captured["prompt"]
    assert len(result["variants"]) == 2
    assert result["used_fallback"] is False
    assert result["error"] is None


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

    result = extract_via_router(
        SAMPLE_PAPER, "KCNH2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )
    assert result["used_fallback"] is True
    assert result["variants"] == []
    assert result["error"] is None
