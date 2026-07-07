"""Tests for pipeline.table_router (table enumeration + router + parser)."""

import json
from pathlib import Path

import pytest

from pipeline.table_router import (
    MarkdownTable,
    _field_header_match,
    _infer_column_mapping_from_headers,
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


def test_enumerate_finds_etable_captions():
    text = """
eTable 1. LQT1 Mutations or Rare Variants

| cDNA | Protein | No. of patients |
| --- | --- | --- |
| c.521G>A | p.R147H | 2 |
"""

    tables = enumerate_markdown_tables(text)

    assert len(tables) == 1
    assert tables[0].caption == "eTable 1. LQT1 Mutations or Rare Variants"


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


def test_parse_routed_table_inherits_blank_gene_cell_for_filter():
    # Markdown rowspan: gene-grouped tables name the gene once, then leave
    # continuation rows blank. Blank cells must inherit the gene above so
    # off-target continuation rows are filtered (regression: PMID 33013630
    # Table 1, where HCN4 Val759Ile leaked into KCNH2 as 870 "carriers").
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Variants identified in long QT probands",
        header_line="| Gene | Protein | Proband | Phenotype |",
        header_cells=["Gene", "Protein", "Proband", "Phenotype"],
        data_lines=[
            "| KCNH2 | p.Ala561Val | P1 | LQT2 |",
            "|  | p.Arg176Trp | P2 | LQT2 |",
            "| SCN5A | p.Arg1623Gln | P3 | LQT3 |",
            "|  | p.Glu1784Lys | P4 | LQT3 |",
        ],
        char_start=0,
        char_end=240,
    )
    mapping = {"gene": 0, "protein": 1, "phenotype": 3}

    variants = parse_routed_table(table, mapping, "KCNH2")

    # Blank-cell continuation row p.Glu1784Lys inherits SCN5A and is dropped;
    # p.Arg176Trp inherits KCNH2 and is kept.
    assert [v["protein_notation"] for v in variants] == [
        "p.Ala561Val",
        "p.Arg176Trp",
    ]


def test_parse_routed_table_skips_population_annotation_table_without_subjects():
    # gnomAD/SIFT/PolyPhen annotation table with a clinical caption but no
    # patient/proband column must not mint one inferred carrier per row
    # (regression: PMID 33013630 Table 1).
    table = MarkdownTable(
        table_id="T1",
        caption=(
            "Table 1. Nonsynonymous variants in cardiac arrhythmia genes "
            "identified in sudden unexpected death in epilepsy (SUDEP)"
        ),
        header_line="| Gene | Variant | gnomAD allele count | SIFT | PolyPhen-2 |",
        header_cells=["Gene", "Variant", "gnomAD allele count", "SIFT", "PolyPhen-2"],
        data_lines=[
            "| KCNH2 | Ile82Thr | 0 | Deleterious | Benign |",
            "|  | Arg176Trp | 44 | Deleterious | Possibly damaging |",
        ],
        char_start=0,
        char_end=240,
    )
    mapping = {"gene": 0, "protein": 1, "patient_count": 2}

    # gnomAD column is rejected as a count, and the table has no subject column,
    # so no carriers are inferred -> zero variants from this annotation table.
    assert parse_routed_table(table, mapping, "KCNH2") == []


def test_parse_routed_table_skips_prediction_score_annotation_table():
    # In-silico predictor table (SIFT/PolyPhen/CADD/REVEL) with no gnomAD column
    # and no subject column must not mint one inferred carrier per row, even
    # though "score" is a clinical-context cue. Generalizes the population-table
    # guard to prediction-score annotation tables (extraction rule 3).
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. In silico pathogenicity predictions",
        header_line="| Gene | Variant | REVEL score | CADD score |",
        header_cells=["Gene", "Variant", "REVEL score", "CADD score"],
        data_lines=[
            "| KCNH2 | Ile82Thr | 0.7 | 24.1 |",
            "|  | Arg176Trp | 0.9 | 27.3 |",
        ],
        char_start=0,
        char_end=200,
    )

    assert parse_routed_table(table, {"gene": 0, "protein": 1}, "KCNH2") == []


def test_parse_routed_table_still_infers_clinical_list_without_annotation_columns():
    # Guard rail for the two guards above: a genuine one-proband-per-row clinical
    # list (mutation + phenotype, no population/prediction columns, no explicit
    # subject header) must still infer one carrier per row.
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Mutations identified in long QT probands",
        header_line="| Gene | Mutation | Phenotype | Age |",
        header_cells=["Gene", "Mutation", "Phenotype", "Age"],
        data_lines=[
            "| KCNH2 | p.Ala561Val | LQT2 | 34 |",
            "|  | p.Arg176Trp | LQT2 | 12 |",
        ],
        char_start=0,
        char_end=200,
    )

    variants = parse_routed_table(
        table, {"gene": 0, "protein": 1, "phenotype": 2}, "KCNH2"
    )

    assert [v["protein_notation"] for v in variants] == ["p.Ala561Val", "p.Arg176Trp"]
    assert all(v["penetrance_data"]["total_carriers_observed"] == 1 for v in variants)


def test_parse_routed_table_drops_off_target_gene_cell_even_when_gene_unmapped():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Large sequencing panel",
        header_line="| Gene | Protein | Carriers |",
        header_cells=["Gene", "Protein", "Carriers"],
        data_lines=[
            "| PIK3CA | p.Gly1049Arg | 40 |",
            "| KCNH2 | p.Arg176Trp | 2 |",
            "| RPS3A | p.Ala538Val | 12 |",
        ],
        char_start=0,
        char_end=220,
    )
    # Simulates a router response that forgot to map the Gene column.
    mapping = {"protein": 1, "patient_count": 2}

    variants = parse_routed_table(table, mapping, "KCNH2")

    assert [v["protein_notation"] for v in variants] == ["p.Arg176Trp"]


def test_parse_routed_table_accepts_target_gene_alias_in_gene_column():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. HERG carriers",
        header_line="| Gene | Protein | Carriers |",
        header_cells=["Gene", "Protein", "Carriers"],
        data_lines=["| HERG | p.Gly572Ser | 3 |"],
        char_start=0,
        char_end=120,
    )
    mapping = {"gene": 0, "protein": 1, "patient_count": 2}

    variants = parse_routed_table(table, mapping, "KCNH2")

    assert len(variants) == 1
    assert variants[0]["protein_notation"] == "p.Gly572Ser"


def test_parse_routed_table_accepts_default_phrase_alias_in_gene_column():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Cardiomyopathy carrier counts",
        header_line="| Gene | Protein | Carriers |",
        header_cells=["Gene", "Protein", "Carriers"],
        data_lines=["| cardiac myosin-binding protein C | p.Arg502Trp | 4 |"],
        char_start=0,
        char_end=140,
    )
    mapping = {"gene": 0, "protein": 1, "patient_count": 2}

    variants = parse_routed_table(table, mapping, "MYBPC3")

    assert len(variants) == 1
    assert variants[0]["protein_notation"] == "p.Arg502Trp"


def test_parse_routed_table_scopes_lqt_caption_without_gene_column():
    table = MarkdownTable(
        table_id="T1",
        caption="eTable 1. LQT1 Mutations or Rare Variants",
        header_line="| cDNA | Protein | No. of patients |",
        header_cells=["cDNA", "Protein", "No. of patients"],
        data_lines=["| c.521G>A | p.R147H | 2 |"],
        char_start=0,
        char_end=120,
    )
    mapping = {"cdna": 0, "protein": 1, "patient_count": 2}

    assert parse_routed_table(table, mapping, "SCN5A") == []
    variants = parse_routed_table(table, mapping, "KCNQ1")
    assert [v["protein_notation"] for v in variants] == ["p.R147H"]


def test_parse_routed_table_splits_parallel_pairs_and_sums_duplicate_rows():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Clinical and genetic summaries of probands.",
        header_line="| Patient | Nucleotide | Amino Acids | Symptom |",
        header_cells=["Patient", "Nucleotide", "Amino Acids", "Symptom"],
        data_lines=[
            "| 1 | exon 3 deletion | N57_G91del35 | Syncope |",
            "| 32 | 14311g>a | V4771I | Syncope |",
            "| 33 | 14311g>a | V4771I | CPA |",
            "| 35 | 14834_14835insTCA | 4944_4945insH | CPA |",
            "| 36 | 9910c>g, 14222c>t | Q3304E, A4741V | CPA |",
        ],
        char_start=0,
        char_end=300,
    )
    mapping = {"cdna": 1, "protein": 2, "patient_count": -1, "phenotype": 3}

    variants = parse_routed_table(table, mapping, "RYR2")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {
        "N57_G91del35",
        "V4771I",
        "4944_4945insH",
        "Q3304E",
        "A4741V",
    }
    assert by_protein["N57_G91del35"]["cdna_notation"] is None
    assert by_protein["V4771I"]["penetrance_data"]["total_carriers_observed"] == 2
    assert by_protein["4944_4945insH"]["cdna_notation"] == "c.14834_14835insTCA"
    assert by_protein["Q3304E"]["cdna_notation"] == "c.9910c>g"
    assert by_protein["A4741V"]["cdna_notation"] == "c.14222c>t"


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


def test_extract_via_router_retries_empty_llm_response_then_reports_error():
    calls = 0

    def stub(**_):
        nonlocal calls
        calls += 1
        return _stub_response("")

    paper = """Table 1. Functional assays

| construct | tail current | activation V1/2 |
|-----------|--------------|-----------------|
| WT | 100% | -32 |
| K897T | 95% | -30 |
"""

    result = extract_via_router(
        paper, "KCNH2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )

    assert calls == 2
    assert result["used_fallback"] is True
    assert result["variants"] == []
    assert "empty LLM response" in result["error"]


# --- B3: cohort-label guard (strict_cohort_labels) -------------------------


def test_field_header_match_guards_ambiguous_cohort_labels_on_assay_tables():
    # Default (flag off): any keyword match maps, as before.
    assert _field_header_match(
        "controls", "unaffected", strict_cohort_labels=False, has_assay_or_gwas_cue=True
    )
    # Ambiguous-only term on an assay/case-control table -> rejected under strict.
    assert not _field_header_match(
        "controls", "unaffected", strict_cohort_labels=True, has_assay_or_gwas_cue=True
    )
    assert not _field_header_match(
        "cases", "affected", strict_cohort_labels=True, has_assay_or_gwas_cue=True
    )
    # Ambiguous term but NOT an assay table -> still mapped (context allows it).
    assert _field_header_match(
        "controls", "unaffected", strict_cohort_labels=True, has_assay_or_gwas_cue=False
    )
    # Unambiguous per-variant terms are always mapped, even strict + assay.
    assert _field_header_match(
        "unaffected",
        "unaffected",
        strict_cohort_labels=True,
        has_assay_or_gwas_cue=True,
    )
    assert _field_header_match(
        "affected", "affected", strict_cohort_labels=True, has_assay_or_gwas_cue=True
    )
    # No keyword at all -> never a match.
    assert not _field_header_match(
        "age", "affected", strict_cohort_labels=True, has_assay_or_gwas_cue=True
    )


def _case_control_table() -> MarkdownTable:
    # 'rsID' makes has_assay_or_gwas_cue True; protein col keeps it from being
    # rejected as a pure-assay table; 'Carriers' is an unambiguous count column.
    return MarkdownTable(
        table_id="T1",
        caption="Table 2. Case-control association",
        header_line="| rsID | Protein | Carriers | Cases | Controls |",
        header_cells=["rsID", "Protein", "Carriers", "Cases", "Controls"],
        data_lines=["| rs1 | p.Arg190Trp | 3 | 5 | 200 |"],
        char_start=0,
        char_end=80,
    )


def test_infer_mapping_default_skips_cohort_columns_when_carrier_present():
    mapping = _infer_column_mapping_from_headers(_case_control_table())
    assert mapping is not None
    assert "affected" not in mapping
    assert "unaffected" not in mapping
    assert mapping.get("patient_count") == 2


def test_infer_mapping_strict_skips_cohort_columns_on_assay_table():
    mapping = _infer_column_mapping_from_headers(
        _case_control_table(), strict_cohort_labels=True
    )
    assert mapping is not None
    assert "affected" not in mapping
    assert "unaffected" not in mapping
    # The unambiguous per-variant count column is still mapped.
    assert "patient_count" in mapping


def test_infer_mapping_strict_keeps_unambiguous_labels():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 3. KCNH2 variant carriers",
        header_line="| Protein | Affected | Unaffected |",
        header_cells=["Protein", "Affected", "Unaffected"],
        data_lines=["| K897T | 2 | 3 |"],
        char_start=0,
        char_end=80,
    )
    mapping = _infer_column_mapping_from_headers(table, strict_cohort_labels=True)
    assert mapping is not None
    # 'Affected'/'Unaffected' are unambiguous -> still mapped even on an assay table.
    assert mapping.get("affected") == 1
    assert mapping.get("unaffected") == 2

    variants = parse_routed_table(table, mapping, "KCNH2")
    assert (
        variants[0]["count_provenance"]["affected_count_type"] == "per_variant_carrier"
    )
    assert (
        variants[0]["count_provenance"]["unaffected_count_type"] == "unaffected_control"
    )


def test_infer_mapping_prefers_carrier_over_total_case_denominator():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. BRCA1 variants in cases",
        header_line="| cDNA | Protein | Total case | Carrier |",
        header_cells=["cDNA", "Protein", "Total case", "Carrier"],
        data_lines=["| c.5467del | p.Ile1824fs | 1505 | 20 |"],
        char_start=0,
        char_end=100,
    )

    mapping = _infer_column_mapping_from_headers(table)

    assert mapping is not None
    assert mapping.get("patient_count") == 3
    assert "affected" not in mapping

    variants = parse_routed_table(table, mapping, "BRCA1")

    assert len(variants) == 1
    pen = variants[0]["penetrance_data"]
    assert pen["total_carriers_observed"] == 20
    assert pen["affected_count"] == 20


def test_parse_routed_table_treats_adult_number_as_row_identifier():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. MYBPC3 adult HCM patients",
        header_line="| Adult number | Diagnosis | Protein |",
        header_cells=["Adult number", "Diagnosis", "Protein"],
        data_lines=["| 172 | HCM | p.Trp916Ter |"],
        char_start=0,
        char_end=100,
    )
    mapping = {"patient_count": 0, "phenotype": 1, "protein": 2}

    variants = parse_routed_table(table, mapping, "MYBPC3")

    assert len(variants) == 1
    assert variants[0]["penetrance_data"]["total_carriers_observed"] == 1
    assert variants[0]["penetrance_data"]["affected_count"] == 1
    assert (
        variants[0]["count_provenance"]["carriers_column_label"]
        == "implicit one carrier per clinical row"
    )


def test_parse_routed_table_rejects_population_occurrence_counts():
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. 1000 Genomes MYBPC3 allele frequencies",
        header_line="| Gene | Protein | MAF (No. of occurrences from n=2184 alleles) |",
        header_cells=[
            "Gene",
            "Protein",
            "MAF (No. of occurrences from n=2184 alleles)",
        ],
        data_lines=["| MYBPC3 | p.Ser236Gly | 0.09 (189) |"],
        char_start=0,
        char_end=120,
    )
    mapping = {"gene": 0, "protein": 1, "patient_count": 2}

    assert parse_routed_table(table, mapping, "MYBPC3") == []
