"""Tests for pipeline.table_router (table enumeration + router + parser)."""

import json
from pathlib import Path

import pytest

from pipeline.table_router import (
    MarkdownTable,
    _clinical_frequency_denominator,
    _count_from_frequency,
    _field_header_match,
    _gene_symbol_tokens,
    _infer_column_mapping_from_headers,
    _normalize_header,
    _normalize_protein,
    _row_has_off_target_gene_without_target,
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


def test_parse_routed_table_keeps_footnoted_legacy_indels():
    """A terminal star is a table footnote, not part of either allele."""
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. KCNQ1 mutations identified in patients",
        header_line="| cDNA | Protein | No. of patients |",
        header_cells=["cDNA", "Protein", "No. of patients"],
        data_lines=[
            "| c.1071_1076dupGAAGCA* | p.360_361dupKQ* | 1 |",
            "| c.1149insT* | p.A384fs79X* | 1 |",
        ],
        char_start=0,
        char_end=160,
    )

    variants = parse_routed_table(
        table, {"cdna": 0, "protein": 1, "patient_count": 2}, "KCNQ1"
    )

    assert [
        (
            v["cdna_notation"],
            v["protein_notation"],
            v["penetrance_data"]["total_carriers_observed"],
        )
        for v in variants
    ] == [
        ("c.1071_1076dupGAAGCA", "p.360_361dupKQ", 1),
        ("c.1149insT", "p.A384fs79X", 1),
    ]


def test_reference_less_protein_dup_rejects_cdna_shaped_payload():
    assert _normalize_protein("p.360_361dupKQ†") == "p.360_361dupKQ"
    assert _normalize_protein("360_361dupKQ") is None
    assert _normalize_protein("p.1071_1076dupGAAGCA") is None
    assert _normalize_protein("1071_1076dupGAAGCA") is None


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


def _one_carrier_row(phenotype: str, functional: str) -> dict:
    """Parse a single clinical row and return its penetrance_data."""
    header = "| Gene | Protein | Phenotype | Functional effect |"
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Mutations identified in long QT probands",
        header_line=header,
        header_cells=[c.strip() for c in header.strip("|").split("|")],
        data_lines=[f"| KCNH2 | p.Ala561Val | {phenotype} | {functional} |"],
        char_start=0,
        char_end=400,
    )
    variants = parse_routed_table(
        table, {"gene": 0, "protein": 1, "phenotype": 2}, "KCNH2"
    )
    assert len(variants) == 1
    return variants[0]["penetrance_data"]


@pytest.mark.parametrize(
    "phenotype,functional,affected,unaffected",
    [
        # The bug: an unrelated functional-effect cell decided the phenotype.
        # "normal"/"control"/"healthy" are ordinary assay English, and a
        # row-wide search on them recorded a symptomatic proband as a confirmed
        # unaffected carrier — fabricating the non-penetrance signal this
        # database exists to measure.
        ("LQT2, syncope", "normal trafficking", 1, None),
        ("LQT2, aborted cardiac arrest", "control-like kinetics", 1, None),
        ("LQT2, syncope", "reduced current", 1, None),
        # Genotype-positive/phenotype-negative must survive: an unambiguous
        # status term outranks a disease token in the same cell.
        ("asymptomatic", "reduced current", 0, 1),
        ("asymptomatic LQT2", "reduced current", 0, 1),
        ("unaffected carrier", "reduced current", 0, 1),
        # Ambiguous words are believed as the whole clinical cell, nowhere else.
        ("normal", "reduced current", 0, 1),
        ("healthy", "reduced current", 0, 1),
        # A denied finding is neither: naming it affected because the cell also
        # mentions a condition would invent the opposite partition.
        ("normal ECG, no syncope", "reduced current", None, None),
        ("LQT2 without arrhythmic events", "reduced current", None, None),
        # ...but an ordinal is not a denial.
        ("Patient No. 3, LQT2", "reduced current", 1, None),
        # Silence stays silent.
        ("unknown", "normal trafficking", None, None),
        ("", "normal trafficking", None, None),
    ],
)
def test_routed_clinical_row_reads_phenotype_only(
    phenotype, functional, affected, unaffected
):
    pdata = _one_carrier_row(phenotype, functional)
    assert pdata["total_carriers_observed"] == 1
    assert pdata.get("affected_count") == affected
    assert pdata.get("unaffected_count") == unaffected


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
    assert affected["unaffected_count"] is None
    assert unaffected["total_carriers_observed"] == 1
    assert unaffected["affected_count"] == 0
    assert unaffected["unaffected_count"] == 1


def test_extract_via_router_deduplicates_repeated_complete_table():
    def stub(**_):
        raise AssertionError("LLM router should not run for clinical tables")

    paper = """Table 1. BMPR2 variants in patients

| Gene | cDNA | Number of patients |
|---|---|---|
| BMPR2 | c.76+1G>T | 2 |

Table 2. Duplicate archive rendering

| Gene | cDNA | Number of patients |
|---|---|---|
| BMPR2 | c.76+1G>T | 2 |
"""

    result = extract_via_router(
        paper, "BMPR2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )

    assert len(result["variants"]) == 1
    assert result["variants"][0]["cdna_notation"] == "c.76+1G>T"
    assert result["variants"][0]["patients"]["count"] == 2
    assert result["duplicate_routed_table_ids"] == ["T2"]


def test_extract_via_router_dedup_does_not_let_off_target_caption_win():
    def stub(**_):
        raise AssertionError("LLM router should not run for deterministic tables")

    paper = """Table 1. BRCA1 variants in patients

| cDNA | Number of patients |
|---|---|
| c.68_69del | 2 |

Table 2. BRCA2 variants in patients

| cDNA | Number of patients |
|---|---|
| c.68_69del | 2 |
"""

    result = extract_via_router(
        paper, "BRCA2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )

    assert [v["cdna_notation"] for v in result["variants"]] == ["c.68_69del"]
    assert result["variants"][0]["source_ref"] == "Table 2. BRCA2 variants in patients"
    assert result["duplicate_routed_table_ids"] == []


def test_extract_via_router_keeps_distinct_tables_sharing_one_variant():
    def stub(**_):
        raise AssertionError("LLM router should not run for deterministic tables")

    paper = """Table 1. BMPR2 variants in patients

| Gene | cDNA | Number of patients |
|---|---|---|
| BMPR2 | c.76+1G>T | 2 |

Table 2. BMPR2 variants in replication patients

| Gene | cDNA | Number of patients |
|---|---|---|
| BMPR2 | c.76+1G>T | 3 |
| BMPR2 | c.350G>A | 1 |
"""

    result = extract_via_router(
        paper, "BMPR2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )

    assert len(result["variants"]) == 3
    assert [v["patients"]["count"] for v in result["variants"]] == [2, 3, 1]
    assert result["duplicate_routed_table_ids"] == []


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
    assert pen["affected_count"] is None
    assert variants[0]["count_provenance"]["affected_column_label"] is None
    assert variants[0]["count_provenance"]["affected_count_type"] is None


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


def test_router_infers_blank_gene_group_and_ignores_family_history_case_counts():
    """The routed and fast paths must enforce the same structural role rules."""
    table = MarkdownTable(
        table_id="T5",
        caption="Table 5. Family and pathological characteristics",
        header_line=(
            "|  | Nucleotide change | AA change | Age of diagnosis | "
            "No. of breast/ovarian cancers in family | Mean age diagnosis | "
            "No. of cases (Breast cancer <= 50) |"
        ),
        header_cells=[
            "",
            "Nucleotide change",
            "AA change",
            "Age of diagnosis",
            "No. of breast/ovarian cancers in family",
            "Mean age diagnosis",
            "No. of cases (Breast cancer <= 50)",
        ],
        data_lines=[
            "| BRCA1 | 180 delA | STOP 22 | 55 | 2 | 43 | 1 |",
            "|  | 185 delAG | STOP 39 | 33 | 1 | 33 | 1 |",
            "| BRCA2 | 490 delCT | STOP 99 | 39 | 2 | 56 | 1 |",
            "|  | 1184 insA | STOP 326 | 34 | 1 | 34 | 1 |",
        ],
        char_start=0,
        char_end=400,
    )

    mapping = _infer_column_mapping_from_headers(table)

    assert mapping is not None
    assert mapping["gene"] == 0
    assert mapping["patient_count"] == -1
    assert "affected" not in mapping
    variants = parse_routed_table(table, mapping, "BRCA2")
    assert [v["cdna_notation"] for v in variants] == [None, None]
    assert [v["legacy_notation"] for v in variants] == [
        "490delCT",
        "1184insA",
    ]
    assert [v["source_notation"] for v in variants] == ["490 delCT", "1184 insA"]
    assert all(v["patients"]["count"] == 1 for v in variants)
    assert all(
        v["count_provenance"]["carriers_column_label"]
        == "implicit one carrier per clinical row"
        for v in variants
    )


def test_router_never_invents_cdna_prefix_for_bic_deleted_base_count():
    table = MarkdownTable(
        table_id="T6",
        caption="Table 6. BRCA2 variants identified in study probands",
        header_line="| Nucleotide change | Proband count |",
        header_cells=["Nucleotide change", "Proband count"],
        data_lines=["| 4184del4 | 7 |", "| c.9117G>A | 2 |"],
        char_start=0,
        char_end=120,
    )
    mapping = {"cdna": 0, "patient_count": 1}

    variants = parse_routed_table(table, mapping, "BRCA2")

    assert len(variants) == 2
    legacy = next(row for row in variants if row["legacy_notation"])
    assert legacy["cdna_notation"] is None
    assert legacy["legacy_notation"] == "4184del4"
    assert legacy["source_notation"] == "4184del4"
    assert legacy["patients"]["count"] == 7
    assert any(
        fact["fact_type"] == "variant_identity" and fact["fact_value"] == "4184del4"
        for fact in legacy["fact_provenance"]
    )

    explicit = next(row for row in variants if row["cdna_notation"])
    assert explicit["cdna_notation"] == "c.9117G>A"
    assert explicit["legacy_notation"] is None


def test_router_treats_non_brca_prefixless_indel_as_omitted_prefix_cdna():
    table = MarkdownTable(
        table_id="T7",
        caption="BMPR2 variants identified in PAH probands",
        header_line="| Nucleotide change | Proband count |",
        header_cells=["Nucleotide change", "Proband count"],
        data_lines=["| 1234delA | 3 |"],
        char_start=0,
        char_end=80,
    )

    (variant,) = parse_routed_table(table, {"cdna": 0, "patient_count": 1}, "BMPR2")

    assert variant["cdna_notation"] == "c.1234delA"
    assert variant["legacy_notation"] is None


@pytest.mark.parametrize(
    "timing_header,family_header,ambiguous_count_header",
    [
        ("Age at onset", "Affected relatives", "Cases before age 40"),
        (
            "Duration of symptoms",
            "Family members with disease",
            "Clinical events",
        ),
        ("Follow-up (years)", "Family history", "Cancer diagnoses"),
    ],
)
def test_router_classifies_clinical_measures_by_role_not_exact_header(
    timing_header: str,
    family_header: str,
    ambiguous_count_header: str,
):
    """Header families, not one paper's wording, keep phenotype values out."""
    table = MarkdownTable(
        table_id="T-clinical",
        caption="Clinical characteristics by variant",
        header_line="",
        header_cells=[
            "Gene",
            "cDNA",
            "Protein",
            timing_header,
            family_header,
            ambiguous_count_header,
        ],
        data_lines=["| PALB2 | c.3113G>A | p.Gly1038Asp | 42 | 3 | 2 |"],
        char_start=0,
        char_end=160,
    )

    mapping = _infer_column_mapping_from_headers(table, target_gene="PALB2")

    assert mapping is not None
    assert mapping["patient_count"] == -1
    assert "affected" not in mapping
    assert "unaffected" not in mapping


def test_router_infers_open_vocabulary_unnamed_gene_groups():
    """Unknown panel genes are scoped from rowspan structure, not a registry."""
    table = MarkdownTable(
        table_id="T-open",
        caption="Grouped hereditary-disease variants",
        header_line="|  | cDNA | Protein | No. of patients |",
        header_cells=["", "cDNA", "Protein", "No. of patients"],
        data_lines=[
            "| PALB2 | c.3113G>A | p.Gly1038Asp | 4 |",
            "|  | c.2257C>T | p.Arg753Ter | 2 |",
            "| ATM | c.7271T>G | p.Val2424Gly | 3 |",
            "|  | c.8147T>C | p.Val2716Ala | 5 |",
        ],
        char_start=0,
        char_end=240,
    )

    mapping = _infer_column_mapping_from_headers(table, target_gene="PALB2")

    assert mapping is not None
    assert mapping["gene"] == 0
    variants = parse_routed_table(table, mapping, "PALB2")
    assert [variant["cdna_notation"] for variant in variants] == [
        "c.3113G>A",
        "c.2257C>T",
    ]


def test_router_does_not_treat_unlabeled_row_ids_as_open_gene_groups():
    """Open-vocabulary inference requires rowspan blanks, not merely ID tokens."""
    table = MarkdownTable(
        table_id="T-ids",
        caption="Variant observations",
        header_line="|  | cDNA | No. of patients |",
        header_cells=["", "cDNA", "No. of patients"],
        data_lines=[
            "| FAM101 | c.3113G>A | 4 |",
            "| FAM102 | c.2257C>T | 2 |",
        ],
        char_start=0,
        char_end=120,
    )

    mapping = _infer_column_mapping_from_headers(table, target_gene="PALB2")

    assert mapping is not None
    assert "gene" not in mapping


# ---------------------------------------------------------------------------
# Router decision observability
#
# The router's column->field mapping is applied deterministically to every row
# of a routed table, so a mapping error scales to hundreds of rows, and a
# decline silently hands the paper to the full-text path. These tests pin the
# decision record onto extraction_metadata so both are auditable after a run.
# ---------------------------------------------------------------------------


FUNCTIONAL_ASSAY_PAPER = """Table 1. Functional assays (no variant data)

| construct | tail current | activation V1/2 |
|-----------|--------------|-----------------|
| WT | 100% | -32 |
| K897T | 95% | -30 |
"""

NARRATIVE_PAPER = (
    "We followed a KCNH2 mutation carrier family with long QT syndrome. "
    "The proband's variant was identified by sequencing. " * 20
)

# extract() has a 500-char circuit breaker, so pad the two-table sample with
# prose (no pipes, so the enumerated table count stays at 2).
SAMPLE_PAPER_PADDED = SAMPLE_PAPER + "\n\n" + NARRATIVE_PAPER


def _extractor():
    from pipeline.extraction import ExpertExtractor

    return ExpertExtractor(models=["test-model"], tier_threshold=1)


def _paper(pmid="12345678", full_text=NARRATIVE_PAPER, gene_symbol="KCNH2"):
    from utils.models import Paper

    return Paper(
        pmid=pmid,
        title="Router observability fixture",
        gene_symbol=gene_symbol,
        full_text=full_text,
    )


def _stub_router_llm(monkeypatch, content):
    """Point the router's default llm_caller at a canned response."""
    calls = []

    def _fake(**kwargs):
        calls.append(kwargs)
        return _stub_response(content)

    monkeypatch.setattr("utils.llm_utils.litellm_completion", _fake)
    return calls


def test_recorded_outcome_codes_are_documented():
    """Pin the wire values downstream scoring will filter on."""
    from pipeline.table_router import ROUTER_OUTCOMES

    for code in (
        "approved",
        "approved_below_threshold",
        "no_usable_tables",
        "llm_error",
        "crashed",
        "disabled",
        "skipped",
    ):
        assert code in ROUTER_OUTCOMES


def test_build_router_decision_metadata_records_mapping_and_declines():
    from pipeline.table_router import RoutedTable, build_router_decision_metadata

    routed = [
        RoutedTable(
            table_id="T1",
            column_mapping={"cdna": 0, "protein": 1, "patient_count": 2},
            confidence=0.9,
            notes="carrier counts in column 2",
        )
    ]

    meta = build_router_decision_metadata(
        model="azure_ai/Kimi-K2.6-1",
        outcome="approved",
        routed=routed,
        enumerated_table_ids=["T1", "T2", "T3"],
        variants_parsed=17,
    )

    assert meta["table_router_model"] == "azure_ai/Kimi-K2.6-1"
    assert meta["table_router_outcome"] == "approved"
    assert meta["table_router_tables_enumerated"] == 3
    assert meta["table_router_tables_routed"] == 1
    assert meta["table_router_variants_parsed"] == 17
    assert meta["table_router_mappings"] == {
        "T1": {"cdna": 0, "protein": 1, "patient_count": 2}
    }
    assert meta["table_router_confidences"] == {"T1": 0.9}
    assert meta["table_router_notes"] == {"T1": "carrier counts in column 2"}
    # The tables the router looked at and passed over are the audit target.
    assert meta["table_router_tables_declined"] == ["T2", "T3"]
    assert "table_router_error" not in meta
    assert "table_router_detail_truncated" not in meta
    # Must survive the extraction-JSON round trip.
    assert json.loads(json.dumps(meta)) == meta


def test_build_router_decision_metadata_caps_per_table_detail():
    from pipeline.table_router import (
        MAX_ROUTER_DETAIL_TABLES,
        RoutedTable,
        build_router_decision_metadata,
    )

    n = MAX_ROUTER_DETAIL_TABLES + 8
    routed = [
        RoutedTable(
            table_id=f"T{i}",
            column_mapping={"cdna": 0},
            confidence=0.5,
            notes="x" * 500,
        )
        for i in range(n)
    ]

    meta = build_router_decision_metadata(
        model="m",
        outcome="approved",
        routed=routed,
        enumerated_table_ids=[f"T{i}" for i in range(n)],
        variants_parsed=3,
    )

    # Counts stay exact; only the per-table detail is capped.
    assert meta["table_router_tables_routed"] == n
    assert meta["table_router_tables_enumerated"] == n
    assert len(meta["table_router_mappings"]) == MAX_ROUTER_DETAIL_TABLES
    assert len(meta["table_router_confidences"]) == MAX_ROUTER_DETAIL_TABLES
    assert all(len(v) <= 200 for v in meta["table_router_notes"].values())
    assert meta["table_router_detail_truncated"] is True


def test_build_router_decision_metadata_records_error():
    from pipeline.table_router import build_router_decision_metadata

    meta = build_router_decision_metadata(
        model="m", outcome="llm_error", error="boom " * 200
    )

    assert meta["table_router_outcome"] == "llm_error"
    assert meta["table_router_tables_enumerated"] == 0
    assert meta["table_router_mappings"] == {}
    assert len(meta["table_router_error"]) <= 300


def test_extract_via_router_reports_enumerated_table_ids():
    def stub(**_):
        return _stub_response(json.dumps({"variant_tables": []}))

    result = extract_via_router(
        SAMPLE_PAPER, "KCNH2", model="azure_ai/Kimi-K2.6-1", llm_caller=stub
    )

    expected = [t.table_id for t in enumerate_markdown_tables(SAMPLE_PAPER)]
    assert result["enumerated_table_ids"] == expected
    assert len(expected) == 2


def test_router_metadata_recorded_on_approved_path(monkeypatch):
    """Approved path persists the model, mapping, confidence and declines."""
    _stub_router_llm(monkeypatch, json.dumps({"variant_tables": []}))
    extractor = _extractor()
    paper = _paper(full_text=SAMPLE_PAPER)

    result = extractor._try_table_router(paper, SAMPLE_PAPER)

    assert result is not None and result.success
    meta = result.extracted_data["extraction_metadata"]
    assert meta["table_router_outcome"] == "approved"
    assert (
        meta["table_router_model"]
        == extractor.last_table_router_decision["table_router_model"]
    )
    assert meta["table_router_tables_enumerated"] == 2
    assert meta["table_router_tables_routed"] == 1
    assert meta["table_router_variants_parsed"] == 2
    # The mapping that was applied to every row of the routed table.
    (mapping,) = meta["table_router_mappings"].values()
    assert mapping["cdna"] == 0
    assert mapping["patient_count"] == 2
    # The functional-assay table was enumerated but not routed.
    assert len(meta["table_router_tables_declined"]) == 1
    # Pre-existing metadata is untouched.
    assert meta["total_variants_found"] == 2
    assert extractor.last_table_router_decision["table_router_outcome"] == "approved"


def test_router_metadata_recorded_on_decline_path(monkeypatch):
    """A decline records what was enumerated, so it can be re-examined later."""
    _stub_router_llm(monkeypatch, json.dumps({"variant_tables": []}))
    extractor = _extractor()
    paper = _paper(full_text=FUNCTIONAL_ASSAY_PAPER)

    result = extractor._try_table_router(paper, FUNCTIONAL_ASSAY_PAPER)

    # Behaviour unchanged: the caller still falls through to full-text Tier 3.
    assert result is None
    decision = extractor.last_table_router_decision
    assert decision["table_router_outcome"] == "no_usable_tables"
    assert decision["table_router_tables_enumerated"] == 1
    assert decision["table_router_tables_routed"] == 0
    assert decision["table_router_variants_parsed"] == 0
    assert decision["table_router_mappings"] == {}
    assert len(decision["table_router_tables_declined"]) == 1


def test_router_metadata_records_crash(monkeypatch):
    import pipeline.table_router as table_router

    def _boom(*_args, **_kwargs):
        raise RuntimeError("router exploded")

    monkeypatch.setattr(table_router, "extract_via_router", _boom)
    extractor = _extractor()

    assert extractor._try_table_router(_paper(), NARRATIVE_PAPER) is None
    decision = extractor.last_table_router_decision
    assert decision["table_router_outcome"] == "crashed"
    assert "router exploded" in decision["table_router_error"]


def test_router_metadata_records_llm_error(monkeypatch):
    """An empty router response is recorded as llm_error, not a clean decline."""
    _stub_router_llm(monkeypatch, "")
    extractor = _extractor()

    assert extractor._try_table_router(_paper(), SAMPLE_PAPER) is None
    decision = extractor.last_table_router_decision
    assert decision["table_router_outcome"] == "llm_error"
    assert decision["table_router_error"]
    assert decision["table_router_tables_enumerated"] == 2


def _stub_full_text_extraction(monkeypatch, extractor):
    """Make the full-text Tier 3 path offline and deterministic."""

    class EmptyScanner:
        variants = []

        def get_hints_for_prompt(self, max_hints):
            return ""

    monkeypatch.setattr(extractor, "_extract_variants_from_tables", lambda *_: [])
    monkeypatch.setattr(
        "pipeline.extraction.scan_document_for_variants",
        lambda *_, **__: EmptyScanner(),
    )
    monkeypatch.setattr(
        extractor,
        "call_llm_json_with_status",
        lambda _prompt, **_kw: (
            {
                "variants": [
                    {"gene_symbol": "KCNH2", "protein_notation": "p.Arg176Trp"}
                ],
                "extraction_metadata": {"total_variants_found": 1},
            },
            False,
            "{}",
        ),
    )


def test_decline_metadata_lands_in_persisted_extraction_metadata(monkeypatch):
    """The 51k-case shape: nothing to route, answer comes from full text."""
    extractor = _extractor()
    _stub_full_text_extraction(monkeypatch, extractor)

    result = extractor.extract(_paper())

    assert result.success
    meta = result.extracted_data["extraction_metadata"]
    assert meta["table_router_outcome"] == "no_usable_tables"
    assert meta["table_router_tables_enumerated"] == 0
    assert meta["table_router_tables_routed"] == 0
    assert "table_router_tables_declined" not in meta


def test_router_disabled_records_disabled_outcome(monkeypatch):
    monkeypatch.setenv("ENABLE_TABLE_ROUTER", "false")
    extractor = _extractor()
    _stub_full_text_extraction(monkeypatch, extractor)

    def _never(*_args, **_kwargs):
        raise AssertionError("router must not run when disabled")

    monkeypatch.setattr(extractor, "_try_table_router", _never)

    result = extractor.extract(_paper())

    assert result.success
    meta = result.extracted_data["extraction_metadata"]
    assert meta["table_router_outcome"] == "disabled"
    assert meta["table_router_tables_enumerated"] == 0


def test_router_skipped_when_gene_symbol_missing(monkeypatch):
    extractor = _extractor()
    _stub_full_text_extraction(monkeypatch, extractor)

    def _never(*_args, **_kwargs):
        raise AssertionError("router needs a gene symbol")

    monkeypatch.setattr(extractor, "_try_table_router", _never)

    result = extractor.extract(_paper(gene_symbol=None))

    assert result.success
    assert (
        result.extracted_data["extraction_metadata"]["table_router_outcome"]
        == "skipped"
    )


def test_approved_metadata_lands_in_persisted_extraction_metadata(monkeypatch):
    """Router output becomes the answer: mapping is persisted with it."""
    from pipeline.extraction import ExpertExtractor

    _stub_router_llm(monkeypatch, json.dumps({"variant_tables": []}))
    # tier_threshold=0 lowers the short-circuit bar to one variant.
    extractor = ExpertExtractor(models=["test-model"], tier_threshold=0)

    result = extractor.extract(_paper(full_text=SAMPLE_PAPER_PADDED))

    assert result.success
    assert result.model_used.startswith("router+")
    meta = result.extracted_data["extraction_metadata"]
    assert meta["table_router_outcome"] == "approved"
    assert meta["table_router_tables_routed"] == 1
    assert list(meta["table_router_mappings"].values())[0]["patient_count"] == 2


def test_approved_below_threshold_recorded_when_output_demoted_to_hints(monkeypatch):
    """Router approved tables but too few variants to skip the full-text path."""
    _stub_router_llm(monkeypatch, json.dumps({"variant_tables": []}))
    extractor = _extractor()
    _stub_full_text_extraction(monkeypatch, extractor)

    result = extractor.extract(_paper(full_text=SAMPLE_PAPER_PADDED))

    assert result.success
    assert not (result.model_used or "").startswith("router+")
    meta = result.extracted_data["extraction_metadata"]
    # "approved" alone would overstate what the run actually used.
    assert meta["table_router_outcome"] == "approved_below_threshold"
    assert meta["table_router_tables_routed"] == 1
    assert meta["table_router_variants_parsed"] == 2


def test_router_decision_does_not_leak_between_papers(monkeypatch):
    """A fresh attempt clears the previous paper's decision."""
    _stub_router_llm(monkeypatch, json.dumps({"variant_tables": []}))
    extractor = _extractor()

    approved = extractor._try_table_router(_paper(full_text=SAMPLE_PAPER), SAMPLE_PAPER)
    assert approved is not None
    assert extractor.last_table_router_decision["table_router_outcome"] == "approved"

    _stub_full_text_extraction(monkeypatch, extractor)
    result = extractor.extract(_paper(pmid="87654321"))

    meta = result.extracted_data["extraction_metadata"]
    assert meta["table_router_outcome"] == "no_usable_tables"
    assert meta["table_router_mappings"] == {}


# ---------------------------------------------------------------------------
# Separator-less ("borderless") pipe tables with wrapped headers.
#
# Office-format supplements (.doc/.docx) converted to text emit pipe-delimited
# rows with no `|---|` separator row anywhere in the document, and wrap long
# header cells across physical lines. Before the borderless branch existed,
# `enumerate_markdown_tables` required a separator, so every such table was
# invisible and the whole paper scored "no usable variant tables" — PMID
# 26496715 carries 100 count-bearing cardiac gold rows in exactly this shape.
# ---------------------------------------------------------------------------

BORDERLESS_WRAPPED = """Table 1.  Summary of putative LQT1-associated mutations in KCNQ1
|Region  |Nucleotide  |Variant          |Mutation Type   |Location  |No. of |
|        |            |                 |                |          |patient|
|        |            |                 |                |          |s      |
|Exon 1  |c.2T>C      |p.Met1Thr        |Missense        |N-Terminal|3      |
|Exon 2  |c.420C>A    |p.Ser140Arg      |Missense        |S1 Domain |1      |
|Exon 3  |c.551A>G    |p.Tyr184Cys      |Missense        |S2-S3     |2      |
"""


def test_borderless_table_is_enumerated_without_a_separator_row():
    tables = enumerate_markdown_tables(BORDERLESS_WRAPPED)
    assert len(tables) == 1, "a separator-less pipe block must still be a table"
    table = tables[0]
    assert table.caption is not None and table.caption.startswith("Table 1.")
    assert len(table.data_lines) == 3


def test_borderless_wrapped_header_is_rejoined_into_one_cell():
    table = enumerate_markdown_tables(BORDERLESS_WRAPPED)[0]
    # "No. of" + "patient" + "s" must land in the SAME column, not become
    # extra data rows. The join is space-free; _normalize_header keeps only
    # alphanumerics, so this is what the count-column matcher sees.
    assert len(table.header_cells) == 6
    assert _normalize_header(table.header_cells[-1]) == "noofpatients"


def test_borderless_wrapped_header_rows_are_not_treated_as_data():
    table = enumerate_markdown_tables(BORDERLESS_WRAPPED)[0]
    # The two continuation lines are mostly-blank; if they leaked into
    # data_lines the parser would emit two variant-less rows.
    for row in table.data_lines:
        assert "Exon" in row


def test_borderless_wrapped_header_yields_a_patient_count_mapping():
    """The point of the fix: the count column becomes bindable."""
    table = enumerate_markdown_tables(BORDERLESS_WRAPPED)[0]
    mapping = _infer_column_mapping_from_headers(table)
    assert mapping is not None, "notation + count columns should both be found"
    assert "patient_count" in mapping
    assert mapping["patient_count"] == 5
    assert mapping.get("protein") == 2


def test_borderless_parse_binds_counts_to_variants():
    table = enumerate_markdown_tables(BORDERLESS_WRAPPED)[0]
    mapping = _infer_column_mapping_from_headers(table)
    variants = parse_routed_table(table, mapping, gene_symbol="KCNQ1")
    counts = {
        v.get("protein_notation"): (v.get("patients") or {}).get("count")
        for v in variants
    }
    assert len(variants) == 3, f"expected 3 variants, got {counts}"
    assert counts == {"p.Met1Thr": 3, "p.Ser140Arg": 1, "p.Tyr184Cys": 2}
    # The rebuilt header must reach count provenance, and the count must be
    # typed as a per-variant carrier count rather than a cohort total.
    prov = variants[0]["count_provenance"]
    assert prov["carriers_count_type"] == "per_variant_carrier"
    assert "patients" in prov["carriers_column_label"].replace(" ", "").lower()
    assert variants[0]["penetrance_data"]["total_carriers_observed"] == 3


def test_structural_annotation_cells_are_not_read_as_other_genes():
    """`Missense` / `N-Terminal` must not read as an off-target gene symbol.

    A mutation-type or protein-domain column with no mapped Gene column used to
    make `_row_has_off_target_gene_without_target` fire on every row, silently
    discarding the entire table.
    """
    for word in (
        "Missense",
        "Nonsense",
        "Frameshift",
        "Splice",
        "Domain",
        "N-Terminal",
        "C-Terminal",
        "Synonymous",
    ):
        assert not _gene_symbol_tokens(word), word

    row = ["Exon 1", "c.2T>C", "p.Met1Thr", "Missense", "N-Terminal"]
    assert not _row_has_off_target_gene_without_target(row, "KCNQ1")


def test_real_off_target_gene_symbols_still_trip_the_guard():
    """The guard must keep working on genuine misrouted panel rows.

    Digit-free symbols (LMNA, TTN, ENG) must keep counting as genes — that is
    the cross-gene contamination defense, so the annotation vocabulary is
    excluded by name rather than by a "must look like a gene" shape rule.
    """
    for symbol in ("SCN5A", "KCNH2", "CACNA1C", "BMPR2", "LMNA", "TTN", "ENG"):
        assert _gene_symbol_tokens(symbol) == {symbol}, symbol
    assert _row_has_off_target_gene_without_target(
        ["SCN5A", "c.100A>G", "p.Lys34Arg"], "KCNQ1"
    )
    assert _row_has_off_target_gene_without_target(
        ["LMNA", "c.100A>G", "p.Lys34Arg"], "KCNQ1"
    )


def test_bordered_table_is_not_emitted_twice_by_the_borderless_branch():
    """A table WITH a separator must still be claimed by the separator path."""
    bordered = """Table 1. Variants
|Variant |Patients |
|---|---|
|p.Met1Thr |3 |
|p.Ser140Arg |1 |
"""
    tables = enumerate_markdown_tables(bordered)
    assert len(tables) == 1
    assert len(tables[0].data_lines) == 2
    assert _normalize_header(tables[0].header_cells[0]) == "variant"


def test_short_pipe_runs_are_not_promoted_to_tables():
    """Prose containing pipes must not become a table."""
    prose = "Some text | with a pipe\n|only one row |here |\nmore prose\n"
    assert enumerate_markdown_tables(prose, only_variant_like=False) == []


def test_caption_scoped_table_ignores_row_level_gene_token_noise():
    """Domain/mutation-type shorthand must not delete rows of a scoped table.

    `_gene_symbol_tokens` accepts any 3-12 char uppercase token, so CNBD, DII,
    DI-S6, FRAME, SHIFT, SITE and even raw nucleotide runs read as "some other
    gene". When the caption already names the target gene, the row-level guard
    is redundant and must not run.
    """
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Summary of putative LQT1-associated mutations in KCNQ1",
        header_line="|Variant |Mutation Type |Location |No. of patients |",
        header_cells=["Variant", "Mutation Type", "Location", "No. of patients"],
        data_lines=[
            "|p.Met1Thr |Missense |N-Terminal |3 |",
            "|p.His105Thrfs |Frame shift/del |CNBD |1 |",
            "|p.Ser140Arg |Splice site |DI-S6 |2 |",
        ],
        char_start=0,
        char_end=300,
    )
    mapping = {"protein": 0, "patient_count": 3}

    variants = parse_routed_table(table, mapping, "KCNQ1")

    assert len(variants) == 3, [v.get("protein_notation") for v in variants]
    assert [v["penetrance_data"]["total_carriers_observed"] for v in variants] == [
        3,
        1,
        2,
    ]


def test_caption_naming_another_gene_still_rejects_the_whole_table():
    """The cross-gene defense must survive the scoped-caption exemption.

    Note the scope check resolves gene names against the built-in roster, so an
    off-roster symbol (LMNA, TTN) yields an EMPTY scope and the table reads as
    unscoped — the row-level guard then still applies. That predates this
    exemption; `test_unscoped_caption_still_runs_the_off_target_row_guard`
    covers that path.
    """
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Summary of SCN5A mutations in Brugada syndrome",
        header_line="|Variant |Mutation Type |No. of patients |",
        header_cells=["Variant", "Mutation Type", "No. of patients"],
        data_lines=["|p.Met1Thr |Missense |3 |"],
        char_start=0,
        char_end=120,
    )
    mapping = {"protein": 0, "patient_count": 2}

    assert parse_routed_table(table, mapping, "KCNQ1") == []


def test_unscoped_caption_still_runs_the_off_target_row_guard():
    """A caption naming no gene keeps the misrouted-panel guard active."""
    table = MarkdownTable(
        table_id="T1",
        caption="Table 1. Variants identified by panel sequencing",
        header_line="|Variant |Notes |No. of patients |",
        header_cells=["Variant", "Notes", "No. of patients"],
        data_lines=[
            "|p.Met1Thr |seen in KCNQ1 |3 |",
            "|p.Arg1623Gln |SCN5A |4 |",
        ],
        char_start=0,
        char_end=200,
    )
    mapping = {"protein": 0, "patient_count": 2}

    variants = parse_routed_table(table, mapping, "KCNQ1")
    proteins = [v["protein_notation"] for v in variants]
    assert "p.Arg1623Gln" not in proteins, "off-target SCN5A row must be dropped"


# --------------------------------------------------------------------------
# Cross-gene attribution (PMID 21232165 regression)
#
# These tests deliberately drive `enumerate_markdown_tables` on raw markdown
# instead of hand-building `MarkdownTable(caption=...)`. Every other caption
# test in this file injects a caption already in "Table N." form, so the whole
# caption-derived guard stack stayed green while the real enumerator captured
# 0 captions out of 330 corpus tables and the cross-gene reject never fired.
# --------------------------------------------------------------------------

PMC_HEADING_PAPER = """# Paper

### Table 1

Mutations in BRCA1 gene in Slovenian population.

| BIC nomenclature | HGVS nomenclature | Protein change | No. of positive families |
|---|---|---|---|
| 300T > G | c.181T > G | C61G | 15 |
| 1806C > T | c.1687C > T | Q563X | 13 |

### Table 2

Mutations in BRCA2 gene in Slovenian population.

| BIC nomenclature | HGVS nomenclature | Protein change | No. of positive families |
|---|---|---|---|
| 8138delC | c.7910delC | R259X | 2 |
"""


MIXED_GENE_PAPER = """# Paper

### Table 3

Polymorphisms in BRCA1 and BRCA2 genes in Slovenian population.

| BIC nomenclature | HGVS nomenclature | Protein change | No. of positive families |
|---|---|---|---|
| BRCA1 |  |  |  |
| 1186A > G | c.1067A > G | Q356R | 25 |
| 4956A > G | c.4837A > G | S1613G | 66 |
| BRCA2 |  |  |  |
| 1093A > C | c.865A > C | N289H | 18 |
| 1342C > A | c.1114C > A | H372N | 160 |
"""


def _protein_calls(paper: str, table_index: int, gene: str) -> set[str]:
    table = enumerate_markdown_tables(paper)[table_index]
    mapping = {"cdna": 1, "protein": 2, "patient_count": 3}
    return {
        (row.get("protein_notation") or "").strip()
        for row in parse_routed_table(table, mapping, gene)
        if (row.get("protein_notation") or "").strip()
    }


def test_enumerate_captures_caption_under_markdown_heading():
    """PMC's XML→markdown emits '### Table 1' then the caption as a paragraph."""
    tables = enumerate_markdown_tables(PMC_HEADING_PAPER)
    captions = [t.caption for t in tables]
    assert captions == [
        "Table 1. Mutations in BRCA1 gene in Slovenian population.",
        "Table 2. Mutations in BRCA2 gene in Slovenian population.",
    ]


def test_caption_naming_other_gene_rejects_whole_table():
    """A BRCA1-captioned table must contribute nothing to a BRCA2 extraction."""
    assert _protein_calls(PMC_HEADING_PAPER, 0, "BRCA2") == set()
    assert _protein_calls(PMC_HEADING_PAPER, 1, "BRCA1") == set()


def test_caption_naming_target_gene_still_extracts():
    """The reject must not cost recall on correctly-scoped tables."""
    assert _protein_calls(PMC_HEADING_PAPER, 0, "BRCA1") == {"C61G", "Q563X"}
    assert _protein_calls(PMC_HEADING_PAPER, 1, "BRCA2") == {"R259X"}


def test_in_band_gene_divider_partitions_mixed_table():
    """A caption naming both genes passes the reject; rows still must split."""
    assert _protein_calls(MIXED_GENE_PAPER, 0, "BRCA1") == {"Q356R", "S1613G"}
    assert _protein_calls(MIXED_GENE_PAPER, 0, "BRCA2") == {"N289H", "H372N"}


def test_brca2_is_a_known_gene_for_caption_scoping():
    """BRCA2 was absent from BUILTIN_GENE_METADATA, so no caption could scope
    to it and every BRCA2-captioned table stayed unguarded."""
    from pipeline.table_router import _caption_gene_scope

    assert _caption_gene_scope("Table 5. Variants in BRCA2 gene") == {"BRCA2"}


def test_ambiguous_query_alias_alone_does_not_scope_a_caption():
    """ "FH" is "family history" far more often than it is LDLR.

    `_caption_gene_scope` drives a whole-table reject, so admitting a gene on an
    ambiguous alias discarded entire cardiac mutation tables.
    """
    from pipeline.table_router import _caption_gene_scope

    assert (
        _caption_gene_scope(
            "Table 2. Mutations in probands with FH of sudden cardiac death"
        )
        == set()
    )
    # ...but an informative query alias must still scope, and a real mention of
    # the gene alongside the ambiguous one must still resolve.
    assert _caption_gene_scope("eTable 1. LQT1 Mutations or Rare Variants") == {"KCNQ1"}
    assert _caption_gene_scope("Table 9. FH cohort LDLR variants") == {"LDLR"}


def test_gene_header_substring_does_not_mint_a_gene_column():
    """ "Generation" contains "gene". A false gene column is worse than none: it
    suppresses the caption reject AND makes the row filter compare variants
    against values like F1/F2, which dropped every row of the table."""
    paper = (
        "# P\n\nTable 1. Mutations in BRCA1 gene\n\n"
        "| Generation | cDNA | Protein | Carriers |\n|---|---|---|---|\n"
        "| F1 | c.181T>G | C61G | 15 |\n| F2 | c.1687C>T | Q563X | 13 |\n"
    )
    table = enumerate_markdown_tables(paper)[0]
    mapping = _infer_column_mapping_from_headers(table)
    assert "gene" not in (mapping or {})
    # the table keeps its own rows, and the caption still rejects the other gene
    assert len(parse_routed_table(table, mapping or {}, "BRCA1")) == 2
    assert parse_routed_table(table, mapping or {}, "BRCA2") == []


def test_real_gene_column_is_still_detected():
    paper = (
        "# P\n\nTable 1. Panel variants\n\n"
        "| Gene | cDNA | Protein | Carriers |\n|---|---|---|---|\n"
        "| TTN | c.43828C>T | R14610X | 12 |\n| MYBPC3 | c.1224-52G>A | R943X | 5 |\n"
    )
    table = enumerate_markdown_tables(paper)[0]
    assert "gene" in (_infer_column_mapping_from_headers(table) or {})


def test_bic_shaped_sample_id_does_not_steal_bmpr2_cdna_column():
    paper = (
        "# P\n\nTable 1. BMPR2 mutation carriers\n\n"
        "| Sample | cDNA | Number of patients |\n|---|---|---|\n"
        "| 4184del4 | c.1459G>A | 2 |\n"
    )
    table = enumerate_markdown_tables(paper)[0]

    mapping = _infer_column_mapping_from_headers(table, target_gene="BMPR2")
    variants = parse_routed_table(table, mapping or {}, "BMPR2")

    assert mapping is not None
    assert mapping["cdna"] == 1
    assert [variant["cdna_notation"] for variant in variants] == ["c.1459G>A"]


def test_borderless_table_captures_heading_caption():
    """The separator-less branch kept the pre-fix lookback, so borderless
    tables stayed captionless and therefore unscoped."""
    from pipeline.table_router import _caption_gene_scope

    paper = (
        "# P\n\n### Table 4\n\nUnclassified sequence variants in BRCA1 gene\n\n"
        "| Nucleotide | Protein | Carriers |\n"
        "| 212C>G | I31M | 3 |\n| 462C>A | P115S | 1 |\n"
    )
    tables = enumerate_markdown_tables(paper)
    assert tables
    assert _caption_gene_scope(tables[0].caption) == {"BRCA1"}


# --- Spreadsheet title-row promotion + self-describing frequency columns -----


def _sheet_table(title, header, continuation, *rows):
    """Render a pandas-style xlsx sheet: merged title row, then the real header."""
    width = len(header)
    cells = [title] + [f"Unnamed: {i}" for i in range(1, width)]
    lines = [
        "| " + " | ".join(cells) + " |",
        "|" + "|".join([":---"] * width) + "|",
        "| " + " | ".join(header) + " |",
    ]
    if continuation:
        lines.append("| " + " | ".join(continuation) + " |")
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


SHEET_TITLE = (
    "Supplementary Data 1: A list of all 1,781 variants and their significance"
)
SHEET_HEADER = [
    "Chr",
    "Pos",
    "Gene",
    "HGVS.c",
    "HGVS.p",
    "Carrier frequency",
    "Carrier frequency",
    "Final clinical significance",
]
SHEET_CONTINUATION = ["", "", "", "", "", "in 7,051 cases", "in 11,241 controls", ""]


def test_spreadsheet_title_row_is_promoted_to_header():
    text = _sheet_table(
        SHEET_TITLE,
        SHEET_HEADER,
        SHEET_CONTINUATION,
        [
            "17",
            "41197708",
            "BRCA1",
            "c.5266dupC",
            "p.Gln1756fs",
            "0.00028",
            "0",
            "Pathogenic",
        ],
        [
            "17",
            "41197712",
            "BRCA1",
            "c.5262G>A",
            "p.Ser1754Ser",
            "0.00014",
            "0.00018",
            "Benign",
        ],
    )
    tables = enumerate_markdown_tables(text)
    assert len(tables) == 1
    table = tables[0]
    assert table.header_cells[:5] == ["Chr", "Pos", "Gene", "HGVS.c", "HGVS.p"]
    # the wrapped continuation is folded into the header it belongs to
    assert table.header_cells[5] == "Carrier frequency in 7,051 cases"
    assert table.header_cells[6] == "Carrier frequency in 11,241 controls"
    # the discarded title becomes the caption, so caption gene-scoping still works
    assert table.caption == SHEET_TITLE
    # continuation row is consumed, not left as a data row
    assert len(table.data_lines) == 2


def test_promoted_sheet_maps_frequency_columns_to_affected_and_unaffected():
    text = _sheet_table(
        SHEET_TITLE,
        SHEET_HEADER,
        SHEET_CONTINUATION,
        [
            "17",
            "41197708",
            "BRCA1",
            "c.5266dupC",
            "p.Gln1756fs",
            "0.00028",
            "0",
            "Pathogenic",
        ],
        [
            "17",
            "41197712",
            "BRCA1",
            "c.5262G>A",
            "p.Ser1754Ser",
            "0.00014",
            "0.00018",
            "Benign",
        ],
    )
    table = enumerate_markdown_tables(text)[0]
    mapping = _infer_column_mapping_from_headers(table, target_gene="BRCA1")
    assert mapping["gene"] == 2
    assert mapping["cdna"] == 3
    assert mapping["protein"] == 4
    assert mapping["affected"] == 5
    assert mapping["unaffected"] == 6
    # the one-carrier-per-row fallback must NOT fire once real counts exist
    assert "patient_count" not in mapping

    variants = parse_routed_table(table, mapping, "BRCA1")
    by_cdna = {v["cdna_notation"]: v for v in variants}
    first = by_cdna["c.5266dupC"]["penetrance_data"]
    assert first["affected_count"] == 2  # round(0.00028 * 7051)
    assert first["unaffected_count"] == 0  # a sourced zero, not an inferred one
    second = by_cdna["c.5262G>A"]["penetrance_data"]
    assert second["affected_count"] == 1  # round(0.00014 * 7051)
    assert second["unaffected_count"] == 2  # round(0.00018 * 11241)


def test_hierarchical_group_header_is_not_promoted():
    """A spanned group label also renders continuation cells as "Unnamed: N".

    Promoting under one would destroy the case/control grouping, so the first
    cell must actually read like a sheet title.
    """
    text = "\n".join(
        [
            "| Cases | Unnamed: 1 | Controls | Unnamed: 3 |",
            "|:---|:---|:---|:---|",
            "| Variant | n | Variant | n |",
            "| c.100A>G | 4 | c.100A>G | 1 |",
            "| c.200C>T | 2 | c.200C>T | 0 |",
        ]
    )
    tables = enumerate_markdown_tables(text, only_variant_like=False)
    assert tables[0].header_cells[0] == "Cases"


def test_population_frequency_headers_are_never_treated_as_counts():
    for header in (
        "Allele frequency in 7,051 cases",  # per-chromosome: denominator is 2N
        "gnomAD frequency in 125,748 individuals",
        "ExAC carrier frequency",
        "Minor allele frequency",
        "MAF",
        "Allele frequency",
        "Frequency in 500 cases with a family history",  # sub-cohort
        "Carrier frequency in 300 individuals",  # side of the split unstated
    ):
        assert _clinical_frequency_denominator(header) is None, header


def test_clinical_frequency_headers_resolve_role_and_denominator():
    assert _clinical_frequency_denominator("Carrier frequency in 7,051 cases") == (
        "affected",
        7051,
        False,
    )
    assert _clinical_frequency_denominator("Carrier frequency in 11,241 controls") == (
        "unaffected",
        11241,
        False,
    )
    assert _clinical_frequency_denominator("Frequency in 1,000 patients (%)") == (
        "affected",
        1000,
        True,
    )


def test_frequency_conversion_abstains_rather_than_fabricating():
    # exact: 0.00667 * 7051 = 47.03, well inside the 5-decimal ulp window
    assert _count_from_frequency("0.00667", 7051, False) == (47, "ok")
    # a stated zero is an observation, not an abstention
    assert _count_from_frequency("0", 11241, False) == (0, "ok")
    # 0.48596 * 7051 = 3426.5 -> no integer is pinned; must NOT return 0
    value, reason = _count_from_frequency("0.48596", 7051, False)
    assert value is None and reason == "ambiguous_rounding"
    # 4 decimals against a 5-figure denominator cannot pin an integer at all
    assert _count_from_frequency("0.5653", 11241, False) == (
        None,
        "precision_too_coarse",
    )
    # a fraction outside [0, 1] is not a frequency
    assert _count_from_frequency("1.5", 100, False) == (None, "out_of_range")


def test_explicit_count_column_wins_over_a_derived_frequency():
    """A stated integer must never be shadowed by a rounded frequency."""
    header = ["Gene", "HGVS.c", "No. of carriers", "Carrier frequency"]
    continuation = ["", "", "", "in 1,000 cases"]
    text = _sheet_table(
        "Supplementary Table 2: variants and carrier counts in the cohort",
        header,
        continuation,
        ["BRCA1", "c.5266dupC", "7", "0.00650"],
    )
    table = enumerate_markdown_tables(text)[0]
    mapping = _infer_column_mapping_from_headers(table, target_gene="BRCA1")
    assert mapping["patient_count"] == 2
    variants = parse_routed_table(table, mapping, "BRCA1")
    # 0.00650 * 1000 would round to 6; the stated 7 must win
    assert variants[0]["patients"]["count"] == 7


def test_external_population_frequency_with_a_cohort_size_is_still_rejected():
    """A stated denominator does not make a reference population a cohort."""
    for header in (
        "Frequency in 60,706 ExAC controls",
        "Frequency in 125,748 gnomAD controls",
        "Frequency in 2,504 1000 Genomes controls",
    ):
        assert _clinical_frequency_denominator(header) is None, header


def test_promoted_sheet_title_naming_a_gene_is_not_adopted_as_caption():
    """A caption scoping to another gene whole-table-rejects downstream."""
    text = _sheet_table(
        "Supplementary Data 2: BRCA2 variants observed across the cohort",
        ["Gene", "HGVS.c", "HGVS.p", "No. of carriers"],
        None,
        ["KCNQ1", "c.1032G>A", "p.Ala344Ala", "3"],
    )
    table = enumerate_markdown_tables(text)[0]
    assert table.header_cells[0] == "Gene"  # promotion still happened
    assert table.caption is None  # but the gene-naming title is not adopted
