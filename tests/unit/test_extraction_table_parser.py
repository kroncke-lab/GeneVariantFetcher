"""Regression tests for deterministic extraction table parsing."""

from pipeline.extraction import ExpertExtractor
from utils.models import ExtractionResult, Paper


def test_deterministic_parser_rejects_gwas_allele_table():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
#### Sheet: GWAS

| Locus | SNV | CHR | BP | EA | AA | EAF | beta | SE | p | n |
|-------|-----|-----|----|----|----|-----|------|----|---|---|
| KCNH2 | rs113843864 | 7 | 150618509 | G | A | 0.752 | 0.059 | 0.006 | 4.10E-23 | 29762 |
| KCNJ2 | rs4399570 | 17 | 68179345 | G | A | 0.699 | 0.142 | 0.006 | 5.30E-143 | 29508 |
"""

    variants = extractor._parse_markdown_table_variants(text, "KCNH2")

    assert variants == []


def test_deterministic_parser_keeps_clinical_variant_table():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 1. KCNH2 variant carriers

| protein | no. of patients | affected |
|---------|-----------------|----------|
| p.Lys897Thr | 7 | 3 |
| A561V | 2 | 2 |
"""

    variants = extractor._parse_markdown_table_variants(text, "KCNH2")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"p.Lys897Thr", "A561V"}
    assert by_protein["p.Lys897Thr"]["penetrance_data"]["total_carriers_observed"] == 7


def test_deterministic_parser_preserves_affected_unaffected_counts():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 1. KCNH2 affected and unaffected carriers

| protein | affected | unaffected |
|---------|----------|------------|
| p.Lys897Thr | 3 | 4 |
| p.Arg1047Leu | 4 | 10 |
"""

    variants = extractor._parse_markdown_table_variants(text, "KCNH2")

    by_protein = {v["protein_notation"]: v for v in variants}
    pen = by_protein["p.Lys897Thr"]["penetrance_data"]
    assert pen["total_carriers_observed"] == 7
    assert pen["affected_count"] == 3
    assert pen["unaffected_count"] == 4

    pen = by_protein["p.Arg1047Leu"]["penetrance_data"]
    assert pen["total_carriers_observed"] == 14
    assert pen["affected_count"] == 4
    assert pen["unaffected_count"] == 10


def test_deterministic_parser_infers_one_carrier_per_patient_row_without_count():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 1. KCNH2 mutation-positive patients

| Patient | Mutation | QTc | Phenotype |
|---------|----------|-----|-----------|
| P1 | p.Arg176Trp | 510 | syncope |
| P2 | p.Asn629Ser | 430 | unaffected |
"""

    variants = extractor._parse_markdown_table_variants(text, "KCNH2")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"p.Arg176Trp", "p.Asn629Ser"}
    pen = by_protein["p.Arg176Trp"]["penetrance_data"]
    assert pen["total_carriers_observed"] == 1
    assert pen["affected_count"] == 1
    assert pen["unaffected_count"] == 0
    pen = by_protein["p.Asn629Ser"]["penetrance_data"]
    assert pen["total_carriers_observed"] == 1
    assert pen["affected_count"] == 0
    assert pen["unaffected_count"] == 1


def test_artifact_filter_removes_malformed_protein_notations():
    extractor = ExpertExtractor(models=["gpt-4"])
    data = {
        "extraction_metadata": {},
        "variants": [
            {"gene_symbol": "KCNH2", "protein_notation": "A"},
            {"gene_symbol": "KCNH2", "protein_notation": "0.734027"},
            {"gene_symbol": "KCNH2", "protein_notation": "p.Lys897Thr"},
            {"gene_symbol": "KCNH2", "cdna_notation": "c.2398+1G>A"},
        ],
    }

    filtered = extractor._filter_extraction_artifacts(data, "KCNH2")

    remaining = {
        v.get("protein_notation") or v.get("cdna_notation")
        for v in filtered["variants"]
    }
    assert remaining == {"p.Lys897Thr", "c.2398+1G>A"}
    assert filtered["extraction_metadata"]["malformed_filtered"] == 2


def test_low_yield_router_result_does_not_short_circuit_full_text(monkeypatch):
    extractor = ExpertExtractor(models=["test-model"], tier_threshold=1)
    paper = Paper(
        pmid="23098067",
        title="Large cohort",
        gene_symbol="KCNH2",
        full_text="KCNH2 mutation variant carrier table. " * 40,
    )

    router_variant = {
        "gene_symbol": "KCNH2",
        "protein_notation": "p.His240His",
        "penetrance_data": {"total_carriers_observed": 1},
        "source_location": "Router table",
    }

    def fake_router(_paper, _text):
        return ExtractionResult(
            pmid=_paper.pmid,
            success=True,
            extracted_data={
                "variants": [router_variant],
                "extraction_metadata": {"total_variants_found": 1},
            },
            model_used="router+stub",
        )

    class EmptyScanner:
        variants = []

        def get_hints_for_prompt(self, max_hints):
            return ""

    monkeypatch.setattr(extractor, "_try_table_router", fake_router)
    monkeypatch.setattr(extractor, "_extract_variants_from_tables", lambda *_: [])
    monkeypatch.setattr(
        "pipeline.extraction.scan_document_for_variants",
        lambda *_, **__: EmptyScanner(),
    )
    monkeypatch.setattr(
        extractor,
        "call_llm_json_with_status",
        lambda _prompt: (
            {
                "variants": [
                    {
                        "gene_symbol": "KCNH2",
                        "protein_notation": "p.Arg176Trp",
                    },
                    {
                        "gene_symbol": "KCNH2",
                        "protein_notation": "p.Leu552Ser",
                    },
                ],
                "extraction_metadata": {"total_variants_found": 2},
            },
            False,
            "{}",
        ),
    )

    result = extractor.extract(paper)

    assert result.success
    assert result.model_used == "test-model"
    proteins = {v.get("protein_notation") for v in result.extracted_data["variants"]}
    assert {"p.Arg176Trp", "p.Leu552Ser", "p.His240His"} <= proteins
