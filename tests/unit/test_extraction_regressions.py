import json

from pipeline.extraction import ExpertExtractor
from pipeline.steps import extract_variants
from utils.models import ExtractionResult


def test_extract_variants_ignores_stale_abstract_cache_when_fulltext_exists(
    tmp_path, monkeypatch
):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    abstract_dir = tmp_path / "abstract_json"
    harvest_dir.mkdir()
    extraction_dir.mkdir()
    abstract_dir.mkdir()

    pmid = "12345678"
    fulltext = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
    fulltext.write_text(
        "# MAIN TEXT\n\n"
        + ("Table S1\nGene\nNucleotide\nAmino acids\nKCNH2\n1601G>A\nR534Q\n" * 120),
        encoding="utf-8",
    )
    abstract_record = abstract_dir / f"{pmid}.json"
    abstract_record.write_text(
        json.dumps({"abstract": "KCNH2 mutation carriers were studied."}),
        encoding="utf-8",
    )
    output_json = extraction_dir / f"KCNH2_PMID_{pmid}.json"
    output_json.write_text(
        json.dumps(
            {
                "variants": [],
                "extraction_metadata": {
                    "total_variants_found": 0,
                    "abstract_only": True,
                    "notes": "Abstract-only extraction; full text not available",
                    "model_used": "old-model",
                },
            }
        ),
        encoding="utf-8",
    )

    calls = []

    def fake_extract(self, paper):
        calls.append((paper.pmid, bool(paper.full_text), bool(paper.abstract)))
        return ExtractionResult(
            pmid=paper.pmid,
            success=True,
            model_used="fake-model",
            extracted_data={
                "variants": [
                    {
                        "gene_symbol": "KCNH2",
                        "cdna_notation": "c.1601G>A",
                        "protein_notation": "R534Q",
                    }
                ],
                "extraction_metadata": {"total_variants_found": 1},
            },
        )

    monkeypatch.setattr(ExpertExtractor, "extract", fake_extract)

    result = extract_variants(
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        gene_symbol="KCNH2",
        abstract_records={pmid: str(abstract_record)},
        abstract_only_pmids=[pmid],
        max_workers=1,
    )

    assert result.success is True
    assert calls == [(pmid, True, False)]
    saved = json.loads(output_json.read_text(encoding="utf-8"))
    metadata = saved["extraction_metadata"]
    assert metadata["abstract_only"] is False
    assert metadata["source_type"] == "fulltext"
    assert metadata["source_file"].endswith(f"{pmid}_FULL_CONTEXT.md")
    assert metadata["source_sha256"]
    assert saved["variants"][0]["protein_notation"] == "R534Q"


def test_extract_variants_does_not_treat_fallback_markdown_as_fulltext(
    tmp_path, monkeypatch
):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    abstract_dir = tmp_path / "abstract_json"
    harvest_dir.mkdir()
    extraction_dir.mkdir()
    abstract_dir.mkdir()

    pmid = "87654321"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# ABSTRACT-ONLY FALLBACK\n\n"
        "> **WARNING:** Full text could not be retrieved for PMID 87654321.\n"
        "> This document contains only the PubMed abstract and metadata.\n",
        encoding="utf-8",
    )
    abstract_record = abstract_dir / f"{pmid}.json"
    abstract_record.write_text(
        json.dumps({"abstract": "KCNH2 c.1601G>A R534Q was reported."}),
        encoding="utf-8",
    )

    calls = []

    def fake_extract(self, paper):
        calls.append((paper.pmid, bool(paper.full_text), bool(paper.abstract)))
        return ExtractionResult(
            pmid=paper.pmid,
            success=True,
            model_used="fake-model",
            extracted_data={
                "variants": [
                    {
                        "gene_symbol": "KCNH2",
                        "cdna_notation": "c.1601G>A",
                        "protein_notation": "R534Q",
                    }
                ],
                "extraction_metadata": {"total_variants_found": 1},
            },
        )

    monkeypatch.setattr(ExpertExtractor, "extract", fake_extract)

    result = extract_variants(
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        gene_symbol="KCNH2",
        abstract_records={pmid: str(abstract_record)},
        abstract_only_pmids=[pmid],
        max_workers=1,
    )

    assert result.success is True
    assert calls == [(pmid, False, True)]
    saved = json.loads((extraction_dir / f"KCNH2_PMID_{pmid}.json").read_text())
    metadata = saved["extraction_metadata"]
    assert metadata["abstract_only"] is True
    assert metadata["source_type"] == "abstract_only"


def test_vertical_gene_table_parser_recovers_stacked_supplement_rows():
    extractor = ExpertExtractor(models=["noop"], tier_threshold=0)
    text = """
Table S1. Subject Clinical and Genetic Characteristics
Patient
Gene
Nucleotide
Amino acids
1
RYR2
1258c>t
R420W
NT
2
KCNH2
3095g>a
R1032Q
NT
3
RYR2
1259 g>a
R420Q
NT
4
RYR2
1259 g>a
R420Q
NT
"""

    variants = extractor._parse_vertical_gene_table_variants(text, "RYR2")

    assert len(variants) == 2
    first, second = variants
    assert first["cdna_notation"] == "c.1258C>T"
    assert first["protein_notation"] == "R420W"
    assert first["patients"]["count"] == 1
    assert second["cdna_notation"] == "c.1259G>A"
    assert second["protein_notation"] == "R420Q"
    assert second["patients"]["count"] == 2
    assert second["penetrance_data"]["affected_count"] == 2


def test_table_hints_are_bounded_for_noisy_large_sources():
    extractor = ExpertExtractor(models=["noop"], tier_threshold=0)
    variants = [
        {
            "gene_symbol": "RYR2",
            "cdna_notation": f"c.{1000 + idx}A>G",
            "protein_notation": f"R{idx}W",
        }
        for idx in range(75)
    ]

    hints = extractor._format_table_hints(variants, max_hints=10)

    assert "c.1000A>G" in hints
    assert "c.1009A>G" in hints
    assert "c.1010A>G" not in hints
    assert "65 additional table-detected variants omitted" in hints
