import json

import pipeline.extraction as extraction_module
from pipeline.extraction import ExpertExtractor
from pipeline.steps import _write_dense_table_overflow_report, extract_variants
from harvesting.supplement_fold import FOLD_BEGIN, FOLD_END
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


def test_extract_variants_restricts_fulltext_to_candidate_pmids(tmp_path, monkeypatch):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    keep_pmid = "22222222"
    stale_pmid = "11111111"
    for pmid in (stale_pmid, keep_pmid):
        (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
            "# MAIN TEXT\n\n"
            + (
                "KCNH2 variant carriers were sequenced. "
                "Table 1 reports c.1601G>A p.R534Q in affected patients.\n" * 120
            ),
            encoding="utf-8",
        )

    calls = []

    def fake_extract(self, paper):
        calls.append(paper.pmid)
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
        candidate_pmids=[keep_pmid],
        max_workers=1,
    )

    assert result.success is True
    assert calls == [keep_pmid]
    assert (extraction_dir / f"KCNH2_PMID_{keep_pmid}.json").exists()
    assert not (extraction_dir / f"KCNH2_PMID_{stale_pmid}.json").exists()


def test_extract_variants_prefers_folded_full_context_over_stale_data_zones(
    tmp_path, monkeypatch
):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "22334455"
    (harvest_dir / f"{pmid}_DATA_ZONES.md").write_text(
        "# DATA ZONES\n\nKCNH2 cohort text without supplement rows.\n" * 100,
        encoding="utf-8",
    )
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        (
            "# FULL CONTEXT\n\n"
            + ("KCNH2 cohort text with clinical and genetic details.\n" * 100)
            + f"{FOLD_BEGIN}\n"
            + "# FOLDED SUPPLEMENTS (re-extraction aid)\n"
            + ("Supplement table reports c.1601G>A p.R534Q.\n" * 40)
            + f"{FOLD_END}\n"
        ),
        encoding="utf-8",
    )

    seen_sources = []

    def fake_extract(self, paper):
        seen_sources.append(paper.full_text)
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
        candidate_pmids=[pmid],
        max_workers=1,
    )

    assert result.success is True
    assert len(seen_sources) == 1
    assert "FOLDED SUPPLEMENTS" in seen_sources[0]
    saved = json.loads((extraction_dir / f"KCNH2_PMID_{pmid}.json").read_text())
    assert saved["extraction_metadata"]["source_file"].endswith(
        f"{pmid}_FULL_CONTEXT.md"
    )


def test_extract_variants_retries_failed_fulltext_by_default(tmp_path, monkeypatch):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "33333333"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# MAIN TEXT\n\n"
        + (
            "KCNH2 variant carriers were sequenced. "
            "Table 1 reports c.1601G>A p.R534Q in affected patients.\n" * 120
        ),
        encoding="utf-8",
    )

    calls = []

    def fake_extract(self, paper):
        calls.append(paper.pmid)
        if len(calls) == 1:
            return ExtractionResult(
                pmid=paper.pmid,
                success=False,
                error="RateLimitReached: upstream provider rate limited the request",
            )
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
        max_workers=1,
        extraction_retry_backoff_seconds=0,
    )

    assert result.success is True
    assert calls == [pmid, pmid]
    assert result.stats["papers_extracted"] == 1
    assert result.stats["initial_failures"] == 1
    assert result.stats["failures"] == 0
    assert result.stats["extraction_retry_attempted"] == 1
    assert result.stats["extraction_retry_succeeded"] == 1
    assert result.stats["extraction_retry_failed"] == 0
    assert (extraction_dir / f"KCNH2_PMID_{pmid}.json").exists()

    failure_csv = tmp_path / "pmid_status" / "extraction_failures.csv"
    assert failure_csv.read_text(encoding="utf-8").strip().splitlines() == [
        "PMID,Error"
    ]
    retry_summary = json.loads(
        (tmp_path / "pmid_status" / "extraction_retry_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert retry_summary["attempts"][0]["pmid"] == pmid
    assert retry_summary["attempts"][0]["status"] == "succeeded"
    assert retry_summary["final_failures"] == []


def test_extract_variants_keeps_final_failure_after_retry(tmp_path, monkeypatch):
    harvest_dir = tmp_path / "pmc_fulltext"
    extraction_dir = tmp_path / "extractions"
    harvest_dir.mkdir()
    extraction_dir.mkdir()

    pmid = "44444444"
    (harvest_dir / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# MAIN TEXT\n\n"
        + (
            "KCNH2 variant carriers were sequenced. "
            "Table 1 reports c.1601G>A p.R534Q in affected patients.\n" * 120
        ),
        encoding="utf-8",
    )

    calls = []

    def fake_extract(self, paper):
        calls.append(paper.pmid)
        return ExtractionResult(
            pmid=paper.pmid,
            success=False,
            error="RateLimitReached: still limited",
        )

    monkeypatch.setattr(ExpertExtractor, "extract", fake_extract)

    result = extract_variants(
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        gene_symbol="KCNH2",
        max_workers=1,
        extraction_retry_backoff_seconds=0,
    )

    assert result.success is True
    assert calls == [pmid, pmid]
    assert result.stats["papers_extracted"] == 0
    assert result.stats["initial_failures"] == 1
    assert result.stats["failures"] == 1
    assert result.stats["extraction_retry_attempted"] == 1
    assert result.stats["extraction_retry_succeeded"] == 0
    assert result.stats["extraction_retry_failed"] == 1

    failure_csv = tmp_path / "pmid_status" / "extraction_failures.csv"
    failure_rows = failure_csv.read_text(encoding="utf-8").strip().splitlines()
    assert failure_rows[0] == "PMID,Error"
    assert failure_rows[1].startswith(f"{pmid},")


def test_table_overflow_merge_chunks_and_records_qc(monkeypatch):
    monkeypatch.setattr(extraction_module, "TABLE_REGEX_MERGE_MAX_VARIANTS", 5)
    monkeypatch.setattr(extraction_module, "TABLE_REGEX_OVERFLOW_MERGE_MAX_VARIANTS", 8)
    monkeypatch.setattr(extraction_module, "TABLE_REGEX_OVERFLOW_CHUNK_SIZE", 3)

    extractor = ExpertExtractor(models=["noop"], tier_threshold=0)
    variants = []
    for idx in range(12):
        variant = {
            "gene_symbol": "KCNH2",
            "cdna_notation": f"c.{100 + idx}A>G",
            "protein_notation": f"p.Arg{100 + idx}Gly",
            "source_location": "router+deterministic table"
            if idx < 6
            else "table regex scan",
        }
        if idx < 6:
            variant["patients"] = {"count": idx + 1}
        variants.append(variant)

    extracted = {"variants": [], "extraction_metadata": {"total_variants_found": 0}}
    merged = extractor._merge_table_variants_with_overflow_qc(
        extracted, variants, "12345678"
    )

    metadata = merged["extraction_metadata"]
    overflow = metadata["table_merge_overflow"]
    assert len(merged["variants"]) == 8
    assert metadata["table_merge_candidate_count"] == 12
    assert overflow["candidate_count"] == 12
    assert overflow["normal_safety_cap"] == 5
    assert overflow["overflow_merge_cap"] == 8
    assert overflow["chunk_size"] == 3
    assert overflow["structured_candidate_count"] == 6
    assert overflow["regex_only_candidate_count"] == 6
    assert overflow["selected_for_merge"] == 8
    assert overflow["omitted_after_dedupe"] == 4
    assert overflow["merged_added"] == 8


def test_dense_table_overflow_report_summarizes_cap_trips_and_supplements(tmp_path):
    harvest_dir = tmp_path / "pmc_fulltext"
    output_dir = tmp_path
    harvest_dir.mkdir()
    pmid = "98765432"
    (harvest_dir / f"{pmid}_artifacts.json").write_text(
        json.dumps(
            {
                "main_text": {
                    "source": "publisher_html",
                    "chars": 10000,
                    "supplement_descriptions_count": 3,
                },
                "summary": {
                    "main_text_chars": 10000,
                    "supplement_count": 1,
                    "supplements_converted": 1,
                    "supplements_total_chars": 2500,
                },
                "supplements": [{"name": "supp1.xlsx"}],
            }
        ),
        encoding="utf-8",
    )
    extraction = ExtractionResult(
        pmid=pmid,
        success=True,
        model_used="fake-model",
        extracted_data={
            "variants": [{"protein_notation": "p.Arg1Gly"}],
            "extraction_metadata": {
                "source_type": "fulltext",
                "scanner_merge_skipped": {
                    "candidate_count": 151,
                    "safety_cap": 150,
                    "reason": "candidate_count_exceeds_safety_cap",
                },
                "table_merge_overflow": {
                    "candidate_count": 600,
                    "deduped_count": 550,
                    "normal_safety_cap": 500,
                    "overflow_merge_cap": 2000,
                    "selected_for_merge": 550,
                    "omitted_after_dedupe": 0,
                    "structured_candidate_count": 500,
                    "regex_only_candidate_count": 50,
                    "merged_added": 25,
                },
            },
        },
    )

    json_path, tsv_path, summary = _write_dense_table_overflow_report(
        extractions=[extraction],
        harvest_dir=harvest_dir,
        output_dir=output_dir,
    )

    assert summary == {
        "records": 1,
        "scanner_cap_trips": 1,
        "table_merge_cap_trips": 1,
        "table_candidates_omitted_after_dedupe": 0,
        "missing_supplement_ref_pmids": 1,
    }
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["records"][0]["pmid"] == pmid
    assert payload["records"][0]["supplement_refs"] == 3
    assert payload["records"][0]["supplements_downloaded"] == 1
    assert (
        tsv_path.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("pmid\tscanner_cap_tripped")
    )


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


def test_markdown_table_parser_prefers_carrier_over_total_case_denominator():
    extractor = ExpertExtractor(models=["noop"], tier_threshold=0)
    text = """
### Table 2

BRCA1 variants identified in familial breast and ovarian cancer patients

| cDNA | Protein | Variation type | Total case | Carrier |
|---|---|---|---|---|
| c.981_982delAT | p.Cys328* | Frameshift | 1142 | 18 |
"""

    variants = extractor._parse_markdown_table_variants(text, "BRCA1")

    assert len(variants) == 1
    variant = variants[0]
    assert variant["patients"]["count"] == 18
    assert variant["penetrance_data"]["total_carriers_observed"] == 18
    assert variant["count_provenance"]["carriers_column_label"] == "Carrier"
    assert variant["count_provenance"]["carriers_count_type"] == "per_variant_carrier"
    assert "Total case" not in variant["count_provenance"]["carriers_column_label"]


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
