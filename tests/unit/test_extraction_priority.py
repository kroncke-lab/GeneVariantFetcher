import json

from pipeline.extraction_priority import prioritize_extraction_sources


def _write_abstract(path, *, title, abstract, year="2024", journal="Test Journal"):
    path.write_text(
        json.dumps(
            {
                "metadata": {
                    "pmid": path.stem,
                    "title": title,
                    "journal": journal,
                    "year": year,
                },
                "abstract": abstract,
            }
        ),
        encoding="utf-8",
    )


def test_prioritization_prefers_original_variant_count_fulltext(tmp_path):
    harvest = tmp_path / "pmc_fulltext"
    abstracts = tmp_path / "abstract_json"
    report_dir = tmp_path / "priority"
    harvest.mkdir()
    abstracts.mkdir()

    strong = harvest / "111_FULL_CONTEXT.md"
    strong.write_text(
        "\n".join(
            [
                "Methods",
                "We sequenced a cohort of 432 Alzheimer disease patients and 500 controls.",
                "Table 1. APOE genotype and carrier counts.",
                "Variant rs429358 carriers n=101; epsilon 4 carriers 101/432.",
            ]
            * 80
        ),
        encoding="utf-8",
    )
    review = harvest / "222_FULL_CONTEXT.md"
    review.write_text(
        "\n".join(
            [
                "This review summarizes APOE biology and Alzheimer disease.",
                "No original cohort or carrier counts are reported.",
            ]
            * 80
        ),
        encoding="utf-8",
    )
    _write_abstract(
        abstracts / "111.json",
        title="APOE variants in an Alzheimer disease cohort",
        abstract="Patients were genotyped for rs429358 with carrier counts.",
    )
    _write_abstract(
        abstracts / "222.json",
        title="APOE and Alzheimer disease: a review",
        abstract="A narrative review of prior studies.",
    )

    result = prioritize_extraction_sources(
        markdown_files=[review, strong],
        abstract_papers=[],
        gene_symbol="APOE",
        disease_terms=["Alzheimer disease"],
        abstract_records={
            "111": str(abstracts / "111.json"),
            "222": str(abstracts / "222.json"),
        },
        top_n=1,
        report_dir=report_dir,
    )

    assert [path.name for path in result.selected_markdown_files] == [
        "111_FULL_CONTEXT.md"
    ]
    assert result.candidates[0].pmid == "111"
    assert result.candidates[0].selected is True
    assert result.candidates[1].pmid == "222"
    assert result.candidates[1].selected is False
    assert (report_dir / "priority_candidates.tsv").exists()
    assert (report_dir / "priority_candidates.json").exists()
    assert (report_dir / "priority_pmids.txt").read_text(encoding="utf-8") == "111\n"


def test_prioritization_keeps_abstract_only_as_low_priority_candidate(tmp_path):
    abstracts = tmp_path / "abstract_json"
    report_dir = tmp_path / "priority"
    abstracts.mkdir()
    _write_abstract(
        abstracts / "333.json",
        title="BRCA1 carrier counts in ovarian cancer",
        abstract="Original cohort with BRCA1 c.68_69delAG carriers n=20.",
    )

    result = prioritize_extraction_sources(
        markdown_files=[],
        abstract_papers=[("333", str(abstracts / "333.json"))],
        gene_symbol="BRCA1",
        disease_terms=["ovarian cancer"],
        top_n=1,
        report_dir=report_dir,
    )

    assert result.selected_markdown_files == []
    assert result.selected_abstract_papers == [("333", str(abstracts / "333.json"))]
    assert result.candidates[0].source_kind == "abstract"
    assert result.candidates[0].selected is True


def test_prioritization_offset_selects_next_ranked_band(tmp_path):
    harvest = tmp_path / "pmc_fulltext"
    abstracts = tmp_path / "abstract_json"
    harvest.mkdir()
    abstracts.mkdir()
    paths = []
    for pmid, variant_count in [("101", 5), ("102", 3), ("103", 1)]:
        path = harvest / f"{pmid}_FULL_CONTEXT.md"
        path.write_text(
            "\n".join(
                [
                    "BRCA1 cohort with breast cancer patients and carrier counts.",
                    *([f"BRCA1 c.{pmid}A>G carriers n=2 Table 1."] * variant_count),
                ]
                * 30
            ),
            encoding="utf-8",
        )
        paths.append(path)
        _write_abstract(
            abstracts / f"{pmid}.json",
            title=f"BRCA1 carrier study {pmid}",
            abstract="BRCA1 variants and carrier counts.",
        )

    result = prioritize_extraction_sources(
        markdown_files=paths,
        abstract_papers=[],
        gene_symbol="BRCA1",
        disease_terms=["breast cancer"],
        abstract_records={p.stem: str(p) for p in abstracts.glob("*.json")},
        top_n=1,
        offset=1,
    )

    assert len(result.selected_markdown_files) == 1
    assert result.selected_candidates[0].rank == 2


def test_prioritization_uses_default_gene_aliases(tmp_path):
    harvest = tmp_path / "pmc_fulltext"
    abstracts = tmp_path / "abstract_json"
    harvest.mkdir()
    abstracts.mkdir()

    alias_hit = harvest / "444_FULL_CONTEXT.md"
    alias_hit.write_text(
        "\n".join(
            [
                "We sequenced hypertrophic cardiomyopathy probands.",
                "Table 1. cMyBP-C variants and carrier counts.",
                "p.Arg502Trp carriers n=4 in affected families.",
            ]
            * 50
        ),
        encoding="utf-8",
    )
    weak = harvest / "555_FULL_CONTEXT.md"
    weak.write_text(
        "Hypertrophic cardiomyopathy review without gene-specific variant counts.\n"
        * 50,
        encoding="utf-8",
    )
    _write_abstract(
        abstracts / "444.json",
        title="cMyBP-C variants in hypertrophic cardiomyopathy families",
        abstract="Original cohort with carrier counts.",
    )
    _write_abstract(
        abstracts / "555.json",
        title="Hypertrophic cardiomyopathy review",
        abstract="Review of prior studies.",
    )

    result = prioritize_extraction_sources(
        markdown_files=[weak, alias_hit],
        abstract_papers=[],
        gene_symbol="MYBPC3",
        disease_terms=["hypertrophic cardiomyopathy"],
        abstract_records={
            "444": str(abstracts / "444.json"),
            "555": str(abstracts / "555.json"),
        },
        top_n=1,
    )

    assert result.selected_markdown_files == [alias_hit]
    assert result.candidates[0].pmid == "444"
    assert result.candidates[0].signals["gene_title_mentions"] > 0
