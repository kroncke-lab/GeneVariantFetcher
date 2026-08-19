"""Source-grounded gene attribution (PMID 21232165 regression).

The bug these pin: `_filter_by_gene` only ever consulted `variant["gene_symbol"]`,
which the extractor stamps with the run's target gene — so the filter was a
tautology and could never reject anything. The BRCA1 and BRCA2 runs for PMID
21232165 produced variant lists that were identical except for that column.

These tests drive the real parser over real markdown rather than asserting on a
hand-built variant dict, because the original defect was precisely that the
guards were verified against synthetic input the pipeline never produces.
"""

import pytest

from pipeline.extraction import ExpertExtractor


TWO_GENE_PAPER = """# Paper

### Table 1

Mutations in BRCA1 gene in Slovenian population.

| BIC nomenclature | HGVS nomenclature | Protein change | No. of positive families |
|---|---|---|---|
| 300T > G | c.181T > G | C61G | 15 |
| 1806C > T | c.1687C > T | Q563X | 13 |

### Table 3

Polymorphisms in BRCA1 and BRCA2 genes in Slovenian population.

| BIC nomenclature | HGVS nomenclature | Protein change | No. of alleles |
|---|---|---|---|
| BRCA1 |  |  |  |
| 1186A > G | c.1067A > G | Q356R | 25 |
| BRCA2 |  |  |  |
| 1093A > C | c.865A > C | N289H | 18 |

### Table 5

Unclassified sequence variants in BRCA2 gene in Slovenian population.

| BIC nomenclature | HGVS nomenclature | Protein change | No. of alleles |
|---|---|---|---|
| 4091delTAA | c.3863delTAA | del288N | 2 |
"""


@pytest.fixture
def extractor():
    return ExpertExtractor.__new__(ExpertExtractor)


def _attribution(extractor):
    return ExpertExtractor._source_gene_attribution(extractor, TWO_GENE_PAPER)


def test_caption_attributes_variants_to_the_gene_it_names(extractor):
    attribution = _attribution(extractor)
    assert attribution["C61G"] == {"BRCA1"}
    assert attribution["Q563X"] == {"BRCA1"}


def test_in_band_divider_splits_a_two_gene_table(extractor):
    attribution = _attribution(extractor)
    assert attribution["Q356R"] == {"BRCA1"}
    assert attribution["N289H"] == {"BRCA2"}


def test_structural_allele_is_fail_open_not_collided(extractor):
    """Known residual, deliberately fail-open.

    The source writes `del288N` and the extractor emits `N288del`, so the two
    never match literally. A position+event key ("288DEL") would bridge them but
    carries no gene, no allele and no coordinate system, so KCNQ1 `c.1032delG`
    and RYR2 `c.1032_1034delAAG` collapse onto one key and the legitimate
    variant is deleted. Keeping an unattributed structural allele is the correct
    trade for a recall-critical filter.
    """
    attribution = _attribution(extractor)
    assert "288DEL" not in attribution
    variants = [{"gene_symbol": "BRCA1", "protein_notation": "N288del"}]
    assert _filter(extractor, variants, "BRCA1") == ["N288del"]


def _filter(extractor, variants, gene):
    data = {"variants": variants}
    out = ExpertExtractor._filter_by_gene(extractor, data, gene, TWO_GENE_PAPER)
    return [v.get("protein_notation") for v in out["variants"]]


def test_filter_drops_variants_the_source_assigns_to_another_gene(extractor):
    # Every variant self-reports the target gene, exactly as the extractor stamps it.
    # Real extractor output carries table provenance, and rejection requires
    # it: attribution is evidence about the table a row came from, not a
    # document-wide blocklist (see test_variant_without_table_provenance_is_kept).
    variants = [
        {"gene_symbol": "BRCA2", "protein_notation": "C61G", "source_table": "Table 1"},
        {
            "gene_symbol": "BRCA2",
            "protein_notation": "Q356R",
            "source_table": "Table 3",
        },
        {
            "gene_symbol": "BRCA2",
            "protein_notation": "N289H",
            "source_table": "Table 3",
        },
        {
            "gene_symbol": "BRCA2",
            "protein_notation": "N288del",
            "source_table": "Table 5",
        },
    ]
    assert _filter(extractor, variants, "BRCA2") == ["N289H", "N288del"]


def test_filter_is_symmetric_for_the_other_gene(extractor):
    variants = [
        {"gene_symbol": "BRCA1", "protein_notation": "C61G", "source_table": "Table 1"},
        {
            "gene_symbol": "BRCA1",
            "protein_notation": "Q356R",
            "source_table": "Table 3",
        },
        {
            "gene_symbol": "BRCA1",
            "protein_notation": "N289H",
            "source_table": "Table 3",
        },
    ]
    assert _filter(extractor, variants, "BRCA1") == ["C61G", "Q356R"]


def test_unscoped_variants_are_never_dropped(extractor):
    """Silence must not reject: only a positive assignment elsewhere rejects."""
    variants = [
        {
            "gene_symbol": "BRCA2",
            "protein_notation": "V9999X",
            "source_table": "Table 1",
        }
    ]
    assert _filter(extractor, variants, "BRCA2") == ["V9999X"]


def test_variant_without_table_provenance_survives_a_mere_mention(extractor):
    """Attribution is evidence about a TABLE, not a document-wide blocklist.

    A row from prose, a figure, or the text scanner carries no source_table.
    Judging it against a document-global map deleted it whenever an unrelated
    table listed the same notation. It is kept as long as ANY table that carries
    the notation attributes it to the target gene.
    """
    mixed = """# Paper

Table 1. Variants in BRCA2 gene

| cDNA | Protein | Carriers |
|---|---|---|
| c.9976A>T | R190W | 7 |

Table 5. Previously reported BRCA1 variants

| cDNA | Protein | Carriers |
|---|---|---|
| c.181T>G | R190W | 3 |
"""
    data = {"variants": [{"gene_symbol": "BRCA2", "protein_notation": "R190W"}]}
    out = ExpertExtractor._filter_by_gene(extractor, data, "BRCA2", mixed)
    assert [v["protein_notation"] for v in out["variants"]] == ["R190W"]


def test_unanimous_source_contradiction_rejects_even_without_provenance(extractor):
    """The text scanner reads a table as prose and labels provenance
    "Text scan (...)", so a row whose source reads
    `| BRCA1 | c.4964C>T | p.S1655F |` was never judged and BRCA1's S1655F and
    R1699Q published under BRCA2 (PMID 19563646).

    When EVERY table carrying the notation attributes it to another gene, the
    paper is unanimous and the row is wrong whichever pass produced it.
    """
    variants = [
        {
            "gene_symbol": "BRCA2",
            "protein_notation": "C61G",
            "source_location": "Text scan (protein_hgvs_short)",
        }
    ]
    assert _filter(extractor, variants, "BRCA2") == []
    # ...and the same row is kept for the gene the paper actually names.
    keep = [
        {
            "gene_symbol": "BRCA1",
            "protein_notation": "C61G",
            "source_location": "Text scan (protein_hgvs_short)",
        }
    ]
    assert _filter(extractor, keep, "BRCA1") == ["C61G"]


def test_variant_is_judged_against_its_own_table(extractor):
    """Two tables listing the same notation must not contaminate each other."""
    same_token = """# Paper

Table 1. Variants in BRCA2 gene

| cDNA | Protein | Carriers |
|---|---|---|
| c.9976A>T | R190W | 7 |

Table 5. Previously reported BRCA1 variants

| cDNA | Protein | Carriers |
|---|---|---|
| c.181T>G | R190W | 3 |
"""
    data = {
        "variants": [
            {
                "gene_symbol": "BRCA2",
                "protein_notation": "R190W",
                "source_table": "Table 1",
            },
            {
                "gene_symbol": "BRCA2",
                "protein_notation": "R190W",
                "source_table": "Table 5",
            },
        ]
    }
    out = ExpertExtractor._filter_by_gene(extractor, data, "BRCA2", same_token)
    kept = [v["source_table"] for v in out["variants"]]
    assert kept == ["Table 1"]


# --------------------------------------------------------------------------
# Cases an external review (grok-4.6) found unpinned. Each of these was a real
# defect at some point in this iteration; none were covered by the tests above,
# which only ever built the shape the code already accepted.
# --------------------------------------------------------------------------

PHENOTYPE_MENTION_PAPER = """# Paper

Table 1. Mutations in KCNQ1

| cDNA | Protein | Phenotype | Notes | Carriers |
|---|---|---|---|---|
| c.1032G>A | G314S | LQT2 | compared to HERG | 7 |
| c.477+1G>A | S140G | LQT1 | also seen in BRCA2 | 3 |
"""


PAIRED_TWO_GENE_PAPER = """# Paper

Table 2. Compound genotypes

| ID | Protein | Region | Second variant | Region |
|---|---|---|---|---|
| K  KCNQ1 | G314A | Pore | SCN5A  D1114N | DII/DIII |
"""


def test_gene_mention_in_phenotype_or_notes_does_not_reassign_the_row(extractor):
    """A cell that MENTIONS a gene is not a cell that ASSIGNS one.

    `LQT2` resolves to KCNH2 and `breast cancer 2` to BRCA2, so scanning every
    cell let a phenotype or notes column claim the row and delete the target
    gene's own variant.
    """
    data = {
        "variants": [
            {"gene_symbol": "KCNQ1", "protein_notation": "G314S"},
            {"gene_symbol": "KCNQ1", "protein_notation": "S140G"},
        ]
    }
    out = ExpertExtractor._filter_by_gene(
        extractor, data, "KCNQ1", PHENOTYPE_MENTION_PAPER
    )
    assert [v["protein_notation"] for v in out["variants"]] == ["G314S", "S140G"]


def test_paired_two_gene_row_attributes_per_cell(extractor):
    """PMID 20541041 pairs two genes on one line; the variant belongs to the
    gene in its OWN cell, not to whichever gene the row mentions first."""
    attribution = ExpertExtractor._source_gene_attribution(
        extractor, PAIRED_TWO_GENE_PAPER
    )
    assert attribution["D1114N"] == {"SCN5A"}
    assert attribution["G314A"] == {"KCNQ1"}


# --------------------------------------------------------------------------
# TASKS.md section 4c — the audit findings not covered by the first fix round.
# --------------------------------------------------------------------------


def test_off_roster_gene_can_partition_a_table(extractor):
    """Only 14 genes are registered, so `| TTN |` after `| BRCA1 |` was not a
    divider: `section_gene` stayed BRCA1, the TTN run lost every row, and BRCA1
    absorbed TTN's variant."""
    from pipeline.table_router import (
        _infer_column_mapping_from_headers,
        enumerate_markdown_tables,
        parse_routed_table,
    )

    paper = """# P

Table 1. Panel variants

| Nucleotide | Protein | Carriers |
|---|---|---|
| BRCA1 |  |  |
| c.181T>G | C61G | 15 |
| TTN |  |  |
| c.43828C>T | R14610X | 12 |
"""
    table = enumerate_markdown_tables(paper)[0]
    mapping = _infer_column_mapping_from_headers(table) or {}
    assert [
        r.get("cdna_notation") for r in parse_routed_table(table, mapping, "TTN")
    ] == ["c.43828C>T"]
    assert [
        r.get("cdna_notation") for r in parse_routed_table(table, mapping, "BRCA1")
    ] == ["c.181T>G"]


def test_borderless_divider_is_not_swallowed_into_the_header(extractor):
    """A divider has one populated cell — the same shape as a wrapped header
    fragment — so the separator-less branch concatenated it into the header
    ("Nucleotide" + "BRCA1") and the table never partitioned."""
    from pipeline.table_router import (
        _infer_column_mapping_from_headers,
        _is_header_continuation,
        enumerate_markdown_tables,
        parse_routed_table,
    )

    assert _is_header_continuation(["BRCA2", "", ""]) is False

    paper = """# P

Table 3. Polymorphisms in BRCA1 and BRCA2 genes

| Nucleotide | Protein | Carriers |
| BRCA1 |  |  |
| 1186A>G | Q356R | 25 |
| BRCA2 |  |  |
| 1093A>C | N289H | 18 |
"""
    got: dict = {}
    for table in enumerate_markdown_tables(paper):
        mapping = _infer_column_mapping_from_headers(table)
        if not mapping:
            continue
        for gene in ("BRCA1", "BRCA2"):
            got.setdefault(gene, []).extend(
                r.get("protein_notation")
                for r in parse_routed_table(table, mapping, gene)
            )
    assert [x for x in got["BRCA1"] if x] == ["Q356R"]
    assert [x for x in got["BRCA2"] if x] == ["N289H"]


def test_vertical_parser_does_not_pair_across_records(extractor):
    """The 8-line lookahead ran past the end of a record and paired this gene's
    cDNA with the NEXT gene's protein change, minting a variant that does not
    exist."""
    paper = """# P

Table 2. Variants

KCNQ1
c.1032G>A
-
KCNH2
c.1841C>T
A614V
"""
    rows = ExpertExtractor._parse_vertical_gene_table_variants(
        extractor, paper, "KCNQ1"
    )
    assert [(v.get("cdna_notation"), v.get("protein_notation")) for v in rows] == [
        ("c.1032G>A", None)
    ]


def test_gene_symbol_is_not_a_protein_change(extractor):
    """The protein regex is IGNORECASE, so "SCN5A" parsed as S-5-A."""
    assert ExpertExtractor._valid_table_protein(extractor, "SCN5A") is None
    assert ExpertExtractor._valid_table_protein(extractor, "p.R190W") == "p.R190W"


def test_fixed_width_caption_without_a_table_label_still_scopes(extractor):
    """PDF-layout captions often carry no "Table N" label, so the gene scope
    stayed empty and BRCA1 founder alleles shipped as BRCA2 with counts."""
    rows = """
c185delAG   12   pGlu23fs   Pathogenic   40
c5382insC    7   pGln1829fs Pathogenic   40
"""
    scoped = "Polymorphisms in the BRCA1 gene among probands\n" + rows
    assert (
        ExpertExtractor._parse_fixed_width_table_variants(extractor, scoped, "BRCA2")
        == []
    )
    assert (
        len(
            ExpertExtractor._parse_fixed_width_table_variants(
                extractor, scoped, "BRCA1"
            )
        )
        == 2
    )
    # A comparison sentence is not a caption and must not reject the table.
    compared = "Carriers compared with BRCA1 probands\n" + rows
    assert (
        len(
            ExpertExtractor._parse_fixed_width_table_variants(
                extractor, compared, "BRCA2"
            )
        )
        == 2
    )


def test_nonhuman_ortholog_paper_contributes_no_variants(extractor):
    """PMID 19944633 is "Single nucleotide variation in exon 11 of canine
    BRCA2". Its residue numbering is the dog's and does not map to human, yet 16
    variants were staged for a clinical browser."""
    canine = (
        "# MAIN TEXT\n\n## Single nucleotide variation in exon 11 of canine "
        "BRCA2 in healthy and cancerous mammary tissue\n\n"
        "| Protein | Carriers |\n|---|---|\n| p.Glu119Lys | 3 |\n"
    )
    data = {"variants": [{"gene_symbol": "BRCA2", "protein_notation": "p.Glu119Lys"}]}
    out = ExpertExtractor._filter_by_gene(extractor, data, "BRCA2", canine)
    assert out["variants"] == []
    assert out["extraction_metadata"]["dropped_nonhuman_ortholog"] == 1


def test_model_system_assay_of_human_variants_is_kept(extractor):
    """A yeast/mouse assay OF human variants is legitimate human data. Only a
    species adjective modifying the GENE means the animal's own ortholog."""
    yeast = (
        "# MAIN TEXT\n\n## Functional Interaction Between BRCA1 and DNA Repair in "
        "Yeast May Uncover a Role of RAD50, RAD51, MRE11A\n\n"
        "| Protein | Carriers |\n|---|---|\n| p.C61G | 1 |\n"
    )
    data = {"variants": [{"gene_symbol": "BRCA1", "protein_notation": "p.C61G"}]}
    out = ExpertExtractor._filter_by_gene(extractor, data, "BRCA1", yeast)
    assert [v["protein_notation"] for v in out["variants"]] == ["p.C61G"]


def test_variant_with_no_identifier_is_dropped(extractor):
    """Nameless rows reach the public browser as blank variants — 30 BMPR2,
    14 BRCA1, 9 BRCA2 in the 150-paper re-extraction."""
    data = {
        "variants": [
            {"gene_symbol": "BRCA2", "protein_notation": "p.R190W"},
            {"gene_symbol": "BRCA2", "additional_notes": "carrier total inferred"},
            # A structural event is legitimately notation-free and must survive.
            {
                "gene_symbol": "BRCA2",
                "variant_class": "exon_duplication",
                "structural_description": "duplication of exons 8-10",
            },
        ]
    }
    out = ExpertExtractor._filter_extraction_artifacts(extractor, data, "BRCA2")
    kept = out["variants"]
    assert [v.get("protein_notation") for v in kept if v.get("protein_notation")] == [
        "p.R190W"
    ]
    assert any(v.get("variant_class") == "exon_duplication" for v in kept)
    assert not any(v.get("additional_notes") for v in kept)
