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
    variants = [
        {"gene_symbol": "BRCA2", "protein_notation": "C61G"},  # Table 1 → BRCA1
        {"gene_symbol": "BRCA2", "protein_notation": "Q356R"},  # Table 3 BRCA1 block
        {"gene_symbol": "BRCA2", "protein_notation": "N289H"},  # Table 3 BRCA2 block
        {"gene_symbol": "BRCA2", "protein_notation": "N288del"},  # Table 5 → BRCA2
    ]
    assert _filter(extractor, variants, "BRCA2") == ["N289H", "N288del"]


def test_filter_is_symmetric_for_the_other_gene(extractor):
    variants = [
        {"gene_symbol": "BRCA1", "protein_notation": "C61G"},
        {"gene_symbol": "BRCA1", "protein_notation": "Q356R"},
        {"gene_symbol": "BRCA1", "protein_notation": "N289H"},
    ]
    assert _filter(extractor, variants, "BRCA1") == ["C61G", "Q356R"]


def test_unscoped_variants_are_never_dropped(extractor):
    """Silence must not reject: only a positive assignment elsewhere rejects."""
    variants = [{"gene_symbol": "BRCA2", "protein_notation": "V9999X"}]
    assert _filter(extractor, variants, "BRCA2") == ["V9999X"]


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
