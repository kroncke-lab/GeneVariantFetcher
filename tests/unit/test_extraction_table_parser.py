"""Regression tests for deterministic extraction table parsing."""

from types import SimpleNamespace

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


def test_deterministic_parser_uses_plain_table_titles_to_filter_gene():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 1. Summary of putative LQT1-associated mutations in KCNQ1
| Region | Nucleotide | Variant | Mutation Type | Location | No. of patients |
|---|---|---|---|---|---|
| Exon 1 | 5 C>T | A2V | Missense | N-Terminal | 1 |

Table 2. Summary of putative LQT2-associated mutations in KCNH2
| Region | Nucleotide | Variant | Mutation Type | Location | No. of patients |
|---|---|---|---|---|---|
| Exon 1 | 47 A>C | D16A* | Missense | N-Terminal | 1 |

Table 3. Summary of putative LQT3-associated mutations in SCN5A
| Region | Nucleotide | Variant | Mutation Type | Location | No. of patients |
|---|---|---|---|---|---|
| Exon 2 | 155 C>A | P52H | Missense | N-Terminal | 1 |
"""

    variants = extractor._parse_markdown_table_variants(text, "SCN5A")

    assert [v["protein_notation"] for v in variants] == ["P52H"]
    assert variants[0]["source_location"].startswith("Table 3.")
    assert variants[0]["penetrance_data"]["total_carriers_observed"] == 1


def test_table_regex_skips_off_target_gene_column_rows():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
| Gene | Protein | cDNA |
|---|---|---|
| KCNH2 | p.(Gly1036Asp) | c.3107G>A |
| SCN5A | p.(Arg18Gln) | c.53G>A |
| SCN5A | p.(Gln1507_Pro1509del) | c.4519_4527delCAGAAGCCC |
"""

    variants = extractor._extract_variants_from_tables(text, "SCN5A")

    proteins = {v.get("protein_notation") for v in variants}
    cdnas = {v.get("cdna_notation") for v in variants}
    assert "R18Q" in proteins
    assert "P.GLN1507_PRO1509DEL" in proteins
    assert "c.53G>A" in cdnas
    assert "c.4519_4527delCAGAAGCCC" in cdnas
    assert "G1036D" not in proteins
    assert "c.3107G>A" not in cdnas


def test_deterministic_parser_filters_generic_noncardiac_table_titles():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 1. BRCA1 variants observed in hereditary cancer probands
| Nucleotide | Variant | No. of patients |
|---|---|---|
| 181 T>G | C61G | 3 |

Table 2. TP53 variants observed in hereditary cancer probands
| Nucleotide | Variant | No. of patients |
|---|---|---|
| 743 G>A | R248Q | 5 |
"""

    variants = extractor._parse_markdown_table_variants(text, "TP53")

    assert [v["protein_notation"] for v in variants] == ["R248Q"]
    assert variants[0]["cdna_notation"] == "c.743G>A"
    assert variants[0]["source_location"].startswith("Table 2.")


def test_deterministic_parser_marks_controls_unaffected_and_continues_to_later_table():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
### Table 2

Control variants found in reference alleles

| Exon | Nucleotide change | Variant | Mutation type | Number | Status |
|---|---|---|---|---|---|
| 2 | 52 C>T | R18W | Missense | 1 | Rare control |
| 2 | 100 C>T | R34C | Missense | 44 | Polymorphism |

### Table 4

Compendium of Brugada syndrome-associated SCN5A mutations

| Exon | Nucleotide change | Coding effect | Mutation type | Number | Center |
|---|---|---|---|---|---|
| Exon 28 | 5350 G>A | E1784K* | Missense | 14 | 1, 2, 5, 6, 7 |
| Exon 12 | 2582_2583delTT | F861WfsX90* | Frame shift | 11 | 3, 7 |
"""

    variants = extractor._parse_markdown_table_variants(text, "SCN5A")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"R18W", "R34C", "E1784K", "F861WfsX90"}
    assert by_protein["R18W"]["clinical_significance"] == "benign"
    assert by_protein["R18W"]["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": 0,
        "unaffected_count": 1,
    }
    assert by_protein["R34C"]["penetrance_data"] == {
        "total_carriers_observed": 44,
        "affected_count": 0,
        "unaffected_count": 44,
    }
    assert by_protein["E1784K"]["cdna_notation"] == "c.5350G>A"
    assert by_protein["E1784K"]["penetrance_data"]["total_carriers_observed"] == 14
    assert by_protein["F861WfsX90"]["source_location"] == "Table 4"


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


def test_markdown_caption_row_scopes_gene_without_gene_column():
    """A caption emitted as a table row (no per-row Gene column) must scope the
    table to the gene it names. Regression for PMID 26669661, where a KCNQ1
    supplemental table leaked all rows into an SCN5A run (1 -> 161 variants)."""
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
| **Supplemental Table 1. *KCNQ1* variants by location and type** | | | | |
| --- | --- | --- | --- | --- |
| **Amino acid change** | **Nucleotide change** | **Location** | **Total families** | **Total carriers** |
| p.(Tyr111Cys) | c.332A>G | N-terminus | 3 | 9 |
| p.(Pro117Leu) | c.350C>T | N-terminus | 1 | 3 |
"""

    # The caption names KCNQ1, so an SCN5A run must claim nothing.
    assert extractor._parse_markdown_table_variants(text, "SCN5A") == []

    # The same table is kept for the gene it actually describes.
    kcnq1 = extractor._parse_markdown_table_variants(text, "KCNQ1")
    assert {v["cdna_notation"] for v in kcnq1} == {"c.332A>G", "c.350C>T"}


def test_markdown_no_caption_no_gene_column_is_not_oversuppressed():
    """A gene-column-less table with no gene named anywhere stays claimable by
    the target gene (single-gene papers must not regress)."""
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
| **Amino acid change** | **Nucleotide change** | **Location** | **Total carriers** |
| --- | --- | --- | --- |
| p.(Tyr111Cys) | c.332A>G | N-terminus | 9 |
"""

    variants = extractor._parse_markdown_table_variants(text, "KCNQ1")
    assert [v["cdna_notation"] for v in variants] == ["c.332A>G"]


def test_markdown_gene_column_filters_rows_in_multigene_table():
    """A table with a per-row Gene column keeps only the target gene's rows."""
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
| **Gene** | **Amino acid change** | **Nucleotide change** | **Total carriers** |
| --- | --- | --- | --- |
| SCN5A | p.(Arg18Gln) | c.53G>A | 4 |
| KCNH2 | p.(Pro1093Leu) | c.3278C>T | 2 |
| KCNQ1 | p.(Tyr111Cys) | c.332A>G | 9 |
"""

    variants = extractor._parse_markdown_table_variants(text, "SCN5A")
    assert [v["cdna_notation"] for v in variants] == ["c.53G>A"]


def test_fixed_width_parser_reads_pdftotext_layout_rows():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 5: SCN5A mutations (442 patients, 445 mutations, 185 unique mutations)

SCN5A mutation (c.)            exon              Mutant (p.)    Functional effect   n
Truncation mutations [n=81 mutations, 44 distinct mutations]
c127C>T                       2         pArg43*                Loss of function     1
c393-1C>T                     4                                Loss of function     1
c5350G>A                      28        pGlu1784Lys            Gain of function     69
c4519-4527del                 26        pGln1507_Pro1509del                        9
"""

    variants = extractor._parse_fixed_width_table_variants(text, "SCN5A")

    by_variant = {v["protein_notation"] or v["cdna_notation"]: v for v in variants}
    assert set(by_variant) == {
        "p.Arg43*",
        "c.393-1C>T",
        "p.Glu1784Lys",
        "p.Gln1507_Pro1509del",
    }
    assert by_variant["p.Glu1784Lys"]["cdna_notation"] == "c.5350G>A"
    assert (
        by_variant["p.Glu1784Lys"]["penetrance_data"]["total_carriers_observed"] == 69
    )
    assert by_variant["c.393-1C>T"]["penetrance_data"]["affected_count"] == 1


def test_fixed_width_parser_reads_wes_gene_column_rows():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Supplemental Table 1: Compendium of variants identified by WES testing

Position          Desig   Gene    Location   Nucleotide         Amino Acid       Zygosity
Chr1: 116311152   VUS     CASQ2   Exon 1     c.11C>T            p.T4I            Het
Chr1: 237433824   VUS     RYR2    Exon 2     c.76G>A            p.A26T           Het
Chr1: 237608827   VUS     RYR2    IVS        c.1292+5T>C        N/A              Het
"""

    variants = extractor._parse_fixed_width_table_variants(text, "RYR2")

    by_variant = {v["protein_notation"] or v["cdna_notation"]: v for v in variants}
    assert set(by_variant) == {"p.A26T", "c.1292+5T>C"}
    assert by_variant["p.A26T"]["cdna_notation"] == "c.76G>A"
    assert by_variant["p.A26T"]["source_location"] == "Supplemental Table 1"
    assert by_variant["c.1292+5T>C"]["penetrance_data"]["affected_count"] == 1


def test_fixed_width_parser_reads_patient_mutation_rows_under_gene_table_title():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Supplementary Table 2. RyR2 mutations in CPVT patients.

  # Patient   RyR2 Mutation   Blood glucose (mg/dl) 2h
        #1       P2404T                 187
        #2       K4392R                 147
        #3       L2534V                 192
"""

    variants = extractor._parse_fixed_width_table_variants(text, "RYR2")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"P2404T", "K4392R", "L2534V"}
    assert by_protein["P2404T"]["source_location"].startswith("Supplementary Table 2.")
    assert by_protein["P2404T"]["penetrance_data"]["total_carriers_observed"] == 1


def test_fixed_width_lqt_etable_rows_are_gene_scoped_and_keep_protein_counts():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
eTable 1. LQT1 Mutations or Rare Variants
                          Mutation    site     Site             N   Female    Proband       Mean QTc     Syncope    CA/VF
              c.153C>A      p.Y51X    N        N/C-term         1        0              1         443          1            0
         c.1022C>T p.A341V        MS   S5-pore-S6   35      20        11        500        23         7
          c.1032 G>A p.A344A splicing       MS   S5-pore-S6   47      26        24        494        23         4
                 * p.A178-G189del     C-loop   non-pore MS      1        1              1         471          0            0
    ****** p.Q361_K362 ins RQ     C    N/C-term      1       1         1        470         1         0

eTable 2. LQT2 Mutations or Rare Variants
                         c.1003C>T     p.Q335X     N    N/C-term      3      2         1        481        2        0

eTable 3. LQT3 Mutations or Rare Variants
            c.5350G>A      p.E1784K C                    50         23       26         481         5         0 BrS(7)
"""

    variants = extractor._parse_fixed_width_table_variants(text, "KCNQ1")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {
        "p.Y51X",
        "p.A341V",
        "p.A344A",
        "p.A178-G189del",
        "p.Q361_K362insRQ",
    }
    assert by_protein["p.Y51X"]["cdna_notation"] == "c.153C>A"
    assert by_protein["p.A341V"]["penetrance_data"] == {
        "total_carriers_observed": 35,
        "affected_count": 23,
        "unaffected_count": 12,
    }
    assert by_protein["p.A344A"]["cdna_notation"] == "c.1032G>A"
    assert all(v["source_location"].startswith("eTable 1.") for v in variants)


def test_pdf_linearized_lqt_etable_reconstructs_rows_for_markdown_parser():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
eTable 1. LQT1 Mutations or Rare Variants
Mutation
site
Site
N
Female
(n)
Proband
(n)
Mean QTc
(proband)
Syncope
(n)
CA/VF
(n)
c.521G>A  p.R147H
MS
non-pore MS
2
1
1
444
1
0
c.1111G>T  p.D317Y  MS
S5-pore-S6
2
1
1
495
2
0
c.1683G>T  p.R561S
C
N/C-term
3
2
0
1
0
c.3093insG  p.G1031fsX
1118
C
N/C-term
2
2
1
525
1
0
N: N-terminus, C: C-terminus, MS: membrane spanning, CA/VF: cardiac arrest
"""

    blocks = extractor._reconstruct_pdf_linearized_tables(text)

    assert len(blocks) == 1
    block = blocks[0]
    assert (
        "| cDNA | Protein | site | Site | No. of patients | Female (n) | "
        "Proband (n) | Mean QTc (proband) | affected | CA/VF (n) |"
    ) in block
    assert (
        "| c.521G>A | p.R147H | MS | non-pore MS | 2 | 1 | 1 | 444 | 1 | 0 |"
    ) in block
    assert (
        "| c.1111G>T | p.D317Y | MS | S5-pore-S6 | 2 | 1 | 1 | 495 | 2 | 0 |"
    ) in block
    assert "| c.1683G>T | p.R561S | C | N/C-term | 3 | 2 | 0 | - | 1 | 0 |" in block
    assert (
        "| c.3093insG | p.G1031fsX1118 | C | N/C-term | 2 | 2 | 1 | 525 | 1 | 0 |"
    ) in block

    variants = extractor._parse_markdown_table_variants(
        extractor._augment_pdf_linearized_tables(text), "KCNQ1"
    )

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {
        "p.R147H",
        "p.D317Y",
        "p.R561S",
        "p.G1031fsX1118",
    }
    assert by_protein["p.R147H"]["cdna_notation"] == "c.521G>A"
    assert by_protein["p.R147H"]["penetrance_data"] == {
        "total_carriers_observed": 2,
        "affected_count": 1,
        "unaffected_count": 1,
    }
    assert by_protein["p.D317Y"]["penetrance_data"] == {
        "total_carriers_observed": 2,
        "affected_count": 2,
        "unaffected_count": 0,
    }
    assert by_protein["p.R561S"]["penetrance_data"] == {
        "total_carriers_observed": 3,
        "affected_count": 1,
        "unaffected_count": 2,
    }
    assert by_protein["p.G1031fsX1118"]["penetrance_data"] == {
        "total_carriers_observed": 2,
        "affected_count": 1,
        "unaffected_count": 1,
    }


def test_pdf_linearized_table_updates_stale_estimate_for_deterministic_short_circuit():
    extractor = ExpertExtractor(models=["test-model"], tier_threshold=0)
    row_blocks = []
    for idx in range(110):
        pos = 100 + idx
        row_blocks.extend(
            [
                f"c.{pos}G>A  p.R{pos}H",
                "MS",
                "non-pore MS",
                "2",
                "1",
                "1",
                "444",
                "1",
                "0",
            ]
        )
    text = (
        """
eTable 1. LQT1 Mutations or Rare Variants
Mutation
site
Site
N
Female
(n)
Proband
(n)
Mean QTc
(proband)
Syncope
(n)
CA/VF
(n)
"""
        + "\n".join(row_blocks)
        + "\nN: N-terminus, C: C-terminus, MS: membrane spanning\n"
        + ("Long QT cohort methods and results. " * 80)
    )

    result = extractor._attempt_extraction(
        Paper(
            pmid="30758498",
            title="LQT1 supplement",
            gene_symbol="KCNQ1",
            full_text=text,
        ),
        "test-model",
        text,
        estimated_variants=0,
    )

    assert result.success
    assert result.model_used == "deterministic-table-parser"
    assert len(result.extracted_data["variants"]) == 110


def test_fixed_width_parser_reads_lqts_compendium_summary_rows():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 2. Compendium Summary of All Mutations
Gene  Exon      Nucleotide Change   Mutation        Mutation Type Location         Number Ethnicity   Status
KCNQ1         1 153 C>G             Y51X            Nonsense          N-Terminal        1W            Case
KCNQ1         1 Dup 160-168         IAP 54-56 dup   Inframe insertion N-Terminal        5 B>H         Polymorphism
KCNQ1         1 Del 211-219         71-73 del AAP   Inframe deletion N-Terminal         1W            Case
KCNQ1 intron 2 477+5G>C             M159sp          Splice site       S2 Domain         2W            Case
KCNQ1         6 Del 826-828         276 del S       Inframe deletion S5 domain          1W            Case
KCNQ1         9 Ins 1177-1179 TGG   392 ins W       Inframe insertion C-Terminal        1W            Case
KCNQ1         16 Del 1842-1844      614 del H       Inframe deletion C-terminus (SAD)   1W            Case
KCNH2         1 77-78 GC>TT         S26I            Missense          N-Terminal        1W            Case
"""

    variants = extractor._parse_fixed_width_table_variants(text, "KCNQ1")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {
        "Y51X",
        "I54ins",
        "P73del",
        "M159X",
        "S276del",
        "W392ins",
        "H614del",
    }
    assert by_protein["Y51X"]["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": 1,
        "unaffected_count": 0,
    }
    assert by_protein["I54ins"]["clinical_significance"] == "benign"
    assert by_protein["I54ins"]["penetrance_data"] == {
        "total_carriers_observed": 5,
        "affected_count": 0,
        "unaffected_count": 5,
    }
    assert all(v["source_location"].startswith("Table 2.") for v in variants)


def test_fixed_width_parser_reads_markitdown_lqts_compendium_rows():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 2. Compendium Summary of All Mutations
Gene Exon Nucleotide Change Mutation Mutation Type Location Number Ethnicity Status
| KCNQ1 | 1153 C>G | Y51X | Nonsense | N-Terminal | 1W | Case |
| ----- | -------- | ---- | -------- | ---------- | -- | ---- |
| KCNQ1 | 6817 C>T | L273F | Missense | S5 domain | 3W>Un | Case |
| KCNQ1 | 6898 G>A | A300T | Missense | Pore | 1A | Rare control |
| KCNH2 | 177-78 GC>TT | S26I | Missense | N-Terminal | 1W | Case |
"""

    variants = extractor._parse_fixed_width_table_variants(text, "KCNQ1")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"Y51X", "L273F", "A300T"}
    assert by_protein["L273F"]["penetrance_data"] == {
        "total_carriers_observed": 3,
        "affected_count": 3,
        "unaffected_count": 0,
    }
    assert by_protein["A300T"]["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": 0,
        "unaffected_count": 1,
    }


def test_fixed_width_parser_reads_vertical_lqts_compendium_rows():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 2. Compendium Summary of All Mutations
Gene
Exon
Nucleotide Change
Mutation
Mutation Type
Location
Number Ethnicity
Status
KCNQ1
1 153 C>G
Y51X
Nonsense
N-Terminal
1 W
Case
KCNQ1
1 Dup 160-168
IAP 54-56 dup
Inframe insertion N-Terminal
5 B>H
Polymorphism
KCNQ1 intron 2
477+5G>C
M159sp
Splice site
S2 Domain
2 W
Case
KCNQ1
6 Del 826-828
276 del S
Inframe deletion S5 domain
1 W
Case
KCNH2
1 77-78 GC>TT
S26I
Missense
N-Terminal
1 W
Case
"""

    variants = extractor._parse_fixed_width_table_variants(text, "KCNQ1")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"Y51X", "I54ins", "M159X", "S276del"}
    assert by_protein["I54ins"]["penetrance_data"] == {
        "total_carriers_observed": 5,
        "affected_count": 0,
        "unaffected_count": 5,
    }
    assert by_protein["M159X"]["penetrance_data"]["affected_count"] == 2


def test_lqts_compendium_short_circuit_does_not_duplicate_vertical_parse():
    extractor = ExpertExtractor(models=["test-model"], tier_threshold=0)
    text = """
Table 2. Compendium Summary of All Mutations
Gene
Exon
Nucleotide Change
Mutation
Mutation Type
Location
Number Ethnicity
Status
KCNQ1
1 153 C>G
Y51X
Nonsense
N-Terminal
1 W
Case
KCNQ1
1 Dup 160-168
IAP 54-56 dup
Inframe insertion N-Terminal
5 B>H
Polymorphism
""" + ("This paragraph provides full article context. " * 20)

    result = extractor._attempt_extraction(
        Paper(pmid="19841300", title="LQTS EPV", gene_symbol="KCNQ1", full_text=text),
        "test-model",
        text,
        estimated_variants=2,
    )

    assert result.success
    assert result.model_used == "deterministic-fixed-width-table-parser"
    assert [v["protein_notation"] for v in result.extracted_data["variants"]] == [
        "Y51X",
        "I54ins",
    ]
    assert (
        result.extracted_data["extraction_metadata"]["deterministic_parser_counts"][
            "vertical"
        ]
        == 0
    )


def test_fixed_width_parser_reads_coding_effect_count_table():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table S1: List of Mutations by Coding Effect, Location, and Frequency in 406
LQT3 Patients.

                      Coding Effect    Location        COUNT
                      V125L            N-term              1
                      K1505_Q1507del   DIII/DIV            9
                      Q1507_P1509del   DIII/DIV           55
                      F1617del         DIV-S3/S4           5
TOTAL                       406
"""

    variants = extractor._parse_fixed_width_table_variants(text, "SCN5A")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {
        "V125L",
        "K1505_Q1507del",
        "Q1507_P1509del",
        "F1617del",
    }
    assert by_protein["Q1507_P1509del"]["penetrance_data"] == {
        "total_carriers_observed": 55,
        "affected_count": 55,
        "unaffected_count": 0,
    }


def test_fixed_width_parser_reads_scn5a_nssnv_case_control_table():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Supplemental Table 1: Properties of SCN5A nsSNVs
                Nucleotid                         BrS      LQT     Control
     Status      e Change    Variant    Region   (2111)   (2888)   (8975)
case nsSNV      2249 A>G   Q750R   DII-S2     0   1   0   Benign
control nsSNV   553 G>A    A185T   DI-S2/3    0   0   1   Benign
Polymorphism    2074 C>A   Q692K   DI/DII     0   0   5   Benign
"""

    variants = extractor._parse_fixed_width_table_variants(text, "SCN5A")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"Q750R", "A185T", "Q692K"}
    assert by_protein["Q750R"]["cdna_notation"] == "c.2249A>G"
    assert by_protein["Q750R"]["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": 1,
        "unaffected_count": 0,
    }
    assert by_protein["A185T"]["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": 0,
        "unaffected_count": 1,
    }
    assert by_protein["Q692K"]["penetrance_data"]["unaffected_count"] == 5


def test_fixed_width_parser_reads_scn5a_functional_table_missing_from_case_table():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Supplemental Table 1: Properties of SCN5A nsSNVs
                Nucleotid                         BrS      LQT     Control
     Status      e Change    Variant    Region   (2111)   (2888)   (8975)
case nsSNV      2249 A>G   Q750R   DII-S2     0   1   0   Benign

Supplemental Table 4: Functionally Characterized SCN5A nsSNVs
              Control
Mutation   Frequency (%)    EP Status                              Electrophysiology Reference
Q750R    0.00%   Abnormal      This study
R800L    0.00%   Abnormal      This study
M1320V   0.00%    Wildtype     This study
A1330T   0.00%   Abnormal      PMID: 16039271
"""

    variants = extractor._parse_fixed_width_table_variants(text, "SCN5A")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"Q750R", "R800L", "M1320V", "A1330T"}
    assert by_protein["Q750R"]["cdna_notation"] == "c.2249A>G"
    assert by_protein["R800L"]["clinical_significance"] == "pathogenic"
    assert by_protein["M1320V"]["clinical_significance"] == "benign"
    assert by_protein["A1330T"]["source_location"].startswith("Supplemental Table 4")


def test_fixed_width_parser_reads_clinical_mutation_table_with_coding_effect():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Table 2. Included SCN5A Mutations and Variants
Nucleotide Change              Coding Effect            Region
163C>T                         Q55X                     N-terminal
407T>C                         L136P (2)                DI-S1
1537delC                       R513VfsTer8              NA
4140_4142 delCAA               del1380N                 NA
"""

    variants = extractor._parse_fixed_width_table_variants(text, "SCN5A")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"Q55X", "L136P", "R513VfsTer8", "N1380del"}
    assert by_protein["Q55X"]["cdna_notation"] == "c.163C>T"
    assert by_protein["L136P"]["penetrance_data"] == {
        "total_carriers_observed": 2,
        "affected_count": 2,
        "unaffected_count": 0,
    }
    assert by_protein["N1380del"]["cdna_notation"] == "c.4140_4142delCAA"


def test_fixed_width_clinical_table_short_circuits_without_gold_standard():
    extractor = ExpertExtractor(models=["test-model"], tier_threshold=1)
    rows = "\n".join(
        f"{100 + idx}A>G                         A{idx + 1}G                    Region"
        for idx in range(20)
    )
    text = f"""
Table 2. Included KCNH2 Mutations and Variants
Nucleotide Change              Coding Effect            Region
{rows}
"""

    result = extractor._attempt_extraction(
        Paper(pmid="123", title="KCNH2 table", gene_symbol="KCNH2", full_text=text),
        "test-model",
        text,
        estimated_variants=20,
    )

    assert result.success
    assert result.model_used == "deterministic-fixed-width-table-parser"
    assert len(result.extracted_data["variants"]) == 20


def test_fixed_width_caption_scopes_arbitrary_noncardiac_gene():
    """A fixed-width caption naming an arbitrary (non-cardiac) gene must scope the
    table to that gene, mirroring the markdown parser's open-vocab gene scope.
    Regression for multi-gene leakage: a BRCA1 supplement must claim nothing in
    an MYH7 run even though BRCA1 is outside the hard-coded cardiac gene set."""
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Supplemental Table 1. BRCA1 mutations in cancer patients

Nucleotide Change              Coding Effect            Region
163C>T                         Q55X                     N-terminal
407T>C                         L136P                    DI-S1
"""

    # Pre-fix behavior: the BRCA1 caption is outside TABLE_LABEL_EXPLICIT_GENES,
    # so the off-target guard never tripped and these two rows leaked into MYH7
    # (yielded 2 variants). After the fix the caption scopes to BRCA1 only.
    assert extractor._parse_fixed_width_table_variants(text, "MYH7") == []

    # The same table is still claimable for the gene it actually describes.
    brca1 = extractor._parse_fixed_width_table_variants(text, "BRCA1")
    by_protein = {v["protein_notation"]: v for v in brca1}
    assert set(by_protein) == {"Q55X", "L136P"}


def test_fixed_width_caption_scopes_contextual_all_alpha_gene():
    """All-letter genes need contextual scoping for novel-gene turnkey runs."""
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Supplemental Table 1. LMNA mutations in cardiomyopathy probands

Nucleotide Change              Coding Effect            Region
163C>T                         Q55X                     N-terminal
407T>C                         L136P                    DI-S1
"""

    assert extractor._parse_fixed_width_table_variants(text, "MYH7") == []

    lmna = extractor._parse_fixed_width_table_variants(text, "LMNA")
    by_protein = {v["protein_notation"]: v for v in lmna}
    assert set(by_protein) == {"Q55X", "L136P"}


def test_fixed_width_no_gene_caption_is_not_oversuppressed():
    """A fixed-width caption that names NO gene must stay claimable by the target
    gene. The noisy open-vocab tokens in a prose caption (COMPENDIUM/TESTING/...)
    must not be mistaken for genes and over-suppress a valid single-gene table."""
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Supplemental Table 1: Compendium of variants identified by WES testing

Position          Desig   Gene    Location   Nucleotide         Amino Acid       Zygosity
Chr1: 116311152   VUS     CASQ2   Exon 1     c.11C>T            p.T4I            Het
Chr1: 237433824   VUS     RYR2    Exon 2     c.76G>A            p.A26T           Het
"""

    variants = extractor._parse_fixed_width_table_variants(text, "RYR2")
    by_variant = {v["protein_notation"] or v["cdna_notation"]: v for v in variants}
    assert "p.A26T" in by_variant
    assert by_variant["p.A26T"]["cdna_notation"] == "c.76G>A"


def test_fixed_width_parser_reads_target_gene_functional_table_without_scn5a_name():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Supplemental Table 2: Functionally Characterized KCNH2 Variants
Mutation   Frequency (%)    EP Status                              Electrophysiology Reference
R56Q    0.00%   Abnormal      This study
K897T   0.10%   Wildtype      PMID: 12345678
"""

    variants = extractor._parse_fixed_width_table_variants(text, "KCNH2")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"R56Q", "K897T"}
    assert by_protein["R56Q"]["clinical_significance"] == "pathogenic"
    assert by_protein["K897T"]["clinical_significance"] == "benign"
    assert "KCNH2" in by_protein["R56Q"]["source_location"]


def test_deterministic_parser_skips_caption_unnamed_header_and_sums_case_columns():
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
#### Sheet: S3

| Gene | CDS | Protein | Vartype | Cases-Europe | Cases-Japan |
|---|---|---|---|---|---|
| SCN5A | c.3G>A | p.Met1? | start_lost | 99 | 99 |

#### Sheet: S4

| Table S4: List of rare protein altering variants in KCNQ1 cases | Unnamed: 1 | Unnamed: 2 | Unnamed: 3 | Unnamed: 4 | Unnamed: 5 |
|---|---|---|---|---|---|
| nan | nan | nan | nan | nan | nan |
| Gene | CDS | Protein | Vartype | Cases-Europe | Cases-Japan |
| KCNQ1 | c.2T>C | p.Met1? | start_lost | 1 | 2 |
| KCNQ1 | c.502G>A | p.Gly168Arg | missense | 9 | 1 |
"""

    variants = extractor._parse_markdown_table_variants(text, "KCNQ1")

    by_protein = {v["protein_notation"]: v for v in variants}
    assert set(by_protein) == {"p.Met1?", "p.Gly168Arg"}
    pen = by_protein["p.Met1?"]["penetrance_data"]
    assert pen["total_carriers_observed"] == 3
    assert pen["affected_count"] == 3
    assert pen["unaffected_count"] == 0
    pen = by_protein["p.Gly168Arg"]["penetrance_data"]
    assert pen["total_carriers_observed"] == 10
    assert pen["affected_count"] == 10
    assert pen["unaffected_count"] == 0


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
    assert by_protein["p.Arg176Trp"]["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": 1,
        "unaffected_count": 0,
    }
    assert by_protein["p.Asn629Ser"]["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": 0,
        "unaffected_count": 1,
    }


def test_suppresses_repeated_study_wide_counts_from_llm_variants():
    extractor = ExpertExtractor(models=["gpt-4"])
    data = {
        "variants": [
            {
                "protein_notation": f"p.Arg{i}Trp",
                "patients": {"count": 43},
                "penetrance_data": {
                    "total_carriers_observed": 43,
                    "affected_count": 28,
                    "unaffected_count": 15,
                },
                "additional_notes": "LLM extraction from cohort summary",
            }
            for i in range(1, 5)
        ],
        "extraction_metadata": {},
    }

    cleaned = extractor._suppress_repeated_study_wide_counts(data)

    assert cleaned["extraction_metadata"]["study_wide_count_suppressed"] == 4
    for variant in cleaned["variants"]:
        assert variant["patients"]["count"] is None
        assert variant["penetrance_data"]["total_carriers_observed"] is None
        assert variant["penetrance_data"]["affected_count"] is None
        assert variant["penetrance_data"]["unaffected_count"] is None
        assert "Cleared repeated cohort-wide count tuple" in variant["additional_notes"]


def test_does_not_suppress_deterministic_table_counts():
    extractor = ExpertExtractor(models=["gpt-4"])
    data = {
        "variants": [
            {
                "protein_notation": f"p.Arg{i}Trp",
                "patients": {"count": 43},
                "penetrance_data": {
                    "total_carriers_observed": 43,
                    "affected_count": 43,
                    "unaffected_count": 0,
                },
                "additional_notes": "Parsed via deterministic fixed-width table parser",
            }
            for i in range(1, 5)
        ],
        "extraction_metadata": {},
    }

    cleaned = extractor._suppress_repeated_study_wide_counts(data)

    assert "study_wide_count_suppressed" not in cleaned["extraction_metadata"]
    assert all(v["patients"]["count"] == 43 for v in cleaned["variants"])


def test_suppresses_repeated_study_wide_context_fields_from_llm_variants():
    extractor = ExpertExtractor(models=["gpt-4"])
    data = {
        "variants": [
            {
                "protein_notation": f"p.Arg{i}Trp",
                "patients": {
                    "phenotype": (
                        "Arrhythmia/ECG phenotypes in 35% of designated "
                        "variant carriers overall"
                    ),
                    "demographics": "median age 61 years",
                },
                "penetrance_data": {},
                "additional_notes": "LLM extraction from cohort summary",
            }
            for i in range(1, 7)
        ],
        "extraction_metadata": {},
    }

    cleaned = extractor._suppress_repeated_study_wide_context_fields(data)

    assert cleaned["extraction_metadata"]["study_wide_context_suppressed"] == 12
    for variant in cleaned["variants"]:
        assert variant["patients"]["phenotype"] is None
        assert variant["patients"]["demographics"] is None
        assert "Cleared repeated cohort-wide phenotype" in variant["additional_notes"]
        assert (
            "Cleared repeated cohort-wide demographics" in variant["additional_notes"]
        )


def test_does_not_suppress_repeated_plain_disease_phenotype():
    extractor = ExpertExtractor(models=["gpt-4"])
    data = {
        "variants": [
            {
                "protein_notation": f"p.Arg{i}Trp",
                "patients": {"phenotype": "Long QT syndrome"},
                "penetrance_data": {},
            }
            for i in range(1, 7)
        ],
        "extraction_metadata": {},
    }

    cleaned = extractor._suppress_repeated_study_wide_context_fields(data)

    assert "study_wide_context_suppressed" not in cleaned["extraction_metadata"]
    assert all(
        v["patients"]["phenotype"] == "Long QT syndrome" for v in cleaned["variants"]
    )


def test_artifact_filter_removes_malformed_protein_notations():
    extractor = ExpertExtractor(models=["gpt-4"])
    data = {
        "extraction_metadata": {},
        "variants": [
            {"gene_symbol": "KCNH2", "protein_notation": "A"},
            {"gene_symbol": "KCNH2", "protein_notation": "0.734027"},
            {"gene_symbol": "KCNH2", "protein_notation": "378"},
            {"gene_symbol": "SCN5A", "protein_notation": "SCN5A"},
            {"gene_symbol": "KCNH2", "protein_notation": "p.Lys897Thr"},
            {"gene_symbol": "KCNH2", "cdna_notation": "c.2398+1G>A"},
            {
                "gene_symbol": "SCN5A",
                "protein_notation": "p.Lys1505_Gln1507del",
                "cdna_notation": "c.4513_4521del",
            },
        ],
    }

    filtered = extractor._filter_extraction_artifacts(data, "SCN5A")

    remaining = {
        v.get("protein_notation") or v.get("cdna_notation")
        for v in filtered["variants"]
    }
    assert remaining == {"p.Lys897Thr", "c.2398+1G>A", "p.Lys1505_Gln1507del"}
    assert filtered["extraction_metadata"]["malformed_filtered"] == 4


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


def test_large_scanner_result_skips_hints_and_merge(monkeypatch):
    import pipeline.extraction as extraction

    extractor = ExpertExtractor(models=["test-model"], tier_threshold=1)
    paper = Paper(
        pmid="26669661",
        title="Multi-gene supplemental table",
        gene_symbol="SCN5A",
        full_text="SCN5A variant supplemental table discussion. " * 120,
    )

    table_variant = {
        "gene_symbol": "SCN5A",
        "protein_notation": "p.Arg18Gln",
        "source_location": "Supplemental Table 2",
    }

    class LargeScanner:
        def __init__(self):
            self.variants = [
                SimpleNamespace(
                    normalized=f"p.Ala{idx}Val",
                    notation_type="protein",
                    confidence=0.95,
                    raw_text=f"A{idx}V",
                    source="PMID_26669661",
                )
                for idx in range(extraction.SCANNER_REGEX_MERGE_MAX_VARIANTS + 1)
            ]

        def get_hints_for_prompt(self, max_hints):
            raise AssertionError("oversized scanner result should not provide hints")

    def fail_scanner_merge(*args, **kwargs):
        raise AssertionError("oversized scanner result should not be merged")

    monkeypatch.setattr(extractor, "_try_table_router", lambda *_: None)
    monkeypatch.setattr(
        extractor, "_extract_variants_from_tables", lambda *_: [table_variant]
    )
    monkeypatch.setattr(
        extraction, "scan_document_for_variants", lambda *_, **__: LargeScanner()
    )
    monkeypatch.setattr(extraction, "merge_scanner_results", fail_scanner_merge)
    monkeypatch.setattr(
        extractor,
        "call_llm_json_with_status",
        lambda _prompt: (
            {
                "variants": [
                    {
                        "gene_symbol": "SCN5A",
                        "protein_notation": "p.Glu1784Lys",
                    }
                ],
                "extraction_metadata": {"total_variants_found": 1},
            },
            False,
            "{}",
        ),
    )

    result = extractor.extract(paper)

    assert result.success
    proteins = {v.get("protein_notation") for v in result.extracted_data["variants"]}
    assert proteins == {"p.Glu1784Lys", "p.Arg18Gln"}
    assert result.extracted_data["extraction_metadata"]["scanner_merge_skipped"] == {
        "candidate_count": extraction.SCANNER_REGEX_MERGE_MAX_VARIANTS + 1,
        "safety_cap": extraction.SCANNER_REGEX_MERGE_MAX_VARIANTS,
        "reason": "candidate_count_exceeds_safety_cap",
    }


def test_pairs_from_reconstructed_blocks_maps_both_directions():
    # The reconstructed table carries an explicit cDNA+protein pairing per row;
    # the map lets a cDNA-only (or protein-only) extracted row recover its partner.
    extractor = ExpertExtractor(models=["gpt-4"])
    block = "\n".join(
        [
            "eTable 1. Mutations",
            "| cDNA | Protein | No. of patients |",
            "| --- | --- | --- |",
            "| c.153C>A | p.Y51X | 1 |",
            "| c.521G>A | p.R147H | 2 |",
        ]
    )
    pairs = extractor._pairs_from_reconstructed_blocks([block])
    assert pairs["c.153c>a"] == "p.Y51X"
    assert pairs["c.521g>a"] == "p.R147H"
    assert pairs["p.y51x"] == "c.153C>A"  # reverse direction


def test_backfill_fills_missing_side_without_overwriting():
    extractor = ExpertExtractor(models=["gpt-4"])
    extractor._linearized_variant_pairs = {
        "c.153c>a": "p.Y51X",
        "p.y51x": "c.153C>A",
    }
    data = {
        "variants": [
            {"cdna_notation": "c.153C>A", "protein_notation": None},  # fill protein
            {"cdna_notation": "", "protein_notation": "p.Y51X"},  # fill cDNA
            {"cdna_notation": "c.999Z>Q", "protein_notation": None},  # absent: stays
            {"cdna_notation": "c.153C>A", "protein_notation": "p.KEEP"},  # keep
        ]
    }
    extractor._backfill_variant_notation_pairs(data)
    variants = data["variants"]
    assert variants[0]["protein_notation"] == "p.Y51X"
    assert variants[1]["cdna_notation"] == "c.153C>A"
    assert variants[2]["protein_notation"] is None
    assert variants[3]["protein_notation"] == "p.KEEP"


def test_backfill_is_noop_when_no_table_reconstructed():
    extractor = ExpertExtractor(models=["gpt-4"])
    extractor._linearized_variant_pairs = {}
    data = {"variants": [{"cdna_notation": "c.1A>T", "protein_notation": None}]}
    extractor._backfill_variant_notation_pairs(data)
    assert data["variants"][0]["protein_notation"] is None


def test_pdf_linearized_reconstruction_generalizes_to_noncardiac_table():
    # Generalization guard against overfitting: the reconstruction must fire on a
    # generic supplement table with NO cardiac-specific columns (no Syncope, QTc,
    # or CA/VF) — only generic Mutation/Region/N/Cases/Affected — so it works on
    # genes we have not manually curated. The detector keys off generic
    # variant + count column names, not the LQTS eTable it was first built on.
    extractor = ExpertExtractor(models=["gpt-4"])
    text = """
Supplementary Table 1. BRCA1 pathogenic variants
Mutation
Region
N
Cases
(n)
Affected
(n)
c.68A>G  p.Glu23Gly
RING
3
2
2
c.181T>G  p.Cys61Gly
BRCT
4
3
3
c.5123C>A  p.Ala1708Glu
BRCT
6
5
5
"""
    blocks = extractor._reconstruct_pdf_linearized_tables(text)
    assert len(blocks) == 1
    block = blocks[0]
    # cDNA and Protein columns were reconstructed for the generic variants.
    assert "| c.68A>G | p.Glu23Gly |" in block
    assert "| c.181T>G | p.Cys61Gly |" in block
    assert "| c.5123C>A | p.Ala1708Glu |" in block
    # And the pairing map (used for the post-extraction backfill) carries them.
    extractor._augment_pdf_linearized_tables(text)
    pairs = extractor._linearized_variant_pairs
    assert pairs.get("c.68a>g") == "p.Glu23Gly"
    assert pairs.get("c.5123c>a") == "p.Ala1708Glu"
