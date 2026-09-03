# Gold-row source presence sweep

Scope: 1533 gene-paper attempts from inventory.tsv (including source_unavailable)

Classification is blind to predictions and reads only the gold row and
the source on disk (article body, converted supplements incl. .docx/.xlsx/
.doc/.pptx/.zip, article PDFs). Nothing under `corpus/` was written.

- Gold rows: **7316** across 1533 gene-paper attempts
- Hard ceiling (source absent or substitution absent from all text, no figures): **1321 rows (18.056%)** — max paper-derived recall if excluded: 81.944%
- Wide ceiling (also figures-present and non-searchable notation): **2236 rows (30.563%)** — max recall if excluded: 69.437%
- Non-zero carrier gold rows: 6442; behind hard ceiling 1134, behind wide ceiling 2005
- Supplement files converted: 626; failed: 17 (papers with failures: 5)

## By class

| class | rows | % |
|---|---:|---:|
| `present_in_body` | 5071 | 69.314 |
| `present_in_supplement_only` | 9 | 0.123 |
| `present_in_pdf_only` | 0 | 0.0 |
| `source_absent` | 649 | 8.871 |
| `text_absent_stub_body` | 137 | 1.873 |
| `text_absent_garbled_body` | 23 | 0.314 |
| `text_absent_figures_present` | 582 | 7.955 |
| `text_absent_notation_inconclusive` | 333 | 4.552 |
| `text_absent_substitution` | 512 | 6.998 |

## By gene

| gene | rows | present_in_body | present_in_supplement_only | present_in_pdf_only | source_absent | text_absent_stub_body | text_absent_garbled_body | text_absent_figures_present | text_absent_notation_inconclusive | text_absent_substitution |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| APOE | 12 | 7 | 0 | 0 | 0 | 0 | 0 | 2 | 2 | 1 |
| BRCA1 | 89 | 63 | 0 | 0 | 0 | 6 | 0 | 15 | 5 | 0 |
| BRCA2 | 110 | 107 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 |
| KCNH2 | 1001 | 582 | 0 | 0 | 57 | 9 | 23 | 153 | 108 | 69 |
| KCNQ1 | 1753 | 1223 | 0 | 0 | 184 | 16 | 0 | 123 | 57 | 150 |
| MYBPC3 | 134 | 129 | 0 | 0 | 0 | 0 | 0 | 4 | 1 | 0 |
| RYR2 | 1013 | 811 | 0 | 0 | 38 | 2 | 0 | 57 | 42 | 63 |
| SCN5A | 3204 | 2149 | 9 | 0 | 370 | 104 | 0 | 225 | 118 | 229 |

## By tranche

| tranche | rows | hard ceiling | wide ceiling | present_in_supplement_only |
|---|---:|---:|---:|---:|
| mixed_gold_01 | 242 | 75 | 93 | 0 |
| mixed_gold_02 | 303 | 24 | 71 | 0 |
| mixed_gold_03 | 399 | 39 | 94 | 0 |
| mixed_gold_04 | 150 | 65 | 79 | 0 |
| mixed_gold_05 | 177 | 21 | 30 | 0 |
| mixed_gold_06 | 511 | 25 | 66 | 0 |
| mixed_gold_07 | 593 | 41 | 152 | 0 |
| mixed_gold_08 | 134 | 27 | 38 | 1 |
| mixed_gold_09 | 247 | 38 | 90 | 0 |
| mixed_gold_10 | 204 | 35 | 52 | 0 |
| mixed_gold_11 | 317 | 52 | 68 | 0 |
| mixed_gold_12 | 153 | 47 | 63 | 0 |
| mixed_gold_13 | 306 | 207 | 222 | 0 |
| mixed_gold_14 | 170 | 15 | 60 | 0 |
| mixed_gold_15 | 110 | 13 | 25 | 0 |
| mixed_gold_16 | 158 | 47 | 69 | 0 |
| mixed_gold_17 | 173 | 15 | 88 | 0 |
| mixed_gold_18 | 138 | 11 | 29 | 0 |
| mixed_gold_19 | 204 | 10 | 46 | 0 |
| mixed_gold_20 | 344 | 28 | 56 | 0 |
| mixed_gold_21 | 130 | 20 | 29 | 0 |
| mixed_gold_22 | 194 | 47 | 76 | 0 |
| mixed_gold_23 | 191 | 38 | 66 | 0 |
| mixed_gold_24 | 276 | 32 | 76 | 0 |
| mixed_gold_25 | 116 | 4 | 16 | 0 |
| mixed_gold_26 | 344 | 55 | 66 | 0 |
| mixed_gold_27 | 198 | 17 | 95 | 0 |
| mixed_gold_28 | 199 | 50 | 65 | 0 |
| mixed_gold_29 | 426 | 24 | 57 | 8 |
| unassigned | 209 | 199 | 199 | 0 |

## By gold provenance

| provenance | rows | hard ceiling | wide ceiling |
|---|---:|---:|---:|
| collaborator_approved_nonexhaustive | 7 | 0 | 2 |
| curated_override_mixed_provenance | 340 | 7 | 37 |
| human_curated_cardiac | 5956 | 1211 | 1995 |
| human_curated_ryr2_pilot | 1013 | 103 | 202 |

## Papers with the most text-absent gold rows (wide ceiling)

| gene | PMID | text-absent rows | gold rows | classes | figures | unfetched links | advertised unresolved supp refs | converted supps | body chars |
|---|---:|---:|---:|---|---:|---:|---|---:|---:|
| SCN5A | 20129283 | 76 | 417 | text_absent_figures_present=76 | 3 | 0 | False | 0 | 59786 |
| KCNQ1 | 23631430 | 69 | 69 | source_absent=69 | 0 | 0 | False | 0 | 0 |
| KCNQ1 | 17192539 | 56 | 57 | text_absent_notation_inconclusive=5, text_absent_substitution=51 | 0 | 0 | False | 0 | 22331 |
| KCNQ1 | 32893267 | 43 | 252 | text_absent_figures_present=43 | 16 | 0 | False | 5 | 8378983 |
| KCNQ1 | 14678125 | 41 | 41 | text_absent_notation_inconclusive=4, text_absent_substitution=37 | 0 | 0 | False | 0 | 13779 |
| SCN5A | 23631430 | 41 | 41 | source_absent=41 | 0 | 0 | False | 0 | 0 |
| RYR2 | 27452199 | 34 | 34 | text_absent_notation_inconclusive=34 | 0 | 0 | False | 0 | 21659 |
| SCN5A | 21273195 | 33 | 34 | text_absent_notation_inconclusive=11, text_absent_substitution=22 | 0 | 0 | False | 0 | 25258 |
| SCN5A | 26746457 | 30 | 32 | text_absent_figures_present=30 | 1 | 0 | False | 1 | 27190 |
| KCNH2 | 14661677 | 29 | 29 | text_absent_figures_present=29 | 2 | 0 | False | 4 | 93136 |
| KCNH2 | 19038855 | 28 | 28 | text_absent_figures_present=28 | 4 | 0 | False | 0 | 35148 |
| SCN5A | 30059973 | 27 | 185 | text_absent_notation_inconclusive=27 | 0 | 0 | False | 0 | 147539 |
| RYR2 | 19398665 | 25 | 27 | text_absent_figures_present=25 | 4 | 0 | False | 0 | 27888 |
| SCN5A | 24631775 | 25 | 27 | text_absent_substitution=25 | 0 | 0 | False | 0 | 64851 |
| KCNH2 | 15840476 | 24 | 89 | text_absent_notation_inconclusive=23, text_absent_substitution=1 | 0 | 0 | False | 0 | 39800 |
| KCNH2 | 24606995 | 24 | 28 | text_absent_figures_present=24 | 5 | 0 | False | 3 | 125064 |
| KCNQ1 | 19716085 | 23 | 182 | text_absent_notation_inconclusive=23 | 0 | 0 | False | 1 | 63004 |
| SCN5A | 11901046 | 23 | 23 | text_absent_substitution=23 | 0 | 0 | False | 0 | 20643 |
| SCN5A | 22885917 | 22 | 22 | text_absent_notation_inconclusive=8, text_absent_substitution=14 | 0 | 0 | False | 0 | 33123 |
| KCNH2 | 29650123 | 21 | 22 | text_absent_notation_inconclusive=6, text_absent_substitution=15 | 0 | 0 | False | 1 | 52691 |
| SCN5A | 25163546 | 20 | 20 | text_absent_stub_body=20 | 0 | 0 | False | 0 | 6330 |
| KCNH2 | 16922724 | 19 | 20 | text_absent_notation_inconclusive=6, text_absent_substitution=13 | 0 | 0 | False | 0 | 24150 |
| SCN5A | 28781330 | 18 | 18 | text_absent_stub_body=18 | 0 | 0 | False | 0 | 3256 |
| KCNH2 | 23098067 | 17 | 22 | text_absent_figures_present=17 | 16 | 0 | False | 6 | 103766 |
| RYR2 | 23595086 | 17 | 18 | text_absent_substitution=17 | 0 | 0 | False | 0 | 26979 |
| KCNH2 | 10973849 | 16 | 60 | text_absent_garbled_body=16 | 5 | 0 | False | 0 | 51783 |
| SCN5A | 12106943 | 16 | 16 | text_absent_figures_present=16 | 2 | 0 | False | 0 | 35954 |
| SCN5A | 29325976 | 16 | 93 | text_absent_notation_inconclusive=8, text_absent_substitution=8 | 0 | 0 | False | 1 | 64703 |
| SCN5A | 24721456 | 14 | 14 | text_absent_substitution=14 | 0 | 0 | False | 0 | 16580 |
| BRCA1 | 10528853 | 13 | 13 | text_absent_figures_present=13 | 3 | 0 | False | 1 | 73412 |
| KCNH2 | 26496715 | 13 | 54 | text_absent_notation_inconclusive=12, text_absent_substitution=1 | 0 | 0 | False | 3 | 40571 |
| KCNQ1 | 16922724 | 13 | 13 | text_absent_notation_inconclusive=2, text_absent_substitution=11 | 0 | 0 | False | 0 | 24150 |
| KCNQ1 | 24606995 | 13 | 21 | text_absent_figures_present=13 | 5 | 0 | False | 3 | 125064 |
| SCN5A | 16643399 | 13 | 13 | text_absent_notation_inconclusive=3, text_absent_substitution=10 | 0 | 0 | False | 0 | 11501 |
| SCN5A | 20031634 | 13 | 13 | source_absent=13 | 0 | 0 | False | 0 | 0 |
| SCN5A | 22360817 | 13 | 21 | text_absent_substitution=13 | 0 | 0 | False | 0 | 20284 |
| KCNH2 | 23631430 | 12 | 12 | source_absent=12 | 0 | 0 | False | 0 | 0 |
| KCNH2 | 29622001 | 12 | 35 | text_absent_figures_present=12 | 10 | 0 | False | 5 | 122596 |
| KCNQ1 | 32470535 | 12 | 12 | text_absent_notation_inconclusive=3, text_absent_substitution=9 | 0 | 0 | False | 0 | 18147 |
| SCN5A | 15996170 | 12 | 12 | source_absent=12 | 0 | 0 | False | 0 | 0 |
