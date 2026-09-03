# Gold-row source presence sweep

Scope: 49 gene-paper attempts from inventory.tsv restricted to tranche_01.tsv

Classification is blind to predictions and reads only the gold row and
the source on disk (article body, converted supplements incl. .docx/.xlsx/
.doc/.pptx/.zip, article PDFs). Nothing under `corpus/` was written.

- Gold rows: **242** across 49 gene-paper attempts
- Hard ceiling (source absent or substitution absent from all text, no figures): **75 rows (30.992%)** — max paper-derived recall if excluded: 69.008%
- Wide ceiling (also figures-present and non-searchable notation): **93 rows (38.43%)** — max recall if excluded: 61.57%
- Non-zero carrier gold rows: 233; behind hard ceiling 74, behind wide ceiling 92
- Supplement files converted: 26; failed: 1 (papers with failures: 1)

## By class

| class | rows | % |
|---|---:|---:|
| `present_in_body` | 149 | 61.57 |
| `present_in_supplement_only` | 0 | 0.0 |
| `present_in_pdf_only` | 0 | 0.0 |
| `source_absent` | 9 | 3.719 |
| `text_absent_stub_body` | 8 | 3.306 |
| `text_absent_garbled_body` | 7 | 2.893 |
| `text_absent_figures_present` | 2 | 0.826 |
| `text_absent_notation_inconclusive` | 16 | 6.612 |
| `text_absent_substitution` | 51 | 21.074 |

## By gene

| gene | rows | present_in_body | present_in_supplement_only | present_in_pdf_only | source_absent | text_absent_stub_body | text_absent_garbled_body | text_absent_figures_present | text_absent_notation_inconclusive | text_absent_substitution |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| KCNH2 | 25 | 12 | 0 | 0 | 0 | 1 | 7 | 0 | 5 | 0 |
| KCNQ1 | 66 | 10 | 0 | 0 | 6 | 0 | 0 | 0 | 5 | 45 |
| MYBPC3 | 74 | 74 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| RYR2 | 23 | 23 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| SCN5A | 54 | 30 | 0 | 0 | 3 | 7 | 0 | 2 | 6 | 6 |

## By tranche

| tranche | rows | hard ceiling | wide ceiling | present_in_supplement_only |
|---|---:|---:|---:|---:|
| mixed_gold_01 | 242 | 75 | 93 | 0 |

## By gold provenance

| provenance | rows | hard ceiling | wide ceiling |
|---|---:|---:|---:|
| curated_override_mixed_provenance | 74 | 0 | 0 |
| human_curated_cardiac | 145 | 75 | 93 |
| human_curated_ryr2_pilot | 23 | 0 | 0 |

## Papers with the most text-absent gold rows (wide ceiling)

| gene | PMID | text-absent rows | gold rows | classes | figures | unfetched links | advertised unresolved supp refs | converted supps | body chars |
|---|---:|---:|---:|---|---:|---:|---|---:|---:|
| KCNQ1 | 14678125 | 41 | 41 | text_absent_notation_inconclusive=4, text_absent_substitution=37 | 0 | 0 | False | 0 | 13779 |
| KCNQ1 | 21131640 | 9 | 9 | text_absent_notation_inconclusive=1, text_absent_substitution=8 | 0 | 0 | False | 0 | 22771 |
| KCNH2 | 15176425 | 7 | 7 | text_absent_garbled_body=7 | 0 | 0 | False | 0 | 225302 |
| KCNQ1 | 15176425 | 6 | 6 | source_absent=6 | 0 | 0 | False | 0 | 0 |
| SCN5A | 15176425 | 5 | 5 | text_absent_stub_body=5 | 0 | 0 | False | 0 | 2967 |
| SCN5A | 12566525 | 4 | 4 | text_absent_notation_inconclusive=1, text_absent_substitution=3 | 0 | 0 | False | 0 | 8976 |
| KCNH2 | 15851119 | 3 | 4 | text_absent_notation_inconclusive=3 | 0 | 0 | False | 0 | 15010 |
| SCN5A | 17081365 | 3 | 3 | source_absent=3 | 0 | 0 | False | 0 | 0 |
| SCN5A | 27232914 | 3 | 5 | text_absent_notation_inconclusive=2, text_absent_substitution=1 | 0 | 0 | False | 0 | 31200 |
| SCN5A | 11804990 | 2 | 2 | text_absent_figures_present=2 | 41 | 0 | False | 2 | 104766 |
| SCN5A | 16764707 | 2 | 10 | text_absent_notation_inconclusive=2 | 0 | 0 | False | 0 | 47221 |
| SCN5A | 28779003 | 2 | 3 | text_absent_substitution=2 | 0 | 0 | False | 0 | 11173 |
| KCNH2 | 15466642 | 1 | 2 | text_absent_notation_inconclusive=1 | 0 | 0 | False | 0 | 48807 |
| KCNH2 | 16470702 | 1 | 1 | text_absent_notation_inconclusive=1 | 0 | 0 | False | 0 | 6168 |
| KCNH2 | 29881912 | 1 | 1 | text_absent_stub_body=1 | 0 | 0 | False | 0 | 2266 |
| SCN5A | 12574983 | 1 | 1 | text_absent_stub_body=1 | 0 | 0 | False | 0 | 3510 |
| SCN5A | 17971661 | 1 | 1 | text_absent_stub_body=1 | 0 | 0 | False | 0 | 2170 |
| SCN5A | 23349452 | 1 | 1 | text_absent_notation_inconclusive=1 | 0 | 0 | False | 0 | 72010 |
