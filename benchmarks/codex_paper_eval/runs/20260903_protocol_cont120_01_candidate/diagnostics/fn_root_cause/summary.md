# FN root cause — `20260903_protocol_cont120_01_candidate`

123 paper-derived false negatives.

| leaf | rows |
|---|---:|
| `acquisition` | 72 |
| `unknown_notation` | 40 |
| `model_missed` | 5 |
| `parser_dropped` | 3 |
| `projection_dropped` | 2 |
| `paper_row_lost_to_linkage_origin` | 1 |

| gene | PMID | FN | leaves |
|---|---:|---:|---|
| SCN5A | 25163546 | 20 | acquisition=20 |
| SCN5A | 20031634 | 13 | acquisition=13 |
| BRCA1 | 10528853 | 13 | unknown_notation=13 |
| KCNH2 | 17905336 | 11 | acquisition=6, unknown_notation=4, parser_dropped=1 |
| KCNH2 | 12402336 | 9 | acquisition=5, unknown_notation=3, model_missed=1 |
| RYR2 | 29447731 | 9 | acquisition=6, unknown_notation=3 |
| KCNH2 | 24103226 | 7 | acquisition=4, unknown_notation=3 |
| KCNQ1 | 22539601 | 6 | acquisition=6 |
| SCN5A | 17905336 | 5 | acquisition=4, unknown_notation=1 |
| RYR2 | 22677073 | 4 | model_missed=3, projection_dropped=1 |
| SCN5A | 17227473 | 3 | acquisition=1, model_missed=1, unknown_notation=1 |
| MYBPC3 | 21302287 | 3 | unknown_notation=3 |
| KCNH2 | 21185499 | 2 | unknown_notation=2 |
| SCN5A | 20451667 | 2 | acquisition=2 |
| KCNH2 | 18774102 | 2 | paper_row_lost_to_linkage_origin=1, projection_dropped=1 |
| KCNH2 | 18452873 | 1 | unknown_notation=1 |
| SCN5A | 21677263 | 1 | acquisition=1 |
| KCNQ1 | 22677073 | 1 | unknown_notation=1 |
| KCNH2 | 16969682 | 1 | unknown_notation=1 |
| KCNQ1 | 26481773 | 1 | parser_dropped=1 |
| KCNQ1 | 20167303 | 1 | acquisition=1 |
| KCNH2 | 10790218 | 1 | unknown_notation=1 |
| RYR2 | 18463390 | 1 | acquisition=1 |
| SCN5A | 27401036 | 1 | parser_dropped=1 |
| KCNQ1 | 30462975 | 1 | acquisition=1 |
| KCNQ1 | 17224687 | 1 | unknown_notation=1 |
| SCN5A | 20219395 | 1 | acquisition=1 |
| SCN5A | 8917568 | 1 | unknown_notation=1 |
| SCN5A | 28018021 | 1 | unknown_notation=1 |

## Rows the reading protocol could have found

| gene | PMID | variant | leaf | db layers | run_text | request | response | extraction | db | paper lane | linkage lane |
|---|---:|---|---|---|---|---|---|---|---|---|---|
| KCNH2 | 12402336 | `D501N` | model_missed | clinvar | True | True | False | False | True | False | True |
| KCNQ1 | 26481773 | `R594Q` | parser_dropped |  | True | True | True | False | False | False | False |
| SCN5A | 17227473 | `R1193Q` | model_missed |  | True | True | False | False | False | False | False |
| RYR2 | 22677073 | `A1136V` | model_missed | clinvar | True | True | False | False | True | False | True |
| RYR2 | 22677073 | `G4936R` | projection_dropped | regex_table | True | True | True | True | True | False | False |
| RYR2 | 22677073 | `Q2958R` | model_missed |  | True | True | False | False | False | False | False |
| RYR2 | 22677073 | `R4037C` | model_missed |  | True | True | False | False | False | False | False |
| KCNH2 | 17905336 | `R534C` | parser_dropped | clinvar | True | True | True | False | True | False | True |
| SCN5A | 27401036 | `F1760A` | parser_dropped | pubtator | True | True | True | False | True | False | True |
| KCNH2 | 18774102 | `EX6-14DEL` | paper_row_lost_to_linkage_origin | figure,pubtator | True | True | True | True | True | False | True |
| KCNH2 | 18774102 | `EX9-14DUP` | projection_dropped | figure | True | True | True | True | True | False | False |

## Notation-unknown rows (probe could not search the notation)

| gene | PMID | variant | sweep class |
|---|---:|---|---|
| KCNH2 | 18452873 | `G189ins` | text_absent_figures_present |
| KCNQ1 | 22677073 | `L162fsX` | text_absent_figures_present |
| KCNH2 | 16969682 | `G601S` | text_absent_figures_present |
| KCNH2 | 12402336 | `G192fsX` | text_absent_notation_inconclusive |
| KCNH2 | 12402336 | `P872fsX` | text_absent_notation_inconclusive |
| KCNH2 | 12402336 | `R823fsX` | text_absent_notation_inconclusive |
| KCNH2 | 21185499 | `EXON1-14DEL` | text_absent_notation_inconclusive |
| KCNH2 | 21185499 | `EXON4-14DEL` | text_absent_notation_inconclusive |
| KCNH2 | 24103226 | `IVS9+5G>T` | text_absent_notation_inconclusive |
| KCNH2 | 24103226 | `P946fsX` | text_absent_notation_inconclusive |
| KCNH2 | 24103226 | `V1038fsX` | text_absent_notation_inconclusive |
| SCN5A | 17227473 | `c.3840+1G>A` | text_absent_notation_inconclusive |
| MYBPC3 | 21302287 | `c.2258_2259insT` | text_absent_figures_present |
| MYBPC3 | 21302287 | `c.3192_3193insC` | text_absent_figures_present |
| MYBPC3 | 21302287 | `c.506-12delC` | text_absent_figures_present |
| KCNH2 | 10790218 | `G189ins` | text_absent_notation_inconclusive |
| KCNH2 | 17905336 | `F627fsX` | text_absent_notation_inconclusive |
| KCNH2 | 17905336 | `G1031fsX` | text_absent_notation_inconclusive |
| KCNH2 | 17905336 | `L171ins` | text_absent_notation_inconclusive |
| KCNH2 | 17905336 | `R328fsX` | text_absent_notation_inconclusive |
| KCNQ1 | 17224687 | `R59RP` | text_absent_notation_inconclusive |
| RYR2 | 29447731 | `169-198_273+823DEL` | text_absent_notation_inconclusive |
| RYR2 | 29447731 | `1827+140_1961+426DEL` | text_absent_notation_inconclusive |
| RYR2 | 29447731 | `4299+1DELG` | text_absent_notation_inconclusive |
| BRCA1 | 10528853 | `c.1186A>G` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.172T>C` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.185delAG` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.2196G>A` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.2201C>T` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.2430T>C` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.2434T>C` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.2731C>T` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.3107A>T` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.3232A>G` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.3827T>G` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.5382insC` | text_absent_figures_present |
| BRCA1 | 10528853 | `c.546G>T` | text_absent_figures_present |
| SCN5A | 8917568 | `P.K1505_Q1507DEL` | text_absent_figures_present |
| SCN5A | 28018021 | `V1279I` | text_absent_figures_present |
| SCN5A | 17905336 | `c.4299_4300insG` | text_absent_notation_inconclusive |
