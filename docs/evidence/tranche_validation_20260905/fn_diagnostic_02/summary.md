# FN root cause — `20260905_protocol_cont120_02_candidate`

126 paper-derived primary false negatives.

| leaf | rows |
|---|---:|
| `acquisition` | 73 |
| `unknown_notation` | 31 |
| `model_missed` | 16 |
| `matcher` | 3 |
| `condensing` | 2 |
| `projection_dropped` | 1 |

| gene | PMID | FN | leaves |
|---|---:|---:|---|
| SCN5A | 22885917 | 22 | acquisition=14, unknown_notation=8 |
| KCNH2 | 29650123 | 21 | acquisition=15, unknown_notation=6 |
| SCN5A | 15996170 | 11 | acquisition=11 |
| KCNH2 | 11854117 | 11 | model_missed=6, unknown_notation=4, acquisition=1 |
| KCNH2 | 21216356 | 8 | acquisition=8 |
| SCN5A | 24144883 | 6 | acquisition=6 |
| RYR2 | 30471092 | 5 | unknown_notation=2, matcher=2, acquisition=1 |
| KCNH2 | 11844290 | 5 | acquisition=5 |
| SCN5A | 27566755 | 4 | unknown_notation=3, condensing=1 |
| KCNH2 | 9693036 | 4 | acquisition=4 |
| MYBPC3 | 20433692 | 3 | model_missed=3 |
| SCN5A | 21216356 | 2 | acquisition=2 |
| KCNH2 | 24057343 | 2 | acquisition=2 |
| KCNQ1 | 34135346 | 2 | unknown_notation=2 |
| SCN5A | 21321465 | 2 | model_missed=2 |
| SCN5A | 28491684 | 1 | unknown_notation=1 |
| KCNQ1 | 18464931 | 1 | condensing=1 |
| KCNH2 | 14676148 | 1 | unknown_notation=1 |
| KCNQ1 | 32830254 | 1 | acquisition=1 |
| KCNH2 | 15364333 | 1 | unknown_notation=1 |
| KCNH2 | 26496715 | 1 | unknown_notation=1 |
| KCNQ1 | 20368164 | 1 | model_missed=1 |
| KCNQ1 | 30244407 | 1 | model_missed=1 |
| SCN5A | 24363796 | 1 | unknown_notation=1 |
| SCN5A | 16039271 | 1 | model_missed=1 |
| KCNQ1 | 32405922 | 1 | unknown_notation=1 |
| SCN5A | 29672598 | 1 | matcher=1 |
| SCN5A | 26496715 | 1 | acquisition=1 |
| SCN5A | 30246897 | 1 | acquisition=1 |
| RYR2 | 18285261 | 1 | model_missed=1 |
| KCNQ1 | 26496715 | 1 | model_missed=1 |
| SCN5A | 25236808 | 1 | projection_dropped=1 |
| KCNQ1 | 18441444 | 1 | acquisition=1 |

## Rows the reading protocol could have found

| gene | PMID | variant | leaf | db layers | run_text | request | response | extraction | db | paper lane | linkage lane |
|---|---:|---|---|---|---|---|---|---|---|---|---|
| KCNQ1 | 18464931 | `Y315S` | condensing |  | True | False | False | False | False | False | False |
| SCN5A | 27566755 | `P.K1505_Q1507DEL` | condensing | regex_table | True | False | False | True | True | False | False |
| MYBPC3 | 20433692 | `IVS11-9G>A` | model_missed |  | True | True | False | False | False | False | False |
| MYBPC3 | 20433692 | `IVS29+5G>A` | model_missed |  | True | True | False | False | False | False | False |
| MYBPC3 | 20433692 | `IVS6+5G>A` | model_missed |  | True | True | False | False | False | False | False |
| KCNQ1 | 20368164 | `A340E` | model_missed | figure | True | True | False | False | True | False | False |
| KCNQ1 | 30244407 | `F617V` | model_missed |  | True | True | False | False | False | False | False |
| SCN5A | 16039271 | `A1330P` | model_missed | pubtator | True | True | False | False | True | False | True |
| RYR2 | 30471092 | `Q2060H` | matcher | llm_table | True | True | True | True | True | True | False |
| RYR2 | 30471092 | `Q293P` | matcher | llm_table | True | True | True | True | True | True | False |
| SCN5A | 21321465 | `P1090L` | model_missed |  | True | True | False | False | False | False | False |
| SCN5A | 21321465 | `R1193Q` | model_missed |  | True | True | False | False | False | False | False |
| SCN5A | 29672598 | `F2004L` | matcher | llm_table | True | True | True | True | True | True | False |
| RYR2 | 18285261 | `K4481R` | model_missed | regex_text | True | True | False | True | True | False | False |
| KCNQ1 | 26496715 | `360_361DUPKQ` | model_missed |  | True | True | False | False | False | False | False |
| SCN5A | 25236808 | `P.F1617DEL` | projection_dropped | llm_text | True | True | True | True | True | False | False |
| KCNH2 | 11854117 | `C44X` | model_missed |  | True | True | False | False | False | False | False |
| KCNH2 | 11854117 | `L799SP` | model_missed |  | True | True | False | False | False | False | False |
| KCNH2 | 11854117 | `Q376SP` | model_missed |  | True | True | False | False | False | False | False |
| KCNH2 | 11854117 | `R744X` | model_missed | clinvar | True | True | False | False | True | False | True |
| KCNH2 | 11854117 | `S428X` | model_missed |  | True | True | False | False | False | False | False |
| KCNH2 | 11854117 | `W1001X` | model_missed | clinvar | True | True | False | False | True | False | True |

## Notation-unknown rows (probe could not search the notation)

| gene | PMID | variant | sweep class |
|---|---:|---|---|
| SCN5A | 28491684 | `EXON23_DELETION` | text_absent_figures_present |
| KCNH2 | 14676148 | `N588K` | text_absent_figures_present |
| SCN5A | 27566755 | `P.F1617DEL` | text_absent_figures_present |
| SCN5A | 27566755 | `P.I1762DEL` | text_absent_figures_present |
| SCN5A | 27566755 | `P.Q1507_P1509DEL` | text_absent_figures_present |
| KCNH2 | 15364333 | `c.1945+6T>C` | text_absent_notation_inconclusive |
| KCNH2 | 26496715 | `T443fsX` | text_absent_notation_inconclusive |
| KCNQ1 | 34135346 | `G148R` | text_absent_figures_present |
| KCNQ1 | 34135346 | `R387Q` | text_absent_figures_present |
| SCN5A | 24363796 | `c.5445_5446insT` | text_absent_figures_present |
| KCNQ1 | 32405922 | `V416M` | text_absent_figures_present |
| KCNH2 | 29650123 | `F617fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `G911fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `L109fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `L987fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `R1035fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `R892fsX` | text_absent_notation_inconclusive |
| RYR2 | 30471092 | `169-?_273+?DEL` | text_absent_notation_inconclusive |
| RYR2 | 30471092 | `c.169-?_c.273+?del;` | text_absent_notation_inconclusive |
| SCN5A | 22885917 | `P.L729DEL` | text_absent_notation_inconclusive |
| SCN5A | 22885917 | `c.1570_1571insG` | text_absent_notation_inconclusive |
| SCN5A | 22885917 | `c.2582_2583delTT` | text_absent_notation_inconclusive |
| SCN5A | 22885917 | `c.3840+1G>A` | text_absent_notation_inconclusive |
| SCN5A | 22885917 | `c.4118delT` | text_absent_notation_inconclusive |
| SCN5A | 22885917 | `c.5280delG` | text_absent_notation_inconclusive |
| SCN5A | 22885917 | `c.704-1G>C` | text_absent_notation_inconclusive |
| SCN5A | 22885917 | `c.934+1G>A` | text_absent_notation_inconclusive |
| KCNH2 | 11854117 | `A83fsX` | text_absent_notation_inconclusive |
| KCNH2 | 11854117 | `G925fsX` | text_absent_notation_inconclusive |
| KCNH2 | 11854117 | `P968fsX` | text_absent_notation_inconclusive |
| KCNH2 | 11854117 | `V295fsX` | text_absent_notation_inconclusive |
