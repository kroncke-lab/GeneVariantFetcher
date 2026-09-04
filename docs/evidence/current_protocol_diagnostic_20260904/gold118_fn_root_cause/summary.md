# FN root cause — `20260902_false_zero_recovery_gold118`

89 legacy trusted-projection false negatives.

| leaf | rows |
|---|---:|
| `acquisition` | 39 |
| `unknown_notation` | 27 |
| `model_missed` | 7 |
| `parser_dropped` | 5 |
| `source_selection` | 4 |
| `projection_dropped` | 4 |
| `condensing` | 2 |
| `matcher` | 1 |

| gene | PMID | FN | leaves |
|---|---:|---:|---|
| KCNH2 | 29650123 | 21 | acquisition=15, unknown_notation=6 |
| SCN5A | 15898185 | 10 | acquisition=9, parser_dropped=1 |
| RYR2 | 27114410 | 9 | acquisition=9 |
| KCNQ1 | 31520628 | 7 | unknown_notation=7 |
| KCNQ1 | 14678125 | 6 | unknown_notation=4, acquisition=2 |
| RYR2 | 19926015 | 5 | unknown_notation=2, model_missed=2, projection_dropped=1 |
| RYR2 | 22677073 | 4 | model_missed=3, projection_dropped=1 |
| KCNH2 | 29121719 | 3 | parser_dropped=3 |
| KCNQ1 | 24667783 | 3 | unknown_notation=3 |
| RYR2 | 29350269 | 3 | source_selection=3 |
| KCNQ1 | 19687231 | 2 | model_missed=2 |
| SCN5A | 22966897 | 2 | acquisition=2 |
| KCNH2 | 10086971 | 1 | unknown_notation=1 |
| KCNH2 | 16155735 | 1 | unknown_notation=1 |
| KCNH2 | 19184172 | 1 | unknown_notation=1 |
| KCNH2 | 29214556 | 1 | source_selection=1 |
| KCNQ1 | 20368164 | 1 | parser_dropped=1 |
| KCNQ1 | 26496715 | 1 | condensing=1 |
| KCNQ1 | 27114410 | 1 | acquisition=1 |
| KCNQ1 | 33082985 | 1 | unknown_notation=1 |
| RYR2 | 18285261 | 1 | projection_dropped=1 |
| RYR2 | 19216760 | 1 | matcher=1 |
| RYR2 | 30835254 | 1 | condensing=1 |
| SCN5A | 17971661 | 1 | acquisition=1 |
| SCN5A | 25236808 | 1 | projection_dropped=1 |
| SCN5A | 28087622 | 1 | unknown_notation=1 |

## Rows the reading protocol could have found

| gene | PMID | variant | leaf | db layers | run_text | request | response | extraction | db | scored lane | linkage lane |
|---|---:|---|---|---|---|---|---|---|---|---|---|
| KCNH2 | 29121719 | `K897T` | parser_dropped |  | True | True | True | False | False | False | False |
| KCNH2 | 29121719 | `M645R` | parser_dropped |  | True | True | True | False | False | False | False |
| KCNH2 | 29121719 | `R92C` | parser_dropped |  | True | True | True | False | False | False | False |
| KCNH2 | 29214556 | `E876X` | source_selection |  | False | False | False | False | False | False | False |
| KCNQ1 | 19687231 | `S338C` | model_missed |  | True | True | False | False | False | False | False |
| KCNQ1 | 19687231 | `F340C` | model_missed |  | True | True | False | False | False | False | False |
| KCNQ1 | 20368164 | `A340E` | parser_dropped | figure | True | True | True | False | True | False | False |
| KCNQ1 | 26496715 | `360_361DUPKQ` | condensing | regex_table | True | False | False | True | True | False | False |
| RYR2 | 18285261 | `K4481R` | projection_dropped | regex_text | True | True | True | True | True | False | False |
| RYR2 | 19216760 | `EXON 3 DELETION` | matcher | figure | True | True | True | True | True | True | False |
| RYR2 | 19926015 | `EXON 3 DELETION` | projection_dropped | llm_table | True | True | True | True | True | False | False |
| RYR2 | 19926015 | `N4736 DEL` | model_missed |  | True | True | False | False | False | False | False |
| RYR2 | 19926015 | `R4822H` | model_missed |  | True | True | False | False | False | False | False |
| RYR2 | 22677073 | `G4936R` | projection_dropped | regex_table | True | True | True | True | True | False | False |
| RYR2 | 22677073 | `Q2958R` | model_missed |  | True | True | False | False | False | False | False |
| RYR2 | 22677073 | `A1136V` | model_missed |  | True | True | False | False | False | False | False |
| RYR2 | 22677073 | `R4037C` | model_missed |  | True | True | False | False | False | False | False |
| RYR2 | 29350269 | `Q4164E` | source_selection |  | False | False | False | False | False | False | False |
| RYR2 | 29350269 | `E1127G` | source_selection |  | False | False | False | False | False | False | False |
| RYR2 | 29350269 | `A3814V` | source_selection |  | False | False | False | False | False | False | False |
| RYR2 | 30835254 | `P1124L` | condensing |  | True | False | False | False | False | False | False |
| SCN5A | 15898185 | `S1103Y` | parser_dropped |  | True | True | True | False | False | False | False |
| SCN5A | 25236808 | `P.F1617DEL` | projection_dropped | llm_text | True | True | True | True | True | False | False |

## Notation-unknown rows (probe could not search the notation)

| gene | PMID | variant | sweep class |
|---|---:|---|---|
| KCNH2 | 10086971 | `C1117fsX` | text_absent_figures_present |
| KCNH2 | 16155735 | `P923fsX` | text_absent_notation_inconclusive |
| KCNH2 | 19184172 | `I782fsX` | text_absent_figures_present |
| KCNH2 | 29650123 | `F617fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `G911fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `L109fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `L987fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `R1035fsX` | text_absent_notation_inconclusive |
| KCNH2 | 29650123 | `R892fsX` | text_absent_notation_inconclusive |
| KCNQ1 | 14678125 | `L151` | text_absent_notation_inconclusive |
| KCNQ1 | 14678125 | `DELF340` | text_absent_notation_inconclusive |
| KCNQ1 | 14678125 | `A344/SP` | text_absent_notation_inconclusive |
| KCNQ1 | 14678125 | `A344A/SPLICE` | text_absent_notation_inconclusive |
| KCNQ1 | 24667783 | `I328_S330DEL` | text_absent_figures_present |
| KCNQ1 | 24667783 | `D454TFS*7` | text_absent_figures_present |
| KCNQ1 | 24667783 | `Y522X` | text_absent_figures_present |
| KCNQ1 | 31520628 | `G589D` | text_absent_figures_present |
| KCNQ1 | 31520628 | `R366W` | text_absent_figures_present |
| KCNQ1 | 31520628 | `Y184S` | text_absent_figures_present |
| KCNQ1 | 31520628 | `Y315N` | text_absent_figures_present |
| KCNQ1 | 31520628 | `V254M` | text_absent_figures_present |
| KCNQ1 | 31520628 | `I235N` | text_absent_figures_present |
| KCNQ1 | 31520628 | `R252fsX` | text_absent_figures_present |
| KCNQ1 | 33082985 | `W248F + L347R` | text_absent_figures_present |
| RYR2 | 19926015 | `L62F` | text_absent_figures_present |
| RYR2 | 19926015 | `M81L` | text_absent_figures_present |
| SCN5A | 28087622 | `P.K1505_Q1507DEL` | text_absent_figures_present |
