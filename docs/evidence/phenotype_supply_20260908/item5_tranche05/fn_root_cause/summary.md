# FN root cause — `20260908_protocol_cont120_05_baseline`

158 paper-derived primary false negatives.

| leaf | rows |
|---|---:|
| `acquisition` | 84 |
| `unknown_notation` | 57 |
| `model_missed` | 15 |
| `parser_dropped` | 2 |

| gene | PMID | FN | leaves |
|---|---:|---:|---|
| RYR2 | 27452199 | 34 | unknown_notation=34 |
| SCN5A | 24721456 | 14 | acquisition=14 |
| SCN5A | 22360817 | 13 | acquisition=13 |
| KCNQ1 | 32470535 | 12 | acquisition=9, unknown_notation=3 |
| SCN5A | 21126620 | 7 | acquisition=7 |
| SCN5A | 26921764 | 6 | unknown_notation=6 |
| SCN5A | 18508782 | 5 | acquisition=5 |
| KCNH2 | 10973849 | 5 | model_missed=4, acquisition=1 |
| SCN5A | 15992732 | 5 | model_missed=5 |
| SCN5A | 16344400 | 4 | acquisition=4 |
| KCNH2 | 18508782 | 4 | acquisition=4 |
| KCNQ1 | 18595190 | 4 | acquisition=4 |
| SCN5A | 25370050 | 4 | model_missed=4 |
| KCNQ1 | 15028050 | 4 | acquisition=4 |
| KCNQ1 | 27041096 | 3 | acquisition=3 |
| SCN5A | 12569154 | 2 | unknown_notation=2 |
| SCN5A | 15277732 | 2 | acquisition=2 |
| SCN5A | 11535580 | 2 | acquisition=2 |
| SCN5A | 10973849 | 2 | unknown_notation=2 |
| SCN5A | 24269159 | 2 | acquisition=2 |
| KCNH2 | 19184172 | 2 | unknown_notation=2 |
| SCN5A | 14716629 | 1 | model_missed=1 |
| SCN5A | 23425522 | 1 | acquisition=1 |
| RYR2 | 18554199 | 1 | acquisition=1 |
| RYR2 | 22334434 | 1 | acquisition=1 |
| KCNQ1 | 27525866 | 1 | acquisition=1 |
| SCN5A | 18362431 | 1 | acquisition=1 |
| KCNQ1 | 28720088 | 1 | unknown_notation=1 |
| KCNQ1 | 20044973 | 1 | acquisition=1 |
| KCNH2 | 30105468 | 1 | unknown_notation=1 |
| SCN5A | 29132927 | 1 | acquisition=1 |
| KCNQ1 | 26019114 | 1 | acquisition=1 |
| KCNQ1 | 19808498 | 1 | model_missed=1 |
| RYR2 | 29668588 | 1 | unknown_notation=1 |
| BRCA1 | 30309222 | 1 | parser_dropped=1 |
| KCNH2 | 15028050 | 1 | parser_dropped=1 |
| KCNH2 | 22727609 | 1 | unknown_notation=1 |
| SCN5A | 12123759 | 1 | acquisition=1 |
| SCN5A | 12193783 | 1 | unknown_notation=1 |
| RYR2 | 34735682 | 1 | unknown_notation=1 |
| SCN5A | 24895456 | 1 | unknown_notation=1 |
| KCNH2 | 30844837 | 1 | unknown_notation=1 |
| SCN5A | 18071069 | 1 | acquisition=1 |

## Rows the reading protocol could have found

| gene | PMID | variant | leaf | db layers | run_text | request | response | extraction | db | paper lane | linkage lane |
|---|---:|---|---|---|---|---|---|---|---|---|---|
| SCN5A | 14716629 | `G1406R` | model_missed | pubtator | True | True | False | False | True | False | True |
| KCNH2 | 10973849 | `D864SP` | model_missed | clinvar | True | True | False | False | True | False | True |
| KCNH2 | 10973849 | `L799SP` | model_missed |  | True | True | False | False | False | False | False |
| KCNH2 | 10973849 | `Q376SP` | model_missed | clinvar | True | True | False | False | True | False | True |
| KCNH2 | 10973849 | `Q376SP` | model_missed | clinvar | True | True | False | False | True | False | True |
| SCN5A | 15992732 | `P1090L` | model_missed | pubtator | True | True | False | False | True | False | True |
| SCN5A | 15992732 | `R1193Q` | model_missed | pubtator | True | True | False | False | True | False | True |
| SCN5A | 15992732 | `R34C` | model_missed | clinvar,pubtator | True | True | False | False | True | False | True |
| SCN5A | 15992732 | `S524Y` | model_missed | clinvar,pubtator | True | True | False | False | True | False | True |
| SCN5A | 15992732 | `V1951L` | model_missed | pubtator | True | True | False | False | True | False | True |
| KCNQ1 | 19808498 | `F296S` | model_missed | clinvar | True | True | False | False | True | False | True |
| BRCA1 | 30309222 | `c.922_923AGCinsT` | parser_dropped |  | True | True | True | False | False | False | False |
| KCNH2 | 15028050 | `S706C` | parser_dropped | clinvar,pubtator | True | True | True | False | True | False | True |
| SCN5A | 25370050 | `K1922E` | model_missed |  | True | True | False | False | False | False | False |
| SCN5A | 25370050 | `R1910E` | model_missed |  | True | True | False | False | False | False | False |
| SCN5A | 25370050 | `R1914A` | model_missed |  | True | True | False | False | False | False | False |
| SCN5A | 25370050 | `R1914E` | model_missed |  | True | True | False | False | False | False | False |

## Notation-unknown rows (probe could not search the notation)

| gene | PMID | variant | sweep class |
|---|---:|---|---|
| SCN5A | 26921764 | `P.Y1795_E1796INSD` | text_absent_notation_inconclusive |
| SCN5A | 26921764 | `c.2559delT` | text_absent_notation_inconclusive |
| SCN5A | 26921764 | `c.2582_2583delTT` | text_absent_notation_inconclusive |
| SCN5A | 26921764 | `c.3667delG` | text_absent_notation_inconclusive |
| SCN5A | 26921764 | `c.3840+1G>A` | text_absent_notation_inconclusive |
| SCN5A | 26921764 | `c.934+1G>A` | text_absent_notation_inconclusive |
| SCN5A | 12569154 | `H558R` | text_absent_figures_present |
| SCN5A | 12569154 | `T512I` | text_absent_figures_present |
| KCNQ1 | 28720088 | `R519X` | text_absent_figures_present |
| KCNH2 | 30105468 | `W585fsX` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `11200 C>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `11590 A>G` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `11812 A>G` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `11812 A>G` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `11836 G>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `11836 G>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `11997 G>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `12066_12067INS CAA G5656A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `1221 A>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `12301 C>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `12331 A>C` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `1244 C>G` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `12475 C>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `12513 G>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `12559 G>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `1258 C>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `1259 G>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `13175 A>G` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `13489 C>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `13759 G>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `14427 G>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `14461 G>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `14593 C>G` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `40-2 A/G SPLICING ERROR` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `4652 A>G` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `5170 G>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `535 G>A` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `6507 G>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `6643 T>C` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `6649 C>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `6737 C>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `6737 C>T` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `718 T>C` | text_absent_notation_inconclusive |
| RYR2 | 27452199 | `7204 T>C` | text_absent_notation_inconclusive |
| KCNQ1 | 32470535 | `E449fsX` | text_absent_notation_inconclusive |
| KCNQ1 | 32470535 | `G245RFSX` | text_absent_notation_inconclusive |
| KCNQ1 | 32470535 | `R518X + G833R` | text_absent_notation_inconclusive |
| RYR2 | 29668588 | `7580 T>C` | text_absent_figures_present |
| SCN5A | 10973849 | `P.K1505_Q1507DEL` | text_absent_figures_present |
| SCN5A | 10973849 | `P.Y1795_E1796INSD` | text_absent_figures_present |
| KCNH2 | 22727609 | `F640del` | text_absent_notation_inconclusive |
| KCNH2 | 19184172 | `E788K` | text_absent_figures_present |
| KCNH2 | 19184172 | `R537W` | text_absent_figures_present |
| SCN5A | 12193783 | `S1103Y` | text_absent_figures_present |
| RYR2 | 34735682 | `G2412R` | text_absent_figures_present |
| SCN5A | 24895456 | `V2016M` | text_absent_figures_present |
| KCNH2 | 30844837 | `D609G` | text_absent_figures_present |
