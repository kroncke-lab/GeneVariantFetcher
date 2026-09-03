# FN root cause — `20260903_protocol_mixed01_candidate`

85 paper-derived false negatives.

| leaf | rows |
|---|---:|
| `acquisition` | 70 |
| `unknown_notation` | 12 |
| `matcher` | 2 |
| `condensing` | 1 |

| gene | PMID | FN | leaves |
|---|---:|---:|---|
| KCNQ1 | 14678125 | 41 | acquisition=37, unknown_notation=4 |
| KCNQ1 | 21131640 | 9 | acquisition=8, unknown_notation=1 |
| KCNH2 | 15176425 | 7 | acquisition=7 |
| KCNQ1 | 15176425 | 6 | acquisition=6 |
| SCN5A | 15176425 | 5 | acquisition=5 |
| SCN5A | 12566525 | 4 | acquisition=3, unknown_notation=1 |
| SCN5A | 27232914 | 2 | acquisition=1, unknown_notation=1 |
| SCN5A | 16764707 | 2 | unknown_notation=2 |
| RYR2 | 17556193 | 2 | matcher=2 |
| SCN5A | 11804990 | 2 | unknown_notation=2 |
| SCN5A | 28779003 | 1 | acquisition=1 |
| SCN5A | 17971661 | 1 | acquisition=1 |
| KCNH2 | 16470702 | 1 | unknown_notation=1 |
| SCN5A | 12574983 | 1 | acquisition=1 |
| MYBPC3 | 33673806 | 1 | condensing=1 |

## Rows the reading protocol could have found

| gene | PMID | variant | leaf | db layers | run_text | request | response | extraction | db | paper lane | linkage lane |
|---|---:|---|---|---|---|---|---|---|---|---|---|
| MYBPC3 | 33673806 | `c.3408C>A` | condensing | regex_table | True | False | False | True | True | True | False |
| RYR2 | 17556193 | `R2267H` | matcher | figure,llm_table,pubtator | True | True | True | True | True | True | True |
| RYR2 | 17556193 | `S4565R` | matcher | figure,llm_table,pubtator | True | True | True | True | True | True | True |

## Notation-unknown rows (probe could not search the notation)

| gene | PMID | variant | sweep class |
|---|---:|---|---|
| SCN5A | 27232914 | `c.3391-1G>A` | text_absent_notation_inconclusive |
| SCN5A | 12566525 | `P.Y1795_E1796INSD` | text_absent_notation_inconclusive |
| KCNQ1 | 14678125 | `A344/SP` | text_absent_notation_inconclusive |
| KCNQ1 | 14678125 | `A344A/SPLICE` | text_absent_notation_inconclusive |
| KCNQ1 | 14678125 | `DELF340` | text_absent_notation_inconclusive |
| KCNQ1 | 14678125 | `L151` | text_absent_notation_inconclusive |
| KCNQ1 | 21131640 | `S333CFS/129X` | text_absent_notation_inconclusive |
| KCNH2 | 16470702 | `7Q34Q36.2DEL` | text_absent_notation_inconclusive |
| SCN5A | 16764707 | `P.1570_F1571INSI` | text_absent_notation_inconclusive |
| SCN5A | 16764707 | `P.Y1795_E1796INSD` | text_absent_notation_inconclusive |
| SCN5A | 11804990 | `D1595N` | text_absent_figures_present |
| SCN5A | 11804990 | `G298S` | text_absent_figures_present |
