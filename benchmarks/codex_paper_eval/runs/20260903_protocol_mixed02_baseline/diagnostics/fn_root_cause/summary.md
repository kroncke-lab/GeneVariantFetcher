# FN root cause — `20260903_protocol_mixed02_baseline`

35 paper-derived false negatives.

| leaf | rows |
|---|---:|
| `unknown_notation` | 20 |
| `acquisition` | 13 |
| `projection_dropped` | 1 |
| `condensing` | 1 |

| gene | PMID | FN | leaves |
|---|---:|---:|---|
| SCN5A | 12106943 | 16 | unknown_notation=16 |
| SCN5A | 25650408 | 12 | acquisition=12 |
| SCN5A | 30059973 | 3 | unknown_notation=2, condensing=1 |
| BRCA1 | 18627636 | 1 | projection_dropped=1 |
| KCNH2 | 26022375 | 1 | unknown_notation=1 |
| KCNH2 | 14642688 | 1 | acquisition=1 |
| SCN5A | 21167004 | 1 | unknown_notation=1 |

## Rows the reading protocol could have found

| gene | PMID | variant | leaf | db layers | run_text | request | response | extraction | db | paper lane | linkage lane |
|---|---:|---|---|---|---|---|---|---|---|---|---|
| BRCA1 | 18627636 | `180delA` | projection_dropped | regex_table | True | True | True | True | True | False | False |
| SCN5A | 30059973 | `P.I1570DUP` | condensing |  | True | False | False | False | False | False | False |

## Notation-unknown rows (probe could not search the notation)

| gene | PMID | variant | sweep class |
|---|---:|---|---|
| SCN5A | 12106943 | `A1924T` | text_absent_figures_present |
| SCN5A | 12106943 | `D951X` | text_absent_figures_present |
| SCN5A | 12106943 | `E1225K` | text_absent_figures_present |
| SCN5A | 12106943 | `E161K` | text_absent_figures_present |
| SCN5A | 12106943 | `G1319V` | text_absent_figures_present |
| SCN5A | 12106943 | `G1406R` | text_absent_figures_present |
| SCN5A | 12106943 | `G1502S` | text_absent_figures_present |
| SCN5A | 12106943 | `G1743E` | text_absent_figures_present |
| SCN5A | 12106943 | `G752R` | text_absent_figures_present |
| SCN5A | 12106943 | `L867X` | text_absent_figures_present |
| SCN5A | 12106943 | `M369K` | text_absent_figures_present |
| SCN5A | 12106943 | `R1512W` | text_absent_figures_present |
| SCN5A | 12106943 | `R367C` | text_absent_figures_present |
| SCN5A | 12106943 | `R535X` | text_absent_figures_present |
| SCN5A | 12106943 | `S1382I` | text_absent_figures_present |
| SCN5A | 12106943 | `V1405L` | text_absent_figures_present |
| SCN5A | 30059973 | `P.G1015DFSX14` | text_absent_notation_inconclusive |
| SCN5A | 30059973 | `P.L1302VFS18` | text_absent_notation_inconclusive |
| KCNH2 | 26022375 | `T1083fsX` | text_absent_figures_present |
| SCN5A | 21167004 | `I1836T` | text_absent_figures_present |
