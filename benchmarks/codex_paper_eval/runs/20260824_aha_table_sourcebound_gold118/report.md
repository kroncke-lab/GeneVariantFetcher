# Codex extraction-blinded paper evaluation — `20260824_aha_table_sourcebound_gold118`

## Technical summary

This hash-locked run evaluated **118 papers** (**per-gene counts {'SCN5A': 30, 'KCNH2': 28, 'KCNQ1': 30, 'RYR2': 30}**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **66.2%**, recall **87.7%**, F1 **75.4%** (554 TP, 283 FP, 78 FN).
- Precision versus counted extras **97.9%** (554 matched rows; 12 extra rows with patient counts). The stricter count-bearing-only diagnostic is **94.6%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **2,730,460 total tokens** (1,951,910 input; 778,550 output).
- Elapsed: **8933.0s wall clock**; 8189.8s summed per-paper route + read time.
- Notation twins merged before scoring: **1** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 118}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `tier2_gold_120.tsv` (118 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility and count-field presence during selection; extraction exported no gold values or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 207 / 632 | 32.8% | 0.193 | 0.742 |
| affected | 56 / 632 | 8.9% | 0.321 | 1.102 |
| unaffected | 28 / 631 | 4.4% | 0.321 | 0.779 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 204 / 472 | 43.2% | 3 / 160 | 1.9% |
| affected | 49 / 406 | 12.1% | 7 / 226 | 3.1% |
| unaffected | 19 / 82 | 23.2% | 9 / 549 | 1.6% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 183 | 110 | 16 | 62.5% | 92.0% | 74.4% | 97.3% | 90.7% | 24.6% / 0.306 / 1.187 | 8.0% / 0.000 / 0.000 | 7.0% / 0.143 / 0.378 |
| KCNH2 | 66 | 49 | 22 | 57.4% | 75.0% | 65.0% | 97.1% | 93.3% | 29.5% / 0.154 / 0.620 | 17.0% / 0.933 / 2.033 | 3.4% / 1.667 / 2.082 |
| KCNQ1 | 146 | 101 | 21 | 59.1% | 87.4% | 70.5% | 96.7% | 92.8% | 38.3% / 0.156 / 0.559 | 2.4% / 0.250 / 0.500 | 0.6% / 0.000 / 0.000 |
| RYR2 | 159 | 23 | 19 | 87.4% | 89.3% | 88.3% | 100.0% | 100.0% | 38.2% / 0.162 / 0.470 | 11.8% / 0.143 / 0.488 | 5.6% / 0.200 / 0.447 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| KCNH2 | 10086971 | text | 3 | 4 | 1 | 42.9% | 75.0% | 54.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 108.9 | 38,198 |
| KCNH2 | 15500450 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 25.8 | 11,594 |
| KCNH2 | 16029385 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 84.0 | 22,473 |
| KCNH2 | 16155735 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 35.8 | 17,188 |
| KCNH2 | 18675227 | text | 3 | 1 | 0 | 75.0% | 100.0% | 85.7% | 66.7% / 0.000 | 66.7% / 1.000 | 0.0% / n/a | 69.0 | 26,422 |
| KCNH2 | 18808722 | text | 8 | 1 | 0 | 88.9% | 100.0% | 94.1% | 25.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 151.8 | 47,767 |
| KCNH2 | 19034806 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 48.8 | 10,791 |
| KCNH2 | 19065538 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 38.2 | 10,073 |
| KCNH2 | 19184172 | text | 3 | 2 | 1 | 60.0% | 75.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 257.2 | 35,219 |
| KCNH2 | 20181576 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 50.0% / 0.000 | 0.0% / n/a | 50.0% / 2.000 | 368.7 | 29,233 |
| KCNH2 | 21308345 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 49.0 | 18,680 |
| KCNH2 | 21779290 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 98.3 | 24,407 |
| KCNH2 | 22052944 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 17.5 | 15,949 |
| KCNH2 | 22067087 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 27.5 | 9,916 |
| KCNH2 | 22314138 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.7 | 16,704 |
| KCNH2 | 22338672 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 50.0% / 0.000 | 100.0% / 3.500 | 50.0% / 0.000 | 144.5 | 27,710 |
| KCNH2 | 22764740 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.2 | 21,130 |
| KCNH2 | 23917959 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.6 | 9,925 |
| KCNH2 | 23936059 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 85.9 | 25,896 |
| KCNH2 | 25819988 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 3.000 | 100.0% / 3.000 | 40.9 | 15,509 |
| KCNH2 | 26118593 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 52.3 | 18,357 |
| KCNH2 | 26746457 | text | 10 | 32 | 0 | 23.8% | 100.0% | 38.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 5.7 | 1,143 |
| KCNH2 | 29016797 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 50.0% / 1.000 | 0.0% / n/a | 93.8 | 25,393 |
| KCNH2 | 29121719 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.2 | 12,373 |
| KCNH2 | 29214556 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.2 | 8,465 |
| KCNH2 | 29650123 | text | 2 | 0 | 20 | 100.0% | 9.1% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 16.7 | 19,065 |
| KCNH2 | 30036649 | text | 5 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 80.0% / 0.000 | 0.0% / n/a | 77.7 | 31,542 |
| KCNH2 | 30246897 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.333 | 0.0% / n/a | 45.4 | 24,152 |
| KCNQ1 | 11351021 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.7 | 13,725 |
| KCNQ1 | 14678125 | text | 35 | 9 | 6 | 79.5% | 85.4% | 82.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 6.3 | 7,718 |
| KCNQ1 | 16155735 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 23.1 | 16,607 |
| KCNQ1 | 18567635 | text | 1 | 9 | 0 | 10.0% | 100.0% | 18.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 107.6 | 35,124 |
| KCNQ1 | 18580685 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 27.8 | 15,010 |
| KCNQ1 | 18808722 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 64.0 | 24,735 |
| KCNQ1 | 19056345 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 30.7 | 13,729 |
| KCNQ1 | 19114714 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.9 | 15,899 |
| KCNQ1 | 19632626 | text | 1 | 63 | 0 | 1.6% | 100.0% | 3.1% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 97.4 | 23,550 |
| KCNQ1 | 19687231 | text | 0 | 1 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.3 | 15,128 |
| KCNQ1 | 20348026 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 47.3 | 20,383 |
| KCNQ1 | 20368164 | text | 1 | 2 | 1 | 33.3% | 50.0% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.9 | 20,501 |
| KCNQ1 | 20959120 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.7 | 16,245 |
| KCNQ1 | 21070882 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 105.1 | 17,689 |
| KCNQ1 | 21129503 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 36.7 | 16,156 |
| KCNQ1 | 21956039 | text | 16 | 0 | 0 | 100.0% | 100.0% | 100.0% | 12.5% / 0.000 | 0.0% / n/a | 0.0% / n/a | 259.2 | 78,082 |
| KCNQ1 | 22613981 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.6 | 6,018 |
| KCNQ1 | 24667783 | text | 14 | 1 | 2 | 93.3% | 87.5% | 90.3% | 62.5% / 0.000 | 6.2% / 0.000 | 0.0% / n/a | 183.3 | 65,016 |
| KCNQ1 | 25471708 | text | 8 | 1 | 0 | 88.9% | 100.0% | 94.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 60.2 | 20,079 |
| KCNQ1 | 26496715 | text | 41 | 6 | 1 | 87.2% | 97.6% | 92.1% | 97.6% / 0.000 | 0.0% / n/a | 0.0% / n/a | 2.1 | 1,559 |
| KCNQ1 | 27114410 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.9 | 5,458 |
| KCNQ1 | 28491547 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 3.000 | 100.0% / 0.000 | 0.0% / n/a | 86.5 | 18,669 |
| KCNQ1 | 28739325 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 34.9 | 24,234 |
| KCNQ1 | 29677589 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 44.4 | 17,243 |
| KCNQ1 | 29851656 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 31.9 | 10,869 |
| KCNQ1 | 29876285 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 64.6 | 23,765 |
| KCNQ1 | 31293497 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 103.1 | 31,696 |
| KCNQ1 | 31520628 | text | 2 | 2 | 7 | 50.0% | 22.2% | 30.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 45.6 | 12,652 |
| KCNQ1 | 32168391 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 65.1 | 24,364 |
| KCNQ1 | 33082985 | text | 2 | 1 | 1 | 66.7% | 66.7% | 66.7% | 66.7% / 1.000 | 33.3% / 1.000 | 33.3% / 0.000 | 91.3 | 30,988 |
| RYR2 | 15466642 | text | 8 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.250 | 100.0% / 0.250 | 0.0% / n/a | 106.1 | 38,926 |
| RYR2 | 16517285 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 27.5 | 8,945 |
| RYR2 | 17161793 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 27.9 | 10,115 |
| RYR2 | 18285261 | text | 4 | 0 | 1 | 100.0% | 80.0% | 88.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 12.9 | 14,514 |
| RYR2 | 19216760 | text | 4 | 0 | 1 | 100.0% | 80.0% | 88.9% | 40.0% / 0.500 | 0.0% / n/a | 0.0% / n/a | 86.6 | 38,877 |
| RYR2 | 19398417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 19.8 | 6,443 |
| RYR2 | 19398665 | text | 27 | 1 | 0 | 96.4% | 100.0% | 98.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 52.0 | 19,771 |
| RYR2 | 19926015 | text | 36 | 6 | 4 | 85.7% | 90.0% | 87.8% | 0.0% / n/a | 0.0% / n/a | 2.5% / 1.000 | 153.6 | 51,216 |
| RYR2 | 20395638 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 30.4 | 8,759 |
| RYR2 | 21652165 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 36.9 | 10,336 |
| RYR2 | 22178870 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 41.5 | 14,257 |
| RYR2 | 22222782 | text | 7 | 0 | 0 | 100.0% | 100.0% | 100.0% | 71.4% / 0.200 | 0.0% / n/a | 0.0% / n/a | 164.5 | 48,888 |
| RYR2 | 22677073 | text | 22 | 1 | 3 | 95.7% | 88.0% | 91.7% | 80.0% / 0.000 | 12.0% / 0.000 | 4.0% / 0.000 | 262.9 | 93,012 |
| RYR2 | 25435091 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 35.4 | 12,354 |
| RYR2 | 25463374 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 82.4 | 31,657 |
| RYR2 | 25500949 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 26.4 | 13,843 |
| RYR2 | 25814417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 19.4 | 20,254 |
| RYR2 | 27114410 | text | 0 | 0 | 9 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.5 | 5,506 |
| RYR2 | 27839804 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 37.7 | 14,205 |
| RYR2 | 28798025 | text | 6 | 1 | 0 | 85.7% | 100.0% | 92.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 67.3 | 42,145 |
| RYR2 | 29350269 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 49.9 | 30,755 |
| RYR2 | 30403697 | text | 21 | 0 | 0 | 100.0% | 100.0% | 100.0% | 61.9% / 0.231 | 23.8% / 0.000 | 15.0% / 0.333 | 281.8 | 91,546 |
| RYR2 | 30546600 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 0.000 | 69.3 | 19,614 |
| RYR2 | 30835254 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 32866913 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 100.0% / 0.000 | 25.2 | 9,458 |
| RYR2 | 33315912 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 35.8 | 23,769 |
| RYR2 | 33536282 | text | 1 | 10 | 0 | 9.1% | 100.0% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 52.5 | 33,059 |
| RYR2 | 33640691 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 16.9 | 9,210 |
| RYR2 | 34661651 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 75.0% / 0.333 | 0.0% / n/a | 0.0% / n/a | 55.5 | 30,393 |
| RYR2 | 35663620 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 37.9 | 14,671 |
| SCN5A | 15671429 | text | 4 | 9 | 0 | 30.8% | 100.0% | 47.1% | 25.0% / 0.000 | 25.0% / 0.000 | 0.0% / n/a | 140.3 | 30,573 |
| SCN5A | 15828879 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 43.1 | 11,963 |
| SCN5A | 15898185 | text | 1 | 0 | 10 | 100.0% | 9.1% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.1 | 6,005 |
| SCN5A | 16301357 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.3 | 7,034 |
| SCN5A | 16929919 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.4 | 7,638 |
| SCN5A | 17675083 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 66.2 | 30,883 |
| SCN5A | 17971661 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 5.4 | 4,680 |
| SCN5A | 19305639 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.4 | 10,888 |
| SCN5A | 19406494 | text | 1 | 8 | 0 | 11.1% | 100.0% | 20.0% | 100.0% / 1.000 | 0.0% / n/a | 100.0% / 1.000 | 90.6 | 14,259 |
| SCN5A | 19549036 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.0 | 15,692 |
| SCN5A | 20038812 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.4 | 20,091 |
| SCN5A | 20539757 | text | 11 | 13 | 0 | 45.8% | 100.0% | 62.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 243.1 | 42,317 |
| SCN5A | 21193062 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 35.7 | 20,515 |
| SCN5A | 21288276 | text | 0 | 2 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.3 | 12,557 |
| SCN5A | 22101522 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.2 | 18,462 |
| SCN5A | 22677073 | text | 12 | 3 | 0 | 80.0% | 100.0% | 88.9% | 91.7% / 0.091 | 16.7% / 0.000 | 8.3% / 0.000 | 194.2 | 70,799 |
| SCN5A | 22685113 | text | 10 | 7 | 0 | 58.8% | 100.0% | 74.1% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 162.8 | 39,524 |
| SCN5A | 22882672 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 7.000 | 0.0% / n/a | 0.0% / n/a | 38.2 | 10,642 |
| SCN5A | 22966897 | text | 0 | 3 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 26.0 | 13,148 |
| SCN5A | 23538678 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 41.4 | 12,312 |
| SCN5A | 24573164 | text | 20 | 8 | 0 | 71.4% | 100.0% | 83.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 107.1 | 27,524 |
| SCN5A | 24581105 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 31.1 | 10,775 |
| SCN5A | 25171853 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 1.000 | 50.0% / 0.000 | 50.0% / 1.000 | 60.4 | 15,384 |
| SCN5A | 25236808 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 42.4 | 14,123 |
| SCN5A | 26820365 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 49.3 | 33,463 |
| SCN5A | 28087622 | text | 6 | 5 | 1 | 54.5% | 85.7% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 193.2 | 26,462 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.6 | 14,448 |
| SCN5A | 29544605 | text | 6 | 1 | 0 | 85.7% | 100.0% | 92.3% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 159.3 | 53,349 |
| SCN5A | 29709101 | text | 11 | 3 | 0 | 78.6% | 100.0% | 88.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 246.6 | 74,300 |
| SCN5A | 32533946 | text | 83 | 44 | 0 | 65.4% | 100.0% | 79.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 165.7 | 75,987 |

## Errors and representation choices

### KCNH2 PMID 10086971

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: C1117fsX
- Extra predictions: c.3108insG, p.Asp1037fs c.3107dup, p.Gly965fs c.2892dup, p.Leu953fs c.2857dup

### KCNH2 PMID 15500450

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 16029385

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 16155735

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 18675227

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Thr1062Ile c.3185C>T
- Count disagreements: p.L955V c.2863C>G affected 1 vs 0 (error +1); p.R954C c.2860C>T affected 1 vs 0 (error +1)

### KCNH2 PMID 18808722

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A490P

### KCNH2 PMID 19034806

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 19065538

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 19184172

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: I782fsX
- Extra predictions: p.Arg176Trp c.526C>T, p.Glu788Asp c.2364G>C

### KCNH2 PMID 20181576

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: P926A
- Count disagreements: K897T unaffected 1 vs 3 (error -2)

### KCNH2 PMID 21308345

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.1600G>T

### KCNH2 PMID 21779290

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Leu564Leu, p.Tyr652Tyr

### KCNH2 PMID 22052944

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 22067087

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 22314138

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 22338672

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A1116V
- Count disagreements: K897T affected 7 vs 0 (error +7)

### KCNH2 PMID 22764740

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 23917959

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 23936059

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 25819988

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.His562Pro c.1685A>C
- Count disagreements: p.H562R c.1685A>G affected 6 vs 9 (error -3); p.H562R c.1685A>G unaffected 7 vs 4 (error +3)

### KCNH2 PMID 26118593

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.F513F c.1539C>T, p.I489I c.1467C>T, p.L564L c.1692A>G

### KCNH2 PMID 26746457

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A175D, A188V, A190T, A9V, D982N, G238S, G873S, G903R, K897T, L1105S, M291R, P1018A, P1030L, P347S, P967L, R1035Q, R1047L, R148W, R181Q, R252W, R326H, R328C, R394H, R885C, R894H, R912W, R948S, S981G, V1063L, V325M, W927L, Y403C

### KCNH2 PMID 29016797

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Asn588Lys c.1764C>G carriers 1 vs 0 (error +1); p.Ser631Ala c.1891T>G carriers 3 vs 0 (error +3); p.Asn588Lys c.1764C>G affected 1 vs 0 (error +1)

### KCNH2 PMID 29121719

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 29214556

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 29650123

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A561V, C49Y, F617fsX, G572S, G911fsX, L109fsX, L779P, N633S, Q1046X, R1035fsX, R148W, R176W, R328C, R534C, R744X, R892fsX, S660L, S818P, W412X, Y43C

### KCNH2 PMID 30036649

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 30246897

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Cys566Arg c.1696T>C affected 1 vs 0 (error +1)

### KCNQ1 PMID 11351021

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 14678125

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: L151, G185S, DELF340, A344/SP, A344A/SPLICE, R518X
- Extra predictions: p.Arg243His c.728G>A, p.Arg243Ser c.727C>A, p.Arg366Leu c.1097G>T, p.Arg366Pro c.1097G>C, p.Asp317Ala c.950A>C, p.Glu284del c.850_852del, p.Gly186Cys c.556G>T, p.Gly186Ser c.556G>A, p.Leu266Arg c.797T>G

### KCNQ1 PMID 16155735

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Pro448Arg c.1343C>G

### KCNQ1 PMID 18567635

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: D242A, E261A, K218A, R228A, R231A, R237A, R243A, R249A, R259A

### KCNQ1 PMID 18580685

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 18808722

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.2020insAG

### KCNQ1 PMID 19056345

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 19114714

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 19632626

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.1033-2del, c.1251+9C>A, c.1394-1G>T, c.1515-4G>A, c.477+5G>A, p.Ala178Thr c.532G>A, p.Ala341Glu c.1022C>A, p.Ala341Val c.1022C>T, p.Ala344= c.1032G>A, p.Ala344Val c.1031C>T, p.Arg190Gln c.569G>A, p.Arg190Leu c.569G>T, p.Arg190Trp c.568C>T, p.Arg243Cys c.727C>T, p.Arg243His c.728G>A, p.Arg380Lys c.1139G>A, p.Arg401fs c.1201dup, p.Arg507Gln c.1520G>A, p.Arg518Ter c.1552C>T, p.Arg539Trp c.1615C>T, p.Arg555Ser c.1663C>A, p.Arg583His c.1748G>A, p.Arg591His c.1772G>A, p.Arg594Gln c.1781G>A, p.Arg632fs c.1893dup, p.Asp202Asn c.604G>A, p.Cys214Tyr c.641G>A, p.Cys642Phe c.1925G>T, p.Gln107His c.321G>C, p.Gln359Ter c.1075C>T, p.Gln376Ter c.1126C>T, p.Gln530Ter c.1588C>T, p.Glu449fs c.1343dup, p.Gly168Arg c.502G>A, p.Gly269Ser c.805G>A, p.Gly306Arg c.916G>A, p.Gly314Ser c.940G>A, p.Gly325Arg c.973G>A, p.Gly626fs c.1875dup, p.His258Arg c.773A>G, p.Ile375fs c.1124_1127del, p.Leu266Pro c.797T>C, p.Leu266fs c.796del, p.Leu496fs c.1486_1487del, p.Lys362Arg c.1085A>G, p.Lys526Glu c.1576A>G, p.Met520Arg c.1559T>G, p.Pro343Ala c.1027C>G, p.Ser277del c.825CTC[1], p.Ser546Leu c.1637C>T, p.Ser566Phe c.1697C>T, p.Thr311Ile c.932C>T, p.Thr312Ile c.935C>T, p.Thr322Met c.965C>T, p.Trp15Ter c.45G>A, p.Trp305Leu c.914G>T, p.Trp305Ter c.914G>A, p.Tyr111Cys c.332A>G, p.Tyr184Ser c.551A>C, p.Val110Ile c.328G>A, p.Val215Met c.643G>A, p.Val254Met c.760G>A, p.Val280Gly c.839T>G
- Count disagreements: S209P carriers 7 vs 6 (error +1)

### KCNQ1 PMID 19687231

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: S338C, F340C
- Extra predictions: I328L

### KCNQ1 PMID 20348026

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 20368164

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A340E
- Extra predictions: p.Ala341Glu, p.Thr312Ile

### KCNQ1 PMID 20959120

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 21070882

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Thr96Arg carriers 1 vs 3 (error -2)

### KCNQ1 PMID 21129503

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 21956039

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 22613981

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 24667783

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: I328_S330DEL, Y522X
- Extra predictions: p.Thr224Met c.671C>T

### KCNQ1 PMID 25471708

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.572_576del

### KCNQ1 PMID 26496715

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: 360_361DUPKQ
- Extra predictions: c.1071_1076dupGAAGCA, c.1251+2T>C, c.477+5G>A, c.605-2G>A, c.921+1G>T, p.Gln359_Arg c.1067AGCAGA[3]

### KCNQ1 PMID 27114410

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G584S

### KCNQ1 PMID 28491547

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg519His c.1556G>A
- Count disagreements: p.Val141Met c.421G>A carriers 3 vs 6 (error -3)

### KCNQ1 PMID 28739325

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Asp242Asn c.724G>A carriers 2 vs 0 (error +2)

### KCNQ1 PMID 29677589

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 29851656

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 29876285

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Trp305Ter c.914G>A

### KCNQ1 PMID 31293497

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.R360_Q361insQKQR, p.Lys362_His c.1077_1088dup

### KCNQ1 PMID 31520628

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G589D, R366W, Y184S, Y315N, V254M, I235N, R252fsX
- Extra predictions: A341V, p.Arg591His c.1772G>A

### KCNQ1 PMID 32168391

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 33082985

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: W248F + L347R
- Extra predictions: p.Trp248Arg c.742T>C
- Count disagreements: p.Leu347Arg c.1040T>G carriers 2 vs 1 (error +1); p.Trp248Phe c.743_744delGGinsTC carriers 2 vs 1 (error +1); p.Trp248Phe c.743_744delGGinsTC affected 1 vs 0 (error +1)

### RYR2 PMID 15466642

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Ile4848Val carriers 6 vs 4 (error +2); p.Ile4848Val affected 6 vs 4 (error +2)

### RYR2 PMID 16517285

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 17161793

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 18285261

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: K4481R

### RYR2 PMID 19216760

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: EXON 3 DELETION
- Count disagreements: p.Asn57_Gly91del deletion of exon 3 (1.1 kb) carriers 3 vs 4 (error -1)

### RYR2 PMID 19398417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 19398665

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: G3946S

### RYR2 PMID 19926015

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: L62F, M81L, EXON 3 DELETION, N4736 DEL
- Extra predictions: G1886S, p.A4282V, p.M2389L, p.R169Q, p.V2321M, p.V4010M
- Count disagreements: p.Y4149S unaffected 1 vs 0 (error +1)

### RYR2 PMID 20395638

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 21652165

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 22178870

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 22222782

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.D1872N carriers 1 vs 2 (error -1)

### RYR2 PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: Q2958R, A1136V, R4037C
- Extra predictions: c.2398+5G>T

### RYR2 PMID 25435091

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 25463374

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R3570W, S4565R
- Count disagreements: H29D carriers 2 vs 4 (error -2)

### RYR2 PMID 25500949

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 25814417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 27114410

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: H240R, T4158P, R420W, S2246L, L3879P, Q3925E, S3959L, W4645R, G4936R

### RYR2 PMID 27839804

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H29D

### RYR2 PMID 28798025

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.34930+2A>G

### RYR2 PMID 29350269

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 30403697

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: R2028H carriers 2 vs 3 (error -1); Y4721C carriers 2 vs 3 (error -1); c.14091-11dupT carriers 1 vs 2 (error -1); Y4721C unaffected 1 vs 2 (error -1)

### RYR2 PMID 30546600

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg56Gln c.167G>A

### RYR2 PMID 30835254

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P1124L

### RYR2 PMID 32866913

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Thr85Ile c.254C>T carriers 3 vs 2 (error +1); p.Thr85Ile c.254C>T affected 2 vs 1 (error +1)

### RYR2 PMID 33315912

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 33536282

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.A4860G, p.D4112N, p.D4646A, p.I2075T, p.I3995V, p.I4855M, p.K4594R, p.Q4879H, p.S4938F, p.T4196I

### RYR2 PMID 33640691

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 34661651

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.S4168P c.12502T>C carriers 1 vs 2 (error -1)

### RYR2 PMID 35663620

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15671429

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.2550-2551insTG, c.2550_2551insTG, c.704-2A>G, p.Asp1274Asn c.3820G>A, p.Asp1594Asn c.4780G>A, p.Asp1594His c.4780G>C, p.Gly1572= c.4716C>T, p.Phe851fs c.2550_2551dup, p.Pro1331fs c.3992del

### SCN5A PMID 15828879

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Asn406Ser c.1217A>G

### SCN5A PMID 15898185

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F1293S, H558R, P1090L, R1512W, R34C, R481W, S1103Y, S1787N, S524Y, V1951L

### SCN5A PMID 16301357

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 16929919

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17675083

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg34Cys
- Count disagreements: H558R carriers 5 vs 1 (error +4)

### SCN5A PMID 17971661

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: E1784K

### SCN5A PMID 19305639

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 19406494

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Ala1949Pro c.5845G>C, p.Trp301* c.903G>A, p.Arg2011Cys c.6031C>T, p.Glu1866Ter c.5596G>T, p.Leu1785fs c.5353_5354del, p.Phe2003Leu c.6007T>C, p.Pro1840fs c.5517dup, p.Trp1712Ter c.5135G>A
- Count disagreements: p.Arg808Cys c.2422C>T carriers 4 vs 3 (error +1); p.Arg808Cys c.2422C>T unaffected 1 vs 0 (error +1)

### SCN5A PMID 19549036

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 20038812

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Phe1485Leu c.4453T>C

### SCN5A PMID 20539757

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H558R, R1632X, T1620M, p.Lys1578fs, p.Phe1617del, c.475_482+3del, p.Arg1622Ter c.4864C>T, p.Arg1631His c.4892G>A, p.Arg1631Ser c.4891C>A, p.Asp1274Asn c.3820G>A, p.Gly1407Arg c.4219G>A, p.Phe1616del c.4844TCT[1], p.Pro1297Leu c.3890C>T

### SCN5A PMID 21193062

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 21288276

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.999-424_1338+81del
- Extra predictions: deletion of exons 9 and 10, p.Trp193fs c.576del

### SCN5A PMID 22101522

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Ile1767Val c.5299A>G, p.Pro2005Ala c.6013C>G, p.Val1950Leu c.5848G>T
- Count disagreements: S216L carriers 2 vs 1 (error +1)

### SCN5A PMID 22685113

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg1625His c.4874G>A, p.Arg1896Trp c.5686C>T, p.Asp1818Asn c.5452G>A, p.Phe1595Ile c.4783T>A, p.Phe2003Leu c.6007T>C, p.Thr1303Met c.3908C>T, p.Val1950Met c.5848G>A

### SCN5A PMID 22882672

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Pro1177Leu carriers 1 vs 8 (error -7)

### SCN5A PMID 22966897

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P1177L, R190W
- Extra predictions: p.Ile1768Val, p.Pro2006Ala, p.Ser1103Tyr

### SCN5A PMID 23538678

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.S231CfsX251 c.692-693delCA carriers 2 vs 1 (error +1)

### SCN5A PMID 24573164

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg223Trp, p.Leu567Gln, p.Ala1112Val c.3335C>T, p.Asp1274Asn c.3820G>A, p.Glu1937Lys c.5809G>A, p.Gly1318Val c.3953G>T, p.Leu1500Val c.4498C>G, p.Ser1139Asn c.3416G>A

### SCN5A PMID 24581105

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 25171853

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Thr220Ile c.659C>T carriers 2 vs 1 (error +1); p.Thr220Ile c.659C>T unaffected 1 vs 0 (error +1)

### SCN5A PMID 25236808

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.F1617DEL
- Extra predictions: p.Phe1616del c.4844TCT[1]

### SCN5A PMID 26820365

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 28087622

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.K1505_Q1507DEL
- Extra predictions: R1918H, Y1795H, p.Arg1912His c.5735G>A, p.Gln1908Arg c.5723A>G, p.Glu1900Gln c.5698G>C

### SCN5A PMID 28339995

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 29544605

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Val1623Ile c.4867G>A

### SCN5A PMID 29709101

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Ile1660Val c.4978A>G, p.Arg1631His c.4892G>A, p.Ile1659Val c.4975A>G

### SCN5A PMID 32533946

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: D1430N, D356N, G1408R, G1712C, G1740R, G1743E, G1743R, G897E, I1660V, L846R, LEU839P, R104Q, R104W, R282H, R878C, R878H, R893H, S1218I, S910L, T187I, p.Ala1679Thr c.5035G>A, p.Ala1745Val c.5234C>T, p.Arg1582Cys c.4744C>T, p.Arg1631His c.4892G>A, p.Arg1631Ser c.4891C>A, p.Arg1897Cys c.5689C>T, p.Arg1957Ter c.5869C>T, p.Asn1721Asp c.5161A>G, p.Asp1242Asn c.3724G>A, p.Asp1369Gly c.4106A>G, p.Glu1224Lys c.3670G>A, p.Glu1573Lys c.4717G>A, p.Glu1783Lys c.5347G>A, p.Gly1261Ser c.3781G>A, p.Gly1419Ala c.4256G>C, p.Gly1419Cys c.4255G>T, p.Gly1641Glu c.4922G>A, p.Gly1739Arg c.5215G>A, p.Thr1708Lys c.5123C>A, p.Thr1708Met c.5123C>T, p.Tyr1448Cys c.4343A>G, p.Val1250Met c.3748G>A, p.Val1278Ile c.3832G>A, p.Val1352Met c.4054G>A

## Scope, method, and limitations

- Population: fixed manifest `tier2_gold_120.tsv` (118 papers); per-gene counts {'SCN5A': 30, 'KCNH2': 28, 'KCNQ1': 30, 'RYR2': 30}; every PMID has downloaded source and at least one gold assertion in each count field.
- Blinding: gold was used only for PMID eligibility and count-field presence during selection; extraction exported no gold values or row counts, and predictions were locked before `score` opened gold.
- Variant metrics are micro-averaged over gold rows. Precision treats unmatched predictions as false positives, although the curated recall packet may omit some real variants.
- Count MAE/RMSE are conditional on a supplied value. Count recall must be read alongside them because abstentions and missed variants are excluded from error magnitude.
- Source acquisition and gold completeness are separate from model reading quality; abstract-only or incomplete source is retained and labeled rather than silently excluded.
- The audited notation score is primary; the preserved raw score bounds sensitivity to post-lock matching adjudication.

## Reproducibility and evidence

- `selection.json`: selected PMIDs, source paths, source hashes, and available representations.
- `predictions.json`: immutable per-paper tools, rationales, extracted variants, counts, evidence quotes, source locations, and telemetry when captured.
- The production `gvf-run` trace manifest for every gene, including its exact call/decision records and write-time digest index, is SHA-256-bound in `predictions.json` and `LOCK.json` and revalidated before scoring.
- Each source run retains its own `llm_traces/<GENE>/<PMID>/` records and `llm_trace_report.html`; the evaluation projection does not copy or relabel those run-scoped records.
- `evidence.csv`: flat evidence ledger for every predicted variant.
- `paper_metrics.csv`: exact per-paper metrics.
- `LOCK.json`: SHA-256 digests proving prediction finalization before scoring.
- `report.json`: complete machine-readable score, errors, timing, and token usage.
- `matcher_adjudication.csv`: post-lock notation-equivalence audit; no extraction was edited.
- `report_raw_matcher.json` and `report_raw_matcher.md`: preserved pre-adjudication score.
- `validation_notes.md`: independent arithmetic, integrity checks, failure concentration, count outliers, and Claude comparison.
- `model_comparison.csv`: compact Codex/Claude comparison with scorer and telemetry caveats.
- `report_queries.sql`: executable DuckDB queries for the bounded analytical report datasets.

## Recommended next steps

1. Adjudicate extra predictions against the paper before treating precision as a production false-positive rate.
2. Review count outliers by source location and distinguish model mistakes from gold disagreements.
3. Add automatic fallback routing for data-rich papers that return zero or very few variants, then repeat with the same lock and count-recall definitions.

## Further questions

- Does table/PDF/OCR routing improve recall enough to justify its additional routing-call tokens?
- How much of the residual error is source incompleteness versus count-role interpretation?
