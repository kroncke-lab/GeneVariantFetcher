# Codex extraction-blinded paper evaluation — `20260820_gold119_noinference`

## Technical summary

This hash-locked run evaluated **119 papers** (**per-gene counts {'SCN5A': 30, 'KCNH2': 29, 'KCNQ1': 30, 'RYR2': 30}**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **45.1%**, recall **86.6%**, F1 **59.3%** (548 TP, 668 FP, 85 FN).
- Precision versus counted extras **97.5%** (548 matched rows; 14 extra rows with patient counts). The stricter count-bearing-only diagnostic is **93.3%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **2,463,525 total tokens** (1,803,614 input; 659,911 output).
- Elapsed: **14723.0s wall clock**; 8739.1s summed per-paper route + read time.
- Notation twins merged before scoring: **0** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 119}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `tier2_gold_120.tsv` (119 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: prediction content was finalized and SHA-256 locked before this external score. The source production workflow may have read registered gold for read-only layer scorecards before the projection lock; those scores did not feed back into extraction. This is not the stricter native lock-before-any-gold-read protocol.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 193 / 633 | 30.5% | 0.425 | 1.734 |
| affected | 48 / 633 | 7.6% | 0.688 | 1.898 |
| unaffected | 32 / 632 | 5.1% | 1.000 | 2.926 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 189 / 473 | 40.0% | 4 / 160 | 2.5% |
| affected | 41 / 407 | 10.1% | 7 / 226 | 3.1% |
| unaffected | 19 / 82 | 23.2% | 13 / 550 | 2.4% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 168 | 162 | 31 | 50.9% | 84.4% | 63.5% | 98.2% | 93.9% | 22.6% / 1.156 / 3.419 | 7.5% / 0.667 / 1.862 | 5.0% / 2.000 / 4.517 |
| KCNH2 | 66 | 90 | 23 | 42.3% | 74.2% | 53.9% | 90.4% | 83.3% | 38.2% / 0.324 / 0.748 | 18.0% / 0.688 / 1.820 | 7.9% / 0.286 / 0.756 |
| KCNQ1 | 148 | 131 | 19 | 53.0% | 88.6% | 66.4% | 97.4% | 94.3% | 39.5% / 0.167 / 0.564 | 4.8% / 1.375 / 2.894 | 3.0% / 1.800 / 3.606 |
| RYR2 | 166 | 285 | 12 | 36.8% | 93.3% | 52.8% | 100.0% | 100.0% | 27.0% / 0.167 / 0.540 | 5.1% / 0.111 / 0.333 | 5.6% / 0.100 / 0.316 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| KCNH2 | 10086971 | text | 3 | 4 | 1 | 42.9% | 75.0% | 54.5% | 75.0% / 0.333 | 0.0% / n/a | 25.0% / 0.000 | 115.4 | 27,588 |
| KCNH2 | 14642689 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 15500450 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 29.1 | 11,535 |
| KCNH2 | 16029385 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 182.9 | 24,871 |
| KCNH2 | 16155735 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 27.9 | 13,628 |
| KCNH2 | 18675227 | text | 3 | 6 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.333 | 66.7% / 1.000 | 0.0% / n/a | 68.3 | 25,087 |
| KCNH2 | 18808722 | text | 8 | 3 | 0 | 72.7% | 100.0% | 84.2% | 37.5% / 1.000 | 0.0% / n/a | 0.0% / n/a | 191.0 | 43,294 |
| KCNH2 | 19034806 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 30.7 | 12,292 |
| KCNH2 | 19065538 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 41.4 | 8,126 |
| KCNH2 | 19184172 | text | 4 | 7 | 0 | 36.4% | 100.0% | 53.3% | 25.0% / 0.000 | 25.0% / 0.000 | 25.0% / 0.000 | 436.0 | 40,198 |
| KCNH2 | 20181576 | text | 2 | 3 | 0 | 40.0% | 100.0% | 57.1% | 50.0% / 0.000 | 0.0% / n/a | 50.0% / 2.000 | 70.0 | 24,081 |
| KCNH2 | 21308345 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.1 | 14,253 |
| KCNH2 | 21779290 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 87.3 | 20,536 |
| KCNH2 | 22052944 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 1.000 | 100.0% / 1.000 | 100.0% / 0.000 | 32.4 | 19,695 |
| KCNH2 | 22067087 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 33.2 | 8,567 |
| KCNH2 | 22314138 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 36.0 | 14,760 |
| KCNH2 | 22338672 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 50.0% / 0.000 | 100.0% / 3.500 | 50.0% / 0.000 | 113.7 | 30,286 |
| KCNH2 | 22764740 | text | 1 | 13 | 0 | 7.1% | 100.0% | 13.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 51.5 | 21,550 |
| KCNH2 | 23917959 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 100.0% / 0.000 | 33.3 | 8,327 |
| KCNH2 | 23936059 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 65.1 | 24,581 |
| KCNH2 | 25819988 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 41.9 | 18,239 |
| KCNH2 | 26118593 | text | 1 | 5 | 0 | 16.7% | 100.0% | 28.6% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 40.5 | 16,249 |
| KCNH2 | 26746457 | text | 10 | 32 | 0 | 23.8% | 100.0% | 38.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 2.3 | 1,144 |
| KCNH2 | 29016797 | text | 2 | 4 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 2.000 | 50.0% / 1.000 | 0.0% / n/a | 64.1 | 25,022 |
| KCNH2 | 29121719 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 31.2 | 11,964 |
| KCNH2 | 29214556 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.5 | 6,916 |
| KCNH2 | 29650123 | text | 2 | 0 | 20 | 100.0% | 9.1% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 17.1 | 18,194 |
| KCNH2 | 30036649 | text | 5 | 3 | 0 | 62.5% | 100.0% | 76.9% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 61.4 | 17,892 |
| KCNH2 | 30246897 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 33.3% / 0.000 | 0.0% / n/a | 83.5 | 25,407 |
| KCNQ1 | 11351021 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 66.5 | 13,578 |
| KCNQ1 | 14678125 | text | 35 | 9 | 6 | 79.5% | 85.4% | 82.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 79.2 | 7,265 |
| KCNQ1 | 16155735 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 45.6 | 16,265 |
| KCNQ1 | 18567635 | text | 1 | 10 | 0 | 9.1% | 100.0% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 159.1 | 40,479 |
| KCNQ1 | 18580685 | text | 1 | 7 | 0 | 12.5% | 100.0% | 22.2% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 29.6 | 14,309 |
| KCNQ1 | 18808722 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 0.000 | 100.0% / 8.000 | 100.0% / 8.000 | 62.8 | 24,100 |
| KCNQ1 | 19056345 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.8 | 13,116 |
| KCNQ1 | 19114714 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 38.2 | 15,669 |
| KCNQ1 | 19632626 | text | 1 | 67 | 0 | 1.5% | 100.0% | 2.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 139.9 | 23,724 |
| KCNQ1 | 19687231 | text | 2 | 4 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 76.0 | 18,466 |
| KCNQ1 | 20348026 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 63.3 | 20,777 |
| KCNQ1 | 20368164 | text | 2 | 3 | 0 | 40.0% | 100.0% | 57.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 55.0 | 19,997 |
| KCNQ1 | 20959120 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.6 | 14,436 |
| KCNQ1 | 21070882 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 52.1 | 15,914 |
| KCNQ1 | 21129503 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 16.9 | 13,134 |
| KCNQ1 | 21956039 | text | 16 | 0 | 0 | 100.0% | 100.0% | 100.0% | 12.5% / 0.000 | 0.0% / n/a | 0.0% / n/a | 290.7 | 80,877 |
| KCNQ1 | 22613981 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.5 | 5,601 |
| KCNQ1 | 24667783 | text | 14 | 1 | 2 | 93.3% | 87.5% | 90.3% | 62.5% / 0.200 | 18.8% / 0.000 | 18.8% / 0.333 | 236.4 | 65,139 |
| KCNQ1 | 25471708 | text | 8 | 1 | 0 | 88.9% | 100.0% | 94.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 74.4 | 17,865 |
| KCNQ1 | 26496715 | text | 40 | 5 | 2 | 88.9% | 95.2% | 92.0% | 95.2% / 0.000 | 0.0% / n/a | 0.0% / n/a | 8.8 | 2,146 |
| KCNQ1 | 27114410 | text | 1 | 2 | 1 | 33.3% | 50.0% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 9.5 | 4,872 |
| KCNQ1 | 28491547 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 3.000 | 100.0% / 0.000 | 0.0% / n/a | 50.3 | 16,798 |
| KCNQ1 | 28739325 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 29.6 | 22,107 |
| KCNQ1 | 29677589 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 63.3 | 16,632 |
| KCNQ1 | 29851656 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 39.1 | 10,373 |
| KCNQ1 | 29876285 | text | 1 | 7 | 0 | 12.5% | 100.0% | 22.2% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 79.0 | 23,423 |
| KCNQ1 | 31293497 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 86.7 | 31,841 |
| KCNQ1 | 31520628 | text | 2 | 2 | 7 | 50.0% | 22.2% | 30.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 16.8 | 12,580 |
| KCNQ1 | 32168391 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 54.3 | 18,675 |
| KCNQ1 | 33082985 | text | 2 | 4 | 1 | 33.3% | 66.7% | 44.4% | 66.7% / 1.000 | 33.3% / 1.000 | 33.3% / 0.000 | 98.7 | 32,070 |
| RYR2 | 15466642 | text | 5 | 4 | 3 | 55.6% | 62.5% | 58.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 86.8 | 8,306 |
| RYR2 | 16517285 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 26.2 | 7,987 |
| RYR2 | 17161793 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 37.7 | 10,354 |
| RYR2 | 18285261 | text | 5 | 1 | 0 | 83.3% | 100.0% | 90.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 61.4 | 14,078 |
| RYR2 | 19216760 | text | 4 | 1 | 1 | 80.0% | 80.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 49.7 | 29,346 |
| RYR2 | 19398417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 0.000 | 42.2 | 8,716 |
| RYR2 | 19398665 | text | 27 | 1 | 0 | 96.4% | 100.0% | 98.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.9 | 17,524 |
| RYR2 | 19926015 | text | 38 | 71 | 2 | 34.9% | 95.0% | 51.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 152.2 | 41,681 |
| RYR2 | 20395638 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 53.2 | 8,938 |
| RYR2 | 21652165 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 39.1 | 10,133 |
| RYR2 | 22178870 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 51.5 | 18,343 |
| RYR2 | 22222782 | text | 7 | 0 | 0 | 100.0% | 100.0% | 100.0% | 85.7% / 0.833 | 0.0% / n/a | 0.0% / n/a | 247.1 | 51,892 |
| RYR2 | 22677073 | text | 24 | 1 | 1 | 96.0% | 96.0% | 96.0% | 80.0% / 0.000 | 4.0% / 0.000 | 12.0% / 0.000 | 363.1 | 94,263 |
| RYR2 | 25435091 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 38.6 | 11,438 |
| RYR2 | 25463374 | text | 1 | 10 | 0 | 9.1% | 100.0% | 16.7% | 100.0% / 1.000 | 0.0% / n/a | 100.0% / 1.000 | 54.9 | 25,418 |
| RYR2 | 25500949 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 36.2 | 13,405 |
| RYR2 | 25814417 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 31.1 | 20,090 |
| RYR2 | 27114410 | text | 4 | 3 | 5 | 57.1% | 44.4% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.1 | 4,861 |
| RYR2 | 27839804 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 50.9 | 12,523 |
| RYR2 | 28798025 | text | 6 | 76 | 0 | 7.3% | 100.0% | 13.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 67.4 | 37,562 |
| RYR2 | 29350269 | text | 3 | 41 | 0 | 6.8% | 100.0% | 12.8% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 49.0 | 29,114 |
| RYR2 | 30403697 | text | 21 | 14 | 0 | 60.0% | 100.0% | 75.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 128.4 | 35,202 |
| RYR2 | 30546600 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 0.000 | 106.4 | 19,934 |
| RYR2 | 30835254 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 2.9 | 973 |
| RYR2 | 32866913 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 100.0% / 0.000 | 48.1 | 8,498 |
| RYR2 | 33315912 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 34.2 | 20,641 |
| RYR2 | 33536282 | text | 1 | 43 | 0 | 2.3% | 100.0% | 4.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 106.0 | 32,869 |
| RYR2 | 33640691 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 312.5 | 11,382 |
| RYR2 | 34661651 | text | 4 | 8 | 0 | 33.3% | 100.0% | 50.0% | 75.0% / 0.333 | 0.0% / n/a | 0.0% / n/a | 85.0 | 30,378 |
| RYR2 | 35663620 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 29.3 | 12,168 |
| SCN5A | 15671429 | text | 4 | 15 | 0 | 21.1% | 100.0% | 34.8% | 100.0% / 5.750 | 100.0% / 1.500 | 100.0% / 4.250 | 219.9 | 28,199 |
| SCN5A | 15828879 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 46.4 | 9,335 |
| SCN5A | 15898185 | text | 2 | 0 | 9 | 100.0% | 18.2% | 30.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.3 | 5,211 |
| SCN5A | 16301357 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 43.6 | 9,188 |
| SCN5A | 16929919 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 36.3 | 7,822 |
| SCN5A | 17675083 | text | 2 | 4 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 2.000 | 100.0% / 2.000 | 0.0% / n/a | 95.6 | 32,301 |
| SCN5A | 17971661 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.8 | 4,248 |
| SCN5A | 19305639 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.1 | 10,392 |
| SCN5A | 19406494 | text | 1 | 8 | 0 | 11.1% | 100.0% | 20.0% | 100.0% / 1.000 | 0.0% / n/a | 100.0% / 1.000 | 65.1 | 17,307 |
| SCN5A | 19549036 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 31.3 | 13,064 |
| SCN5A | 20038812 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 18.5 | 19,017 |
| SCN5A | 20539757 | text | 11 | 14 | 0 | 44.0% | 100.0% | 61.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 102.5 | 25,460 |
| SCN5A | 21193062 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 36.9 | 18,602 |
| SCN5A | 21288276 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.7 | 10,337 |
| SCN5A | 22101522 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 35.3 | 18,342 |
| SCN5A | 22677073 | text | 12 | 5 | 0 | 70.6% | 100.0% | 82.8% | 100.0% / 0.083 | 16.7% / 0.000 | 8.3% / 0.000 | 242.7 | 70,882 |
| SCN5A | 22685113 | text | 0 | 0 | 10 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 22882672 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 7.000 | 0.0% / n/a | 0.0% / n/a | 42.3 | 8,896 |
| SCN5A | 22966897 | text | 0 | 3 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.7 | 12,168 |
| SCN5A | 23538678 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 28.6 | 10,700 |
| SCN5A | 24573164 | text | 20 | 9 | 0 | 69.0% | 100.0% | 81.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 104.8 | 22,350 |
| SCN5A | 24581105 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 44.9 | 11,503 |
| SCN5A | 25171853 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 1.000 | 50.0% / 0.000 | 0.0% / n/a | 87.2 | 14,888 |
| SCN5A | 25236808 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 2.000 | 100.0% / 0.000 | 100.0% / 2.000 | 40.1 | 8,702 |
| SCN5A | 26820365 | text | 1 | 10 | 0 | 9.1% | 100.0% | 16.7% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 63.1 | 32,819 |
| SCN5A | 28087622 | text | 6 | 8 | 1 | 42.9% | 85.7% | 57.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 82.1 | 27,510 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 47.4 | 16,063 |
| SCN5A | 29544605 | text | 6 | 27 | 0 | 18.2% | 100.0% | 30.8% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 88.5 | 41,584 |
| SCN5A | 29709101 | text | 11 | 5 | 0 | 68.8% | 100.0% | 81.5% | 90.9% / 1.100 | 27.3% / 0.000 | 27.3% / 0.000 | 232.6 | 66,344 |
| SCN5A | 32533946 | text | 77 | 45 | 6 | 63.1% | 92.8% | 75.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 242.5 | 75,764 |

## Errors and representation choices

### KCNH2 PMID 10086971

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: C1117fsX
- Extra predictions: c.3108insG, p.Asp1037fs c.3107dup, p.Gly965fs c.2892dup, p.Leu953fs c.2857dup
- Count disagreements: c.2592+1G>A carriers 6 vs 5 (error +1)

### KCNH2 PMID 14642689

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A561T

### KCNH2 PMID 15500450

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 16029385

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G601S, T65P

### KCNH2 PMID 16155735

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P923fsX

### KCNH2 PMID 18675227

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1116V, K897T, L955C, Q725X, R1014X, T1062I
- Count disagreements: p.Leu955Val c.2863C>G carriers 1 vs 2 (error -1); p.Arg954Cys c.2860C>T affected 1 vs 0 (error +1); p.Leu955Val c.2863C>G affected 1 vs 0 (error +1)

### KCNH2 PMID 18808722

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1116V, A490P, G189R
- Count disagreements: A490T carriers 5 vs 7 (error -2); K897T carriers 6 vs 7 (error -1)

### KCNH2 PMID 19034806

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 19065538

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 19184172

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: K897T, p.Ala536Val c.1607C>T, p.Arg1068Ter c.3202C>T, p.Arg176Trp c.526C>T, p.Glu538Lys c.1612G>A, p.Glu786Lys c.2356G>A, p.Glu788Asp c.2364G>C

### KCNH2 PMID 20181576

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: K897K, P926A, p.Ala926ProfsTer14
- Count disagreements: K897T unaffected 1 vs 3 (error -2)

### KCNH2 PMID 21308345

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.1600G>T

### KCNH2 PMID 21779290

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Leu564Leu, p.Tyr652Tyr

### KCNH2 PMID 22052944

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A614V
- Count disagreements: p.Arg176Trp carriers 1 vs 2 (error -1); p.Arg176Trp affected 0 vs 1 (error -1)

### KCNH2 PMID 22067087

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 22314138

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: R744X, R774P

### KCNH2 PMID 22338672

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1116V
- Count disagreements: K897T affected 7 vs 0 (error +7)

### KCNH2 PMID 22764740

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A490T, A915V, G572R, G601S, I593R, K897T, N470D, P596R, R1047L, T474I, T613M, V612L, Y611H

### KCNH2 PMID 23917959

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Gly603Asp c.1808G>A carriers 1 vs 2 (error -1)

### KCNH2 PMID 23936059

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 25819988

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: H562P

### KCNH2 PMID 26118593

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: F29L, c.2229_2230insC, p.F513F c.1539C>T, p.I489I c.1467C>T, p.L564L c.1692A>G

### KCNH2 PMID 26746457

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A175D, A188V, A190T, A9V, D982N, G238S, G873S, G903R, K897T, L1105S, M291R, P1018A, P1030L, P347S, P967L, R1035Q, R1047L, R148W, R181Q, R252W, R326H, R328C, R394H, R885C, R894H, R912W, R948S, S981G, V1063L, V325M, W927L, Y403C

### KCNH2 PMID 29016797

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.2264G>C, c.5010G>T, c.5605C>T, c.631C>T
- Count disagreements: p.Asn588Lys c.1764C>G carriers 1 vs 0 (error +1); p.Ser631Ala c.1891T>G carriers 3 vs 0 (error +3); p.Asn588Lys c.1764C>G affected 1 vs 0 (error +1)

### KCNH2 PMID 29121719

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 29214556

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 29650123

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A561V, C49Y, F617fsX, G572S, G911fsX, L109fsX, L779P, N633S, Q1046X, R1035fsX, R148W, R176W, R328C, R534C, R744X, R892fsX, S660L, S818P, W412X, Y43C

### KCNH2 PMID 30036649

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G314F, G314S, G628D

### KCNH2 PMID 30246897

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 11351021

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 14678125

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: L151, G185S, DELF340, A344/SP, A344A/SPLICE, R518X
- Extra predictions: p.Arg243His c.728G>A, p.Arg243Ser c.727C>A, p.Arg366Leu c.1097G>T, p.Arg366Pro c.1097G>C, p.Asp317Ala c.950A>C, p.Glu284del c.850_852del, p.Gly186Cys c.556G>T, p.Gly186Ser c.556G>A, p.Leu266Arg c.797T>G

### KCNQ1 PMID 16155735

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Pro448Arg c.1343C>G

### KCNQ1 PMID 18567635

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: D242A, E261A, K218A, R228A, R231A, R237A, R243A, R249A, R259A, p.Pro475Asp

### KCNQ1 PMID 18580685

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: D202H, I204M, Q357R, R294H, R555C, S209F, V215M

### KCNQ1 PMID 18808722

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G189R, R190Q, c.2020insAG
- Count disagreements: L187P affected 5 vs 13 (error -8); L187P unaffected 26 vs 18 (error +8)

### KCNQ1 PMID 19056345

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 19114714

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 19632626

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: S140G, S140P, S209F, V141M, c.1033-2del, c.1251+9C>A, c.1394-1G>T, c.1515-4G>A, c.477+5G>A, p.Ala178Thr c.532G>A, p.Ala341Glu c.1022C>A, p.Ala341Val c.1022C>T, p.Ala344= c.1032G>A, p.Ala344Val c.1031C>T, p.Arg190Gln c.569G>A, p.Arg190Leu c.569G>T, p.Arg190Trp c.568C>T, p.Arg243Cys c.727C>T, p.Arg243His c.728G>A, p.Arg380Lys c.1139G>A, p.Arg401fs c.1201dup, p.Arg507Gln c.1520G>A, p.Arg518Ter c.1552C>T, p.Arg539Trp c.1615C>T, p.Arg555Ser c.1663C>A, p.Arg583His c.1748G>A, p.Arg591His c.1772G>A, p.Arg594Gln c.1781G>A, p.Arg632fs c.1893dup, p.Asp202Asn c.604G>A, p.Cys214Tyr c.641G>A, p.Cys642Phe c.1925G>T, p.Gln107His c.321G>C, p.Gln359Ter c.1075C>T, p.Gln376Ter c.1126C>T, p.Gln530Ter c.1588C>T, p.Glu449fs c.1343dup, p.Gly168Arg c.502G>A, p.Gly269Ser c.805G>A, p.Gly306Arg c.916G>A, p.Gly314Ser c.940G>A, p.Gly325Arg c.973G>A, p.Gly626fs c.1875dup, p.His258Arg c.773A>G, p.Ile375fs c.1124_1127del, p.Leu266Pro c.797T>C, p.Leu266fs c.796del, p.Leu496fs c.1486_1487del, p.Lys362Arg c.1085A>G, p.Lys526Glu c.1576A>G, p.Met520Arg c.1559T>G, p.Pro343Ala c.1027C>G, p.Ser277del c.825CTC[1], p.Ser546Leu c.1637C>T, p.Ser566Phe c.1697C>T, p.Thr311Ile c.932C>T, p.Thr312Ile c.935C>T, p.Thr322Met c.965C>T, p.Trp15Ter c.45G>A, p.Trp305Leu c.914G>T, p.Trp305Ter c.914G>A, p.Tyr111Cys c.332A>G, p.Tyr184Ser c.551A>C, p.Val110Ile c.328G>A, p.Val215Met c.643G>A, p.Val254Met c.760G>A, p.Val280Gly c.839T>G

### KCNQ1 PMID 19687231

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Ile328Leu, p.Lys326Arg, p.Thr327Val, p.Val324Leu

### KCNQ1 PMID 20348026

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G628S

### KCNQ1 PMID 20368164

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A340I, p.Ala341Glu, p.Thr312Ile

### KCNQ1 PMID 20959120

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 21070882

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Thr96Arg carriers 1 vs 3 (error -2)

### KCNQ1 PMID 21129503

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A341V

### KCNQ1 PMID 21956039

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 22613981

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 24667783

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: I328_S330DEL, Y522X
- Extra predictions: p.Thr224Met c.671C>T
- Count disagreements: p.Asp454Thrfs*7 carriers 2 vs 3 (error -1); p.Pro73Thr c.217C>A carriers 1 vs 2 (error -1); p.Asp454Thrfs*7 unaffected 1 vs 2 (error -1)

### KCNQ1 PMID 25471708

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.572_576del

### KCNQ1 PMID 26496715

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: 360_361DUPKQ, A384fsX
- Extra predictions: c.1251+2T>C, c.477+5G>A, c.605-2G>A, c.921+1G>T, p.Gln359_Arg c.1067AGCAGA[3]

### KCNQ1 PMID 27114410

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: G584S
- Extra predictions: N634FS, N98S

### KCNQ1 PMID 28491547

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Val141Met c.421G>A carriers 3 vs 6 (error -3)

### KCNQ1 PMID 28739325

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: D242A
- Count disagreements: p.D242N c.724G>A carriers 2 vs 0 (error +2)

### KCNQ1 PMID 29677589

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.R190W c.568C>T affected 1 vs 0 (error +1); p.R594Q c.1781G>A affected 1 vs 0 (error +1)

### KCNQ1 PMID 29851656

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 29876285

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A300T, A302V, T311I, V319L, W305L, W305S, W305X

### KCNQ1 PMID 31293497

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.R360_Q361insQKQR, p.Lys362_His c.1077_1088dup

### KCNQ1 PMID 31520628

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: G589D, R366W, Y184S, Y315N, V254M, I235N, R252fsX
- Extra predictions: A341V, p.Arg591His c.1772G>A

### KCNQ1 PMID 32168391

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 33082985

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: W248F + L347R
- Extra predictions: L347P, T372M, W248C, W248R
- Count disagreements: Leu347Arg carriers 2 vs 1 (error +1); Trp248Phe carriers 2 vs 1 (error +1); Trp248Phe affected 1 vs 0 (error +1)

### RYR2 PMID 15466642

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P164S, F4499C, G4671R
- Extra predictions: p.Ala2403Val c.7208C>T, p.Ala4510Val c.13529C>T, p.Arg414Cys c.1240C>T, p.Arg414His c.1241G>A

### RYR2 PMID 16517285

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 17161793

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 18285261

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: S453S

### RYR2 PMID 19216760

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: EXON 3 DELETION
- Extra predictions: c.(?_169

### RYR2 PMID 19398417

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 19398665

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G3946S

### RYR2 PMID 19926015

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: M81L, EXON 3 DELETION
- Extra predictions: A4282V, G1885E, G1886S, M2389L, Q2958R, R169Q, V2321M, V4010M, c.(?_169, c.169-198_273+820del, p.Ala2213Val c.6638C>T, p.Ala2387Thr c.7159G>A, p.Ala4510Ser c.13528G>T, p.Ala4510Thr c.13528G>A, p.Ala77Val c.230C>T, p.Arg1051His c.3152G>A, p.Arg176Gln c.527G>A, p.Arg176Leu c.527G>T, p.Arg2401Cys c.7201C>T, p.Arg2401His c.7202G>A, p.Arg3570Trp c.10708C>T, p.Arg414Cys c.1240C>T, p.Arg4307Cys c.12919C>T, p.Arg4959Gln c.14876G>A, p.Arg838Leu c.2513G>T, p.Asn3308Ser c.9923A>G, p.Asp2431Tyr c.7291G>T, p.Asp4001Asn c.12001G>A, p.Gln3811Pro c.11432A>C, p.Gln4936Arg c.14807A>G, p.Glu1724Lys c.5170G>A, p.Glu4431Lys c.13291G>A, p.Gly172Glu c.515G>A, p.Gly2949Val c.8846G>T, p.Gly3946Ser c.11836G>A, p.Gly4662Ser c.13984G>A, p.Gly4864Cys c.14590G>T, p.Gly4935Arg c.14803G>A, p.Gly809Glu c.2426G>A, p.His202Gln c.606C>A, p.His240Arg c.719A>G, p.Ile2009Met c.6027T>G, p.Ile419Ser c.1256T>G, p.Ile4857Val c.14569A>G, p.Leu1459Ser c.4376T>C, p.Leu2123Phe c.6367C>T, p.Leu4105Phe c.12313C>T, p.Leu4698Pro c.14093T>C, p.Lys3717Arg c.11150A>G, p.Met4109Thr c.12326T>C, p.Phe2307Leu c.6921C>G, p.Pro1395Ala c.4183C>G, p.Pro2328Leu c.6983C>T, p.Pro466Arg c.1397C>G, p.Pro990Gln c.2969C>A, p.Ser221Gly c.661A>G, p.Ser2246Leu c.6737C>T, p.Ser2312Gly c.6934A>G, p.Ser2393Pro c.7177T>C, p.Ser2653Tyr c.7958C>A, p.Ser4155Thr c.12463T>A, p.Ser4667Ile c.14000G>T, p.Thr1093Pro c.3277A>C, p.Thr2510Ala c.7528A>G, p.Thr4755Ile c.14264C>T, p.Tyr2392Cys c.7175A>G, p.Val2113Met c.6337G>A, p.Val2306Ile c.6916G>A, p.Val4771Ile c.14311G>A, p.Val4771Phe c.14311G>T, p.Val4778Leu c.14332G>T

### RYR2 PMID 20395638

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Gly4935Arg c.14803G>A

### RYR2 PMID 21652165

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: R4496C

### RYR2 PMID 22178870

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 22222782

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.D1872N carriers 1 vs 2 (error -1); p.H4579Y carriers 1 vs 4 (error -3); p.Q486H carriers 1 vs 2 (error -1)

### RYR2 PMID 22677073

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: Q2958R
- Extra predictions: c.2398+5G>T

### RYR2 PMID 25435091

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: R420Q

### RYR2 PMID 25463374

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A77V, G230C, P2328S, Q4201R, R2474S, R3570W, R4497C, S2226L, S4565R, V4653F
- Count disagreements: H29D carriers 3 vs 4 (error -1); H29D unaffected 0 vs 1 (error -1)

### RYR2 PMID 25500949

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1744S, A189T, D4301N

### RYR2 PMID 25814417

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: V2475F

### RYR2 PMID 27114410

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: T4158P, L3879P, Q3925E, S3959L, G4936R
- Extra predictions: N634FS, N98S, p.Gly4935Arg c.14803G>A

### RYR2 PMID 27839804

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: H29D

### RYR2 PMID 28798025

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A119S, A129V, A1462T, A177V, A2019V, A2025V, A222T, A281G, A449V, A4919T, A856T, C132X, C4829Y, C81Y, D145E, D62N, D660Y, D888N, E102K, E26Q, G542R, H117H, H802N, I4286L, I615N, K133I, K1617T, K319N, K3674N, L1253V, L1256L, L1320P, L3057V, N797K, P1750L, P225L, P2777H, P407R, P638L, Q1024X, Q232R, R1042C, R1314Q, R138W, R1468Q, R199C, R208H, R2189Q, R274K, R286W, R370L, R408H, R438FS, R552X, R688X, R955W, S111F, S1195Y, S455L, T1929M, T2063M, T250I, T2664N, T266M, T277M, T286A, T429M, T46I, V39M, V532M, V927G, c.34930+2A>G, c.37448_37467del, c.46521delT, c.81261delC, c.98972_98975del

### RYR2 PMID 29350269

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A3039P, D548E, E2013K, G550R, L22P, N2045S, N3982D, P736L, Q2432H, R1250H, R1278Q, R2639Q, R272C, R420C, R66W, S160T, T116M, T240M, Y62N, c.1190T>C, c.1234G>A, c.1297+5G>T, c.1555G>A, c.184A>T, c.19505T>C, c.196G>A, c.1978+3C>A, c.2207G>A, c.2246G>A, c.2928+5C>T, c.2953G>A, c.331C>A, c.347G>A, c.3749C>T, c.478A>T, c.554C>T, c.65A>G, c.707A>G, c.7296C>A, c.8605A>G, c.9115C>G

### RYR2 PMID 30403697

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A103V, C1274Y, D4642H, H1042R, I882V, IVS5+1G>C, K296Q, K4742R, L714R, M0781T, Q584R, R1458G, R4954K, S2774G

### RYR2 PMID 30546600

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 30835254

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Thr1107Met c.3320C>T

### RYR2 PMID 32866913

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Thr85Ile c.254C>T carriers 3 vs 2 (error +1); p.Thr85Ile c.254C>T affected 2 vs 1 (error +1)

### RYR2 PMID 33315912

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G357S

### RYR2 PMID 33536282

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: L3995V, Q3776L, c.11983A>C, p.Ala4860Gly, p.Arg4496Cys, p.Arg4594Lys, p.Asp4112Asn, p.Asp4122Asn, p.Asp4646Ala, p.Gln4879His, p.Gln4897His, p.Ile2074Thr, p.Ile2075Thr, p.Ile3995Val, p.Ile4855Met, p.Lys4594Arg, p.Met1045Val, p.Ser11Gly, p.Ser4938Phe, p.Ser88*, p.Thr4196Ile, c.848+1G>A, p.Ala2387Thr c.7159G>A, p.Ala2403Val c.7208C>T, p.Ala2633Thr c.7897G>A, p.Arg1013Gln c.3038G>A, p.Arg1089Cys c.3265C>T, p.Arg1671Gln c.5012G>A, p.Arg2144Cys c.6430C>T, p.Arg76Trp c.226C>T, p.Glu3716Lys c.11146G>A, p.Glu4076Lys c.12226G>A, p.Glu4182Gln c.12544G>C, p.Gly2284Asp c.6851G>A, p.His4579Tyr c.13735C>T, p.His4908Asn c.14722C>A, p.Lys304Glu c.910A>G, p.Pro2328Ser c.6982C>T, p.Pro916Leu c.2747C>T, p.Ser4153Arg c.12457A>C, p.Thr2188Ala c.6562A>G, p.Tyr2392Cys c.7175A>G, p.Tyr4725Cys c.14174A>G

### RYR2 PMID 33640691

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 34661651

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A4860G, I4855M, K4750Q, R169Q, R2267H, R4496C, S4565R, c.14151+1G>A
- Count disagreements: p.Ser4168Pro c.12502T>C carriers 1 vs 2 (error -1)

### RYR2 PMID 35663620

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1357V

### SCN5A PMID 15671429

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: C659T, D1595N, R225W, S1103Y, V1397X, W156X, c.2551insTG, c.2550_2551insTG, c.704-2A>G, p.Asp1274Asn c.3820G>A, p.Asp1594Asn c.4780G>A, p.Asp1594His c.4780G>C, p.Gly1572= c.4716C>T, p.Phe851fs c.2550_2551dup, p.Pro1331fs c.3992del
- Count disagreements: D1275N carriers 2 vs 22 (error -20); D1595H carriers 2 vs 1 (error +1); R814W carriers 2 vs 1 (error +1); T220I carriers 2 vs 1 (error +1); D1275N affected 1 vs 7 (error -6); D1275N unaffected 1 vs 15 (error -14); D1595H unaffected 1 vs 0 (error +1); R814W unaffected 1 vs 0 (error +1); T220I unaffected 1 vs 0 (error +1)

### SCN5A PMID 15828879

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Asn406Ser c.1217A>G

### SCN5A PMID 15898185

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: F1293S, H558R, P1090L, R1512W, R34C, R481W, S1787N, S524Y, V1951L

### SCN5A PMID 16301357

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Leu1825Pro carriers 1 vs 0 (error +1)

### SCN5A PMID 16929919

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17675083

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: M1766L, R282H, T512I, p.Arg34Cys
- Count disagreements: H558R carriers 5 vs 1 (error +4); H558R affected 5 vs 1 (error +4)

### SCN5A PMID 17971661

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: E1784K

### SCN5A PMID 19305639

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 19406494

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Ala1949Pro c.5845G>C, p.Trp301* c.903G>A, p.Arg2011Cys c.6031C>T, p.Glu1866Ter c.5596G>T, p.Leu1785fs c.5353_5354del, p.Phe2003Leu c.6007T>C, p.Pro1840fs c.5517dup, p.Trp1712Ter c.5135G>A
- Count disagreements: p.Arg808Cys c.2422C>T carriers 4 vs 3 (error +1); p.Arg808Cys c.2422C>T unaffected 1 vs 0 (error +1)

### SCN5A PMID 19549036

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Arg1193Gln, p.Ser1103Tyr

### SCN5A PMID 20038812

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: I1461T, R1448P, p.Phe1485Leu c.4453T>C

### SCN5A PMID 20539757

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: F1617del, H1298L, H558R, K1578FS, R1623H, R1632X, c.475_482+3del, p.Arg1622Ter c.4864C>T, p.Arg1631His c.4892G>A, p.Arg1631Ser c.4891C>A, p.Asp1274Asn c.3820G>A, p.Gly1407Arg c.4219G>A, p.Phe1616del c.4844TCT[1], p.Pro1297Leu c.3890C>T

### SCN5A PMID 21193062

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 21288276

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.999-424_1338+81del
- Extra predictions: p.Trp193fs c.576del

### SCN5A PMID 22101522

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: S1503A

### SCN5A PMID 22677073

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G38S, S1787N, p.Ile1767Val c.5299A>G, p.Pro2005Ala c.6013C>G, p.Val1950Leu c.5848G>T
- Count disagreements: S216L carriers 2 vs 1 (error +1)

### SCN5A PMID 22685113

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R1626H, D1819N, R340Q, R1897W, V1951M, F1596I, F2004L, S216L, T1304M, T220I

### SCN5A PMID 22882672

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Pro1177Leu carriers 1 vs 8 (error -7)

### SCN5A PMID 22966897

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P1177L, R190W
- Extra predictions: p.Ile1768Val, p.Pro2006Ala, p.Ser1103Tyr

### SCN5A PMID 23538678

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.S231CfsX251
- Count disagreements: c.693delCA carriers 2 vs 1 (error +1)

### SCN5A PMID 24573164

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: C373Y, R223W, p.Leu567Gln, p.Ala1112Val c.3335C>T, p.Asp1274Asn c.3820G>A, p.Glu1937Lys c.5809G>A, p.Gly1318Val c.3953G>T, p.Leu1500Val c.4498C>G, p.Ser1139Asn c.3416G>A

### SCN5A PMID 24581105

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 25171853

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.T220I c.659C>T carriers 2 vs 1 (error +1)

### SCN5A PMID 25236808

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Phe1616del c.4844TCT[1]
- Count disagreements: F1617DEL carriers 2 vs 4 (error -2); F1617DEL unaffected 1 vs 3 (error -2)

### SCN5A PMID 26820365

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: D561A, G582S, G737R, I376K, I376T, P1204L, Q854R, R442C c.1324C>T, T1127C, T286T

### SCN5A PMID 28087622

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P.K1505_Q1507DEL
- Extra predictions: C373Y, H1853R, R1918H, Y1795H, Y371C, p.Arg1912His c.5735G>A, p.Gln1908Arg c.5723A>G, p.Glu1900Gln c.5698G>C

### SCN5A PMID 28339995

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 29544605

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: K416FS*0, P.F33LFS, P.K60RFS, R547X, T372M, c.1244_1245insT, c.124G>A, c.13774G>T, c.1378G>A, c.1560+1G>A, c.179_180delAA, c.2148+1G>A, c.26011C>T, c.325-2A>G, c.36_38delAAG, c.37432_37433insGTGGTTACTACAGCCTC, c.37437_37438insG, c.415-1G>T, c.47653delA, c.523+2T>C, c.546delT, c.650-2A>G, c.6540delG, c.843+1G>T, c.861C>G, c.98_107delTTCACTTCAC, p.Val1623Ile c.4867G>A

### SCN5A PMID 29709101

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: F756L, F756fsX, I1660V, p.Arg1631His c.4892G>A, p.Ile1659Val c.4975A>G
- Count disagreements: H886Q carriers 4 vs 7 (error -3); R376H carriers 6 vs 12 (error -6); W822* carriers 4 vs 6 (error -2)

### SCN5A PMID 32533946

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: D1243N, L1346P, S1382I, E1574K, R1583C, N1722D
- Extra predictions: A1632H, D1430N, D356N, G1406R, G1408R, G1740R, G752A, G897E, I1660V, L846R, M1766L, P.LEU839P, R104Q, R104W, R282H, R878C, R878H, R893H, S1218I, S910L, T187I, p.Ala1679Thr c.5035G>A, p.Ala1745Val c.5234C>T, p.Arg1582Cys c.4744C>T, p.Arg1631His c.4892G>A, p.Arg1631Ser c.4891C>A, p.Arg1897Cys c.5689C>T, p.Arg1957Ter c.5869C>T, p.Asn1721Asp c.5161A>G, p.Asp1242Asn c.3724G>A, p.Asp1369Gly c.4106A>G, p.Glu1224Lys c.3670G>A, p.Glu1573Lys c.4717G>A, p.Glu1783Lys c.5347G>A, p.Gly1261Ser c.3781G>A, p.Gly1419Ala c.4256G>C, p.Gly1419Cys c.4255G>T, p.Gly1641Glu c.4922G>A, p.Gly1739Arg c.5215G>A, p.Thr1708Lys c.5123C>A, p.Thr1708Met c.5123C>T, p.Tyr1448Cys c.4343A>G, p.Val1250Met c.3748G>A, p.Val1278Ile c.3832G>A, p.Val1352Met c.4054G>A

## Scope, method, and limitations

- Population: fixed manifest `tier2_gold_120.tsv` (119 papers); per-gene counts {'SCN5A': 30, 'KCNH2': 29, 'KCNQ1': 30, 'RYR2': 30}; every PMID has downloaded source and at least one gold assertion in each count field.
- Blinding: prediction content was finalized and SHA-256 locked before this external score. The source production workflow may have read registered gold for read-only layer scorecards before the projection lock; those scores did not feed back into extraction. This is not the stricter native lock-before-any-gold-read protocol.
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
