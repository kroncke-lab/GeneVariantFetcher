# Codex extraction-blinded paper evaluation — `20260902_false_zero_recovery_gold118`

## Technical summary

This hash-locked run evaluated **118 papers** (**per-gene counts {'SCN5A': 30, 'KCNH2': 28, 'KCNQ1': 30, 'RYR2': 30}**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **66.4%**, recall **85.9%**, F1 **74.9%** (543 TP, 275 FP, 89 FN).
- Precision versus counted extras **98.2%** (543 matched rows; 10 extra rows with patient counts). The stricter count-bearing-only diagnostic is **95.4%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **2,746,411 total tokens** (1,981,017 input; 765,394 output).
- Elapsed: **10822.0s wall clock**; 9928.8s summed per-paper route + read time.
- Notation twins merged before scoring: **1** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 118}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `tier2_gold_120.tsv` (118 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility and count-field presence during selection; extraction exported no gold values or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 207 / 632 | 32.8% | 0.179 | 0.692 |
| affected | 120 / 632 | 19.0% | 0.508 | 2.365 |
| unaffected | 37 / 631 | 5.9% | 2.081 | 7.427 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 204 / 472 | 43.2% | 3 / 160 | 1.9% |
| affected | 105 / 406 | 25.9% | 15 / 226 | 6.6% |
| unaffected | 25 / 82 | 30.5% | 12 / 549 | 2.2% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 184 | 110 | 15 | 62.6% | 92.5% | 74.6% | 98.4% | 93.3% | 21.1% / 0.286 / 1.155 | 17.1% / 0.118 / 0.485 | 7.5% / 0.400 / 0.816 |
| KCNH2 | 60 | 41 | 28 | 59.4% | 68.2% | 63.5% | 100.0% | 100.0% | 31.8% / 0.179 / 0.627 | 21.6% / 0.474 / 1.051 | 9.1% / 1.250 / 1.936 |
| KCNQ1 | 145 | 100 | 22 | 59.2% | 86.8% | 70.4% | 96.0% | 93.8% | 54.5% / 0.110 / 0.445 | 18.6% / 0.516 / 1.244 | 3.0% / 2.000 / 2.757 |
| RYR2 | 154 | 24 | 24 | 86.5% | 86.5% | 86.5% | 99.4% | 97.9% | 25.8% / 0.217 / 0.552 | 20.2% / 0.889 / 4.062 | 5.1% / 5.667 / 14.769 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| KCNH2 | 10086971 | text | 3 | 1 | 1 | 75.0% | 75.0% | 75.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 162.1 | 31,797 |
| KCNH2 | 15500450 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 39.3 | 12,943 |
| KCNH2 | 16029385 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 103.6 | 22,320 |
| KCNH2 | 16155735 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 66.2 | 12,269 |
| KCNH2 | 18675227 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 67.4 | 19,148 |
| KCNH2 | 18808722 | text | 8 | 1 | 0 | 88.9% | 100.0% | 94.1% | 62.5% / 0.000 | 25.0% / 0.000 | 0.0% / n/a | 236.1 | 57,432 |
| KCNH2 | 19034806 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 40.3 | 10,720 |
| KCNH2 | 19065538 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 43.9 | 10,130 |
| KCNH2 | 19184172 | text | 3 | 0 | 1 | 100.0% | 75.0% | 85.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 224.6 | 21,821 |
| KCNH2 | 20181576 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 0.000 | 100.0% / 2.000 | 100.0% / 3.000 | 105.6 | 26,092 |
| KCNH2 | 21308345 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 73.1 | 16,992 |
| KCNH2 | 21779290 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 114.8 | 24,796 |
| KCNH2 | 22052944 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 49.6 | 17,532 |
| KCNH2 | 22067087 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 35.7 | 9,942 |
| KCNH2 | 22314138 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 53.0 | 17,169 |
| KCNH2 | 22338672 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 50.0% / 0.000 | 50.0% / 0.000 | 50.0% / 0.000 | 94.5 | 30,600 |
| KCNH2 | 22764740 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 42.4 | 19,916 |
| KCNH2 | 23917959 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 0.000 | 0.0% / n/a | 62.4 | 11,134 |
| KCNH2 | 23936059 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 116.9 | 28,973 |
| KCNH2 | 25819988 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 3.000 | 100.0% / 3.000 | 51.1 | 17,585 |
| KCNH2 | 26118593 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 72.1 | 18,951 |
| KCNH2 | 26746457 | text | 10 | 32 | 0 | 23.8% | 100.0% | 38.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 2.1 | 1,143 |
| KCNH2 | 29016797 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 50.0% / 1.000 | 0.0% / n/a | 85.6 | 28,780 |
| KCNH2 | 29121719 | text | 0 | 0 | 3 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.9 | 12,222 |
| KCNH2 | 29214556 | text | 1 | 3 | 1 | 25.0% | 50.0% | 33.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 150.2 | 14,520 |
| KCNH2 | 29650123 | text | 1 | 0 | 21 | 100.0% | 4.5% | 8.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 26.0 | 20,741 |
| KCNH2 | 30036649 | text | 5 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 20.0% / 0.000 | 126.6 | 35,704 |
| KCNH2 | 30246897 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.333 | 100.0% / 0.333 | 56.7 | 18,212 |
| KCNQ1 | 11351021 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 54.6 | 15,797 |
| KCNQ1 | 14678125 | text | 35 | 9 | 6 | 79.5% | 85.4% | 82.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 13.2 | 5,803 |
| KCNQ1 | 16155735 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 44.7 | 10,478 |
| KCNQ1 | 18567635 | text | 1 | 9 | 0 | 10.0% | 100.0% | 18.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 238.9 | 36,062 |
| KCNQ1 | 18580685 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 64.3 | 17,840 |
| KCNQ1 | 18808722 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 75.3 | 25,501 |
| KCNQ1 | 19056345 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 65.5 | 14,463 |
| KCNQ1 | 19114714 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 43.3 | 16,875 |
| KCNQ1 | 19632626 | text | 1 | 63 | 0 | 1.6% | 100.0% | 3.1% | 100.0% / 1.000 | 100.0% / 6.000 | 100.0% / 5.000 | 92.8 | 22,502 |
| KCNQ1 | 19687231 | text | 0 | 1 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.0 | 15,525 |
| KCNQ1 | 20348026 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 36.9 | 20,739 |
| KCNQ1 | 20368164 | text | 1 | 2 | 1 | 33.3% | 50.0% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 63.9 | 21,873 |
| KCNQ1 | 20959120 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 36.5 | 16,667 |
| KCNQ1 | 21070882 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 74.6 | 19,514 |
| KCNQ1 | 21129503 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 35.7 | 14,909 |
| KCNQ1 | 21956039 | text | 16 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 68.8% / 0.000 | 0.0% / n/a | 180.6 | 41,606 |
| KCNQ1 | 22613981 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.5 | 7,514 |
| KCNQ1 | 24667783 | text | 13 | 1 | 3 | 92.9% | 81.2% | 86.7% | 75.0% / 0.083 | 43.8% / 0.000 | 0.0% / n/a | 182.0 | 25,682 |
| KCNQ1 | 25471708 | text | 8 | 1 | 0 | 88.9% | 100.0% | 94.1% | 87.5% / 0.143 | 50.0% / 1.000 | 0.0% / n/a | 108.0 | 21,056 |
| KCNQ1 | 26496715 | text | 41 | 6 | 1 | 87.2% | 97.6% | 92.1% | 97.6% / 0.000 | 0.0% / n/a | 0.0% / n/a | 10.1 | 1,605 |
| KCNQ1 | 27114410 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 25.8 | 6,266 |
| KCNQ1 | 28491547 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 3.000 | 100.0% / 0.000 | 100.0% / 3.000 | 62.4 | 17,605 |
| KCNQ1 | 28739325 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 55.3 | 24,870 |
| KCNQ1 | 29677589 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 87.4 | 20,638 |
| KCNQ1 | 29851656 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 43.3 | 11,671 |
| KCNQ1 | 29876285 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 2.000 | 71.2 | 23,440 |
| KCNQ1 | 31293497 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 90.8 | 36,681 |
| KCNQ1 | 31520628 | text | 2 | 2 | 7 | 50.0% | 22.2% | 30.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 25.4 | 14,433 |
| KCNQ1 | 32168391 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 82.0 | 25,999 |
| KCNQ1 | 33082985 | text | 2 | 1 | 1 | 66.7% | 66.7% | 66.7% | 66.7% / 1.000 | 66.7% / 1.500 | 33.3% / 0.000 | 102.3 | 33,587 |
| RYR2 | 15466642 | text | 8 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.250 | 100.0% / 0.250 | 0.0% / n/a | 163.9 | 43,981 |
| RYR2 | 16517285 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 54.7 | 10,382 |
| RYR2 | 17161793 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 42.6 | 11,027 |
| RYR2 | 18285261 | text | 4 | 0 | 1 | 100.0% | 80.0% | 88.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.8 | 14,658 |
| RYR2 | 19216760 | text | 4 | 0 | 1 | 100.0% | 80.0% | 88.9% | 40.0% / 0.500 | 40.0% / 1.500 | 0.0% / n/a | 173.7 | 50,361 |
| RYR2 | 19398417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 29.7 | 7,393 |
| RYR2 | 19398665 | text | 27 | 1 | 0 | 96.4% | 100.0% | 98.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.5 | 19,740 |
| RYR2 | 19926015 | text | 35 | 4 | 5 | 89.7% | 87.5% | 88.6% | 2.5% / 1.000 | 2.5% / 0.000 | 2.5% / 1.000 | 82.8 | 42,748 |
| RYR2 | 20395638 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 67.4 | 11,268 |
| RYR2 | 21652165 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 47.2 | 11,660 |
| RYR2 | 22178870 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 72.4 | 18,962 |
| RYR2 | 22222782 | text | 7 | 0 | 0 | 100.0% | 100.0% | 100.0% | 28.6% / 0.500 | 14.3% / 0.000 | 0.0% / n/a | 206.7 | 36,987 |
| RYR2 | 22677073 | text | 21 | 5 | 4 | 80.8% | 84.0% | 82.4% | 12.0% / 0.000 | 8.0% / 0.000 | 0.0% / n/a | 354.9 | 108,205 |
| RYR2 | 25435091 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 77.8 | 16,195 |
| RYR2 | 25463374 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 66.4 | 30,745 |
| RYR2 | 25500949 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 44.6 | 15,568 |
| RYR2 | 25814417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 24.000 | 100.0% / 44.000 | 49.7 | 22,912 |
| RYR2 | 27114410 | text | 0 | 0 | 9 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.8 | 6,254 |
| RYR2 | 27839804 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 43.9 | 12,006 |
| RYR2 | 28798025 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 66.7% / 0.000 | 66.7% / 0.000 | 0.0% / n/a | 231.2 | 69,822 |
| RYR2 | 29350269 | text | 0 | 0 | 3 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 16.4 | 5,924 |
| RYR2 | 30403697 | text | 21 | 0 | 0 | 100.0% | 100.0% | 100.0% | 52.4% / 0.364 | 28.6% / 0.000 | 10.0% / 0.500 | 357.3 | 110,669 |
| RYR2 | 30546600 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 126.8 | 22,593 |
| RYR2 | 30835254 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 32866913 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 100.0% / 0.000 | 60.2 | 9,970 |
| RYR2 | 33315912 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 2.000 | 100.0% / 5.000 | 67.6 | 27,004 |
| RYR2 | 33536282 | text | 1 | 10 | 0 | 9.1% | 100.0% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 103.2 | 33,209 |
| RYR2 | 33640691 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.0 | 10,100 |
| RYR2 | 34661651 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 0.000 | 25.0% / 0.000 | 0.0% / n/a | 111.8 | 36,155 |
| RYR2 | 35663620 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 43.0 | 14,160 |
| SCN5A | 15671429 | text | 4 | 8 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 160.1 | 28,399 |
| SCN5A | 15828879 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 2.000 | 100.0% / 2.000 | 43.5 | 8,881 |
| SCN5A | 15898185 | text | 1 | 0 | 10 | 100.0% | 9.1% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.6 | 6,263 |
| SCN5A | 16301357 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 67.7 | 11,801 |
| SCN5A | 16929919 | text | 2 | 2 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 103.4 | 15,145 |
| SCN5A | 17675083 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 50.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 52.0 | 25,303 |
| SCN5A | 17971661 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 13.5 | 5,452 |
| SCN5A | 19305639 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 40.8 | 13,357 |
| SCN5A | 19406494 | text | 1 | 8 | 0 | 11.1% | 100.0% | 20.0% | 100.0% / 1.000 | 100.0% / 0.000 | 100.0% / 1.000 | 86.8 | 20,346 |
| SCN5A | 19549036 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 34.6 | 17,012 |
| SCN5A | 20038812 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.9 | 20,596 |
| SCN5A | 20539757 | text | 11 | 14 | 0 | 44.0% | 100.0% | 61.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 138.6 | 29,181 |
| SCN5A | 21193062 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 42.2 | 20,091 |
| SCN5A | 21288276 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 64.9 | 15,370 |
| SCN5A | 22101522 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 41.8 | 20,607 |
| SCN5A | 22677073 | text | 12 | 3 | 0 | 80.0% | 100.0% | 88.9% | 50.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 176.5 | 58,263 |
| SCN5A | 22685113 | text | 10 | 7 | 0 | 58.8% | 100.0% | 74.1% | 90.0% / 0.000 | 10.0% / 0.000 | 0.0% / n/a | 223.8 | 46,963 |
| SCN5A | 22882672 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 7.000 | 100.0% / 2.000 | 0.0% / n/a | 50.9 | 10,855 |
| SCN5A | 22966897 | text | 0 | 3 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.7 | 11,986 |
| SCN5A | 23538678 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 0.000 | 0.0% / n/a | 56.3 | 13,619 |
| SCN5A | 24573164 | text | 20 | 7 | 0 | 74.1% | 100.0% | 85.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 169.2 | 26,175 |
| SCN5A | 24581105 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 69.9 | 14,748 |
| SCN5A | 25171853 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.500 | 100.0% / 0.000 | 100.0% / 1.500 | 90.3 | 12,682 |
| SCN5A | 25236808 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 69.0 | 14,537 |
| SCN5A | 26820365 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 99.5 | 39,840 |
| SCN5A | 28087622 | text | 6 | 5 | 1 | 54.5% | 85.7% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 94.1 | 29,666 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 90.9 | 17,374 |
| SCN5A | 29544605 | text | 6 | 1 | 0 | 85.7% | 100.0% | 92.3% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 132.1 | 62,891 |
| SCN5A | 29709101 | text | 11 | 3 | 0 | 78.6% | 100.0% | 88.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 222.9 | 77,980 |
| SCN5A | 32533946 | text | 83 | 44 | 0 | 65.4% | 100.0% | 79.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 181.8 | 63,585 |

## Errors and representation choices

### KCNH2 PMID 10086971

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: C1117fsX
- Extra predictions: c.3108insG

### KCNH2 PMID 15500450

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 16029385

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 16155735

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P923fsX

### KCNH2 PMID 18675227

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

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

### KCNH2 PMID 20181576

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: P926A
- Count disagreements: K897T affected 1 vs 0 (error +1); P926AfsX14 c.2775insG affected 3 vs 0 (error +3); K897T unaffected 1 vs 3 (error -2); P926AfsX14 c.2775insG unaffected 0 vs 4 (error -4)

### KCNH2 PMID 21308345

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

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

### KCNH2 PMID 22764740

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 23917959

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Gly603Asp c.1808G>A carriers 1 vs 2 (error -1)

### KCNH2 PMID 23936059

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 25819988

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.H562R c.1685A>G affected 6 vs 9 (error -3); p.H562R c.1685A>G unaffected 7 vs 4 (error +3)

### KCNH2 PMID 26118593

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 26746457

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A175D, A188V, A190T, A9V, D982N, G238S, G873S, G903R, L1105S, M291R, P1018A, P1030L, P347S, P967L, R1035Q, R1047L, R148W, R181Q, R252W, R326H, R328C, R394H, R885C, R894H, R912W, R948S, S981G, V1063L, V325M, W927L, Y403C, p.Lys897Thr

### KCNH2 PMID 29016797

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Asn588Lys c.1764C>G carriers 1 vs 0 (error +1); p.Ser631Ala c.1891T>G carriers 3 vs 0 (error +3); p.Asn588Lys c.1764C>G affected 1 vs 0 (error +1)

### KCNH2 PMID 29121719

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: K897T, M645R, R92C

### KCNH2 PMID 29214556

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: E876X
- Extra predictions: R1047L, p.Asp572Asn, p.His558Arg

### KCNH2 PMID 29650123

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A561V, C49Y, F617fsX, G572S, G911fsX, L109fsX, L779P, L987fsX, N633S, Q1046X, R1035fsX, R148W, R176W, R328C, R534C, R744X, R892fsX, S660L, S818P, W412X, Y43C

### KCNH2 PMID 30036649

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 30246897

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Cys566Arg c.1696T>C affected 1 vs 0 (error +1); p.Cys566Arg c.1696T>C unaffected 0 vs 1 (error -1)

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

- Extra predictions: c.1033-2del, c.1251+9C>A, c.1394-1G>T, c.1515-4G>A, c.477+5G>A, p.Ala178Thr c.532G>A, p.Ala341Glu c.1022C>A, p.Ala341Val c.1022C>T, p.Ala344= c.1032G>A, p.Ala344Val c.1031C>T, p.Arg190Gln c.569G>A, p.Arg190Leu c.569G>T, p.Arg243Cys c.727C>T, p.Arg243His c.728G>A, p.Arg380Lys c.1139G>A, p.Arg401fs c.1201dup, p.Arg507Gln c.1520G>A, p.Arg518X c.1552C>T, p.Arg539Trp c.1615C>T, p.Arg555Ser c.1663C>A, p.Arg583His c.1748G>A, p.Arg591His c.1772G>A, p.Arg594Gln c.1781G>A, p.Arg632fs c.1893dup, p.Asp202Asn c.604G>A, p.Cys214Tyr c.641G>A, p.Cys642Phe c.1925G>T, p.Gln107His c.321G>C, p.Gln359Ter c.1075C>T, p.Gln376Ter c.1126C>T, p.Glu449fs c.1343dup, p.Gly168Arg c.502G>A, p.Gly269Ser c.805G>A, p.Gly306Arg c.916G>A, p.Gly314Ser c.940G>A, p.Gly325Arg c.973G>A, p.Gly626fs c.1875dup, p.His258Arg c.773A>G, p.Ile375fs c.1124_1127del, p.Leu266Pro c.797T>C, p.Leu266fs c.796del, p.Leu496fs c.1486_1487del, p.Lys362Arg c.1085A>G, p.Lys526Glu c.1576A>G, p.Met520Arg c.1559T>G, p.Pro343Ala c.1027C>G, p.Q530* c.1588C>T, p.R190W c.568C>T, p.Ser277del c.825CTC[1], p.Ser546Leu c.1637C>T, p.Ser566Phe c.1697C>T, p.Thr311Ile c.932C>T, p.Thr312Ile c.935C>T, p.Thr322Met c.965C>T, p.Trp15Ter c.45G>A, p.Trp305Leu c.914G>T, p.Trp305Ter c.914G>A, p.Tyr111Cys c.332A>G, p.Tyr184Ser c.551A>C, p.Val110Ile c.328G>A, p.Val215Met c.643G>A, p.Val254Met c.760G>A, p.Val280Gly c.839T>G
- Count disagreements: S209P carriers 7 vs 6 (error +1); S209P affected 6 vs 0 (error +6); S209P unaffected 1 vs 6 (error -5)

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

- No scored variant or count disagreement.

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

- Missed gold variants: I328_S330DEL, D454TFS*7, Y522X
- Extra predictions: p.Thr224Met c.671C>T
- Count disagreements: p.Pro73Thr c.217C>A carriers 1 vs 2 (error -1)

### KCNQ1 PMID 25471708

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.572_576del
- Count disagreements: p.Tyr111Cys c.332A>G carriers 1 vs 2 (error -1); p.R190W c.568C>T affected 1 vs 0 (error +1); p.S227del c.828_830del affected 1 vs 0 (error +1); p.S349W c.1046C>G affected 1 vs 0 (error +1); p.Tyr111Cys c.332A>G affected 1 vs 0 (error +1)

### KCNQ1 PMID 26496715

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: 360_361DUPKQ
- Extra predictions: c.1071_1076dupGAAGCA, c.1251+2T>C, c.477+5G>A, c.605-2G>A, c.921+1G>T, p.Gln359_Arg c.1067AGCAGA[3]

### KCNQ1 PMID 27114410

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G584S

### KCNQ1 PMID 28491547

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Val141Met c.421G>A carriers 3 vs 6 (error -3); p.Val141Met c.421G>A unaffected 0 vs 3 (error -3)

### KCNQ1 PMID 28739325

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Asp242Asn c.724G>A carriers 2 vs 0 (error +2)

### KCNQ1 PMID 29677589

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg594Gln c.1781G>A affected 1 vs 0 (error +1); p.R190W c.568C>T affected 1 vs 0 (error +1)

### KCNQ1 PMID 29851656

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.R192L c.575G>T affected 1 vs 0 (error +1)

### KCNQ1 PMID 29876285

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Trp305Ter c.914G>A
- Count disagreements: A300S unaffected 2 vs 4 (error -2)

### KCNQ1 PMID 31293497

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.R360_Q361insQKQR, p.Lys362_His c.1077_1088dup

### KCNQ1 PMID 31520628

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G589D, R366W, Y184S, Y315N, V254M, I235N, R252fsX
- Extra predictions: p.Ala341Val, p.Arg591His c.1772G>A

### KCNQ1 PMID 32168391

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 33082985

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: W248F + L347R
- Extra predictions: p.Trp248Arg c.742T>C
- Count disagreements: p.Leu347Arg c.1040T>G carriers 2 vs 1 (error +1); p.Trp248Phe c.743_744delGGinsTC carriers 2 vs 1 (error +1); p.Leu347Arg c.1040T>G affected 2 vs 0 (error +2); p.Trp248Phe c.743_744delGGinsTC affected 1 vs 0 (error +1)

### RYR2 PMID 15466642

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Ile4848Val carriers 2 vs 4 (error -2); p.Ile4848Val affected 2 vs 4 (error -2)

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
- Count disagreements: p.Asn3308Ser carriers 5 vs 4 (error +1); p.Asn3308Ser affected 1 vs 4 (error -3)

### RYR2 PMID 19398417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 19398665

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: G3946S

### RYR2 PMID 19926015

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: L62F, M81L, EXON 3 DELETION, N4736 DEL, R4822H
- Extra predictions: A4282V, G1886S, M2389L, V4010M
- Count disagreements: p.Y4149S carriers 2 vs 1 (error +1); p.Y4149S unaffected 1 vs 0 (error +1)

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

- Missed gold variants: G4936R, Q2958R, A1136V, R4037C
- Extra predictions: c.13610C>T, c.14803G>A, c.2398+5G>T, c.3321C>T, c.6739C>T

### RYR2 PMID 25435091

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 25463374

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R3570W, S4565R

### RYR2 PMID 25500949

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 25814417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Gly357Ser c.1069G>A affected 97 vs 73 (error +24); p.Gly357Ser c.1069G>A unaffected 62 vs 106 (error -44)

### RYR2 PMID 27114410

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: H240R, T4158P, R420W, S2246L, L3879P, Q3925E, S3959L, W4645R, G4936R

### RYR2 PMID 27839804

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H29D

### RYR2 PMID 28798025

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 29350269

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: Q4164E, E1127G, A3814V

### RYR2 PMID 30403697

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: R2028H carriers 1 vs 3 (error -2); Y4721C carriers 2 vs 3 (error -1); c.14091-11dupT carriers 1 vs 2 (error -1); Y4721C unaffected 1 vs 2 (error -1)

### RYR2 PMID 30546600

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Pro56= c.168G>A

### RYR2 PMID 30835254

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P1124L

### RYR2 PMID 32866913

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Thr85Ile c.254C>T carriers 3 vs 2 (error +1); p.Thr85Ile c.254C>T affected 2 vs 1 (error +1)

### RYR2 PMID 33315912

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Pro2328Ser c.6982C>T affected 17 vs 15 (error +2); p.Pro2328Ser c.6982C>T unaffected 42 vs 47 (error -5)

### RYR2 PMID 33536282

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.A4860G, p.D4112N, p.D4646A, p.I2075T, p.I3995V, p.I4855M, p.K4594R, p.Q4879H, p.S4938F, p.T4196I

### RYR2 PMID 33640691

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 34661651

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 35663620

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15671429

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.2550_2551insTG, c.704-2A>G, p.Asp1274Asn c.3820G>A, p.Asp1594Asn c.4780G>A, p.Asp1594His c.4780G>C, p.Gly1572= c.4716C>T, p.Phe851fs c.2550_2551dup, p.Pro1331fs c.3992del

### SCN5A PMID 15828879

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Asn406Ser c.1217A>G
- Count disagreements: p.Arg282His affected 1 vs 3 (error -2); p.Arg282His unaffected 2 vs 0 (error +2)

### SCN5A PMID 15898185

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F1293S, H558R, P1090L, R1512W, R34C, R481W, S1103Y, S1787N, S524Y, V1951L

### SCN5A PMID 16301357

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 16929919

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: F1760C, V1763M

### SCN5A PMID 17675083

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R34C

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

- Extra predictions: H558R, R1623H, R1632X, T1620M, p.Lys1578fs, p.Phe1617del, c.475_482+3del, p.Arg1622Ter c.4864C>T, p.Arg1631His c.4892G>A, p.Arg1631Ser c.4891C>A, p.Asp1274Asn c.3820G>A, p.Gly1407Arg c.4219G>A, p.Phe1616del c.4844TCT[1], p.Pro1297Leu c.3890C>T

### SCN5A PMID 21193062

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 21288276

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Trp193fs c.576del

### SCN5A PMID 22101522

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Ile1767Val c.5299A>G, p.Pro2005Ala c.6013C>G, p.Val1950Leu c.5848G>T

### SCN5A PMID 22685113

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg1625His c.4874G>A, p.Arg1896Trp c.5686C>T, p.Asp1818Asn c.5452G>A, p.Phe1595Ile c.4783T>A, p.Phe2003Leu c.6007T>C, p.Thr1303Met c.3908C>T, p.Val1950Met c.5848G>A

### SCN5A PMID 22882672

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Pro1177Leu carriers 1 vs 8 (error -7); p.Pro1177Leu affected 1 vs 3 (error -2)

### SCN5A PMID 22966897

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P1177L, R190W
- Extra predictions: I1768V, P2006A, S1103Y

### SCN5A PMID 23538678

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Ser231CysfsTer251 c.692-693delCA carriers 2 vs 1 (error +1)

### SCN5A PMID 24573164

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Leu567Gln, p.Ala1112Val c.3335C>T, p.Asp1274Asn c.3820G>A, p.Glu1937Lys c.5809G>A, p.Gly1318Val c.3953G>T, p.Leu1500Val c.4498C>G, p.Ser1139Asn c.3416G>A

### SCN5A PMID 24581105

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 25171853

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.P1048SfsX97 carriers 3 vs 1 (error +2); p.Thr220Ile c.659C>T carriers 2 vs 1 (error +1); p.P1048SfsX97 unaffected 2 vs 0 (error +2); p.Thr220Ile c.659C>T unaffected 1 vs 0 (error +1)

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
- Extra predictions: R1918H, p.Tyr1795His, p.Arg1912His c.5735G>A, p.Gln1908Arg c.5723A>G, p.Glu1900Gln c.5698G>C

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

- Extra predictions: D1430N, D356N, G1712C, G1740R, G1743E, G1743R, G897E, I1660V, L846R, LEU839P, R104Q, R104W, R878H, R893H, S1218I, S910L, p.Arg282His, p.Arg878Cys, p.Gly1408Arg, p.Thr187Ile, p.Ala1679Thr c.5035G>A, p.Ala1745Val c.5234C>T, p.Arg1582Cys c.4744C>T, p.Arg1631His c.4892G>A, p.Arg1631Ser c.4891C>A, p.Arg1897Cys c.5689C>T, p.Arg1957Ter c.5869C>T, p.Asn1721Asp c.5161A>G, p.Asp1242Asn c.3724G>A, p.Asp1369Gly c.4106A>G, p.Glu1224Lys c.3670G>A, p.Glu1573Lys c.4717G>A, p.Glu1783Lys c.5347G>A, p.Gly1261Ser c.3781G>A, p.Gly1419Ala c.4256G>C, p.Gly1419Cys c.4255G>T, p.Gly1641Glu c.4922G>A, p.Gly1739Arg c.5215G>A, p.Thr1708Lys c.5123C>A, p.Thr1708Met c.5123C>T, p.Tyr1448Cys c.4343A>G, p.Val1250Met c.3748G>A, p.Val1278Ile c.3832G>A, p.Val1352Met c.4054G>A

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
- Phenotype-count recovery figure and inspectable source data: `figures/phenotype_count_recovery.svg`, `figures/phenotype_count_recovery.png`, `figures/phenotype_count_recovery.pdf`, `figures/data/phenotype_count_recovery.csv`, `figures/data/phenotype_count_recovery.json`.
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
