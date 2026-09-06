# Codex extraction-blinded paper evaluation — `20260905_protocol_cont120_03_baseline`

## Technical summary

This hash-locked run evaluated **120 papers** (**per-gene counts {'SCN5A': 57, 'KCNH2': 21, 'KCNQ1': 25, 'RYR2': 14, 'BRCA1': 1, 'BRCA2': 1, 'MYBPC3': 1}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **85.3%**, recall **80.6%**, F1 **82.9%** (773 TP, 133 FP, 186 FN).
- Precision versus counted extras **96.4%** (773 matched rows; 29 extra rows with patient counts). The stricter count-bearing-only diagnostic is **95.3%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **2,333,356 total tokens** (1,700,231 input; 633,125 output).
- Elapsed: **6947.0s wall clock**; 6302.6s summed per-paper route + read time.
- Notation twins merged before scoring: **1** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 120}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 773 | 133 | 186 | 85.3% | 80.6% | 82.9% |
| `linkage_assisted` | secondary_diagnostic | 803 | 610 | 156 | 56.8% | 83.7% | 67.7% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `tranche_03.tsv` (120 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 579 / 959 | 60.4% | 0.316 | 1.173 |
| affected | 134 / 959 | 14.0% | 0.910 | 3.473 |
| unaffected | 139 / 959 | 14.5% | 0.101 | 0.550 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 575 / 813 | 70.7% | 4 / 146 | 2.7% |
| affected | 117 / 710 | 16.5% | 17 / 249 | 6.8% |
| unaffected | 87 / 139 | 62.6% | 52 / 820 | 6.3% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 565 | 77 | 115 | 88.0% | 83.1% | 85.5% | 97.9% | 97.3% | 62.6% / 0.312 / 1.156 | 4.1% / 0.464 / 1.268 | 9.3% / 0.111 / 0.577 |
| KCNH2 | 41 | 18 | 6 | 69.5% | 87.2% | 77.4% | 95.3% | 91.7% | 44.7% / 0.095 / 0.436 | 31.9% / 0.533 / 0.816 | 2.1% / 0.000 / 0.000 |
| KCNQ1 | 127 | 14 | 22 | 90.1% | 85.2% | 87.6% | 96.2% | 95.7% | 73.2% / 0.128 / 0.650 | 56.4% / 1.095 / 4.279 | 46.3% / 0.087 / 0.538 |
| RYR2 | 29 | 24 | 22 | 54.7% | 56.9% | 55.8% | 74.4% | 68.8% | 41.2% / 1.619 / 2.911 | 13.7% / 1.286 / 1.813 | 11.8% / 0.167 / 0.408 |
| BRCA1 | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | n/a | n/a | 0.0% / n/a / n/a | 0.0% / n/a / n/a | 0.0% / n/a / n/a |
| BRCA2 | 9 | 0 | 18 | 100.0% | 33.3% | 50.0% | 100.0% | 100.0% | 3.7% / 0.000 / 0.000 | 0.0% / n/a / n/a | 0.0% / n/a / n/a |
| MYBPC3 | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 100.0% | 100.0% | 33.3% / 0.000 / 0.000 | 0.0% / n/a / n/a | 0.0% / n/a / n/a |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 18822425 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.3 | 13,152 |
| KCNH2 | 19843919 | text | 5 | 1 | 0 | 83.3% | 100.0% | 90.9% | 100.0% / 0.000 | 80.0% / 1.000 | 0.0% / n/a | 106.8 | 43,029 |
| SCN5A | 18060054 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 45.1 | 28,945 |
| KCNQ1 | 24269949 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 100.0% / 4.000 | 100.0% / 15.000 | 0.0% / n/a | 62.2 | 29,315 |
| KCNQ1 | 10220144 | text | 4 | 0 | 3 | 100.0% | 57.1% | 72.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 30.1 | 9,499 |
| KCNQ1 | 25230101 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.4 | 10,825 |
| SCN5A | 15689442 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 12.000 | 100.0% / 1.000 | 0.0% / n/a | 89.4 | 22,893 |
| KCNQ1 | 22885918 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 105.0 | 28,390 |
| KCNH2 | 21499742 | text | 0 | 4 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 73.7 | 12,447 |
| KCNQ1 | 21895724 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 67.5 | 27,610 |
| RYR2 | 23498838 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 46.6 | 11,413 |
| SCN5A | 29017927 | text | 16 | 4 | 5 | 80.0% | 76.2% | 78.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 100.0 | 34,675 |
| KCNQ1 | 34130155 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 20.1 | 13,082 |
| SCN5A | 17697823 | text | 1 | 0 | 6 | 100.0% | 14.3% | 25.0% | 14.3% / 0.000 | 14.3% / 0.000 | 0.0% / n/a | 37.2 | 14,004 |
| SCN5A | 27871843 | text | 1 | 0 | 4 | 100.0% | 20.0% | 33.3% | 20.0% / 0.000 | 20.0% / 0.000 | 0.0% / n/a | 25.7 | 18,997 |
| KCNH2 | 10735633 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 152.3 | 8,718 |
| RYR2 | 24113177 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 2.000 | 100.0% / 2.000 | 100.0% / 0.000 | 90.4 | 28,249 |
| SCN5A | 25483584 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.9 | 13,605 |
| SCN5A | 26412604 | text | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.4 | 5,762 |
| SCN5A | 25261036 | text | 1 | 6 | 0 | 14.3% | 100.0% | 25.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 35.4 | 24,598 |
| SCN5A | 30419068 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 17.7 | 12,760 |
| SCN5A | 19026623 | text | 5 | 0 | 0 | 100.0% | 100.0% | 100.0% | 80.0% / 0.000 | 80.0% / 0.000 | 0.0% / n/a | 67.6 | 30,166 |
| KCNQ1 | 25634836 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 35.1 | 16,204 |
| RYR2 | 21659649 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 51.9 | 24,624 |
| SCN5A | 23995044 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 34.9 | 14,380 |
| SCN5A | 28493952 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 8.000 | 100.0% / 4.000 | 100.0% / 4.000 | 89.5 | 37,607 |
| SCN5A | 25401102 | text | 3 | 6 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 1.000 | 66.7% / 0.500 | 100.0% / 1.000 | 40.2 | 18,311 |
| SCN5A | 12796143 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.9 | 9,308 |
| RYR2 | 25463374 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 67.3 | 31,551 |
| KCNH2 | 25914329 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 53.5 | 27,730 |
| KCNQ1 | 31520628 | text | 0 | 1 | 9 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.0 | 14,341 |
| SCN5A | 29572929 | text | 0 | 1 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 30.6 | 22,951 |
| KCNQ1 | 30337886 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 31.2 | 30,996 |
| SCN5A | 28011106 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 36.6 | 13,204 |
| KCNQ1 | 21459285 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 42.6 | 15,864 |
| SCN5A | 21051419 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.3 | 5,920 |
| KCNH2 | 12062363 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 100.0% / 2.000 | 0.0% / n/a | 38.3 | 11,961 |
| KCNH2 | 19157587 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 19781797 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 0.000 | 50.0% / 4.000 | 0.0% / n/a | 73.9 | 22,166 |
| SCN5A | 22885918 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 24.0 | 14,471 |
| KCNQ1 | 19843919 | text | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 6.4 | 5,538 |
| SCN5A | 28104484 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 33.3% / 0.000 | 33.3% / 0.000 | 0.0% / n/a | 63.2 | 38,698 |
| SCN5A | 28341781 | text | 52 | 0 | 3 | 100.0% | 94.5% | 97.2% | 94.5% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 24372464 | text | 5 | 1 | 1 | 83.3% | 83.3% | 83.3% | 66.7% / 0.500 | 83.3% / 0.800 | 33.3% / 1.000 | 160.1 | 45,354 |
| SCN5A | 22966897 | text | 0 | 3 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.8 | 12,314 |
| KCNH2 | 15500450 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 38.2 | 12,443 |
| KCNH2 | 23995044 | text | 1 | 1 | 1 | 50.0% | 50.0% | 50.0% | 50.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 64.5 | 29,418 |
| SCN5A | 26733869 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 21.0 | 17,770 |
| SCN5A | 15338453 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 26.3 | 8,595 |
| MYBPC3 | 16858239 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 33.3% / 0.000 | 0.0% / n/a | 0.0% / n/a | 50.7 | 12,077 |
| KCNH2 | 14642691 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| BRCA1 | 10441573 | text | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 45.9 | 26,683 |
| KCNH2 | 18675227 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 66.7% / 0.000 | 100.0% / 0.667 | 0.0% / n/a | 79.9 | 30,155 |
| SCN5A | 29540853 | text | 1 | 1 | 2 | 50.0% | 33.3% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 48.1 | 22,993 |
| RYR2 | 28237968 | text | 13 | 10 | 5 | 56.5% | 72.2% | 63.4% | 72.2% / 2.154 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 25616976 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 88.1 | 23,098 |
| SCN5A | 22402334 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 33.2 | 17,010 |
| KCNQ1 | 31565860 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 158.7 | 39,224 |
| KCNH2 | 18551196 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 99.3 | 30,918 |
| KCNH2 | 22885918 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 45.9 | 18,984 |
| SCN5A | 28294644 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 25.5 | 8,958 |
| SCN5A | 16945804 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 43.2 | 13,826 |
| SCN5A | 28781849 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.6 | 12,697 |
| RYR2 | 16918210 | text | 0 | 0 | 4 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 26412604 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 10.4 | 5,594 |
| SCN5A | 20129283 | text | 339 | 7 | 78 | 98.0% | 81.3% | 88.9% | 81.1% / 0.281 | 0.0% / n/a | 12.9% / 0.000 | 8.7 | 2,551 |
| SCN5A | 12569159 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 50.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 118.2 | 30,157 |
| SCN5A | 19377070 | text | 1 | 7 | 0 | 12.5% | 100.0% | 22.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 30.5 | 20,082 |
| SCN5A | 24871449 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.8 | 8,464 |
| RYR2 | 32866913 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 100.0% / 0.000 | 30.2 | 8,842 |
| RYR2 | 27756708 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 25.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 77.3 | 32,455 |
| KCNH2 | 11997281 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 43.9 | 15,123 |
| KCNH2 | 22402334 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 79.2 | 28,874 |
| SCN5A | 23998552 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 24.8 | 23,587 |
| RYR2 | 27225049 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 36.9 | 16,267 |
| KCNH2 | 27871843 | text | 8 | 3 | 1 | 72.7% | 88.9% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 53.8 | 28,920 |
| KCNH2 | 10220144 | text | 5 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 35.9 | 10,224 |
| SCN5A | 30354299 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 10.1 | 6,159 |
| SCN5A | 18451998 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 0.000 | 50.0% / 0.000 | 50.0% / 0.000 | 57.0 | 19,867 |
| SCN5A | 26279430 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 22.6 | 18,486 |
| SCN5A | 14990510 | text | 1 | 6 | 0 | 14.3% | 100.0% | 25.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 93.9 | 26,368 |
| SCN5A | 24059039 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 34.8 | 8,996 |
| KCNQ1 | 26423924 | text | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 10.5 | 6,118 |
| KCNQ1 | 31751991 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.1 | 9,303 |
| SCN5A | 29506689 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 1.000 | 100.0% / 1.000 | 100.0% / 0.000 | 19.7 | 18,039 |
| SCN5A | 32533946 | text | 83 | 20 | 0 | 80.6% | 100.0% | 89.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 95.0 | 63,782 |
| SCN5A | 21705349 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 43.9 | 10,510 |
| SCN5A | 15520322 | text | 3 | 0 | 1 | 100.0% | 75.0% | 85.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 108.5 | 38,193 |
| BRCA2 | 15365993 | text | 9 | 0 | 18 | 100.0% | 33.3% | 50.0% | 3.7% / 0.000 | 0.0% / n/a | 0.0% / n/a | 172.5 | 55,091 |
| RYR2 | 25435091 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 1.000 | 49.5 | 12,295 |
| KCNQ1 | 19841300 | text | 87 | 4 | 3 | 95.6% | 96.7% | 96.1% | 96.7% / 0.046 | 74.4% / 0.000 | 73.3% / 0.061 | 19.9 | 10,819 |
| KCNH2 | 18984536 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 74.0 | 29,549 |
| SCN5A | 11997281 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.2 | 6,980 |
| KCNQ1 | 33552729 | text | 3 | 1 | 0 | 75.0% | 100.0% | 85.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 137.1 | 30,642 |
| SCN5A | 22519808 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 100.0% / 0.000 | 27.7 | 13,999 |
| KCNQ1 | 18266681 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 41.3 | 22,558 |
| KCNQ1 | 28249770 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 30.2 | 10,978 |
| KCNH2 | 25987402 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 42.3 | 14,349 |
| KCNH2 | 26412604 | text | 1 | 5 | 0 | 16.7% | 100.0% | 28.6% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 86.8 | 25,084 |
| SCN5A | 25065297 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.9 | 11,864 |
| SCN5A | 16239976 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 31.2 | 9,728 |
| KCNQ1 | 28491547 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 3.000 | 100.0% / 0.000 | 0.0% / n/a | 48.2 | 17,593 |
| SCN5A | 26667357 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 24.6 | 9,471 |
| KCNH2 | 12808265 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 31.9 | 14,787 |
| SCN5A | 28069705 | text | 6 | 2 | 0 | 75.0% | 100.0% | 85.7% | 83.3% / 0.000 | 83.3% / 0.000 | 0.0% / n/a | 77.8 | 47,632 |
| KCNQ1 | 18713323 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 83.3% / 0.000 | 83.3% / 14.600 | 0.0% / n/a | 126.8 | 42,509 |
| KCNH2 | 19136169 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 35.6 | 9,549 |
| RYR2 | 33536282 | text | 1 | 9 | 0 | 10.0% | 100.0% | 18.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 98.4 | 43,056 |
| SCN5A | 28584071 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 6.000 | 0.0% / n/a | 0.0% / n/a | 29.4 | 15,553 |
| SCN5A | 21552533 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.3 | 18,949 |
| KCNQ1 | 28814790 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 54.8 | 30,344 |
| RYR2 | 12093772 | text | 2 | 0 | 12 | 100.0% | 14.3% | 25.0% | 14.3% / 1.000 | 14.3% / 1.000 | 14.3% / 0.000 | 51.3 | 18,466 |
| SCN5A | 25179549 | text | 1 | 2 | 1 | 33.3% | 50.0% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 91.4 | 36,490 |
| KCNQ1 | 26279191 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 0.000 | 0.0% / n/a | 32.6 | 11,926 |
| SCN5A | 20102920 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 49.8 | 13,946 |
| RYR2 | 31970460 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 235.3 | 15,943 |
| SCN5A | 18156160 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 5.000 | 100.0% / 5.000 | 0.0% / n/a | 229.1 | 16,271 |
| SCN5A | 28087622 | text | 6 | 2 | 1 | 75.0% | 85.7% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 71.3 | 31,221 |
| SCN5A | 19318916 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 38.9 | 12,877 |
| SCN5A | 22155680 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 18.3 | 15,400 |

## Errors and representation choices

### SCN5A PMID 18822425

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 19843919

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Gly604Ser
- Count disagreements: p.Ala614Val affected 1 vs 0 (error +1); p.Asp342Val affected 1 vs 0 (error +1); p.His492Tyr affected 1 vs 0 (error +1); p.Met756Val affected 1 vs 0 (error +1)

### SCN5A PMID 18060054

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 24269949

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: S27A, S27D, T322A, Y315C
- Count disagreements: I235N carriers 19 vs 15 (error +4); I235N affected 15 vs 0 (error +15)

### KCNQ1 PMID 10220144

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G314S, G345R, Y315S

### KCNQ1 PMID 25230101

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: V241F

### SCN5A PMID 15689442

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg1193Gln carriers 12 vs 0 (error +12); p.Arg1193Gln affected 1 vs 0 (error +1)

### KCNQ1 PMID 22885918

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 21499742

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: C566G, H562R
- Extra predictions: G601S, K897fs, R1014*, R176W

### KCNQ1 PMID 21895724

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 23498838

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 29017927

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.F1617DEL, P.K1505_Q1507DEL, P.Y1795_E1796INSD, S216L, T1304M
- Extra predictions: p.Leu409Pro, p.Phe1760Ala, p.Thr304Met, p.Tyr1767Ala

### KCNQ1 PMID 34130155

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17697823

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A735V, L136P, L276Q, T1709M, c.5157delC, c.5290delG

### SCN5A PMID 27871843

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A385T, L1988R, R1629X, R504T

### KCNH2 PMID 10735633

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 24113177

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Met4109Arg c.12326T>G
- Count disagreements: p.Glu2311Asp c.6933G>C carriers 2 vs 4 (error -2); p.Glu2311Asp c.6933G>C affected 1 vs 3 (error -2)

### SCN5A PMID 25483584

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26412604

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P2005A, Q1033R

### SCN5A PMID 25261036

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Asp1275Asn, p.Gly1743Arg, p.Gly292Ser, p.Lys317Asn, p.Phe1473Ser, p.Val294Met

### SCN5A PMID 30419068

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 19026623

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 25634836

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 21659649

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: P2328S

### SCN5A PMID 23995044

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 28493952

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: R1193Q c.3578G>A carriers 8 vs 0 (error +8); R1193Q c.3578G>A affected 4 vs 0 (error +4); R1193Q c.3578G>A unaffected 4 vs 0 (error +4)

### SCN5A PMID 25401102

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A1818G, E1061E, F1750F, L1749L, p.A29A, p.D1818D
- Count disagreements: p.K974D carriers 4 vs 1 (error +3); p.K974D affected 2 vs 1 (error +1); p.K974D unaffected 2 vs 0 (error +2); p.S321Y unaffected 1 vs 0 (error +1)

### SCN5A PMID 12796143

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: S1503A

### RYR2 PMID 25463374

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: S4565R, p.Arg3570Trp

### KCNH2 PMID 25914329

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 31520628

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G589D, I235N, R192fsX, R252fsX, R366W, R591C, V254M, Y184S, Y315N
- Extra predictions: p.Ala341Val

### SCN5A PMID 29572929

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: IVS21+5G>A, V1281X
- Extra predictions: c.3840+5G>A

### KCNQ1 PMID 30337886

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 28011106

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 21459285

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 21051419

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F1596I

### KCNH2 PMID 12062363

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Glu637Lys c.1909G>A carriers 1 vs 3 (error -2); p.Glu637Lys c.1909G>A affected 1 vs 3 (error -2)

### KCNH2 PMID 19157587

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: Q738X

### RYR2 PMID 19781797

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg3570Trp affected 2 vs 6 (error -4)

### SCN5A PMID 22885918

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 19843919

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R231C, R243H

### SCN5A PMID 28104484

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 28341781

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.3840+1G>A, c.4245+1G>A, c.4299+1delG

### KCNQ1 PMID 24372464

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: T377del
- Extra predictions: p.Thr377_Val417del c.1251+1G>T
- Count disagreements: p.Ile567Thr c.1700T>C carriers 5 vs 7 (error -2); p.Asp202Asn c.604G>A affected 1 vs 0 (error +1); p.Ile567Thr c.1700T>C affected 2 vs 0 (error +2); p.Leu273Phe c.817C>T affected 3 vs 2 (error +1); p.Ile567Thr c.1700T>C unaffected 3 vs 5 (error -2)

### SCN5A PMID 22966897

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P1177L, R190W
- Extra predictions: I1768V, P2006A, S1103Y

### KCNH2 PMID 15500450

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 23995044

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: T613L
- Extra predictions: p.Thr612Leu

### SCN5A PMID 26733869

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15338453

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Gly1262Ser c.3934G>A carriers 4 vs 2 (error +2)

### MYBPC3 PMID 16858239

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.1065insC

### KCNH2 PMID 14642691

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R1014X

### BRCA1 PMID 10441573

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.1135insA, c.1675delA

### KCNH2 PMID 18675227

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.L955V c.2863C>G affected 1 vs 0 (error +1); p.R954C c.2860C>T affected 1 vs 0 (error +1)

### SCN5A PMID 29540853

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: L1825P, Y1795insD
- Extra predictions: K1505E

### RYR2 PMID 28237968

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: E1083K, E3716Q, G2866del, M3972I, T4630C
- Extra predictions: c.11147A>G, c.1191G>A, c.12470G>A, c.1258C>T, c.13890G>A, c.14553C>A, c.3407C>T, c.5170G>A, c.5654G>A, c.8598del
- Count disagreements: c.11983A>G carriers 1 vs 7 (error -6); c.13352del carriers 1 vs 3 (error -2); c.13735C>T carriers 1 vs 11 (error -10); c.5614G>A carriers 1 vs 2 (error -1); c.7159G>A carriers 1 vs 3 (error -2); c.7160C>T carriers 1 vs 4 (error -3); c.94C>A carriers 1 vs 5 (error -4)

### KCNQ1 PMID 25616976

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg533Trp

### SCN5A PMID 22402334

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 31565860

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 18551196

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: F640V, p.Ala561Pro c.1681G>C

### KCNH2 PMID 22885918

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 28294644

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 16945804

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 28781849

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 16918210

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P164S, S3938R, T4196A, V186M

### KCNQ1 PMID 26412604

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R555H

### SCN5A PMID 20129283

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1680T, A735V, A735V, A735V, D1243N, D1243N, D1275N, E1225K, E1225K, E1225K, E161K, E746K, E746K, G1319V, G1408R, G1408R, G1661R, G1743E, G1743R, G1743R, G1743R, G752R, I1660V, L1393X, L1393X, L868X, N927S, P.F2004DUP, P.I137_C139DUP, P.I1570DUP, P.K1493DEL, P.Y1795_E1796INSD, P336L, Q1695X, R104Q, R104Q, R104W, R121Q, R1232W, R1232W, R1583C, R1623X, R1638X, R225W, R367H, R367H, R367H, R376H, R376H, R526H, R535X, R535X, R535X, R878H, R878H, R878H, R878H, R893H, R893H, R965C, R965C, S1672Y, T1620M, T1709M, T220I, T632M, V1405M, V232I, c.1890+5G>A, c.1936delC, c.1936delC, c.2435_2436+3delTGGTAinsCGCCT, c.2549_2550insTG, c.2602delC, c.2602delC, c.3667delG, c.3840+1G>A, c.4437+5G>A
- Extra predictions: A1223PfsX7 c.3666delG, F2004dup c.6010_6012dupTTC, I137_C139dup c.410_418dupTCATGTGCA, I1570dup c.4708_4710dupATC, V1525M c.4573G>A, c.2435_24363delTGGTAinsCGCCT, c.5387_5388insTGA
- Count disagreements: A1680T c.5038G>A carriers 2 vs 1 (error +1); A226V carriers 1 vs 2 (error -1); A735V c.2204C>T carriers 4 vs 1 (error +3); D1275N c.3823G>A carriers 3 vs 1 (error +2); E1225K c.3673G>A carriers 4 vs 1 (error +3); E161K c.481G>A carriers 3 vs 1 (error +2); E746K c.2236G>A carriers 3 vs 1 (error +2); G1319V c.3956G>T carriers 5 vs 1 (error +4); G1408R c.4222G>A carriers 7 vs 1 (error +6); G1743E c.5228G>A carriers 6 vs 1 (error +5); G1743R c.5227G>A carriers 5 vs 1 (error +4); G752R c.2254G>A carriers 5 vs 1 (error +4); I1660V carriers 1 vs 4 (error -3); K1493X c.4477A>T carriers 1 vs 2 (error -1); K1493del c.4477_4479delAAG carriers 2 vs 1 (error +1); L1393X c.4178T>A carriers 3 vs 1 (error +2); L868X c.2602delC carriers 2 vs 1 (error +1); N927S c.2780A>G carriers 3 vs 1 (error +2); P336L c.1007C>T carriers 2 vs 1 (error +1); Q1695X c.5083C>T carriers 2 vs 1 (error +1); Q646RfsX5 c.1936delC carriers 3 vs 1 (error +2); R104Q c.311G>A carriers 3 vs 1 (error +2); R104W c.310C>T carriers 2 vs 1 (error +1); R121Q c.362G>A carriers 2 vs 1 (error +1); R1232W c.3694C>T carriers 3 vs 1 (error +2); R1623X c.4867C>T carriers 2 vs 1 (error +1); R1638X c.4912C>T carriers 3 vs 1 (error +2); R225W c.673C>T carriers 3 vs 2 (error +1); R367H c.1100G>A carriers 6 vs 1 (error +5); R376H c.1127G>A carriers 4 vs 1 (error +3); R526H c.1577G>A carriers 2 vs 1 (error +1); R535X c.1603C>T carriers 4 vs 1 (error +3); R878H c.2633G>A carriers 5 vs 1 (error +4); R893H c.2678G>A carriers 3 vs 1 (error +2); R965C c.2893C>T carriers 3 vs 1 (error +2); S1672Y c.5015C>A carriers 2 vs 1 (error +1); T1620M c.4859C>T carriers 2 vs 1 (error +1); T1709M c.5126C>T carriers 2 vs 1 (error +1); T220I c.659C>T carriers 2 vs 1 (error +1); T632M c.1895C>T carriers 2 vs 1 (error +1); V1405M c.4213G>A carriers 2 vs 1 (error +1); V232I c.694G>A carriers 2 vs 1 (error +1); c.1890+5G>A carriers 2 vs 1 (error +1); c.3840+1G>A carriers 6 vs 1 (error +5); c.4437+5G>A carriers 2 vs 1 (error +1)

### SCN5A PMID 12569159

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H588R

### SCN5A PMID 19377070

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: Q1118X, Q55X c.163C>T, R1623X c.4867C>T, R535X c.1603C>T, S1812X c.5435C>A, W1421X, W156X c.468G>A

### SCN5A PMID 24871449

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 32866913

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Thr85Ile c.254C>T carriers 3 vs 2 (error +1); p.Thr85Ile c.254C>T affected 2 vs 1 (error +1)

### RYR2 PMID 27756708

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Ser4938Phe c.14813C>T carriers 1 vs 2 (error -1)

### KCNH2 PMID 11997281

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 22402334

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 23998552

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: N1325S

### RYR2 PMID 27225049

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 27871843

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R823T
- Extra predictions: c.1539C>T, c.1956T>C, c.3069C>G

### KCNH2 PMID 10220144

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 30354299

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.1795INSD

### SCN5A PMID 18451998

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26279430

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 14990510

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: C1850A, F1794C, Y1795A, Y1795F, Y1795S, p.Tyr1795His

### SCN5A PMID 24059039

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Phe1617del

### KCNQ1 PMID 26423924

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A46T, R583H

### KCNQ1 PMID 31751991

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 29506689

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R222Q, p.Ile141Val
- Count disagreements: p.Gly213Asp c.638G>A carriers 14 vs 13 (error +1); p.Gly213Asp c.638G>A affected 13 vs 12 (error +1)

### SCN5A PMID 32533946

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: D1430N, D356N, G1408R, G1712C, G1740R, G897E, I1660V, L846R, LEU839P, R104Q, R104W, R282H, R878C, R878H, R893H, S1218I, S910L, T187I, p.Gly1743Arg, p.Gly1743Glu

### SCN5A PMID 21705349

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15520322

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.4813+5insTGGG

### BRCA2 PMID 15365993

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.1-26G>A, c.2138-51G>T, c.2138-73T>C, c.653+67A>C, c.7070-9T>C, c.7663+10G>C, c.8034-14T>C, c.909+54G>C, c.909+56C>T, c.9877-19G>A, p.A565, p.H743, p.K1132, p.K454, p.N830, p.S2414, p.S455, p.V1269

### RYR2 PMID 25435091

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Cys2277Arg unaffected 0 vs 1 (error -1)

### KCNQ1 PMID 19841300

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: 392INSW, 73-73 DEL AAP, Y315F
- Extra predictions: F335L, P73del, W392ins, Y315S
- Count disagreements: I54ins carriers 5 vs 1 (error +4); I54ins unaffected 5 vs 1 (error +4)

### KCNH2 PMID 18984536

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: G965X, R1014X

### SCN5A PMID 11997281

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F1250L, G615E, H558R, L618F, R34C

### KCNQ1 PMID 33552729

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.477+1G>A

### SCN5A PMID 22519808

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: R1193Q c.3578G>A carriers 2 vs 3 (error -1); R1193Q c.3578G>A affected 1 vs 2 (error -1)

### KCNQ1 PMID 18266681

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: I313M

### KCNQ1 PMID 28249770

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 25987402

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 26412604

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.F513F, p.I489I, p.K897T, p.L564L, p.Y652Y

### SCN5A PMID 25065297

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 16239976

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 28491547

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Val141Met c.421G>A carriers 3 vs 6 (error -3)

### SCN5A PMID 26667357

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 12808265

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 28069705

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: L729D, R190Q

### KCNQ1 PMID 18713323

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: A341V affected 20 vs 0 (error +20); G269S affected 25 vs 0 (error +25); G314S affected 8 vs 0 (error +8); R243C affected 13 vs 0 (error +13); p.Phe340del affected 7 vs 0 (error +7)

### KCNH2 PMID 19136169

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 33536282

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.A4860G, p.D4112N, p.D4646A, p.I3995V, p.I4855M, p.K4594R, p.Q4879H, p.S4938F, p.T4196I

### SCN5A PMID 28584071

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H558R c.1673A>G, R34C, S1103Y
- Count disagreements: R1193Q c.3578G>A carriers 6 vs 0 (error +6)

### SCN5A PMID 21552533

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 28814790

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Val141Met

### RYR2 PMID 12093772

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A4860G, E2311D, E4950K, G3946S, I4867M, L3778F, N4104K, N4895D, R2474S, R4497C, S2246L, V4771I
- Count disagreements: p.Gly3946Ser carriers 2 vs 1 (error +1); p.Ser2246Leu carriers 2 vs 1 (error +1); p.Gly3946Ser affected 2 vs 1 (error +1); p.Ser2246Leu affected 2 vs 1 (error +1)

### SCN5A PMID 25179549

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G1319V
- Extra predictions: G1318A, G1318V

### KCNQ1 PMID 26279191

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Val141Met carriers 1 vs 2 (error -1)

### SCN5A PMID 20102920

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Lys1505_Gln1507del

### RYR2 PMID 31970460

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: D1220E
- Extra predictions: p.Asp3291Val

### SCN5A PMID 18156160

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H558R
- Count disagreements: p.Pro1438Leu c.4313C>T carriers 5 vs 0 (error +5); p.Pro1438Leu c.4313C>T affected 5 vs 0 (error +5)

### SCN5A PMID 28087622

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.K1505_Q1507DEL
- Extra predictions: R1918H, p.Tyr1795His

### SCN5A PMID 19318916

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 22155680

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

## Scope, method, and limitations

- Population: fixed manifest `tranche_03.tsv` (120 papers); per-gene counts {'SCN5A': 57, 'KCNH2': 21, 'KCNQ1': 25, 'RYR2': 14, 'BRCA1': 1, 'BRCA2': 1, 'MYBPC3': 1}; every PMID has downloaded source and at least one gold assertion in each count field.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.
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
