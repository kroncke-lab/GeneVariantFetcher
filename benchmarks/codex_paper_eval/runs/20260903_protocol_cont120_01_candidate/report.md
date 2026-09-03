# Codex extraction-blinded paper evaluation — `20260903_protocol_cont120_01_candidate`

## Technical summary

This hash-locked run evaluated **120 papers** (**per-gene counts {'SCN5A': 57, 'KCNH2': 21, 'KCNQ1': 25, 'RYR2': 14, 'APOE': 1, 'BRCA1': 1, 'MYBPC3': 1}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **64.9%**, recall **68.0%**, F1 **66.4%** (261 TP, 141 FP, 123 FN).
- Precision versus counted extras **93.2%** (261 matched rows; 19 extra rows with patient counts). The stricter count-bearing-only diagnostic is **87.5%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **2,829,181 total tokens** (2,086,657 input; 742,524 output).
- Elapsed: **9690.0s wall clock**; 8714.0s summed per-paper route + read time.
- Notation twins merged before scoring: **0** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 120}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 261 | 141 | 123 | 64.9% | 68.0% | 66.4% |
| `linkage_assisted` | secondary_diagnostic | 292 | 273 | 92 | 51.7% | 76.0% | 61.5% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `tranche_01.tsv` (120 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 130 / 384 | 33.9% | 0.231 | 0.679 |
| affected | 91 / 383 | 23.8% | 0.253 | 0.503 |
| unaffected | 22 / 383 | 5.7% | 0.227 | 0.477 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 127 / 332 | 38.3% | 3 / 52 | 5.8% |
| affected | 76 / 291 | 26.1% | 15 / 92 | 16.3% |
| unaffected | 10 / 52 | 19.2% | 12 / 331 | 3.6% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 109 | 89 | 48 | 55.1% | 69.4% | 61.4% | 94.0% | 89.9% | 38.2% / 0.317 / 0.806 | 29.3% / 0.196 / 0.442 | 7.0% / 0.091 / 0.302 |
| KCNH2 | 25 | 20 | 34 | 55.6% | 42.4% | 48.1% | 86.2% | 80.0% | 27.1% / 0.062 / 0.250 | 10.3% / 0.500 / 0.707 | 5.2% / 0.667 / 0.816 |
| KCNQ1 | 37 | 16 | 11 | 69.8% | 77.1% | 73.3% | 90.2% | 84.6% | 43.8% / 0.143 / 0.488 | 33.3% / 0.625 / 0.791 | 10.4% / 0.200 / 0.447 |
| RYR2 | 63 | 11 | 14 | 85.1% | 81.8% | 83.4% | 98.4% | 93.8% | 19.5% / 0.467 / 1.000 | 14.3% / 0.091 / 0.302 | 2.6% / 0.500 / 0.707 |
| APOE | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 / 0.000 | 50.0% / 0.000 / 0.000 | 50.0% / 0.000 / 0.000 |
| BRCA1 | 0 | 0 | 13 | 0.0% | 0.0% | 0.0% | n/a | n/a | 0.0% / n/a / n/a | 0.0% / n/a / n/a | 0.0% / n/a / n/a |
| MYBPC3 | 25 | 5 | 3 | 83.3% | 89.3% | 86.2% | 89.3% | 84.2% | 57.1% / 0.000 / 0.000 | 39.3% / 0.000 / 0.000 | 0.0% / n/a / n/a |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 27684715 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.7 | 19,146 |
| SCN5A | 23810369 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 20.1 | 10,809 |
| SCN5A | 26131924 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 17.7 | 18,751 |
| KCNH2 | 18452873 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 63.0 | 28,511 |
| KCNH2 | 19070294 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 24.5 | 7,869 |
| SCN5A | 24096004 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 13.7 | 16,807 |
| RYR2 | 17875969 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 3.000 | 0.0% / n/a | 0.0% / n/a | 24.4 | 19,813 |
| SCN5A | 12068034 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 80.4 | 26,017 |
| KCNQ1 | 9702906 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 25.6 | 8,413 |
| SCN5A | 24972929 | text | 1 | 17 | 0 | 5.6% | 100.0% | 10.5% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 545.9 | 51,517 |
| SCN5A | 16155735 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 100.0% / 1.000 | 0.0% / n/a | 21.7 | 16,530 |
| SCN5A | 22407026 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 66.7% / 1.000 | 0.0% / n/a | 76.4 | 24,565 |
| SCN5A | 21677263 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 18.5 | 7,359 |
| KCNH2 | 22407026 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.6 | 10,422 |
| KCNQ1 | 22677073 | text | 10 | 2 | 1 | 83.3% | 90.9% | 87.0% | 81.8% / 0.000 | 81.8% / 1.000 | 0.0% / n/a | 143.9 | 72,463 |
| SCN5A | 25923670 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 49.2 | 18,316 |
| APOE | 23990795 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 50.0% / 0.000 | 50.0% / 0.000 | 181.1 | 46,555 |
| SCN5A | 12454206 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.0 | 8,664 |
| SCN5A | 20031634 | text | 0 | 0 | 13 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.8 | 5,780 |
| SCN5A | 22129298 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 20.9 | 9,884 |
| KCNH2 | 16969682 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.4 | 15,714 |
| SCN5A | 30218094 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 51.6 | 31,132 |
| RYR2 | 24978818 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.333 | 100.0% / 0.000 | 33.3% / 0.000 | 124.9 | 37,314 |
| KCNQ1 | 31293497 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 77.1 | 36,508 |
| SCN5A | 29514831 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 26.1 | 27,928 |
| SCN5A | 23503384 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.250 | 0.0% / n/a | 0.0% / n/a | 67.1 | 23,059 |
| KCNQ1 | 25444851 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.2 | 19,561 |
| KCNH2 | 12402336 | text | 0 | 2 | 9 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 53.1 | 35,599 |
| SCN5A | 14500339 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.0 | 9,525 |
| KCNH2 | 21185499 | text | 0 | 1 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 65.2 | 18,199 |
| KCNQ1 | 15528464 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.7 | 40,041 |
| KCNQ1 | 22539601 | text | 1 | 0 | 6 | 100.0% | 14.3% | 25.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.2 | 6,849 |
| SCN5A | 18452873 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 57.3 | 26,537 |
| SCN5A | 24581105 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 1236.2 | 13,353 |
| KCNQ1 | 19880070 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 56.9 | 28,729 |
| SCN5A | 19632629 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 100.0% / 1.000 | 100.0% / 0.000 | 100.0% / 1.000 | 66.7 | 21,742 |
| KCNH2 | 14714110 | text | 1 | 7 | 0 | 12.5% | 100.0% | 22.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 163.0 | 20,659 |
| KCNQ1 | 26481773 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 55.1 | 13,885 |
| KCNH2 | 24103226 | text | 0 | 0 | 7 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.9 | 10,871 |
| RYR2 | 23094885 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 1.000 | 79.6 | 22,028 |
| SCN5A | 17227473 | text | 1 | 3 | 3 | 25.0% | 25.0% | 25.0% | 25.0% / 0.000 | 25.0% / 0.000 | 0.0% / n/a | 102.0 | 33,212 |
| KCNH2 | 27555138 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 30.9 | 12,588 |
| SCN5A | 27816319 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.4 | 13,807 |
| SCN5A | 26401487 | text | 2 | 2 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 35.4 | 25,277 |
| SCN5A | 19719504 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 28.5 | 13,598 |
| SCN5A | 17016421 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 25.8 | 10,331 |
| RYR2 | 22677073 | text | 21 | 5 | 4 | 80.8% | 84.0% | 82.4% | 12.0% / 0.000 | 8.0% / 0.000 | 0.0% / n/a | 274.8 | 109,203 |
| SCN5A | 23683917 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 0.500 | 50.0% / 0.000 | 0.0% / n/a | 60.8 | 26,830 |
| SCN5A | 24951663 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.1 | 15,898 |
| KCNQ1 | 20167303 | text | 1 | 2 | 1 | 33.3% | 50.0% | 40.0% | 50.0% / 1.000 | 50.0% / 0.000 | 50.0% / 0.000 | 68.3 | 30,460 |
| KCNQ1 | 29876285 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 1.000 | 59.2 | 24,414 |
| RYR2 | 19398665 | text | 27 | 1 | 0 | 96.4% | 100.0% | 98.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.1 | 21,638 |
| KCNH2 | 12442276 | text | 3 | 3 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 66.7% / 1.000 | 66.7% / 1.000 | 44.8 | 18,069 |
| SCN5A | 23805106 | text | 9 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 65.2 | 33,663 |
| SCN5A | 22675453 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 59.2 | 25,075 |
| KCNQ1 | 18580685 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 38.2 | 17,336 |
| SCN5A | 17445919 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.8 | 10,511 |
| SCN5A | 20451667 | text | 1 | 0 | 2 | 100.0% | 33.3% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.4 | 7,469 |
| MYBPC3 | 21302287 | text | 25 | 5 | 3 | 83.3% | 89.3% | 86.2% | 57.1% / 0.000 | 39.3% / 0.000 | 0.0% / n/a | 266.8 | 105,488 |
| KCNQ1 | 16155735 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 30.8 | 16,679 |
| KCNH2 | 27816319 | text | 3 | 1 | 0 | 75.0% | 100.0% | 85.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 82.9 | 22,311 |
| SCN5A | 21498565 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.3 | 14,884 |
| SCN5A | 26776555 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 28.2 | 25,732 |
| KCNQ1 | 19056345 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 52.4 | 15,121 |
| KCNQ1 | 25825456 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 37.4 | 8,860 |
| KCNH2 | 10790218 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.8 | 11,595 |
| SCN5A | 25163546 | text | 0 | 0 | 20 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 5.9 | 6,261 |
| KCNH2 | 17905336 | text | 0 | 0 | 11 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 26.0 | 10,535 |
| SCN5A | 28491738 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 32.2 | 13,106 |
| SCN5A | 20950709 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 39.7 | 13,668 |
| KCNH2 | 17531263 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 50.8 | 22,382 |
| SCN5A | 17675083 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 2.000 | 100.0% / 0.000 | 0.0% / n/a | 120.5 | 35,006 |
| RYR2 | 30296944 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 75.0 | 16,630 |
| RYR2 | 18463390 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.0 | 7,249 |
| SCN5A | 27401036 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.7 | 17,930 |
| RYR2 | 16843546 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 75.0% / 1.000 | 50.0% / 0.500 | 0.0% / n/a | 72.1 | 25,927 |
| KCNH2 | 16155735 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 57.0 | 22,474 |
| SCN5A | 18599870 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 73.2 | 20,975 |
| SCN5A | 27243970 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 3.000 | 0.0% / n/a | 0.0% / n/a | 16.2 | 7,237 |
| SCN5A | 27707468 | text | 7 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 203.0 | 86,567 |
| SCN5A | 24439875 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 30.8 | 20,207 |
| KCNH2 | 20167303 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 25.0% / 1.000 | 25.0% / 0.000 | 61.3 | 31,065 |
| KCNQ1 | 30462975 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 64.5 | 17,874 |
| KCNQ1 | 21952006 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 31.3 | 31,667 |
| KCNQ1 | 17224687 | text | 2 | 3 | 1 | 40.0% | 66.7% | 50.0% | 66.7% / 0.000 | 66.7% / 0.000 | 66.7% / 0.000 | 40.7 | 16,788 |
| RYR2 | 29447731 | text | 0 | 0 | 9 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 16.2 | 19,455 |
| KCNH2 | 15670565 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 40.0 | 12,771 |
| SCN5A | 22677073 | text | 12 | 0 | 0 | 100.0% | 100.0% | 100.0% | 91.7% / 0.091 | 91.7% / 0.091 | 58.3% / 0.000 | 229.1 | 83,910 |
| KCNH2 | 18774102 | text | 0 | 2 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 43.6 | 19,242 |
| RYR2 | 26018045 | text | 1 | 5 | 0 | 16.7% | 100.0% | 28.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 45.4 | 19,704 |
| SCN5A | 19808432 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 55.2 | 27,305 |
| SCN5A | 29709244 | text | 1 | 33 | 0 | 2.9% | 100.0% | 5.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 102.4 | 32,948 |
| KCNQ1 | 11021476 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 32.9 | 11,278 |
| KCNQ1 | 19160088 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 82.0 | 17,708 |
| KCNH2 | 19160088 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 37.1 | 14,091 |
| SCN5A | 23612926 | text | 2 | 2 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 139.1 | 32,327 |
| SCN5A | 29759671 | text | 11 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 36.4% / 0.000 | 0.0% / n/a | 173.8 | 63,777 |
| KCNQ1 | 27810088 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 78.0 | 29,635 |
| BRCA1 | 10528853 | text | 0 | 0 | 13 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 49.1 | 36,301 |
| SCN5A | 20219395 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.4 | 10,515 |
| RYR2 | 33640691 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 43.8 | 15,281 |
| KCNQ1 | 20981542 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 137.0 | 14,237 |
| RYR2 | 17199967 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 32.9 | 16,422 |
| KCNH2 | 26173150 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 109.0 | 16,809 |
| SCN5A | 8917568 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 55.0 | 36,128 |
| SCN5A | 15051636 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 41.9 | 16,890 |
| KCNQ1 | 29978770 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.5 | 18,552 |
| KCNQ1 | 31484877 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 50.0% / 0.000 | 0.0% / n/a | 58.0 | 21,406 |
| KCNH2 | 19490267 | text | 2 | 3 | 0 | 40.0% | 100.0% | 57.1% | 100.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 73.9 | 27,399 |
| SCN5A | 28018021 | text | 1 | 6 | 1 | 14.3% | 50.0% | 22.2% | 50.0% / 0.000 | 50.0% / 0.000 | 50.0% / 0.000 | 55.8 | 37,880 |
| SCN5A | 17905336 | text | 1 | 0 | 5 | 100.0% | 16.7% | 28.6% | 16.7% / 1.000 | 16.7% / 1.000 | 0.0% / n/a | 35.0 | 13,448 |
| SCN5A | 23791817 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.5 | 20,232 |
| SCN5A | 28217227 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 76.6 | 21,793 |
| KCNQ1 | 19348785 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.0 | 10,075 |
| RYR2 | 24793461 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 24.5 | 11,759 |
| SCN5A | 10532948 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.5 | 16,544 |
| SCN5A | 15671429 | text | 4 | 1 | 0 | 80.0% | 100.0% | 88.9% | 50.0% / 0.500 | 50.0% / 0.000 | 0.0% / n/a | 160.2 | 38,206 |
| RYR2 | 30308486 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 33.3 | 17,938 |
| KCNQ1 | 29488358 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 187.4 | 32,197 |
| SCN5A | 19669871 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 26.7 | 13,997 |

## Errors and representation choices

### SCN5A PMID 27684715

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 23810369

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26131924

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: F1344S

### KCNH2 PMID 18452873

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G189ins
- Extra predictions: c.569_570insGCGCGGGCG

### KCNH2 PMID 19070294

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 24096004

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 17875969

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Asn57_Gly91del deletion of exon 3 (part of intron 2, exon 3, part of intron 3; Alu-mediated) carriers 16 vs 13 (error +3)

### SCN5A PMID 12068034

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: I1756C, Y1767C

### KCNQ1 PMID 9702906

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 24972929

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A180A c.540G>A, D331D c.993C>T, E1053K, G1743R, His558Arg c.38645420T>C, L145L c.435C>T, PHE109DEL c.324_326del, R1193Q, S21G c.61A>G, T113T c.339C>T, T353I, V1951L, V245V c.735G>A, Y308C c.923A>G, c.1-23G>C, c.238-11C>G, c.353-16C>A

### SCN5A PMID 16155735

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Ala551Thr carriers 1 vs 3 (error -2); p.Ala551Thr affected 1 vs 2 (error -1)

### SCN5A PMID 22407026

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg1193Gln affected 1 vs 0 (error +1); p.P1090L affected 1 vs 0 (error +1)

### SCN5A PMID 21677263

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: N1325S

### KCNH2 PMID 22407026

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: L162fsX
- Extra predictions: G643S, p.L163fsX236
- Count disagreements: A302V affected 1 vs 0 (error +1); A341E affected 1 vs 0 (error +1); G269S affected 1 vs 0 (error +1); G584S affected 1 vs 0 (error +1); I198V affected 1 vs 0 (error +1); T312I affected 1 vs 0 (error +1); T587M affected 1 vs 0 (error +1); W379R affected 1 vs 0 (error +1); Y111C affected 1 vs 0 (error +1)

### SCN5A PMID 25923670

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### APOE PMID 23990795

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 12454206

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 20031634

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: E161K, G1408R, G752R, L839P, N1722D, P.1570_F1571INSI, R225W, R535X, S1382I, c.1983_1993dupGGCCCTCAGCG, c.3667delG, c.3963+2T>C, c.612-2A>G

### SCN5A PMID 22129298

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 16969682

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G601S

### SCN5A PMID 30218094

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R222Q, R225P, R225W, p.R814W
- Count disagreements: R219H carriers 1 vs 0 (error +1); R219H affected 1 vs 0 (error +1)

### RYR2 PMID 24978818

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.V1810L c.5428G>C carriers 3 vs 2 (error +1)

### KCNQ1 PMID 31293497

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.R360_Q361insQKQR

### SCN5A PMID 29514831

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 23503384

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.E553* carriers 1 vs 0 (error +1)

### KCNQ1 PMID 25444851

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.I204F, p.S209F, p.Val205Met

### KCNH2 PMID 12402336

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: D501N, E698X, G192fsX, L87P, P347S, P872fsX, R823fsX, Y447X, Y99S
- Extra predictions: G1006fsX, p.A188Gfs*143

### SCN5A PMID 14500339

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Gln1077del

### KCNH2 PMID 21185499

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: EXON1-14DEL, EXON4-14DEL
- Extra predictions: 650 kb heterozygous deletion of exons 4-15 of KCNH2 including 19 other genes (including ABP1); delimited by probes A_16_P18164054 to A_16_P01835269

### KCNQ1 PMID 15528464

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 22539601

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: Q530X, R190W, R192CFS91X, S277del, S349W, Y111C

### SCN5A PMID 18452873

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.A572D c.1715C>A
- Count disagreements: p.E1784K c.5350G>A carriers 1 vs 2 (error -1); p.N275K c.825C>A carriers 1 vs 2 (error -1); p.E1784K c.5350G>A affected 1 vs 2 (error -1); p.N275K c.825C>A affected 1 vs 2 (error -1)

### SCN5A PMID 24581105

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 19880070

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 19632629

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H558R, R34FS, W1421X, p.Arg1193Gln
- Count disagreements: R1195H carriers 2 vs 1 (error +1); R1195H unaffected 1 vs 0 (error +1)

### KCNH2 PMID 14714110

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: N861I, R176W c.526C>T, R56Q, R823W, S818P, V822M, Y667X

### KCNQ1 PMID 26481773

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R594Q

### KCNH2 PMID 24103226

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A57P, G53R, G604S, I571V, IVS9+5G>T, P946fsX, V1038fsX

### RYR2 PMID 23094885

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: R420Q unaffected 6 vs 5 (error +1)

### SCN5A PMID 17227473

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P1090L, R1193Q, c.3840+1G>A
- Extra predictions: deletion of exon 6, p.Arg1192Gln, p.Pro1089Leu

### KCNH2 PMID 27555138

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 27816319

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26401487

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Ala29Ala, p.Glu1061Glu

### SCN5A PMID 19719504

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17016421

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1136V, G4936R, Q2958R, R4037C
- Extra predictions: c.13610C>T, c.14803G>A, c.2398+5G>T, c.3321C>T, c.6739C>T

### SCN5A PMID 23683917

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H558R
- Count disagreements: p.D1741Y carriers 2 vs 1 (error +1)

### SCN5A PMID 24951663

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Ala29Ala, p.Asp1818Asp, p.Ser1102Tyr

### KCNQ1 PMID 20167303

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: H445Y
- Extra predictions: p.Glu146Lys c.436G>A, p.His455Tyr c.1363C>T
- Count disagreements: p.Arg243Cys c.727C>T carriers 2 vs 3 (error -1)

### KCNQ1 PMID 29876285

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: A300S unaffected 3 vs 4 (error -1)

### RYR2 PMID 19398665

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: G3946S

### KCNH2 PMID 12442276

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: G306V, I489I, V630L
- Count disagreements: p.L413P c.1421T>C affected 1 vs 2 (error -1); p.L559H c.1676T>A affected 1 vs 2 (error -1); p.L413P c.1421T>C unaffected 1 vs 0 (error +1); p.L559H c.1676T>A unaffected 1 vs 0 (error +1)

### SCN5A PMID 23805106

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 22675453

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: C373F

### KCNQ1 PMID 18580685

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17445919

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 20451667

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: H558R, S1103Y

### MYBPC3 PMID 21302287

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.2258_2259insT, c.3192_3193insC, c.506-12delC
- Extra predictions: F244F c.732C>T, K754E, K754EfsX78, L517M c.1549C>A, p.K1065QfsX11 c.3192-3193InsC

### KCNQ1 PMID 16155735

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 27816319

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg518*

### SCN5A PMID 21498565

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26776555

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: L812Q, R811H, R814Q
- Count disagreements: p.Lys817Glu c.2449A>G carriers 1 vs 2 (error -1); p.Lys817Glu c.2449A>G affected 1 vs 2 (error -1)

### KCNQ1 PMID 19056345

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 25825456

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 10790218

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G189ins

### SCN5A PMID 25163546

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A572D, D1114N, F1596I, F2004L, G213D, I1448L, L618F, N70K, P2006A, Q573X, R1898H, R18W, S216L, T1620M, T220I, T220N, V1340I, V294M, Y205X, c.704-2A>G

### KCNH2 PMID 17905336

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: D774Y, F627fsX, G1031fsX, G572S, I42N, L171ins, N629T, P1075L, R328fsX, R534C, W568C

### SCN5A PMID 28491738

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 20950709

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 17531263

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Pro347Ser c.1039C>T carriers 10 vs 9 (error +1)

### SCN5A PMID 17675083

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg34Cys
- Count disagreements: H558R carriers 5 vs 1 (error +4)

### RYR2 PMID 30296944

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 18463390

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P2404T

### SCN5A PMID 27401036

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F1760A

### RYR2 PMID 16843546

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Phe4851Cys carriers 1 vs 3 (error -2); p.Ser2246Leu carriers 1 vs 2 (error -1); p.Ser2246Leu affected 1 vs 2 (error -1)

### KCNH2 PMID 16155735

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 18599870

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 27243970

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: A572D carriers 3 vs 0 (error +3)

### SCN5A PMID 27707468

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 24439875

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Gly1297Glyfs*22

### KCNH2 PMID 20167303

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Q81H c.243G>C affected 1 vs 0 (error +1)

### KCNQ1 PMID 30462975

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G45S
- Extra predictions: p.R16R c.47G>A

### KCNQ1 PMID 21952006

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 17224687

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R59RP
- Extra predictions: S546S, Y662Y, p.Arg594Pro c.1781G>C

### RYR2 PMID 29447731

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: 169-198_273+823DEL, 1827+140_1961+426DEL, 4299+1DELG, E4431K, F4016Y, K2716I, M3695T, R4114Q, R4608W

### KCNH2 PMID 15670565

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: S216L carriers 2 vs 1 (error +1); S216L affected 2 vs 1 (error +1)

### KCNH2 PMID 18774102

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: EX6-14DEL, EX9-14DUP
- Extra predictions: deletion of exons 6-14, duplication of exons 9-14

### RYR2 PMID 26018045

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: L433P, R176Q, R420W, T2504M, p.Ala77Val

### SCN5A PMID 19808432

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.R1623Q

### SCN5A PMID 29709244

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A1656P, D1243N, D501G, E161K, E19X, G1262S, G1319V, G1743E, H445D, I1643N, L1346P, L1373X, L1582P, L618F, L842R, M1701T, N109K, N927S, R1232Q, R1306C, R1629Q, R1897W, R367C, R367H, R808C, R811L, V1667I, W879R, p.Arg1193Gln, p.D1275N, p.R18Q, p.R893H, p.T220I

### KCNQ1 PMID 11021476

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg259Cys affected 1 vs 0 (error +1)

### KCNQ1 PMID 19160088

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 19160088

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 23612926

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg1414Cys c.4240C>T, p.Arg965Cys c.2893C>T

### SCN5A PMID 29759671

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 27810088

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### BRCA1 PMID 10528853

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.1186A>G, c.172T>C, c.185delAG, c.2196G>A, c.2201C>T, c.2430T>C, c.2434T>C, c.2731C>T, c.3107A>T, c.3232A>G, c.3827T>G, c.5382insC, c.546G>T

### SCN5A PMID 20219395

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R282H

### RYR2 PMID 33640691

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 20981542

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Gly316Glu carriers 10 vs 8 (error +2)

### RYR2 PMID 17199967

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 26173150

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 8917568

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.K1505_Q1507DEL

### SCN5A PMID 15051636

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 29978770

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 31484877

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 19490267

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: T474I, V630L, p.Ala614Val

### SCN5A PMID 28018021

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: V1279I
- Extra predictions: R222Q, R225W, c.1292A>G, c.335T>C, p.R219H, p.Val1278Ile c.3832G>A

### SCN5A PMID 17905336

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: I1768V, N1325S, Q692K, R190Q, c.4299_4300insG
- Count disagreements: p.Arg1193Gln carriers 2 vs 3 (error -1); p.Arg1193Gln affected 2 vs 3 (error -1)

### SCN5A PMID 23791817

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 28217227

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 19348785

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 24793461

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 10532948

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg1232Trp

### SCN5A PMID 15671429

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.2550-2551insTG
- Count disagreements: T220I carriers 2 vs 1 (error +1)

### RYR2 PMID 30308486

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 29488358

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg231His, p.Ser140Gly, p.Val241Phe, p.Val307Leu

### SCN5A PMID 19669871

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

## Scope, method, and limitations

- Population: fixed manifest `tranche_01.tsv` (120 papers); per-gene counts {'SCN5A': 57, 'KCNH2': 21, 'KCNQ1': 25, 'RYR2': 14, 'APOE': 1, 'BRCA1': 1, 'MYBPC3': 1}; every PMID has downloaded source and at least one gold assertion in each count field.
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
