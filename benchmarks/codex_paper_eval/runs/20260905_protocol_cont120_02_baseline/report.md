# Codex extraction-blinded paper evaluation — `20260905_protocol_cont120_02_baseline`

## Technical summary

This hash-locked run evaluated **120 papers** (**per-gene counts {'SCN5A': 57, 'KCNH2': 21, 'KCNQ1': 25, 'RYR2': 14, 'APOE': 1, 'BRCA2': 1, 'MYBPC3': 1}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **70.8%**, recall **77.2%**, F1 **73.9%** (441 TP, 182 FP, 130 FN).
- Precision versus counted extras **96.7%** (441 matched rows; 15 extra rows with patient counts). The stricter count-bearing-only diagnostic is **95.0%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **2,883,141 total tokens** (2,078,165 input; 804,976 output).
- Elapsed: **8376.0s wall clock**; 7424.8s summed per-paper route + read time.
- Notation twins merged before scoring: **1** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 120}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 441 | 182 | 130 | 70.8% | 77.2% | 73.9% |
| `linkage_assisted` | secondary_diagnostic | 463 | 368 | 108 | 55.7% | 81.1% | 66.0% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `tranche_02.tsv` (120 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 282 / 571 | 49.4% | 0.121 | 0.526 |
| affected | 107 / 571 | 18.7% | 0.542 | 2.550 |
| unaffected | 34 / 571 | 6.0% | 1.765 | 7.627 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 279 / 519 | 53.8% | 3 / 52 | 5.8% |
| affected | 92 / 460 | 20.0% | 15 / 111 | 13.5% |
| unaffected | 21 / 88 | 23.9% | 13 / 483 | 2.7% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 199 | 46 | 56 | 81.2% | 78.0% | 79.6% | 96.1% | 93.5% | 44.7% / 0.158 / 0.592 | 22.4% / 0.088 / 0.296 | 7.5% / 0.211 / 0.725 |
| KCNH2 | 108 | 4 | 55 | 96.4% | 66.3% | 78.5% | 100.0% | 100.0% | 54.0% / 0.080 / 0.489 | 6.7% / 0.091 / 0.302 | 0.0% / n/a / n/a |
| KCNQ1 | 81 | 18 | 10 | 81.8% | 89.0% | 85.3% | 93.1% | 91.5% | 67.0% / 0.115 / 0.496 | 29.7% / 1.000 / 2.046 | 7.7% / 1.429 / 2.070 |
| RYR2 | 36 | 32 | 6 | 52.9% | 85.7% | 65.5% | 97.3% | 92.9% | 31.0% / 0.154 / 0.392 | 19.0% / 3.125 / 8.493 | 9.5% / 11.250 / 22.006 |
| APOE | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 / 0.000 | 100.0% / 0.000 / 0.000 | 0.0% / n/a / n/a |
| BRCA2 | 3 | 82 | 0 | 3.5% | 100.0% | 6.8% | 100.0% | n/a | 0.0% / n/a / n/a | 0.0% / n/a / n/a | 0.0% / n/a / n/a |
| MYBPC3 | 13 | 0 | 3 | 100.0% | 81.2% | 89.7% | 100.0% | 100.0% | 31.2% / 0.000 / 0.000 | 18.8% / 0.000 / 0.000 | 25.0% / 0.250 / 0.500 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 21216356 | text | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.2 | 5,528 |
| SCN5A | 28491684 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 52.2 | 17,239 |
| KCNH2 | 29331839 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 57.9 | 21,250 |
| SCN5A | 11748104 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 89.5 | 19,898 |
| RYR2 | 32218223 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.667 | 66.7% / 0.500 | 66.7% / 0.500 | 67.8 | 21,307 |
| RYR2 | 35439358 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 70.5 | 29,552 |
| SCN5A | 17897635 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 41.2 | 23,021 |
| SCN5A | 30371189 | text | 5 | 0 | 0 | 100.0% | 100.0% | 100.0% | 80.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 136.6 | 49,160 |
| SCN5A | 18803136 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 20.9 | 6,900 |
| KCNQ1 | 29672598 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 75.0% / 0.000 | 75.0% / 1.000 | 25.0% / 4.000 | 91.2 | 40,095 |
| KCNQ1 | 29372044 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 68.3 | 17,681 |
| SCN5A | 15851440 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 51.4 | 20,128 |
| RYR2 | 21292648 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 31.2 | 13,448 |
| SCN5A | 17205354 | text | 12 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 111.7 | 18,326 |
| SCN5A | 20038812 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 20.1 | 20,197 |
| KCNQ1 | 18464931 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 14676148 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 10.8 | 10,084 |
| SCN5A | 19843921 | text | 12 | 6 | 0 | 66.7% | 100.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 79.3 | 24,788 |
| SCN5A | 11304498 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 30.3 | 10,587 |
| RYR2 | 28100344 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 35.4 | 11,330 |
| SCN5A | 11274952 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 65.4 | 15,614 |
| SCN5A | 28739862 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 57.7 | 30,177 |
| KCNQ1 | 31899541 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 62.7 | 26,821 |
| KCNQ1 | 32830254 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 0.0% / n/a | 66.7% / 1.000 | 0.0% / n/a | 67.7 | 16,934 |
| KCNQ1 | 29037160 | text | 3 | 2 | 0 | 60.0% | 100.0% | 75.0% | 100.0% / 1.333 | 100.0% / 1.667 | 100.0% / 1.667 | 220.8 | 24,170 |
| SCN5A | 27566755 | text | 47 | 0 | 4 | 100.0% | 92.2% | 95.9% | 88.2% / 0.000 | 0.0% / n/a | 0.0% / n/a | 2.7 | 1,653 |
| KCNQ1 | 25139741 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 54.9 | 7,249 |
| SCN5A | 11150514 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 29.3 | 12,621 |
| RYR2 | 16272262 | text | 12 | 31 | 0 | 27.9% | 100.0% | 43.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 67.5 | 32,368 |
| KCNQ1 | 17038145 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 0.000 | 24.1 | 17,259 |
| SCN5A | 10690282 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 44.6 | 14,481 |
| KCNQ1 | 25786344 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 94.5 | 34,462 |
| APOE | 27108409 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 31.9 | 17,660 |
| SCN5A | 26036855 | text | 2 | 3 | 0 | 40.0% | 100.0% | 57.1% | 50.0% / 1.000 | 50.0% / 1.000 | 0.0% / n/a | 76.9 | 35,958 |
| SCN5A | 12051963 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 40.8 | 15,941 |
| RYR2 | 25814417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 24.000 | 100.0% / 44.000 | 24.2 | 22,291 |
| KCNH2 | 15364333 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.0 | 16,913 |
| KCNH2 | 26496715 | text | 53 | 0 | 1 | 100.0% | 98.1% | 99.1% | 98.1% / 0.000 | 0.0% / n/a | 0.0% / n/a | 2.5 | 1,835 |
| KCNQ1 | 30878014 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 26.3 | 14,892 |
| MYBPC3 | 20433692 | text | 13 | 0 | 3 | 100.0% | 81.2% | 89.7% | 31.2% / 0.000 | 18.8% / 0.000 | 25.0% / 0.250 | 310.5 | 92,009 |
| SCN5A | 28837624 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 88.4 | 37,806 |
| SCN5A | 21596231 | text | 4 | 5 | 0 | 44.4% | 100.0% | 61.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 62.0 | 23,343 |
| KCNQ1 | 20368164 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 30.7 | 19,941 |
| KCNH2 | 21130771 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 4.000 | 0.0% / n/a | 0.0% / n/a | 24.2 | 18,047 |
| SCN5A | 26636822 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 24.1 | 32,548 |
| RYR2 | 28158428 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 56.9 | 17,726 |
| SCN5A | 24112685 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.0 | 9,125 |
| KCNQ1 | 26022593 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 47.7 | 13,052 |
| KCNH2 | 30036649 | text | 5 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 52.2 | 20,370 |
| KCNQ1 | 30244407 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 50.0% / 0.000 | 50.0% / 0.000 | 50.0% / 0.000 | 51.9 | 29,340 |
| RYR2 | 31875585 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 36.7 | 11,473 |
| SCN5A | 24613995 | text | 13 | 9 | 0 | 59.1% | 100.0% | 74.3% | 53.8% / 0.000 | 38.5% / 0.000 | 23.1% / 0.000 | 216.5 | 85,021 |
| SCN5A | 23276942 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.1 | 7,278 |
| KCNQ1 | 23400408 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 100.0% / 1.000 | 0.0% / n/a | 43.7 | 12,229 |
| KCNH2 | 9693036 | text | 0 | 0 | 4 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 13.0 | 9,666 |
| SCN5A | 15851227 | text | 11 | 3 | 0 | 78.6% | 100.0% | 88.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 213.4 | 67,712 |
| KCNH2 | 26847485 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 39.7 | 22,090 |
| KCNH2 | 23237912 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 64.8 | 16,364 |
| BRCA2 | 26848529 | text | 3 | 82 | 0 | 3.5% | 100.0% | 6.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 177.2 | 76,216 |
| SCN5A | 19083750 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 25.4 | 8,228 |
| KCNH2 | 24057343 | text | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.8 | 8,356 |
| KCNH2 | 27761169 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 49.9 | 25,682 |
| SCN5A | 24144883 | text | 2 | 0 | 6 | 100.0% | 25.0% | 40.0% | 12.5% / 0.000 | 0.0% / n/a | 0.0% / n/a | 45.0 | 17,310 |
| SCN5A | 29709101 | text | 11 | 1 | 0 | 91.7% | 100.0% | 95.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 159.7 | 77,200 |
| KCNH2 | 24112685 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 37.4 | 22,394 |
| SCN5A | 19808664 | text | 5 | 1 | 0 | 83.3% | 100.0% | 90.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 51.3 | 24,885 |
| KCNQ1 | 23350853 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 9.000 | 0.0% / n/a | 113.6 | 15,817 |
| KCNQ1 | 16627448 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 2.000 | 100.0% / 1.000 | 28.2 | 11,191 |
| KCNQ1 | 30036649 | text | 3 | 2 | 0 | 60.0% | 100.0% | 75.0% | 100.0% / 0.667 | 33.3% / 0.000 | 0.0% / n/a | 50.9 | 22,374 |
| SCN5A | 21609529 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.6 | 7,382 |
| KCNQ1 | 34135346 | text | 8 | 6 | 2 | 57.1% | 80.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 315.3 | 87,237 |
| SCN5A | 24363796 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 66.7 | 29,483 |
| RYR2 | 30157307 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 33.0 | 11,929 |
| SCN5A | 30036649 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 44.9 | 22,170 |
| RYR2 | 35245853 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 31.7 | 15,714 |
| KCNH2 | 18776039 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 87.7 | 34,463 |
| RYR2 | 18929323 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 41.1 | 16,709 |
| SCN5A | 17442746 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 3.000 | 100.0% / 0.000 | 100.0% / 3.000 | 32.6 | 14,437 |
| SCN5A | 27676163 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 29.4 | 15,055 |
| SCN5A | 15161528 | text | 1 | 0 | 2 | 100.0% | 33.3% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 51.8 | 12,990 |
| SCN5A | 16039271 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.1 | 7,694 |
| SCN5A | 27082542 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 2.000 | 0.0% / n/a | 100.0% / 1.000 | 33.3 | 18,535 |
| KCNQ1 | 32405922 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.4 | 11,660 |
| KCNQ1 | 22250012 | text | 2 | 3 | 0 | 40.0% | 100.0% | 57.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 45.1 | 24,159 |
| KCNQ1 | 18400097 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 110.5 | 27,963 |
| KCNH2 | 29650123 | text | 1 | 0 | 21 | 100.0% | 4.5% | 8.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 26.6 | 21,695 |
| KCNQ1 | 30170673 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 40.1 | 24,751 |
| KCNH2 | 21216356 | text | 0 | 0 | 8 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 10.0 | 5,404 |
| SCN5A | 17198989 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 36.9 | 24,292 |
| RYR2 | 30471092 | text | 3 | 0 | 5 | 100.0% | 37.5% | 54.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.0 | 27,912 |
| RYR2 | 14571276 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 73.2 | 18,452 |
| SCN5A | 21321465 | text | 8 | 2 | 2 | 80.0% | 80.0% | 80.0% | 70.0% / 0.000 | 70.0% / 0.000 | 0.0% / n/a | 118.0 | 46,890 |
| KCNH2 | 30244407 | text | 4 | 1 | 1 | 80.0% | 80.0% | 80.0% | 40.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 95.5 | 41,680 |
| SCN5A | 29672598 | text | 11 | 0 | 1 | 100.0% | 91.7% | 95.7% | 91.7% / 0.273 | 75.0% / 0.222 | 8.3% / 0.000 | 230.3 | 98,303 |
| SCN5A | 22885917 | text | 0 | 0 | 22 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.3 | 11,514 |
| SCN5A | 26496715 | text | 3 | 1 | 1 | 75.0% | 75.0% | 75.0% | 50.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 63.4 | 34,467 |
| SCN5A | 17331104 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.5 | 16,522 |
| SCN5A | 29791480 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 54.9 | 27,852 |
| KCNQ1 | 28491650 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 53.2 | 16,254 |
| SCN5A | 16980337 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 24.7 | 8,267 |
| SCN5A | 24445991 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 1.000 | 100.0% / 1.000 | 100.0% / 0.000 | 33.3 | 13,561 |
| KCNH2 | 30246897 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.333 | 0.0% / n/a | 57.2 | 27,353 |
| SCN5A | 25757662 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 16.7% / 0.000 | 16.7% / 0.000 | 0.0% / n/a | 106.4 | 41,871 |
| SCN5A | 24349418 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 63.8 | 33,738 |
| SCN5A | 23237912 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 42.3 | 16,068 |
| SCN5A | 20137763 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 31.0 | 10,566 |
| SCN5A | 15996170 | text | 1 | 2 | 11 | 33.3% | 8.3% | 13.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.2 | 9,996 |
| SCN5A | 30246897 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.4 | 5,493 |
| KCNH2 | 11844290 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.7 | 5,519 |
| RYR2 | 18285261 | text | 4 | 0 | 1 | 100.0% | 80.0% | 88.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 12.4 | 14,991 |
| KCNQ1 | 26496715 | text | 41 | 5 | 1 | 89.1% | 97.6% | 93.2% | 97.6% / 0.000 | 19.0% / 0.000 | 0.0% / n/a | 282.1 | 111,363 |
| SCN5A | 25236808 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 49.3 | 15,457 |
| KCNH2 | 17171344 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 26.3 | 8,725 |
| SCN5A | 20022821 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 31.0 | 20,831 |
| KCNH2 | 22764740 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 88.4 | 19,601 |
| SCN5A | 19762097 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 9.9 | 11,559 |
| SCN5A | 18452876 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.9 | 9,107 |
| KCNH2 | 11854117 | text | 33 | 0 | 11 | 100.0% | 75.0% | 85.7% | 40.9% / 0.000 | 0.0% / n/a | 0.0% / n/a | 226.4 | 93,697 |
| SCN5A | 15121794 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 4.000 | 100.0% / 0.000 | 0.0% / n/a | 115.2 | 28,033 |
| KCNQ1 | 18441444 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.4 | 9,166 |

## Errors and representation choices

### SCN5A PMID 21216356

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G1329S, H558R

### SCN5A PMID 28491684

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: EXON23_DELETION
- Extra predictions: deletion of exon 23 (bases 3841-3963, NM_198056.2; codons 1281-1321); affects S3-S4/S5 of domain III; likely nonsense-mediated decay and haploinsufficiency

### KCNH2 PMID 29331839

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Pro963Thr carriers 4 vs 2 (error +2)

### SCN5A PMID 11748104

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 32218223

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Leu2432Phe carriers 2 vs 3 (error -1); p.Met4002Val carriers 2 vs 1 (error +1); p.Leu2432Phe affected 1 vs 2 (error -1); p.Met4002Val unaffected 1 vs 0 (error +1)

### RYR2 PMID 35439358

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17897635

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.5464_5467delTCTG
- Extra predictions: H558R

### SCN5A PMID 30371189

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 18803136

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 29672598

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg243Cys c.727C>T affected 1 vs 0 (error +1); p.His455Tyr c.1363C>T affected 1 vs 0 (error +1); p.Thr96Arg c.287C>G affected 1 vs 0 (error +1); p.Glu146Lys c.436G>A unaffected 4 vs 0 (error +4)

### KCNQ1 PMID 29372044

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A512fsX

### SCN5A PMID 15851440

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: R1193Q carriers 2 vs 3 (error -1); R1193Q affected 1 vs 2 (error -1)

### RYR2 PMID 21292648

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17205354

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 20038812

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F1486L

### KCNQ1 PMID 18464931

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: Y315S

### KCNH2 PMID 14676148

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: N588K

### SCN5A PMID 19843921

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.3142_3143insTG, c.3840+1G>A, c.4719C>T, c.934+1G>A, p.Phe861fs*90, p.Trp774fs*28

### SCN5A PMID 11304498

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 28100344

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 11274952

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Val1667Ile carriers 9 vs 11 (error -2)

### SCN5A PMID 28739862

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: L409P

### KCNQ1 PMID 31899541

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 32830254

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: W392X + I145SFSX
- Count disagreements: p.I145Sfs*92 c.431delC affected 2 vs 1 (error +1); p.W392* c.1175G>A affected 1 vs 0 (error +1)

### KCNQ1 PMID 29037160

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R555C, c.477+1G>A
- Count disagreements: R366Q c.1097G>A carriers 8 vs 7 (error +1); p.Arg243His c.728G>A carriers 1 vs 4 (error -3); R366Q c.1097G>A affected 4 vs 0 (error +4); p.Arg243His c.728G>A affected 1 vs 0 (error +1); p.Arg174Cys c.520C>T unaffected 2 vs 0 (error +2); p.Arg243His c.728G>A unaffected 0 vs 3 (error -3)

### SCN5A PMID 27566755

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.F1617DEL, P.I1762DEL, P.K1505_Q1507DEL, P.Q1507_P1509DEL

### KCNQ1 PMID 25139741

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 11150514

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 16272262

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A2387P, A2403T, A4510T, A4607P, A4860G, E2311D, E4146K, F4499C, G3946S, G4671R, I419F, I4848V, I4867M, L3778F, L433P, N2386I, N4097S, N4104K, P164S, P2328S, Q4201R, R176Q, R2474S, R414L, R420W, R4497C, S2246L, T2504M, T4158P, V4653F, Y2392C

### KCNQ1 PMID 17038145

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 10690282

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 25786344

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.A302V c.905C>T affected 1 vs 0 (error +1); p.A46T c.136G>A affected 1 vs 0 (error +1); p.R195W c.583C>T affected 1 vs 0 (error +1); p.R670K c.2009G>A affected 1 vs 0 (error +1)

### APOE PMID 27108409

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26036855

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: GLN1301DEL, c.3901_3903delCAG, p.His558Arg c.1673A>G
- Count disagreements: P1730H carriers 1 vs 2 (error -1); P1730H affected 1 vs 2 (error -1)

### SCN5A PMID 12051963

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H558R

### RYR2 PMID 25814417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.G357S c.1069G>A affected 97 vs 73 (error +24); p.G357S c.1069G>A unaffected 62 vs 106 (error -44)

### KCNH2 PMID 15364333

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.1945+6T>C

### KCNH2 PMID 26496715

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: T443fsX

### KCNQ1 PMID 30878014

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### MYBPC3 PMID 20433692

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: IVS11-9G>A, IVS29+5G>A, IVS6+5G>A
- Count disagreements: D75N unaffected 1 vs 2 (error -1)

### SCN5A PMID 28837624

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 21596231

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: D1275N, D1595H, I1835T, R814W, T220I

### KCNQ1 PMID 20368164

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A340E

### KCNH2 PMID 21130771

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: N588K
- Count disagreements: p.T618I c.1853C>T carriers 4 vs 0 (error +4)

### SCN5A PMID 26636822

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R1623X, R475S c.1425A>C, c.3890_3891insA, p.His558Arg c.1673A>G

### RYR2 PMID 28158428

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Ile4867Val

### SCN5A PMID 24112685

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.T290fsX53

### KCNQ1 PMID 26022593

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.R243H c.128G>A carriers 4 vs 3 (error +1); p.R243H c.128G>A affected 1 vs 0 (error +1)

### KCNH2 PMID 30036649

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 30244407

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F617V

### RYR2 PMID 31875585

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 24613995

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A572D, p.F2004L, p.P1090L, p.P2006A, p.R1193Q, p.S1103Y, p.S216L, p.S524Y, p.V1951L

### SCN5A PMID 23276942

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 23400408

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg243His c.728G>A affected 1 vs 0 (error +1)

### KCNH2 PMID 9693036

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A614V, G572C, N588D, V630A

### SCN5A PMID 15851227

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.F2004L, p.P2006A, p.S216L

### KCNH2 PMID 26847485

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 23237912

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### BRCA2 PMID 26848529

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A634fsX c.1901delC, ALA2603_ARG2659DEL c.7976+1G>A, ASN2101_VAL2102DEL c.6301_6306delAATGTA, D687X c.2059_2063delGATTA, D974N c.2920G>A, E2198fsX c.6591_6592delTG, E3096X c.9286G>T, F2801fsX, G2901D c.8702G>A, G3003X c.9007G>T, H2021P c.6062_6063delinsCA, I1724fsX c.5171delT, I2986fsX c.8956_8957insAA, I332fsX c.993_994dupA, I605fsX c.1813delA, I770fsX c.2307delT, K157fsX c.470_474delAGTCA, K2150fsX c.6449_6450insTA, K3263fsX c.9788delA, K3326X c.9976A>T, K585R c.1754A>G, K585fsX c.1754_1754delA, L1908fsX c.5722_5723delCT, L2080X c.6239T>G, L3119X c.9356_9357delTAinsG, M2393fsX c.7178_7179delTG, M815fsX c.2442delC, N1287fsX c.3860delA, N1473fsX c.4416_4417delGA, N1822fsX c.5465delA, N2051X c.6150_6151insT, N361fsX c.1082delA, N372N c.1114C>A, N863fsX c.2588dupA, N991D c.2971A>G, P1702fsX c.5103dupA, P628fsX c.1881delA, P999L c.2996C>T, Q1056fsX c.3167_3170delAAAA, Q1291X c.3871C>T, Q2655fsX c.7963delC, Q2829X c.8485C>T, Q2941fsX c.8820_8823del, Q3034fsX c.9098_9099insA, Q3036fsX c.9105dup, R2108C c.6322C>T, R2336H c.7007G>A, R435fsX c.1303dupA, S1248fsX c.3744_3747delTGAG, S1722fsX c.5164_5165delAG, S1882X c.5645C>A, S1900X c.5699C>G, S1955X c.5864C>G, S2012fsX c.6033_6034delTT, S2120X c.6359C>G, S2267X c.6800C>A, S2984fsX c.8950delT, T2471fsX c.7409dupT, T2746fsX c.8234dupT, T2867S c.8599A>T, T3085fsX c.9253delA, T598A c.1792A>G, V220fsX c.658_659delGT, W2725fsX c.8172delG, Y1894X c.5682C>G, Y2215X, Y2839X c.8517C>A, c.5332+1G>C, c.5332+1delG, c.5681dupA, c.6645_6648delCTCC, c.7435+53C>T, c.8820_8823delACAA, c.9105dupT, c.994_995dupA, p.Ala938Profs*21 c.2808_2811delACAA, p.Arg155Lysfs*26 c.464_468delGAGAT, p.Ile1859Lysfs*3 c.5576_5579delTTAA, p.Leu88Alafs*12 c.262_263delCT, p.Lys467* c.1399A>T, p.Pro2381Hisfs*13 c.7142delC, p.Ser611* c.1832C>A

### SCN5A PMID 19083750

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 24057343

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: H302fsX, H492Y

### KCNH2 PMID 27761169

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A78G, A78V

### SCN5A PMID 24144883

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F1596I, R1897W, R340Q, T1304M, T220I, V1951M

### SCN5A PMID 29709101

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Ile1660Val c.4978A>G

### KCNH2 PMID 24112685

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 19808664

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: F1486A

### KCNQ1 PMID 23350853

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg231His affected 11 vs 2 (error +9)

### KCNQ1 PMID 16627448

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Ala341Val affected 5 vs 3 (error +2); p.Ala341Val unaffected 1 vs 2 (error -1)

### KCNQ1 PMID 30036649

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.940G>T, c.941G>T
- Count disagreements: p.Ala344= c.1032G>C carriers 4 vs 6 (error -2)

### SCN5A PMID 21609529

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 34135346

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G148R, R387Q
- Extra predictions: c.1394-1G>T, c.1590+1G>A, c.477+5G>A, c.683+5G>A, p.Arg397Gln c.1190G>A, p.Gly168Arg c.502G>A

### SCN5A PMID 24363796

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.5445_5446insT

### RYR2 PMID 30157307

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 30036649

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 35245853

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 18776039

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 18929323

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17442746

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg814Gln c.2441G>A carriers 7 vs 10 (error -3); p.Arg814Gln c.2441G>A unaffected 4 vs 7 (error -3)

### SCN5A PMID 27676163

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15161528

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A29A, D1819D

### SCN5A PMID 16039271

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 27082542

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A29A, E1061E
- Count disagreements: p.R1632C c.4894C>T carriers 2 vs 0 (error +2); p.R1632C c.4894C>T unaffected 1 vs 0 (error +1)

### KCNQ1 PMID 32405922

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: V416M

### KCNQ1 PMID 22250012

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: I145C, S140C, V141C

### KCNQ1 PMID 18400097

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 29650123

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A561V, C49Y, F617fsX, G572S, G911fsX, L109fsX, L779P, L987fsX, N633S, Q1046X, R1035fsX, R148W, R176W, R328C, R534C, R744X, R892fsX, S660L, S818P, W412X, Y43C

### KCNQ1 PMID 30170673

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 21216356

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A193T, E637G, G21D, G572S, G628R, R181W, S1028del, V115M

### SCN5A PMID 17198989

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 30471092

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: 169-?_273+?DEL, H4646D, Q2060H, Q293P, c.169-?_c.273+?del;

### RYR2 PMID 14571276

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 21321465

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P1090L, R1193Q
- Extra predictions: H558R, c.393-1c>t

### KCNH2 PMID 30244407

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F617V
- Extra predictions: D46F

### SCN5A PMID 29672598

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F2004L
- Count disagreements: p.His558Arg c.1673A>G carriers 1 vs 2 (error -1); p.Phe2004Leu c.6010T>C carriers 2 vs 1 (error +1); p.R1193Q c.3578G>A carriers 1 vs 2 (error -1); p.Phe2004Leu c.6010T>C affected 2 vs 1 (error +1); p.R1193Q c.3578G>A affected 1 vs 2 (error -1)

### SCN5A PMID 22885917

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: D1243N, D1741Y, D356N, F1293S, G1319V, G1743E, G514C, G752R, M1335R, P.L729DEL, R1623X, R225W, R367C, T630T, W156X, c.1570_1571insG, c.2582_2583delTT, c.3840+1G>A, c.4118delT, c.5280delG, c.704-1G>C, c.934+1G>A

### SCN5A PMID 26496715

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R1193Q
- Extra predictions: p.Arg1139Gln c.3416G>A

### SCN5A PMID 17331104

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 29791480

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 28491650

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 16980337

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 24445991

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: D1818D
- Count disagreements: A1427S carriers 1 vs 0 (error +1); A1427S affected 1 vs 0 (error +1)

### KCNH2 PMID 30246897

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Cys566Arg c.1696T>C affected 1 vs 0 (error +1)

### SCN5A PMID 25757662

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 24349418

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H558R

### SCN5A PMID 23237912

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 20137763

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15996170

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1148T, A1186T, A1932V, E428K, F532C, H1200Y, P701L, R1739Q, R1913C, R689H, V1667I
- Extra predictions: H558R, c.703+130G>A

### SCN5A PMID 30246897

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G292S

### KCNH2 PMID 11844290

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F805C, M124R, R752W, V822M, W1001X

### RYR2 PMID 18285261

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: K4481R

### KCNQ1 PMID 26496715

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: 360_361DUPKQ
- Extra predictions: c.1071_1076dupGAAGCA, c.1251+2T>C, c.477+5G>A, c.605-2G>A, c.921+1G>T

### SCN5A PMID 25236808

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 17171344

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Gly604Ser c.1810G>A carriers 11 vs 10 (error +1)

### SCN5A PMID 20022821

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Gly752Arg c.2254G>A carriers 3 vs 2 (error +1)

### KCNH2 PMID 22764740

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 19762097

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 18452876

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 11854117

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A83fsX, C44X, G925fsX, I593X, L799SP, P968fsX, Q376SP, R744X, S428X, V295fsX, W1001X

### SCN5A PMID 15121794

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.R1193Q c.3578G>A carriers 5 vs 1 (error +4)

### KCNQ1 PMID 18441444

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: W248F

## Scope, method, and limitations

- Population: fixed manifest `tranche_02.tsv` (120 papers); per-gene counts {'SCN5A': 57, 'KCNH2': 21, 'KCNQ1': 25, 'RYR2': 14, 'APOE': 1, 'BRCA2': 1, 'MYBPC3': 1}; every PMID has downloaded source and at least one gold assertion in each count field.
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
