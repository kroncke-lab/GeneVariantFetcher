# Codex extraction-blinded paper evaluation — `20260905_protocol_cont120_02_candidate`

## Technical summary

This hash-locked run evaluated **120 papers** (**per-gene counts {'SCN5A': 57, 'KCNH2': 21, 'KCNQ1': 25, 'RYR2': 14, 'APOE': 1, 'BRCA2': 1, 'MYBPC3': 1}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **71.2%**, recall **77.9%**, F1 **74.4%** (445 TP, 180 FP, 126 FN).
- Precision versus counted extras **83.8%** (445 matched rows; 86 extra rows with patient counts). The stricter count-bearing-only diagnostic is **77.6%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **2,720,749 total tokens** (1,960,361 input; 760,388 output).
- Elapsed: **7845.0s wall clock**; 7172.4s summed per-paper route + read time.
- Notation twins merged before scoring: **2** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 120}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 445 | 180 | 126 | 71.2% | 77.9% | 74.4% |
| `linkage_assisted` | secondary_diagnostic | 467 | 365 | 104 | 56.1% | 81.8% | 66.6% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `tranche_02.tsv` (120 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 293 / 571 | 51.3% | 0.113 | 0.551 |
| affected | 122 / 571 | 21.4% | 0.377 | 2.214 |
| unaffected | 41 / 571 | 7.2% | 1.415 | 6.921 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 290 / 519 | 55.9% | 3 / 52 | 5.8% |
| affected | 106 / 460 | 23.0% | 16 / 111 | 14.4% |
| unaffected | 26 / 88 | 29.5% | 15 / 483 | 3.1% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 201 | 49 | 54 | 80.4% | 78.8% | 79.6% | 97.1% | 95.4% | 48.6% / 0.113 / 0.492 | 27.5% / 0.086 / 0.293 | 8.6% / 0.136 / 0.477 |
| KCNH2 | 109 | 4 | 54 | 96.5% | 66.9% | 79.0% | 100.0% | 100.0% | 54.6% / 0.079 / 0.486 | 8.0% / 0.154 / 0.392 | 1.8% / 0.000 / 0.000 |
| KCNQ1 | 82 | 21 | 9 | 79.6% | 90.1% | 84.5% | 94.3% | 92.9% | 67.0% / 0.115 / 0.587 | 31.9% / 0.448 / 0.670 | 9.9% / 1.000 / 1.528 |
| RYR2 | 36 | 32 | 6 | 52.9% | 85.7% | 65.5% | 97.3% | 92.3% | 28.6% / 0.083 / 0.289 | 21.4% / 2.778 / 8.007 | 11.9% / 9.000 / 19.682 |
| APOE | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 / 0.000 | 100.0% / 0.000 / 0.000 | 0.0% / n/a / n/a |
| BRCA2 | 3 | 74 | 0 | 3.9% | 100.0% | 7.5% | 3.9% | 3.9% | 100.0% / 1.333 / 2.309 | 0.0% / n/a / n/a | 0.0% / n/a / n/a |
| MYBPC3 | 13 | 0 | 3 | 100.0% | 81.2% | 89.7% | 100.0% | 100.0% | 18.8% / 0.000 / 0.000 | 0.0% / n/a / n/a | 12.5% / 0.500 / 0.707 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 21216356 | text | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.9 | 5,814 |
| SCN5A | 28491684 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 48.5 | 17,247 |
| KCNH2 | 29331839 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 0.0% / n/a | 100.0% / 0.000 | 40.1 | 17,766 |
| SCN5A | 11748104 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 35.8 | 16,543 |
| RYR2 | 32218223 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.333 | 66.7% / 0.500 | 66.7% / 0.000 | 72.8 | 20,111 |
| RYR2 | 35439358 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 1.000 | 58.5 | 25,281 |
| SCN5A | 17897635 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 52.0 | 25,114 |
| SCN5A | 30371189 | text | 5 | 0 | 0 | 100.0% | 100.0% | 100.0% | 80.0% / 0.000 | 80.0% / 0.000 | 0.0% / n/a | 147.7 | 49,536 |
| SCN5A | 18803136 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 20.2 | 7,257 |
| KCNQ1 | 29672598 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 25.0% / 4.000 | 108.2 | 41,562 |
| KCNQ1 | 29372044 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 45.9 | 15,584 |
| SCN5A | 15851440 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 56.6 | 20,491 |
| RYR2 | 21292648 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 22.9 | 12,922 |
| SCN5A | 17205354 | text | 12 | 1 | 0 | 92.3% | 100.0% | 96.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 96.3 | 17,392 |
| SCN5A | 20038812 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 29.6 | 20,956 |
| KCNQ1 | 18464931 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 14676148 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 41.4 | 9,817 |
| SCN5A | 19843921 | text | 12 | 6 | 0 | 66.7% | 100.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 86.1 | 25,143 |
| SCN5A | 11304498 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 33.7 | 11,019 |
| RYR2 | 28100344 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 32.7 | 11,246 |
| SCN5A | 11274952 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 100.0% / 0.000 | 100.0% / 2.000 | 43.5 | 14,037 |
| SCN5A | 28739862 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 45.7 | 29,538 |
| KCNQ1 | 31899541 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 50.6 | 23,223 |
| KCNQ1 | 32830254 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 0.0% / n/a | 66.7% / 0.500 | 0.0% / n/a | 54.1 | 16,121 |
| KCNQ1 | 29037160 | text | 3 | 2 | 0 | 60.0% | 100.0% | 75.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 151.7 | 39,795 |
| SCN5A | 27566755 | text | 47 | 0 | 4 | 100.0% | 92.2% | 95.9% | 92.2% / 0.000 | 0.0% / n/a | 0.0% / n/a | 3.4 | 1,653 |
| KCNQ1 | 25139741 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 44.1 | 7,845 |
| SCN5A | 11150514 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.7 | 12,700 |
| RYR2 | 16272262 | text | 12 | 31 | 0 | 27.9% | 100.0% | 43.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 70.0 | 34,156 |
| KCNQ1 | 17038145 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 0.000 | 20.2 | 15,924 |
| SCN5A | 10690282 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 35.8 | 9,352 |
| KCNQ1 | 25786344 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 100.0% / 1.000 | 396.4 | 35,463 |
| APOE | 27108409 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 30.5 | 17,789 |
| SCN5A | 26036855 | text | 2 | 4 | 0 | 33.3% | 100.0% | 50.0% | 50.0% / 1.000 | 50.0% / 1.000 | 0.0% / n/a | 72.1 | 39,164 |
| SCN5A | 12051963 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 133.1 | 15,420 |
| RYR2 | 25814417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 24.000 | 100.0% / 44.000 | 23.2 | 22,623 |
| KCNH2 | 15364333 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.3 | 14,979 |
| KCNH2 | 26496715 | text | 53 | 0 | 1 | 100.0% | 98.1% | 99.1% | 98.1% / 0.000 | 0.0% / n/a | 0.0% / n/a | 4.2 | 2,195 |
| KCNQ1 | 30878014 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 28.7 | 15,081 |
| MYBPC3 | 20433692 | text | 13 | 0 | 3 | 100.0% | 81.2% | 89.7% | 18.8% / 0.000 | 0.0% / n/a | 12.5% / 0.500 | 281.3 | 51,116 |
| SCN5A | 28837624 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 66.7% / 0.000 | 0.0% / n/a | 87.8 | 38,410 |
| SCN5A | 21596231 | text | 4 | 5 | 0 | 44.4% | 100.0% | 61.5% | 50.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 90.7 | 32,732 |
| KCNQ1 | 20368164 | text | 1 | 2 | 1 | 33.3% | 50.0% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 51.9 | 21,268 |
| KCNH2 | 21130771 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 4.000 | 0.0% / n/a | 0.0% / n/a | 29.2 | 18,159 |
| SCN5A | 26636822 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 26.9 | 33,171 |
| RYR2 | 28158428 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 34.8 | 18,353 |
| SCN5A | 24112685 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 75.8 | 22,237 |
| KCNQ1 | 26022593 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 37.6 | 12,547 |
| KCNH2 | 30036649 | text | 5 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 84.9 | 36,129 |
| KCNQ1 | 30244407 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 50.0% / 0.000 | 50.0% / 0.000 | 50.0% / 0.000 | 32.2 | 23,993 |
| RYR2 | 31875585 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 62.4 | 10,877 |
| SCN5A | 24613995 | text | 13 | 9 | 0 | 59.1% | 100.0% | 74.3% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 75.0 | 38,030 |
| SCN5A | 23276942 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.0 | 7,774 |
| KCNQ1 | 23400408 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 45.3 | 12,037 |
| KCNH2 | 9693036 | text | 0 | 0 | 4 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 10.9 | 9,469 |
| SCN5A | 15851227 | text | 11 | 3 | 0 | 78.6% | 100.0% | 88.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 203.3 | 51,584 |
| KCNH2 | 26847485 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 35.1 | 22,129 |
| KCNH2 | 23237912 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 99.4 | 18,660 |
| BRCA2 | 26848529 | text | 3 | 74 | 0 | 3.9% | 100.0% | 7.5% | 100.0% / 1.333 | 0.0% / n/a | 0.0% / n/a | 39.3 | 13,220 |
| SCN5A | 19083750 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 17.2 | 7,020 |
| KCNH2 | 24057343 | text | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.7 | 10,068 |
| KCNH2 | 27761169 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 46.5 | 25,646 |
| SCN5A | 24144883 | text | 2 | 0 | 6 | 100.0% | 25.0% | 40.0% | 12.5% / 0.000 | 0.0% / n/a | 0.0% / n/a | 57.3 | 19,162 |
| SCN5A | 29709101 | text | 11 | 1 | 0 | 91.7% | 100.0% | 95.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 146.4 | 76,034 |
| KCNH2 | 24112685 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 30.2 | 21,144 |
| SCN5A | 19808664 | text | 5 | 1 | 0 | 83.3% | 100.0% | 90.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 71.7 | 27,141 |
| KCNQ1 | 23350853 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 77.1 | 16,291 |
| KCNQ1 | 16627448 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 100.0% / 1.000 | 42.9 | 14,541 |
| KCNQ1 | 30036649 | text | 3 | 2 | 0 | 60.0% | 100.0% | 75.0% | 100.0% / 0.667 | 33.3% / 0.000 | 33.3% / 0.000 | 35.3 | 20,583 |
| SCN5A | 21609529 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 21.0 | 7,647 |
| KCNQ1 | 34135346 | text | 8 | 6 | 2 | 57.1% | 80.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 229.6 | 81,339 |
| SCN5A | 24363796 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 66.1 | 29,397 |
| RYR2 | 30157307 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 28.9 | 12,808 |
| SCN5A | 30036649 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 38.2 | 21,235 |
| RYR2 | 35245853 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 44.1 | 15,243 |
| KCNH2 | 18776039 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 45.9 | 34,422 |
| RYR2 | 18929323 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.0 | 15,990 |
| SCN5A | 17442746 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 34.4 | 13,338 |
| SCN5A | 27676163 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 27.1 | 15,437 |
| SCN5A | 15161528 | text | 3 | 1 | 0 | 75.0% | 100.0% | 85.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 36.4 | 10,826 |
| SCN5A | 16039271 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.4 | 8,689 |
| SCN5A | 27082542 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 2.000 | 0.0% / n/a | 100.0% / 1.000 | 41.4 | 19,460 |
| KCNQ1 | 32405922 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 12.1 | 12,006 |
| KCNQ1 | 22250012 | text | 2 | 3 | 0 | 40.0% | 100.0% | 57.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 48.3 | 24,145 |
| KCNQ1 | 18400097 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 112.6 | 27,568 |
| KCNH2 | 29650123 | text | 1 | 0 | 21 | 100.0% | 4.5% | 8.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.2 | 20,495 |
| KCNQ1 | 30170673 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.1 | 23,694 |
| KCNH2 | 21216356 | text | 0 | 0 | 8 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.5 | 5,620 |
| SCN5A | 17198989 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 45.7 | 25,910 |
| RYR2 | 30471092 | text | 3 | 0 | 5 | 100.0% | 37.5% | 54.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.2 | 27,477 |
| RYR2 | 14571276 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 87.6 | 18,662 |
| SCN5A | 21321465 | text | 8 | 2 | 2 | 80.0% | 80.0% | 80.0% | 70.0% / 0.000 | 60.0% / 0.000 | 10.0% / 0.000 | 139.9 | 48,981 |
| KCNH2 | 30244407 | text | 5 | 1 | 0 | 83.3% | 100.0% | 90.9% | 40.0% / 0.000 | 20.0% / 0.000 | 20.0% / 0.000 | 77.4 | 42,519 |
| SCN5A | 29672598 | text | 11 | 0 | 1 | 100.0% | 91.7% | 95.7% | 83.3% / 0.300 | 91.7% / 0.273 | 8.3% / 0.000 | 257.5 | 102,001 |
| SCN5A | 22885917 | text | 0 | 0 | 22 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.0 | 11,724 |
| SCN5A | 26496715 | text | 3 | 1 | 1 | 75.0% | 75.0% | 75.0% | 50.0% / 0.000 | 25.0% / 0.000 | 0.0% / n/a | 81.1 | 36,030 |
| SCN5A | 17331104 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 24.5 | 17,637 |
| SCN5A | 29791480 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 51.5 | 21,080 |
| KCNQ1 | 28491650 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 100.0% / 1.000 | 0.0% / n/a | 44.3 | 16,572 |
| SCN5A | 16980337 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.2 | 7,542 |
| SCN5A | 24445991 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 32.2 | 13,893 |
| KCNH2 | 30246897 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.333 | 0.0% / n/a | 67.2 | 29,182 |
| SCN5A | 25757662 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 44.6 | 21,714 |
| SCN5A | 24349418 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 54.4 | 28,099 |
| SCN5A | 23237912 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 41.9 | 15,368 |
| SCN5A | 20137763 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 30.0 | 10,430 |
| SCN5A | 15996170 | text | 1 | 2 | 11 | 33.3% | 8.3% | 13.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 31.0 | 9,469 |
| SCN5A | 30246897 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 9.8 | 5,789 |
| KCNH2 | 11844290 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.4 | 5,671 |
| RYR2 | 18285261 | text | 4 | 0 | 1 | 100.0% | 80.0% | 88.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 12.3 | 14,843 |
| KCNQ1 | 26496715 | text | 41 | 5 | 1 | 89.1% | 97.6% | 93.2% | 97.6% / 0.000 | 26.2% / 0.000 | 0.0% / n/a | 245.0 | 111,534 |
| SCN5A | 25236808 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.9 | 11,635 |
| KCNH2 | 17171344 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 38.3 | 9,934 |
| SCN5A | 20022821 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 38.8 | 22,057 |
| KCNH2 | 22764740 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 62.3 | 19,995 |
| SCN5A | 19762097 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 12.0 | 11,747 |
| SCN5A | 18452876 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 16.4 | 8,422 |
| KCNH2 | 11854117 | text | 33 | 0 | 11 | 100.0% | 75.0% | 85.7% | 40.9% / 0.000 | 0.0% / n/a | 0.0% / n/a | 208.3 | 91,752 |
| SCN5A | 15121794 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 4.000 | 100.0% / 0.000 | 0.0% / n/a | 85.2 | 27,590 |
| KCNQ1 | 18441444 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 18.4 | 8,747 |

## Errors and representation choices

### SCN5A PMID 21216356

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G1329S, H558R

### SCN5A PMID 28491684

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: EXON23_DELETION
- Extra predictions: deletion of exon 23 (bases 3841-3963, NM_198056.2; codons 1281-1321; affects S3-S4/S5 of domain III)

### KCNH2 PMID 29331839

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Pro963Thr carriers 4 vs 2 (error +2)

### SCN5A PMID 11748104

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 32218223

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Leu2432Phe carriers 2 vs 3 (error -1); p.Leu2432Phe affected 1 vs 2 (error -1)

### RYR2 PMID 35439358

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: R4790Ter unaffected 3 vs 4 (error -1)

### SCN5A PMID 17897635

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H558R

### SCN5A PMID 30371189

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 18803136

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 29672598

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Glu146Lys c.436G>A carriers 5 vs 1 (error +4); p.Arg243Cys c.727C>T affected 1 vs 0 (error +1); p.Glu146Lys c.436G>A affected 1 vs 0 (error +1); p.His455Tyr c.1363C>T affected 1 vs 0 (error +1); p.Thr96Arg c.287C>G affected 1 vs 0 (error +1); p.Glu146Lys c.436G>A unaffected 4 vs 0 (error +4)

### KCNQ1 PMID 29372044

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.A512Pfs*81 c.1532_1534delG affected 1 vs 0 (error +1)

### SCN5A PMID 15851440

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: R1193Q carriers 2 vs 3 (error -1); R1193Q affected 1 vs 2 (error -1)

### RYR2 PMID 21292648

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17205354

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: F1760A

### SCN5A PMID 20038812

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 18464931

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: Y315S

### KCNH2 PMID 14676148

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: N588K

### SCN5A PMID 19843921

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.3142_3143insTG, c.3840+1G>A, c.4719C>T, c.934+1G>A, p.Phe861fs, p.Trp774fs

### SCN5A PMID 11304498

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 28100344

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 11274952

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: V1667I carriers 9 vs 11 (error -2); V1667I unaffected 7 vs 9 (error -2)

### SCN5A PMID 28739862

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: L409P

### KCNQ1 PMID 31899541

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 32830254

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: W392X + I145SFSX
- Count disagreements: p.W392* c.1175G>A affected 1 vs 0 (error +1)

### KCNQ1 PMID 29037160

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R555C, c.477+1G>A

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

- Count disagreements: p.Ala302Val c.905C>T affected 1 vs 0 (error +1); p.Ala46Thr c.136G>A affected 1 vs 0 (error +1); p.Arg195Trp c.583C>T affected 1 vs 0 (error +1); p.Arg670Lys c.2009G>A affected 1 vs 0 (error +1); p.Ala302Val c.905C>T unaffected 0 vs 1 (error -1); p.Ala46Thr c.136G>A unaffected 0 vs 1 (error -1); p.Arg195Trp c.583C>T unaffected 0 vs 1 (error -1); p.Arg670Lys c.2009G>A unaffected 0 vs 1 (error -1)

### APOE PMID 27108409

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26036855

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: GLN1301DEL, c.1141-3C>A, c.3901_3903delCAG, p.His558Arg c.1673A>G
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
- Extra predictions: p.Ala341Glu, p.Thr312Ile

### KCNH2 PMID 21130771

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: N588K
- Count disagreements: p.T618I c.1853C>T carriers 4 vs 0 (error +4)

### SCN5A PMID 26636822

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R1623X, R475S c.1425A>C, c.3890_3891insA, p.His558Arg c.1673A>G

### RYR2 PMID 28158428

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.I4867V

### SCN5A PMID 24112685

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.T290fsX53

### KCNQ1 PMID 26022593

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg243His c.128G>A carriers 4 vs 3 (error +1); p.Arg243His c.128G>A affected 1 vs 0 (error +1)

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

- Extra predictions: A572D, F2004L, P1090L, P2006A, S524Y, p.Arg1193Gln, p.Ser1103Tyr, p.Ser216Leu, p.Val1951Leu

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

- Extra predictions: F2004L, P2006A, p.Ser216Leu

### KCNH2 PMID 26847485

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 23237912

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### BRCA2 PMID 26848529

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.7435+53C>T, c.7976+1G>A, p.Ala634Valfs*18 c.1901delC, p.Ala938Profs*21 c.2808_2811delACAA, p.Arg155Lysfs*26 c.464_468delGAGAT, p.Arg2108Cys c.6322C>T, p.Arg2336His c.7007G>A, p.Arg435Lysfs*17 c.1303dupA, p.Asn1287Ilefs*6 c.3860delA, p.Asn1473Glnfs*8 c.4416_4417delGA, p.Asn1822Ilefs*18 c.5465delA, p.Asn2051* c.6150_6151insT, p.Asn2101_Val2102del c.6301_6306delAATGTA, p.Asn361Metfs*6 c.1082delA, p.Asn863Lysfs*18 c.2588dupA, p.Asn991Asp c.2971A>G, p.Asp687* c.2059_2063delGATTA, p.Asp974Asn c.2920G>A, p.Gln1056Argfs*3 c.3167_3170delAAAA, p.Gln1291* c.3871C>T, p.Gln2655Asnfs*2 c.7963delC, p.Gln2829* c.8485C>T, p.Gln2941Leufs*34 c.8820_8823del, p.Gln3034Serfs*10 c.9098_9099insA, p.Gln3036Serfs*8 c.9105dup, p.Glu2198Asnfs*4 c.6591_6592delTG, p.Glu3096* c.9286G>T, p.Glu3263Argfs*12 c.9788delA, p.Gly2901Asp c.8702G>A, p.Gly3003* c.9007G>T, p.His2021Pro c.6062_6063delinsCA, p.His372Asn c.1114C>A, p.Ile1724Lysfs*17 c.5171delT, p.Ile1859Lysfs*3 c.5576_5579delTTAA, p.Ile2986Lysfs*3 c.8956_8957insAA, p.Ile332Lysfs*18 c.993_994dupA, p.Ile605Tyrfs*9 c.1813delA, p.Ile770Phefs*2 c.2307delT, p.Leu1908Argfs*2 c.5722_5723delCT, p.Leu2080* c.6239T>G, p.Leu3119* c.9356_9357delTAinsG, p.Leu88Alafs*12 c.262_263delCT, p.Lys157Serfs*24 c.470_474delAGTCA, p.Lys2150Asnfs*19 c.6449_6450insTA, p.Lys3326* c.9976A>T, p.Lys467* c.1399A>T, p.Lys585Arg c.1754A>G, p.Met2393Lysfs*18 c.7178_7179delTG, p.Met815Trpfs*10 c.2442delC, p.Phe2801Leufs*10, p.Pro1702Thrfs*16 c.5103dupA, p.Pro2381Hisfs*13 c.7142delC, p.Pro628Hisfs*16 c.1881delA, p.Pro999Leu c.2996C>T, p.Ser1248Argfs*10 c.3744_3747delTGAG, p.Ser1722Tyrfs*4 c.5164_5165delAG, p.Ser1882* c.5645C>A, p.Ser1900* c.5699C>G, p.Ser1955* c.5864C>G, p.Ser2012Glnfs*5 c.6033_6034delTT, p.Ser2120* c.6359C>G, p.Ser2267* c.6800C>A, p.Ser2984Glnfs*4 c.8950delT, p.Ser611* c.1832C>A, p.Thr2471Hisfs*4 c.7409dupT, p.Thr2746Aspfs*18 c.8234dupT, p.Thr2867Ser c.8599A>T, p.Thr3085Glnfs*19 c.9253delA, p.Thr598Ala c.1792A>G, p.Trp2725Glyfs*8 c.8172delG, p.Tyr1894* c.5681dupA, p.Tyr2215*, p.Tyr2839* c.8517C>A, p.Val220Ilefs*4 c.658_659delGT
- Count disagreements: c.5576_5580delTTAAA carriers 1 vs 5 (error -4)

### SCN5A PMID 19083750

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 24057343

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: H302fsX, H492Y

### KCNH2 PMID 27761169

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A78G, A78V
- Count disagreements: A78T affected 1 vs 0 (error +1)

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

- No scored variant or count disagreement.

### KCNQ1 PMID 16627448

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Ala341Val unaffected 1 vs 2 (error -1)

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

- No scored variant or count disagreement.

### SCN5A PMID 27676163

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15161528

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.4245+82A>G

### SCN5A PMID 16039271

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1330P

### SCN5A PMID 27082542

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A29A, E1061E
- Count disagreements: p.R1632C c.4894C>T carriers 2 vs 0 (error +2); p.R1632C c.4894C>T unaffected 1 vs 0 (error +1)

### KCNQ1 PMID 32405922

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: V416M

### KCNQ1 PMID 22250012

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: S140C, V141C, p.Ile145Cys

### KCNQ1 PMID 18400097

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Thr587Met c.1760C>T

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
- Extra predictions: H558R, c.393-1C>T

### KCNH2 PMID 30244407

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: D46F

### SCN5A PMID 29672598

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F2004L
- Count disagreements: p.Arg1193Gln c.3578G>A carriers 1 vs 2 (error -1); p.His558Arg c.1673A>G carriers 1 vs 2 (error -1); p.Phe2004Leu c.6010T>C carriers 2 vs 1 (error +1); p.Arg1193Gln c.3578G>A affected 1 vs 2 (error -1); p.His558Arg c.1673A>G affected 1 vs 2 (error -1); p.Phe2004Leu c.6010T>C affected 2 vs 1 (error +1)

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

- Count disagreements: p.Ser349* affected 1 vs 2 (error -1)

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

- Missed gold variants: P.F1617DEL

### KCNH2 PMID 17171344

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Gly604Ser c.1810G>A carriers 11 vs 10 (error +1)

### SCN5A PMID 20022821

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

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

- Count disagreements: p.Arg1193Gln c.3578G>A carriers 5 vs 1 (error +4)

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
