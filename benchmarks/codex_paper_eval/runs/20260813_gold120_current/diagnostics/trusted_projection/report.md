# Codex extraction-blinded paper evaluation — `20260813_gold120_current`

## Technical summary

This hash-locked run evaluated **120 papers** (**30 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **37.1%**, recall **84.1%**, F1 **51.5%** (534 TP, 906 FP, 101 FN).
- Precision versus counted extras **90.8%** (534 matched rows; 54 extra rows with patient counts). The stricter count-bearing-only diagnostic is **72.7%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **2,363,592 total tokens** (1,732,043 input; 631,549 output).
- Elapsed: **8360.6s wall clock**; 7559.9s summed per-paper route + read time.
- Representation choices: {'text': 120}.

## Blinding and scorer audit

- Paper selection used a seeded sample of 120 papers (30 per cardiac gene; seed 2026081301) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: prediction content was finalized and SHA-256 locked before this external score. The source production workflow may have read registered gold for read-only layer scorecards before the projection lock; those scores did not feed back into extraction. This is not the stricter native lock-before-any-gold-read protocol.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 140 / 635 | 22.0% | 0.243 | 0.878 |
| affected | 60 / 635 | 9.4% | 0.783 | 1.794 |
| unaffected | 44 / 634 | 6.9% | 0.523 | 1.340 |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 185 | 219 | 14 | 45.8% | 93.0% | 61.4% | 20.1% / 0.375 / 1.313 | 3.5% / 0.000 / 0.000 | 2.5% / 0.400 / 0.632 |
| KCNH2 | 65 | 131 | 24 | 33.2% | 73.0% | 45.6% | 41.6% / 0.189 / 0.593 | 29.2% / 1.500 / 2.638 | 11.2% / 1.500 / 2.550 |
| KCNQ1 | 130 | 277 | 37 | 31.9% | 77.8% | 45.3% | 13.8% / 0.304 / 0.808 | 10.8% / 0.444 / 0.816 | 7.2% / 0.500 / 1.000 |
| RYR2 | 154 | 279 | 26 | 35.6% | 85.6% | 50.2% | 22.2% / 0.125 / 0.524 | 5.0% / 0.000 / 0.000 | 9.5% / 0.000 / 0.000 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| RYR2 | 19926015 | text | 38 | 78 | 2 | 32.8% | 95.0% | 48.7% | 2.5% / 0.000 | 2.5% / 0.000 | 2.5% / 0.000 | 124.8 | 40,320 |
| KCNH2 | 18675227 | text | 3 | 9 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 0.333 | 100.0% / 0.667 | 66.7% / 0.500 | 33.0 | 16,591 |
| RYR2 | 33640691 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 27.1 | 11,330 |
| RYR2 | 30403697 | text | 20 | 22 | 1 | 47.6% | 95.2% | 63.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 98.7 | 35,485 |
| RYR2 | 27839804 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 34.2 | 10,694 |
| RYR2 | 15466642 | text | 5 | 4 | 3 | 55.6% | 62.5% | 58.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 6.0 | 8,346 |
| KCNH2 | 10086972 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 22222782 | text | 7 | 1 | 0 | 87.5% | 100.0% | 93.3% | 71.4% / 0.800 | 0.0% / n/a | 0.0% / n/a | 202.7 | 50,954 |
| SCN5A | 21193062 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 21.5 | 17,852 |
| RYR2 | 22178870 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 45.4 | 17,321 |
| KCNH2 | 22314138 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.3 | 14,076 |
| SCN5A | 32533946 | text | 83 | 52 | 0 | 61.5% | 100.0% | 76.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 195.0 | 75,162 |
| KCNQ1 | 29677589 | text | 2 | 2 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 1.000 | 50.0% / 0.000 | 31.0 | 15,454 |
| SCN5A | 28087622 | text | 7 | 10 | 0 | 41.2% | 100.0% | 58.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 78.7 | 26,981 |
| RYR2 | 19398665 | text | 26 | 7 | 1 | 78.8% | 96.3% | 86.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 48.9 | 18,453 |
| RYR2 | 25463374 | text | 1 | 10 | 0 | 9.1% | 100.0% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 47.8 | 27,356 |
| SCN5A | 22677073 | text | 12 | 10 | 0 | 54.5% | 100.0% | 70.6% | 91.7% / 0.091 | 0.0% / n/a | 16.7% / 0.000 | 170.2 | 66,219 |
| KCNQ1 | 20959120 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.2 | 15,052 |
| KCNQ1 | 28491547 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 3.000 | 100.0% / 0.000 | 100.0% / 3.000 | 40.6 | 16,914 |
| RYR2 | 18285261 | text | 5 | 1 | 0 | 83.3% | 100.0% | 90.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 18.2 | 13,964 |
| KCNQ1 | 19632626 | text | 1 | 67 | 0 | 1.5% | 100.0% | 2.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 188.3 | 21,391 |
| KCNQ1 | 22613981 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 12.4 | 5,773 |
| RYR2 | 27114410 | text | 4 | 4 | 5 | 50.0% | 44.4% | 47.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 9.5 | 4,933 |
| KCNH2 | 29016797 | text | 2 | 5 | 0 | 28.6% | 100.0% | 44.4% | 100.0% / 2.000 | 100.0% / 2.000 | 50.0% / 0.000 | 44.3 | 24,161 |
| KCNH2 | 16029385 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 0.000 | 100.0% / 7.000 | 100.0% / 7.000 | 76.3 | 19,913 |
| RYR2 | 16517285 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 50.0% / 0.000 | 50.0% / 0.000 | 50.0% / 0.000 | 25.9 | 9,017 |
| KCNQ1 | 19114714 | text | 3 | 2 | 0 | 60.0% | 100.0% | 75.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.6 | 15,173 |
| KCNQ1 | 16155735 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 18.8 | 15,838 |
| KCNQ1 | 21129503 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 25.1 | 16,516 |
| KCNQ1 | 20348026 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.1 | 19,513 |
| KCNQ1 | 19056345 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.7 | 13,064 |
| RYR2 | 19398417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 23.3 | 6,551 |
| KCNQ1 | 21956039 | text | 5 | 0 | 11 | 100.0% | 31.2% | 47.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 54.7 | 38,764 |
| KCNH2 | 22067087 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 25.4 | 8,725 |
| SCN5A | 20038812 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 17.7 | 18,896 |
| SCN5A | 16301357 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 31.0 | 8,314 |
| RYR2 | 35663620 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 33.3 | 12,638 |
| RYR2 | 17161793 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 27.5 | 10,432 |
| KCNQ1 | 20368164 | text | 2 | 3 | 0 | 40.0% | 100.0% | 57.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.8 | 20,131 |
| KCNQ1 | 18567635 | text | 1 | 10 | 0 | 9.1% | 100.0% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 125.1 | 37,112 |
| KCNH2 | 14642689 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 30036649 | text | 5 | 4 | 0 | 55.6% | 100.0% | 71.4% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 70.3 | 30,626 |
| KCNH2 | 26118593 | text | 1 | 11 | 0 | 8.3% | 100.0% | 15.4% | 0.0% / n/a | 100.0% / 0.000 | 0.0% / n/a | 53.9 | 22,834 |
| SCN5A | 15671429 | text | 4 | 16 | 0 | 20.0% | 100.0% | 33.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 84.6 | 23,366 |
| KCNQ1 | 18808722 | text | 1 | 6 | 0 | 14.3% | 100.0% | 25.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 67.2 | 24,699 |
| SCN5A | 29709101 | text | 10 | 17 | 1 | 37.0% | 90.9% | 52.6% | 27.3% / 0.000 | 0.0% / n/a | 0.0% / n/a | 66.0 | 31,489 |
| SCN5A | 25171853 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 1.000 | 100.0% / 0.000 | 50.0% / 1.000 | 65.1 | 15,091 |
| SCN5A | 15898185 | text | 2 | 0 | 9 | 100.0% | 18.2% | 30.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 21.7 | 7,447 |
| KCNH2 | 22764740 | text | 1 | 14 | 0 | 6.7% | 100.0% | 12.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 24.1 | 21,540 |
| KCNQ1 | 33082985 | text | 2 | 6 | 1 | 25.0% | 66.7% | 36.4% | 66.7% / 1.000 | 66.7% / 1.500 | 66.7% / 0.500 | 54.9 | 29,939 |
| KCNQ1 | 21070882 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 43.8 | 16,845 |
| SCN5A | 17675083 | text | 2 | 4 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 67.5 | 30,050 |
| SCN5A | 20539757 | text | 11 | 15 | 0 | 42.3% | 100.0% | 59.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 105.9 | 23,777 |
| RYR2 | 25435091 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 24.4 | 10,966 |
| RYR2 | 30546600 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 113.0 | 19,628 |
| SCN5A | 19406494 | text | 1 | 9 | 0 | 10.0% | 100.0% | 18.2% | 100.0% / 1.000 | 100.0% / 0.000 | 100.0% / 1.000 | 73.5 | 18,461 |
| KCNH2 | 10086971 | text | 2 | 6 | 1 | 25.0% | 66.7% | 36.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 102.1 | 27,051 |
| KCNQ1 | 28739325 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 2.000 | 0.0% / n/a | 0.0% / n/a | 36.6 | 23,428 |
| RYR2 | 25814417 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 1228.2 | 20,510 |
| KCNQ1 | 18580685 | text | 1 | 7 | 0 | 12.5% | 100.0% | 22.2% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 22.0 | 14,101 |
| SCN5A | 22882672 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 7.000 | 0.0% / n/a | 0.0% / n/a | 33.4 | 9,090 |
| KCNQ1 | 24667783 | text | 14 | 1 | 2 | 93.3% | 87.5% | 90.3% | 68.8% / 0.000 | 68.8% / 0.273 | 37.5% / 0.167 | 157.2 | 66,937 |
| KCNH2 | 19184172 | text | 3 | 5 | 1 | 37.5% | 75.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 383.5 | 39,269 |
| RYR2 | 20395638 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 34.1 | 9,060 |
| KCNH2 | 21308345 | text | 2 | 2 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 21.0 | 14,004 |
| KCNH2 | 21779290 | text | 1 | 5 | 0 | 16.7% | 100.0% | 28.6% | 100.0% / 0.000 | 100.0% / 1.000 | 100.0% / 1.000 | 69.0 | 20,310 |
| RYR2 | 30835254 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.9 | 22,136 |
| KCNH2 | 22338672 | text | 2 | 2 | 0 | 50.0% | 100.0% | 66.7% | 50.0% / 0.000 | 100.0% / 5.000 | 0.0% / n/a | 188.9 | 27,736 |
| RYR2 | 32866913 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 17.8 | 7,759 |
| RYR2 | 34661651 | text | 4 | 8 | 0 | 33.3% | 100.0% | 50.0% | 75.0% / 0.333 | 0.0% / n/a | 0.0% / n/a | 55.6 | 28,555 |
| RYR2 | 33536282 | text | 1 | 45 | 0 | 2.2% | 100.0% | 4.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 53.7 | 31,346 |
| SCN5A | 24581105 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 26.5 | 9,129 |
| KCNH2 | 30246897 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.333 | 33.3% / 0.000 | 51.8 | 24,019 |
| KCNQ1 | 27114410 | text | 1 | 3 | 1 | 25.0% | 50.0% | 33.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.7 | 4,998 |
| KCNH2 | 29121719 | text | 3 | 3 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.0 | 11,742 |
| KCNH2 | 23936059 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 50.0% / 5.000 | 0.0% / n/a | 53.9 | 22,268 |
| SCN5A | 22685113 | text | 10 | 7 | 0 | 58.8% | 100.0% | 74.1% | 100.0% / 0.000 | 10.0% / 0.000 | 0.0% / n/a | 167.4 | 38,970 |
| KCNQ1 | 26496715 | text | 38 | 131 | 4 | 22.5% | 90.5% | 36.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 107.0 | 33,873 |
| RYR2 | 21652165 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 23.8 | 9,905 |
| SCN5A | 19549036 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.1 | 12,944 |
| KCNQ1 | 11351021 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.9 | 14,285 |
| KCNH2 | 26746457 | text | 10 | 34 | 0 | 22.7% | 100.0% | 37.0% | 100.0% / 0.100 | 0.0% / n/a | 0.0% / n/a | 2.1 | 1,063 |
| KCNH2 | 20181576 | text | 2 | 3 | 0 | 40.0% | 100.0% | 57.1% | 50.0% / 0.000 | 50.0% / 1.000 | 50.0% / 2.000 | 78.7 | 23,549 |
| KCNQ1 | 29851656 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 24.0 | 9,683 |
| SCN5A | 22101522 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.5 | 17,962 |
| KCNH2 | 25819988 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 0.000 | 100.0% / 3.000 | 100.0% / 3.000 | 23.0 | 15,009 |
| KCNH2 | 29650123 | text | 2 | 0 | 20 | 100.0% | 9.1% | 16.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.3 | 18,403 |
| KCNH2 | 22052944 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 18.8 | 16,369 |
| KCNQ1 | 29876285 | text | 1 | 7 | 0 | 12.5% | 100.0% | 22.2% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 1.000 | 58.5 | 23,239 |
| SCN5A | 21288276 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 34.5 | 11,138 |
| KCNH2 | 23917959 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 1.000 | 100.0% / 0.000 | 100.0% / 1.000 | 26.8 | 8,984 |
| SCN5A | 22966897 | text | 0 | 4 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 21.8 | 12,583 |
| KCNQ1 | 14678125 | text | 35 | 9 | 6 | 79.5% | 85.4% | 82.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 4.7 | 7,227 |
| KCNQ1 | 31520628 | text | 2 | 2 | 7 | 50.0% | 22.2% | 30.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.7 | 12,241 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 26.4 | 16,170 |
| RYR2 | 22677073 | text | 24 | 9 | 2 | 72.7% | 92.3% | 81.4% | 76.9% / 0.000 | 0.0% / n/a | 38.5% / 0.000 | 273.2 | 94,548 |
| RYR2 | 25500949 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 27.6 | 13,443 |
| KCNQ1 | 31293497 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 38.3 | 7,128 |
| RYR2 | 29350269 | text | 3 | 69 | 0 | 4.2% | 100.0% | 8.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 40.7 | 28,568 |
| SCN5A | 24573164 | text | 20 | 10 | 0 | 66.7% | 100.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 85.9 | 23,055 |
| KCNH2 | 19065538 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 16.2 | 6,398 |
| SCN5A | 16929919 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 21.9 | 6,387 |
| KCNQ1 | 32168391 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 45.2 | 17,936 |
| KCNQ1 | 25471708 | text | 4 | 5 | 4 | 44.4% | 50.0% | 47.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 48.1 | 17,831 |
| SCN5A | 26820365 | text | 1 | 16 | 0 | 5.9% | 100.0% | 11.1% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 48.9 | 33,099 |
| RYR2 | 19216760 | text | 1 | 0 | 4 | 100.0% | 20.0% | 33.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 20.6 | 6,755 |
| KCNH2 | 16155735 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 50.0% / 0.000 | 50.0% / 0.000 | 0.0% / n/a | 28.4 | 16,261 |
| RYR2 | 28798025 | text | 0 | 0 | 6 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 2.9 | 1,064 |
| KCNH2 | 15500450 | text | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 22.3 | 10,677 |
| SCN5A | 23538678 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 44.9 | 12,465 |
| KCNH2 | 29214556 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.3 | 7,422 |
| SCN5A | 29544605 | text | 6 | 33 | 0 | 15.4% | 100.0% | 26.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 66.4 | 41,239 |
| RYR2 | 33315912 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 30.0 | 20,574 |
| KCNQ1 | 19687231 | text | 2 | 4 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 35.6 | 19,251 |
| SCN5A | 17971661 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 9.0 | 4,453 |
| KCNH2 | 19034806 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 30.6 | 12,973 |
| SCN5A | 19305639 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 13.5 | 9,979 |
| SCN5A | 15828879 | text | 1 | 6 | 0 | 14.3% | 100.0% | 25.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 37.2 | 12,160 |
| SCN5A | 25236808 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 100.0% / 0.000 | 0.0% / n/a | 43.7 | 11,780 |
| KCNH2 | 18808722 | text | 8 | 11 | 0 | 42.1% | 100.0% | 59.3% | 25.0% / 0.000 | 12.5% / 5.000 | 0.0% / n/a | 138.6 | 38,964 |

## Errors and representation choices

### RYR2 PMID 19926015

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: M81L, EXON 3 DELETION
- Extra predictions: A4282V, G1885E, G1886S, M2389L, Q2958R, R169Q, V2321M, V4010M, c.(?_169, c.169-198_273+820del, p.Ala2213Val c.6638C>T, p.Ala2387Thr c.7159G>A, p.Ala4510Ser c.13528G>T, p.Ala4510Thr c.13528G>A, p.Ala77Val c.230C>T, p.Arg1051His c.3152G>A, p.Arg169Gly c.505C>G, p.Arg176Gln c.527G>A, p.Arg176Leu c.527G>T, p.Arg2401Cys c.7201C>T, p.Arg2401His c.7202G>A, p.Arg332Gln c.995G>A, p.Arg3570Trp c.10708C>T, p.Arg414Cys c.1240C>T, p.Arg420Trp c.1258C>T, p.Arg4307Cys c.12919C>T, p.Arg4822Gly c.14464C>G, p.Arg4959Gln c.14876G>A, p.Arg838Leu c.2513G>T, p.Asn3308Ser c.9923A>G, p.Asp2431Tyr c.7291G>T, p.Asp4001Asn c.12001G>A, p.Gln3811Pro c.11432A>C, p.Gln4936Arg c.14807A>G, p.Glu1724Lys c.5170G>A, p.Glu4431Lys c.13291G>A, p.Gly172Glu c.515G>A, p.Gly2949Val c.8846G>T, p.Gly357Val c.1070G>T, p.Gly3946Ser c.11836G>A, p.Gly4662Ser c.13984G>A, p.Gly4864Cys c.14590G>T, p.Gly4935Arg c.14803G>A, p.Gly809Glu c.2426G>A, p.His202Gln c.606C>A, p.His240Arg c.719A>G, p.Ile2009Met c.6027T>G, p.Ile419Ser c.1256T>G, p.Ile4857Val c.14569A>G, p.Leu1459Ser c.4376T>C, p.Leu2123Phe c.6367C>T, p.Leu4105Phe c.12313C>T, p.Leu4698Pro c.14093T>C, p.Lys3717Arg c.11150A>G, p.Met4109Thr c.12326T>C, p.Phe2307Leu c.6921C>G, p.Pro1395Ala c.4183C>G, p.Pro2328Leu c.6983C>T, p.Pro466Arg c.1397C>G, p.Pro990Gln c.2969C>A, p.Ser221Gly c.661A>G, p.Ser2246Leu c.6737C>T, p.Ser2312Gly c.6934A>G, p.Ser2393Pro c.7177T>C, p.Ser2653Tyr c.7958C>A, p.Ser4155Thr c.12463T>A, p.Ser4667Ile c.14000G>T, p.Thr1093Pro c.3277A>C, p.Thr2510Ala c.7528A>G, p.Thr415Ile c.1244C>T, p.Thr4755Ile c.14264C>T, p.Tyr2392Cys c.7175A>G, p.Tyr4149Cys c.12446A>G, p.Val2113Met c.6337G>A, p.Val2306Ile c.6916G>A, p.Val4771Ile c.14311G>A, p.Val4771Phe c.14311G>T, p.Val4778Leu c.14332G>T

### KCNH2 PMID 18675227

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1116V, A341V, G1034D, K897T, L955C, L995V, Q725X, R1014X, T1062I
- Count disagreements: p.Leu955Val c.2863C>G carriers 1 vs 2 (error -1); p.Arg954Cys c.2860C>T affected 1 vs 0 (error +1); p.Leu955Val c.2863C>G affected 1 vs 0 (error +1); p.Arg954Cys c.2860C>T unaffected 1 vs 2 (error -1)

### RYR2 PMID 33640691

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 30403697

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: G4722S
- Extra predictions: A103V, C1274Y, D4642H, D85N, H1042R, I882V, IVS5+1G>C, K296Q, K4742R, L714R, M0781T, Q584R, Q692K, R1458G, R4954K, S2774G, V288I, c.13781A>G, c.6224T>C, p.G4772S, p.Arg2401Leu c.7202G>T, p.Arg2474Gly c.7420A>G

### RYR2 PMID 27839804

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: H29D

### RYR2 PMID 15466642

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P164S, F4499C, G4671R
- Extra predictions: p.Ala2403Val c.7208C>T, p.Ala4510Val c.13529C>T, p.Arg414Cys c.1240C>T, p.Arg414His c.1241G>A

### KCNH2 PMID 10086972

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.2592+1G->A

### RYR2 PMID 22222782

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Leu4767= c.14301C>T
- Count disagreements: p.D1872N carriers 1 vs 2 (error -1); p.H4579Y carriers 1 vs 4 (error -3)

### SCN5A PMID 21193062

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 22178870

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 22314138

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: R744X, R774P, p.Arg744Gln c.2231G>A, p.Arg744Gly c.2230C>G

### SCN5A PMID 32533946

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1632H, D1430N, D356N, G1406R, G1408R, G1712C, G1740R, G1743E, G1743R, G752A, G897E, I1660V, L846R, M1766L, P.LEU839P, R104Q, R104W, R282H, R878C, R878H, R893H, S1218I, S910L, T187I, p.Ala1679Thr c.5035G>A, p.Ala1745Val c.5234C>T, p.Arg1582Cys c.4744C>T, p.Arg1631His c.4892G>A, p.Arg1631Ser c.4891C>A, p.Arg1897Cys c.5689C>T, p.Arg1957Ter c.5869C>T, p.Asn1721Asp c.5161A>G, p.Asp1242Asn c.3724G>A, p.Asp1369Gly c.4106A>G, p.Asp349Gly c.1046A>G, p.Glu1224Lys c.3670G>A, p.Glu1573Lys c.4717G>A, p.Glu1783Lys c.5347G>A, p.Gly1261Ser c.3781G>A, p.Gly1405Glu c.4214G>A, p.Gly1419Ala c.4256G>C, p.Gly1419Cys c.4255G>T, p.Gly1641Glu c.4922G>A, p.Gly1660Arg c.4978G>A, p.Gly1660Glu c.4979G>A, p.Gly1739Arg c.5215G>A, p.Thr1708Lys c.5123C>A, p.Thr1708Met c.5123C>T, p.Tyr1448Cys c.4343A>G, p.Val1250Met c.3748G>A, p.Val1278Ile c.3832G>A, p.Val1352Met c.4054G>A

### KCNQ1 PMID 29677589

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.1781G>A, c.568C>T
- Count disagreements: p.R190W c.568 C>T affected 1 vs 0 (error +1); p.R594Q c.1781 G>A affected 1 vs 0 (error +1)

### SCN5A PMID 28087622

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: C373Y, D1840G, H1853R, R1918H, Y1795H, Y371C, Y373C, p.Arg1912His c.5735G>A, p.Gln1908Arg c.5723A>G, p.Glu1900Gln c.5698G>C

### RYR2 PMID 19398665

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: N4101I
- Extra predictions: G3946S, N4104I, p.Ala2387Val c.7160C>T, p.Ala4091Arg c.12271_12272delinsAG, p.Ala4091Gly c.12272C>G, p.Asn4178Ser c.12533A>G, p.Gly3946Arg c.11836G>C

### RYR2 PMID 25463374

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A77V, G230C, P2328S, Q4201R, R2474S, R3570W, R4497C, S2226L, S4565R, V4653F

### SCN5A PMID 22677073

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1136V, C108R, G38S, R1047L, R420W, S1787N, W379R, p.Ile1767Val c.5299A>G, p.Pro2005Ala c.6013C>G, p.Val1950Leu c.5848G>T
- Count disagreements: S216L carriers 2 vs 1 (error +1)

### KCNQ1 PMID 20959120

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 28491547

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Val141Met c.421G>A carriers 3 vs 6 (error -3); p.Val141Met c.421G>A unaffected 0 vs 3 (error -3)

### RYR2 PMID 18285261

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: S453S

### KCNQ1 PMID 19632626

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: S140G, S140P, S209F, V141M, c.1033-2del, c.1251+9C>A, c.1394-1G>T, c.1515-4G>A, c.477+5G>A, p.Ala178Thr c.532G>A, p.Ala341Glu c.1022C>A, p.Ala341Val c.1022C>T, p.Ala344= c.1032G>A, p.Ala344Val c.1031C>T, p.Arg190Gln c.569G>A, p.Arg190Leu c.569G>T, p.Arg190Trp c.568C>T, p.Arg243Cys c.727C>T, p.Arg243His c.728G>A, p.Arg380Lys c.1139G>A, p.Arg401fs c.1201dup, p.Arg507Gln c.1520G>A, p.Arg518Ter c.1552C>T, p.Arg539Trp c.1615C>T, p.Arg555Ser c.1663C>A, p.Arg583His c.1748G>A, p.Arg591His c.1772G>A, p.Arg594Gln c.1781G>A, p.Arg632fs c.1893dup, p.Asp202Asn c.604G>A, p.Cys214Tyr c.641G>A, p.Cys642Phe c.1925G>T, p.Gln107His c.321G>C, p.Gln359Ter c.1075C>T, p.Gln376Ter c.1126C>T, p.Gln530Ter c.1588C>T, p.Glu449fs c.1343dup, p.Gly168Arg c.502G>A, p.Gly269Ser c.805G>A, p.Gly306Arg c.916G>A, p.Gly314Ser c.940G>A, p.Gly325Arg c.973G>A, p.Gly626fs c.1875dup, p.His258Arg c.773A>G, p.Ile375fs c.1124_1127del, p.Leu266Pro c.797T>C, p.Leu266fs c.796del, p.Leu496fs c.1486_1487del, p.Lys362Arg c.1085A>G, p.Lys526Glu c.1576A>G, p.Met520Arg c.1559T>G, p.Pro343Ala c.1027C>G, p.Ser277del c.825CTC[1], p.Ser546Leu c.1637C>T, p.Ser566Phe c.1697C>T, p.Thr311Ile c.932C>T, p.Thr312Ile c.935C>T, p.Thr322Met c.965C>T, p.Trp15Ter c.45G>A, p.Trp305Leu c.914G>T, p.Trp305Ter c.914G>A, p.Tyr111Cys c.332A>G, p.Tyr184Ser c.551A>C, p.Val110Ile c.328G>A, p.Val215Met c.643G>A, p.Val254Met c.760G>A, p.Val280Gly c.839T>G

### KCNQ1 PMID 22613981

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Arg231Leu c.692G>T

### RYR2 PMID 27114410

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: T4158P, L3879P, Q3925E, S3959L, G4936R
- Extra predictions: F90L, N634FS, N98S, p.Gly4935Arg c.14803G>A

### KCNH2 PMID 29016797

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: S755T, c.2264G>C, c.5010G>T, c.5605C>T, c.631C>T
- Count disagreements: p.Asn588Lys c.1764C>G carriers 1 vs 0 (error +1); p.Ser631Ala c.1891T>G carriers 3 vs 0 (error +3); p.Asn588Lys c.1764C>G affected 1 vs 0 (error +1); p.Ser631Ala c.1891T>G affected 3 vs 0 (error +3)

### KCNH2 PMID 16029385

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G601S, T65P, p.Val822Met c.2464G>A
- Count disagreements: M124R affected 16 vs 9 (error +7); M124R unaffected 0 vs 7 (error -7)

### RYR2 PMID 16517285

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R169Q

### KCNQ1 PMID 19114714

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Pro117Thr c.349C>A, p.Tyr111Ser c.332A>C

### KCNQ1 PMID 16155735

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Pro448Arg c.1343C>G

### KCNQ1 PMID 21129503

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A341V, p.Tyr111Ser c.332A>C

### KCNQ1 PMID 20348026

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G601S, G628S, M124T, N598Q

### KCNQ1 PMID 19056345

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 19398417

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 21956039

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: I235N, G168R, A46T, A344A, S349X, P400FS/62X, R594Q, S546L, R360M, Q530X, G635R

### KCNH2 PMID 22067087

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 20038812

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: I1461T, I1477T, R1448P, p.Phe1485Leu c.4453T>C

### SCN5A PMID 16301357

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Leu1825Pro carriers 1 vs 0 (error +1)

### RYR2 PMID 35663620

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1357V

### RYR2 PMID 17161793

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 20368164

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A340I, p.Ala341Glu, p.Thr312Ile

### KCNQ1 PMID 18567635

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: D242A, E261A, K218A, P475D, R228A, R231A, R237A, R243A, R249A, R259A

### KCNH2 PMID 14642689

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A561T

### KCNH2 PMID 30036649

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G314F, G314S, G628D, I239V

### KCNH2 PMID 26118593

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A29A, E1061E, F29L, G643S, H558R, I145I, S546S, T309I, p.F513F c.1539C>T, p.I489I c.1467C>T, p.L564L c.1692A>G

### SCN5A PMID 15671429

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: C659T, D1595N, R225W, S1103Y, V1397X, W156X, c.2551insTG, c.2550_2551insTG, c.704-2A>G, p.Arg814Gln c.2441G>A, p.Asp1274Asn c.3820G>A, p.Asp1594Asn c.4780G>A, p.Asp1594His c.4780G>C, p.Gly1572= c.4716C>T, p.Phe851fs c.2550_2551dup, p.Pro1331fs c.3992del

### KCNQ1 PMID 18808722

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A490P, A490T, G189R, R190Q, c.2020insAG, p.Leu187Arg c.560T>G

### SCN5A PMID 29709101

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.4813+3_4813_6dup
- Extra predictions: F756L, F756fsX, I1660V, c.1045G>A, c.1100G>A, c.1127G>A, c.1312A>T, c.2465G>A, c.2466G>A, c.2658T>A, c.4895G>A, c.4978A>G, c.4813+3_4813+6dup, p.Arg1631His c.4892G>A, p.Asp349Gly c.1046A>G, p.Ile1659Val c.4975A>G, p.Met438Thr c.1313T>C

### SCN5A PMID 25171853

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.T220I c.659C>T carriers 2 vs 1 (error +1); p.T220I c.659C>T unaffected 1 vs 0 (error +1)

### SCN5A PMID 15898185

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: F1293S, H558R, P1090L, R1512W, R34C, R481W, S1787N, S524Y, V1951L

### KCNH2 PMID 22764740

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A490T, A915V, G572R, G601S, I593R, K897T, N470D, P596R, R1047L, T474I, T613M, V612L, Y611H, p.Gly487Ser c.1459G>A

### KCNQ1 PMID 33082985

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: W248F + L347R
- Extra predictions: L347P, T372M, W248C, W248R, c.1040T>C, c.1115C>T
- Count disagreements: p.Leu347Arg c.1040T>G carriers 2 vs 1 (error +1); p.Trp248Phe c.743_744delGGinsTC carriers 2 vs 1 (error +1); p.Leu347Arg c.1040T>G affected 2 vs 0 (error +2); p.Trp248Phe c.743_744delGGinsTC affected 1 vs 0 (error +1); p.Leu347Arg c.1040T>G unaffected 0 vs 1 (error -1)

### KCNQ1 PMID 21070882

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: R67H, R98W

### SCN5A PMID 17675083

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: M1766L, R282H, T512I, p.Arg34Cys
- Count disagreements: H558R carriers 5 vs 1 (error +4)

### SCN5A PMID 20539757

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: F1617del, H558R, K1578FS, R1632X, T1620M, c.475_482+3del, p.Arg1622Ter c.4864C>T, p.Arg1631His c.4892G>A, p.Arg1631Ser c.4891C>A, p.Asp1274Asn c.3820G>A, p.Gly1407Arg c.4219G>A, p.Phe1616del c.4844TCT[1], p.Pro1297Leu c.3890C>T, p.Thr187Ala c.559A>G, p.Thr187Ser c.560C>G

### RYR2 PMID 25435091

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: R420Q

### RYR2 PMID 30546600

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.1259G>A, p.Arg56Gln c.167G>A

### SCN5A PMID 19406494

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Ala1949Pro c.5845G>C, p.Trp301* c.903G>A, p.Arg2011Cys c.6031C>T, p.Arg808His c.2423G>A, p.Glu1866Ter c.5596G>T, p.Leu1785fs c.5353_5354del, p.Phe2003Leu c.6007T>C, p.Pro1840fs c.5517dup, p.Trp1712Ter c.5135G>A
- Count disagreements: p.Arg808Cys c.2422C>T carriers 4 vs 3 (error +1); p.Arg808Cys c.2422C>T unaffected 1 vs 0 (error +1)

### KCNH2 PMID 10086971

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: C1117fsX
- Extra predictions: c.3108insG, c.2592+1G>A, p.Asp1037fs c.3107dup, p.Gly965fs c.2892dup, p.Leu953fs c.2857dup, p.Ser818Ala c.2452T>G

### KCNQ1 PMID 28739325

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: D242A, D242W
- Count disagreements: p.D242N c.724G>A carriers 2 vs 0 (error +2)

### RYR2 PMID 25814417

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: V2475F, p.Gly357Asp c.1070G>A

### KCNQ1 PMID 18580685

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: D202H, I204M, Q357R, R294H, R555C, S209F, V215M

### SCN5A PMID 22882672

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Pro1177Leu carriers 1 vs 8 (error -7)

### KCNQ1 PMID 24667783

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: I328_S330DEL, Y522X
- Extra predictions: p.Thr224Met c.671C>T
- Count disagreements: p.Ala300Thr c.898G>A affected 2 vs 0 (error +2); p.Gly168Arg c.502G>A affected 2 vs 1 (error +1); p.Gly168Arg c.502G>A unaffected 0 vs 1 (error -1)

### KCNH2 PMID 19184172

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: Q1070X
- Extra predictions: G537V c.1610G>T, K897T, p.Arg1068Ter c.3202C>T, p.Q1070* c.3208C>T, p.Glu788Asp c.2364G>C

### RYR2 PMID 20395638

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Gly4935Arg c.14803G>A

### KCNH2 PMID 21308345

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.1600G>T, p.Arg534Leu c.1601G>T

### KCNH2 PMID 21779290

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A29A, H558R, p.L564L, p.Y652Y, p.Tyr652= c.1956T>C
- Count disagreements: p.K897T affected 1 vs 0 (error +1); p.K897T unaffected 0 vs 1 (error -1)

### RYR2 PMID 30835254

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1107M, G357S, R403Q, T1107M

### KCNH2 PMID 22338672

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1116V, E446K
- Count disagreements: K897T affected 7 vs 0 (error +7); R744X affected 4 vs 1 (error +3)

### RYR2 PMID 32866913

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R85I
- Extra predictions: p.Thr85Ile c.254C>T

### RYR2 PMID 34661651

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A4860G, I4855M, K4750Q, R169Q, R2267H, R4496C, S4565R, c.14151+1G>A
- Count disagreements: p.S4168P c.12502T>C carriers 1 vs 2 (error -1)

### RYR2 PMID 33536282

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: L3995V, M1045V, Q3776L, R4594K, S11G, S88X, c.11321A>T, c.11983A>C, p.A4860G, p.D4112N, p.D4122N, p.D4646A, p.I2074T, p.I2075T, p.I3995V, p.I4855M, p.K4594R, p.Q4879H, p.Q4897H, p.R4496C, p.S4938F, p.T4196I, c.848+1G>A, p.Ala2387Thr c.7159G>A, p.Ala2403Val c.7208C>T, p.Ala2633Thr c.7897G>A, p.Arg1013Gln c.3038G>A, p.Arg1089Cys c.3265C>T, p.Arg1671Gln c.5012G>A, p.Arg2144Cys c.6430C>T, p.Arg76Trp c.226C>T, p.Asp4646Glu c.13938C>A, p.Glu3716Lys c.11146G>A, p.Glu4076Lys c.12226G>A, p.Glu4182Gln c.12544G>C, p.Gly2284Asp c.6851G>A, p.His4579Tyr c.13735C>T, p.His4908Asn c.14722C>A, p.Lys304Glu c.910A>G, p.Pro2328Ser c.6982C>T, p.Pro916Leu c.2747C>T, p.Ser4153Arg c.12457A>C, p.Thr2188Ala c.6562A>G, p.Tyr2392Cys c.7175A>G, p.Tyr4725Cys c.14174A>G

### SCN5A PMID 24581105

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 30246897

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Cys566Arg c.1696T>C affected 1 vs 0 (error +1)

### KCNQ1 PMID 27114410

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: G584S
- Extra predictions: F90L, N634FS, N98S

### KCNH2 PMID 29121719

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: H558R, R1023C, R659W

### KCNH2 PMID 23936059

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Arg1047Leu c.3140G>T affected 5 vs 0 (error +5)

### SCN5A PMID 22685113

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Arg1625His c.4874G>A, p.Arg1896Trp c.5686C>T, p.Asp1818Asn c.5452G>A, p.Phe1595Ile c.4783T>A, p.Phe2003Leu c.6007T>C, p.Thr1303Met c.3908C>T, p.Val1950Met c.5848G>A

### KCNQ1 PMID 26496715

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A344A, 360_361DUPKQ, S373fsX, R562S
- Extra predictions: A169G, A32T, A384F, A561T, A561V, A614V, D102H, D609G, E480V, E544A, E575K, G569A, G572S, G601S, G604S, G626S, G628S, G71W, H105T, L413P, L457P, L559H, N629S, N633S, P.ALA344S, P.PHE494DEL, P191R, P596S, P72T, Q247X, R366X, R534C, S373W, S620I, T443P, T443R, T613M, T634P, V411M, V440A, W154X, Y475C, c.1020T>G, c.1022C>T, c.1032G>A, c.1032_1117dup, c.1033G>C, c.1034G>C, c.1040T>C, c.1071_1076dupGAAGCA, c.1096C>T, c.1138A>G, c.1140G>T, c.1208A>C, c.1231G>A, c.1238T>C, c.1251+2T>C, c.1326delT, c.1328_1329delC, c.1424A>G, c.1439A>T, c.1480_1482delT, c.1552C>T, c.1556G>A, c.1631A>C, c.1676T>A, c.1681G>A, c.1686G>T, c.1714G>A, c.1723G>A, c.1748G>A, c.1760C>T, c.1765G>A, c.1781G>A, c.1786C>T, c.1801G>A, c.1810G>A, c.1826A>G, c.1859G>T, c.1876G>A, c.1886A>G, c.1898A>G, c.1900A>C, c.214C>A, c.2345T>A, c.2389G>A, c.2414T>C, c.2453C>G, c.2453C>T, c.2467C>T, c.2509G>T, c.2587C>T, c.2616delC, c.2682_2685delC, c.2734C>T, c.2775_2776insG, c.2776C>T, c.2944_2948dupG, c.2T>C, c.3048delC, c.304G>C, c.3094_3095insC, c.3102_3111delC, c.3416G>A, c.420C>A, c.461G>A, c.477+5G>A, c.502G>C, c.506C>G, c.551A>G, c.572delC, c.574C>T, c.605-2G>A, c.674C>T, c.727C>T, c.739C>T, c.752T>A, c.760G>A, c.775C>T, c.815G>A, c.830C>T, c.905C>T, c.914G>T, c.917G>T, c.921+1G>T, c.940G>A, c.94G>A, c.958C>T, c.965C>T, c.973G>A, p.Gln359_Arg c.1067AGCAGA[3]

### RYR2 PMID 21652165

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: R4496C, p.Ser4153Ile c.12458G>T, p.Ser4153Thr c.12458G>C

### SCN5A PMID 19549036

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.R1193Q, p.S1103Y

### KCNQ1 PMID 11351021

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 26746457

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A175D, A188V, A190T, A9V, D982N, G238S, G873S, G903R, K897T, L1105S, M291R, P1018A, P1030L, P347S, P967L, R1035Q, R1047L, R148W, R181Q, R252W, R326H, R328C, R394H, R885C, R894H, R912W, R948S, S981G, V1063L, V325M, W927L, Y403C, p.Arg744Gly c.2230C>G, p.Gly924Arg c.2770G>A
- Count disagreements: A913V carriers 1 vs 2 (error -1)

### KCNH2 PMID 20181576

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: H558R, K897K, P926A
- Count disagreements: K897T affected 1 vs 0 (error +1); K897T unaffected 1 vs 3 (error -2)

### KCNQ1 PMID 29851656

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 22101522

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 25819988

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: H562P, p.His562Asn c.1684C>A, p.His562Tyr c.1684C>T
- Count disagreements: p.H562R c.1685A>G affected 6 vs 9 (error -3); p.H562R c.1685A>G unaffected 7 vs 4 (error +3)

### KCNH2 PMID 29650123

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A561V, C49Y, F617fsX, G572S, G911fsX, L109fsX, L779P, N633S, Q1046X, R1035fsX, R148W, R176W, R328C, R534C, R744X, R892fsX, S660L, S818P, W412X, Y43C

### KCNH2 PMID 22052944

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A614V

### KCNQ1 PMID 29876285

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A300T, A302V, T311I, V319L, W305L, W305S, W305X
- Count disagreements: A300S unaffected 3 vs 4 (error -1)

### SCN5A PMID 21288276

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.999-424_1338+81del
- Extra predictions: F92F

### KCNH2 PMID 23917959

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Gly603Arg c.1807G>C
- Count disagreements: p.Gly603Asp c.1808G>A carriers 1 vs 2 (error -1); p.Gly603Asp c.1808G>A unaffected 0 vs 1 (error -1)

### SCN5A PMID 22966897

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P1177L, R190W
- Extra predictions: R176W, p.Ile1768Val, p.Pro2006Ala, p.Ser1103Tyr

### KCNQ1 PMID 14678125

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: L151, G185S, DELF340, A344/SP, A344A/SPLICE, R518X
- Extra predictions: p.Arg243His c.728G>A, p.Arg243Ser c.727C>A, p.Arg366Leu c.1097G>T, p.Arg366Pro c.1097G>C, p.Asp317Ala c.950A>C, p.Glu284del c.850_852del, p.Gly186Cys c.556G>T, p.Gly186Ser c.556G>A, p.Leu266Arg c.797T>G

### KCNQ1 PMID 31520628

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: G589D, R366W, Y184S, Y315N, V254M, I235N, R252fsX
- Extra predictions: A341V, p.Arg591His c.1772G>A

### SCN5A PMID 28339995

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 22677073

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: V2113M, Q2958R
- Extra predictions: A302V, D85N, F2004L, G643S, P2006A, R1047L, V1951L, W379R, c.2398+5G>T

### RYR2 PMID 25500949

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1744S, A189T, D4301N

### KCNQ1 PMID 31293497

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R360_Q361DUPQKQR

### RYR2 PMID 29350269

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A3039P, A735E, D111Y, D548E, E2013K, G550R, G749E, I2869V, I397T, L22P, N2045S, N236S, N3982D, P736L, P864L, Q2432H, R1250H, R1278Q, R185H, R2639Q, R272C, R420C, R519C, R66W, R985W, S160T, T116M, T137P, T240M, V412M, Y62N, c.11441C>T, c.1190T>C, c.11944A>G, c.1234G>A, c.12490C>G, c.1258C>T, c.1297+5G>T, c.1555G>A, c.1644T>G, c.1648G>A, c.184A>T, c.19505T>C, c.196G>A, c.1978+3C>A, c.2204G>T, c.2207G>A, c.2246G>A, c.2591C>T, c.2928+5C>T, c.2953G>A, c.331C>A, c.3380A>G, c.347G>A, c.3749C>T, c.3833G>A, c.409A>C, c.478A>T, c.554C>T, c.6037G>A, c.6134A>G, c.65A>G, c.707A>G, c.719C>T, c.7296C>A, c.7916G>A, c.814C>T, c.8605A>G, c.9115C>G

### SCN5A PMID 24573164

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: C373Y, p.L567Q, p.R223W, p.Ala1112Val c.3335C>T, p.Arg965Cys c.2893C>T, p.Asp1274Asn c.3820G>A, p.Glu1937Lys c.5809G>A, p.Gly1318Val c.3953G>T, p.Leu1500Val c.4498C>G, p.Ser1139Asn c.3416G>A

### KCNH2 PMID 19065538

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 16929919

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 32168391

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 25471708

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R518X, M159X, Q530X, R192CFS91X
- Extra predictions: c.477+1G>A, c.572_576del, p.Q530* c.1588C>T, p.R518* c.1522C>T, p.Arg192fs c.573_577del

### SCN5A PMID 26820365

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A432T, D561A, G582S, G737R, G844D, I376K, I376T, P1204L, Q854R, R252H, R442C c.1324C>T, R892C, R892H c.2675G>A, T1127C, T286T, c.1682A>C

### RYR2 PMID 19216760

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: S616L, R1051P, EXON 3 DELETION, EXON 3 DELETION

### KCNH2 PMID 16155735

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Asn633Ser c.1898A>G

### RYR2 PMID 28798025

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A4556P, V3557A, R2303G, L3061V, S453T, G1885E

### KCNH2 PMID 15500450

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: D611Y, p.Asp609Ala c.1826A>C, p.Asp609Asn c.1825G>A

### SCN5A PMID 23538678

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.693delCA

### KCNH2 PMID 29214556

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 29544605

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: E42K, G460S, K416FS*0, P.F33LFS, P.K60RFS, R14DEL, R547X, S287R, T372M, c.1115C>T, c.1244_1245insT, c.124G>A, c.13774G>T, c.1378G>A, c.1560+1G>A, c.1639C>T, c.2148+1G>A, c.26011C>T, c.325-2A>G, c.36_38delAAG, c.37432_37433insGTGGTTACTACAGCCTC, c.37437_37438insG, c.415-1G>T, c.47653delA, c.523+2T>C, c.546delT, c.650-2A>G, c.6540delG, c.817C>T, c.843+1G>T, c.861C>G, p.Arg1027Gln c.3080G>A, p.Val1623Ile c.4867G>A

### RYR2 PMID 33315912

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G357S, c.6982C>T

### KCNQ1 PMID 19687231

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Ile328Leu, p.Lys326Arg, p.Thr327Val, p.Val324Leu

### SCN5A PMID 17971661

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: E1784K

### KCNH2 PMID 19034806

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: F99F

### SCN5A PMID 19305639

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15828879

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: R1232W, R1432G, T1620M, p.Arg282Cys c.844C>T, p.Arg282Gly c.844C>G, p.Asn406Ser c.1217A>G

### SCN5A PMID 25236808

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: R518X, p.Phe1616del c.4844TCT[1]

### KCNH2 PMID 18808722

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1116V, A490P, G189R, L187P, R190Q, p.Ala490Val c.1469C>T, p.Ala614Ser c.1840G>T, p.Asn629Ile c.1886A>T, p.Asn629Thr c.1886A>C, p.Asp609Glu c.1827C>A, p.Asp609Gly c.1826A>G
- Count disagreements: A490T affected 7 vs 2 (error +5)

## Scope, method, and limitations

- Population: seeded sample of 120 papers (30 per cardiac gene; seed 2026081301); 30 per cardiac gene; every PMID has downloaded source and at least one gold assertion in each count field.
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
