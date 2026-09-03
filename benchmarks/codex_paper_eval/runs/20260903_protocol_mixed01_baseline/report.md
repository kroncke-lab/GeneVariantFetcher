# Codex extraction-blinded paper evaluation — `20260903_protocol_mixed01_baseline`

## Technical summary

This hash-locked run evaluated **49 papers** (**per-gene counts {'SCN5A': 23, 'KCNH2': 9, 'KCNQ1': 10, 'RYR2': 6, 'MYBPC3': 1}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **71.8%**, recall **64.0%**, F1 **67.7%** (155 TP, 61 FP, 87 FN).
- Precision versus counted extras **82.0%** (155 matched rows; 34 extra rows with patient counts). The stricter count-bearing-only diagnostic is **78.2%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **1,078,934 total tokens** (792,800 input; 286,134 output).
- Elapsed: **3769.0s wall clock**; 3366.4s summed per-paper route + read time.
- Notation twins merged before scoring: **1** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 49}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 155 | 61 | 87 | 71.8% | 64.0% | 67.7% |
| `linkage_assisted` | secondary_diagnostic | 209 | 120 | 33 | 63.5% | 86.4% | 73.2% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `tranche_01.tsv` (49 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 48 / 242 | 19.8% | 0.812 | 3.961 |
| affected | 101 / 242 | 41.7% | 0.079 | 0.345 |
| unaffected | 7 / 242 | 2.9% | 0.286 | 0.756 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 44 / 233 | 18.9% | 4 / 9 | 44.4% |
| affected | 97 / 222 | 43.7% | 4 / 20 | 20.0% |
| unaffected | 4 / 13 | 30.8% | 3 / 229 | 1.3% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 34 | 35 | 20 | 49.3% | 63.0% | 55.3% | 59.6% | 43.9% | 33.3% / 1.611 / 6.373 | 9.3% / 0.400 / 0.632 | 5.6% / 0.667 / 1.155 |
| KCNH2 | 17 | 11 | 8 | 60.7% | 68.0% | 64.2% | 65.4% | 50.0% | 32.0% / 0.625 / 1.275 | 20.0% / 0.000 / 0.000 | 0.0% / n/a / n/a |
| KCNQ1 | 10 | 9 | 56 | 52.6% | 15.2% | 23.5% | 90.9% | 83.3% | 7.6% / 0.000 / 0.000 | 4.5% / 0.667 / 0.816 | 1.5% / 0.000 / 0.000 |
| RYR2 | 21 | 5 | 2 | 80.8% | 91.3% | 85.7% | 100.0% | 100.0% | 73.9% / 0.294 / 0.728 | 65.2% / 0.267 / 0.730 | 13.0% / 0.000 / 0.000 |
| MYBPC3 | 73 | 1 | 1 | 98.6% | 98.6% | 98.6% | 98.6% | 98.6% | 0.0% / n/a / n/a | 98.6% / 0.000 / 0.000 | 0.0% / n/a / n/a |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| KCNH2 | 15466642 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.2 | 21,563 |
| SCN5A | 27232914 | text | 3 | 1 | 2 | 75.0% | 60.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.4 | 16,909 |
| SCN5A | 23283979 | text | 2 | 2 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 101.8 | 34,434 |
| SCN5A | 28638671 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 35.5 | 13,260 |
| KCNQ1 | 29857160 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 54.6 | 22,804 |
| SCN5A | 28779003 | text | 2 | 1 | 1 | 66.7% | 66.7% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 42.0 | 12,195 |
| SCN5A | 28469501 | text | 0 | 2 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 51.8 | 19,529 |
| KCNH2 | 29881912 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 18.8 | 7,401 |
| KCNH2 | 9544837 | text | 5 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 51.9 | 9,949 |
| SCN5A | 12566525 | text | 0 | 0 | 4 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.7 | 6,555 |
| KCNQ1 | 15176425 | text | 0 | 0 | 6 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.2 | 36,716 |
| SCN5A | 28637969 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.7 | 7,193 |
| KCNQ1 | 21138517 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 0.000 | 46.4 | 12,390 |
| SCN5A | 17081365 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 38.9 | 9,084 |
| KCNQ1 | 29677589 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 53.1 | 18,132 |
| KCNH2 | 22067087 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 36.9 | 10,595 |
| KCNH2 | 23936059 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 88.5 | 29,121 |
| SCN5A | 17971661 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 10.0 | 5,378 |
| RYR2 | 17062961 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 103.1 | 14,613 |
| KCNQ1 | 14678125 | text | 0 | 0 | 41 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 13.5 | 5,799 |
| SCN5A | 29062695 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 100.0% / 2.000 | 19.1 | 10,550 |
| RYR2 | 26230511 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 120.5 | 42,989 |
| SCN5A | 26498160 | text | 2 | 2 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 82.2 | 45,711 |
| SCN5A | 19272188 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 113.4 | 27,573 |
| KCNQ1 | 33664273 | text | 1 | 4 | 0 | 20.0% | 100.0% | 33.3% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 146.3 | 42,704 |
| KCNQ1 | 21131640 | text | 0 | 0 | 9 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.3 | 6,644 |
| SCN5A | 15579534 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 26.6 | 18,794 |
| KCNH2 | 28491588 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 2.500 | 0.0% / n/a | 0.0% / n/a | 67.2 | 17,092 |
| SCN5A | 23414114 | text | 1 | 21 | 0 | 4.5% | 100.0% | 8.7% | 100.0% / 27.000 | 0.0% / n/a | 0.0% / n/a | 373.5 | 98,106 |
| KCNH2 | 15851119 | text | 4 | 9 | 0 | 30.8% | 100.0% | 47.1% | 75.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 365.2 | 65,907 |
| SCN5A | 23200271 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 51.9 | 11,986 |
| RYR2 | 22222782 | text | 7 | 4 | 0 | 63.6% | 100.0% | 77.8% | 42.9% / 0.333 | 14.3% / 0.000 | 0.0% / n/a | 240.6 | 38,869 |
| KCNQ1 | 12710526 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.6 | 7,556 |
| RYR2 | 15466642 | text | 8 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.500 | 100.0% / 0.500 | 0.0% / n/a | 177.7 | 45,044 |
| SCN5A | 23936059 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 73.4 | 26,403 |
| KCNH2 | 16470702 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 11.5 | 6,898 |
| SCN5A | 23349452 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 95.8 | 41,603 |
| SCN5A | 16764707 | text | 7 | 0 | 3 | 100.0% | 70.0% | 82.4% | 70.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 82.9 | 26,866 |
| SCN5A | 12574983 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 13.0 | 5,953 |
| MYBPC3 | 33673806 | text | 73 | 1 | 1 | 98.6% | 98.6% | 98.6% | 0.0% / n/a | 98.6% / 0.000 | 0.0% / n/a | 6.7 | 3,360 |
| SCN5A | 15176425 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.2 | 5,873 |
| SCN5A | 14654377 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 40.1 | 15,662 |
| KCNQ1 | 30508507 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 40.6 | 15,459 |
| SCN5A | 10727653 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 87.1 | 19,244 |
| KCNH2 | 15176425 | text | 0 | 0 | 7 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 13.1 | 36,444 |
| RYR2 | 27839804 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 35.6 | 12,196 |
| RYR2 | 17556193 | text | 2 | 0 | 2 | 100.0% | 50.0% | 66.7% | 50.0% / 0.000 | 50.0% / 0.000 | 25.0% / 0.000 | 56.7 | 26,191 |
| KCNQ1 | 23092362 | text | 3 | 3 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 92.6 | 23,976 |
| SCN5A | 11804990 | text | 0 | 0 | 2 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.1 | 19,661 |

## Errors and representation choices

### KCNH2 PMID 15466642

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: I419F

### SCN5A PMID 27232914

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R1629G, c.3391-1G>A
- Extra predictions: p.Arg1629Lys

### SCN5A PMID 23283979

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: C373Y, H886Q

### SCN5A PMID 28638671

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg1193Gln carriers 1 vs 0 (error +1); p.Arg1193Gln affected 1 vs 0 (error +1)

### KCNQ1 PMID 29857160

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 28779003

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: I1448N
- Extra predictions: c.4340T>A

### SCN5A PMID 28469501

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G289S
- Extra predictions: p.Lys1573Glu c.4717A>G, p.Ser502Leu c.1505C>T

### KCNH2 PMID 29881912

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 9544837

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 12566525

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1330T, I1768V, P.Y1795_E1796INSD, Q692K

### KCNQ1 PMID 15176425

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: 1685DELAG+1DELG, 277DELS, D202H, R231C, S547L, W248C

### SCN5A PMID 28637969

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 21138517

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17081365

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 29677589

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg190Trp c.568C>T affected 1 vs 0 (error +1); p.Arg594Gln c.1781G>A affected 1 vs 0 (error +1)

### KCNH2 PMID 22067087

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 23936059

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 17971661

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: E1784K

### RYR2 PMID 17062961

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 14678125

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A178P, A300T, A341E, A341V, A344/SP, A344A/SPLICE, D317G, DELF340, E284K, G168R, G185S, G269D, G269S, G306R, G345R, I567S, K393N, L151, L266P, L273F, Q357H, Q530X, R190Q, R243C, R366Q, R366W, R518X, R539W, R594Q, S225L, S349W, S566F, T144A, T311I, T312I, V172M, V254L, V254M, V524G, Y171X, Y315C

### SCN5A PMID 29062695

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Phe532Cys unaffected 2 vs 0 (error +2)

### RYR2 PMID 26230511

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26498160

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Gly1275Ser, p.Leu1626Pro

### SCN5A PMID 19272188

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.5626insATG, p.I848SfsX33 c.2540delC

### KCNQ1 PMID 33664273

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.1686-9T>C, c.1733-357C>T, c.1733-362_1733-361delAGinsGGA, c.1733-402A>C

### KCNQ1 PMID 21131640

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G168R, G316E, L266P, R243C, R360G, R518X, S333CFS/129X, S546L, S566F

### SCN5A PMID 15579534

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Glu1053Lys c.3157G>A carriers 1 vs 0 (error +1); p.Glu1053Lys c.3157G>A affected 1 vs 0 (error +1)

### KCNH2 PMID 28491588

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: N588K
- Count disagreements: p.K897T carriers 3 vs 1 (error +2); p.T618I c.1853C>T carriers 3 vs 0 (error +3)

### SCN5A PMID 23414114

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.A1924T c.5770G>A, p.D1243N c.3727G>A, p.E1240Q c.3718G>C, p.E161K c.481G>A, p.F1293S c.3878T>C, p.G1319V c.3956G>T, p.G35S c.103G>A, p.G615E c.1844G>A, p.I1968S c.5903T>G, p.L1308F c.3922C>T, p.L619F c.1855C>T, p.P717L c.2150C>T, p.Q1832E c.5494C>G, p.R376H c.1127G>A, p.R526H c.1577G>A, p.R661W c.1981C>T, p.S216L c.647C>T, p.T220I c.659C>T, p.V1525M c.4573G>A, p.V1951L c.5851G>T, p.V232I c.694G>A
- Count disagreements: p.F2004L c.6010T>C carriers 27 vs 0 (error +27)

### KCNH2 PMID 15851119

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: Y475del, p.Arg1033fs, p.Asn588Asp, p.Asn633Ser, p.Asn996Ile, p.Asp837Gly, p.Gly572Ser, p.Thr65Pro, p.Val822Met

### SCN5A PMID 23200271

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 22222782

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Asp4646Ala, p.Gln4879His, p.Glu4146Lys, p.Ile4855Met
- Count disagreements: p.Q486H carriers 1 vs 2 (error -1)

### KCNQ1 PMID 12710526

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 15466642

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Ile419Phe carriers 1 vs 3 (error -2); p.Ile4848Val carriers 6 vs 4 (error +2); p.Ile419Phe affected 1 vs 3 (error -2); p.Ile4848Val affected 6 vs 4 (error +2)

### SCN5A PMID 23936059

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.2738C>T

### KCNH2 PMID 16470702

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: 7Q34Q36.2DEL

### SCN5A PMID 23349452

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: Arg454Trp, Thr24Ile

### SCN5A PMID 16764707

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: L860fsX, P.1570_F1571INSI, P.Y1795_E1796INSD

### SCN5A PMID 12574983

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R1623Q

### MYBPC3 PMID 33673806

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: p.Glu843_Arg845del
- Extra predictions: c.2528_2536del

### SCN5A PMID 15176425

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A185T, A691T, I239V, R190G, R340Q

### SCN5A PMID 14654377

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: E1295K

### KCNQ1 PMID 30508507

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.1253_1259delACTCTGG, c.848C>A

### SCN5A PMID 10727653

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 15176425

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: C39X, D460fsX, E1053fsX, G572S, IVS13+1G>A, R273X, W497X

### RYR2 PMID 27839804

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H29D

### RYR2 PMID 17556193

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R2267H, S4565R

### KCNQ1 PMID 23092362

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Gly88Arg, p.Thr85Ala, p.Thr85Met

### SCN5A PMID 11804990

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: D1595N, G298S

## Scope, method, and limitations

- Population: fixed manifest `tranche_01.tsv` (49 papers); per-gene counts {'SCN5A': 23, 'KCNH2': 9, 'KCNQ1': 10, 'RYR2': 6, 'MYBPC3': 1}; every PMID has downloaded source and at least one gold assertion in each count field.
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
