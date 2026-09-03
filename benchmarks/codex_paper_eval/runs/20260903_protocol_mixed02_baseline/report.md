# Codex extraction-blinded paper evaluation — `20260903_protocol_mixed02_baseline`

## Technical summary

This hash-locked run evaluated **49 papers** (**per-gene counts {'SCN5A': 24, 'KCNH2': 8, 'KCNQ1': 10, 'RYR2': 6, 'BRCA1': 1}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **83.2%**, recall **88.4%**, F1 **85.8%** (268 TP, 54 FP, 35 FN).
- Precision versus counted extras **97.5%** (268 matched rows; 7 extra rows with patient counts). The stricter count-bearing-only diagnostic is **97.1%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **1,124,321 total tokens** (842,737 input; 281,584 output).
- Elapsed: **3275.0s wall clock**; 2617.6s summed per-paper route + read time.
- Notation twins merged before scoring: **0** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 49}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 268 | 54 | 35 | 83.2% | 88.4% | 85.8% |
| `linkage_assisted` | secondary_diagnostic | 278 | 150 | 25 | 65.0% | 91.7% | 76.1% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `tranche_02.tsv` (49 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 232 / 303 | 76.6% | 0.026 | 0.161 |
| affected | 35 / 303 | 11.6% | 0.200 | 0.447 |
| unaffected | 6 / 303 | 2.0% | 0.333 | 0.577 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 231 / 285 | 81.1% | 1 / 18 | 5.6% |
| affected | 29 / 276 | 10.5% | 6 / 27 | 22.2% |
| unaffected | 3 / 16 | 18.8% | 3 / 287 | 1.0% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 213 | 6 | 32 | 97.3% | 86.9% | 91.8% | 98.6% | 98.5% | 81.2% / 0.015 / 0.123 | 4.9% / 0.000 / 0.000 | 1.2% / 0.333 / 0.577 |
| KCNH2 | 14 | 8 | 2 | 63.6% | 87.5% | 73.7% | 100.0% | 100.0% | 18.8% / 0.000 / 0.000 | 18.8% / 0.667 / 0.816 | 0.0% / n/a / n/a |
| KCNQ1 | 18 | 22 | 0 | 45.0% | 100.0% | 62.1% | 90.0% | 77.8% | 38.9% / 0.143 / 0.378 | 33.3% / 0.833 / 0.913 | 16.7% / 0.333 / 0.577 |
| RYR2 | 23 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% / 0.087 / 0.295 | 60.9% / 0.000 / 0.000 | 0.0% / n/a / n/a |
| BRCA1 | 0 | 18 | 1 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% / n/a / n/a | 0.0% / n/a / n/a | 0.0% / n/a / n/a |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 27784737 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 21.4 | 7,972 |
| SCN5A | 22835975 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 40.9 | 11,368 |
| BRCA1 | 18627636 | text | 0 | 18 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 262.9 | 100,469 |
| KCNQ1 | 32508908 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 90.1 | 26,524 |
| SCN5A | 23424222 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 77.5 | 20,104 |
| SCN5A | 15057319 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 57.5 | 15,957 |
| KCNH2 | 29330128 | text | 1 | 6 | 0 | 14.3% | 100.0% | 25.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 51.4 | 16,514 |
| SCN5A | 23168001 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 25.6 | 10,992 |
| SCN5A | 26129877 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 35.3 | 26,002 |
| KCNQ1 | 33498651 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 100.0% / 0.000 | 80.1 | 28,930 |
| RYR2 | 21954897 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 25.0% / 0.000 | 0.0% / n/a | 40.0 | 20,768 |
| SCN5A | 26803770 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 39.3 | 24,097 |
| SCN5A | 25341504 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 141.9 | 18,846 |
| SCN5A | 12106943 | text | 0 | 0 | 16 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 21.3 | 94,078 |
| SCN5A | 28815794 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 37.7 | 20,280 |
| KCNH2 | 19215240 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.6 | 8,746 |
| RYR2 | 21652165 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 25.4 | 11,536 |
| RYR2 | 29544605 | text | 12 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 159.8 | 71,591 |
| KCNQ1 | 27868350 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 1.000 | 100.0% / 0.000 | 50.1 | 19,511 |
| SCN5A | 24599044 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 16.8 | 10,348 |
| KCNQ1 | 29544605 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 39.5 | 31,907 |
| SCN5A | 30116708 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 21.4 | 8,943 |
| SCN5A | 11463728 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 48.8 | 25,303 |
| SCN5A | 20025708 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 48.3 | 26,724 |
| KCNH2 | 29752375 | text | 8 | 1 | 0 | 88.9% | 100.0% | 94.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 69.7 | 23,080 |
| KCNQ1 | 29922582 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 100.0% / 1.000 | 43.5 | 15,103 |
| KCNH2 | 26831020 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.5 | 14,416 |
| SCN5A | 30059973 | text | 182 | 3 | 3 | 98.4% | 98.4% | 98.4% | 97.8% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 22710484 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 1.000 | 100.0% / 0.000 | 100.0% / 1.000 | 24.5 | 15,525 |
| KCNH2 | 15572050 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 25.7 | 18,634 |
| RYR2 | 27491078 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 64.4 | 33,782 |
| SCN5A | 27281089 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 37.4 | 14,627 |
| SCN5A | 25650408 | text | 0 | 0 | 12 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.6 | 6,236 |
| KCNQ1 | 29330128 | text | 1 | 5 | 0 | 16.7% | 100.0% | 28.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 61.3 | 17,321 |
| SCN5A | 18242854 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 35.8 | 26,022 |
| SCN5A | 10200053 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 31.5 | 10,279 |
| KCNH2 | 26022375 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.2 | 6,971 |
| KCNH2 | 14642688 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 6.1 | 7,634 |
| SCN5A | 29544605 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 105.5 | 59,184 |
| KCNH2 | 26129877 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 0.0% / n/a | 48.3 | 31,353 |
| SCN5A | 21167004 | text | 2 | 1 | 1 | 66.7% | 66.7% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 111.4 | 35,642 |
| RYR2 | 27646203 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 0.0% / n/a | 0.0% / n/a | 32.8 | 14,324 |
| KCNQ1 | 20421371 | text | 7 | 15 | 0 | 31.8% | 100.0% | 48.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 126.7 | 38,436 |
| SCN5A | 20395683 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 27.6 | 8,052 |
| KCNQ1 | 14527360 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 56.9 | 9,220 |
| SCN5A | 15028074 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 26.8 | 8,225 |
| KCNQ1 | 26168993 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 0.000 | 0.0% / n/a | 28.3 | 10,394 |
| RYR2 | 22519458 | text | 4 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.250 | 0.0% / n/a | 0.0% / n/a | 107.9 | 32,562 |
| KCNQ1 | 20138589 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 29.6 | 9,789 |

## Errors and representation choices

### SCN5A PMID 27784737

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 22835975

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: H558R carriers 1 vs 0 (error +1)

### BRCA1 PMID 18627636

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: 180delA
- Extra predictions: E1221A c.3781A>C, E352X c.1173G>T, E402X c.1323G>T, I124I c.491C>A, N1268S c.3922A>G, N909I c.2845A>T, N913K c.2858T>A, P346S c.1155C>T, P977P c.3050A>G, Q1420X c.4377C>T, R1751X c.5370C>T, R1835Q c.5623G>A, R252C c.873C>T, R762S c.2405A>T, S1631N c.5011G>A, S265S c.914T>C, V191I c.690G>A, Y856H c.2685T>C

### KCNQ1 PMID 32508908

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: c.477+5G>A

### SCN5A PMID 23424222

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15057319

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 29330128

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: A57P, G53D, G53R, G53S, R56Q, Y54H

### SCN5A PMID 23168001

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R1023H

### SCN5A PMID 26129877

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 33498651

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: P631fs*20 affected 1 vs 0 (error +1)

### RYR2 PMID 21954897

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26803770

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 25341504

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 12106943

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1924T, D951X, E1225K, E161K, G1319V, G1406R, G1502S, G1743E, G752R, L867X, M369K, R1512W, R367C, R535X, S1382I, V1405L

### SCN5A PMID 28815794

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 19215240

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 21652165

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 29544605

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 27868350

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Lys258fs
- Count disagreements: p.Arg259Gly c.775C>G affected 1 vs 0 (error +1)

### SCN5A PMID 24599044

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 29544605

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.R518Q c.1553G>A affected 1 vs 0 (error +1); p.S77F c.230C>T affected 1 vs 0 (error +1)

### SCN5A PMID 30116708

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 11463728

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 20025708

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: P1008S carriers 3 vs 4 (error -1)

### KCNH2 PMID 29752375

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R1047L

### KCNQ1 PMID 29922582

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Arg518* affected 2 vs 1 (error +1); p.Arg518* unaffected 1 vs 2 (error -1)

### KCNH2 PMID 26831020

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R863X

### SCN5A PMID 30059973

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.G1015DFSX14, P.I1570DUP, P.L1302VFS18
- Extra predictions: Leu1302Valfs18 c.3900_3903dup, c.3045_3046del, p.Ile1570dup c.4708_4710dup

### SCN5A PMID 22710484

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: H558R
- Count disagreements: p.Arg222Gln c.665G>A carriers 3 vs 2 (error +1); p.Arg222Gln c.665G>A unaffected 1 vs 0 (error +1)

### KCNH2 PMID 15572050

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 27491078

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 27281089

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 25650408

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1416E, C335S, D1041N, E1780G, E312K, H445D, I1660S, I397V, M704T, M734V, R179X, R475K

### KCNQ1 PMID 29330128

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: E160K, E160V, F157C, G168R, T169R

### SCN5A PMID 18242854

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 10200053

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 26022375

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: T1083fsX

### KCNH2 PMID 14642688

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: E682X

### SCN5A PMID 29544605

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNH2 PMID 26129877

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Thr436Met affected 1 vs 0 (error +1); p.Thr895Met affected 1 vs 0 (error +1)

### SCN5A PMID 21167004

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: I1836T
- Extra predictions: I1835T

### RYR2 PMID 27646203

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Ile4855Met carriers 3 vs 2 (error +1)

### KCNQ1 PMID 20421371

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: D202A, D202E, D202K, D202M, I204A, I204N, S209A, S209L, S209M, V205A, V205L, V207M, V215A, V215L, V215P

### SCN5A PMID 20395683

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 14527360

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 15028074

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 26168993

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Phe279Ile carriers 1 vs 2 (error -1)

### RYR2 PMID 22519458

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.M4107V c.12643A>G carriers 1 vs 2 (error -1)

### KCNQ1 PMID 20138589

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

## Scope, method, and limitations

- Population: fixed manifest `tranche_02.tsv` (49 papers); per-gene counts {'SCN5A': 24, 'KCNH2': 8, 'KCNQ1': 10, 'RYR2': 6, 'BRCA1': 1}; every PMID has downloaded source and at least one gold assertion in each count field.
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
