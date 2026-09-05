# Codex extraction-blinded paper evaluation — `20260905_mechanism10_candidate`

## Technical summary

This hash-locked run evaluated **10 papers** (**per-gene counts {'SCN5A': 4, 'KCNQ1': 2, 'RYR2': 3, 'MYBPC3': 1}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **89.5%**, recall **96.7%**, F1 **92.9%** (348 TP, 41 FP, 12 FN).
- Precision versus counted extras **97.5%** (348 matched rows; 9 extra rows with patient counts). The stricter count-bearing-only diagnostic is **96.4%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **514,748 total tokens** (385,191 input; 129,557 output).
- Elapsed: **1551.0s wall clock**; 1271.3s summed per-paper route + read time.
- Notation twins merged before scoring: **0** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 10}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 348 | 41 | 12 | 89.5% | 96.7% | 92.9% |
| `linkage_assisted` | secondary_diagnostic | 348 | 89 | 12 | 79.6% | 96.7% | 87.3% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `mechanism_panel.tsv` (10 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 240 / 360 | 66.7% | 0.017 | 0.129 |
| affected | 56 / 360 | 15.6% | 0.161 | 0.401 |
| unaffected | 15 / 360 | 4.2% | 0.000 | 0.000 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 240 / 276 | 87.0% | 0 / 84 | 0.0% |
| affected | 47 / 265 | 17.7% | 9 / 95 | 9.5% |
| unaffected | 9 / 15 | 60.0% | 6 / 345 | 1.7% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 287 | 27 | 6 | 91.4% | 98.0% | 94.6% | 98.6% | 98.0% | 66.2% / 0.021 / 0.144 | 4.1% / 0.000 / 0.000 | 2.4% / 0.000 / 0.000 |
| KCNQ1 | 10 | 5 | 2 | 66.7% | 83.3% | 74.1% | 100.0% | 100.0% | 75.0% / 0.000 / 0.000 | 75.0% / 1.000 / 1.000 | 50.0% / 0.000 / 0.000 |
| RYR2 | 25 | 2 | 2 | 92.6% | 92.6% | 92.6% | 92.6% | 91.7% | 81.5% / 0.000 / 0.000 | 81.5% / 0.000 / 0.000 | 7.4% / 0.000 / 0.000 |
| MYBPC3 | 26 | 7 | 2 | 78.8% | 92.9% | 85.2% | 89.7% | 83.3% | 53.6% / 0.000 / 0.000 | 46.4% / 0.000 / 0.000 | 0.0% / n/a / n/a |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 20031634 | text | 10 | 2 | 3 | 83.3% | 76.9% | 80.0% | 53.8% / 0.571 | 53.8% / 0.000 | 53.8% / 0.000 | 193.0 | 66,968 |
| KCNQ1 | 26481773 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 22.5 | 13,986 |
| RYR2 | 22677073 | text | 23 | 2 | 2 | 92.0% | 92.0% | 92.0% | 80.0% / 0.000 | 80.0% / 0.000 | 0.0% / n/a | 256.8 | 105,979 |
| KCNQ1 | 22677073 | text | 10 | 5 | 1 | 66.7% | 90.9% | 76.9% | 81.8% / 0.000 | 81.8% / 1.000 | 54.5% / 0.000 | 215.9 | 78,214 |
| SCN5A | 22677073 | text | 12 | 2 | 0 | 85.7% | 100.0% | 92.3% | 41.7% / 0.000 | 41.7% / 0.000 | 0.0% / n/a | 115.7 | 54,702 |
| MYBPC3 | 21302287 | text | 26 | 7 | 2 | 78.8% | 92.9% | 85.2% | 53.6% / 0.000 | 46.4% / 0.000 | 0.0% / n/a | 316.6 | 109,599 |
| SCN5A | 30059973 | text | 182 | 3 | 3 | 98.4% | 98.4% | 98.4% | 98.4% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 32533946 | text | 83 | 20 | 0 | 80.6% | 100.0% | 89.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 87.7 | 63,548 |
| RYR2 | 19398417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 24.0 | 8,453 |
| RYR2 | 25435091 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 39.1 | 13,299 |

## Errors and representation choices

### SCN5A PMID 20031634

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.1983_1993dupGGCCCTCAGCG, c.3963+2T>C, P.1570_F1571INSI
- Extra predictions: c.3963+12T>C, p.Ala665GlyfsX16 c.1983-1993dup
- Count disagreements: p.Arg225Trp c.673C>T carriers 11 vs 10 (error +1); p.Asn1722Asp c.5164A>G carriers 9 vs 8 (error +1); p.Gly1408Arg c.4222G>A carriers 14 vs 13 (error +1); p.Ser1382Ile c.4145G>T carriers 9 vs 8 (error +1)

### KCNQ1 PMID 26481773

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: R594Q

### RYR2 PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: Q2958R, R4037C
- Extra predictions: c.2398+5G>T, c.6739C>T

### KCNQ1 PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: L162fsX
- Extra predictions: G643S, K393N, P448R, V648I, p.L163fs*236
- Count disagreements: A302V affected 1 vs 0 (error +1); A341E affected 1 vs 0 (error +1); G269S affected 1 vs 0 (error +1); G584S affected 1 vs 0 (error +1); I198V affected 1 vs 0 (error +1); T312I affected 1 vs 0 (error +1); T587M affected 1 vs 0 (error +1); W379R affected 1 vs 0 (error +1); Y111C affected 1 vs 0 (error +1)

### SCN5A PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: R18W, S1787N

### MYBPC3 PMID 21302287

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.2258_2259insT, c.506-12delC
- Extra predictions: F244F c.732C>T, K754E, K754EfsX78, L517M c.1549C>A, c.2258-2259InsT, c.2846-2847InsT, c.3192-3193InsC

### SCN5A PMID 30059973

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.G1015DFSX14, P.I1570DUP, P.L1302VFS18
- Extra predictions: Leu1302Valfs18 c.3900_3903dup, c.3045_3046del, p.Ile1570dup c.4708_4710dup

### SCN5A PMID 32533946

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: D1430N, D356N, G1408R, G1712C, G1740R, G1743E, G1743R, G897E, I1660V, L846R, LEU839P, R104Q, R104W, R282H, R878C, R878H, R893H, S1218I, S910L, T187I

### RYR2 PMID 19398417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 25435091

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

## Scope, method, and limitations

- Population: fixed manifest `mechanism_panel.tsv` (10 papers); per-gene counts {'SCN5A': 4, 'KCNQ1': 2, 'RYR2': 3, 'MYBPC3': 1}; every PMID has downloaded source and at least one gold assertion in each count field.
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
