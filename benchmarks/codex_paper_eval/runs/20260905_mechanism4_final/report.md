# Codex extraction-blinded paper evaluation — `20260905_mechanism4_final`

## Technical summary

This hash-locked run evaluated **4 papers** (**per-gene counts {'SCN5A': 2, 'KCNQ1': 1, 'RYR2': 1}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **90.2%**, recall **90.2%**, F1 **90.2%** (55 TP, 6 FP, 6 FN).
- Precision versus counted extras **93.2%** (55 matched rows; 4 extra rows with patient counts). The stricter count-bearing-only diagnostic is **92.5%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **335,739 total tokens** (254,218 input; 81,521 output).
- Elapsed: **1148.0s wall clock**; 857.0s summed per-paper route + read time.
- Notation twins merged before scoring: **0** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 4}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 55 | 6 | 6 | 90.2% | 90.2% | 90.2% |
| `linkage_assisted` | secondary_diagnostic | 56 | 16 | 5 | 77.8% | 91.8% | 84.2% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `mechanism_ablation.tsv` (4 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 49 / 61 | 80.3% | 0.122 | 0.350 |
| affected | 49 / 61 | 80.3% | 0.224 | 0.474 |
| unaffected | 30 / 61 | 49.2% | 0.000 | 0.000 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 49 / 61 | 80.3% | 0 / 0 | n/a |
| affected | 39 / 50 | 78.0% | 10 / 11 | 90.9% |
| unaffected | 8 / 13 | 61.5% | 22 / 48 | 45.8% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 23 | 2 | 2 | 92.0% | 92.0% | 92.0% | 95.8% | 95.0% | 76.0% / 0.316 / 0.562 | 76.0% / 0.053 / 0.229 | 52.0% / 0.000 / 0.000 |
| KCNQ1 | 10 | 2 | 1 | 83.3% | 90.9% | 87.0% | 90.9% | 90.9% | 90.9% / 0.000 / 0.000 | 90.9% / 1.000 / 1.000 | 0.0% / n/a / n/a |
| RYR2 | 22 | 2 | 3 | 91.7% | 88.0% | 89.8% | 91.7% | 90.9% | 80.0% / 0.000 / 0.000 | 80.0% / 0.000 / 0.000 | 68.0% / 0.000 / 0.000 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 20031634 | text | 11 | 2 | 2 | 84.6% | 84.6% | 84.6% | 61.5% / 0.625 | 61.5% / 0.000 | 61.5% / 0.000 | 155.4 | 65,324 |
| KCNQ1 | 22677073 | text | 10 | 2 | 1 | 83.3% | 90.9% | 87.0% | 90.9% / 0.000 | 90.9% / 1.000 | 0.0% / n/a | 153.2 | 76,118 |
| RYR2 | 22677073 | text | 22 | 2 | 3 | 91.7% | 88.0% | 89.8% | 80.0% / 0.000 | 80.0% / 0.000 | 68.0% / 0.000 | 308.0 | 111,328 |
| SCN5A | 22677073 | text | 12 | 0 | 0 | 100.0% | 100.0% | 100.0% | 91.7% / 0.091 | 91.7% / 0.091 | 41.7% / 0.000 | 240.4 | 82,969 |

## Errors and representation choices

### SCN5A PMID 20031634

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.1983_1993dupGGCCCTCAGCG, c.3963+2T>C
- Extra predictions: c.3963+12T>C, p.Ala665GlyfsX16 c.1983_1993dup
- Count disagreements: p.1570_F1571insI carriers 11 vs 10 (error +1); p.Arg225Trp c.673C>T carriers 11 vs 10 (error +1); p.Asn1722Asp c.5164A>G carriers 9 vs 8 (error +1); p.Gly1408Arg c.4222G>A carriers 14 vs 13 (error +1); p.Ser1382Ile c.4145G>T carriers 9 vs 8 (error +1)

### KCNQ1 PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: L162fsX
- Extra predictions: G643S, p.L163fs*236 c.488delT
- Count disagreements: A302V affected 1 vs 0 (error +1); A341E affected 1 vs 0 (error +1); G269S affected 1 vs 0 (error +1); G584S affected 1 vs 0 (error +1); I198V affected 1 vs 0 (error +1); T312I affected 1 vs 0 (error +1); T587M affected 1 vs 0 (error +1); W379R affected 1 vs 0 (error +1); Y111C affected 1 vs 0 (error +1); p.L191fs*90 affected 1 vs 0 (error +1)

### RYR2 PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: Q2958R, A1136V, R4037C
- Extra predictions: c.2398+5G>T, c.6739C>T

### SCN5A PMID 22677073

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: S216L carriers 2 vs 1 (error +1); S216L affected 2 vs 1 (error +1)

## Scope, method, and limitations

- Population: fixed manifest `mechanism_ablation.tsv` (4 papers); per-gene counts {'SCN5A': 2, 'KCNQ1': 1, 'RYR2': 1}; every PMID has downloaded source and at least one gold assertion in each count field.
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
