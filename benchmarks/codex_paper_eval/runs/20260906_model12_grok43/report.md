# Codex extraction-blinded paper evaluation — `20260906_model12_grok43`

## Technical summary

This hash-locked run evaluated **12 papers** (**per-gene counts {'SCN5A': 5, 'RYR2': 5, 'MYBPC3': 2}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **95.1%**, recall **88.7%**, F1 **91.8%** (699 TP, 36 FP, 89 FN).
- Precision versus counted extras **98.2%** (699 matched rows; 13 extra rows with patient counts). The stricter count-bearing-only diagnostic is **97.7%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **447,788 total tokens** (310,939 input; 136,849 output).
- Elapsed: **1216.0s wall clock**; 975.2s summed per-paper route + read time.
- Notation twins merged before scoring: **1** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 12}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 699 | 36 | 89 | 95.1% | 88.7% | 91.8% |
| `linkage_assisted` | secondary_diagnostic | 699 | 368 | 89 | 65.5% | 88.7% | 75.4% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `paper_manifest.tsv` (12 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 553 / 788 | 70.2% | 0.179 | 0.723 |
| affected | 24 / 788 | 3.0% | 1.250 | 5.050 |
| unaffected | 69 / 787 | 8.8% | 0.739 | 5.364 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 553 / 705 | 78.4% | 0 / 83 | 0.0% |
| affected | 24 / 650 | 3.7% | 0 / 138 | 0.0% |
| unaffected | 69 / 86 | 80.2% | 0 / 701 | 0.0% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 635 | 31 | 83 | 95.3% | 88.4% | 91.8% | 98.3% | 98.0% | 74.0% / 0.186 / 0.738 | 1.4% / 0.600 / 1.897 | 9.1% / 0.108 / 0.868 |
| RYR2 | 26 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 19.2% / 0.000 / 0.000 | 11.5% / 8.000 / 13.856 | 12.0% / 14.667 / 25.403 |
| MYBPC3 | 38 | 5 | 6 | 88.4% | 86.4% | 87.4% | 95.0% | 90.0% | 38.6% / 0.000 / 0.000 | 25.0% / 0.000 / 0.000 | 2.3% / 0.000 / 0.000 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| RYR2 | 18929323 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 24.6 | 16,701 |
| SCN5A | 30059973 | text | 182 | 3 | 3 | 98.4% | 98.4% | 98.4% | 98.4% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| MYBPC3 | 20433692 | text | 13 | 0 | 3 | 100.0% | 81.2% | 89.7% | 6.2% / 0.000 | 6.2% / 0.000 | 6.2% / 0.000 | 234.0 | 75,380 |
| MYBPC3 | 21302287 | text | 25 | 5 | 3 | 83.3% | 89.3% | 86.2% | 57.1% / 0.000 | 35.7% / 0.000 | 0.0% / n/a | 208.6 | 109,479 |
| RYR2 | 30403697 | text | 21 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 120.8 | 44,214 |
| SCN5A | 20031634 | text | 11 | 1 | 2 | 91.7% | 84.6% | 88.0% | 76.9% / 0.500 | 76.9% / 0.600 | 76.9% / 0.700 | 113.8 | 63,338 |
| RYR2 | 19398417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 22.3 | 7,859 |
| RYR2 | 25435091 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 32.1 | 13,324 |
| SCN5A | 32533946 | text | 83 | 20 | 0 | 80.6% | 100.0% | 89.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 92.8 | 63,545 |
| SCN5A | 25163546 | text | 20 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 92.3 | 27,711 |
| RYR2 | 25814417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 24.000 | 100.0% / 44.000 | 24.7 | 23,603 |
| SCN5A | 20129283 | text | 339 | 7 | 78 | 98.0% | 81.3% | 88.9% | 81.3% / 0.277 | 0.0% / n/a | 13.2% / 0.000 | 9.0 | 2,634 |

## Errors and representation choices

### RYR2 PMID 18929323

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 30059973

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.G1015DFSX14, P.I1570DUP, P.L1302VFS18
- Extra predictions: Leu1302Valfs18 c.3900_3903dup, c.3045_3046del, p.Ile1570dup c.4708_4710dup

### MYBPC3 PMID 20433692

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: IVS11-9G>A, IVS29+5G>A, IVS6+5G>A

### MYBPC3 PMID 21302287

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: c.2258_2259insT, c.3192_3193insC, c.506-12delC
- Extra predictions: F244F c.732C>T, K754E, K754EfsX78, L517M c.1549C>A, p.K1065QfsX11 c.3192-3193InsC

### RYR2 PMID 30403697

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 20031634

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.1570_F1571INSI, c.1983_1993dupGGCCCTCAGCG
- Extra predictions: p.Ala665GlyfsX16 c.1983_1993dup
- Count disagreements: c.3963+2T>C carriers 10 vs 9 (error +1); p.Arg225Trp c.673C>T carriers 11 vs 10 (error +1); p.Asn1722Asp c.5164A>G carriers 9 vs 8 (error +1); p.Gly1408Arg c.4222G>A carriers 14 vs 13 (error +1); p.Ser1382Ile c.4145G>T carriers 9 vs 8 (error +1); c.3963+2T>C affected 2 vs 8 (error -6); c.3963+2T>C unaffected 8 vs 1 (error +7)

### RYR2 PMID 19398417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 25435091

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 32533946

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: D1430N, D356N, G1408R, G1712C, G1740R, G1743E, G1743R, G897E, I1660V, L846R, LEU839P, R104Q, R104W, R282H, R878C, R878H, R893H, S1218I, S910L, T187I

### SCN5A PMID 25163546

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 25814417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Count disagreements: p.Gly357Ser c.1069G>A affected 97 vs 73 (error +24); p.Gly357Ser c.1069G>A unaffected 62 vs 106 (error -44)

### SCN5A PMID 20129283

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1680T, A735V, A735V, A735V, D1243N, D1243N, D1275N, E1225K, E1225K, E1225K, E161K, E746K, E746K, G1319V, G1408R, G1408R, G1661R, G1743E, G1743R, G1743R, G1743R, G752R, I1660V, L1393X, L1393X, L868X, N927S, P.F2004DUP, P.I137_C139DUP, P.I1570DUP, P.K1493DEL, P.Y1795_E1796INSD, P336L, Q1695X, R104Q, R104Q, R104W, R121Q, R1232W, R1232W, R1583C, R1623X, R1638X, R225W, R367H, R367H, R367H, R376H, R376H, R526H, R535X, R535X, R535X, R878H, R878H, R878H, R878H, R893H, R893H, R965C, R965C, S1672Y, T1620M, T1709M, T220I, T632M, V1405M, V232I, c.1890+5G>A, c.1936delC, c.1936delC, c.2435_2436+3delTGGTAinsCGCCT, c.2549_2550insTG, c.2602delC, c.2602delC, c.3667delG, c.3840+1G>A, c.4437+5G>A
- Extra predictions: A1223PfsX7 c.3666delG, F2004dup c.6010_6012dupTTC, I137_C139dup c.410_418dupTCATGTGCA, I1570dup c.4708_4710dupATC, V1525M c.4573G>A, c.2435_24363delTGGTAinsCGCCT, c.5387_5388insTGA
- Count disagreements: A1680T c.5038G>A carriers 2 vs 1 (error +1); A226V carriers 1 vs 2 (error -1); A735V c.2204C>T carriers 4 vs 1 (error +3); D1275N c.3823G>A carriers 3 vs 1 (error +2); E1225K c.3673G>A carriers 4 vs 1 (error +3); E746K c.2236G>A carriers 3 vs 1 (error +2); G1319V c.3956G>T carriers 5 vs 1 (error +4); G1743E c.5228G>A carriers 6 vs 1 (error +5); G1743R c.5227G>A carriers 5 vs 1 (error +4); I1660V carriers 1 vs 4 (error -3); K1493X c.4477A>T carriers 1 vs 2 (error -1); K1493del c.4477_4479delAAG carriers 2 vs 1 (error +1); L1393X c.4178T>A carriers 3 vs 1 (error +2); L868X c.2602delC carriers 2 vs 1 (error +1); N927S c.2780A>G carriers 3 vs 1 (error +2); P336L c.1007C>T carriers 2 vs 1 (error +1); Q1695X c.5083C>T carriers 2 vs 1 (error +1); Q646RfsX5 c.1936delC carriers 3 vs 1 (error +2); R104Q c.311G>A carriers 3 vs 1 (error +2); R104W c.310C>T carriers 2 vs 1 (error +1); R121Q c.362G>A carriers 2 vs 1 (error +1); R1232W c.3694C>T carriers 3 vs 1 (error +2); R1623X c.4867C>T carriers 2 vs 1 (error +1); R1638X c.4912C>T carriers 3 vs 1 (error +2); R367H c.1100G>A carriers 6 vs 1 (error +5); R376H c.1127G>A carriers 4 vs 1 (error +3); R526H c.1577G>A carriers 2 vs 1 (error +1); R878H c.2633G>A carriers 5 vs 1 (error +4); R893H c.2678G>A carriers 3 vs 1 (error +2); R965C c.2893C>T carriers 3 vs 1 (error +2); S1672Y c.5015C>A carriers 2 vs 1 (error +1); T1709M c.5126C>T carriers 2 vs 1 (error +1); T220I c.659C>T carriers 2 vs 1 (error +1); T632M c.1895C>T carriers 2 vs 1 (error +1); V1405M c.4213G>A carriers 2 vs 1 (error +1); V232I c.694G>A carriers 2 vs 1 (error +1); c.1890+5G>A carriers 2 vs 1 (error +1); c.3840+1G>A carriers 6 vs 1 (error +5); c.4437+5G>A carriers 2 vs 1 (error +1); p.Arg225Trp c.673C>T carriers 3 vs 2 (error +1); p.Arg535* c.1603C>T carriers 4 vs 1 (error +3); p.Glu161Lys c.481G>A carriers 3 vs 1 (error +2); p.Gly1408Arg c.4222G>A carriers 7 vs 1 (error +6); p.Gly752Arg c.2254G>A carriers 5 vs 1 (error +4)

## Scope, method, and limitations

- Population: fixed manifest `paper_manifest.tsv` (12 papers); per-gene counts {'SCN5A': 5, 'RYR2': 5, 'MYBPC3': 2}; every PMID has downloaded source and at least one gold assertion in each count field.
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
