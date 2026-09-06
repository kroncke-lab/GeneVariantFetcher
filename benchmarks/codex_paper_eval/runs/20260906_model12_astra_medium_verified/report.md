# Codex extraction-blinded paper evaluation — `20260906_model12_astra_medium_verified`

## Technical summary

This hash-locked run evaluated **12 papers** (**per-gene counts {'SCN5A': 5, 'RYR2': 5, 'MYBPC3': 2}**) after selecting only PMIDs with downloaded source and at least one named, non-excluded gold variant. Codex predictions were finalized before scoring.

- Variant precision **96.7%**, recall **77.2%**, F1 **85.8%** (608 TP, 21 FP, 180 FN).
- Precision versus counted extras **97.6%** (608 matched rows; 15 extra rows with patient counts). The stricter count-bearing-only diagnostic is **97.4%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Complete API token telemetry is unavailable for this run. Traced failures may have unknown usage; any recorded token sum covers known usage only and must not be interpreted as total cost.
- Elapsed: **19999.0s wall clock**; 11715.5s summed per-paper route + read time.
- Notation twins merged before scoring: **1** same-paper prediction rows that were the same variant in another notation (equivalent-allele identity only; ambiguous or count-conflicting rows were left separate).
- Representation choices: {'text': 12}.

## Provenance-separated identity scores

The paper-derived lane is primary. ClinVar/PubTator citation linkage is retained as a secondary enrichment diagnostic and does not count as finding a variant in the paper.

| Lane | Role | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| `paper_derived` | primary | 608 | 21 | 180 | 96.7% | 77.2% | 85.8% |
| `linkage_assisted` | secondary_diagnostic | 608 | 324 | 180 | 65.2% | 77.2% | 70.7% |

## Blinding and scorer audit

- Paper selection used the fixed manifest `paper_manifest.tsv` (12 papers) from the downloaded-source, named-variant-gold-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: gold was used only for PMID eligibility under the recorded `variant` rule; extraction exported no gold identities, values, or row counts, and predictions were locked before `score` opened gold.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 559 / 788 | 70.9% | 0.175 | 0.718 |
| affected | 37 / 788 | 4.7% | 0.000 | 0.000 |
| unaffected | 65 / 787 | 8.3% | 0.000 | 0.000 |

Gold encodes "no such individuals reported" as an explicit 0 while the pipeline deliberately abstains with NULL, so pooled count recall mixes that convention gap with real attribution misses. The stratified view separates them; the non-zero column is the actionable attribution number.

| field | non-zero gold: supplied / asserted | non-zero recall | zero gold: supplied / asserted | zero recall |
|---|---:|---:|---:|---:|
| carriers | 559 / 705 | 79.3% | 0 / 83 | 0.0% |
| affected | 37 / 650 | 5.7% | 0 / 138 | 0.0% |
| unaffected | 65 / 86 | 75.6% | 0 / 701 | 0.0% |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 551 | 13 | 167 | 97.7% | 76.7% | 86.0% | 98.2% | 98.1% | 73.8% / 0.185 / 0.737 | 1.3% / 0.000 / 0.000 | 8.9% / 0.000 / 0.000 |
| RYR2 | 17 | 1 | 9 | 94.4% | 65.4% | 77.3% | 100.0% | 100.0% | 15.4% / 0.000 / 0.000 | 11.5% / 0.000 / 0.000 | 4.0% / 0.000 / 0.000 |
| MYBPC3 | 40 | 7 | 4 | 85.1% | 90.9% | 87.9% | 88.9% | 83.3% | 56.8% / 0.000 / 0.000 | 56.8% / 0.000 / 0.000 | 0.0% / n/a / n/a |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| RYR2 | 18929323 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | n/a | n/a |
| SCN5A | 30059973 | text | 182 | 3 | 3 | 98.4% | 98.4% | 98.4% | 98.4% / 0.000 | 0.0% / n/a | 0.0% / n/a | n/a | n/a |
| MYBPC3 | 20433692 | text | 13 | 0 | 3 | 100.0% | 81.2% | 89.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | n/a | n/a |
| MYBPC3 | 21302287 | text | 27 | 7 | 1 | 79.4% | 96.4% | 87.1% | 89.3% / 0.000 | 89.3% / 0.000 | 0.0% / n/a | n/a | n/a |
| RYR2 | 30403697 | text | 13 | 0 | 8 | 100.0% | 61.9% | 76.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | n/a | n/a |
| SCN5A | 20031634 | text | 10 | 3 | 3 | 76.9% | 76.9% | 76.9% | 69.2% / 0.444 | 69.2% / 0.000 | 69.2% / 0.000 | n/a | n/a |
| RYR2 | 19398417 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | n/a | n/a |
| RYR2 | 25435091 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | n/a | n/a |
| SCN5A | 32533946 | text | 0 | 0 | 83 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | n/a | n/a |
| SCN5A | 25163546 | text | 20 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | n/a | n/a |
| RYR2 | 25814417 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | n/a | n/a |
| SCN5A | 20129283 | text | 339 | 7 | 78 | 98.0% | 81.3% | 88.9% | 81.3% / 0.277 | 0.0% / n/a | 13.2% / 0.000 | n/a | n/a |

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

- Missed gold variants: c.2258_2259insT
- Extra predictions: F244F c.732C>T, K754E, K754EfsX78, L517M c.1549C>A, c.2258-2259InsT, c.2846-2847InsT, c.3192-3193InsC

### RYR2 PMID 30403697

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: F3496L, G1885E, G1886S, G4772S, H2464D, S3938R, c.14091-11dupT, c.3599-9delT

### SCN5A PMID 20031634

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: P.1570_F1571INSI, c.1983_1993dupGGCCCTCAGCG, c.3963+2T>C
- Extra predictions: c.3816delG, p.Ala665GlyfsTer16 c.1983_1993dup, p.Trp156*
- Count disagreements: p.Arg225Trp c.673C>T carriers 11 vs 10 (error +1); p.Asn1722Asp c.5164A>G carriers 9 vs 8 (error +1); p.Gly1408Arg c.4222G>A carriers 14 vs 13 (error +1); p.Ser1382Ile c.4145G>T carriers 9 vs 8 (error +1)

### RYR2 PMID 19398417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 25435091

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Extra predictions: p.Arg420Gln

### SCN5A PMID 32533946

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1428V, A166T, A1680T, A1746V, A178G, A286S, A735E, A735T, C335R, D1243N, D1370G, D349N, D772N, D785N, D84N, E1225K, E1574K, E1784K, E346D, E746K, E901K, F851L, F892I, F93S, G1262S, G1406E, G1420R, G1420V, G1642E, G1661R, G386R, G752R, G833R, G9V, K175N, L1346P, L136P, L276Q, L299M, L839P, L928P, M369K, M734V, M764K, N109K, N1325S, N1380K, N1722D, N927S, P1014S, P1730H, P773S, R121W, R1432G, R1583C, R1632H, R1898C, R1958X, R282C, R367C, R367L, R620H, R808C, R814Q, S1382I, S1672Y, T1461S, T1709M, T220I, T353I, V1251M, V1279I, V1281F, V1353M, V1405L, V1405M, V223L, V396L, V714A, V924I, W1345C, W879R, Y1449C

### SCN5A PMID 25163546

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 25814417

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: G357S

### SCN5A PMID 20129283

**text** — Production gvf-run calibrated-manifest strategy using the current deterministic and model-backed extraction, verification, trust, and recovery layers. Exact code, prompt, model, source, and call-trace provenance is bound separately by the run and evaluation manifests. Mapped to the harness 'text' route because no single harness route describes the multi-stage production pipeline.

- Missed gold variants: A1680T, A735V, A735V, A735V, D1243N, D1243N, D1275N, E1225K, E1225K, E1225K, E161K, E746K, E746K, G1319V, G1408R, G1408R, G1661R, G1743E, G1743R, G1743R, G1743R, G752R, I1660V, L1393X, L1393X, L868X, N927S, P.F2004DUP, P.I137_C139DUP, P.I1570DUP, P.K1493DEL, P.Y1795_E1796INSD, P336L, Q1695X, R104Q, R104Q, R104W, R121Q, R1232W, R1232W, R1583C, R1623X, R1638X, R225W, R367H, R367H, R367H, R376H, R376H, R526H, R535X, R535X, R535X, R878H, R878H, R878H, R878H, R893H, R893H, R965C, R965C, S1672Y, T1620M, T1709M, T220I, T632M, V1405M, V232I, c.1890+5G>A, c.1936delC, c.1936delC, c.2435_2436+3delTGGTAinsCGCCT, c.2549_2550insTG, c.2602delC, c.2602delC, c.3667delG, c.3840+1G>A, c.4437+5G>A
- Extra predictions: A1223PfsX7 c.3666delG, F2004dup c.6010_6012dupTTC, I137_C139dup c.410_418dupTCATGTGCA, I1570dup c.4708_4710dupATC, V1525M c.4573G>A, c.2435_24363delTGGTAinsCGCCT, c.5387_5388insTGA
- Count disagreements: A1680T c.5038G>A carriers 2 vs 1 (error +1); A226V carriers 1 vs 2 (error -1); A735V c.2204C>T carriers 4 vs 1 (error +3); D1275N c.3823G>A carriers 3 vs 1 (error +2); E1225K c.3673G>A carriers 4 vs 1 (error +3); E746K c.2236G>A carriers 3 vs 1 (error +2); G1319V c.3956G>T carriers 5 vs 1 (error +4); G1408R carriers 7 vs 1 (error +6); G1743E c.5228G>A carriers 6 vs 1 (error +5); G1743R c.5227G>A carriers 5 vs 1 (error +4); I1660V carriers 1 vs 4 (error -3); K1493X c.4477A>T carriers 1 vs 2 (error -1); K1493del c.4477_4479delAAG carriers 2 vs 1 (error +1); L1393X c.4178T>A carriers 3 vs 1 (error +2); L868X c.2602delC carriers 2 vs 1 (error +1); N927S c.2780A>G carriers 3 vs 1 (error +2); P336L c.1007C>T carriers 2 vs 1 (error +1); Q1695X c.5083C>T carriers 2 vs 1 (error +1); Q646RfsX5 c.1936delC carriers 3 vs 1 (error +2); R104Q c.311G>A carriers 3 vs 1 (error +2); R104W c.310C>T carriers 2 vs 1 (error +1); R121Q c.362G>A carriers 2 vs 1 (error +1); R1232W c.3694C>T carriers 3 vs 1 (error +2); R1623X c.4867C>T carriers 2 vs 1 (error +1); R1638X c.4912C>T carriers 3 vs 1 (error +2); R367H c.1100G>A carriers 6 vs 1 (error +5); R376H c.1127G>A carriers 4 vs 1 (error +3); R526H c.1577G>A carriers 2 vs 1 (error +1); R878H c.2633G>A carriers 5 vs 1 (error +4); R893H c.2678G>A carriers 3 vs 1 (error +2); R965C c.2893C>T carriers 3 vs 1 (error +2); S1672Y c.5015C>A carriers 2 vs 1 (error +1); T1709M c.5126C>T carriers 2 vs 1 (error +1); T220I c.659C>T carriers 2 vs 1 (error +1); T632M c.1895C>T carriers 2 vs 1 (error +1); V1405M c.4213G>A carriers 2 vs 1 (error +1); V232I c.694G>A carriers 2 vs 1 (error +1); c.1890+5G>A carriers 2 vs 1 (error +1); c.3840+1G>A carriers 6 vs 1 (error +5); c.4437+5G>A carriers 2 vs 1 (error +1); p.Arg225Trp c.673C>T carriers 3 vs 2 (error +1); p.Arg535* c.1603C>T carriers 4 vs 1 (error +3); p.Glu161Lys c.481G>A carriers 3 vs 1 (error +2); p.Gly752Arg c.2254G>A carriers 5 vs 1 (error +4)

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
