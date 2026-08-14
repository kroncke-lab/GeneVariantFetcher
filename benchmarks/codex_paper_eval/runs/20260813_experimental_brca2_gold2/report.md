# Codex extraction-blinded paper evaluation — `20260813_experimental_brca2_gold2`

## Technical summary

This hash-locked run evaluated **2 papers** (**2 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **3.4%**, recall **100.0%**, F1 **6.6%** (7 TP, 198 FP, 0 FN).
- Precision versus counted extras **5.7%** (7 matched rows; 116 extra rows with patient counts). The stricter count-bearing-only diagnostic is **2.5%** and has a different numerator; it is not comparable to the repository's counted-extra precision floor.
- Exact API telemetry: **1,770,488 total tokens** (1,150,331 input; 620,157 output).
- Elapsed: **5869.0s wall clock**; 57.3s summed per-paper route + read time.
- Representation choices: {'text': 2}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `brca2_2_collaborator_reviewed_20260811.tsv` (2 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- Blinding: prediction content was finalized and SHA-256 locked before this external score. The source production workflow may have read registered gold for read-only layer scorecards before the projection lock; those scores did not feed back into extraction. This is not the stricter native lock-before-any-gold-read protocol.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 3 / 7 | 42.9% | 1.333 | 2.309 |
| affected | 0 / 7 | 0.0% | n/a | n/a |
| unaffected | 0 / 7 | 0.0% | n/a | n/a |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| BRCA2 | 7 | 198 | 0 | 3.4% | 100.0% | 6.6% | 5.7% | 2.5% | 42.9% / 1.333 / 2.309 | 0.0% / n/a / n/a | 0.0% / n/a / n/a |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| BRCA2 | 26833046 | text | 4 | 80 | 0 | 4.8% | 100.0% | 9.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 48.3 | 15,111 |
| BRCA2 | 26848529 | text | 3 | 118 | 0 | 2.5% | 100.0% | 4.8% | 100.0% / 1.333 | 0.0% / n/a | 0.0% / n/a | 9.1 | 4,013 |

## Errors and representation choices

### BRCA2 PMID 26833046

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G1738E, c.130T>A, c.1310_1313delAAGA, c.145G>T, c.1556delA, c.1687C>T, c.1803delC, c.1813delA, c.181T>G, c.1877_1878insTAGT, c.190T>C, c.2231C>G, c.2450delA, c.2475delC, c.2523dupG, c.2787dupA, c.2808_2811delACAA, c.2812delGCAA, c.2830A>T, c.2943_2944insTGAGA, c.301+1G>C, c.3052_3053insTGAGA, c.3319G>T, c.3400G>T, c.3477_3479delAAAinsC, c.3599_3600delGT, c.3607C>T, c.3710delT, c.3718C>T, c.3765delT, c.3847_3848delGT, c.3860delA, c.3874delT, c.3904G>T, c.4035delA, c.4258delG, c.427G>T, c.4284_4285insT, c.4285_4286insT, c.469_470delAA, c.5089T>C, c.5096G>A, c.5130_5133delTGTA, c.5143A>C, c.5164_5165delAG, c.5213G>A, c.5219delT, c.5263insC, c.5266dupC, c.5341delG, c.5351delA, c.5352delC, c.5382insC, c.5503C>T, c.5722_5723delCT, c.583delT, c.5857G>T, c.6082_6086delGAAGA, c.6244G>T, c.6408_6414delAAATGTT, c.6443_6444delCT, c.6486_6489delACAA, c.6490_6492del, c.7069_7070delCT, c.7757G>A, c.7878G>C, c.7913_7917delTTCCT, c.7980T>G, c.8364G>A, c.8536G>T, c.8575delC, c.8754+3G>C, c.8953+1G>T, c.9097delA, c.9106C>T, c.7008-1G>C, p.Ala2825fs c.8474del, p.E1107*, p.Lys157fs c.469_470del, p.Thr2125fs c.6373del

### BRCA2 PMID 26848529

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.212+1G>T, c.4986+1G>A, c.7435+53C>T, c.7976+1G>A, p.Ala634Valfs*18 c.1901delC, p.Ala938Profs*21 c.2808_2811delACAA, p.Arg1203* c.3607C>T, p.Arg155Lysfs*26 c.464_468delGAGAT, p.Arg1835* c.5503C>T, p.Arg2108Cys c.6322C>T, p.Arg2336His c.7007G>A, p.Arg435Lysfs*17 c.1303dupA, p.Arg691Serfs*10 c.2073delA, p.Asn1121Lysfs*10 c.3363_3367delTACAG, p.Asn1287Ilefs*6 c.3860delA, p.Asn1355Lysfs*10 c.4065_4068delTCAA, p.Asn1473Glnfs*8 c.4416_4417delGA, p.Asn1822Ilefs*18 c.5465delA, p.Asn2051* c.6150_6151insT, p.Asn2101_Val2102del c.6301_6306delAATGTA, p.Asn361Metfs*6 c.1082delA, p.Asn704Cysfs*7 c.2110_2111delAA, p.Asn863Lysfs*18 c.2588dupA, p.Asn991Asp c.2971A>G, p.Asp1151Metfs*4 c.3450delT, p.Asp687* c.2059_2063delGATTA, p.Asp693Asn c.2077G>A, p.Asp974Asn c.2920G>A, p.Cys328* c.981_982delAT, p.Gln1056Argfs*3 c.3167_3170delAAAA, p.Gln1144His c.3432G>C, p.Gln1291* c.3871C>T, p.Gln2655Asnfs*2 c.7963delC, p.Gln2829* c.8485C>T, p.Gln2941Leufs*34 c.8820_8823del, p.Gln3034Serfs*10 c.9098_9099insA, p.Gln3036Serfs*8 c.9105dupT, p.Gln759* c.2275C>T, p.Gln858* c.2572C>T, p.Glu1038Gly c.3113A>G, p.Glu1214* c.3640G>T, p.Glu1257Glyfs*9 c.3770_3771delAG, p.Glu2198Asnfs*4 c.6591_6592delTG, p.Glu3096* c.9286G>T, p.Glu3263Argfs*12 c.9788delA, p.Glu554* c.1660G>T, p.Gly1788Val c.5363G>T, p.Gly2901Asp c.8702G>A, p.Gly3003* c.9007G>T, p.Gly933Alafs*4 c.2798_2799delGT, p.His2021Pro, p.His372Asn c.1114C>A, p.Ile1170Phe c.3508A>T, p.Ile1318Leu c.3952A>C, p.Ile1724Lysfs*17 c.5171delT, p.Ile1824Aspfs*3 c.5470_5477delATTGGGCA, p.Ile2986Lysfs*3 c.8956_8957insAA, p.Ile332Lysfs*18 c.993_994dupA, p.Ile605Tyrfs*9 c.1813delA, p.Ile770Phefs*2 c.2307delT, p.Ile980Lys c.2939T>A, p.Leu1306Aspfs*23 c.3916_3917delTT, p.Leu1908Argfs*2 c.5722_5723delCT, p.Leu2080* c.6239T>G, p.Leu3119* c.9356_9357delTAinsG, p.Leu750Valfs*10 c.2248_2252delCTCAT, p.Leu88Alafs*12 c.262_263delCT, p.Lys1183Arg c.3548A>G, p.Lys2150Asnfs*19 c.6449_6450insTA, p.Lys307Glu c.919A>G, p.Lys3326* c.9976A>T, p.Lys355Arg c.1064A>G, p.Lys467* c.1399A>T, p.Lys519Argfs*13 c.1556delA, p.Lys585Arg c.1754A>G, p.Met2393Lysfs*18 c.7178_7179delTG, p.Met815Trpfs*10 c.2442delC, p.Phe1177Leufs*33 c.3531delT, p.Phe2801Leufs*10, p.Pro1099Leufs*10 c.3294delT, p.Pro1702Thrfs*16 c.5103dupA, p.Pro2381Hisfs*13 c.7142delC, p.Pro628Hisfs*16 c.1881delA, p.Pro981Ala c.2941C>G, p.Pro999Leu c.2996C>T, p.Ser1040Asn c.3119G>A, p.Ser1041* c.3122C>G, p.Ser1248Argfs*10 c.3744_3747delTGAG, p.Ser1613Gly c.4837A>G, p.Ser1722Tyrfs*4 c.5164_5165delAG, p.Ser1841Valfs*2 c.5521delA, p.Ser1882* c.5645C>A, p.Ser1900* c.5699C>G, p.Ser1955* c.5864C>G, p.Ser2012Glnfs*5 c.6033_6034delTT, p.Ser2120* c.6359C>G, p.Ser2267* c.6800C>A, p.Ser2984Glnfs*4 c.8950delT, p.Ser611* c.1832C>A, p.Ser645Tyr c.1934C>A, p.Ser713* c.2138C>G, p.Ser868* c.2603C>A, p.Thr2471Hisfs*4 c.7409dupT, p.Thr2746Aspfs*18 c.8234dupT, p.Thr2867Ser c.8599A>T, p.Thr3085Glnfs*19 c.9253delA, p.Thr598Ala c.1792A>G, p.Trp1718* c.5154G>A, p.Trp2725Glyfs*8 c.8172delG, p.Trp372* c.1116G>A, p.Tyr1894* c.5681dupA, p.Tyr2215*, p.Tyr2839* c.8517C>A, p.Val1120Aspfs*11 c.3359_3363delTTAAT, p.Val1145Phefs*10 c.3433delG, p.Val220Ilefs*4 c.658_659delGT, c.5576_5579delTTAA, p.Ile1859fs c.5576_5579del
- Count disagreements: c.5576_5580delTTAAA carriers 1 vs 5 (error -4)

## Scope, method, and limitations

- Population: fixed manifest `brca2_2_collaborator_reviewed_20260811.tsv` (2 papers); 2 per cardiac gene; every PMID has downloaded source and at least one gold assertion in each count field.
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
