# Codex extraction-blinded paper evaluation — `20260726_fixed48_production`

## Technical summary

This hash-locked run evaluated **48 papers** (**12 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **53.3%**, recall **66.9%**, F1 **59.3%** (670 TP, 588 FP, 331 FN).
- Exact API telemetry: **0 total tokens** (0 input; 0 output).
- Elapsed: **0.0s wall clock**; 0.0s summed per-paper route + read time.
- Representation choices: {'text': 48}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `highcarrier_48_papers_20260723.tsv` (48 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- `selection.json` contains source metadata and hashes but no gold values or gold row counts. `predictions.json` was made read-only and SHA-256 locked before scoring first opened the gold CSVs.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 328 / 1001 | 32.8% | 1.424 | 12.891 |
| affected | 293 / 1001 | 29.3% | 0.137 | 1.993 |
| unaffected | 232 / 1001 | 23.2% | 0.358 | 4.663 |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 180 | 75 | 40 | 70.6% | 81.8% | 75.8% | 64.1% / 0.979 / 6.911 | 61.8% / 0.022 / 0.149 | 50.9% / 0.107 / 1.134 |
| KCNH2 | 229 | 234 | 63 | 49.5% | 78.4% | 60.7% | 25.3% / 1.324 / 11.162 | 20.2% / 0.000 / 0.000 | 20.5% / 0.000 / 0.000 |
| KCNQ1 | 116 | 174 | 155 | 40.0% | 42.8% | 41.4% | 5.5% / 15.267 / 50.692 | 3.7% / 3.400 / 10.752 | 0.4% / 70.000 / 70.000 |
| RYR2 | 145 | 105 | 73 | 58.0% | 66.5% | 62.0% | 45.0% / 0.020 / 0.202 | 40.4% / 0.034 / 0.238 | 27.1% / 0.017 / 0.130 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 27566755 | text | 51 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 0.0 | 0 |
| SCN5A | 26669661 | text | 26 | 28 | 1 | 48.1% | 96.3% | 64.2% | 3.7% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 20470418 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 59.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 55.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 29709101 | text | 10 | 12 | 1 | 45.5% | 90.9% | 60.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 28341781 | text | 52 | 1 | 3 | 98.1% | 94.5% | 96.3% | 94.5% / 0.000 | 94.5% / 0.000 | 94.5% / 0.000 | 0.0 | 0 |
| SCN5A | 18451998 | text | 2 | 1 | 0 | 66.7% | 100.0% | 80.0% | 50.0% / 9.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 26921764 | text | 24 | 3 | 3 | 88.9% | 88.9% | 88.9% | 88.9% / 0.083 | 88.9% / 0.083 | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 27554632 | text | 9 | 10 | 0 | 47.4% | 100.0% | 64.3% | 88.9% / 0.125 | 77.8% / 0.143 | 77.8% / 0.000 | 0.0 | 0 |
| SCN5A | 10590249 | text | 0 | 2 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 25051102 | text | 2 | 16 | 1 | 11.1% | 66.7% | 19.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 26746457 | text | 2 | 1 | 30 | 66.7% | 6.2% | 11.4% | 6.2% / 6.000 | 6.2% / 0.000 | 6.2% / 6.000 | 0.0 | 0 |
| KCNH2 | 29622001 | text | 33 | 23 | 2 | 58.9% | 94.3% | 72.5% | 51.4% / 0.000 | 14.3% / 0.000 | 14.3% / 0.000 | 0.0 | 0 |
| KCNH2 | 11854117 | text | 37 | 0 | 7 | 100.0% | 84.1% | 91.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 14661677 | text | 0 | 0 | 29 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 19160088 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 66.7% / 48.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 26496715 | text | 51 | 1 | 3 | 98.1% | 94.4% | 96.2% | 94.4% / 0.000 | 94.4% / 0.000 | 94.4% / 0.000 | 0.0 | 0 |
| KCNH2 | 11844290 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 10973849 | text | 55 | 192 | 5 | 22.3% | 91.7% | 35.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 10862094 | text | 6 | 6 | 2 | 50.0% | 75.0% | 60.0% | 25.0% / 0.000 | 25.0% / 0.000 | 25.0% / 0.000 | 0.0 | 0 |
| KCNH2 | 10841244 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 100.0% / 0.000 | 100.0% / 0.000 | 0.0 | 0 |
| KCNH2 | 23864605 | text | 2 | 5 | 0 | 28.6% | 100.0% | 44.4% | 0.0% / n/a | 0.0% / n/a | 50.0% / 0.000 | 0.0 | 0 |
| KCNH2 | 24667783 | text | 17 | 0 | 6 | 100.0% | 73.9% | 85.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 19038855 | text | 25 | 7 | 3 | 78.1% | 89.3% | 83.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 19490272 | text | 54 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 23153844 | text | 21 | 12 | 0 | 63.6% | 100.0% | 77.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 17470695 | text | 0 | 1 | 56 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 14678125 | text | 0 | 0 | 41 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 28720088 | text | 1 | 5 | 1 | 16.7% | 50.0% | 25.0% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 21129503 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 25087618 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 17192539 | text | 1 | 0 | 56 | 100.0% | 1.8% | 3.4% | 1.8% / 193.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 24052033 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 18713323 | text | 6 | 2 | 0 | 75.0% | 100.0% | 85.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 29197658 | text | 29 | 7 | 0 | 80.6% | 100.0% | 89.2% | 31.0% / 0.000 | 31.0% / 0.000 | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 33141630 | text | 1 | 146 | 0 | 0.7% | 100.0% | 1.4% | 100.0% / 36.000 | 100.0% / 34.000 | 100.0% / 70.000 | 0.0 | 0 |
| RYR2 | 25814417 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 29925740 | text | 49 | 1 | 2 | 98.0% | 96.1% | 97.0% | 96.1% / 0.000 | 96.1% / 0.000 | 96.1% / 0.000 | 0.0 | 0 |
| RYR2 | 33315912 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 16272262 | text | 12 | 37 | 0 | 24.5% | 100.0% | 39.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 34202968 | text | 1 | 9 | 0 | 10.0% | 100.0% | 18.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 19926015 | text | 0 | 0 | 40 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 33686871 | text | 1 | 1 | 7 | 50.0% | 12.5% | 20.0% | 0.0% / n/a | 12.5% / 2.000 | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 28237968 | text | 1 | 22 | 17 | 4.3% | 5.6% | 4.9% | 5.6% / 2.000 | 5.6% / 1.000 | 5.6% / 1.000 | 0.0 | 0 |
| RYR2 | 12106942 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 22677073 | text | 23 | 9 | 3 | 71.9% | 88.5% | 79.3% | 76.9% / 0.000 | 34.6% / 0.000 | 34.6% / 0.000 | 0.0 | 0 |
| RYR2 | 30403697 | text | 20 | 21 | 1 | 48.8% | 95.2% | 64.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 33606749 | text | 30 | 2 | 3 | 93.8% | 90.9% | 92.3% | 84.8% / 0.000 | 84.8% / 0.000 | 0.0% / n/a | 0.0 | 0 |

## Errors and representation choices

### SCN5A PMID 27566755

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### SCN5A PMID 26669661

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P.Y1795_E1796INSD
- Extra predictions: A344sp, P.TYR1795_GLU1796INS, c.1231G>A, c.1844G>A, c.2102C>T, c.2390G>T, c.2441G>A, c.2923C>T, c.3094G>A, c.3098A>G, c.3236C>T, c.3285G>T, c.3342C>A, c.3556G>A, c.3662C>T, c.3833T>A, c.3911C>T, c.3988G>A, c.4519_4527delCAGAAGCCC, c.4850_4852del, c.4859C>T, c.4931G>A, c.5236G>A, c.5302A>G, c.5350G>A, c.5369A>G, c.5385_5387dupTGA, c.53G>A

### SCN5A PMID 20470418

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: S1102Y
- Count disagreements: p.Ser1103Tyr c.3308C>A carriers 26 vs 85 (error -59)

### SCN5A PMID 28339995

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.D1790G carriers 30 vs 85 (error -55)

### SCN5A PMID 29709101

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.4813+3_4813_6dup
- Extra predictions: F756L, F756fsX, I1660V, c.1045G>A, c.1100G>A, c.1127G>A, c.1312A>T, c.2465G>A, c.2466G>A, c.2658T>A, c.4895G>A, c.4978A>G

### SCN5A PMID 28341781

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.3840+1G>A, c.4245+1G>A, c.4299+1delG
- Extra predictions: IVS21+1G>A

### SCN5A PMID 18451998

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Val1098Leu
- Count disagreements: E1784K carriers 41 vs 50 (error -9)

### SCN5A PMID 26921764

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: L1373X, P.Y1795_E1796INSD, R535X
- Extra predictions: c.5387_5388insTGA, p.L1373* c.4118T>A, p.R535* c.1603C>T
- Count disagreements: p.M369K c.1106T>A carriers 2 vs 3 (error -1); p.R1193Q c.3578G>A carriers 2 vs 3 (error -1); p.M369K c.1106T>A affected 2 vs 3 (error -1); p.R1193Q c.3578G>A affected 2 vs 3 (error -1)

### SCN5A PMID 27554632

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1180V, D1275N, T750N, V1279I, p.D1182D c.3546C>T, p.E547E c.1641A>G, p.H558R c.1673A>G, p.L561L c.1681C>T, p.S525S c.1575C>T, p.Y1261Y c.3783C>T
- Count disagreements: p.T1247I c.3740C>T carriers 18 vs 19 (error -1); p.T1247I c.3740C>T affected 18 vs 19 (error -1)

### SCN5A PMID 10590249

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P.Y1795_E1796INSD
- Extra predictions: D1790G, E1784K

### SCN5A PMID 25051102

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.2788-6C>T
- Extra predictions: G400A, IVS16-6C>T, Q371X, R1193E, R680H, S1103Y, S216L, V168M, c.1199G>C, c.287C>T, c.3183G>A, c.4509C>T, c.5457T>A, c.630G>A, c.87G>A, p.His558Arg c.1673A>G

### SCN5A PMID 26746457

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A1180V, A1270S, A1680T, A551V, C982R, D1275N, E1053K, E48K, F1596I, F532C, G1935S, G552W, G615E, I1836T, I848F, L1194M, L1704H, P717L, R1023H, R1195S, R1512W, R1739Q, R1898C, R2012C, R282C, R689H, S1904L, T1304M, V1353M, V1532I
- Extra predictions: T983I
- Count disagreements: p.Arg1193Gln carriers 19 vs 7 (error +12); p.Arg1193Gln unaffected 12 vs 0 (error +12)

### KCNH2 PMID 29622001

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.1129-58_1320Del, P1034fsX
- Extra predictions: A341V, A525T, D202H, D317N, G269S, G325R, G589D, R366W, R518X, R561G, R594Q, S277DEL, S277L, S546L, T169K, T311I, W248C, Y171X, Y315C, c.1032G>A, c.1129-2A>G, c.3093_3106del, c.683+5G>A

### KCNH2 PMID 11854117

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A83fsX, G925fsX, I593X, L799SP, P968fsX, Q376SP, V295fsX

### KCNH2 PMID 14661677

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A1058E, A190T, A203T, A915V, C723R, G187del, G187S, G873S, H254Q, K897T, K897T, K897T, K897T, L1023del, N257H, N33T, P251A, P347S, P910L, P917L, P967L, Q1068R, R1035W, R1047L, R1047L, R176W, R181Q, T367S, V215G

### KCNH2 PMID 19160088

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R176W
- Count disagreements: p.Arg176Trp carriers 16 vs 112 (error -96)

### KCNH2 PMID 26496715

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P191fsX, T443fsX, A797T
- Extra predictions: p.Arg797Thr

### KCNH2 PMID 11844290

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R752W, F805C, M124R, V822M, W1001X

### KCNH2 PMID 10973849

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: D864SP, K638del, L799SP, Q376SP, Q376SP
- Extra predictions: A150fsX, A172C, A178P, A178T, A194P, A209G, A300T, A332G, A341E, A341V, A344V, A371T, A525T, A551C, A558P, A614V, A842G, A944C, A944G, A98C, C1022A, C1022T, C1024T, C1031T, C1039T, C1046G, C1066T, C1096T, C132A, C20T, C215A, C221T, C241T, C520T, C674T, C727T, C817T, C87A, C926G, C932T, C934T, C935T, C939G, C958G, D1114N, D242N, D317N, D76N, E160K, E261K, E543fsX, F157C, F167W, F640L, G1032A, G1032C, G1033C, G1034A, G1036fsX, G1097A, G1097C, G1111A, G1128A, G139T, G140T, G157C, G167A, G168R, G179S, G189fsX, G226A, G232C, G269D, G269S, G306R, G314S, G325R, G345E, G345R, G478A, G502A, G521A, G532A, G532C, G535A, G569A, G572C, G572R, G580C, G601S, G724A, G728A, G760A, G781A, G805A, G806A, G898A, G914C, G916A, G928A, G940A, G949A, G954C, G95A, G973A, I313M, I593G, K318N, L191fsX, L200fsX, L250H, L266P, L273F, L342F, L353P, L51H, L650fsX, N470D, N588D, N629D, N629K, N629S, N633S, P320A, P448R, P630fsX, P631fsX, Q356X, Q530X, Q725X, R174C, R174H, R190Q, R243C, R243H, R32H, R366P, R366Q, R366W, R518X, R534C, R539W, R555C, R582C, R583C, R591H, R594Q, S225L, S349W, S373P, S428X, S566F, S74L, S818L, T1058C, T1117C, T196G, T257G, T259C, T309R, T311I, T312I, T391I, T436M, T470G, T474I, T556fsX, T587M, T742C, T749A, T797C, V254M, V310I, V47F, V612L, V630A, V630L, V796fsX, V822M, W248R, W305S, W392R, W87R, Y111C, Y184S, Y281C, Y315C, Y315S, Y420fsX, Y493X, Y611H, Y611X

### KCNH2 PMID 10862094

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P151fsX, S543fsX
- Extra predictions: K897K, T897T, c.1631_1632delAG, c.1631delAG, c.453delC, p.Lys897Thr c.2690A>C

### KCNH2 PMID 10841244

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Leu552Ser c.1655T>C carriers 42 vs 44 (error -2)

### KCNH2 PMID 23864605

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: F805C, N470D, S25N, Y652T, p.Asn629Asp

### KCNH2 PMID 24667783

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A172V, E807X, L343fsX, S855R, S890C, W705fsX

### KCNH2 PMID 19038855

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: E698X, I560fsX, I642del
- Extra predictions: E689X c.2062C>T, F640V c.1918T>G, I642_V644del c.1926_1934del, N633S c.1898A>G, Q376sp c.1128G>A, Q901fs/71 c.2705delC, S649P c.1945T>C

### KCNQ1 PMID 19490272

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 23153844

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: P400FS, p.A341E, p.G314S, p.M159*, p.Q357R, p.Q530*, p.R243C, p.S225L, p.S349W, p.T312I, p.W305S, p.Y171*

### KCNQ1 PMID 17470695

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: M1V, T144A, A150FS/133[DELCT451-452], E160K, G168R, Y171X[513 C>G], R174H, A178P, Y184S, G185S, G189E, R190Q, L191FS/90[DELTGCGC572-576], R195FS/40[DELG585], S225L, A226V, R237P, D242N, R243C, V254M, R258C, R259C, L266P, G269S, L273F, I274V, S277L, G292D, F296S, G306R, T312I, G314S, Y315C, Y315S, P320H, T322M, G325R, DELF340[DELCTT1017-1019], A341E, A341V, P343S, S349W, S373P, P400FS/62[INSC1201-1202], I517T, R518X[1552C>T], M520R, V524G, Q530X[1588C>T], R562M, S566F, S571FS/20[DELC1714), R591H, R594Q, D611Y, A636FS/28[DELC1909]
- Extra predictions: A344A

### KCNQ1 PMID 14678125

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: T144A, L151, G168R, Y171X, V172M, A178P, G185S, R190Q, S225L, R243C, V254L, V254M, L266P, G269D, G269S, L273F, E284K, A300T, G306R, T311I, T312I, Y315C, D317G, DELF340, A341E, A341V, A344/SP, A344A/SPLICE, G345R, S349W, Q357H, R366Q, R366W, K393N, R518X, V524G, Q530X, R539W, S566F, I567S, R594Q

### KCNQ1 PMID 28720088

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R519X
- Extra predictions: A341V, G589D, L353L, R518A, p.Arg518*

### KCNQ1 PMID 21129503

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A341V

### KCNQ1 PMID 25087618

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A341V

### KCNQ1 PMID 17192539

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P117L, G168R, R174C, R174H, L175fsX, A178T, Y184S, R190W, R190Q, R190L, S225L, R231C, R243C, V254M, H258N, H258P, E261K, G269D, S277W, V280A, V280E, V288fsX, G306R, T309R, T311I, T312I, G314S, Y315S, Y315C, G316E, P320H, P320A, G325R, F339S, A341V, L342F, A344V, SP/A344A, Q359-K362DEL, R366W, A371T, S373P, W379S, K422fsX, M476L, M520R, R539W, R555C, R555H, S566P, I567T, T587M, A590T, R591H, R594C, Q604X
- Count disagreements: p.Gly589Asp carriers 243 vs 50 (error +193)

### KCNQ1 PMID 24052033

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### KCNQ1 PMID 18713323

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G643S, P488R

### KCNQ1 PMID 29197658

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: F127L, L619M, P477L, p.Ala2Val, p.Lys398Arg, p.Phe479Leu, p.Ser389Pro

### KCNQ1 PMID 33141630

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A208V, A223A, A226A, A226V, A99V, C214C, C214Y, D202G, D202H, D202N, D75G, D75H, D75N, G216G, G219E, G219G, I100L, I204F, I204I, I204M, I227L, I77F, I77M, K218E, K91E, L203P, L76P, M210I, M210T, M83T, Q220K, S209F, S209P, S209S, S217F, S225L, S225S, S82F, S82P, S90F, S98L, T224T, T97M, V205M, V206L, V206V, V207L, V207M, V212V, V215G, V215M, V215V, V221V, V78M, V80M, V88G, V88M, c.223G>A, c.223G>C, c.224A>G, c.227T>C, c.229A>T, c.231C>G, c.232G>A, c.238G>A, c.244T>C, c.245C>T, c.248T>C, c.262G>A, c.263T>G, c.269C>T, c.271A>G, c.290C>T, c.293C>T, c.296C>T, c.298A>C, c.604+12C>T, c.604+16C>G, c.604+16C>T, c.604+17C>T, c.604+19G>A, c.604+1G>A, c.604+22G>A, c.604+24G>T, c.604+28A>C, c.604+29C>T, c.604+30G>A, c.604+33C>G, c.604+40C>G, c.604+41C>T, c.604+42A>T, c.604+43G>A, c.604+49G>A, c.604+49G>T, c.604+62C>A, c.604+72C>G, c.604+7C>T, c.604+9T>C, c.604G>A, c.605-10C>T, c.605-11G>A, c.605-16C>T, c.605-18C>A, c.605-18C>G, c.605-18C>T, c.605-24C>T, c.605-28A>G, c.605-30G>A, c.605-31C>T, c.605-33C>A, c.605-36G>A, c.605-37C>T, c.605-38G>A, c.605-41G>A, c.605-46T>G, c.605-48C>A, c.605-52G>T, c.605-65A>G, c.612C>T, c.613G>A, c.616G>C, c.618C>T, c.619G>A, c.619G>C, c.619G>T, c.623C>T, c.627C>T, c.630G>A, c.636C>G, c.641G>A, c.642C>T, c.643G>A, c.644T>G, c.645G>A, c.648C>A, c.648C>T, c.656G>A, c.657G>A, c.658C>A, c.663G>A, c.669C>T, c.672G>A, c.674C>T, c.675G>A, c.678C>G, c.678C>T
- Count disagreements: T224M carriers 124 vs 88 (error +36); T224M affected 0 vs 34 (error -34); T224M unaffected 124 vs 54 (error +70)

### RYR2 PMID 25814417

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: V2475F

### RYR2 PMID 29925740

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: EXON 3 DELETION, N4168S
- Extra predictions: N4178S c.12533A>G

### RYR2 PMID 33315912

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G357S, c.6982C>T

### RYR2 PMID 16272262

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A2387P, A2403T, A2428T, A4510T, A4607P, A4860G, E2311D, E4146K, E4950K, F4499C, G3946S, G4671R, H4833Y, I419F, I4848V, I4867M, L3778F, L433P, M450I, N2386I, N4097S, N4104K, N4895D, P164S, P2328S, Q4201R, R176Q, R2474S, R414L, R420W, R4497C, S2246L, T2504M, T4158P, V4653F, V4880A, Y2392C

### RYR2 PMID 34202968

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A4860G, D3330G, D4646A, G2866DEL, G3037D, G357S, G4897D, N3308S, R3190Q

### RYR2 PMID 19926015

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: L62F, M81L, E243K, F329L, R332W, G357S, V377M, T415R, R420Q, V507I, A549V, R739H, R1013Q, T1107M, A1136V, E1837K, E2045G, Y2156C, H2168Q, E2183V, D2216V, EXON 3 DELETION, E2296Q, R2420W, M3972I, D3973H, S4124G, Y4149S, R4157Q, Q4159P, N4178S, E4187Q, G4315E, K4650E, N4736 DEL, R4790Q, K4805R, R4822H, L3974Q, K3997E

### RYR2 PMID 33686871

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: I3995V, D4112N, T4196I, D4646A, Q4879H, K4594R, I2075T
- Extra predictions: G357S
- Count disagreements: p.G3118R affected 4 vs 6 (error -2)

### RYR2 PMID 28237968

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: T4630C, G4749E, A2387T, D2300A, I3995V, H4579Y, D1872N, E3716Q, I4867T, G2866del, L4670F, A2387V, Q32K, E1083K, M3972I, V2174F, R3190Q
- Extra predictions: c.11147A>G, c.1191G>A, c.11983A>G, c.12470G>A, c.1258C>T, c.13735C>T, c.13890G>A, c.14008C>T, c.14246G>A, c.14553C>A, c.14600T>C, c.3407C>T, c.5170G>A, c.5614G>A, c.5654G>A, c.6520G>T, c.6899A>C, c.7159G>A, c.7160C>T, c.8598del, c.94C>A, c.9569G>A
- Count disagreements: c.13352del carriers 1 vs 3 (error -2); c.13352del affected 1 vs 2 (error -1); c.13352del unaffected 0 vs 1 (error -1)

### RYR2 PMID 12106942

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### RYR2 PMID 22677073

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: V2113M, Q2958R, A1136V
- Extra predictions: A302V, D85N, F2004L, G643S, P2006A, R1047L, V1951L, W379R, c.2398+5G>T

### RYR2 PMID 30403697

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: G4722S
- Extra predictions: A103V, C1274Y, D4642H, D85N, H1042R, I882V, IVS5+1G>C, K296Q, K4742R, L714R, M0781T, Q584R, Q692K, R1458G, R485C, R4954K, S2774G, V288I, c.13781A>G, c.6224T>C, p.G4772S

### RYR2 PMID 33606749

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: N57_G91DEL35, Q3861H, 4944_4945INSH
- Extra predictions: G380R, c.14834_14835insTCA

## Scope, method, and limitations

- Population: fixed manifest `highcarrier_48_papers_20260723.tsv` (48 papers); 12 per cardiac gene; every PMID has downloaded source and at least one gold assertion in each count field.
- Blinding: gold was used only for PMID eligibility and count-field presence during selection; extraction exported no gold values or row counts, and predictions were made read-only and SHA-256 locked before `score` opened gold.
- Variant metrics are micro-averaged over gold rows. Precision treats unmatched predictions as false positives, although the curated recall packet may omit some real variants.
- Count MAE/RMSE are conditional on a supplied value. Count recall must be read alongside them because abstentions and missed variants are excluded from error magnitude.
- Source acquisition and gold completeness are separate from model reading quality; abstract-only or incomplete source is retained and labeled rather than silently excluded.
- The audited notation score is primary; the preserved raw score bounds sensitivity to post-lock matching adjudication.

## Reproducibility and evidence

- `selection.json`: selected PMIDs, source paths, source hashes, and available representations.
- `predictions.json`: immutable per-paper tools, rationales, extracted variants, counts, evidence quotes, source locations, elapsed time, and token telemetry.
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
