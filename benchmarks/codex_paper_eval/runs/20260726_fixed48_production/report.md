# Codex extraction-blinded paper evaluation — `20260726_fixed48_production`

## Technical summary

This hash-locked run evaluated **48 papers** (**12 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **44.5%**, recall **78.8%**, F1 **56.9%** (789 TP, 985 FP, 212 FN).
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
| SCN5A | 184 | 151 | 36 | 54.9% | 83.6% | 66.3% | 64.1% / 0.979 / 6.911 | 61.8% / 0.022 / 0.149 | 50.9% / 0.107 / 1.134 |
| KCNH2 | 257 | 492 | 35 | 34.3% | 88.0% | 49.4% | 25.3% / 1.324 / 11.162 | 20.2% / 0.000 / 0.000 | 20.5% / 0.000 / 0.000 |
| KCNQ1 | 200 | 226 | 71 | 46.9% | 73.8% | 57.4% | 5.5% / 15.267 / 50.692 | 3.7% / 3.400 / 10.752 | 0.4% / 70.000 / 70.000 |
| RYR2 | 148 | 116 | 70 | 56.1% | 67.9% | 61.4% | 45.0% / 0.020 / 0.202 | 40.4% / 0.034 / 0.238 | 27.1% / 0.017 / 0.130 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 27566755 | text | 51 | 11 | 0 | 82.3% | 100.0% | 90.3% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 0.0 | 0 |
| SCN5A | 26669661 | text | 26 | 33 | 1 | 44.1% | 96.3% | 60.5% | 3.7% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 20470418 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 59.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 55.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 29709101 | text | 10 | 17 | 1 | 37.0% | 90.9% | 52.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 28341781 | text | 52 | 26 | 3 | 66.7% | 94.5% | 78.2% | 94.5% / 0.000 | 94.5% / 0.000 | 94.5% / 0.000 | 0.0 | 0 |
| SCN5A | 18451998 | text | 2 | 4 | 0 | 33.3% | 100.0% | 50.0% | 50.0% / 9.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 26921764 | text | 24 | 3 | 3 | 88.9% | 88.9% | 88.9% | 88.9% / 0.083 | 88.9% / 0.083 | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 27554632 | text | 9 | 14 | 0 | 39.1% | 100.0% | 56.2% | 88.9% / 0.125 | 77.8% / 0.143 | 77.8% / 0.000 | 0.0 | 0 |
| SCN5A | 10590249 | text | 0 | 3 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 25051102 | text | 2 | 20 | 1 | 9.1% | 66.7% | 16.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| SCN5A | 26746457 | text | 6 | 19 | 26 | 24.0% | 18.8% | 21.1% | 6.2% / 6.000 | 6.2% / 0.000 | 6.2% / 6.000 | 0.0 | 0 |
| KCNH2 | 29622001 | text | 33 | 25 | 2 | 56.9% | 94.3% | 71.0% | 51.4% / 0.000 | 14.3% / 0.000 | 14.3% / 0.000 | 0.0 | 0 |
| KCNH2 | 11854117 | text | 37 | 13 | 7 | 74.0% | 84.1% | 78.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 14661677 | text | 23 | 2 | 6 | 92.0% | 79.3% | 85.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 19160088 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 66.7% / 48.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 26496715 | text | 52 | 5 | 2 | 91.2% | 96.3% | 93.7% | 94.4% / 0.000 | 94.4% / 0.000 | 94.4% / 0.000 | 0.0 | 0 |
| KCNH2 | 11844290 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 10973849 | text | 56 | 424 | 4 | 11.7% | 93.3% | 20.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 10862094 | text | 6 | 10 | 2 | 37.5% | 75.0% | 50.0% | 25.0% / 0.000 | 25.0% / 0.000 | 25.0% / 0.000 | 0.0 | 0 |
| KCNH2 | 10841244 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 100.0% / 0.000 | 100.0% / 0.000 | 0.0 | 0 |
| KCNH2 | 23864605 | text | 2 | 5 | 0 | 28.6% | 100.0% | 44.4% | 0.0% / n/a | 0.0% / n/a | 50.0% / 0.000 | 0.0 | 0 |
| KCNH2 | 24667783 | text | 20 | 0 | 3 | 100.0% | 87.0% | 93.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNH2 | 19038855 | text | 25 | 8 | 3 | 75.8% | 89.3% | 82.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 19490272 | text | 54 | 7 | 0 | 88.5% | 100.0% | 93.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 23153844 | text | 21 | 12 | 0 | 63.6% | 100.0% | 77.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 17470695 | text | 46 | 35 | 10 | 56.8% | 82.1% | 67.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 14678125 | text | 35 | 9 | 6 | 79.5% | 85.4% | 82.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 28720088 | text | 1 | 5 | 1 | 16.7% | 50.0% | 25.0% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 21129503 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 25087618 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 17192539 | text | 4 | 0 | 53 | 100.0% | 7.0% | 13.1% | 1.8% / 193.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 24052033 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 18713323 | text | 6 | 2 | 0 | 75.0% | 100.0% | 85.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 29197658 | text | 29 | 7 | 0 | 80.6% | 100.0% | 89.2% | 31.0% / 0.000 | 31.0% / 0.000 | 0.0% / n/a | 0.0 | 0 |
| KCNQ1 | 33141630 | text | 1 | 147 | 0 | 0.7% | 100.0% | 1.3% | 100.0% / 36.000 | 100.0% / 34.000 | 100.0% / 70.000 | 0.0 | 0 |
| RYR2 | 25814417 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 29925740 | text | 49 | 1 | 2 | 98.0% | 96.1% | 97.0% | 96.1% / 0.000 | 96.1% / 0.000 | 96.1% / 0.000 | 0.0 | 0 |
| RYR2 | 33315912 | text | 1 | 2 | 0 | 33.3% | 100.0% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 16272262 | text | 12 | 41 | 0 | 22.6% | 100.0% | 36.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 34202968 | text | 1 | 9 | 0 | 10.0% | 100.0% | 18.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 19926015 | text | 0 | 0 | 40 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 33686871 | text | 1 | 1 | 7 | 50.0% | 12.5% | 20.0% | 0.0% / n/a | 12.5% / 2.000 | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 28237968 | text | 3 | 25 | 15 | 10.7% | 16.7% | 13.0% | 5.6% / 2.000 | 5.6% / 1.000 | 5.6% / 1.000 | 0.0 | 0 |
| RYR2 | 12106942 | text | 6 | 1 | 0 | 85.7% | 100.0% | 92.3% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 22677073 | text | 24 | 9 | 2 | 72.7% | 92.3% | 81.4% | 76.9% / 0.000 | 34.6% / 0.000 | 34.6% / 0.000 | 0.0 | 0 |
| RYR2 | 30403697 | text | 20 | 23 | 1 | 46.5% | 95.2% | 62.5% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| RYR2 | 33606749 | text | 30 | 2 | 3 | 93.8% | 90.9% | 92.3% | 84.8% / 0.000 | 84.8% / 0.000 | 0.0% / n/a | 0.0 | 0 |

## Errors and representation choices

### SCN5A PMID 27566755

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Arg1643His c.4928G>A, p.Asn1268Ser c.3803A>G, p.Asn1324Ser c.3971A>G, p.Glu1207Lys c.3619G>A, p.Glu1780Gly c.5339A>G, p.Glu1783Lys c.5347G>A, p.Ile1447Leu c.4339A>C, p.Leu1500Val c.4498C>G, p.Met1765Val c.5293A>G, p.Thr1778Met c.5333C>T, p.Val1666Ile c.4996G>A

### SCN5A PMID 26669661

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P.Y1795_E1796INSD
- Extra predictions: A344sp, P.TYR1795_GLU1796INS, c.1231G>A, c.1844G>A, c.2102C>T, c.2390G>T, c.2441G>A, c.2923C>T, c.3094G>A, c.3098A>G, c.3236C>T, c.3285G>T, c.3342C>A, c.3556G>A, c.3662C>T, c.3833T>A, c.3911C>T, c.3988G>A, c.4519_4527delCAGAAGCCC, c.4850_4852del, c.4859C>T, c.4931G>A, c.5236G>A, c.5302A>G, c.5350G>A, c.5369A>G, c.5385_5387dupTGA, c.53G>A, p.Ala1220Val c.3659C>T, p.Arg53Gln c.158G>A, p.Asp1113Glu c.3339C>A, p.Ser1078Phe c.3233C>T, p.Trp1094Cys c.3282G>T

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
- Extra predictions: F756L, F756fsX, I1660V, c.1045G>A, c.1100G>A, c.1127G>A, c.1312A>T, c.2465G>A, c.2466G>A, c.2658T>A, c.4895G>A, c.4978A>G, c.4813+3_4813+6dup, p.Arg1631His c.4892G>A, p.Asp349Gly c.1046A>G, p.Ile1659Val c.4975A>G, p.Met438Thr c.1313T>C

### SCN5A PMID 28341781

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.3840+1G>A, c.4245+1G>A, c.4299+1delG
- Extra predictions: IVS21+1G>A, c.3960+1G>A, c.3960+1G>T, c.4242+1G>C, c.4242+2T>A, p.Ala1427Ser c.4279G>T, p.Arg1194His c.3581G>A, p.Arg1431Ser c.4293G>C, p.Arg1622Ter c.4864C>T, p.Arg1643Cys c.4927C>T, p.Arg1643His c.4928G>A, p.Arg1912Cys c.5734C>T, p.Arg1918His c.5753G>A, p.Arg376Cys c.1126C>T, p.Arg878His c.2633G>A, p.Arg893Pro c.2678G>C, p.Asn1379del c.4134CAA[1], p.Asp1422_Arg c.4265_4294del, p.Glu1224Gln c.3670G>C, p.Glu1224Lys c.3670G>A, p.Gly1742Arg c.5224G>A, p.Leu1895fs c.5684_5685del, p.Thr1246Ile c.3737C>T, p.Thr1708Met c.5123C>T, p.Trp1094Ter c.3281G>A, p.Val1404Met c.4210G>A

### SCN5A PMID 18451998

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Val1098Leu, p.Glu1783Lys c.5347G>A, p.Thr1303Met c.3908C>T, p.Val1097Leu c.3289G>T
- Count disagreements: E1784K carriers 41 vs 50 (error -9)

### SCN5A PMID 26921764

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: L1373X, P.Y1795_E1796INSD, R535X
- Extra predictions: c.5387_5388insTGA, p.L1373* c.4118T>A, p.R535* c.1603C>T
- Count disagreements: p.M369K c.1106T>A carriers 2 vs 3 (error -1); p.R1193Q c.3578G>A carriers 2 vs 3 (error -1); p.M369K c.1106T>A affected 2 vs 3 (error -1); p.R1193Q c.3578G>A affected 2 vs 3 (error -1)

### SCN5A PMID 27554632

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A1180V, D1275N, T750N, V1279I, p.D1182D c.3546C>T, p.E547E c.1641A>G, p.H558R c.1673A>G, p.L561L c.1681C>T, p.S525S c.1575C>T, p.Y1261Y c.3783C>T, p.Gly1261Cys c.3781G>T, p.Gly1261Ser c.3781G>A, p.Met1244Ile c.3732G>A, p.Thr1246Ile c.3737C>T
- Count disagreements: p.T1247I c.3740C>T carriers 18 vs 19 (error -1); p.T1247I c.3740C>T affected 18 vs 19 (error -1)

### SCN5A PMID 10590249

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P.Y1795_E1796INSD
- Extra predictions: D1790G, E1784K, p.Tyr1794_Glu c.5382_5384dup

### SCN5A PMID 25051102

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.2788-6C>T
- Extra predictions: G400A, IVS16-6C>T, Q371X, R1193E, R680H, S1103Y, S216L, V168M, c.1199G>C, c.287C>T, c.3183G>A, c.4509C>T, c.5457T>A, c.630G>A, c.87G>A, p.His558Arg c.1673A>G, c.IVS16-6C>T, p.Arg1192Gln c.3575G>A, p.Phe1595Ile c.4783T>A, p.Ser1502= c.4506C>T

### SCN5A PMID 26746457

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A1180V, A1270S, A1680T, A551V, D1275N, E1053K, E48K, F1596I, F532C, G1935S, G615E, I1836T, I848F, L1194M, L1704H, R1023H, R1195S, R1512W, R1739Q, R1898C, R2012C, R689H, S1904L, T1304M, V1353M, V1532I
- Extra predictions: T983I, p.Ala997Thr c.2989G>A, p.Arg1023Cys c.3067C>T, p.Arg1023Pro c.3068G>C, p.Arg1027Gln c.3080G>A, p.Arg1115Trp c.3343C>T, p.Arg1738Gln c.5213G>A, p.Arg1928His c.5783G>A, p.Arg535Gln c.1604G>A, p.Arg975Gln c.2924G>A, p.Asp1242Asn c.3724G>A, p.Ile1835Thr c.5504T>C, p.Leu1193Met c.3577T>A, p.Leu1703His c.5108T>A, p.Leu299Met c.895T>A, p.Pro2005Ala c.6013C>G, p.Ser1903Leu c.5708C>T, p.Thr1016Met c.3047C>T, p.Thr1303Met c.3908C>T
- Count disagreements: p.Arg1193Gln carriers 19 vs 7 (error +12); p.Arg1193Gln unaffected 12 vs 0 (error +12)

### KCNH2 PMID 29622001

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.1129-58_1320Del, P1034fsX
- Extra predictions: A341V, A525T, D202H, D317N, G269S, G325R, G589D, R366W, R518X, R561G, R594Q, S277DEL, S277L, S546L, T169K, T311I, W248C, Y171X, Y315C, c.1032G>A, c.1129-2A>G, c.3093_3106del, c.683+5G>A, p.Ala558Val c.1673C>T, p.Thr152fs c.453del

### KCNH2 PMID 11854117

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A83fsX, G925fsX, I593X, L799SP, P968fsX, Q376SP, V295fsX
- Extra predictions: p.Arg56Leu c.167G>T, p.Arg823Gln c.2468G>A, p.Arg823Leu c.2468G>T, p.Asn629Ile c.1886A>T, p.Asp609Glu c.1827C>A, p.Asp609Gly c.1826A>G, p.Cys66Arg c.196T>C, p.Cys66Trp c.198C>G, p.Cys66Tyr c.197G>A, p.Gln376= c.1128G>A, p.Gly604Cys c.1810G>T, p.Phe805Ser c.2414T>C, p.Pro72Arg c.215C>G

### KCNH2 PMID 14661677

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: G187del, K897T, K897T, K897T, L1023del, R1047L
- Extra predictions: c.551GCGCGGGCG[1], p.Ala913Val c.2738C>T

### KCNH2 PMID 19160088

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R176W
- Count disagreements: p.Arg176Trp carriers 16 vs 112 (error -96)

### KCNH2 PMID 26496715

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: T443fsX, A797T
- Extra predictions: p.Arg797Thr, c.1557+1G>C, p.Gly71Arg c.211G>A, p.Gly71Glu c.212G>A, p.Pro72Gln c.215C>A

### KCNH2 PMID 11844290

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R752W, F805C, M124R, V822M, W1001X

### KCNH2 PMID 10973849

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: D864SP, L799SP, Q376SP, Q376SP
- Extra predictions: A150fsX, A172C, A178P, A178T, A194P, A209G, A300T, A332G, A341E, A341V, A344V, A371T, A525T, A551C, A558P, A614V, A842G, A944C, A944G, A98C, C1022A, C1022T, C1024T, C1031T, C1039T, C1046G, C1066T, C1096T, C132A, C20T, C215A, C221T, C241T, C520T, C674T, C727T, C817T, C87A, C926G, C932T, C934T, C935T, C939G, C958G, D1114N, D242N, D317N, D76N, E160K, E261K, E543fsX, F157C, F167W, F640L, G1032A, G1032C, G1033C, G1034A, G1036fsX, G1097A, G1097C, G1111A, G1128A, G139T, G140T, G157C, G167A, G168R, G179S, G189fsX, G226A, G232C, G269D, G269S, G306R, G314S, G325R, G345E, G345R, G478A, G502A, G521A, G532A, G532C, G535A, G569A, G572C, G572R, G580C, G601S, G724A, G728A, G760A, G781A, G805A, G806A, G898A, G914C, G916A, G928A, G940A, G949A, G954C, G95A, G973A, I313M, I593G, K318N, L191fsX, L200fsX, L250H, L266P, L273F, L342F, L353P, L51H, L650fsX, N470D, N588D, N629D, N629K, N629S, N633S, P320A, P448R, P630fsX, P631fsX, Q356X, Q530X, Q725X, R174C, R174H, R190Q, R243C, R243H, R32H, R366P, R366Q, R366W, R518X, R534C, R539W, R555C, R582C, R583C, R591H, R594Q, S225L, S349W, S373P, S428X, S566F, S74L, S818L, T1058C, T1117C, T196G, T257G, T259C, T309R, T311I, T312I, T391I, T436M, T470G, T474I, T556fsX, T587M, T742C, T749A, T797C, V254M, V310I, V47F, V612L, V630A, V630L, V796fsX, V822M, W248R, W305S, W392R, W87R, Y111C, Y184S, Y281C, Y315C, Y315S, Y420fsX, Y493X, Y611H, Y611X, c.1128+1G>A, c.1129-1G>A, c.1129-2A>G, c.1557+1G>C, c.1558-1G>C, c.1720_1945+102del, c.1945+1G>A, c.1946-1G>C, c.1946-2A>C, c.1946-4_1948del, c.2145+1G>A, c.2146-1G>C, c.2146-2A>G, c.2146-50_2167del, c.2398+1G>T, c.2398_2398+21del, c.2399-108_2488del, c.2399-2A>G, c.2592+1G>A, c.2593-1G>C, c.2966-2A>G, c.307+1del, c.307_307+1delinsTT, c.308-2A>G, c.3152+1G>T, c.3153-1G>C, c.472+1G>A, c.473-2A>C, c.473-2A>G, c.76+1G>A, c.76+2T>A, c.77-1G>A, c.916+1G>A, c.917-2A>G, p.Ala185fs c.544_551dup, p.Ala223fs c.668_671del, p.Ala228fs c.678del, p.Ala253fs c.757del, p.Ala285fs c.853_859del, p.Ala34fs c.100del, p.Ala429fs c.1285del, p.Ala505fs c.1513del, p.Ala527fs c.1578del, p.Ala715fs c.2144_2145del, p.Ala753fs c.2257del, p.Ala78Thr c.232G>A, p.Ala79fs c.234_241del, p.Ala824fs c.2470del, p.Ala915fs c.2743del, p.Arg1032Trp c.3094C>T, p.Arg1032fs c.3094del, p.Arg1033fs c.3094_3095dup, p.Arg1035fs c.3099del, p.Arg269fs c.805del, p.Arg271fs c.812_813del, p.Arg273Ter c.817C>T, p.Arg366Ter c.1096C>T, p.Arg394fs c.1181del, p.Arg531Trp c.1591C>T, p.Arg56Leu c.167G>T, p.Arg56fs c.166del, p.Arg62Ter c.184C>T, p.Arg694fs c.2080_2081insA, p.Arg744Ter c.2230C>T, p.Arg752Gln c.2255G>A, p.Arg814fs c.2442_2451del, p.Arg823Gln c.2468G>A, p.Arg823Leu c.2468G>T, p.Arg823Pro c.2468G>C, p.Arg863Ter c.2587C>T, p.Arg885fs c.2648_2652dup, p.Arg894fs c.2675_2679dup, p.Arg912fs c.2734_2738del, p.Arg922Gln c.2765G>A, p.Arg922fs c.2762dup, p.Asn220fs c.659del, p.Asn257fs c.767dup, p.Asn339fs c.1014del, p.Asp102fs c.303del, p.Asp16fs c.46del, p.Asp259fs c.774_780del, p.Asp286fs c.855del, p.Asp287fs c.853_859dup, p.Asp460fs c.1379del, p.Asp580fs c.1739del, p.Asp609Glu c.1827C>A, p.Asp609Gly c.1826A>G, p.Asp727_Ala c.2180_2317del, p.Asp896fs c.2682_2685dup, p.Cys276fs c.826del, p.Cys52Ter c.156C>A, p.Cys566fs c.1695del, p.Cys66Arg c.196T>C, p.Cys66Ter c.198C>A, p.Cys66Trp c.198C>G, p.Cys66Tyr c.197G>A, p.Cys977Ter c.2931C>A, p.Gln1008Ter c.3022C>T, p.Gln1010Ter c.3028C>T, p.Gln1065Ter c.3193C>T, p.Gln1070Ter c.3208C>T, p.Gln11Ter c.31C>T, p.Gln335Ter c.1003C>T, p.Gln376= c.1128G>A, p.Gln450Ter c.1348C>T, p.Gln450fs c.1328_1349dup, p.Gln592Ter c.1774C>T, p.Gln664Ter c.1990C>T, p.Gln676Ter c.2026C>T, p.Gln688fs c.2063_2069del, p.Gln695Ter c.2083C>T, p.Gln702Ter c.2104C>T, p.Gln75fs c.222_223del, p.Gln81fs c.234_241dup, p.Gln84fs c.234_250dup, p.Gln884Ter c.2650C>T, p.Glu177Ter c.529G>T, p.Glu229Ter c.685G>T, p.Glu230Ter c.688G>T, p.Glu372fs c.1112_1113dup, p.Glu444Ter c.1330G>T, p.Glu50Ter c.148G>T, p.Glu544fs c.1629dup, p.Glu58Ter c.172G>T, p.Glu807fs c.2419del, p.Glu876Ter c.2626G>T, p.Glu900fs c.2698dup, p.Glu90Ter c.268G>T, p.Glu929fs c.2785del, p.Glu971fs c.2911del, p.Gly1006fs c.3017del, p.Gly183fs c.548_560del, p.Gly236fs c.707_708del, p.Gly306fs c.913_914dup, p.Gly53Cys c.157G>T, p.Gly584Cys c.1750G>T, p.Gly604Cys c.1810G>T, p.Gly626fs c.1877del, p.Gly785fs c.2354del, p.Gly873fs c.2616dup, p.Gly879fs c.2634_2643del, p.Gly911fs c.2732_2766del, p.Gly921fs c.2762del, p.Gly924fs c.2767_2770dup, p.Gly947fs c.2840_2852del, p.Gly965fs c.2892dup, p.Gly969fs c.2906del, p.His302fs c.904del, p.His762fs c.2285_2286del, p.Ile571fs c.1707_1711dup, p.Ile583fs c.1745_1746dup, p.Ile593Lys c.1778T>A, p.Ile593Thr c.1778T>C, p.Ile662fs c.1983del, p.Leu155fs c.464del, p.Leu170fs c.507dup, p.Leu296fs c.885del, p.Leu343fs c.1027del, p.Leu380fs c.1139del, p.Leu433fs c.1298_1304del, p.Leu552Ter c.1655T>A, p.Leu559fs c.1676_1682del, p.Leu589fs c.1762_1765dup, p.Leu666fs c.1996del, p.Leu769fs c.2304_2305dup, p.Leu838fs c.2512del, p.Leu953fs c.2857del, p.Leu973fs c.2918del, p.Leu987fs c.2959_2960del, p.Lys135fs c.402del, p.Lys21fs c.60dup, p.Lys28Ter c.81dup, p.Lys373fs c.1113_1114del, p.Lys897fs c.2690delinsCGACAC, p.Met651_Tyr c.1956del, p.Phe125fs c.373_374insGTGG, p.Phe163fs c.486del, p.Phe720fs c.2156dup, p.Phe852fs c.2555del, p.Phe891fs c.2662_2669dup, p.Pro1030fs c.3086_3087dup, p.Pro1034fs c.3097_3098dup, p.Pro298fs c.893del, p.Pro299fs c.893dup, p.Pro632Ala c.1894C>G, p.Pro721_Glu c.2159dup, p.Pro721fs c.2162_2163del, p.Pro72Leu c.215C>T, p.Pro72Ser c.214C>T, p.Pro764fs c.2289_2290dup, p.Pro902fs c.2705del, p.Pro910fs c.2724_2728dup, p.Pro917fs c.2748_2766del, p.Pro926fs c.2775del, p.Ser182fs c.544_556del, p.Ser284fs c.850_868del, p.Ser304fs c.906_910dup, p.Ser331fs c.991del, p.Ser581Ter c.1742C>A, p.Ser668Ter c.2003C>A, p.Ser919fs c.2754_2755insG, p.Ser981fs c.2931_2941dup, p.Thr152fs c.453del, p.Thr270fs c.809_821del, p.Thr319fs c.955_956del, p.Thr421fs c.1261del, p.Thr443fs c.1326_1347del, p.Thr613Ala c.1837A>G, p.Thr613Pro c.1837A>C, p.Thr708fs c.2118_2121dup, p.Thr74fs c.220_233del, p.Thr859fs c.2574del, p.Trp1001Ter c.3002G>A, p.Trp154Ter c.461G>A, p.Trp398Ter c.1193G>A, p.Trp563Ter c.1688G>A, p.Trp568Ter c.1703G>A, p.Trp585Leu c.1754G>T, p.Trp705Ter c.2114G>A, p.Trp853Ter c.2559G>A, p.Trp927Ter c.2780G>A, p.Tyr420Ter c.1260C>G, p.Tyr542Ter c.1626C>A, p.Tyr652Ter c.1956T>A, p.Tyr673fs c.2017_2020del, p.Tyr812Ter c.2436T>A, p.Tyr845_Pro c.2536_2537del, p.Tyr845fs c.2534_2541del, p.Tyr99Ter c.297C>G, p.Val1038fs c.3108_3112dup, p.Val222fs c.665_669delinsC, p.Val279fs c.834_849del

### KCNH2 PMID 10862094

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P151fsX, S543fsX
- Extra predictions: K897K, T897T, c.1631_1632delAG, c.1631delAG, c.453delC, p.Lys897Thr c.2690A>C, p.Gly584Cys c.1750G>T, p.Thr152fs c.453del, p.Thr613Ala c.1837A>G, p.Thr613Pro c.1837A>C

### KCNH2 PMID 10841244

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Count disagreements: p.Leu552Ser c.1655T>C carriers 42 vs 44 (error -2)

### KCNH2 PMID 23864605

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: F805C, N470D, S25N, Y652T, p.Asn629Asp

### KCNH2 PMID 24667783

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A172V, E807X, W705fsX

### KCNH2 PMID 19038855

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: E698X, I560fsX, I642del
- Extra predictions: E689X c.2062C>T, F640V c.1918T>G, I642_V644del c.1926_1934del, N633S c.1898A>G, Q376sp c.1128G>A, Q901fs/71 c.2705delC, S649P c.1945T>C, p.Arg582Pro c.1745G>C

### KCNQ1 PMID 19490272

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Ala344Glu c.1031C>A, p.Ala344Gly c.1031C>G, p.Asp242Glu c.726C>A, p.Asp317Ala c.950A>C, p.Gly168Glu c.503G>A, p.Trp305Arg c.913T>A, p.Trp305Leu c.914G>T

### KCNQ1 PMID 23153844

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: P400FS, p.A341E, p.G314S, p.M159*, p.Q357R, p.Q530*, p.R243C, p.S225L, p.S349W, p.T312I, p.W305S, p.Y171*

### KCNQ1 PMID 17470695

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A150FS/133[DELCT451-452], G185S, L191FS/90[DELTGCGC572-576], R195FS/40[DELG585], R237P, R258C, P400FS/62[INSC1201-1202], R518X[1552C>T], Q530X[1588C>T], S571FS/20[DELC1714)
- Extra predictions: A344A, c.1032+2T>C, c.1032+5G>T, c.683+5G>A, p.Ala178Thr c.532G>A, p.Ala344Glu c.1031C>A, p.Ala344Gly c.1031C>G, p.Ala344Val c.1031C>T, p.Arg190Trp c.568C>T, p.Arg243His c.728G>A, p.Arg360Gly c.1078A>G, p.Arg397Trp c.1189C>T, p.Arg401fs c.1201dup, p.Arg591Cys c.1771C>T, p.Asp317Gly c.950A>G, p.Glu284Lys c.850G>A, p.Glu284del c.850_852del, p.Glu449fs c.1345dup, p.Gly186Ser c.556G>A, p.Gly189Arg c.565G>A, p.Gly269Asp c.806G>A, p.Gly57Val c.170G>T, p.Ile567Ser c.1700T>G, p.Ile567Thr c.1700T>C, p.Leu266Arg c.797T>G, p.Leu353Pro c.1058T>C, p.Lys393Asn c.1179G>T, p.Met1Leu c.1A>T, p.Pro320Leu c.959C>T, p.Pro320Ser c.958C>T, p.Pro343Ala c.1027C>G, p.Thr322Ala c.964A>G, p.Thr322Arg c.965C>G, p.Tyr278His c.832T>C, p.Val310Ile c.928G>A

### KCNQ1 PMID 14678125

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: L151, G185S, DELF340, A344/SP, A344A/SPLICE, R518X
- Extra predictions: p.Arg243His c.728G>A, p.Arg243Ser c.727C>A, p.Arg366Leu c.1097G>T, p.Arg366Pro c.1097G>C, p.Asp317Ala c.950A>C, p.Glu284del c.850_852del, p.Gly186Cys c.556G>T, p.Gly186Ser c.556G>A, p.Leu266Arg c.797T>G

### KCNQ1 PMID 28720088

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: R519X
- Extra predictions: A341V, G589D, L353L, R518A, p.Arg518*

### KCNQ1 PMID 21129503

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A341V, p.Tyr111Ser c.332A>C

### KCNQ1 PMID 25087618

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: A341V

### KCNQ1 PMID 17192539

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: P117L, G168R, R174C, R174H, L175fsX, A178T, Y184S, R190W, R190Q, R190L, S225L, R231C, R243C, V254M, H258N, H258P, E261K, G269D, S277W, V280A, V280E, V288fsX, G306R, T309R, T311I, T312I, G314S, Y315S, Y315C, G316E, P320H, P320A, G325R, F339S, A341V, L342F, A344V, SP/A344A, Q359-K362DEL, A371T, S373P, W379S, K422fsX, M476L, R539W, R555H, S566P, I567T, T587M, A590T, R591H, R594C, Q604X
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

- Extra predictions: A208V, A223A, A226A, A226V, A99V, C214C, C214Y, D202G, D202H, D202N, D75G, D75H, D75N, G216G, G219E, G219G, I100L, I204F, I204I, I204M, I227L, I77F, I77M, K218E, K91E, L203P, L76P, M210I, M210T, M83T, Q220K, S209F, S209P, S209S, S217F, S225L, S225S, S82F, S82P, S90F, S98L, T224T, T97M, V205M, V206L, V206V, V207L, V207M, V212V, V215G, V215M, V215V, V221V, V78M, V80M, V88G, V88M, c.223G>A, c.223G>C, c.224A>G, c.227T>C, c.229A>T, c.231C>G, c.232G>A, c.238G>A, c.244T>C, c.245C>T, c.248T>C, c.262G>A, c.263T>G, c.269C>T, c.271A>G, c.290C>T, c.293C>T, c.296C>T, c.298A>C, c.604+12C>T, c.604+16C>G, c.604+16C>T, c.604+17C>T, c.604+19G>A, c.604+1G>A, c.604+22G>A, c.604+24G>T, c.604+28A>C, c.604+29C>T, c.604+30G>A, c.604+33C>G, c.604+40C>G, c.604+41C>T, c.604+42A>T, c.604+43G>A, c.604+49G>A, c.604+49G>T, c.604+62C>A, c.604+72C>G, c.604+7C>T, c.604+9T>C, c.604G>A, c.605-10C>T, c.605-11G>A, c.605-16C>T, c.605-18C>A, c.605-18C>G, c.605-18C>T, c.605-24C>T, c.605-28A>G, c.605-30G>A, c.605-31C>T, c.605-33C>A, c.605-36G>A, c.605-37C>T, c.605-38G>A, c.605-41G>A, c.605-46T>G, c.605-48C>A, c.605-52G>T, c.605-65A>G, c.612C>T, c.613G>A, c.616G>C, c.618C>T, c.619G>A, c.619G>C, c.619G>T, c.623C>T, c.627C>T, c.630G>A, c.636C>G, c.641G>A, c.642C>T, c.643G>A, c.644T>G, c.645G>A, c.648C>A, c.648C>T, c.656G>A, c.657G>A, c.658C>A, c.663G>A, c.669C>T, c.672G>A, c.674C>T, c.675G>A, c.678C>G, c.678C>T, c.671C>T
- Count disagreements: T224M carriers 124 vs 88 (error +36); T224M affected 0 vs 34 (error -34); T224M unaffected 124 vs 54 (error +70)

### RYR2 PMID 25814417

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: V2475F, p.Gly357Asp c.1070G>A

### RYR2 PMID 29925740

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: EXON 3 DELETION, N4168S
- Extra predictions: N4178S c.12533A>G

### RYR2 PMID 33315912

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G357S, c.6982C>T

### RYR2 PMID 16272262

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A2387P, A2403T, A2428T, A4510T, A4607P, A4860G, E2311D, E4146K, E4950K, F4499C, G3946S, G4671R, H4833Y, I419F, I4848V, I4867M, L3778F, L433P, M450I, N2386I, N4097S, N4104K, N4895D, P164S, P2328S, Q4201R, R176Q, R2474S, R414L, R420W, R4497C, S2246L, T2504M, T4158P, V4653F, V4880A, Y2392C, p.Arg332Trp c.994C>T, p.Glu1724Val c.5171A>T, p.His4108Asp c.12322C>G, p.His4108Tyr c.12322C>T

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

- Missed gold variants: T4630C, G4749E, A2387T, D2300A, I3995V, H4579Y, D1872N, E3716Q, I4867T, L4670F, A2387V, Q32K, E1083K, V2174F, R3190Q
- Extra predictions: c.11147A>G, c.1191G>A, c.11983A>G, c.12470G>A, c.1258C>T, c.13735C>T, c.13890G>A, c.14008C>T, c.14246G>A, c.14553C>A, c.14600T>C, c.3407C>T, c.5170G>A, c.5614G>A, c.5654G>A, c.6520G>T, c.6899A>C, c.7159G>A, c.7160C>T, c.8598del, c.94C>A, c.9569G>A, p.Glu1724Val c.5171A>T, p.Glu3716Lys c.11146G>A, p.His4579Gln c.13737C>A
- Count disagreements: c.13352del carriers 1 vs 3 (error -2); c.13352del affected 1 vs 2 (error -1); c.13352del unaffected 0 vs 1 (error -1)

### RYR2 PMID 12106942

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: p.Asn2386Lys c.7158C>A

### RYR2 PMID 22677073

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: V2113M, Q2958R
- Extra predictions: A302V, D85N, F2004L, G643S, P2006A, R1047L, V1951L, W379R, c.2398+5G>T

### RYR2 PMID 30403697

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: G4722S
- Extra predictions: A103V, C1274Y, D4642H, D85N, H1042R, I882V, IVS5+1G>C, K296Q, K4742R, L714R, M0781T, Q584R, Q692K, R1458G, R485C, R4954K, S2774G, V288I, c.13781A>G, c.6224T>C, p.G4772S, p.Arg2401Leu c.7202G>T, p.Arg2474Gly c.7420A>G

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
