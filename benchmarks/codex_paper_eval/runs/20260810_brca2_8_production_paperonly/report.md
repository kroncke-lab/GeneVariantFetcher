# Codex extraction-blinded paper evaluation — `20260810_brca2_8_production`

## Technical summary

This hash-locked run evaluated **8 papers** (**8 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **17.4%**, recall **73.6%**, F1 **28.2%** (81 TP, 384 FP, 29 FN).
- Exact API telemetry: **0 total tokens** (0 input; 0 output).
- Elapsed: **0.0s wall clock**; 0.0s summed per-paper route + read time.
- Representation choices: {'text': 8}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `brca2_8_papers_20260810.tsv` (8 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- `selection.json` contains source metadata and hashes but no gold values or gold row counts. `predictions.json` was made read-only and SHA-256 locked before scoring first opened the gold CSVs.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 45 / 110 | 40.9% | 0.311 | 1.606 |
| affected | 37 / 110 | 33.6% | 1.378 | 2.736 |
| unaffected | 37 / 110 | 33.6% | 0.324 | 1.660 |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| BRCA2 | 81 | 384 | 29 | 17.4% | 73.6% | 28.2% | 40.9% / 0.311 / 1.606 | 33.6% / 1.378 / 2.736 | 33.6% / 0.324 / 1.660 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| BRCA2 | 10398279 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| BRCA2 | 15365993 | text | 19 | 105 | 8 | 15.3% | 70.4% | 25.2% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| BRCA2 | 18489799 | text | 21 | 50 | 20 | 29.6% | 51.2% | 37.5% | 51.2% / 0.000 | 51.2% / 2.143 | 51.2% / 0.000 | 0.0 | 0 |
| BRCA2 | 21356067 | text | 6 | 26 | 0 | 18.8% | 100.0% | 31.6% | 33.3% / 0.000 | 33.3% / 1.000 | 33.3% / 1.000 | 0.0 | 0 |
| BRCA2 | 22655046 | text | 10 | 3 | 0 | 76.9% | 100.0% | 87.0% | 30.0% / 3.333 | 20.0% / 0.000 | 20.0% / 5.000 | 0.0 | 0 |
| BRCA2 | 25802882 | text | 17 | 9 | 1 | 65.4% | 94.4% | 77.3% | 83.3% / 0.000 | 50.0% / 0.000 | 50.0% / 0.000 | 0.0 | 0 |
| BRCA2 | 26833046 | text | 4 | 75 | 0 | 5.1% | 100.0% | 9.6% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 0.0 | 0 |
| BRCA2 | 26848529 | text | 3 | 116 | 0 | 2.5% | 100.0% | 4.9% | 100.0% / 1.333 | 100.0% / 1.333 | 100.0% / 0.000 | 0.0 | 0 |

## Errors and representation choices

### BRCA2 PMID 10398279

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- No scored variant or count disagreement.

### BRCA2 PMID 15365993

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: p.K454, p.S455, p.A565, p.H743, p.N830, p.K1132, p.V1269, p.S2414
- Extra predictions: E1038G, IVS10-51G>T, IVS10-73T>C, IVS11+44A>G, IVS11-9T>C, IVS12+41T>C, IVS14+10G>C, IVS14+14A>G, IVS16-14T>C, IVS16-68A>G, IVS18+66G>A, IVS2-95G>T, IVS26-19G>A, IVS3+84G>C, IVS3+91T>C, IVS4+67A>C, IVS8+13delG, IVS8+54G>C, IVS8+56C>T, IVS8-64delT, IVS9-34T>C, K1183R, K1690Q, M1628T, P871L, S1713G, S308X, W1815X, Y856H, c.10234A>G, c.1041delAGCinsT, c.10462A>G, c.1093A>C, c.1114A>C, c.114G>A, c.134+84G>C, c.134+91T>C, c.1342A>C, c.1362A>G, c.1365A>G, c.1590A>G, c.1593A>G, c.1695C>A, c.1923C>A, c.203G>A, c.2082C>T, c.2201C>T, c.2229T>C, c.2311T>C, c.233G>A, c.2350A>G, c.2430T>C, c.2457T>C, c.2487delT, c.2490C>T, c.2566T>C, c.2578A>G, c.2612C>T, c.2685T>C, c.2718C>T, c.2731C>T, c.2971A>G, c.3113A>G, c.3199A>G, c.3232A>G, c.3396A>G, c.3548A>G, c.3624A>G, c.3667A>G, c.3741T>C, c.3807T>C, c.3860T>C, c.3975G>A, c.4035T>C, c.4094G>A, c.4096+44A>G, c.4158C>T, c.4185+41T>C, c.4277C>T, c.4308T>C, c.4334A>C, c.4427T>C, c.4484+14A>G, c.4562A>C, c.4837A>G, c.4883T>C, c.4956A>G, c.4987-68A>G, c.5002T>C, c.5068A>C, c.5152+66G>A, c.5187A>C, c.5214A>G, c.5333A>G, c.5445G>A, c.5564G>A, c.5785A>G, c.594-34T>C, c.6013A>G, c.6091A>G, c.6319A>G, c.7242A>G, c.7470A>G, c.81-95G>T, c.865A>C

### BRCA2 PMID 18489799

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.462_463del2, c.1796_1800del5, c.2808_2811del4, c.3076A>T, c.3109C>T, c.3744_3747del4, c.5042_c.5043del2, c.5645C>A, c.6405_6409del5, c.7471C>T, c.7913_7917del5, c.8169_8172dup4, c.8363G>A, c.9435_9436del2, c.475G>A, c.476-2A>G, c.7007G>A, c.8755-1G>A, c.9117+2T>A, c.9118-2A>G
- Extra predictions: p.Ala938ProfsX21, p.Arg504ValfsX28, p.Asn1215fsX3, p.Asn2135LysfsX3, p.Asp156X, p.Cys24Tyr, p.Cys39Arg, p.Cys61Gly, p.Cys64Tyr, p.Gln1037X, p.Gln1388GlufsX2, p.Gln1447X, p.Gln1756ProfsX74, p.Gln19X, p.Gln2491X, p.Gln356HisfsX15, p.Gln534X, p.Gln563X, p.Gln804LeufsX5, p.Gln921ArgfsX79, p.Glu23ValfsX17, p.Glu402SerfsX8, p.Glu720ArgfsX6, p.Glu732ArgfsX3, p.Glu755X, p.Gly1055AlafsX4, p.Gly2919fsX3, p.Leu1351X, p.Leu347ArgfsX27, p.Leu833X, p.Lys1026X, p.Lys468ArgfsX7, p.Phe143GlyfsX26, p.Phe2638fsX, p.Phe559X, p.Ser1248ArgfsX10, p.Ser1253ArgfsX10, p.Ser1389X, p.Ser1882X, p.Ser282TyrfsX15, p.Ser3147CysfsX2, p.Trp1782X, p.Trp1837X, p.Trp2725MetfsX3, p.Trp2788X, p.Val1234GlnfsX8, p.Val1681GlufsX7, p.Val2985GlyfsX4, p.Val3040MetfsX20, p.Val340GlyfsX6
- Count disagreements: p.Ala1922CysfsX2 affected 1 vs 0 (error +1); p.Asp252ValfsX24 affected 1 vs 0 (error +1); p.Gln2157IlefsX18 affected 1 vs 0 (error +1); p.Gln2384ArgfsX7 affected 1 vs 0 (error +1); p.Glu2846GlyfsX22 affected 14 vs 0 (error +14); p.Ile605AsnfsX11 affected 2 vs 0 (error +2); p.Leu103IlefsX10 affected 1 vs 0 (error +1); p.Leu1616LysfsX2 affected 1 vs 0 (error +1); p.Leu1908ArgfsX2 affected 1 vs 0 (error +1); p.Leu1930TyrfsX33 affected 1 vs 0 (error +1); p.Leu3135PhefsX28 affected 4 vs 0 (error +4); p.Lys1881GlnfsX27 affected 1 vs 0 (error +1); p.Lys2150SerfsX25 affected 3 vs 0 (error +3); p.Phe3155GlufsX9 affected 1 vs 0 (error +1); p.Ser2213X affected 2 vs 0 (error +2); p.Ser2252PhefsX9 affected 1 vs 0 (error +1); p.Thr1738IlefsX2 affected 1 vs 0 (error +1); p.Thr2681CysfsX11 affected 3 vs 0 (error +3); p.Thr3033AsnfsX11 affected 1 vs 0 (error +1); p.Val1283LysfsX2 affected 2 vs 0 (error +2); p.Val464GlyfsX3 affected 2 vs 0 (error +2)

### BRCA2 PMID 21356067

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: A692A, D438D, D515G, E272K, E879E, G878G, G998E, H430H, H632H, I157T, L195P, L337S, M1628T, P735P, P919S, Q559R, R1699W, S1436S, S1613G, S2414S, V455I, V932M, Y1137Y, Y334D, c.1100delC, c.1592delT
- Count disagreements: c.10234A>G affected 1 vs 0 (error +1); c.8182G>A affected 1 vs 0 (error +1); c.10234A>G unaffected 8 vs 9 (error -1); c.8182G>A unaffected 1 vs 2 (error -1)

### BRCA2 PMID 22655046

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.2808_2811delACAA, c.3124_3133delAGCAATATTA, c.9382C>T
- Count disagreements: c.5114_5117delTAAA carriers 1 vs 11 (error -10); c.5114_5117delTAAA unaffected 0 vs 10 (error -10)

### BRCA2 PMID 25802882

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Missed gold variants: c.8023A>G
- Extra predictions: I2675V, N2113fs, R2318X, S1882X, c.188T>A, c.1952_1953insG, c.2800C>T, c.5802delTTAA, p.Ser1915Ter c.5744C>A

### BRCA2 PMID 26833046

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: G1738E, c.130T>A, c.1310_1313delAAGA, c.145G>T, c.1556delA, c.1687C>T, c.1803delC, c.1813delA, c.181T>G, c.1877_1878insTAGT, c.190T>C, c.2231C>G, c.2450delA, c.2475delC, c.2523dupG, c.2787dupA, c.2808_2811delACAA, c.2812delGCAA, c.2830A>T, c.2943_2944insTGAGA, c.301+1G>C, c.3052_3053insTGAGA, c.3319G>T, c.3400G>T, c.3477_3479delAAAinsC, c.3599_3600delGT, c.3607C>T, c.3710delT, c.3718C>T, c.3765delT, c.3847_3848delGT, c.3860delA, c.3874delT, c.3904G>T, c.4035delA, c.4258delG, c.427G>T, c.4284_4285insT, c.4285_4286insT, c.469_470delAA, c.5089T>C, c.5096G>A, c.5130_5133delTGTA, c.5143A>C, c.5164_5165delAG, c.5213G>A, c.5219delT, c.5263insC, c.5266dupC, c.5341delG, c.5351delA, c.5352delC, c.5382insC, c.5503C>T, c.5722_5723delCT, c.583delT, c.5857G>T, c.6082_6086delGAAGA, c.6244G>T, c.6408_6414delAAATGTT, c.6443_6444delCT, c.6486_6489delACAA, c.6490_6492del, c.7069_7070delCT, c.7757G>A, c.7878G>C, c.7913_7917delTTCCT, c.7980T>G, c.8364G>A, c.8536G>T, c.8575delC, c.8754+3G>C, c.8953+1G>T, c.9097delA, c.9106C>T

### BRCA2 PMID 26848529

**text** — Production gvf-run strategy: deterministic regex/table sweep over the untruncated source, Kimi table routing, grok-4.3 Tier 3 extraction on a 60k-char gene-focused prompt, gpt-5.6-sol claim verification on high-risk papers, then recovery layers. It reads the same corpus/ markdown the single-model runs read; --no-source-recovery kept the source identical. Mapped to the harness 'text' route because no single harness route describes a multi-source pipeline.

- Extra predictions: c.212+1G>T, c.4986+1G>A, c.7435+53C>T, c.7976+1G>A, p.Ala634Valfs*18 c.1901delC, p.Ala938Profs*21 c.2808_2811delACAA, p.Arg1203* c.3607C>T, p.Arg155Lysfs*26 c.464_468delGAGAT, p.Arg1835* c.5503C>T, p.Arg2108Cys c.6322C>T, p.Arg2336His c.7007G>A, p.Arg435Lysfs*17 c.1303dupA, p.Arg691Serfs*10 c.2073delA, p.Asn1121Lysfs*10 c.3363_3367delTACAG, p.Asn1287Ilefs*6 c.3860delA, p.Asn1355Lysfs*10 c.4065_4068delTCAA, p.Asn1473Glnfs*8 c.4416_4417delGA, p.Asn1822Ilefs*18 c.5465delA, p.Asn2051* c.6150_6151insT, p.Asn2101_Val2102del c.6301_6306delAATGTA, p.Asn361Metfs*6 c.1082delA, p.Asn704Cysfs*7 c.2110_2111delAA, p.Asn863Lysfs*18 c.2588dupA, p.Asn991Asp c.2971A>G, p.Asp1151Metfs*4 c.3450delT, p.Asp687* c.2059_2063delGATTA, p.Asp693Asn c.2077G>A, p.Asp974Asn c.2920G>A, p.Cys328* c.981_982delAT, p.Gln1056Argfs*3 c.3167_3170delAAAA, p.Gln1144His c.3432G>C, p.Gln1291* c.3871C>T, p.Gln2655Asnfs*2 c.7963delC, p.Gln2829* c.8485C>T, p.Gln2941Leufs*34 c.8820_8823del, p.Gln3034Serfs*10 c.9098_9099insA, p.Gln3036Serfs*8 c.9105dupT, p.Gln759* c.2275C>T, p.Gln858* c.2572C>T, p.Glu1038Gly c.3113A>G, p.Glu1214* c.3640G>T, p.Glu1257Glyfs*9 c.3770_3771delAG, p.Glu2198Asnfs*4 c.6591_6592delTG, p.Glu3096* c.9286G>T, p.Glu3263Argfs*12 c.9788delA, p.Glu554* c.1660G>T, p.Gly1788Val c.5363G>T, p.Gly2901Asp c.8702G>A, p.Gly3003* c.9007G>T, p.Gly933Alafs*4 c.2798_2799delGT, p.His2021Pro, p.His372Asn c.1114C>A, p.Ile1170Phe c.3508A>T, p.Ile1318Leu c.3952A>C, p.Ile1724Lysfs*17 c.5171delT, p.Ile1824Aspfs*3 c.5470_5477delATTGGGCA, p.Ile2986Lysfs*3 c.8956_8957insAA, p.Ile332Lysfs*18 c.993_994dupA, p.Ile605Tyrfs*9 c.1813delA, p.Ile770Phefs*2 c.2307delT, p.Ile980Lys c.2939T>A, p.Leu1306Aspfs*23 c.3916_3917delTT, p.Leu1908Argfs*2 c.5722_5723delCT, p.Leu2080* c.6239T>G, p.Leu3119* c.9356_9357delTAinsG, p.Leu750Valfs*10 c.2248_2252delCTCAT, p.Leu88Alafs*12 c.262_263delCT, p.Lys1183Arg c.3548A>G, p.Lys2150Asnfs*19 c.6449_6450insTA, p.Lys307Glu c.919A>G, p.Lys3326* c.9976A>T, p.Lys355Arg c.1064A>G, p.Lys467* c.1399A>T, p.Lys519Argfs*13 c.1556delA, p.Lys585Arg c.1754A>G, p.Met2393Lysfs*18 c.7178_7179delTG, p.Met815Trpfs*10 c.2442delC, p.Phe1177Leufs*33 c.3531delT, p.Phe2801Leufs*10, p.Pro1099Leufs*10 c.3294delT, p.Pro1702Thrfs*16 c.5103dupA, p.Pro2381Hisfs*13 c.7142delC, p.Pro628Hisfs*16 c.1881delA, p.Pro981Ala c.2941C>G, p.Pro999Leu c.2996C>T, p.Ser1040Asn c.3119G>A, p.Ser1041* c.3122C>G, p.Ser1248Argfs*10 c.3744_3747delTGAG, p.Ser1613Gly c.4837A>G, p.Ser1722Tyrfs*4 c.5164_5165delAG, p.Ser1841Valfs*2 c.5521delA, p.Ser1882* c.5645C>A, p.Ser1900* c.5699C>G, p.Ser1955* c.5864C>G, p.Ser2012Glnfs*5 c.6033_6034delTT, p.Ser2120* c.6359C>G, p.Ser2267* c.6800C>A, p.Ser2984Glnfs*4 c.8950delT, p.Ser611* c.1832C>A, p.Ser645Tyr c.1934C>A, p.Ser713* c.2138C>G, p.Ser868* c.2603C>A, p.Thr2471Hisfs*4 c.7409dupT, p.Thr2746Aspfs*18 c.8234dupT, p.Thr2867Ser c.8599A>T, p.Thr3085Glnfs*19 c.9253delA, p.Thr598Ala c.1792A>G, p.Trp1718* c.5154G>A, p.Trp2725Glyfs*8 c.8172delG, p.Trp372* c.1116G>A, p.Tyr1894* c.5681dupA, p.Tyr2215*, p.Tyr2839* c.8517C>A, p.Val1120Aspfs*11 c.3359_3363delTTAAT, p.Val1145Phefs*10 c.3433delG, p.Val220Ilefs*4 c.658_659delGT
- Count disagreements: c.5576_5580delTTAAA carriers 1 vs 5 (error -4); c.5576_5580delTTAAA affected 1 vs 5 (error -4)

## Scope, method, and limitations

- Population: fixed manifest `brca2_8_papers_20260810.tsv` (8 papers); 8 per cardiac gene; every PMID has downloaded source and at least one gold assertion in each count field.
- Blinding: gold was used only for PMID eligibility and count-field presence during selection; extraction exported no gold values or row counts, and predictions were made read-only and SHA-256 locked before `score` opened gold.
- Variant metrics are micro-averaged over gold rows. Precision treats unmatched predictions as false positives, although the curated recall packet may omit some real variants.
- Count MAE/RMSE are conditional on a supplied value. Count recall must be read alongside them because abstentions and missed variants are excluded from error magnitude.
- Source acquisition and gold completeness are separate from model reading quality; abstract-only or incomplete source is retained and labeled rather than silently excluded.
- The audited notation score is primary; the preserved raw score bounds sensitivity to post-lock matching adjudication.

## Reproducibility and evidence

- `selection.json`: selected PMIDs, source paths, source hashes, and available representations.
- `predictions.json`: immutable per-paper tools, rationales, extracted variants, counts, evidence quotes, source locations, elapsed time, and token telemetry.
- `llm_traces/<GENE>/<PMID>/`: exact textual requests, safe parameters, raw provider response envelopes, parse attempts, and explicit route/final-selection events for each model call.
- `llm_traces/trace_manifest.json`: SHA-256 inventory locked before gold scoring; provider-returned reasoning summaries are retained, while private hidden chain-of-thought is not available.
- `llm_trace_report.html`: self-contained per-paper browser view of the locked trace timeline, prompts, responses, rationales, retries, and integrity state.
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
