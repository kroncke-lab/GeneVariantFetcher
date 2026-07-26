# Codex extraction-blinded paper evaluation — `20260724_legacy48_grok`

## Technical summary

This hash-locked run evaluated **48 papers** (**12 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **92.1%**, recall **52.5%**, F1 **66.9%** (526 TP, 45 FP, 475 FN).
- Exact API telemetry: **847,139 total tokens** (609,102 input; 238,037 output).
- Elapsed: **2441.3s wall clock**; 2440.2s summed per-paper route + read time.
- Representation choices: {'table': 20, 'text': 26, 'pdf': 2}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `highcarrier_48_papers_20260723.tsv` (48 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- `selection.json` contains source metadata and hashes but no gold values or gold row counts. `predictions.json` was made read-only and SHA-256 locked before scoring first opened the gold CSVs.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 406 / 1001 | 40.6% | 1.172 | 11.587 |
| affected | 236 / 1001 | 23.6% | 0.581 | 4.681 |
| unaffected | 82 / 1001 | 8.2% | 0.427 | 3.225 |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 128 | 15 | 92 | 89.5% | 58.2% | 70.5% | 58.2% / 1.086 / 7.233 | 41.4% / 1.286 / 7.381 | 25.0% / 0.527 / 3.910 |
| KCNH2 | 149 | 7 | 143 | 95.5% | 51.0% | 66.5% | 17.1% / 1.920 / 13.576 | 1.4% / 4.250 / 7.159 | 0.7% / 0.000 / 0.000 |
| KCNQ1 | 88 | 14 | 183 | 86.3% | 32.5% | 47.2% | 31.7% / 2.663 / 21.171 | 0.0% / n/a / n/a | 0.0% / n/a / n/a |
| RYR2 | 161 | 9 | 57 | 94.7% | 73.9% | 83.0% | 65.1% / 0.085 / 0.581 | 64.7% / 0.021 / 0.253 | 11.5% / 0.240 / 0.693 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 27566755 | table | 0 | 0 | 51 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.3 | 6,835 |
| SCN5A | 26669661 | text | 26 | 1 | 1 | 96.3% | 96.3% | 96.3% | 96.3% / 0.000 | 0.0% / n/a | 0.0% / n/a | 55.9 | 57,754 |
| SCN5A | 20470418 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 59.000 | 100.0% / 22.000 | 0.0% / n/a | 42.5 | 10,764 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 55.000 | 100.0% / 55.000 | 100.0% / 0.000 | 45.2 | 16,062 |
| SCN5A | 29709101 | table | 10 | 2 | 1 | 83.3% | 90.9% | 87.0% | 90.9% / 0.000 | 0.0% / n/a | 0.0% / n/a | 51.0 | 14,855 |
| SCN5A | 28341781 | text | 52 | 3 | 3 | 94.5% | 94.5% | 94.5% | 94.5% / 0.000 | 94.5% / 0.000 | 94.5% / 0.000 | 81.5 | 23,005 |
| SCN5A | 18451998 | text | 1 | 0 | 1 | 100.0% | 50.0% | 66.7% | 50.0% / 9.000 | 50.0% / 38.000 | 50.0% / 29.000 | 49.7 | 11,821 |
| SCN5A | 26921764 | table | 26 | 2 | 1 | 92.9% | 96.3% | 94.5% | 96.3% / 0.077 | 96.3% / 0.077 | 0.0% / n/a | 79.3 | 16,181 |
| SCN5A | 27554632 | text | 9 | 6 | 0 | 60.0% | 100.0% | 75.0% | 100.0% / 0.444 | 100.0% / 0.000 | 0.0% / n/a | 47.2 | 17,313 |
| SCN5A | 10590249 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 10.000 | 0.0% / n/a | 0.0% / n/a | 41.3 | 14,919 |
| SCN5A | 25051102 | table | 1 | 1 | 2 | 50.0% | 33.3% | 40.0% | 33.3% / 0.000 | 33.3% / 0.000 | 33.3% / 0.000 | 51.5 | 18,341 |
| SCN5A | 26746457 | table | 0 | 0 | 32 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 27.5 | 7,113 |
| KCNH2 | 29622001 | table | 34 | 1 | 1 | 97.1% | 97.1% | 97.1% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 84.4 | 23,778 |
| KCNH2 | 11854117 | table | 44 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 57.0 | 15,838 |
| KCNH2 | 14661677 | text | 0 | 0 | 29 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 33.3 | 21,672 |
| KCNH2 | 19160088 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 66.7% / 48.000 | 0.0% / n/a | 0.0% / n/a | 40.3 | 11,674 |
| KCNH2 | 26496715 | table | 0 | 0 | 54 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 20.8 | 18,873 |
| KCNH2 | 11844290 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 24.6 | 6,508 |
| KCNH2 | 10973849 | text | 60 | 4 | 0 | 93.8% | 100.0% | 96.8% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 84.6 | 25,635 |
| KCNH2 | 10862094 | text | 6 | 2 | 2 | 75.0% | 75.0% | 75.0% | 25.0% / 0.000 | 25.0% / 0.000 | 25.0% / 0.000 | 61.0 | 14,060 |
| KCNH2 | 10841244 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 69.7 | 15,469 |
| KCNH2 | 23864605 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 8.500 | 0.0% / n/a | 32.5 | 18,776 |
| KCNH2 | 24667783 | text | 0 | 0 | 23 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.7 | 20,894 |
| KCNH2 | 19038855 | table | 0 | 0 | 28 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 25.1 | 8,394 |
| KCNQ1 | 19490272 | table | 54 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 97.1 | 16,084 |
| KCNQ1 | 23153844 | table | 21 | 13 | 0 | 61.8% | 100.0% | 76.4% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 63.1 | 14,012 |
| KCNQ1 | 17470695 | table | 0 | 0 | 56 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 15.5 | 7,425 |
| KCNQ1 | 14678125 | text | 0 | 0 | 41 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 20.9 | 8,051 |
| KCNQ1 | 28720088 | table | 1 | 1 | 1 | 50.0% | 50.0% | 50.0% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 38.0 | 9,672 |
| KCNQ1 | 21129503 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 32.9 | 12,618 |
| KCNQ1 | 25087618 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 31.5 | 14,678 |
| KCNQ1 | 17192539 | text | 1 | 0 | 56 | 100.0% | 1.8% | 3.4% | 1.8% / 193.000 | 0.0% / n/a | 0.0% / n/a | 33.2 | 9,782 |
| KCNQ1 | 24052033 | table | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 20.2 | 8,245 |
| KCNQ1 | 18713323 | table | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 38.4 | 10,992 |
| KCNQ1 | 29197658 | text | 0 | 0 | 29 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 32.4 | 14,842 |
| KCNQ1 | 33141630 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 36.000 | 0.0% / n/a | 0.0% / n/a | 47.0 | 30,085 |
| RYR2 | 25814417 | table | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 6.000 | 0.0% / n/a | 0.0% / n/a | 35.2 | 11,907 |
| RYR2 | 29925740 | text | 50 | 1 | 1 | 98.0% | 98.0% | 98.0% | 98.0% / 0.000 | 98.0% / 0.000 | 0.0% / n/a | 192.5 | 42,542 |
| RYR2 | 33315912 | table | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 3.000 | 0.0% / n/a | 33.5 | 9,994 |
| RYR2 | 16272262 | pdf | 12 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 53.8 | 21,533 |
| RYR2 | 34202968 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 48.2 | 22,487 |
| RYR2 | 19926015 | text | 2 | 1 | 38 | 66.7% | 5.0% | 9.3% | 5.0% / 0.500 | 5.0% / 0.000 | 5.0% / 0.500 | 54.3 | 28,227 |
| RYR2 | 33686871 | text | 1 | 0 | 7 | 100.0% | 12.5% | 22.2% | 0.0% / n/a | 12.5% / 0.000 | 0.0% / n/a | 31.8 | 29,084 |
| RYR2 | 28237968 | pdf | 17 | 6 | 1 | 73.9% | 94.4% | 82.9% | 94.4% / 0.176 | 94.4% / 0.000 | 94.4% / 0.176 | 105.2 | 27,864 |
| RYR2 | 12106942 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 16.7% / 0.000 | 0.0% / n/a | 0.0% / n/a | 45.7 | 18,507 |
| RYR2 | 22677073 | table | 19 | 0 | 7 | 100.0% | 73.1% | 84.4% | 73.1% / 0.000 | 73.1% / 0.000 | 0.0% / n/a | 57.1 | 16,377 |
| RYR2 | 30403697 | table | 18 | 1 | 3 | 94.7% | 85.7% | 90.0% | 85.7% / 0.111 | 85.7% / 0.000 | 28.6% / 0.333 | 59.2 | 28,194 |
| RYR2 | 33606749 | table | 33 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 82.0 | 17,403 |

## Errors and representation choices

### SCN5A PMID 27566755

**table** — structured table rows carry the variant-level person counts for SCN5A mutations in the 406 LQT3 patients, preserving carrier/affected/unaffected evidence lost in running text aggregates

- Missed gold variants: A1330D, A1330T, A1746T, A413T, D1114E, D1790G, E1208K, E1781G, E1784K, F1596I, G1631D, G615E, I1278N, I1448L, I1768V, L1501V, L1560F, L1786Q, L604V, M1766L, M1766V, N1269S, N1325S, N406K, P.F1617DEL, P.I1762DEL, P.K1505_Q1507DEL, P.Q1507_P1509DEL, P1008S, P1021S, P627L, Q245K, Q692K, R1623L, R1623Q, R1644H, R1991Q, R340W, R971C, S941N, T1069M, T1304M, T1779M, T370M, V125L, V1667I, V1763M, V1777M, V411M, Y1767C, Y1795C

### SCN5A PMID 26669661

**text** — Text preview supplies SCN5A/LQT3-specific transmission ratios and aggregates (57% maternal origin) directly matching target gene; table preview is restricted to KCNQ1 variant-level carriers/transmissions only and omits SCN5A entirely.

- Missed gold variants: c.703+1G>A
- Extra predictions: p.0?

### SCN5A PMID 20470418

**text** — Full text provides the explicit carrier counts (26 heterozygotes) and expression-imbalance findings (13/26) along with affected/unaffected subgroup comparisons for the sole variant p.Ser1103Tyr; table preview explicitly records table_captions_count=0 and no structured variant-level rows.

- Count disagreements: p.Ser1103Tyr carriers 26 vs 85 (error -59); p.Ser1103Tyr affected 17 vs 39 (error -22)

### SCN5A PMID 28339995

**text** — Running full text contains explicit variant-level counts for the D1790G SCN5A mutation (20 carriers, 30% event rate after stopping therapy, zero events while compliant) plus kindred context; table preview supplies only metadata with no structured rows or person-level data.

- Count disagreements: D1790G carriers 30 vs 85 (error -55); D1790G affected 30 vs 85 (error -55)

### SCN5A PMID 29709101

**table** — structured table rows carry the variant-level person counts for SCN5A LoF carriers, affected status, and unaffected evidence; 5 table captions are present and running text only aggregates to 79/39 carriers without per-variant detail

- Missed gold variants: c.4813+3_4813_6dup
- Extra predictions: c.4813+3_4813+6dup, c.4978A>G (p.(Ile1660Val))

### SCN5A PMID 28341781

**text** — Only 'text' representation is available per catalog; running full text supplies the clearest continuous source of SCN5A proband registry data, genotype-phenotype descriptions, and any variant-level counts or carrier details present in the article.

- Missed gold variants: c.3840+1G>A, c.4245+1G>A, c.4732_4733dupAA
- Extra predictions: L1579SfsTer53, IVS21+1 G>A, IVS23+1 G>A

### SCN5A PMID 18451998

**text** — Pedigree images (Figure 1 panels A-O) supply the only complete per-individual variant-level carrier/affected/unaffected evidence (symptoms, QTc, BrS pattern, SND, device status) that text summaries and figure captions aggregate or omit. Requested ocr, which was unavailable; used text.

- Missed gold variants: T1304M
- Count disagreements: E1784K carriers 41 vs 50 (error -9); E1784K affected 9 vs 47 (error -38); E1784K unaffected 32 vs 3 (error +29)

### SCN5A PMID 26921764

**table** — The parsed table supplies explicit per-variant patient counts (carrier numbers) across multiple SCN5A mutations with their molecular details; running text supplies only aggregate SCN5A detection rates and clinical correlations without listing individual variants.

- Missed gold variants: c.2559delT
- Extra predictions: F853LfsX16, variant notation unavailable
- Count disagreements: G1743E carriers 4 vs 2 (error +2); G1743E affected 4 vs 2 (error +2)

### SCN5A PMID 27554632

**text** — Text preview supplies running full-text results with explicit variant details (T1247I, A1260D, G1262S in exon 21 plus 7 other missense changes), cohort sizes (27 DCM + 12 HCM + 16 controls), and pathogenicity calls. Table preview contains zero tables or variant-row counts (table_captions_count=0), only repeated header fragments.

- Extra predictions: S525S, E547E, H558R, L561L, D1182D, Y1261Y
- Count disagreements: T570N carriers 1 vs 2 (error -1); M1245I carriers 7 vs 9 (error -2); A1260D carriers 3 vs 4 (error -1)

### SCN5A PMID 10590249

**text** — Running full text supplies the complete clinical/genetic descriptions, patient selection criteria, ECG/phenotype definitions, obligate/proven carrier assignments, and sudden-death classifications for the 1795insD SCN5A variant across the multi-generation family; the supplied table preview confirms zero tables exist.

- Count disagreements: 1795insD carriers 53 vs 43 (error +10)

### SCN5A PMID 25051102

**table** — Table representation supplies structured rows with patient/control counts, variant primers for SCN5A typing, and group characteristics (AMI/VF+ n=49 etc.) that directly preserve variant-level carrier/affected/unaffected evidence.

- Missed gold variants: c.2788-6C>T, R1193Q
- Extra predictions: 3578G>A (Arg1193Glu)

### SCN5A PMID 26746457

**table** — Structured table rows carry the variant-level person counts for SCN5A (carrier/affected/unaffected evidence), unlike the aggregate summary statistics shown in text previews.

- Missed gold variants: A1180V, A1270S, A1680T, A551V, C982R, D1243N, D1275N, E1053K, E48K, F1596I, F532C, G1935S, G552W, G615E, I1836T, I848F, L1194M, L1704H, P717L, R1023H, R1193Q, R1195S, R1512W, R1739Q, R1898C, R2012C, R282C, R689H, S1904L, T1304M, V1353M, V1532I

### KCNH2 PMID 29622001

**table** — Table rows supply explicit per-variant carrier counts (n, n(%), families) for KCNH2 mutations (e.g., p.R176W n=86, p.L552S n=73) including mutation type/location, directly preserving variant-level carrier evidence; text supplies only aggregated HRs without per-variant person counts.

- Missed gold variants: P1034fsX
- Extra predictions: c.3093_3106del

### KCNH2 PMID 11854117

**table** — Structured table rows directly supply variant-level No. of Subjects counts (e.g., L552S=4, A561V=6, subtotals 35/58) for every KCNH2/HERG allele, best preserving carrier/person evidence across pore vs non-pore regions.

- No scored variant or count disagreement.

### KCNH2 PMID 14661677

**text** — Text recovers no usable article body (supplement-only mismatch). Table/PDF previews are unrelated (impact factors, reviewer instructions) and omit variant data. OCR of supplied figure images (fig_p1_1.png, fig_p1_2.png) needed for pedigree/genotype/phenotype evidence of KCNH2 variant carriers, affected/unaffected status. Requested ocr, which was unavailable; used text.

- Missed gold variants: A1058E, A190T, A203T, A915V, C723R, G187del, G187S, G873S, H254Q, K897T, K897T, K897T, K897T, L1023del, N257H, N33T, P251A, P347S, P910L, P917L, P967L, Q1068R, R1035W, R1047L, R1047L, R176W, R181Q, T367S, V215G

### KCNH2 PMID 19160088

**text** — Text preview alone supplies explicit variant-level carrier counts (KCNH2 R176W n=16) plus QT phenotype associations; table preview contains zero variant/person data despite metadata noting a caption.

- Missed gold variants: R176W
- Count disagreements: KCNH2 R176W carriers 16 vs 112 (error -96)

### KCNH2 PMID 26496715

**table** — Structured table rows preserve variant-level patient counts (No. of patients) with explicit carrier/affected evidence per row; text recovers only fragmented supplement without usable main-body context for KCNH2.

- Missed gold variants: A1017fsX, R1032fsX, D896fsX, C984fsX, G873fsX, L987fsX, P191fsX, T443fsX, T443fsX, P1034fsX, P926fsX, c.1557+1G>A, A169G, A32T, A561T, A561V, A614V, D102H, D609G, D837Y, E480V, E544A, E575K, F494del, F805S, G572S, G601S, G604S, G626S, G628S, G71W, I782N, L413P, L457P, L559H, N629S, N633S, P596S, P72T, P926S, Q247X, R366X, R534C, A797T, R823W, R863X, R912W, S620I, S818L, S818W, T613M, T634P, W154X, Y475C

### KCNH2 PMID 11844290

**text** — Only text preview is supplied in catalog; it supplies the running main text of the paper on variable LQT expression among carriers of five HERG (KCNH2) mutations and therefore contains the narrative evidence on variant-level carriers, affected and unaffected individuals.

- Missed gold variants: R752W, F805C, M124R, V822M, W1001X

### KCNH2 PMID 10973849

**text** — Full running text (548 lines, 51668 chars) supplies complete mutational screening results, gene-specific mutation totals (HERG/KCNH2 45%), type/domain breakdowns, and per-family/individual occurrence statements for the 177 variants; table preview contains only 1836 chars of main text plus schematic figures lacking any carrier/affected/unaffected counts per variant.

- Extra predictions: G1672C (A558P), C1841T (A614V), G1801A (G601S), C1920A (F640L)

### KCNH2 PMID 10862094

**text** — Full text is the clearest authoritative source; the supplied REHARVESTED FULL TEXT (426 lines, 20k+ chars from unpaywall_pdf) directly enumerates all eight HERG/KCNH2 variants (P451L, Y569H, 1631delAG, G584S, G601S, T613M, 453delC, R176W) plus phenotypic notes and polymorphism data in the 39 Finnish LQTS cohort.

- Missed gold variants: P151fsX, S543fsX
- Extra predictions: 1631delAG, 453delC

### KCNH2 PMID 10841244

**text** — Running full text supplies explicit variant-level evidence: L552S homozygotes (two siblings) with severe phenotypes (neonatal 2:1 AV block, torsades, hypoglycemia, death), heterozygotes (parents plus 38 subjects in six families) including quantitative QTc stratification (symptomatic carriers 500±59 ms, asymptomatic 452±34 ms, non-carriers 412±23 ms) and functional studies; no tables or images are supplied to capture these counts or genotype-phenotype mappings.

- No scored variant or count disagreement.

### KCNH2 PMID 23864605

**text** — Full running text contains the only variant descriptions (G601S and other LQT2 mutations) plus all functional/trafficking data; table preview supplies only metadata with no variant-level carrier/affected/unaffected counts or genotype-phenotype rows.

- Count disagreements: G601S affected 3 vs 0 (error +3); A614V affected 14 vs 0 (error +14)

### KCNH2 PMID 24667783

**text** — Running text describes family segregation analysis for KCNH2 mutations (plus pathogenicity calls using segregation + allele frequency + conservation), directly supplying variant-level carrier/affected/unaffected counts that aggregate tables omit.

- Missed gold variants: A561V, T613M, G262fsX, G785fsX, L559H, M645V, N629S, P596L, P72R, R176W, R534C, R582C, S428X, S818L, A172V, E807X, G238R, G880V, L343fsX, R356H, S855R, S890C, W705fsX

### KCNH2 PMID 19038855

**table** — structured table rows carry the variant-level person counts

- Missed gold variants: A193fsX, A558P, C64Y, D456Y, D501H, E698X, E876X, F640L, G1036fsX, G306W, G572S, G925fsX, I30T, I560fsX, I642del, M645L, N629I, P241L, P872fsX, Q676fsX, R176W, R252fsX, R366X, R534C, R582C, R73fsX, T613M, Y99S

### KCNQ1 PMID 19490272

**table** — Table supplies explicit variant-level person counts (n) for each listed KCNQ1 mutation together with aggregate carrier/event statistics (cardiac events, ACA, SCD) stratified by tertile; text supplies only narrative summary and overall N without per-variant breakdown.

- No scored variant or count disagreement.

### KCNQ1 PMID 23153844

**table** — Table directly lists 34 KCNQ1 mutations with per-variant patient counts (e.g., V254M n=118, G168R n=87), preserving variant-level carrier evidence and enabling linkage to clinical event rates (syncope/ACA/SCD) shown in adjacent phenotype tables.

- Extra predictions: M1V, M159sp, Y171X, L191fs/90, S225L, R243C, W305S, T312I, G314S, A341E, S349W, P400fs/62, Q530X

### KCNQ1 PMID 17470695

**table** — Structured table rows explicitly preserve variant-level carrier, affected, and unaffected counts (with table_captions_count=5 and main_text reference to 581 subjects across 74 mutations), outperforming running full-text descriptions for precise person-level evidence.

- Missed gold variants: M1V, T144A, A150FS/133[DELCT451-452], E160K, G168R, Y171X[513 C>G], R174H, A178P, Y184S, G185S, G189E, R190Q, L191FS/90[DELTGCGC572-576], R195FS/40[DELG585], S225L, A226V, R237P, D242N, R243C, V254M, R258C, R259C, L266P, G269S, L273F, I274V, S277L, G292D, F296S, G306R, T312I, G314S, Y315C, Y315S, P320H, T322M, G325R, DELF340[DELCTT1017-1019], A341E, A341V, P343S, S349W, S373P, P400FS/62[INSC1201-1202], I517T, R518X[1552C>T], M520R, V524G, Q530X[1588C>T], R562M, S566F, S571FS/20[DELC1714), R591H, R594Q, D611Y, A636FS/28[DELC1909]

### KCNQ1 PMID 14678125

**text** — Text preview supplies the full methods/results narrative with aggregate KCNQ1-region patient counts (164 pre-pore, 101 pore, 29 post-pore), QTc values, and cardiac-event rates through age 40; no tables exist (table_captions_count=0) and PDF/OCR previews are empty, so running text is the only source of phenotypic evidence.

- Missed gold variants: T144A, L151, G168R, Y171X, V172M, A178P, G185S, R190Q, S225L, R243C, V254L, V254M, L266P, G269D, G269S, L273F, E284K, A300T, G306R, T311I, T312I, Y315C, D317G, DELF340, A341E, A341V, A344/SP, A344A/SPLICE, G345R, S349W, Q357H, R366Q, R366W, K393N, R518X, V524G, Q530X, R539W, S566F, I567S, R594Q

### KCNQ1 PMID 28720088

**table** — structured table rows carry the variant-level person counts for the two KCNQ1 founder mutations (Y111C n=148, R518* n=79) plus genotype-negatives, with carrier/affected/unaffected stratification by QTc and NOS1AP genotypes

- Missed gold variants: R519X
- Extra predictions: R518*

### KCNQ1 PMID 21129503

**text** — Only text representation is supplied and it alone contains explicit variant-level evidence including 170 carriers of Y111C/KCNQ1 across 37 families plus founder-couple and haplotype details that preserve carrier counts and genealogic context.

- No scored variant or count disagreement.

### KCNQ1 PMID 25087618

**text** — Full text supplies the explicit cohort counts (168 KCNQ1 A341V carriers + 181 non-carriers), genotype/haplotype-level risk effects on cardiac events/QTc/disease severity, and phenotype associations that best preserve variant-level carrier/affected/unaffected evidence.

- No scored variant or count disagreement.

### KCNQ1 PMID 17192539

**text** — Running full text is the clearest source and alone preserves all available evidence on nuclear families, parental/descendant carriers (480 type-1/KCNQ1), clinical symptoms, and affected/unaffected phenotype data for the 142 mutations.

- Missed gold variants: P117L, G168R, R174C, R174H, L175fsX, A178T, Y184S, R190W, R190Q, R190L, S225L, R231C, R243C, V254M, H258N, H258P, E261K, G269D, S277W, V280A, V280E, V288fsX, G306R, T309R, T311I, T312I, G314S, Y315S, Y315C, G316E, P320H, P320A, G325R, F339S, A341V, L342F, A344V, SP/A344A, Q359-K362DEL, R366W, A371T, S373P, W379S, K422fsX, M476L, M520R, R539W, R555C, R555H, S566P, I567T, T587M, A590T, R591H, R594C, Q604X
- Count disagreements: KCNQ1 G589D carriers 243 vs 50 (error +193)

### KCNQ1 PMID 24052033

**table** — Tables (explicitly captioned as Table 1 Clinical data, Table 2 Age/sex, Table 3/4 vectorcardiographic comparisons by LQT1 mutation) contain the structured variant-level counts and clinical/phenotype breakdowns for the KCNQ1 founder mutations Y111C and R518X, including carrier/affected/unaffected evidence that running text only summarizes narratively.

- No scored variant or count disagreement.

### KCNQ1 PMID 18713323

**table** — Structured table rows carry the variant-level person counts and matched KCNQ1 mutation details (carrier/affected/unaffected evidence) for the six identical mutations across ethnic groups; running text provides only aggregates.

- No scored variant or count disagreement.

### KCNQ1 PMID 29197658

**text** — Full-text narrative supplies complete counts (244 MVs, 29 gnomAD-common, 157 gnomAD-absent, 7 low in-silico) plus functional data for F127L/P477L/L619M; supplied table preview contains only metadata and figure captions, not variant rows or clinical carrier/affected status.

- Missed gold variants: P73T, V110I, T153M, V172M, R195Q, P197L, I274V, A287E, G292D, R293C, A300T, K362R, A370V, K393N, R397W, A399S, P408A, D446E, P448R, R452Q, R452W, G460S, V576I, G589D, T600M, D611N, G629S, G643S, V648I

### KCNQ1 PMID 33141630

**text** — Main text directly supplies variant identification, exact carrier counts (124/5521), frequency, and carrier vs non-carrier QTc phenotype comparison needed for variant-level affected/unaffected evidence.

- Count disagreements: c.671C>T (p.T224M) carriers 124 vs 88 (error +36)

### RYR2 PMID 25814417

**table** — Structured table rows preserve variant-level patient details (gender, age at diagnosis, previous symptoms, ventricular arrhythmias in basal tests) for p.G357S_RyR2 carriers, directly supplying carrier/affected/unaffected counts and phenotype evidence that running text summarizes only at aggregate level (179 carriers, 6 SCD).

- Count disagreements: p.G357S carriers 185 vs 179 (error +6)

### RYR2 PMID 29925740

**text** — Text preview alone contains the explicit aggregate carrier counts for RYR2 (9/104 LQTS-diagnosed patients) plus clinical context on misdiagnosis reasons and modified Schwartz scoring that directly references genotype-phenotype linkage; table preview supplies only metadata with zero variant rows, person counts, or affected/unaffected evidence.

- Missed gold variants: N4168S
- Extra predictions: p.N4178S (c.12533a>g)

### RYR2 PMID 33315912

**table** — Table rows directly supply structured variant-level counts for the RYR2 P2328S carriers (e.g., 41 No CE vs 18 CE, proband status, family history, β-blocker use, follow-up intervals) that best preserve carrier/affected/unaffected person-level evidence; running text supplies only summary aggregates.

- Count disagreements: P2328S affected 18 vs 15 (error +3)

### RYR2 PMID 16272262

**pdf** — PDF preserves complete article layout and any embedded tables (e.g., Table 2 referenced alongside pedigrees/fig 1) that list variant-specific carrier counts, affected/unaffected status, and clinical data for the 13 RYR2 mutations; text is limited to repeated abstract snippets only, while parsed table fragments are incomplete and fragmented.

- No scored variant or count disagreement.

### RYR2 PMID 34202968

**text** — Pedigree/figure images supply the complete variant-level genotype-phenotype evidence (carrier status, affected/unaffected individuals, homozygotes, segregation) across the three families that textual summaries only partially describe. Requested ocr, which was unavailable; used text.

- No scored variant or count disagreement.

### RYR2 PMID 19926015

**text** — Running full text supplies the core mutational analysis results (155 patients screened, 73 mutation-positive, 63 distinct variants with 47% prevalence, 34 novel, exon distribution) plus study context that directly supports variant-level carrier/affected counts and phenotype correlations; table and pdf previews contain only species-conservation alignments for already-identified variants and lack person-level evidence.

- Missed gold variants: L62F, M81L, E243K, F329L, R332W, G357S, V377M, T415R, R420Q, V507I, A549V, R739H, R1013Q, T1107M, A1136V, E1837K, E2045G, Y2156C, H2168Q, E2183V, D2216V, E2296Q, R2420W, M3972I, D3973H, S4124G, R4157Q, Q4159P, N4178S, E4187Q, G4315E, K4650E, N4736 DEL, R4790Q, K4805R, R4822H, L3974Q, K3997E
- Extra predictions: exon 3 deletion (1.1 kb)
- Count disagreements: Y4149S carriers 2 vs 1 (error +1); Y4149S unaffected 1 vs 0 (error +1)

### RYR2 PMID 33686871

**text** — figure/pedigree images are necessary because textual representations omit the genotype/phenotype evidence Requested ocr, which was unavailable; used text.

- Missed gold variants: I3995V, D4112N, T4196I, D4646A, Q4879H, K4594R, I2075T

### RYR2 PMID 28237968

**pdf** — PDF-layout extraction supplies the most complete running text of results (probands/relatives counts, mutation descriptions, and any supplemental variant annotations) without the truncation or repetition visible in HTML previews; no tables exist and figures are only survival/schematic plots, so PDF text best retains carrier/affected/unaffected evidence.

- Missed gold variants: E1083K
- Extra predictions: p.Glu1724Lys (c.5170G>A), p.Arg4157Gln (c.12470G>A), p.Phe4851Leu (c.14553C>A), p.Arg420Trp (c.1258C>T), p.Ala1136Val (c.3407C>T), p.Gly1886Ser (c.5656G>A), p.Ala1136Val (c.3407C>T), p.Gly1885Glu (c.5654G>A), p.Gly1886Ser (c.5656G>A)
- Count disagreements: p.Ala2387Val (c.7160C>T) carriers 1 vs 4 (error -3); p.Ala2387Val (c.7160C>T) unaffected 0 vs 3 (error -3)

### RYR2 PMID 12106942

**text** — Main text supplies the sole detailed narrative on the six distinct RyR2 mutations, exact carrier counts (43/81), symptomatic status (28/43 affected), sudden-death events, and per-family phenotypic correlations; table preview explicitly reports zero tables and zero table captions, while PDF/OCR options are unavailable.

- No scored variant or count disagreement.

### RYR2 PMID 22677073

**table** — Table preserves structured variant-level rows with case-specific details (gene, exon, nucleotide/amino-acid change, SUD event, sentinel event, family history, sex/age/ethnicity) for the 20 CPVT1/RYR2-positive cases plus carrier/polymorphism flags, directly supplying affected/unaffected evidence absent from narrative summaries.

- Missed gold variants: T1107M, V2113M, G1885E, G1886S, Q2958R, A1136V, R4037C

### RYR2 PMID 30403697

**table** — Table rows enumerate per-patient RYR2 variants with explicit carrier status, SCA/phenotype, and family segregation (affected/unaffected counts) that narrative summaries in text/pdf omit.

- Missed gold variants: G4722S, c.3599-9delT, c.14091-11dupT
- Extra predictions: RYR2-p.G4772S
- Count disagreements: RYR2-p.R2028H carriers 2 vs 3 (error -1); RYR2-p.Y4721C carriers 2 vs 3 (error -1); RYR2-p.R2028H unaffected 1 vs 2 (error -1); RYR2-p.Y4721C unaffected 1 vs 2 (error -1)

### RYR2 PMID 33606749

**table** — Table rows supply per-proband RYR2 variant details (nucleotide/amino-acid changes, location), inheritance origin (de novo/maternal/paternal), genotyped family members, and explicit parent phenotypes (none/syncope/AF/CPA), directly preserving variant-level carrier, affected, and unaffected counts that text aggregates omit.

- No scored variant or count disagreement.

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
