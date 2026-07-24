# Codex blinded paper evaluation — `codex_20260723_highcarrier48_routed`

## Technical summary

This hash-locked run evaluated **48 papers** (**5 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **81.0%**, recall **61.0%**, F1 **69.6%** (611 TP, 143 FP, 390 FN).
- Exact API telemetry: **856,768 total tokens** (628,595 input; 228,173 output).
- Elapsed: **3490.6s wall clock**; 3489.4s summed per-paper route + read time.
- Representation choices: {'text': 22, 'table': 19, 'ocr': 4, 'pdf': 3}.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 511 / 1001 | 51.0% | 0.523 | 5.309 |
| affected | 300 / 1001 | 30.0% | 0.383 | 3.019 |
| unaffected | 202 / 1001 | 20.2% | 0.584 | 2.856 |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 136 | 67 | 84 | 67.0% | 61.8% | 64.3% | 61.4% / 0.637 / 5.247 | 30.0% / 0.909 / 5.405 | 28.6% / 0.714 / 3.966 |
| KCNH2 | 163 | 44 | 129 | 78.7% | 55.8% | 65.3% | 37.0% / 1.000 / 9.281 | 18.2% / 0.396 / 2.014 | 3.1% / 3.333 / 5.249 |
| KCNQ1 | 117 | 17 | 154 | 87.3% | 43.2% | 57.8% | 42.8% / 0.310 / 3.343 | 11.1% / 0.000 / 0.000 | 0.0% / n/a / n/a |
| RYR2 | 195 | 15 | 23 | 92.9% | 89.4% | 91.1% | 69.7% / 0.243 / 0.774 | 69.3% / 0.225 / 1.980 | 59.6% / 0.331 / 1.774 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 27566755 | text | 51 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 11.8% / 0.000 | 11.8% / 0.000 | 164.2 | 20,160 |
| SCN5A | 26669661 | table | 1 | 26 | 26 | 3.7% | 3.7% | 3.7% | 3.7% / 0.000 | 0.0% / n/a | 0.0% / n/a | 78.1 | 69,473 |
| SCN5A | 20470418 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 59.000 | 100.0% / 22.000 | 0.0% / n/a | 13.1 | 8,240 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 46.2 | 16,188 |
| SCN5A | 29709101 | table | 2 | 10 | 9 | 16.7% | 18.2% | 17.4% | 18.2% / 0.000 | 18.2% / 0.000 | 0.0% / n/a | 70.4 | 14,651 |
| SCN5A | 28341781 | text | 45 | 10 | 10 | 81.8% | 81.8% | 81.8% | 81.8% / 0.000 | 81.8% / 0.000 | 81.8% / 0.000 | 204.3 | 24,841 |
| SCN5A | 18451998 | table | 1 | 1 | 1 | 50.0% | 50.0% | 50.0% | 50.0% / 9.000 | 50.0% / 38.000 | 50.0% / 29.000 | 126.4 | 16,374 |
| SCN5A | 26921764 | table | 23 | 4 | 4 | 85.2% | 85.2% | 85.2% | 85.2% / 0.087 | 0.0% / n/a | 0.0% / n/a | 103.4 | 16,079 |
| SCN5A | 27554632 | text | 9 | 6 | 0 | 60.0% | 100.0% | 75.0% | 88.9% / 0.500 | 88.9% / 0.000 | 100.0% / 0.444 | 51.0 | 17,219 |
| SCN5A | 10590249 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 58.0 | 15,921 |
| SCN5A | 25051102 | table | 0 | 9 | 3 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 46.3 | 17,157 |
| SCN5A | 26746457 | text | 2 | 0 | 30 | 100.0% | 6.2% | 11.8% | 6.2% / 6.000 | 6.2% / 0.000 | 6.2% / 6.000 | 30.9 | 11,122 |
| KCNH2 | 29622001 | table | 23 | 12 | 12 | 65.7% | 65.7% | 65.7% | 65.7% / 0.043 | 20.0% / 0.000 | 0.0% / n/a | 122.6 | 24,449 |
| KCNH2 | 11854117 | table | 37 | 7 | 7 | 84.1% | 84.1% | 84.1% | 84.1% / 0.000 | 0.0% / n/a | 0.0% / n/a | 68.8 | 15,100 |
| KCNH2 | 14661677 | ocr | 0 | 0 | 29 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.6 | 3,157 |
| KCNH2 | 19160088 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 66.7% / 48.000 | 0.0% / n/a | 0.0% / n/a | 54.8 | 12,792 |
| KCNH2 | 26496715 | table | 37 | 17 | 17 | 68.5% | 68.5% | 68.5% | 68.5% / 0.000 | 68.5% / 0.000 | 0.0% / n/a | 140.9 | 24,928 |
| KCNH2 | 11844290 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 9.6 | 4,921 |
| KCNH2 | 10973849 | text | 55 | 5 | 5 | 91.7% | 91.7% | 91.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 202.3 | 31,941 |
| KCNH2 | 10862094 | text | 6 | 3 | 2 | 66.7% | 75.0% | 70.6% | 75.0% / 1.500 | 75.0% / 0.667 | 75.0% / 2.167 | 152.3 | 21,171 |
| KCNH2 | 10841244 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 100.0% / 0.000 | 100.0% / 0.000 | 48.4 | 11,713 |
| KCNH2 | 23864605 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 8.500 | 100.0% / 8.500 | 79.8 | 17,966 |
| KCNH2 | 24667783 | table | 0 | 0 | 23 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.2 | 6,080 |
| KCNH2 | 19038855 | text | 0 | 0 | 28 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 68.3 | 17,214 |
| KCNQ1 | 19490272 | table | 54 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 150.5 | 15,051 |
| KCNQ1 | 23153844 | table | 21 | 13 | 0 | 61.8% | 100.0% | 76.4% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 69.0 | 12,151 |
| KCNQ1 | 17470695 | table | 0 | 0 | 56 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.3 | 6,923 |
| KCNQ1 | 14678125 | text | 0 | 0 | 41 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.0 | 6,694 |
| KCNQ1 | 28720088 | text | 1 | 1 | 1 | 50.0% | 50.0% | 50.0% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 19.0 | 37,180 |
| KCNQ1 | 21129503 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 20.5 | 11,341 |
| KCNQ1 | 25087618 | ocr | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 11.3 | 7,766 |
| KCNQ1 | 17192539 | text | 1 | 0 | 56 | 100.0% | 1.8% | 3.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 28.2 | 8,392 |
| KCNQ1 | 24052033 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 20.9 | 11,950 |
| KCNQ1 | 18713323 | table | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 42.0 | 8,735 |
| KCNQ1 | 29197658 | table | 29 | 3 | 0 | 90.6% | 100.0% | 95.1% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 84.3 | 13,042 |
| KCNQ1 | 33141630 | pdf | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 36.000 | 100.0% / 0.000 | 0.0% / n/a | 59.1 | 36,548 |
| RYR2 | 25814417 | table | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 6.000 | 100.0% / 24.000 | 100.0% / 19.000 | 82.8 | 15,945 |
| RYR2 | 29925740 | text | 50 | 1 | 1 | 98.0% | 98.0% | 98.0% | 98.0% / 0.000 | 98.0% / 0.000 | 98.0% / 0.000 | 140.4 | 22,882 |
| RYR2 | 33315912 | table | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 3.000 | 0.0% / n/a | 20.0 | 8,071 |
| RYR2 | 16272262 | pdf | 12 | 0 | 0 | 100.0% | 100.0% | 100.0% | 50.0% / 0.833 | 50.0% / 0.167 | 75.0% / 0.444 | 137.2 | 26,883 |
| RYR2 | 34202968 | ocr | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 72.7 | 14,629 |
| RYR2 | 19926015 | text | 34 | 1 | 6 | 97.1% | 85.0% | 90.7% | 5.0% / 0.500 | 5.0% / 0.000 | 2.5% / 1.000 | 119.7 | 32,866 |
| RYR2 | 33686871 | ocr | 1 | 0 | 7 | 100.0% | 12.5% | 22.2% | 12.5% / 0.000 | 12.5% / 0.000 | 12.5% / 0.000 | 55.9 | 17,967 |
| RYR2 | 28237968 | pdf | 16 | 9 | 2 | 64.0% | 88.9% | 74.4% | 88.9% / 0.188 | 88.9% / 0.000 | 88.9% / 0.188 | 123.8 | 26,097 |
| RYR2 | 12106942 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 33.3% / 2.000 | 16.7% / 0.000 | 16.7% / 4.000 | 32.2 | 16,914 |
| RYR2 | 22677073 | table | 22 | 1 | 4 | 95.7% | 84.6% | 89.8% | 84.6% / 0.045 | 84.6% / 0.045 | 0.0% / n/a | 65.7 | 16,328 |
| RYR2 | 30403697 | table | 20 | 1 | 1 | 95.2% | 95.2% | 95.2% | 95.2% / 0.200 | 95.2% / 0.000 | 95.2% / 0.200 | 54.3 | 27,175 |
| RYR2 | 33606749 | table | 31 | 2 | 2 | 93.9% | 93.9% | 93.9% | 93.9% / 0.419 | 93.9% / 0.161 | 93.9% / 0.258 | 101.8 | 16,381 |

## Errors and representation choices

### SCN5A PMID 27566755

**text** — The running full text is the most complete available representation, including the SCN5A cohort composition, mutation context, and clinical outcomes. The table artifact contains no preserved variant-level rows, while OCR is limited to two figures and is unlikely to provide carrier, affected, or unaffected counts.

- No scored variant or count disagreement.

### SCN5A PMID 26669661

**table** — The supplemental tables provide structured variant-level rows with family, total-carrier, and maternal/paternal transmission counts. This is more suitable for SCN5A carrier evidence than the running text, which presents only aggregate LQT3 results; the figures are not necessary.

- Missed gold variants: A1186T, A1221V, A1330T, A1746T, D1114E, D1790G, E1032K, E1784K, G615E, G797V, I1278N, I1768V, P.F1617DEL, P.Q1507_P1509DEL, P.Y1795_E1796INSD, P701L, Q1033R, R1644H, R18Q, R814Q, R975W, S1079F, T1304M, T1620M, V411M, W1095C
- Extra predictions: SCN5A c.53G>A (p.(Arg18Gln)), SCN5A c.1231G>A (p.(Val411Met)), SCN5A c.1844G>A (p.(Gly615Glu)), SCN5A c.2102C>T (p.(Pro701Leu)), SCN5A c.2390G>T (p.(Gly797Val)), SCN5A c.2441G>A (p.(Arg814Gln)), SCN5A c.2923C>T (p.(Arg975Trp)), SCN5A c.3094G>A (p.(Glu1032Lys)), SCN5A c.3098A>G (p.(Gln1033Arg)), SCN5A c.3236C>T (p.(Ser1079Phe)), SCN5A c.3285G>T (p.(Trp1095Cys)), SCN5A c.3342C>A (p.(Asp1114Glu)), SCN5A c.3556G>A (p.(Ala1186Thr)), SCN5A c.3662C>T (p.(Ala1221Val)), SCN5A c.3833T>A (p.(Ile1278Asn)), SCN5A c.3911C>T (p.(Thr1304Met)), SCN5A c.3988G>A (p.(Ala1330Thr)), SCN5A c.4519_4527delCAGAAGCCC (p.(Gln1507_Pro1509del)), SCN5A c.4850_4852del (p.(Phe1617del)), SCN5A c.4859C>T (p.(Thr1620Met)), SCN5A c.4931G>A (p.(Arg1644His)), SCN5A c.5236G>A (p.(Ala1746Thr)), SCN5A c.5302A>G (p.(Ile1768Val)), SCN5A c.5350G>A (p.(Glu1784Lys)), SCN5A c.5369A>G (p.(Asp1790Gly)), SCN5A c.5385_5387dupTGA (p.(Tyr1795_Glu1796insAsp))

### SCN5A PMID 20470418

**text** — The running full text directly reports variant-level carrier and phenotype counts, including 26 heterozygous infant hearts and SIDS versus other-cause groups. No structured tables are present, and the single image does not appear necessary to recover genotype–phenotype evidence.

- Count disagreements: SCN5A c.3308C>A (p.Ser1103Tyr; also designated p.Ser1102Tyr) carriers 26 vs 85 (error -59); SCN5A c.3308C>A (p.Ser1103Tyr; also designated p.Ser1102Tyr) affected 17 vs 39 (error -22)

### SCN5A PMID 28339995

**text** — The publisher full text contains the SCN5A D1790G carrier cohort and cardiac-event counts, whereas the table representation is only artifact metadata and does not provide variant-level person rows.

- No scored variant or count disagreement.

### SCN5A PMID 29709101

**table** — The PMC XML artifact includes five structured tables, which are most likely to preserve variant-level SCN5A carrier and phenotype counts. The running-text preview is truncated and repetitive, while the figures concern ECG and risk-marker analyses rather than pedigree genotype evidence.

- Missed gold variants: c.2268_2271del, c.4813+3_4813_6dup, D349N, H886Q, M438L, R1632H, R367H, R376H, W822X
- Extra predictions: c.1045G>A (p.(Asp349Asn)), c.1100G>A (p.(Arg367His)), c.1127G>A (p.(Arg376His)), c.1312A>T (p.(Met438Leu)), c.2268_2271del (p.(Phe756Leufs*8)), c.2466G>A (p.(Trp822*)), c.2658T>A (p.(His886Gln)), c.4813+3_4813+6dup, c.4895G>A (p.(Arg1632His)), c.4978A>G (p.(Ile1660Val))

### SCN5A PMID 28341781

**text** — Text is the only available representation and contains the original research article’s cohort, genotype, phenotype, and follow-up content; no table, PDF, or OCR representation is supplied.

- Missed gold variants: c.1537delC, c.3840+1G>A, c.4245+1G>A, c.4299+1delG, c.4389_4396delCCTCTTTA, c.4732_4733dupAA, c.5157delC, c.5290delG, c.870delC, P.N1380DEL
- Extra predictions: 870delC (N291TfsTer52), 1537delC (R513VfsTer8), 4140_4142delCAA (del1380N), 4389_4396delCCTCTTTA (L1464WfsTer5), 4732_4733dupAA (L1579SfsTer53), 5157delC (I1720SfsTer67), 5290delG (V1764SfsTer23), IVS21+1 G>A, IVS23+1 G>A, IVS24+1 delG

### SCN5A PMID 18451998

**table** — The parsed artifact preserves exact E1784K person counts from the pedigree caption: 41 mutation carriers, including 9 symptomatic and 32 asymptomatic carriers, with additional phenotype counts. This is clearer for variant-level carrier, affected, and unaffected curation than the abstract percentages or OCR of pedigree symbols.

- Missed gold variants: T1304M
- Extra predictions: V1098L
- Count disagreements: E1784K carriers 41 vs 50 (error -9); E1784K affected 9 vs 47 (error -38); E1784K unaffected 32 vs 3 (error +29)

### SCN5A PMID 26921764

**table** — The preserved structured table directly links each SCN5A nucleotide and coding variant to the number of patients carrying it. The running text mainly provides aggregate genotype and outcome counts, so the table is the clearest authoritative source for variant-level person counts, although affected versus unaffected status is not resolved per variant in the preview.

- Missed gold variants: c.2559delT, c.2582_2583delTT, c.3667delG, P.Y1795_E1796INSD
- Extra predictions: c.2559delT (p.F853LfsX16), c.2582_2583delTT (p.F861WfsX90), c.3666delG (p.A1223PfsX7), c.5387_5388insTGA (p.1795_1796insD)
- Count disagreements: c.5228G>A (p.G1743E) carriers 4 vs 2 (error +2)

### SCN5A PMID 27554632

**text** — The running full text is the clearest available representation for variant and cohort evidence. The table artifact appears to be a malformed layout conversion without reliably structured variant-level rows or captions.

- Extra predictions: S525S (c.1575C>T; rs370684004), E547E (c.1641A>G; rs749857735), H558R (c.1673A>G; rs1805124), L561L (c.1681C>T; rs45522138), D1182D (c.3546C>T; rs752243560), Y1261Y (c.3783C>T; rs750678689)
- Count disagreements: T570N (c.1709C>A; rs780825701) carriers 1 vs 2 (error -1); M1245I (c.3735G>A; rs753677814) carriers 7 vs 9 (error -2); A1260D (c.3779C>A; rs199540275) carriers 3 vs 4 (error -1); T570N (c.1709C>A; rs780825701) unaffected 0 vs 1 (error -1); M1245I (c.3735G>A; rs753677814) unaffected 0 vs 2 (error -2); A1260D (c.3779C>A; rs199540275) unaffected 0 vs 1 (error -1)

### SCN5A PMID 10590249

**text** — The running full text is the most complete available representation and includes the clinical/genetic methods and family-level evidence. The table representation contains only figure metadata and captions, with no structured variant-level person counts or pedigree data.

- Missed gold variants: P.Y1795_E1796INSD
- Extra predictions: 1795insD (5537insTGA)

### SCN5A PMID 25051102

**table** — Structured article tables are preserved and are the clearest source for SCN5A variant-level cohort and genotype counts; the text preview is noisy and duplicated, while the available PDF and OCR concern an unrelated library form.

- Missed gold variants: c.2788-6C>T, F1596I, R1193Q
- Extra predictions: 87G>A, 630G>A, 1673A>G (His558Arg), IVS16-6C>T, 3183G>A, 3578G>A, 4509C>T, 4786T>A (Phe1596Ile), 5457T>A

### SCN5A PMID 26746457

**text** — The PMC full text is available and contains the carrier and phenotype counts. The table representation provides only artifact metadata with no variant-level rows, while OCR contains a single figure image without indicated genotype/phenotype evidence.

- Missed gold variants: A1180V, A1270S, A1680T, A551V, C982R, D1275N, E1053K, E48K, F1596I, F532C, G1935S, G552W, G615E, I1836T, I848F, L1194M, L1704H, P717L, R1023H, R1195S, R1512W, R1739Q, R1898C, R2012C, R282C, R689H, S1904L, T1304M, V1353M, V1532I
- Count disagreements: SCN5A p.R1193Q carriers 19 vs 7 (error +12); SCN5A p.R1193Q unaffected 12 vs 0 (error +12)

### KCNH2 PMID 29622001

**table** — Structured table rows provide KCNH2 variant-level carrier and family counts, including p.R176W and p.L552S, and are the clearest available source for person-level aggregation.

- Missed gold variants: D460fsX, E544fsX, Q725fsX, T74fsX, L987fsX, G1006fsX, P1034fsX, T1054fsX, T152fsX, V215X, R281fsX, A285fsX
- Extra predictions: c.221_242del, c.453delC, c.643delG, c.842dupG, c.853_859dupGCCGACG, c.1379delA, c.1631_1632delAG, c.2173_2196del, c.2959_2960delCT, c.3017delG, c.3093_3106del, c.3160dupA
- Count disagreements: p.L552S carriers 73 vs 74 (error -1)

### KCNH2 PMID 11854117

**table** — The structured table directly preserves KCNH2 variant-level subject counts, family annotations, exon, and mutation type. The running-text preview is fragmented and repeatedly embeds the same table, while no PDF or OCR content is available.

- Missed gold variants: A83fsX, G925fsX, I593X, L799SP, P968fsX, Q376SP, V295fsX
- Extra predictions: 593stop, fsA83/37aa, fsV295/63aa, Q376splice, L799splice, fsGly925, fsP968/4aa

### KCNH2 PMID 14661677

**ocr** — The text has no usable article body, while the table and PDF previews contain only unrelated journal materials. The available figure images are therefore the only representation that may preserve genotype, pedigree, and variant-level carrier/phenotype evidence.

- Missed gold variants: A1058E, A190T, A203T, A915V, C723R, G187del, G187S, G873S, H254Q, K897T, K897T, K897T, K897T, L1023del, N257H, N33T, P251A, P347S, P910L, P917L, P967L, Q1068R, R1035W, R1047L, R1047L, R176W, R181Q, T367S, V215G

### KCNH2 PMID 19160088

**text** — The full text directly reports variant-level evidence for KCNH2 R176W, including 16 carriers and the associated adjusted QT prolongation. The table representation contains only artifact metadata with no preserved rows, and OCR is not needed for the genotype–phenotype counts.

- Missed gold variants: R176W
- Count disagreements: KCNH2 R176W carriers 16 vs 112 (error -96)

### KCNH2 PMID 26496715

**table** — The supplement’s structured mutation tables are the clearest available source for variant-level patient counts; the running text has no usable article body and only repeats truncated table content.

- Missed gold variants: A1017fsX, R1032fsX, D896fsX, C984fsX, G873fsX, L987fsX, P191fsX, T443fsX, T443fsX, P1034fsX, P926fsX, F494del, Q247X, R366X, A797T, R863X, W154X
- Extra predictions: c.461G>A (p.Trp154X), c.572delC (p.Pro191ArgX10), c.739C>T (p.Gln247X), c.1096C>T (p.Arg366X), c.1326delT (p.Thr443ProfsX78), c.1328_1329delCC/insG (p.Thr443ArgfsX77), c.1480_1482delTTC (p.Phe494del), c.2389G>A (p.Arg797Thr), c.2587C>T (p.Arg863X), c.2616delC (p.Gly873AlafsX5), c.2682_2685delCACG (p.Asp896ArgfsX77), c.2775_2776insG (p.Pro926AlaX14), c.2944_2948dupGACAC (p.Cys984ThrfsX75), c.2959-2960delCT (p.Leu987Valfs131X), c.3048delC (p.Ala1017ProfsX40), c.3094_3095insC (p.Arg1032AlafsX86), c.3102_3111delCCGGGGCGAC (p.Pro1034ArgfsX20)

### KCNH2 PMID 11844290

**text** — The running full text is the only available representation and includes the article body needed to capture family-level genotype and phenotype evidence for the five HERG/KCNH2 mutations; no separate table, PDF, or OCR representation is available.

- Missed gold variants: R752W, F805C, M124R, V822M, W1001X

### KCNH2 PMID 10973849

**text** — The reharvested full text is substantially complete and preserves the study narrative and embedded mutation tables. The parsed table artifact contains only limited HTML text and figure metadata, while OCR covers mutation-topology figures rather than variant-level carrier, affected, or unaffected counts.

- Missed gold variants: D864SP, K638del, L799SP, Q376SP, Q376SP
- Extra predictions: G1128A (Q376sp), A1129-2G (Q376sp), del1913–1915 (delK638), G2398+1C (L799sp), G2592+1A (D864sp)

### KCNH2 PMID 10862094

**text** — Text is the only available representation and is labeled as reharvested full text; it includes the KCNH2/HERG variant survey and associated cohort information. No table, PDF, or OCR representation is available to provide better-preserved variant-level carrier or phenotype counts.

- Missed gold variants: P151fsX, S543fsX
- Extra predictions: 453delC, 1631delAG, K897T (2690A>C)
- Count disagreements: G584S (1750G>A) carriers 19 vs 10 (error +9); R176W (526C>T) affected 1 vs 4 (error -3); G601S (1801G>A) affected 2 vs 1 (error +1); R176W (526C>T) unaffected 5 vs 2 (error +3); G584S (1750G>A) unaffected 15 vs 6 (error +9); G601S (1801G>A) unaffected 1 vs 2 (error -1)

### KCNH2 PMID 10841244

**text** — Text is the only available representation and includes variant-level evidence for HERG/KCNH2 L552S, including two homozygous affected siblings, heterozygous parents, 38 additional carriers, and symptomatic versus asymptomatic carrier groups.

- Count disagreements: T1655C (L552S; HERG-Fin) carriers 42 vs 44 (error -2)

### KCNH2 PMID 23864605

**text** — The complete PMC running text is available and is the clearest authoritative source for this cell-based functional study. The artifact preview contains experimental figure metadata rather than variant-level carrier, affected, or unaffected person counts, and OCR is unnecessary because the figures do not appear to contain pedigree evidence.

- Count disagreements: G601S affected 3 vs 0 (error +3); A614V affected 14 vs 0 (error +14); G601S unaffected 6 vs 0 (error +6); A614V unaffected 11 vs 0 (error +11)

### KCNH2 PMID 24667783

**table** — The structured patient-level table is the clearest available source for linking KCNH2 variants to individual carriers and associated phenotype data. The running text is fragmented and primarily provides cohort-level summaries, while the OCR inventory does not establish that figures contain necessary genotype–phenotype counts.

- Missed gold variants: A561V, T613M, G262fsX, G785fsX, L559H, M645V, N629S, P596L, P72R, R176W, R534C, R582C, S428X, S818L, A172V, E807X, G238R, G880V, L343fsX, R356H, S855R, S890C, W705fsX

### KCNH2 PMID 19038855

**text** — The PMC running full text is complete and directly reports KCNH2/LQT2 cohort counts and seizure phenotypes. The table preview does not expose structured variant-level rows or carrier/affected/unaffected counts, while OCR contains only aggregate figures and no necessary pedigree evidence.

- Missed gold variants: A193fsX, A558P, C64Y, D456Y, D501H, E698X, E876X, F640L, G1036fsX, G306W, G572S, G925fsX, I30T, I560fsX, I642del, M645L, N629I, P241L, P872fsX, Q676fsX, R176W, R252fsX, R366X, R534C, R582C, R73fsX, T613M, Y99S

### KCNQ1 PMID 19490272

**table** — The structured mutation table directly preserves variant-level KCNQ1 carrier counts (n) for the 54 missense mutations, while accompanying rows provide cardiac-event counts by conservation tertile. This is clearer and less ambiguous than the running text for person-count curation.

- No scored variant or count disagreement.

### KCNQ1 PMID 23153844

**table** — The preserved structured mutation table provides variant-level KCNQ1 patient counts for all 34 mutations, making it clearer and more reliable than running text for carrier-count curation.

- Extra predictions: M1V, M159sp, Y171X, L191fs/90, S225L, R243C, W305S, T312I, G314S, A341E, S349W, P400fs/62, Q530X

### KCNQ1 PMID 17470695

**table** — The parsed PMC artifact includes full text plus five structured tables and a supplement description, making it the best available representation for preserving variant-level subject counts. The figures provide only mutation locations and coarse frequency bins, not exact carrier, affected, and unaffected counts.

- Missed gold variants: M1V, T144A, A150FS/133[DELCT451-452], E160K, G168R, Y171X[513 C>G], R174H, A178P, Y184S, G185S, G189E, R190Q, L191FS/90[DELTGCGC572-576], R195FS/40[DELG585], S225L, A226V, R237P, D242N, R243C, V254M, R258C, R259C, L266P, G269S, L273F, I274V, S277L, G292D, F296S, G306R, T312I, G314S, Y315C, Y315S, P320H, T322M, G325R, DELF340[DELCTT1017-1019], A341E, A341V, P343S, S349W, S373P, P400FS/62[INSC1201-1202], I517T, R518X[1552C>T], M520R, V524G, Q530X[1588C>T], R562M, S566F, S571FS/20[DELC1714), R591H, R594Q, D611Y, A636FS/28[DELC1909]

### KCNQ1 PMID 14678125

**text** — The running full text contains the study population and KCNQ1 mutation-region carrier counts, whereas the table representation contains only artifact metadata and no preserved tables or variant-level rows.

- Missed gold variants: T144A, L151, G168R, Y171X, V172M, A178P, G185S, R190Q, S225L, R243C, V254L, V254M, L266P, G269D, G269S, L273F, E284K, A300T, G306R, T311I, T312I, Y315C, D317G, DELF340, A341E, A341V, A344/SP, A344A/SPLICE, G345R, S349W, Q357H, R366Q, R366W, K393N, R518X, V524G, Q530X, R539W, S566F, I567S, R594Q

### KCNQ1 PMID 28720088

**text** — The full PMC text explicitly reports KCNQ1 variant-level carrier counts for Y111C and R518* and genotype-negative individuals. The table preview contains only artifact metadata with no visible rows, while no pedigree or figure evidence requires OCR.

- Missed gold variants: R519X
- Extra predictions: KCNQ1 p.R518*

### KCNQ1 PMID 21129503

**text** — Text is the only available representation and explicitly reports the variant-level cohort count of 170 Y111C/KCNQ1 carriers across 37 proband families. No table, PDF, or OCR content is available to provide more complete carrier or phenotype counts.

- No scored variant or count disagreement.

### KCNQ1 PMID 25087618

**ocr** — Figure 1 contains the participant flow and numerical breakdown by KCNQ1 A341V carrier status and cardiac-event phenotype; its caption identifies these categories but omits the actual counts shown in the image. OCR is therefore necessary to preserve carrier, affected, and unaffected evidence.

- No scored variant or count disagreement.

### KCNQ1 PMID 17192539

**text** — Text is the only available representation and contains the full article narrative, including genotyped family structure, carrier cohorts, phenotyping criteria, and mutation-level study context.

- Missed gold variants: P117L, G168R, R174C, R174H, L175fsX, A178T, Y184S, R190W, R190Q, R190L, S225L, R231C, R243C, V254M, H258N, H258P, E261K, G269D, S277W, V280A, V280E, V288fsX, G306R, T309R, T311I, T312I, G314S, Y315S, Y315C, G316E, P320H, P320A, G325R, F339S, A341V, L342F, A344V, SP/A344A, Q359-K362DEL, R366W, A371T, S373P, W379S, K422fsX, M476L, M520R, R539W, R555C, R555H, S566P, I567T, T587M, A590T, R591H, R594C, Q604X

### KCNQ1 PMID 24052033

**text** — The running full text is available and includes the study methods and cohort descriptions for the KCNQ1 Y111C and R518X carriers. The table representation preserves only table captions, not the rows needed for variant-level carrier, affected, or unaffected counts.

- No scored variant or count disagreement.

### KCNQ1 PMID 18713323

**table** — The structured full-text artifacts include three tables and are most likely to preserve exact person counts across the six matched KCNQ1 variants and ethnic groups. The OCR figures are Kaplan-Meier curves and are less suitable for exact carrier or event counts.

- No scored variant or count disagreement.

### KCNQ1 PMID 29197658

**table** — Structured tables and supplementary artifacts are the best available source for the many KCNQ1 variant-level records and gnomAD individual/allele counts; the running text mainly reports aggregate results, while the figures summarize categories rather than person-level phenotype evidence.

- Extra predictions: F127L, P477L, L619M

### KCNQ1 PMID 33141630

**pdf** — The supplemental PDF contains the participant protocol and supplemental tables with carrier-level and phenotype-count evidence for KCNQ1 c.671C>T (p.T224M), including all 124 identified carriers and QTc classifications. These details are more complete than the running-text preview, while the table preview exposes mainly artifact metadata rather than table rows.

- Count disagreements: KCNQ1 c.671C>T (p.Thr224Met; p.T224M; rs199472706) carriers 124 vs 88 (error +36)

### RYR2 PMID 25814417

**table** — The preserved structured tables contain patient-level clinical fields and are best suited to retaining variant-level carrier and phenotype evidence for the RYR2 p.G357S founder variant.

- Count disagreements: p.G357S carriers 185 vs 179 (error +6); p.G357S affected 97 vs 73 (error +24); p.G357S unaffected 87 vs 106 (error -19)

### RYR2 PMID 29925740

**text** — The text contains the only substantive RYR2 carrier evidence, including cohort-level counts and nine RYR2-positive patients; the table representation contains metadata but no preserved variant-level rows or person counts.

- Missed gold variants: N4168S
- Extra predictions: c.12533A>G (p.N4178S)

### RYR2 PMID 33315912

**table** — The preserved structured table most clearly reports variant-carrier phenotype counts, including carriers with cardiac events (n=18) and without cardiac events (n=41), while the text provides broader cohort context.

- Count disagreements: P2328S affected 18 vs 15 (error +3)

### RYR2 PMID 16272262

**pdf** — The PDF contains the complete article and preserves the layout of the mutation and clinical tables, whereas the running text is largely abstract-only and the parsed table representation is severely fragmented. It is therefore the clearest source for variant-level carrier and phenotype evidence.

- Count disagreements: G4662S carriers 2 vs 5 (error -3); H4762P carriers 4 vs 5 (error -1); V4771I carriers 1 vs 2 (error -1); V4771I affected 1 vs 2 (error -1); G4662S unaffected 1 vs 4 (error -3); H4762P unaffected 3 vs 4 (error -1)

### RYR2 PMID 34202968

**ocr** — Variant-level carrier, affected, and unaffected evidence is likely encoded in the family pedigree figures (g001-g003), while the text provides only aggregate statements and the table artifact contains no person-level rows.

- No scored variant or count disagreement.

### RYR2 PMID 19926015

**text** — The running full text contains the complete cohort methods, variant findings, and clinical context. The table and PDF previews are limited mainly to a supplemental cross-species conservation table and do not preserve variant-level carrier, affected, or unaffected person counts; the figures do not appear necessary for genotype-phenotype evidence.

- Missed gold variants: L62F, M81L, V377M, Y2156C, E2183V, G4315E
- Extra predictions: exon 3 deletion (1.1 kb)
- Count disagreements: Y4149S carriers 2 vs 1 (error +1); Y4149S unaffected 1 vs 0 (error +1)

### RYR2 PMID 33686871

**ocr** — The variant evidence is family-segregation based and likely encoded in a pedigree figure. Running text reports four affected homozygotes and asymptomatic heterozygous relatives but does not preserve exact individual-level carrier and unaffected counts; no structured tables are present. OCR of the figure images is therefore the best source for genotype-phenotype counts.

- Missed gold variants: I3995V, D4112N, T4196I, D4646A, Q4879H, K4594R, I2075T

### RYR2 PMID 28237968

**pdf** — The PDF is the most complete authoritative article representation and preserves page-layout content, including mutation-level tables or supplementary references that the HTML/text extraction may omit. The parsed table artifact contains no actual table rows, and the OCR inventory provides no demonstrated pedigree-based genotype/phenotype evidence.

- Missed gold variants: E4451del, G2866del
- Extra predictions: NM_001035.2:c.5170G>A (p.Glu1724Lys), NM_001035.2:c.13350-13352del (p.Glu4451del), NM_001035.2:c.8596-8598del (p.Gly2866del), NM_001035.2:c.12470G>A (p.Arg4157Gln), NM_001035.2:c.14553C>A (p.Phe4851Leu), NM_001035.2:c.1258C>T (p.Arg420Trp), NM_001035.2:c.3407C>T (p.Ala1136Val), NM_001035.2:c.5656G>A (p.Gly1886Ser), NM_001035.2:c.5654G>A (p.Gly1885Glu)
- Count disagreements: NM_001035.2:c.7160C>T (p.Ala2387Val) carriers 1 vs 4 (error -3); NM_001035.2:c.7160C>T (p.Ala2387Val) unaffected 0 vs 3 (error -3)

### RYR2 PMID 12106942

**text** — The running full text is available and substantial, while the table representation contains only artifact metadata and no preserved table rows. Text therefore provides the clearest available source for variant-level carrier and phenotype counts.

- Count disagreements: R420W carriers 4 vs 8 (error -4); R420W unaffected 2 vs 6 (error -4)

### RYR2 PMID 22677073

**table** — The structured case-level table preserves gene, nucleotide and amino-acid changes, individual demographics, SUD event, sentinel phenotype, and family history, making it the clearest source for RYR2 variant-level carrier and phenotype evidence. OCR is unnecessary, and the running text provides mainly aggregate results.

- Missed gold variants: V2113M, Q2958R, A1136V, R4037C
- Extra predictions: RYR2 c.6739C>T (reported p.S2246L)
- Count disagreements: RYR2 c.6737C>T (p.S2246L) carriers 1 vs 2 (error -1); RYR2 c.6737C>T (p.S2246L) affected 1 vs 2 (error -1)

### RYR2 PMID 30403697

**table** — The PMC XML artifact includes two structured tables and the full article context; these tables are most likely to preserve patient-by-patient RYR2 variants and associated phenotype evidence more reliably than running text, PDF extraction, or structural figures.

- Missed gold variants: G4722S
- Extra predictions: p.G4772S
- Count disagreements: p.R2028H carriers 2 vs 3 (error -1); p.Y4721C carriers 2 vs 3 (error -1); c.3599-9delT carriers 1 vs 2 (error -1); c.14091-11dupT carriers 1 vs 2 (error -1); p.R2028H unaffected 1 vs 2 (error -1); p.Y4721C unaffected 1 vs 2 (error -1); c.3599-9delT unaffected 0 vs 1 (error -1); c.14091-11dupT unaffected 0 vs 1 (error -1)

### RYR2 PMID 33606749

**table** — The corrected structured Table 1 directly links each RYR2 variant to the proband, genotyped relatives, inheritance, and parental phenotypes, best preserving variant-level carrier and affected/unaffected evidence.

- Missed gold variants: N57_G91DEL35, 4944_4945INSH
- Extra predictions: exon 3 deletion (N57_G91del35), 14834_14835insTCA (4944_4945insH)
- Count disagreements: 1259g>a (R420Q) carriers 2 vs 1 (error +1); 3667a>g (T1223A) carriers 2 vs 1 (error +1); 3766c>a (P1256T) carriers 2 vs 1 (error +1); 4552c>t (L1518F) carriers 2 vs 1 (error +1); 5170g>a (E1724K) carriers 2 vs 1 (error +1); 6574a>t (M2192L) carriers 2 vs 1 (error +1); 7024g>a (G2342R) carriers 2 vs 1 (error +1); 7169c>t (T2390I) carriers 2 vs 1 (error +1); 11917g>a (D3973N) carriers 2 vs 1 (error +1); 12371g>a (S4124N) carriers 2 vs 1 (error +1); 14311g>a (V4771I) carriers 3 vs 2 (error +1); 9910c>g (Q3304E) carriers 2 vs 1 (error +1); 14222c>t (A4741V) carriers 2 vs 1 (error +1); 1259g>a (R420Q) affected 2 vs 1 (error +1); 5170g>a (E1724K) affected 2 vs 1 (error +1); 14311g>a (V4771I) affected 3 vs 2 (error +1); 9910c>g (Q3304E) affected 2 vs 1 (error +1); 14222c>t (A4741V) affected 2 vs 1 (error +1); 3667a>g (T1223A) unaffected 1 vs 0 (error +1); 3766c>a (P1256T) unaffected 1 vs 0 (error +1); 4552c>t (L1518F) unaffected 1 vs 0 (error +1); 6574a>t (M2192L) unaffected 1 vs 0 (error +1); 7024g>a (G2342R) unaffected 1 vs 0 (error +1); 7169c>t (T2390I) unaffected 1 vs 0 (error +1); 11917g>a (D3973N) unaffected 1 vs 0 (error +1); 12371g>a (S4124N) unaffected 1 vs 0 (error +1)

## Scope, method, and limitations

- Population: the fixed 20-paper comparison set, five papers per gene; every PMID has downloaded source and at least one gold assertion in each count field.
- Blinding: selection exported no gold values or row counts; each paper was routed and read from local source only; predictions were made read-only and SHA-256 locked before `score` opened gold.
- Variant metrics are micro-averaged over gold rows. Precision treats unmatched predictions as false positives, although the curated recall packet may omit some real variants.
- Count MAE/RMSE are conditional on a supplied value. Count recall must be read alongside them because abstentions and missed variants are excluded from error magnitude.
- Source acquisition and gold completeness are separate from model reading quality; abstract-only or incomplete source is retained and labeled rather than silently excluded.

## Reproducibility and evidence

- `selection.json`: selected PMIDs, source paths, source hashes, and available representations.
- `predictions.json`: immutable per-paper tools, rationales, extracted variants, counts, evidence quotes, source locations, elapsed time, and token telemetry.
- `evidence.csv`: flat evidence ledger for every predicted variant.
- `paper_metrics.csv`: exact per-paper metrics.
- `LOCK.json`: SHA-256 digests proving prediction finalization before scoring.
- `report.json`: complete machine-readable score, errors, timing, and token usage.

## Recommended next steps

1. Adjudicate extra predictions against the paper before treating precision as a production false-positive rate.
2. Review count outliers by source location and distinguish model mistakes from gold disagreements.
3. Repeat on the high-carrier 48-paper stress set only after preserving the same lock and count-recall definitions.

## Further questions

- Does table/PDF/OCR routing improve recall enough to justify its additional routing-call tokens?
- How much of the residual error is source incompleteness versus count-role interpretation?
