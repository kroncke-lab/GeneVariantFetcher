# Codex extraction-blinded paper evaluation — `20260724_fixed48_grok`

## Technical summary

This hash-locked run evaluated **48 papers** (**12 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **90.8%**, recall **59.8%**, F1 **72.1%** (599 TP, 61 FP, 402 FN).
- Exact API telemetry: **945,426 total tokens** (687,231 input; 258,195 output).
- Elapsed: **1807.1s wall clock**; 2703.5s summed per-paper route + read time.
- Representation choices: {'text': 24, 'table': 22, 'pdf': 2}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `highcarrier_48_papers_20260723.tsv` (48 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- `selection.json` contains source metadata and hashes but no gold values or gold row counts. `predictions.json` was made read-only and SHA-256 locked before scoring first opened the gold CSVs.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 520 / 1001 | 51.9% | 0.612 | 5.860 |
| affected | 211 / 1001 | 21.1% | 0.706 | 4.958 |
| unaffected | 68 / 1001 | 6.8% | 0.985 | 3.939 |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 129 | 9 | 91 | 93.5% | 58.6% | 72.1% | 57.3% / 1.095 / 7.312 | 17.3% / 3.079 / 11.421 | 2.7% / 4.833 / 11.839 |
| KCNH2 | 141 | 3 | 151 | 97.9% | 48.3% | 64.7% | 29.1% / 1.447 / 10.736 | 1.7% / 3.400 / 6.403 | 1.7% / 3.400 / 5.604 |
| KCNQ1 | 170 | 38 | 101 | 81.7% | 62.7% | 71.0% | 62.4% / 0.213 / 2.769 | 10.7% / 0.000 / 0.000 | 0.0% / n/a / n/a |
| RYR2 | 159 | 11 | 59 | 93.5% | 72.9% | 82.0% | 64.2% / 0.150 / 0.455 | 63.8% / 0.108 / 0.424 | 26.1% / 0.368 / 1.000 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 27566755 | text | 51 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 84.8 | 22,163 |
| SCN5A | 26669661 | text | 26 | 1 | 1 | 96.3% | 96.3% | 96.3% | 96.3% / 0.000 | 0.0% / n/a | 0.0% / n/a | 58.8 | 58,178 |
| SCN5A | 20470418 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 59.000 | 100.0% / 22.000 | 0.0% / n/a | 38.4 | 10,410 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 55.000 | 100.0% / 55.000 | 0.0% / n/a | 34.5 | 16,913 |
| SCN5A | 29709101 | table | 10 | 2 | 1 | 83.3% | 90.9% | 87.0% | 90.9% / 0.000 | 0.0% / n/a | 0.0% / n/a | 46.4 | 16,752 |
| SCN5A | 28341781 | text | 4 | 0 | 51 | 100.0% | 7.3% | 13.6% | 7.3% / 0.000 | 7.3% / 0.000 | 7.3% / 0.000 | 57.5 | 20,896 |
| SCN5A | 18451998 | text | 1 | 1 | 1 | 50.0% | 50.0% | 50.0% | 50.0% / 9.000 | 50.0% / 38.000 | 50.0% / 29.000 | 41.3 | 14,428 |
| SCN5A | 26921764 | table | 27 | 1 | 0 | 96.4% | 100.0% | 98.2% | 100.0% / 0.074 | 100.0% / 0.074 | 0.0% / n/a | 81.6 | 16,194 |
| SCN5A | 27554632 | table | 3 | 0 | 6 | 100.0% | 33.3% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 34.5 | 12,696 |
| SCN5A | 10590249 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 47.6 | 17,923 |
| SCN5A | 25051102 | table | 2 | 4 | 1 | 33.3% | 66.7% | 44.4% | 66.7% / 0.500 | 66.7% / 0.000 | 0.0% / n/a | 56.7 | 19,423 |
| SCN5A | 26746457 | text | 2 | 0 | 30 | 100.0% | 6.2% | 11.8% | 6.2% / 6.000 | 6.2% / 0.000 | 3.1% / 0.000 | 48.1 | 12,607 |
| KCNH2 | 29622001 | table | 34 | 1 | 1 | 97.1% | 97.1% | 97.1% | 97.1% / 0.735 | 0.0% / n/a | 0.0% / n/a | 81.5 | 23,228 |
| KCNH2 | 11854117 | table | 44 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 82.0 | 15,964 |
| KCNH2 | 14661677 | text | 0 | 0 | 29 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.6 | 20,726 |
| KCNH2 | 19160088 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 66.7% / 48.000 | 0.0% / n/a | 0.0% / n/a | 41.4 | 11,928 |
| KCNH2 | 26496715 | table | 0 | 0 | 54 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 23.3 | 19,047 |
| KCNH2 | 11844290 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 19.1 | 5,888 |
| KCNH2 | 10973849 | table | 52 | 0 | 8 | 100.0% | 86.7% | 92.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 197.3 | 80,492 |
| KCNH2 | 10862094 | text | 6 | 2 | 2 | 75.0% | 75.0% | 75.0% | 25.0% / 0.000 | 25.0% / 0.000 | 25.0% / 0.000 | 71.8 | 15,927 |
| KCNH2 | 10841244 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 100.0% / 0.000 | 100.0% / 0.000 | 43.2 | 14,410 |
| KCNH2 | 23864605 | table | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 8.500 | 100.0% / 8.500 | 30.0 | 12,351 |
| KCNH2 | 24667783 | table | 0 | 0 | 23 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 26.1 | 7,951 |
| KCNH2 | 19038855 | text | 0 | 0 | 28 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 18.0 | 15,064 |
| KCNQ1 | 19490272 | table | 54 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 106.2 | 15,554 |
| KCNQ1 | 23153844 | table | 21 | 13 | 0 | 61.8% | 100.0% | 76.4% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 58.9 | 14,567 |
| KCNQ1 | 17470695 | table | 53 | 24 | 3 | 68.8% | 94.6% | 79.7% | 94.6% / 0.000 | 0.0% / n/a | 0.0% / n/a | 198.0 | 31,987 |
| KCNQ1 | 14678125 | text | 0 | 0 | 41 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 17.6 | 8,255 |
| KCNQ1 | 28720088 | pdf | 1 | 1 | 1 | 50.0% | 50.0% | 50.0% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 38.1 | 20,440 |
| KCNQ1 | 21129503 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 25.3 | 12,504 |
| KCNQ1 | 25087618 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 31.5 | 15,770 |
| KCNQ1 | 17192539 | text | 1 | 0 | 56 | 100.0% | 1.8% | 3.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.0 | 9,721 |
| KCNQ1 | 24052033 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 46.1 | 16,611 |
| KCNQ1 | 18713323 | table | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 44.9 | 14,119 |
| KCNQ1 | 29197658 | table | 29 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 73.5 | 16,605 |
| KCNQ1 | 33141630 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 36.000 | 0.0% / n/a | 0.0% / n/a | 45.8 | 36,025 |
| RYR2 | 25814417 | table | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 49.9 | 12,616 |
| RYR2 | 29925740 | text | 50 | 1 | 1 | 98.0% | 98.0% | 98.0% | 98.0% / 0.000 | 98.0% / 0.000 | 0.0% / n/a | 94.2 | 23,000 |
| RYR2 | 33315912 | table | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 3.000 | 100.0% / 6.000 | 47.8 | 11,277 |
| RYR2 | 16272262 | table | 12 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 45.0 | 17,409 |
| RYR2 | 34202968 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 30.7 | 22,365 |
| RYR2 | 19926015 | table | 0 | 0 | 40 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 42.7 | 13,831 |
| RYR2 | 33686871 | pdf | 1 | 0 | 7 | 100.0% | 12.5% | 22.2% | 0.0% / n/a | 12.5% / 2.000 | 0.0% / n/a | 36.8 | 18,930 |
| RYR2 | 28237968 | text | 18 | 9 | 0 | 66.7% | 100.0% | 80.0% | 100.0% / 0.167 | 100.0% / 0.111 | 100.0% / 0.278 | 91.0 | 48,646 |
| RYR2 | 12106942 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 16.7% / 0.000 | 0.0% / n/a | 0.0% / n/a | 38.4 | 17,439 |
| RYR2 | 22677073 | table | 19 | 0 | 7 | 100.0% | 73.1% | 84.4% | 73.1% / 0.000 | 73.1% / 0.000 | 0.0% / n/a | 73.1 | 18,061 |
| RYR2 | 30403697 | table | 18 | 1 | 3 | 94.7% | 85.7% | 90.0% | 85.7% / 0.111 | 85.7% / 0.000 | 28.6% / 0.333 | 72.2 | 31,139 |
| RYR2 | 33606749 | table | 32 | 0 | 1 | 100.0% | 97.0% | 98.5% | 97.0% / 0.500 | 97.0% / 0.250 | 97.0% / 0.250 | 75.2 | 16,996 |

## Errors and representation choices

### SCN5A PMID 27566755

**text** — Running full text supplies the clearest narrative description of the SCN5A/LQT3 cohort, patient counts, mutation locations, QTc-event relationships, and beta-blocker outcomes that support blinded variant-level carrier/affected/unaffected evidence extraction.

- No scored variant or count disagreement.

### SCN5A PMID 26669661

**text** — Text preview contains all mentions of SCN5A/LQT3 with aggregated maternal/paternal transmission ratios and family-level counts; available table preview is restricted to KCNQ1 variants and supplies zero SCN5A variant-level carrier or transmission data.

- Missed gold variants: c.703+1G>A
- Extra predictions: p.0?

### SCN5A PMID 20470418

**text** — Full main_text (PMC XML, 16025 chars) contains the complete methods, results, and counts for the single SCN5A variant p.Ser1103Tyr (26 heterozygotes, 13 with allelic imbalance, bimodal log2 ratios, SIDS vs. other-cause controls), directly preserving all carrier/affected/unaffected evidence at the variant level; no tables exist, PDF is unavailable, and the supplement jpg is a non-pedigree expression plot unnecessary for genotype/phenotype curation.

- Count disagreements: p.Ser1103Tyr carriers 26 vs 85 (error -59); p.Ser1103Tyr affected 17 vs 39 (error -22)

### SCN5A PMID 28339995

**text** — Running full text supplies the clearest narrative on D1790G SCN5A carriers (20 total), compliance status, and counts of cardiac events in affected/unaffected individuals; table/pdf/ocr previews are empty and unavailable.

- Count disagreements: D1790G carriers 30 vs 85 (error -55); D1790G affected 30 vs 85 (error -55)

### SCN5A PMID 29709101

**table** — Table supplies explicit variant-level rows with BrS+ (affected) vs BrS- (carrier/unaffected) patient counts plus ACMG classifications; text supplies only aggregate totals and no per-variant breakdown.

- Missed gold variants: c.4813+3_4813_6dup
- Extra predictions: c.4813+3_4813+6dup, c.4978A>G (p.(Ile1660Val))

### SCN5A PMID 28341781

**text** — Running full text is the clearest source and the only representation supplied; it contains the complete proband registry narrative on SCN5A mutations, cardiac events, and genotype-phenotype data without reliance on absent tables, PDFs or figures.

- Missed gold variants: A1428S, A735V, c.1537delC, c.3840+1G>A, c.4245+1G>A, c.4299+1delG, c.4389_4396delCCTCTTTA, c.4732_4733dupAA, c.5157delC, c.5290delG, c.870delC, E428K, F532C, G1408R, G1420R, G1743R, H1200Y, H1923D, K1527R, L276Q, L846R, N406S, P.N1380DEL, Q55X, Q779K, R104Q, R1195H, R1432S, R1623X, R1644C, R1644H, R179X, R1913C, R1919H, R367H, R376H, R689H, R693C, R814Q, R878C, R893C, R893H, R988Q, T1247I, T1709M, V1328M, V1405M, W1095X, W904X, Y1409X, Y416C

### SCN5A PMID 18451998

**text** — Figure 1 pedigrees supply per-individual variant-level carrier/affected/unaffected symbols, QTc values, sinus node dysfunction, BrS ST elevation, ICD/PPM, and SCD data that textual summaries explicitly omit. Requested ocr, which was unavailable; used text.

- Missed gold variants: T1304M
- Extra predictions: V1098L
- Count disagreements: E1784K carriers 41 vs 50 (error -9); E1784K affected 9 vs 47 (error -38); E1784K unaffected 32 vs 3 (error +29)

### SCN5A PMID 26921764

**table** — Table supplies explicit per-variant patient counts and lists of SCN5A mutations (e.g., 9 carriers of E1784K, 3 of S1672Y) that directly encode carrier-level evidence; narrative text only reports aggregate genotype-positive totals without variant resolution or individual carrier/phenotype mapping.

- Extra predictions: variant notation unavailable
- Count disagreements: 5228 G>A (G1743E) carriers 4 vs 2 (error +2); 5228 G>A (G1743E) affected 4 vs 2 (error +2)

### SCN5A PMID 27554632

**table** — Table rows provide structured variant-level counts of carriers across DCM/HCM patients and controls (affected/unaffected), which text only summarizes narratively; this directly preserves the required carrier/phenotype evidence for SCN5A variants like T1247I/A1260D/G1262S.

- Missed gold variants: M1245I, R513C, R513H, R526H, T512I, T570N

### SCN5A PMID 10590249

**text** — Pedigree figure (F1) plus variant confirmation figure (F2) supply the only complete variant-level mapping of 1795insD carriers to affected/unaffected status, obligate carriers, nocturnal deaths, and ECG phenotypes; running text and captions alone omit these person-level genotype/phenotype linkages. Requested ocr, which was unavailable; used text.

- No scored variant or count disagreement.

### SCN5A PMID 25051102

**table** — Table previews directly list SCN5A variants (e.g., 87G>A, G400A, H558R, R1193Q) with PCR-SSP typing details and tie to baseline patient/control groupings (AMI/VF+ n=49, AMI/VF- n=74, controls n=480) that contain the carrier/affected/unaffected counts.

- Missed gold variants: c.2788-6C>T
- Extra predictions: 87G>A, 1673A>G (His558Arg), IVS16-6C>T, 5457T>A
- Count disagreements: 3578G>A (Arg1193Gln) carriers 1 vs 2 (error -1)

### SCN5A PMID 26746457

**text** — Running full text supplies the clearest source of cohort counts, carrier vs non-carrier phenotype comparisons, and ECG/ICD-9 arrhythmia data for SCN5A variant carriers.

- Missed gold variants: A1180V, A1270S, A1680T, A551V, C982R, D1275N, E1053K, E48K, F1596I, F532C, G1935S, G552W, G615E, I1836T, I848F, L1194M, L1704H, P717L, R1023H, R1195S, R1512W, R1739Q, R1898C, R2012C, R282C, R689H, S1904L, T1304M, V1353M, V1532I
- Count disagreements: SCN5A-R1193Q carriers 19 vs 7 (error +12)

### KCNH2 PMID 29622001

**table** — Table rows directly supply variant-level carrier counts (n, families) for specific KCNH2 mutations (R176W, L552S) plus mutation type/location, which best preserves the needed person-level evidence over narrative risk associations in text.

- Missed gold variants: P1034fsX
- Extra predictions: c.3093_3106del
- Count disagreements: p.L552S carriers 73 vs 74 (error -1); c.453delC carriers 24 vs 0 (error +24)

### KCNH2 PMID 11854117

**table** — structured table rows carry the variant-level person counts (No. of Subjects) with mutation-specific entries, subtotals for pore/non-pore regions, and direct variant-level carrier evidence

- No scored variant or count disagreement.

### KCNH2 PMID 14661677

**text** — Text supplies no usable article body (only supplement recovery). Table and PDF previews contain only unrelated impact-factor/reviewer metadata with zero variant, carrier, affected, or unaffected counts. OCR targets the figure/pedigree images (fig_p1_1.png, fig_p1_2.png, Figure 8) that alone retain genotype-phenotype evidence required for variant-level curation. Requested ocr, which was unavailable; used text.

- Missed gold variants: A1058E, A190T, A203T, A915V, C723R, G187del, G187S, G873S, H254Q, K897T, K897T, K897T, K897T, L1023del, N257H, N33T, P251A, P347S, P910L, P917L, P967L, Q1068R, R1035W, R1047L, R1047L, R176W, R181Q, T367S, V215G

### KCNH2 PMID 19160088

**text** — Text preview alone supplies explicit variant-level carrier counts (KCNH2 R176W n=16 plus QT prolongation statistics) while table/PDF/OCR previews are empty or unavailable; main-text extraction is therefore the sole source preserving the required carrier/affected evidence.

- Missed gold variants: R176W
- Count disagreements: KCNH2 R176W carriers 16 vs 112 (error -96)

### KCNH2 PMID 26496715

**table** — Table supplies structured variant rows with explicit No. of patients counts per allele; text preview is limited to repeated, truncated supplement fragments lacking main-body context and offers no additional carrier/affected/unaffected breakdowns.

- Missed gold variants: A1017fsX, R1032fsX, D896fsX, C984fsX, G873fsX, L987fsX, P191fsX, T443fsX, T443fsX, P1034fsX, P926fsX, c.1557+1G>A, A169G, A32T, A561T, A561V, A614V, D102H, D609G, D837Y, E480V, E544A, E575K, F494del, F805S, G572S, G601S, G604S, G626S, G628S, G71W, I782N, L413P, L457P, L559H, N629S, N633S, P596S, P72T, P926S, Q247X, R366X, R534C, A797T, R823W, R863X, R912W, S620I, S818L, S818W, T613M, T634P, W154X, Y475C

### KCNH2 PMID 11844290

**text** — running full text is the clearest source and supplies the variant-level carrier/affected/unaffected evidence across the five HERG (KCNH2) mutations and families described in the title and body

- Missed gold variants: R752W, F805C, M124R, V822M, W1001X

### KCNH2 PMID 10973849

**table** — Structured table rows with nucleotide changes, coding effects, No.of counts, plus TABLE1 genotype rows showing age/sex/QTc/symptoms % per gene (incl. HERG/KCNH2) directly preserve variant-level carrier and affected/unaffected phenotypic evidence; running text disperses this data.

- Missed gold variants: A561T, A561V, D864SP, G604S, G628S, L799SP, Q376SP, T613M

### KCNH2 PMID 10862094

**text** — Only available representation; full-text extraction (426 lines) contains complete mutation survey with evolutionarily conserved HERG regions, specific variants (P451L/Y569H/etc.), and LQTS patient counts needed for variant-level carrier/affected/unaffected evidence.

- Missed gold variants: P151fsX, S543fsX
- Extra predictions: 453delC, 1631delAG

### KCNH2 PMID 10841244

**text** — Text provides the only available representation and directly preserves variant-level evidence including L552S homozygotes (2 affected siblings with specific phenotypes: 2:1 AV block, torsades, death), heterozygotes (parents + 38 additional subjects across 6 families), plus aggregate counts and stats by symptomatic carriers (mean QTc 500 ms), asymptomatic carriers (452 ms), and noncarriers (412 ms).

- Count disagreements: L552S carriers 42 vs 44 (error -2)

### KCNH2 PMID 23864605

**table** — Table preview directly supplies variant-level counts (families/subjects/symptomatic) for G601S and A614V in KCNH2, preserving carrier/affected evidence lost from textual biochemical assays.

- Count disagreements: G601S affected 3 vs 0 (error +3); A614V affected 14 vs 0 (error +14); G601S unaffected 6 vs 0 (error +6); A614V unaffected 11 vs 0 (error +11)

### KCNH2 PMID 24667783

**table** — structured table rows carry the variant-level person counts with gene/exon/protein region, zygosity, and phenotype stats (QTc, syncope, SCD) that preserve carrier/affected/unaffected evidence for KCNH2

- Missed gold variants: A561V, T613M, G262fsX, G785fsX, L559H, M645V, N629S, P596L, P72R, R176W, R534C, R582C, S428X, S818L, A172V, E807X, G238R, G880V, L343fsX, R356H, S855R, S890C, W705fsX

### KCNH2 PMID 19038855

**text** — Full-text extraction supplies complete methods/results on the LQT2 seizure-phenotype aggregate (77 LQT2 subjects, 30/77 positive phenotype, p-values vs other genotypes) plus background genotype-negative rates; no per-variant carrier/affected counts or pedigrees exist in any supplied preview.

- Missed gold variants: A193fsX, A558P, C64Y, D456Y, D501H, E698X, E876X, F640L, G1036fsX, G306W, G572S, G925fsX, I30T, I560fsX, I642del, M645L, N629I, P241L, P872fsX, Q676fsX, R176W, R252fsX, R366X, R534C, R582C, R73fsX, T613M, Y99S

### KCNQ1 PMID 19490272

**table** — Table rows enumerate each KCNQ1 missense mutation with its exact patient count (n), enabling direct variant-level carrier, affected, and unaffected tallies by conservation tertile; running text aggregates only at tertile level and omits per-variant n.

- No scored variant or count disagreement.

### KCNQ1 PMID 23153844

**table** — Structured table rows explicitly list 34 KCNQ1 mutations with per-variant patient counts (e.g., V254M:118, G168R:87), directly supplying the variant-level carrier/person evidence needed for curation; main text only references an Online Table and supplies context without these counts.

- Extra predictions: M1V, M159sp, Y171X, L191fs/90, S225L, R243C, W305S, T312I, G314S, A341E, S349W, P400fs/62, Q530X

### KCNQ1 PMID 17470695

**table** — Table rows explicitly list each KCNQ1 variant with per-variant subject counts (No. of Subjects), mutation type, and functional effect; this directly supplies the required variant-level carrier evidence that textual narrative summaries omit.

- Missed gold variants: V254M, T322M, R562M
- Extra predictions: G57V, W120C, G189R, V254 mol/L, G269D, Y278H, E284K, V310I, D317G, T322 mol/L, A344A/sp [1032 G>A], A344V, L353P, Q357H, R360G, K393N, R397W, P448fs/13 [ins G 1344-1345], R562 mol/L, I567S, R591C, IVS2+1 G>A, IVS4+5 G>A, IVS7+5 G>A

### KCNQ1 PMID 14678125

**text** — Main text supplies aggregate counts and phenotypic data for all 294 KCNQ1-mutation carriers (164 pre-pore, 101 pore, 29 post-pore) together with QTc values and cardiac-event follow-up; no tables exist and figure captions are uninformative for genotype-phenotype linkage.

- Missed gold variants: T144A, L151, G168R, Y171X, V172M, A178P, G185S, R190Q, S225L, R243C, V254L, V254M, L266P, G269D, G269S, L273F, E284K, A300T, G306R, T311I, T312I, Y315C, D317G, DELF340, A341E, A341V, A344/SP, A344A/SPLICE, G345R, S349W, Q357H, R366Q, R366W, K393N, R518X, V524G, Q530X, R539W, S566F, I567S, R594Q

### KCNQ1 PMID 28720088

**pdf** — PDF preserves table structure and variant-level counts (Y111C n=148, R518* n=79, genotype-negative cohort) that markdown text extraction omits despite 4 table captions; best matches carrier/affected/unaffected evidence preference.

- Missed gold variants: R519X
- Extra predictions: R518*

### KCNQ1 PMID 21129503

**text** — Text preview supplies main-text abstract with explicit variant-level carrier counts (170 Y111C/KCNQ1 carriers in 37 families), genealogic tracing to founder couple, and haplotype allele-sharing data that directly capture carrier/affected evidence; table/PDF/OCR previews absent.

- No scored variant or count disagreement.

### KCNQ1 PMID 25087618

**text** — Text supplies the explicit study population counts (181 non-carriers, 168 KCNQ1-A341V carriers) plus genotype/haplotype-level associations with QTc, cardiac events and severity; no tables, PDFs or figure-derived genotype/phenotype counts are supplied in the catalog.

- No scored variant or count disagreement.

### KCNQ1 PMID 17192539

**text** — Running full text is the clearest source and preserves all details on 484 nuclear families, 142 mutations, 480 KCNQ1 carriers, phenotype-genotype correlations, clinical symptoms, and transmission in pedigrees with affected/unaffected evidence.

- Missed gold variants: P117L, G168R, R174C, R174H, L175fsX, A178T, Y184S, R190W, R190Q, R190L, S225L, R231C, R243C, V254M, H258N, H258P, E261K, G269D, S277W, V280A, V280E, V288fsX, G306R, T309R, T311I, T312I, G314S, Y315S, Y315C, G316E, P320H, P320A, G325R, F339S, A341V, L342F, A344V, SP/A344A, Q359-K362DEL, R366W, A371T, S373P, W379S, K422fsX, M476L, M520R, R539W, R555C, R555H, S566P, I567T, T587M, A590T, R591H, R594C, Q604X

### KCNQ1 PMID 24052033

**text** — running full text is the clearest source and directly describes the 150 individuals, 22 LQTS families, Y111C and R518X carrier phenotypes, and affected/unaffected comparisons without reliance on omitted figure or table content

- No scored variant or count disagreement.

### KCNQ1 PMID 18713323

**table** — Table rows supply per-mutation carrier counts (e.g., A341V: 4 Caucasian/16 Japanese) plus biophysical classification, directly preserving variant-level person-level evidence needed for carrier/affected/unaffected mapping; running text supplies only aggregate totals.

- No scored variant or count disagreement.

### KCNQ1 PMID 29197658

**table** — Table rows explicitly list per-variant 'No. of cases' together with gnomAD allele counts, directly preserving variant-level carrier/person counts for affected cases; narrative text supplies only aggregate summaries and functional results for three variants.

- No scored variant or count disagreement.

### KCNQ1 PMID 33141630

**text** — Main text directly states 124 carriers identified among 5521 Amish participants plus quantitative QTc association (20.2 ms shift, boxplot comparison of carriers vs non-carriers), supplying the core variant-level carrier and phenotype counts; table preview contains only external gnomAD variant data, not study-specific carrier/affected counts; pdf preview is limited to supplemental methods and recruitment protocol without the numerical evidence.

- Count disagreements: p.T224M (c.671C>T) carriers 124 vs 88 (error +36)

### RYR2 PMID 25814417

**table** — structured table rows carry the variant-level person counts and list individual mutation carriers with phenotype data (symptoms, VA/CVA test results) that best preserve carrier/affected/unaffected evidence for p.G357S_RyR2

- No scored variant or count disagreement.

### RYR2 PMID 29925740

**text** — Text is the only available representation and contains patient counts (9 of 104 LQTS-diagnosed cases) plus RYR2 mutation context for carrier/affected status; no tables, figures, or PDF content supplied.

- Missed gold variants: N4168S
- Extra predictions: c.12533A>G (p.N4178S)

### RYR2 PMID 33315912

**table** — Structured table rows directly preserve variant-level carrier counts (62 total RYR2 P2328S carriers) stratified by affected (CE n=18) vs unaffected (No CE n=41) status, with explicit person-level data on probands, ACA/syncope events, family history, β-blocker/ICD use, and exercise-stress PVC findings that text abstracts omit.

- Count disagreements: p.P2328S affected 18 vs 15 (error +3); p.P2328S unaffected 41 vs 47 (error -6)

### RYR2 PMID 16272262

**table** — Structured table rows (via table 2 / genotype-phenotype section) directly preserve variant-level carrier, affected, and unaffected counts plus low-penetrance silent carriers, as referenced in snippets naming fig 1/table 2 and explicitly stating 'Eleven RYR2 mutation carriers were considered to be phenotypically unaffected'.

- No scored variant or count disagreement.

### RYR2 PMID 34202968

**text** — Main text preview supplies full running-text details on the p.Asp3291Val variant across three families, including eight sudden deaths before age 30, all affected subjects carrying at least one copy, and three homozygous affected sisters, which directly supplies variant-level carrier/affected counts and phenotypes.

- No scored variant or count disagreement.

### RYR2 PMID 19926015

**table** — structured table rows carry the variant-level person counts

- Missed gold variants: L62F, M81L, E243K, F329L, R332W, G357S, V377M, T415R, R420Q, V507I, A549V, R739H, R1013Q, T1107M, A1136V, E1837K, E2045G, Y2156C, H2168Q, E2183V, D2216V, EXON 3 DELETION, E2296Q, R2420W, M3972I, D3973H, S4124G, Y4149S, R4157Q, Q4159P, N4178S, E4187Q, G4315E, K4650E, N4736 DEL, R4790Q, K4805R, R4822H, L3974Q, K3997E

### RYR2 PMID 33686871

**pdf** — PDF-layout text is more complete and preserves family pedigree/genotype details (homozygous affected vs heterozygous unaffected) that running-text extraction truncates or fragments.

- Missed gold variants: I3995V, D4112N, T4196I, D4646A, Q4879H, K4594R, I2075T
- Count disagreements: G3118R affected 4 vs 6 (error -2)

### RYR2 PMID 28237968

**text** — Running full text (41k chars) supplies all aggregate and per-group counts for probands vs relatives (affected/unaffected) plus mutation descriptions; no tables exist to preserve and figure captions contain only Kaplan-Meier or schematic data without carrier counts.

- Extra predictions: c.5170G>A (p.Glu1724Lys), c.12470G>A (p.Arg4157Gln), c.14553C>A (p.Phe4851Leu), c.1258C>T (p.Arg420Trp), c.3407C>T (p.Ala1136Val), c.5656G>A (p.Gly1886Ser), c.3407C>T (p.Ala1136Val) [repeat count already aggregated], c.5654G>A (p.Gly1885Glu), c.5656G>A (p.Gly1886Ser) [repeat count already aggregated]
- Count disagreements: c.7160C>T (p.Ala2387Val) carriers 1 vs 4 (error -3); c.94C>A (p.Gln32Lys) affected 4 vs 5 (error -1); c.3247G>A (p.Glu1083Lys) affected 4 vs 5 (error -1); c.7160C>T (p.Ala2387Val) unaffected 0 vs 3 (error -3); c.94C>A (p.Gln32Lys) unaffected 1 vs 0 (error +1); c.3247G>A (p.Glu1083Lys) unaffected 1 vs 0 (error +1)

### RYR2 PMID 12106942

**text** — Full main text (56701 chars) supplies aggregate carrier counts (43 mutation-positive subjects), symptomatic status (28 affected), sudden death outcomes, and family-level phenotype correlations for the six RYR2 mutations; no tables, figures, or OCR material exist to supply superior variant-level carrier/affected/unaffected breakdowns.

- No scored variant or count disagreement.

### RYR2 PMID 22677073

**table** — Table rows explicitly list individual cases with RYR2/CPVT1 variants (20 positive cases), plus per-case fields for sex/age/ethnicity/SUD event/family history that directly encode carrier status, affected/unaffected evidence, and genotype-phenotype links; running text supplies only aggregate yields and screening methods without these variant-level details.

- Missed gold variants: T1107M, V2113M, G1885E, G1886S, Q2958R, A1136V, R4037C

### RYR2 PMID 30403697

**table** — Table rows explicitly enumerate each of the 15 multi-variant subjects with their exact RYR2 (or RYR2+CASQ2) alleles, cis/trans phase, inheritance, family-history carrier status (e.g., asymptomatic heterozygous parents/siblings), and phenotype (symptoms/SCA/affected vs. asymptomatic/unaffected), directly supplying the variant-level carrier/affected/unaffected counts required for curation.

- Missed gold variants: G4722S, c.3599-9delT, c.14091-11dupT
- Extra predictions: RYR2-p.G4772S
- Count disagreements: RYR2-p.R2028H carriers 2 vs 3 (error -1); RYR2-p.Y4721C carriers 2 vs 3 (error -1); RYR2-p.R2028H unaffected 1 vs 2 (error -1); RYR2-p.Y4721C unaffected 1 vs 2 (error -1)

### RYR2 PMID 33606749

**table** — Table rows directly enumerate each proband's specific RYR2 variant (nucleotide/amino-acid change), inheritance origin (de novo/maternal/paternal), and explicit father/mother phenotypes (none/syncope/AF), supplying the precise variant-level carrier/affected/unaffected counts required; text supplies only aggregate statistics, PDF is limited to a correction notice, and OCR images are unnecessary when the structured table already preserves the evidence.

- Missed gold variants: Q3861H
- Count disagreements: N57_G91del35 carriers 4 vs 2 (error +2); R420Q carriers 2 vs 1 (error +1); T1223A carriers 2 vs 1 (error +1); P1256T carriers 2 vs 1 (error +1); L1518F carriers 2 vs 1 (error +1); E1724K carriers 2 vs 1 (error +1); M2192L carriers 2 vs 1 (error +1); G2342R carriers 2 vs 1 (error +1); T2390I carriers 2 vs 1 (error +1); Q3861H carriers 2 vs 1 (error +1); D3973N carriers 2 vs 1 (error +1); S4124N carriers 2 vs 1 (error +1); V4771I carriers 3 vs 2 (error +1); Q3304E carriers 2 vs 1 (error +1); A4741V carriers 2 vs 1 (error +1); N57_G91del35 affected 4 vs 2 (error +2); R420Q affected 2 vs 1 (error +1); E1724K affected 2 vs 1 (error +1); Q3861H affected 2 vs 1 (error +1); V4771I affected 3 vs 2 (error +1); Q3304E affected 2 vs 1 (error +1); A4741V affected 2 vs 1 (error +1); T1223A unaffected 1 vs 0 (error +1); P1256T unaffected 1 vs 0 (error +1); L1518F unaffected 1 vs 0 (error +1); M2192L unaffected 1 vs 0 (error +1); G2342R unaffected 1 vs 0 (error +1); T2390I unaffected 1 vs 0 (error +1); D3973N unaffected 1 vs 0 (error +1); S4124N unaffected 1 vs 0 (error +1)

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
