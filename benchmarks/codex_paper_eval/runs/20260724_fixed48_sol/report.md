# Codex extraction-blinded paper evaluation — `20260724_fixed48_sol`

## Technical summary

This hash-locked run evaluated **48 papers** (**12 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **87.7%**, recall **72.3%**, F1 **79.3%** (724 TP, 102 FP, 277 FN).
- Exact API telemetry: **881,425 total tokens** (644,348 input; 237,077 output).
- Elapsed: **4554.9s wall clock**; 4553.9s summed per-paper route + read time.
- Representation choices: {'text': 19, 'table': 20, 'ocr': 6, 'pdf': 3}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `highcarrier_48_papers_20260723.tsv` (48 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- `selection.json` contains source metadata and hashes but no gold values or gold row counts. `predictions.json` was made read-only and SHA-256 locked before scoring first opened the gold CSVs.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 642 / 1001 | 64.1% | 0.481 | 4.881 |
| affected | 329 / 1001 | 32.9% | 0.687 | 6.008 |
| unaffected | 169 / 1001 | 16.9% | 0.852 | 3.784 |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 179 | 16 | 41 | 91.8% | 81.4% | 86.3% | 79.5% / 0.526 / 4.666 | 34.5% / 0.921 / 5.166 | 13.6% / 1.367 / 5.730 |
| KCNH2 | 206 | 39 | 86 | 84.1% | 70.5% | 76.7% | 50.0% / 1.007 / 8.319 | 25.3% / 0.351 / 1.755 | 4.5% / 2.923 / 5.262 |
| KCNQ1 | 173 | 35 | 98 | 83.2% | 63.8% | 72.2% | 63.5% / 0.209 / 2.745 | 11.4% / 3.065 / 17.063 | 0.4% / 26.000 / 26.000 |
| RYR2 | 166 | 12 | 52 | 93.3% | 76.1% | 83.8% | 68.3% / 0.228 / 0.751 | 67.9% / 0.236 / 2.005 | 57.3% / 0.312 / 1.787 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 27566755 | text | 51 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 11.8% / 0.000 | 11.8% / 0.000 | 178.6 | 21,452 |
| SCN5A | 26669661 | table | 27 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 60.3 | 68,733 |
| SCN5A | 20470418 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 59.000 | 100.0% / 22.000 | 0.0% / n/a | 16.5 | 9,193 |
| SCN5A | 28339995 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 77.4 | 19,945 |
| SCN5A | 29709101 | table | 10 | 2 | 1 | 83.3% | 90.9% | 87.0% | 90.9% / 0.000 | 90.9% / 0.000 | 27.3% / 0.000 | 75.3 | 19,082 |
| SCN5A | 28341781 | text | 53 | 2 | 2 | 96.4% | 96.4% | 96.4% | 96.4% / 0.000 | 96.4% / 0.000 | 29.1% / 0.000 | 401.6 | 25,552 |
| SCN5A | 18451998 | ocr | 2 | 4 | 0 | 33.3% | 100.0% | 50.0% | 50.0% / 9.000 | 50.0% / 38.000 | 50.0% / 29.000 | 63.9 | 16,371 |
| SCN5A | 26921764 | table | 27 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.074 | 0.0% / n/a | 0.0% / n/a | 111.0 | 15,432 |
| SCN5A | 27554632 | table | 3 | 0 | 6 | 100.0% | 33.3% | 50.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 512.9 | 12,267 |
| SCN5A | 10590249 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 10.000 | 100.0% / 10.000 | 100.0% / 0.000 | 38.4 | 15,665 |
| SCN5A | 25051102 | table | 1 | 8 | 2 | 11.1% | 33.3% | 16.7% | 33.3% / 0.000 | 33.3% / 0.000 | 0.0% / n/a | 54.3 | 18,959 |
| SCN5A | 26746457 | text | 2 | 0 | 30 | 100.0% | 6.2% | 11.8% | 6.2% / 6.000 | 6.2% / 0.000 | 6.2% / 6.000 | 20.1 | 10,509 |
| KCNH2 | 29622001 | table | 34 | 1 | 1 | 97.1% | 97.1% | 97.1% | 97.1% / 0.735 | 22.9% / 0.000 | 0.0% / n/a | 261.5 | 26,211 |
| KCNH2 | 11854117 | table | 44 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 103.5 | 15,593 |
| KCNH2 | 14661677 | ocr | 0 | 0 | 29 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.4 | 3,267 |
| KCNH2 | 19160088 | text | 2 | 0 | 1 | 100.0% | 66.7% | 80.0% | 66.7% / 48.000 | 0.0% / n/a | 0.0% / n/a | 73.1 | 14,590 |
| KCNH2 | 26496715 | table | 53 | 1 | 1 | 98.1% | 98.1% | 98.1% | 98.1% / 0.000 | 98.1% / 0.000 | 0.0% / n/a | 163.2 | 24,847 |
| KCNH2 | 11844290 | text | 0 | 0 | 5 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 10.1 | 4,894 |
| KCNH2 | 10973849 | table | 60 | 34 | 0 | 63.8% | 100.0% | 77.9% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 377.7 | 32,709 |
| KCNH2 | 10862094 | text | 6 | 3 | 2 | 66.7% | 75.0% | 70.6% | 75.0% / 3.833 | 75.0% / 1.333 | 75.0% / 3.500 | 157.9 | 20,745 |
| KCNH2 | 10841244 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 2.000 | 100.0% / 0.000 | 100.0% / 0.000 | 13.0 | 11,491 |
| KCNH2 | 23864605 | table | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 8.500 | 100.0% / 8.500 | 11.8 | 10,691 |
| KCNH2 | 24667783 | ocr | 4 | 0 | 19 | 100.0% | 17.4% | 29.6% | 17.4% / 0.250 | 17.4% / 0.250 | 17.4% / 0.000 | 53.1 | 9,767 |
| KCNH2 | 19038855 | text | 0 | 0 | 28 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 61.2 | 17,337 |
| KCNQ1 | 19490272 | table | 54 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 145.2 | 14,420 |
| KCNQ1 | 23153844 | table | 21 | 13 | 0 | 61.8% | 100.0% | 76.4% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 47.8 | 13,557 |
| KCNQ1 | 17470695 | table | 56 | 21 | 0 | 72.7% | 100.0% | 84.2% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 178.5 | 19,340 |
| KCNQ1 | 14678125 | text | 0 | 0 | 41 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 8.1 | 7,371 |
| KCNQ1 | 28720088 | pdf | 1 | 1 | 1 | 50.0% | 50.0% | 50.0% | 50.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 21.3 | 18,934 |
| KCNQ1 | 21129503 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 22.3 | 11,682 |
| KCNQ1 | 25087618 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 95.000 | 100.0% / 26.000 | 40.6 | 16,332 |
| KCNQ1 | 17192539 | text | 1 | 0 | 56 | 100.0% | 1.8% | 3.4% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 26.4 | 8,333 |
| KCNQ1 | 24052033 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 14.1 | 13,384 |
| KCNQ1 | 18713323 | table | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 0.0% / n/a | 0.0% / n/a | 24.9 | 11,906 |
| KCNQ1 | 29197658 | table | 29 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 64.6 | 14,957 |
| KCNQ1 | 33141630 | pdf | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 36.000 | 100.0% / 0.000 | 0.0% / n/a | 42.9 | 41,531 |
| RYR2 | 25814417 | table | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 6.000 | 100.0% / 24.000 | 100.0% / 19.000 | 100.5 | 17,536 |
| RYR2 | 29925740 | text | 50 | 1 | 1 | 98.0% | 98.0% | 98.0% | 98.0% / 0.000 | 98.0% / 0.000 | 98.0% / 0.000 | 152.2 | 22,601 |
| RYR2 | 33315912 | table | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 3.000 | 0.0% / n/a | 42.4 | 10,737 |
| RYR2 | 16272262 | ocr | 12 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 96.6 | 16,321 |
| RYR2 | 34202968 | ocr | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 108.2 | 18,620 |
| RYR2 | 19926015 | text | 1 | 2 | 39 | 33.3% | 2.5% | 4.7% | 2.5% / 1.000 | 2.5% / 0.000 | 2.5% / 1.000 | 110.5 | 32,861 |
| RYR2 | 33686871 | ocr | 1 | 0 | 7 | 100.0% | 12.5% | 22.2% | 12.5% / 0.000 | 12.5% / 0.000 | 12.5% / 0.000 | 42.7 | 17,454 |
| RYR2 | 28237968 | pdf | 18 | 7 | 0 | 72.0% | 100.0% | 83.7% | 100.0% / 0.167 | 100.0% / 0.000 | 100.0% / 0.167 | 110.9 | 27,555 |
| RYR2 | 12106942 | text | 6 | 0 | 0 | 100.0% | 100.0% | 100.0% | 33.3% / 2.000 | 16.7% / 0.000 | 16.7% / 4.000 | 29.6 | 17,031 |
| RYR2 | 22677073 | table | 22 | 1 | 4 | 95.7% | 84.6% | 89.8% | 84.6% / 0.045 | 84.6% / 0.045 | 0.0% / n/a | 54.2 | 15,662 |
| RYR2 | 30403697 | table | 20 | 1 | 1 | 95.2% | 95.2% | 95.2% | 95.2% / 0.200 | 95.2% / 0.000 | 95.2% / 0.200 | 84.8 | 31,332 |
| RYR2 | 33606749 | table | 33 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.455 | 100.0% / 0.212 | 100.0% / 0.242 | 80.6 | 16,666 |

## Errors and representation choices

### SCN5A PMID 27566755

**text** — The PMC full text is substantially complete and includes the cohort composition, SCN5A mutation context, and four table captions, whereas the available OCR images appear limited to figures and provide no indication of pedigrees or necessary genotype/phenotype carrier evidence.

- No scored variant or count disagreement.

### SCN5A PMID 26669661

**table** — The supplemental tables are the only representation preserving structured variant-level counts, including total families, total carriers, and maternal/paternal transmissions. This is more suitable for SCN5A variant-level curation than the aggregate abstract or unlabeled figure images.

- No scored variant or count disagreement.

### SCN5A PMID 20470418

**text** — The PMC full text explicitly reports variant-level carrier and phenotype-group evidence, including 26 heterozygous p.Ser1103Tyr infant hearts and comparisons between SIDS and other-cause deaths. The lone OCR image is not indicated to contain necessary pedigree or genotype/phenotype counts.

- Count disagreements: SCN5A c.3308C>A (p.Ser1103Tyr; also designated p.Ser1102Tyr) carriers 26 vs 85 (error -59); SCN5A c.3308C>A (p.Ser1103Tyr; also designated p.Ser1102Tyr) affected 17 vs 39 (error -22)

### SCN5A PMID 28339995

**text** — Text is the only available authoritative representation and contains publisher-sourced full article content with explicit SCN5A D1790G carrier and cardiac-event counts. No table, PDF, or OCR representation is available to preserve additional variant-level evidence.

- No scored variant or count disagreement.

### SCN5A PMID 29709101

**table** — The structured variant table directly links each SCN5A DNA/protein change to patient counts stratified by BrS-positive and BrS-negative status, best preserving variant-level carrier and phenotype evidence. The figures contain ECG examples rather than necessary genotype–phenotype counts, while the text preview is truncated.

- Missed gold variants: c.4813+3_4813_6dup
- Extra predictions: c.4813+3_4813+6dup, c.4978A>G (p.Ile1660Val)

### SCN5A PMID 28341781

**text** — Text is the only available representation and contains the running original research article, including cohort methods and results relevant to SCN5A probands. No table, PDF, or OCR representation is supplied to preserve additional variant-level evidence.

- Missed gold variants: c.3840+1G>A, c.4245+1G>A
- Extra predictions: IVS21+1 G>A, IVS23+1 G>A

### SCN5A PMID 18451998

**ocr** — Figure 1 contains pedigrees for 15 E1784K families and visually encodes carrier status, symptomatic versus asymptomatic carriers, QT prolongation, sudden cardiac death, ST elevation, and sinus node dysfunction. OCR of the pedigree images is therefore necessary to preserve person-level genotype/phenotype evidence that the truncated text and captions summarize but do not fully enumerate.

- Extra predictions: V1098L, 1795insD, ΔKPQ, ΔK1500
- Count disagreements: E1784K carriers 41 vs 50 (error -9); E1784K affected 9 vs 47 (error -38); E1784K unaffected 32 vs 3 (error +29)

### SCN5A PMID 26921764

**table** — The structured table directly links each SCN5A nucleotide and coding variant to the number of patients carrying it, preserving variant-level carrier counts more clearly than the running text, which mainly reports aggregate genotype and outcome totals.

- Count disagreements: 5228 G>A (G1743E) carriers 4 vs 2 (error +2)

### SCN5A PMID 27554632

**table** — The article contains a structured variant table reporting SCN5A variants across DCM, HCM, and control groups, making it the clearest source for variant-level affected and unaffected carrier counts.

- Missed gold variants: M1245I, R513C, R513H, R526H, T512I, T570N

### SCN5A PMID 10590249

**text** — Running full text is the only available representation and includes the clinical and genetic methods, mutation-carrier definitions, and pedigree figure captions relevant to variant-level evidence. No table, PDF, or OCR representation is available.

- Count disagreements: 1795insD (TGA insertion at nucleotide 5537) carriers 53 vs 43 (error +10); 1795insD (TGA insertion at nucleotide 5537) affected 53 vs 43 (error +10)

### SCN5A PMID 25051102

**table** — Structured article tables are the clearest representation for SCN5A variant-level genotype and cohort counts, preserving comparisons among AMI/VF-positive patients, AMI/VF-negative patients, and controls. The PDF is an unrelated library form, and no pedigree or figure evidence requires OCR.

- Missed gold variants: c.2788-6C>T, R1193Q
- Extra predictions: 87G>A, 630G>A, 1673A>G (His558Arg), IVS16-6C>T, 3183G>A, 3578G>A, 4509C>T, 5457T>A

### SCN5A PMID 26746457

**text** — The PMC full text is available and contains the cohort-level carrier and phenotype counts, whereas the sole OCR image is not indicated to contain necessary genotype–phenotype or pedigree evidence. No usable structured table or PDF representation is available.

- Missed gold variants: A1180V, A1270S, A1680T, A551V, C982R, D1275N, E1053K, E48K, F1596I, F532C, G1935S, G552W, G615E, I1836T, I848F, L1194M, L1704H, P717L, R1023H, R1195S, R1512W, R1739Q, R1898C, R2012C, R282C, R689H, S1904L, T1304M, V1353M, V1532I
- Count disagreements: SCN5A p.R1193Q carriers 19 vs 7 (error +12); SCN5A p.R1193Q unaffected 12 vs 0 (error +12)

### KCNH2 PMID 29622001

**table** — Structured tables provide KCNH2 mutation-specific carrier counts and family counts, including p.R176W and p.L552S, more clearly than the running text. The OCR inventory appears limited to figures and is not necessary for the genotype/person-count evidence.

- Missed gold variants: P1034fsX
- Extra predictions: c.3093_3106del
- Count disagreements: p.L552S carriers 73 vs 74 (error -1); c.453delC carriers 24 vs 0 (error +24)

### KCNH2 PMID 11854117

**table** — The structured table most clearly preserves KCNH2 variant-level subject counts, along with mutation location, exon, type, and family footnotes. The running-text preview is fragmented and repeatedly embeds the same table, while no PDF or OCR content is available.

- No scored variant or count disagreement.

### KCNH2 PMID 14661677

**ocr** — The running text has no usable article body, while the table and PDF previews are unrelated administrative or journal-metric material. The OCR inventory is the only representation that may preserve figure or pedigree-based genotype/phenotype evidence omitted from the textual previews.

- Missed gold variants: A1058E, A190T, A203T, A915V, C723R, G187del, G187S, G873S, H254Q, K897T, K897T, K897T, K897T, L1023del, N257H, N33T, P251A, P347S, P910L, P917L, P967L, Q1068R, R1035W, R1047L, R1047L, R176W, R181Q, T367S, V215G

### KCNH2 PMID 19160088

**text** — The PMC full text directly reports variant-level evidence for KCNH2 R176W, including 16 carriers and the associated QT-interval effect. The available OCR is limited to one image and is not needed to recover genotype or phenotype counts.

- Missed gold variants: R176W
- Count disagreements: KCNH2 R176W carriers 16 vs 112 (error -96)

### KCNH2 PMID 26496715

**table** — The recovered content is supplement-only, and the structured mutation tables preserve variant-level patient counts more clearly than the duplicated and truncated text representation.

- Missed gold variants: A797T
- Extra predictions: c.2389G>A (p.Arg797Thr)

### KCNH2 PMID 11844290

**text** — Text is the only available representation and contains the article's running full text; no table, PDF, or OCR representation is supplied to better preserve variant-level carrier and phenotype evidence.

- Missed gold variants: R752W, F805C, M124R, V822M, W1001X

### KCNH2 PMID 10973849

**table** — The structured tables include HERG (KCNH2) mutation-level fields such as nucleotide change, coding effect, and number of individuals, making them the clearest source for variant-level carrier counts. The figures are topology schematics and do not preserve person-level phenotype evidence.

- Extra predictions: dup558–600 (L200fs/144), del1261 (Y420fs/12), C1283A (S428X), C1307T (T436M), A1408G (N470D), C1421T (T474I), C1479G (Y493X), del1498–1524 (del500–508), C1600T (R534C), delT1671 (T556fs/7), G1672C (A558P), G1714C (G572R), G1714T (G572C), C1744T (R582C), A1762G (N588D), T1778G (I593G), G1801A (G601S), T1831C (Y611H), T1833(A or G) (Y611X), G1834T (V612L), C1841T (A614V), A1885G (N629D), A1886G (N629S), C1887A (N629K), G1888C (V630L), T1889C (V630A), A1898G (N633S), C1920A (F640L), del1951–1952 (L650fs/2), C2173T (Q725X), dup2356–2386 (V796fs/22), C2453T (S818L), G2464A (V822M), insG3107–3108 (G1036fs/82)

### KCNH2 PMID 10862094

**text** — The reharvested running full text is the only available representation and includes the KCNH2/HERG variant survey and associated cohort evidence; no table, PDF-layout, or OCR representation is available.

- Missed gold variants: P151fsX, S543fsX
- Extra predictions: c.453delC (453delC), c.1631_1632delAG (1631delAG), c.2690A>C (p.Lys897Thr; K897T)
- Count disagreements: c.1801G>A (p.Gly601Ser; G601S) carriers 19 vs 3 (error +16); c.1750G>A (p.Gly584Ser; G584S) carriers 3 vs 10 (error -7); c.1801G>A (p.Gly601Ser; G601S) affected 4 vs 1 (error +3); c.526C>T (p.Arg176Trp; R176W) affected 1 vs 4 (error -3); c.1750G>A (p.Gly584Ser; G584S) affected 2 vs 4 (error -2); c.1801G>A (p.Gly601Ser; G601S) unaffected 15 vs 2 (error +13); c.526C>T (p.Arg176Trp; R176W) unaffected 5 vs 2 (error +3); c.1750G>A (p.Gly584Ser; G584S) unaffected 1 vs 6 (error -5)

### KCNH2 PMID 10841244

**text** — Text is the only available representation and includes variant-level genotype and phenotype evidence for HERG/KCNH2 L552S, including two affected homozygous siblings, heterozygous parents, 38 additional carriers, and symptomatic, asymptomatic, and noncarrier groups.

- Count disagreements: c.1655T>C (p.Leu552Ser; L552S, HERG-Fin) carriers 42 vs 44 (error -2)

### KCNH2 PMID 23864605

**table** — The structured table directly preserves variant-level genotype-positive subject counts and symptomatic counts for KCNH2 G601S and A614V. The figures are cell-biological experiments rather than pedigrees, and the running-text preview is less direct for carrier and phenotype counts.

- Count disagreements: G601S affected 3 vs 0 (error +3); A614V affected 14 vs 0 (error +14); G601S unaffected 6 vs 0 (error +6); A614V unaffected 11 vs 0 (error +11)

### KCNH2 PMID 24667783

**ocr** — The textual preview omits detailed segregation evidence, and the table preview primarily shows index-patient variant fields and aggregate phenotype summaries rather than carrier-level affected and unaffected relatives. The figure images are therefore most likely to preserve the pedigrees needed for variant-level genotype–phenotype counts.

- Missed gold variants: G262fsX, G785fsX, L559H, M645V, N629S, P596L, P72R, R176W, R534C, R582C, S428X, S818L, A172V, E807X, G238R, L343fsX, S855R, S890C, W705fsX
- Count disagreements: p.Arg356His carriers 2 vs 1 (error +1); p.Arg356His affected 2 vs 1 (error +1)

### KCNH2 PMID 19038855

**text** — The PMC running full text is complete and includes the study’s cohort/genotype results and two table-associated sections. The available images show aggregate LQTS-subtype seizure prevalence rather than necessary variant-level carrier, affected, or unaffected evidence, so OCR is not authoritative.

- Missed gold variants: A193fsX, A558P, C64Y, D456Y, D501H, E698X, E876X, F640L, G1036fsX, G306W, G572S, G925fsX, I30T, I560fsX, I642del, M645L, N629I, P241L, P872fsX, Q676fsX, R176W, R252fsX, R366X, R534C, R582C, R73fsX, T613M, Y99S

### KCNQ1 PMID 19490272

**table** — The structured mutation table provides variant-level carrier counts for individual KCNQ1 missense variants, while the accompanying table preserves affected-event counts by conservation tertile. This is clearer and less ambiguous than the running text for count extraction.

- No scored variant or count disagreement.

### KCNQ1 PMID 23153844

**table** — The structured mutation table provides variant-level KCNQ1 carrier counts for all 34 mutations, whereas the running text primarily gives aggregate cohort and outcome data without preserving per-variant person counts.

- Extra predictions: M1V, M159sp, Y171X, L191fs/90, S225L, R243C, W305S, T312I, G314S, A341E, S349W, P400fs/62, Q530X

### KCNQ1 PMID 17470695

**table** — The structured table directly maps individual KCNQ1 variants to subject counts, mutation types, and functional effects, providing the clearest available variant-level carrier evidence. The figures only show approximate frequency categories and do not preserve affected/unaffected counts.

- Extra predictions: G57V, W120C, G189R, G269D, Y278H, E284K, V310I, D317G, A344A/sp [1032 G>A], A344V, L353P, Q357H, R360G, K393N, R397W, P448fs/13 [ins G 1344-1345], I567S, R591C, IVS2+1 G>A, IVS4+5 G>A, IVS7+5 G>A

### KCNQ1 PMID 14678125

**text** — The running text is the only available representation and provides the KCNQ1 carrier cohort size, mutation-region subgroup counts, and cardiac-event phenotype definitions and comparisons. No table, PDF, or OCR content is available to preserve additional variant-level evidence.

- Missed gold variants: T144A, L151, G168R, Y171X, V172M, A178P, G185S, R190Q, S225L, R243C, V254L, V254M, L266P, G269D, G269S, L273F, E284K, A300T, G306R, T311I, T312I, Y315C, D317G, DELF340, A341E, A341V, A344/SP, A344A/SPLICE, G345R, S349W, Q357H, R366Q, R366W, K393N, R518X, V524G, Q530X, R539W, S566F, I567S, R594Q

### KCNQ1 PMID 28720088

**pdf** — The article contains four tables, but no structured table representation is available. The PDF is most likely to preserve the table layout and variant-level counts for KCNQ1 Y111C and R518* carriers and their phenotype groupings that may be degraded in running text.

- Missed gold variants: R519X
- Extra predictions: KCNQ1 R518*

### KCNQ1 PMID 21129503

**text** — Text is the only available representation and directly reports variant-level evidence for Y111C/KCNQ1, including 170 carriers across 37 proband families. No table, PDF, or OCR content is available to provide more complete affected/unaffected counts.

- No scored variant or count disagreement.

### KCNQ1 PMID 25087618

**text** — The PMC full text directly reports KCNQ1 p.Ala341Val carrier and non-carrier counts and describes phenotype and cardiac-event analyses. The OCR inventory contains analysis-flow and modeled-result figures rather than necessary genotype/phenotype evidence.

- Count disagreements: p.Ala341Val (A341V) affected 95 vs 0 (error +95); p.Ala341Val (A341V) unaffected 26 vs 0 (error +26)

### KCNQ1 PMID 17192539

**text** — Text is the only available representation and contains the full article narrative, including genotyped family cohorts, carrier counts, and transmission/phenotype evidence relevant to KCNQ1.

- Missed gold variants: P117L, G168R, R174C, R174H, L175fsX, A178T, Y184S, R190W, R190Q, R190L, S225L, R231C, R243C, V254M, H258N, H258P, E261K, G269D, S277W, V280A, V280E, V288fsX, G306R, T309R, T311I, T312I, G314S, Y315S, Y315C, G316E, P320H, P320A, G325R, F339S, A341V, L342F, A344V, SP/A344A, Q359-K362DEL, R366W, A371T, S373P, W379S, K422fsX, M476L, M520R, R539W, R555C, R555H, S566P, I567T, T587M, A590T, R591H, R594C, Q604X

### KCNQ1 PMID 24052033

**text** — Running full text is the only available representation and includes the study methods, cohort description, and KCNQ1 variant groups Y111C and R518X. No separate table, PDF, or OCR representation is available.

- No scored variant or count disagreement.

### KCNQ1 PMID 18713323

**table** — Structured table rows provide exact person counts for each of the six KCNQ1 variants, stratified by ethnicity. The figures contain mutation-specific Kaplan–Meier curves but do not preserve exact variant-level affected and unaffected counts as clearly as the tables.

- No scored variant or count disagreement.

### KCNQ1 PMID 29197658

**table** — The structured table directly links each KCNQ1 variant to its reported number of cases and gnomAD count, preserving variant-level person counts more clearly than the narrative text. The figures do not provide necessary pedigree or genotype–phenotype evidence.

- No scored variant or count disagreement.

### KCNQ1 PMID 33141630

**pdf** — The supplemental PDF contains the detailed KCNQ1 p.T224M carrier cohort and phenotype information, including the 124 identified carriers and supplemental tables that are not adequately preserved in the malformed table preview or abbreviated text preview.

- Count disagreements: KCNQ1 c.671C>T (p.Thr224Met; p.T224M; rs199472706; hg38 chr11:2571391C>T) carriers 124 vs 88 (error +36)

### RYR2 PMID 25814417

**table** — Structured patient-level rows for p.G357S RYR2-positive individuals preserve carrier status and phenotype evidence, including symptoms and ventricular arrhythmia findings, more clearly than the running text.

- Count disagreements: RYR2 p.G357S carriers 185 vs 179 (error +6); RYR2 p.G357S affected 97 vs 73 (error +24); RYR2 p.G357S unaffected 87 vs 106 (error -19)

### RYR2 PMID 29925740

**text** — Text is the only available representation. It provides cohort-level RYR2 carrier information, while no table, PDF, or OCR content is available to preserve additional variant-level counts.

- Missed gold variants: N4168S
- Extra predictions: c.12533A>G (p.N4178S)

### RYR2 PMID 33315912

**table** — The structured table most clearly preserves variant-carrier phenotype counts, explicitly separating 41 carriers without cardiac events from 18 with cardiac events and providing associated clinical findings.

- Count disagreements: RYR2 p.Pro2328Ser (P2328S) affected 18 vs 15 (error +3)

### RYR2 PMID 16272262

**ocr** — The pedigree figure is necessary to preserve family-level genotype and phenotype relationships, including variant-specific affected and unaffected carriers, which are omitted or mangled in the text and table representations.

- No scored variant or count disagreement.

### RYR2 PMID 34202968

**ocr** — The study reports segregation across three families, and variant-level carrier, affected, and unaffected evidence is likely encoded in pedigree figures. OCR best preserves the individual genotype and phenotype relationships that running text may only summarize.

- No scored variant or count disagreement.

### RYR2 PMID 19926015

**text** — The PMC running full text is the most complete representation of the study cohort and RYR2 mutation findings. The table and PDF previews primarily contain cross-species conservation data rather than variant-level carrier or phenotype counts, and the figures do not appear necessary for genotype–phenotype evidence.

- Missed gold variants: L62F, M81L, E243K, F329L, R332W, G357S, V377M, T415R, R420Q, V507I, A549V, R739H, R1013Q, T1107M, A1136V, E1837K, E2045G, Y2156C, H2168Q, E2183V, D2216V, EXON 3 DELETION, E2296Q, R2420W, M3972I, D3973H, S4124G, R4157Q, Q4159P, N4178S, E4187Q, G4315E, K4650E, N4736 DEL, R4790Q, K4805R, R4822H, L3974Q, K3997E
- Extra predictions: 3.6-kb deletion involving exon 3, 1.1-kb deletion involving exon 3
- Count disagreements: Y4149S carriers 2 vs 1 (error +1); Y4149S unaffected 1 vs 0 (error +1)

### RYR2 PMID 33686871

**ocr** — The study's variant-level segregation evidence is family-based, and the pedigree figure is needed to preserve individual genotypes and affected versus unaffected status. No structured tables are available, while running text only summarizes four affected homozygotes and asymptomatic heterozygous relatives.

- Missed gold variants: I3995V, D4112N, T4196I, D4646A, Q4879H, K4594R, I2075T

### RYR2 PMID 28237968

**pdf** — The PDF preserves the complete article layout, including likely variant-level tables and person counts that are absent from the HTML text extraction; OCR figures mainly show mutation locations rather than carrier and phenotype counts.

- Extra predictions: NM_001035.2:c.5170G>A (p.Glu1724Lys), NM_001035.2:c.12470G>A (p.Arg4157Gln), NM_001035.2:c.14553C>A (p.Phe4851Leu), NM_001035.2:c.1258C>T (p.Arg420Trp), NM_001035.2:c.3407C>T (p.Ala1136Val), NM_001035.2:c.5656G>A (p.Gly1886Ser), NM_001035.2:c.5654G>A (p.Gly1885Glu)
- Count disagreements: NM_001035.2:c.7160C>T (p.Ala2387Val) carriers 1 vs 4 (error -3); NM_001035.2:c.7160C>T (p.Ala2387Val) unaffected 0 vs 3 (error -3)

### RYR2 PMID 12106942

**text** — Text is the only available representation and contains substantial Elsevier full text, including aggregate carrier and phenotype counts across the RYR2-mutated families. No table, PDF, or OCR representation is available to preserve additional variant-level evidence.

- Count disagreements: R420W carriers 4 vs 8 (error -4); R420W unaffected 2 vs 6 (error -4)

### RYR2 PMID 22677073

**table** — The structured case-level table links each variant and gene, including RYR2, to an individual’s demographics, sudden-death event, sentinel phenotype, and family history. It therefore best preserves variant-level carrier and phenotype evidence; the running text is truncated, and OCR offers only an uncaptioned image inventory.

- Missed gold variants: V2113M, Q2958R, A1136V, R4037C
- Extra predictions: 6739C>T (S2246L)
- Count disagreements: 6737C>T (S2246L) carriers 1 vs 2 (error -1); 6737C>T (S2246L) affected 1 vs 2 (error -1)

### RYR2 PMID 30403697

**table** — The structured subject-level rows directly link each RYR2 variant pair to carrier status, symptoms, asymptomatic relatives, inheritance, phase, and family relationships, best preserving variant-level affected and unaffected evidence.

- Missed gold variants: G4722S
- Extra predictions: p.G4772S
- Count disagreements: p.R2028H carriers 2 vs 3 (error -1); p.Y4721C carriers 2 vs 3 (error -1); c.3599-9delT carriers 1 vs 2 (error -1); c.14091-11dupT carriers 1 vs 2 (error -1); p.R2028H unaffected 1 vs 2 (error -1); p.Y4721C unaffected 1 vs 2 (error -1); c.3599-9delT unaffected 0 vs 1 (error -1); c.14091-11dupT unaffected 0 vs 1 (error -1)

### RYR2 PMID 33606749

**table** — The structured corrected table preserves patient-level RYR2 variants, genotyped relatives, inheritance, and parental phenotypes, providing the clearest variant-level carrier and affected/unaffected evidence.

- Count disagreements: exon 3 deletion (N57_G91del35) carriers 4 vs 2 (error +2); 1259g>a (R420Q) carriers 2 vs 1 (error +1); 3667a>g (T1223A) carriers 2 vs 1 (error +1); 3766c>a (P1256T) carriers 2 vs 1 (error +1); 4552c>t (L1518F) carriers 2 vs 1 (error +1); 5170g>a (E1724K) carriers 2 vs 1 (error +1); 6574a>t (M2192L) carriers 2 vs 1 (error +1); 7024g>a (G2342R) carriers 2 vs 1 (error +1); 7169c>t (T2390I) carriers 2 vs 1 (error +1); 11917g>a (D3973N) carriers 2 vs 1 (error +1); 12371g>a (S4124N) carriers 2 vs 1 (error +1); 14311g>a (V4771I) carriers 3 vs 2 (error +1); 9910c>g (Q3304E) carriers 2 vs 1 (error +1); 14222c>t (A4741V) carriers 2 vs 1 (error +1); exon 3 deletion (N57_G91del35) affected 4 vs 2 (error +2); 1259g>a (R420Q) affected 2 vs 1 (error +1); 5170g>a (E1724K) affected 2 vs 1 (error +1); 14311g>a (V4771I) affected 3 vs 2 (error +1); 9910c>g (Q3304E) affected 2 vs 1 (error +1); 14222c>t (A4741V) affected 2 vs 1 (error +1); 3667a>g (T1223A) unaffected 1 vs 0 (error +1); 3766c>a (P1256T) unaffected 1 vs 0 (error +1); 4552c>t (L1518F) unaffected 1 vs 0 (error +1); 6574a>t (M2192L) unaffected 1 vs 0 (error +1); 7024g>a (G2342R) unaffected 1 vs 0 (error +1); 7169c>t (T2390I) unaffected 1 vs 0 (error +1); 11917g>a (D3973N) unaffected 1 vs 0 (error +1); 12371g>a (S4124N) unaffected 1 vs 0 (error +1)

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
