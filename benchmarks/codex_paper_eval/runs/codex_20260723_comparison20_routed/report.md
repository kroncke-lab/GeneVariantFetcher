# Codex extraction-blinded paper evaluation — `codex_20260723_comparison20_routed`

## Technical summary

This hash-locked run evaluated **20 papers** (**5 per cardiac gene**) after selecting only PMIDs with downloaded source and gold assertions for carriers, affected, and unaffected. Codex predictions were finalized before scoring.

- Variant precision **84.0%**, recall **67.5%**, F1 **74.9%** (79 TP, 15 FP, 38 FN).
- Exact API telemetry: **217,579 total tokens** (180,415 input; 37,164 output).
- Elapsed: **622.1s wall clock**; 621.9s summed per-paper route + read time.
- Representation choices: {'text': 15, 'table': 3, 'pdf': 1, 'ocr': 1}.

## Blinding and scorer audit

- Paper selection used the fixed manifest `comparison_papers_20260723.tsv` (20 papers) from the downloaded-source, gold-count-eligible pool. Routing, extraction, counts, evidence, and source locations were gold-value-blind.
- `selection.json` contains source metadata and hashes but no gold values or gold row counts. `predictions.json` was made read-only and SHA-256 locked before scoring first opened the gold CSVs.

## Count fidelity

Count recall is the share of all gold count assertions for which the locked prediction supplied a value; MAE/RMSE are computed only where both gold and prediction supplied a value.

| field | supplied / gold assertions | count recall | MAE | RMSE |
|---|---:|---:|---:|---:|
| carriers | 63 / 117 | 53.8% | 0.127 | 0.563 |
| affected | 63 / 117 | 53.8% | 0.143 | 0.577 |
| unaffected | 63 / 117 | 53.8% | 0.016 | 0.126 |

## Per-gene results

| gene | TP | FP | FN | precision | recall | F1 | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| SCN5A | 7 | 4 | 1 | 63.6% | 87.5% | 73.7% | 12.5% / 4.000 / 4.000 | 12.5% / 4.000 / 4.000 | 12.5% / 0.000 / 0.000 |
| KCNH2 | 10 | 9 | 1 | 52.6% | 90.9% | 66.7% | 63.6% / 0.286 / 0.535 | 63.6% / 0.429 / 0.655 | 45.5% / 0.200 / 0.447 |
| KCNQ1 | 6 | 1 | 5 | 85.7% | 54.5% | 66.7% | 18.2% / 0.000 / 0.000 | 18.2% / 0.000 / 0.000 | 9.1% / 0.000 / 0.000 |
| RYR2 | 56 | 1 | 31 | 98.2% | 64.4% | 77.8% | 60.9% / 0.038 / 0.194 | 60.9% / 0.038 / 0.194 | 64.4% / 0.000 / 0.000 |

## Per-paper results

| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|
| SCN5A | 16054936 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 10.1 | 19,598 |
| SCN5A | 16929919 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.1 | 3,433 |
| SCN5A | 8917568 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.0 | 10,926 |
| SCN5A | 17544529 | text | 1 | 1 | 0 | 50.0% | 100.0% | 66.7% | 100.0% / 4.000 | 100.0% / 4.000 | 100.0% / 0.000 | 19.8 | 6,724 |
| SCN5A | 23414114 | table | 1 | 3 | 0 | 25.0% | 100.0% | 40.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 37.7 | 12,524 |
| KCNH2 | 18452873 | table | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 1.000 | 100.0% / 1.000 | 23.1 | 8,987 |
| KCNH2 | 16043162 | text | 2 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 1.000 | 100.0% / 1.000 | 0.0% / n/a | 8.3 | 4,554 |
| KCNH2 | 29121719 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.5 | 11,381 |
| KCNH2 | 30844837 | pdf | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 5.6 | 2,700 |
| KCNH2 | 15851119 | table | 4 | 9 | 0 | 30.8% | 100.0% | 47.1% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 59.2 | 11,696 |
| KCNQ1 | 27000522 | text | 0 | 1 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 14.0 | 3,234 |
| KCNQ1 | 21138517 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 13.5 | 8,579 |
| KCNQ1 | 10220144 | text | 4 | 0 | 3 | 100.0% | 57.1% | 72.7% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 21.4 | 9,063 |
| KCNQ1 | 31245483 | ocr | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 0.0% / n/a | 54.5 | 14,717 |
| KCNQ1 | 18567635 | text | 0 | 0 | 1 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.6 | 31,359 |
| RYR2 | 12093772 | text | 2 | 0 | 12 | 100.0% | 14.3% | 25.0% | 14.3% / 1.000 | 14.3% / 1.000 | 14.3% / 0.000 | 28.6 | 10,771 |
| RYR2 | 14571276 | text | 3 | 0 | 0 | 100.0% | 100.0% | 100.0% | 0.0% / n/a | 0.0% / n/a | 100.0% / 0.000 | 11.9 | 3,972 |
| RYR2 | 29925740 | text | 50 | 1 | 1 | 98.0% | 98.0% | 98.0% | 98.0% / 0.000 | 98.0% / 0.000 | 98.0% / 0.000 | 253.8 | 24,697 |
| RYR2 | 23595086 | text | 0 | 0 | 18 | 0.0% | 0.0% | 0.0% | 0.0% / n/a | 0.0% / n/a | 0.0% / n/a | 7.7 | 10,653 |
| RYR2 | 27225049 | text | 1 | 0 | 0 | 100.0% | 100.0% | 100.0% | 100.0% / 0.000 | 100.0% / 0.000 | 100.0% / 0.000 | 9.5 | 8,011 |

## Errors and representation choices

### SCN5A PMID 16054936

**text** — The running full text most clearly establishes that the human familial variant is in SCN1A, while SCN5A/Gln1478Lys was an engineered homolog used only for functional assays. The preserved tables and PDF webtable describe SCN1A mutation carriers or electrophysiology, not SCN5A carrier, affected, or unaffected person counts.

- Missed gold variants: Q1476K

### SCN5A PMID 16929919

**text** — The text is the only substantive representation and identifies the SCN5A Y1795C and Y1795H variants with their associated phenotypes. The table representation contains only metadata and no variant-level carrier, affected, or unaffected counts.

- No scored variant or count disagreement.

### SCN5A PMID 8917568

**text** — The running full text directly identifies the three SCN5A variants and describes a heterologous functional study. The available images appear to be experimental figures rather than pedigrees needed for carrier or phenotype counts, so OCR is not authoritative for human genotype-phenotype evidence.

- No scored variant or count disagreement.

### SCN5A PMID 17544529

**text** — The running full text is the only substantive representation and reports the SCN5A variants and cohort-level carrier findings. The table representation contains metadata only and no preserved variant-level rows or person counts.

- Extra predictions: D1823D
- Count disagreements: H558R carriers 5 vs 1 (error +4); H558R affected 5 vs 1 (error +4)

### SCN5A PMID 23414114

**table** — The structured table directly reports SCN5A variants with genotype counts, including heterozygous and homozygous carriers and major-allele noncarriers across ancestry groups and overall. This preserves variant-level person counts more clearly than the repetitive running-text preview.

- Extra predictions: SCN5A c.103G>A (p.G35S), SCN5A c.647C>T (p.S216L), SCN5A c.1127G>A (p.R376H)

### KCNH2 PMID 18452873

**table** — The structured patient-level table directly links the KCNH2 insertion to patient 6 and preserves genotype, cohort, QTc, and age at AF onset. This is clearer for variant-level carrier and phenotype evidence than the fragmented running-text preview; no image-based evidence is required.

- Count disagreements: KCNH2 INS "GCGCGGGCG" 569–570 (INS GAG 189–190) affected 1 vs 0 (error +1); KCNH2 INS "GCGCGGGCG" 569–570 (INS GAG 189–190) unaffected 0 vs 1 (error -1)

### KCNH2 PMID 16043162

**text** — The text explicitly identifies KCNH2 variants G873S and N985S in two unrelated affected patients. The table representation contains only artifact metadata and no variant-level carrier or phenotype counts.

- Count disagreements: G873S carriers 1 vs 0 (error +1); N985S carriers 1 vs 0 (error +1); G873S affected 1 vs 0 (error +1); N985S affected 1 vs 0 (error +1)

### KCNH2 PMID 29121719

**text** — The running full text is the most complete representation: it names the KCNH2 variants (R92C, M645R, and K897T) and provides cohort-level carrier/phenotype context. The table preview contains only aggregate long-QT-gene counts and its artifact extraction failed content validation, so it does not preserve KCNH2 variant-level person counts.

- No scored variant or count disagreement.

### KCNH2 PMID 30844837

**pdf** — No usable article body or variant-level table is available. The PDF is the most complete representation of the recovered artifact and preserves its layout, although the preview appears to contain only publisher permission material rather than KCNH2 carrier or phenotype evidence.

- Missed gold variants: D609G

### KCNH2 PMID 15851119

**table** — The preserved table directly links each KCNH2 variant to case or family-member counts, cardiac-event phenotype, event timing, and genotype, making it the clearest source for variant-level carrier and affected-person evidence.

- Extra predictions: T65P, Y475del, G572S, N588D, N633S, V822M, D837G, N996I, R1033fs/23

### KCNQ1 PMID 27000522

**text** — Text is the only available representation. It is an abstract-only fallback and reports one KCNQ1-positive case but lacks variant-level carrier, affected, and unaffected counts; no table, PDF, or OCR evidence is available.

- Missed gold variants: L266P
- Extra predictions: unspecified likely pathogenic KCNQ1 variant

### KCNQ1 PMID 21138517

**text** — The running full text is available from PMC XML and is the clearest source for the KCNQ1 carrier cohort and phenotype information. No structured tables are present, and the figures show aggregate QT responses rather than variant-level carrier, affected, or unaffected counts.

- No scored variant or count disagreement.

### KCNQ1 PMID 10220144

**text** — The running full text contains the mutation and family-segregation narrative, whereas the table representation contains only artifact metadata and no preserved variant-level rows.

- Missed gold variants: G314S, Y315S, G345R

### KCNQ1 PMID 31245483

**ocr** — Figure 1B contains the family pedigree, and the caption states that members with both M437V and QT prolongation are marked visually in red. The pedigree image is therefore necessary to recover variant-level carrier and phenotype counts that are not enumerated in the text or tables.

- No scored variant or count disagreement.

### KCNQ1 PMID 18567635

**text** — The full running text is the clearest authoritative source for the KCNQ1 variants and functional findings. No structured tables or person-level carrier, affected, or unaffected counts are indicated, and OCR figures are not needed to recover genotype-phenotype evidence.

- Missed gold variants: F340W

### RYR2 PMID 12093772

**text** — The running full text is the most complete available representation and includes the clinical, genetic, and family-member descriptions. The table artifact preview contains only metadata and captions, without structured variant-level carrier, affected, or unaffected counts.

- Missed gold variants: R4497C, R2474S, N4104K, V4771I, E2311D, N4895D, L3778F, I4867M, G3946S, E4950K, S2246L, A4860G
- Count disagreements: S2246L carriers 2 vs 1 (error +1); G3946S carriers 2 vs 1 (error +1); S2246L affected 2 vs 1 (error +1); G3946S affected 2 vs 1 (error +1)

### RYR2 PMID 14571276

**text** — The text contains the only substantive variant-level evidence, naming RYR2 V2306I, P4902L, and R4959Q and stating that they were absent in unaffected and control individuals. The table representation contains only artifact metadata and no preserved variant or person-count rows.

- No scored variant or count disagreement.

### RYR2 PMID 29925740

**text** — The running text contains the relevant RYR2 carrier information, including that 9 of 104 initially diagnosed LQTS patients carried RYR2 mutations. The table representation contains only extraction metadata and no variant-level person-count rows.

- Missed gold variants: N4168S
- Extra predictions: c.12533A>G (p.N4178S)

### RYR2 PMID 23595086

**text** — The table representation contains only extraction metadata and no structured variant or person-level rows. The running main text is therefore the most complete available source for RYR2 variant, carrier, affected, and unaffected evidence.

- Missed gold variants: R407S, N1551S, M2192L, A2387V, G2400T, R2474G, G2628E, D3638A, Q3861H, D3876E, M4002I, K4392R, S4124R, I4587V, Y4725C, K4750Q, L4865V, L4919S

### RYR2 PMID 27225049

**text** — The running full text contains the case-level phenotype, the RYR2 c.7169C>T (p.T2390I) result, and familial inheritance evidence involving the father. No table, PDF, or OCR representation is available.

- No scored variant or count disagreement.

## Scope, method, and limitations

- Population: fixed manifest `comparison_papers_20260723.tsv` (20 papers); 5 per cardiac gene; every PMID has downloaded source and at least one gold assertion in each count field.
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
