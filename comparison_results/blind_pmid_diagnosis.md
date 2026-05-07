# Blind PMID Diagnostic — 46 PMIDs / 412 Missing Variants

**Generated:** 2026-05-07
**Source run:** `results/KCNH2/20260503_161203/` (DB used by `comparison_results/summary.json`)
**Comparison input:** `comparison_results/missing_in_sqlite.csv` (816 missing rows across 164 PMIDs)
**Scope:** PMIDs that were harvested into the `papers` table but produced **zero** rows in `variant_papers`, AND appear in the missing-in-sqlite list with at least one gold-standard variant.

## Headline

**The blind-PMID problem is overwhelmingly a harvest-stage problem, not an extraction-stage problem.** Of the 46 blind PMIDs, **40 (87%) have no full-text harvested at all** — only an abstract — and the abstracts genuinely do not contain the variant notations the gold standard expects (the regex scanner recovers 0/412 missing variants from the harvested content). The remaining 6 PMIDs have small (<14 KB) FULL_CONTEXT shims that explicitly reference tables, appendices, or supplements whose contents never reached the LLM.

A separate but stackable issue: 10 of the 46 PMIDs (114 missing variants) hit Azure rate-limit errors during extraction. Re-running these would not by itself recover any variants because the harvested content is still abstract-only; but they should be re-run after a successful re-harvest so the stat does not mask harvest failures.

## Reconciling counts

The task brief cited 59 blind PMIDs / 417 missing variants. The actual numbers from the `20260503_161203` run that produced `comparison_results/summary.json` are **46 PMIDs / 412 missing variants** under the strict definition above. The discrepancy is small enough to be a different snapshot or a slightly different definition (e.g., including PMIDs with non-zero variants but missing the specific gold-standard variant). The conclusions are unaffected.

## Failure categories

| Code | Category | PMIDs | Missing variants | Mechanism | Fix |
|------|----------|------:|-----------------:|-----------|-----|
| **A** | Abstract-only, LLM ran cleanly | 31 | 287 | Harvester returned only the PubMed abstract; LLM correctly returned 0 variants because abstracts of these large cohort papers do not enumerate individual variants | **Tier 3.5 browser fallback** to fetch full text + supplement tables |
| **B** | Abstract-only, LLM rate-limited | 9 | 114 | Same harvest gap as A. The Azure rate-limit error is a red herring — even a successful LLM run on the abstract would have returned 0 | Re-harvest with Tier 3.5; re-run extraction. Rate-limit alone is not the problem |
| **C** | Full-text harvested, LLM ran cleanly but extracted 0 | 5 | 7 | FULL_CONTEXT files are 2–14 KB shims; extraction notes explicitly cite "Table 1 / Appendix A / supplemental ZIP referenced but content not included" | Improve table+appendix harvesting (HTML rendering may have stripped them); add supplement ZIP unpacker |
| **D** | Full-text harvested, LLM rate-limited | 1 | 4 | Single PMID (19322600) with ~3 KB shim that hit the rate limiter | Re-run extraction; investigate whether shim contains the variants |
| | **Total** | **46** | **412** | | |

### Why the regex scanner found 0 of 412 missing variants

The scanner ran on whatever text was available (FULL_CONTEXT for 6 PMIDs, raw abstract JSON for 40). It did not recover a single missing gold variant. This is consistent with the diagnosis: the variants live in **tables, supplements, or full-text body that was never harvested**. Building a smarter scanner is not the lever that closes this gap; getting the right text into the pipeline is.

## Top 10 highest-impact blind PMIDs

These 10 PMIDs alone account for **325 of 412 (79%)** of the blind-PMID gap.

| Rank | PMID | Missing | Category | Title |
|----:|------|--------:|----------|-------|
| 1 | [15840476](https://pubmed.ncbi.nlm.nih.gov/15840476/) | 86 | A | Compendium of cardiac channel mutations in 541 consecutive unrelated patients referred for long QT syndrome genetic testing. |
| 2 | [10973849](https://pubmed.ncbi.nlm.nih.gov/10973849/) | 59 | A | Spectrum of mutations in long-QT syndrome genes. KVLQT1, HERG, SCN5A, KCNE1, and KCNE2. |
| 3 | [26496715](https://pubmed.ncbi.nlm.nih.gov/26496715/) | 53 | B | Common Genotypes of Long QT Syndrome in China and the Role of ECG Prediction. |
| 4 | [11854117](https://pubmed.ncbi.nlm.nih.gov/11854117/) | 44 | A | Increased risk of arrhythmic events in long-QT syndrome with mutations in the pore region of the human ether-a-go-go-related gene potassium channel. |
| 5 | [16922724](https://pubmed.ncbi.nlm.nih.gov/16922724/) | 20 | B | Spectrum of pathogenic mutations and associated polymorphisms in a cohort of 44 unrelated patients with long QT syndrome. |
| 6 | [23631430](https://pubmed.ncbi.nlm.nih.gov/23631430/) | 12 | B | Results of genetic testing in 855 consecutive unrelated patients referred for long QT syndrome in a clinical laboratory. |
| 7 | [17905336](https://pubmed.ncbi.nlm.nih.gov/17905336/) | 11 | A | Long QT and Brugada syndrome gene mutations in New Zealand. |
| 8 | [19996378](https://pubmed.ncbi.nlm.nih.gov/19996378/) | 9 | A | Clinical characteristics and genetic background of congenital long-QT syndrome diagnosed in fetal, neonatal, and infantile life: a nationwide questionnaire surv |
| 9 | [12402336](https://pubmed.ncbi.nlm.nih.gov/12402336/) | 9 | A | DHPLC analysis of potassium ion channel genes in congenital long QT syndrome. |
| 10 | [21216356](https://pubmed.ncbi.nlm.nih.gov/21216356/) | 8 | B | Genetic testing of patients with long QT syndrome. |

**Pattern:** every one of the top 10 is a multi-gene LQTS cohort/screening study where individual variant tables live in the published full text or supplement, not the abstract. PMID 15840476 (Tester et al., *Heart Rhythm* 2005) is the single largest blind paper at 86 missing variants — a 541-patient compendium with the variant table almost certainly in the body.

## Per-PMID diagnostic table (all 46)

Columns: `pmid` | `miss` (gold-variant rows missing) | `full_kb` (FULL_CONTEXT.md size) | `figs` | `sup` (supplement files) | `status` (extraction status) | `scan_rec` (regex scanner gold variants recovered) | `title`

| PMID | Miss | Full KB | Figs | Sup | Status | Scan✓ | Title |
|------|----:|--------:|----:|----:|--------|------:|-------|
| 15840476 | 86 | 0.0 | 0 | 0 | abstract_only | 0 | Compendium of cardiac channel mutations in 541 consecutive unrelated patients referred for |
| 10973849 | 59 | 0.0 | 0 | 0 | abstract_only | 0 | Spectrum of mutations in long-QT syndrome genes. KVLQT1, HERG, SCN5A, KCNE1, and KCNE2. |
| 26496715 | 53 | 0.0 | 0 | 0 | rate-limited | 0 | Common Genotypes of Long QT Syndrome in China and the Role of ECG Prediction. |
| 11854117 | 44 | 0.0 | 0 | 0 | abstract_only | 0 | Increased risk of arrhythmic events in long-QT syndrome with mutations in the pore region  |
| 16922724 | 20 | 0.0 | 0 | 0 | rate-limited | 0 | Spectrum of pathogenic mutations and associated polymorphisms in a cohort of 44 unrelated  |
| 23631430 | 12 | 0.0 | 0 | 0 | rate-limited | 0 | Results of genetic testing in 855 consecutive unrelated patients referred for long QT synd |
| 17905336 | 11 | 0.0 | 0 | 0 | abstract_only | 0 | Long QT and Brugada syndrome gene mutations in New Zealand. |
| 19996378 | 9 | 0.0 | 0 | 0 | abstract_only | 0 | Clinical characteristics and genetic background of congenital long-QT syndrome diagnosed i |
| 12402336 | 9 | 0.0 | 0 | 0 | abstract_only | 0 | DHPLC analysis of potassium ion channel genes in congenital long QT syndrome. |
| 21216356 | 8 | 0.0 | 0 | 0 | rate-limited | 0 | Genetic testing of patients with long QT syndrome. |
| 22727609 | 8 | 0.0 | 0 | 0 | rate-limited | 0 | Clinical characteristics of 30 Czech families with long QT syndrome and KCNQ1 and KCNH2 ge |
| 26823142 | 6 | 0.0 | 0 | 0 | abstract_only | 0 | Pediatric Cohort With Long QT Syndrome　- KCNH2 Mutation Carriers Present Late Onset But Se |
| 28532774 | 6 | 0.0 | 0 | 0 | abstract_only | 0 | Genotype Positive Long QT Syndrome in Patients With Coexisting Congenital Heart Disease. |
| 17210839 | 6 | 0.0 | 0 | 0 | abstract_only | 0 | Prevalence of long-QT syndrome gene variants in sudden infant death syndrome. |
| 19843919 | 5 | 0.0 | 0 | 0 | rate-limited | 0 | Latent genetic backgrounds and molecular pathogenesis in drug-induced long-QT syndrome. |
| 16379539 | 4 | 0.0 | 0 | 0 | abstract_only | 0 | Gene sequencing in neonates and infants with the long QT syndrome. |
| 29020304 | 4 | 0.0 | 0 | 0 | abstract_only | 0 | Identification of a targeted and testable antiarrhythmic therapy for long-QT syndrome type |
| 26063740 | 4 | 0.0 | 0 | 0 | rate-limited | 0 | Follow-up of 316 molecularly defined pediatric long-QT syndrome patients: clinical course, |
| 18508782 | 4 | 0.0 | 0 | 0 | abstract_only | 0 | Sudden arrhythmic death syndrome: familial evaluation identifies inheritable heart disease |
| 22429796 | 4 | 0.0 | 0 | 0 | abstract_only | 0 | End-recovery QTc: a useful metric for assessing genetic variants of unknown significance i |
| 15851119 | 4 | 0.0 | 0 | 0 | abstract_only | 0 | Identification of a common genetic substrate underlying postpartum cardiac events in conge |
| 19322600 | 4 | 2.9 | 0 | 0 | rate-limited | 0 | Contribution of Long-QT Syndrome Genetic Variants in Sudden Infant Death Syndrome |
| 20975234 | 3 | 0.0 | 0 | 0 | abstract_only | 0 | Atrioventricular block-induced Torsades de Pointes with clinical and molecular backgrounds |
| 20197117 | 3 | 0.0 | 0 | 0 | abstract_only | 0 | Congenital long QT syndrome and 2:1 atrioventricular block: an optimistic outcome in the c |
| 28438721 | 3 | 0.0 | 0 | 0 | abstract_only | 0 | Clinical profile and mutation spectrum of long QT syndrome in Saudi Arabia: The impact of  |
| 30246897 | 3 | 0.0 | 0 | 0 | abstract_only | 0 | Patients diagnosed with long QT syndrome after repair of congenital heart disease. |
| 12808265 | 3 | 0.0 | 0 | 0 | abstract_only | 0 | Q-T peak dispersion in congenital long QT syndrome: possible marker of mutation of HERG. |
| 29121719 | 3 | 0.0 | 0 | 0 | abstract_only | 0 | Implantable cardioverter defibrillator therapy in repaired tetralogy of Fallot after pulmo |
| 19352046 | 2 | 3.2 | 0 | 0 | ran_full | 0 | Paper 19352046 |
| 19127321 | 2 | 0.0 | 3 | 0 | abstract_only | 0 | Pregnancy and the risk of torsades de pointes in congenital long-QT syndrome. |
| 14642687 | 2 | 13.7 | 0 | 0 | ran_full | 0 | Paper 14642687 |
| 29016797 | 2 | 0.0 | 0 | 0 | abstract_only | 0 | Multiple clinical profiles of families with the short QT syndrome. |
| 26129877 | 2 | 0.0 | 0 | 0 | rate-limited | 0 | Functional Characterization of Rare Variants Implicated in Susceptibility to Lone Atrial F |
| 15466642 | 2 | 0.0 | 0 | 0 | rate-limited | 0 | Spectrum and frequency of cardiac channel defects in swimming-triggered arrhythmia syndrom |
| 18348270 | 1 | 0.0 | 0 | 0 | abstract_only | 0 | Delineation of the phenotype associated with 7q36.1q36.2 deletion: long QT syndrome, renal |
| 17595376 | 1 | 0.0 | 0 | 0 | abstract_only | 0 | Neonatal long QT syndrome due to a de novo dominant negative hERG mutation. |
| 27555138 | 1 | 0.0 | 0 | 0 | abstract_only | 0 | Sevoflurane-associated torsade de pointes in a patient with congenital long QT syndrome ge |
| 30844837 | 1 | 0.0 | 0 | 0 | abstract_only | 0 | Arrhythmogenic Ventricular Cardiomyopathy Associated With Fibromuscular Dysplasia of Ostia |
| 18452873 | 1 | 2.0 | 0 | 0 | ran_full | 0 | Paper 18452873 |
| 16969682 | 1 | 0.0 | 0 | 0 | abstract_only | 0 | [Long QT syndrome causing grand mal epilepsy: case report, pedigree, therapeutic options,  |
| 20636320 | 1 | 0.0 | 0 | 0 | abstract_only | 0 | QTc prolongation and family history of sudden death in a patient with desmin cardiomyopath |
| 23237912 | 1 | 0.0 | 0 | 0 | abstract_only | 0 | A case of long QT syndrome with triple gene abnormalities: digenic mutations in KCNH2 and  |
| 26412604 | 1 | 3.1 | 0 | 1 | ran_full | 0 | Paper 26412604 |
| 11332568 | 1 | 0.0 | 0 | 0 | abstract_only | 0 | Electrocardiographic prediction of abnormal genotype in congenital long QT syndrome: exper |
| 26937396 | 1 | 12.9 | 3 | 1 | ran_full | 0 | Paper 26937396 |
| 22382559 | 1 | 0.0 | 0 | 0 | abstract_only | 0 | Peripartum cardiomyopathy presenting with syncope due to Torsades de pointes: a case of lo |

## Category A details — abstract-only, clean LLM run (31 PMIDs / 287 variants)

These are the highest-leverage targets: harvest succeeded only at the abstract level; nothing about the extraction prompt or LLM model would have recovered the variants. Almost all are mid-to-large LQTS cohort papers where the variant list is in a full-text table.

Sample LLM diagnostics (these are *correct* responses given the input):
- PMID 11854117 (44 missing): abstract names "the pore region of HERG" but lists no individual variants — full text contains the per-patient table.
- PMID 17905336 (11 missing): "Long QT and Brugada syndrome gene mutations in New Zealand" — abstract gives counts only.
- PMID 17210839 (6 missing): SIDS/LQTS prevalence study — variant table not in abstract.

**Recommendation A:** Route every Category-A PMID through Tier 3.5 browser fallback. Prioritize by missing count (top of the list above).

## Category B details — abstract-only, rate-limited (9 PMIDs / 114 variants)

Same harvest gap as A, plus a transient Azure rate-limit. The rate-limit is not load-bearing: the harvested input was abstract-only, so even an unlimited-quota LLM call would not have recovered the missing variants.

- **26496715** (53 missing): Common Genotypes of Long QT Syndrome in China and the Role of ECG Prediction.
- **16922724** (20 missing): Spectrum of pathogenic mutations and associated polymorphisms in a cohort of 44 unrelated patients w
- **23631430** (12 missing): Results of genetic testing in 855 consecutive unrelated patients referred for long QT syndrome in a
- **21216356** (8 missing): Genetic testing of patients with long QT syndrome.
- **22727609** (8 missing): Clinical characteristics of 30 Czech families with long QT syndrome and KCNQ1 and KCNH2 gene mutatio
- **19843919** (5 missing): Latent genetic backgrounds and molecular pathogenesis in drug-induced long-QT syndrome.
- **26063740** (4 missing): Follow-up of 316 molecularly defined pediatric long-QT syndrome patients: clinical course, treatment
- **26129877** (2 missing): Functional Characterization of Rare Variants Implicated in Susceptibility to Lone Atrial Fibrillatio
- **15466642** (2 missing): Spectrum and frequency of cardiac channel defects in swimming-triggered arrhythmia syndromes.

**Recommendation B:** Same as A — Tier 3.5 browser fallback. Once a real full-text harvest exists, re-run extraction with current rate-limit-respecting settings.

## Category C details — full text but zero variants extracted (5 PMIDs / 7 variants)

These are the **only true extraction-stage failures** in the blind set. In each case the LLM clearly diagnosed the problem in its own `challenges` / `notes` fields:

### PMID 19352046 — 2 missing | 3278 bytes | 0 figs | 0 sup
Missing variants: `['A561V', 'D501N']`

LLM challenges:
- No specific KCNH2 variant details (e.g. c. or p. notation) provided in the available text, abstract, or tables
- Only generic reference to 'the KCNH2 mutation' in two cases

### PMID 14642687 — 2 missing | 14052 bytes | 0 figs | 0 sup
Missing variants: `['G572C', 'N470D']`

LLM challenges:
- Specific variant notations (e.g. c. or p.) not listed in provided full text
- Tables are described but actual mutation data rows not included
- Only generic reference to 'six distinct HERG mutations' and 'eight subjects'

### PMID 18452873 — 1 missing | 2026 bytes | 0 figs | 0 sup
Missing variants: `['G189Ins']`

LLM challenges:
- No individual KCNH2 variants or mutation notations are provided in the main text or abstract
- Only high-level genotype category counts (e.g., 174 LQT2 patients) are given without variant details

### PMID 26412604 — 1 missing | 3162 bytes | 0 figs | 1 sup
Missing variants: `['R1047L']`

LLM challenges:
- No specific variant details or notations provided for KCNH2 (only gene-level count of n=3 pathogenic mutations)
- Supplemental file referenced but no content available in input

### PMID 26937396 — 1 missing | 13179 bytes | 3 figs | 1 sup
Missing variants: `['R835fsX']`

LLM challenges:
- Specific nucleotide or amino acid change for the KCNH2 variant is not stated in the provided full text
- Table 1 and Appendix A are referenced but their content is not included in the prompt
- Supplemental ZIP file not accessible

**Pattern:** in every case the FULL_CONTEXT shim mentions Table 1 / Appendix A / a supplemental file by reference, but the table rows or supplement contents are not in the prompt. The HTML or PDF-to-Markdown converter is stripping the variant tables.

**Recommendation C:**
1. Re-harvest with the supplement fetcher (these papers have non-zero supplement directories for 2 of 5).
2. Audit the full-text format converter (likely `harvesting/format_converters.py`) on these specific PMIDs to see why the table content is missing — this is a small enough set (5 papers) to inspect manually.
3. For PMID 26412604 and 26937396, the supplement was fetched (`sup=1`) but apparently not included in FULL_CONTEXT — verify the supplement-to-prompt path.

## Category D details — full text, rate-limited (1 PMID / 4 variants)

PMID 19322600 (4 missing: K897T, P408A, R1047L, R148W). FULL_CONTEXT is 2.9 KB — likely a stub. Extraction was rate-limited but a fallback "abstract-only" extraction returned `confidence=low`, `abstract_only=False` and 0 variants. The harvested text is too small to contain the per-patient variant data anyway, so this is effectively another harvest gap.

**Recommendation D:** Treat as Category B — Tier 3.5 browser fallback first; the rate-limit is secondary.

## Tier 3.5 browser fallback — projected impact

Of the 46 blind PMIDs / 412 missing variants:

| Cohort | PMIDs | Missing variants | Browser fallback impact |
|--------|------:|----:|-------------------------|
| Category A (abstract-only, clean) | 31 | 287 | **Direct fix** — full text very likely available via publisher site, OA repository, or PMC; some may be paywalled |
| Category B (abstract-only, rate-limited) | 9 | 114 | **Direct fix** — same as A; rate-limit is incidental |
| Category D (rate-limited stub) | 1 | 4 | **Likely fix** — current 3 KB harvest is a stub |
| Category C (full text, table parser failed) | 5 | 7 | **Indirect** — needs format-converter / supplement-extractor work, not browser |
| **Total addressable by Tier 3.5** | **41** | **405** | |

**Top candidates for browser fallback (by impact):**

| Rank | PMID | Missing | Category | Why browser fallback should help |
|----:|------|--------:|----------|----------------------------------|
| 1 | [15840476](https://pubmed.ncbi.nlm.nih.gov/15840476/) | 86 | A | Cohort/screening paper; variant table in main body |
| 2 | [10973849](https://pubmed.ncbi.nlm.nih.gov/10973849/) | 59 | A | Cohort/screening paper; variant table in main body |
| 3 | [26496715](https://pubmed.ncbi.nlm.nih.gov/26496715/) | 53 | B | Cohort/screening paper; variant table in main body |
| 4 | [11854117](https://pubmed.ncbi.nlm.nih.gov/11854117/) | 44 | A | Cohort/screening paper; variant table in main body |
| 5 | [16922724](https://pubmed.ncbi.nlm.nih.gov/16922724/) | 20 | B | Cohort/screening paper; variant table in main body |
| 6 | [23631430](https://pubmed.ncbi.nlm.nih.gov/23631430/) | 12 | B | Cohort/screening paper; variant table in main body |
| 7 | [17905336](https://pubmed.ncbi.nlm.nih.gov/17905336/) | 11 | A | Cohort/screening paper; variant table in main body |
| 8 | [19996378](https://pubmed.ncbi.nlm.nih.gov/19996378/) | 9 | A | Cohort/screening paper; variant table in main body |
| 9 | [12402336](https://pubmed.ncbi.nlm.nih.gov/12402336/) | 9 | A | Cohort/screening paper; variant table in main body |
| 10 | [21216356](https://pubmed.ncbi.nlm.nih.gov/21216356/) | 8 | B | Cohort/screening paper; variant table in main body |

## Summary recommendations (priority-ordered)

1. **Run Tier 3.5 browser fallback against the 40 abstract-only blind PMIDs (Categories A + B + D = 405 missing variants).** Order by missing count. The top 10 alone are worth ~325 variants. This is the single largest closeable chunk of the recall gap (~9–10 percentage points if every fetch succeeds).

2. **Audit the format converter for Category C (5 PMIDs).** Run the existing harvest on these 5 PMIDs in isolation and dump the intermediate XML/HTML before markdown conversion to see where the variant tables drop out. Specific PMIDs: 14642687, 19352046, 18452873, 26412604, 26937396. Especially check why `_supplements/` content for 26412604 and 26937396 is not in FULL_CONTEXT.

3. **Re-extract Category B + D after re-harvest.** The 10 rate-limited PMIDs are not a rate-limit problem alone — they need full text first. Once Tier 3.5 has produced full text, re-run extraction.

4. **(Optional) Add a "harvest grade" column to the run manifest.** Tag each PMID as `abstract_only` / `partial_fulltext` (<5 KB) / `fulltext` so this diagnostic can be reproduced from the manifest in future runs without cross-referencing the DB.

## Methodology

1. Loaded `comparison_results/missing_in_sqlite.csv` (816 missing rows, 164 unique PMIDs).
2. Pulled `papers.pmid` (586) and `DISTINCT variant_papers.pmid` (397) from `results/KCNH2/20260503_161203/KCNH2.db`.
3. Defined "blind" as: `pmid IN papers AND pmid NOT IN variant_papers AND pmid IN missing_in_sqlite`. Result: 46 PMIDs / 412 missing variants.
4. For each blind PMID, inspected: `pmc_fulltext/{pmid}_FULL_CONTEXT.md`, `_figures/`, `_supplements/`, `abstract_json/{pmid}.json`, `extractions/KCNH2_PMID_{pmid}.json` (`extraction_metadata.abstract_only`, `challenges`, `notes`), `pmid_status/extraction_failures.csv`.
5. Ran `utils.variant_scanner.scan_document_for_variants` against the harvested content (FULL_CONTEXT preferred, abstract fallback) and counted gold-standard variant overlap.
