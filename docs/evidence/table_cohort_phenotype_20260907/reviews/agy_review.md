REVIEW_VERDICT

Your core diagnosis is correct: the bottleneck is translating parsed per-variant frequencies into explicit phenotype counts. However, your proposed classifier risks violating the strict semantic boundaries required by the curation contract. Here is my adversarial breakdown.

### Q1: Contract Defensibility and Tier 1 vs Tier 2
Tier 1 is defensible. If a table structurally defines its columns as "No. of Patients" or its caption as "Mutations in LQT1 Patients," extracting those counts as `affected` is a direct read of the authors' structural claim, not a manufactured partition. You are mapping the explicitly stated cohort role to the lab's nomenclature.

Tier 2 crosses the line. Projecting a paper-level ascertainment (e.g., from the abstract) onto an isolated supplementary table that just lists "COUNT" requires a biological and methodological assumption that the table strictly represents the proband cohort and not, for example, a mixed validation cohort or allele frequency list. I would not ship Tier 2; it introduces high risk for low-quality recall.

### Q2: Adversarial Shapes and Exclusion Tuning
**Shapes that will fool Tier 1:**
*   **Mixed trait tables:** "Table 3: Clinical characteristics of KCNQ1 variant carriers" (if "carriers" is missed, it maps everyone to affected, missing asymptomatic relatives).
*   **In vitro/assay tables:** Tables where the "N" refers to transfected cells (HEK293), patch-clamp experiments, or simulated structures, not human patients.
*   **Molecular autopsy:** "Variants identified in the SUD cohort" (Sudden Unexplained Death). If SUD is not in your target disease list, it fails; if it is, it might include unaffected individuals who died of other causes.

**Exclusion words:**
*   **Missing:** "in vitro", "cells", "assays", "HEK", "Xenopus", "cohort" (often implies a mix), "SCD/SUD", "relatives", "family members", "sub-clinical".
*   **Too aggressive:** "phenotype" (excluding "Phenotypes of BrS patients" breaks a perfect Tier 1 candidate) and "genotype-positive" (in clinical series, genotype-positive individuals are literally the cases).

### Q3: Controls and Alleles
Emitting `affected = 0` for explicit control tables is logically sound and aligns with the lab's convention, provided the source text explicitly defines the cohort as healthy or normal.

However, mapping "reference alleles" directly to human counts is mathematically flawed. An allele count is not a person count. If a table reports 44 reference alleles, that could be 22 people (homozygous) or 44 people (heterozygous). In dominant cardiac traits, variants are typically heterozygous, but asserting `carriers = unaffected = N` from an *allele* count introduces factual inaccuracy unless the parser explicitly handles zygosity. You should leave allele counts as NULL for people counts.

### Q4: Measurement Plan Critique
Your within-arm ablation (LLM off/on for the deterministic module) is excellent. It perfectly isolates the deterministic logic from the LLM's run-to-run stochasticity.

**Flaws in the plan:**
1.  **Wasteful baseline:** If the only difference is a code-owned deterministic module applied *post-extraction*, paying for a separate baseline LLM arm is a waste of money. Execute one extraction run, and score it twice (once with the derivation module bypassed, once with it applied).
2.  **Loose exactness threshold:** A conditional exactness drop from 92% to 85% means you are willing to inject up to a 7% error rate into previously pristine supplied data just to boost recall. I would require conditional exactness to remain >= 90%. If it drops below 90%, your exclusions are failing.

### Q5: Highest-Leverage Next Step
The single highest-leverage fix is addressing the **fixed-width identity fast path** that currently bypasses phenotype tables (like in PMID 30059973).

If your pipeline identifies 185 variants but intentionally skips reading the tables containing their actual clinical counts, you are mathematically guaranteeing 0% recall for those rows. Fixing a routing bug to ensure identified variants are actually passed to the table parser is lower risk and higher reward than building complex semantic classifiers. The 232 identity misses are also critical; you cannot count what you cannot find.

### Q6: Unproven Claims
**Claim:** "conditional exactness stays around 90% (the 45 multi-center duplicate rows become wrong affected because their carriers are already 'wrong' against gold's row-splitting)."

**Challenge:** You are assuming the *only* major source of new errors will be those 45 known edge-case duplicates. This is unproven. It is highly likely that Tier 1 (and especially Tier 2) will misclassify tables containing mixed proband/relative cohorts as pure cases, injecting novel extraction errors.

**How to check:** During your free calibration run on the open locks, do not just look at the aggregate 90% number. You must manually audit *every single newly supplied affected count* that the scorer flags as wrong. If they are all the 45 duplicates, your claim holds. If you find newly misclassified tables, your exclusion logic requires tightening before spending the API budget.
