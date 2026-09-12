# Grok statistical review — 2026-09-12

Scope: user-supplied aggregate table and narrative plus generic mathematical questions only. Grok did not inspect internal source or CSVs. The external model review is advisory and does not establish the implementation or biological validity.

**The table cannot support the narrative’s causal story, and it cannot show a missing spike at background risk.** It only describes the distribution of *some estimated quantity* over *literature/warehouse variants already selected into the set*. BRCA2 is shifted left relative to the other four genes; that is all the aggregates establish.

**What the numbers actually show.** Variant counts differ by more than an order of magnitude (HNF1A 233 vs BRCA2 5962). Medians: GCK/LDLR 0.94, HNF1A 0.84, KCNQ1 0.71, BRCA2 0.10. Share below 0.10: BRCA2 49%, HNF1A 16%, others 2–4%. Share ≥0.50 is high for four genes (66–87%) and low for BRCA2 (11%). Share ≥0.90 is high only for GCK/LDLR (62–63%); HNF1A 16%, KCNQ1 4%, BRCA2 0%. Those contrasts are real *as descriptions of this selected set*. They are not penetrance calibration, not a variant census, and not evidence about gnomAD.

---

**(1) A histogram of literature-observed variants cannot represent all variants.**
The sampling frame is variants that were sequenced, interpreted, and ingested. That frame over-represents affected-proband alleles, ClinVar-submitted alleles, and gene-specific testing panels. It under-represents benign, synonymous, deep-intronic, and genotype-first alleles. Most true variants in any gene sit near population risk for the linked disease; they never enter this histogram. BRCA2’s larger *n* is more literature density, not a complete catalog. The feeling that there “should be more ~0 variants” is the expected shape of an *unascertained* variant set. This table cannot test that expectation. **Hypothesis, not shown:** adding genotype-first / population-only alleles would create a large near-background mass in HNF1A, GCK, LDLR, and KCNQ1.

**(2) A large *summed* population allele count does not imply large *per-variant* denominators.**
Gene-level AC is dominated by a few common alleles. Distinct variants are mostly rare. BRCA2 having many variants plus a large gene-level carrier sum (an external claim, not in the table) is compatible with *most* variants remaining singletons/doubletons. The sentence “most BRCA2 variants appear a handful of times against a large population denominator” conflates gene-level or common-allele denominators with the denominator of the typical variant. **Not established.** The left shift in BRCA2 could equally be prior mass, endpoint definition, missense mix, or many tiny denominators plus a low prior—not “two million carriers doing the work.”

**(3) Phenotype-unknown gnomAD samples are not confirmed unaffected. Allele count is not distinct carriers.**
Unknown phenotype ≠ unaffected. Treating unknowns as non-cases biases estimated penetrance downward (denominator inflation, possible case misclassification). AC ≠ people: homozygotes, relatedness, overlapping cohorts, sex chromosomes, and multi-allelic sites break that identity. “Carriers” in a narrative total is not a statistical unit the table provides. **Hypothesis:** BRCA2’s mass below 0.10 is partly this coding choice. **Unknowns:** per-variant AC, person-level carrier counts, phenotype filters, relatedness handling.

**(4) Yes — a count-weighted feature prior can move common-allele influence onto unanchored rare variants.**
If truncating / domain / functional classes are pooled by allele or carrier counts, high-AC alleles dominate the fitted class mean. Rare alleles then inherit that mean with little of their own data. That is shrinkage to a *count-weighted* class, not to the typical rare allele of the class. The narrative’s own caveat (grey bars, truncating prior) is this mechanism stated informally. The table cannot show that the warehouse uses count weighting; it is a **design risk**, not a finding.

**(5) Generic Beta mean \((10 p_{\text{prior}} + a)/(10 + N)\), one affected, mathematical only.**
Prior weight 10 means ten pseudo-observations at \(p_{\text{prior}}\). One real affected is weak.

- \(N=1\) (proband only): \(p=0.94 \to 10.4/11 \approx 0.945\) (prior barely moves). \(p=0.003 \to 1.03/11 \approx 0.094\).
- \(N=6\) (proband + 5 coded unaffected): \(0.94 \to 10.4/16=0.65\); \(0.003 \to \approx 0.064\).
- \(N=101\): \(0.94 \to \approx 0.094\); \(0.003 \to \approx 0.009\).

So a high class prior plus small \(N\) parks variants near 1 regardless of sparse counts. A near-null prior plus one proband parks near ~0.09 even with no population data, and near 0 as \(N\) grows if extras are coded unaffected. **Implication for the table, as illustration only:** GCK/LDLR medians of 0.94 are the *expected parking place* of a high prior with small \(N\), not proof of near-complete penetrance. BRCA2’s median 0.10 is numerically close to “one affected, null-ish prior, \(N=1\)” *or* “high prior washed by a large unaffected-coded \(N\)”—the table cannot tell which. Do not treat this formula as a verified description of the warehouse.

**(6) Near zero is not background absolute risk, and not comparable across genes.**
Lifetime/age/sex-specific risk differs by endpoint: hyperglycemia vs familial hypercholesterolemia vs breast/ovarian cancer vs arrhythmic events. A BRCA2 value of 0.10 is not “population breast-cancer risk.” A GCK value of 0.94 is not the same kind of “disease” as LDLR CAD or KCNQ1 sudden death. Age truncation, screening, male/female mix, and competing mortality are unidentified. **Below 0.10 ≠ population rate.** Cross-gene comparison of these histograms is not a comparison of comparable absolute risks.

**(7) Visual spread and ClinVar do not validate absolute penetrance.**
Spread can be mixed priors, mixed \(N\), mixed endpoints, missing ClinVar, and stacking. ClinVar is an assertion layer correlated with the same literature that feeds the estimator (circularity). KCNQ1 is *not* shown to be a calibrated control: median 0.71, 66% ≥0.50, 4% ≥0.90 is still upper-half pile-up, not an empirical penetrance law. “Most resembles a genuine distribution” is a **visual hypothesis**. Grey-bar dominance at the high end, if present in the plots, would argue the opposite of validation: position tracking the prior, as the caveat already says.

**LDLR “spike at 0.94 = no-ClinVar frameshifts the warehouse does not enumerate” is not in the table.** Neither is ClinVar-class composition, nor that grey = truncating.

---

**(8) Ranked audit before changing a standing protocol**
Do not retune priors, denominators, or display rules from this table. Ordered checks; each must pass with pre-registered acceptance, or the protocol stays.

1. **Define the estimand.** Absolute penetrance of a named endpoint, by age/sex, in a named population — or a relative enrichment score. Acceptance: written estimand; if it is not absolute penetrance, stop calling the axis penetrance.
2. **Sampling frame.** List inclusion rules. Report the fraction of ClinVar / gnomAD / panel variants *absent* from the histogram. Acceptance: frame documented; no claim about “all variants.”
3. **Unit of analysis.** Per-variant \(N\) = distinct people, not AC. Split homozygotes, relatedness, overlapping datasets. Acceptance: AC-vs-carrier table; median and IQR of per-variant \(N\) by gene (the missing piece for claims (2)–(3)).
4. **Phenotype coding.** Sensitivity: gnomAD as unknown (excluded from \(N\)) vs unaffected vs multiple-imputation. Acceptance: BRCA2 and GCK location of mass must be stable in direction under the unknown-not-unaffected rule, or the left/right contrast is coded, not biological.
5. **Prior influence.** For each variant, report posterior, likelihood-only estimate, prior mean, and effective prior weight. Acceptance: variants with \(N \le 5\) labeled prior-dominated; truncating/no-ClinVar mass near 0.94 must shrink when prior weight is halved, or the spike is the prior.
6. **Count-weighted class priors.** Re-fit class priors unweighted and inverse-frequency-weighted. Acceptance: rare-variant posteriors do not move by more than a pre-set amount (e.g. 0.10 median shift) when common alleles are down-weighted; else common alleles are leaking (question 4).
7. **Background risk.** Overlay age/sex-specific population incidence for *that endpoint*, not 0. Acceptance: “near 0” is never claimed unless it matches that overlay; BRCA2 0.10 is not called population rate.
8. **External calibration, not ClinVar shape.** Restricted to variants with published segregation or genotype-first cohorts. Acceptance: calibration slope/intercept with CIs; ClinVar stack and visual spread are not acceptance criteria.
9. **Only then** consider protocol change, with a locked analysis plan and a holdout gene.

---

**Hypotheses vs table.**
Supported: BRCA2’s estimated-quantity distribution is left-shifted; GCK/LDLR concentrate at the top bin; variant *n* is largest for BRCA2.
Unsupported: large per-variant denominators; gnomAD-as-unaffected; grey = truncating prior; LDLR 0.94 mechanism; KCNQ1 as gold-standard penetrance; near-0 = population risk; literature set = all variants.

**On the user’s feeling:** For a *complete* variant list one should see a large near-background spike. This display is not that list. Absence of ~0 mass in GCK/LDLR/KCNQ1 is the expected ascertainment pattern (and/or high prior + small \(N\)), not evidence that almost every variant is highly penetrant. BRCA2’s 49% below 0.10 is the only gene that even *looks* like a low-risk pile, and 0.10 is not shown to be population incidence. Treat the histograms as maps of an estimator over an ascertained set until the audit above is done.
