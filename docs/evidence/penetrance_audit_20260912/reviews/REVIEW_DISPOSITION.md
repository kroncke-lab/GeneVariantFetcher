# Grok review disposition — 2026-09-12

Grok `grok-4.6-build` completed the bounded conceptual review in
[`grok_review.md`](grok_review.md). The successful call used CLI
`grok 1.0.13 (5e9a58528b76) [stable]`, requested `grok-4.6`, reasoning `xhigh`,
tools disabled, web and subagents disabled, verbatim mode and a two-turn limit
(completed in one turn). The exact minimal payload is
[`grok_minimal_prompt.txt`](grok_minimal_prompt.txt); the returned structured
response and usage are in [`grok_minimal_retry_raw.json`](grok_minimal_retry_raw.json).
Two minimal-payload calls cost $0.02447728 according to CLI-reported usage;
the first stopped at a plan and was not a substantive review.

The original attempt could not create a session or reach the service in the
sandbox. Automatic approval review rejected an escalation that would have
sent internal source paths and research data. No full-context review was
executed. The approved safer alternative sent only the user's supplied
aggregate table/narrative plus generic mathematical questions. Grok did not
inspect the local CSVs or implementation. Local source findings below are our
independent verification, not discoveries established by Grok.

## Accepted and independently checked

- **Selected variant universe:** the upstream observation builder filters to
  usable literature phenotype counts; the warehouse is then joined onto this
  selected set. The protocol script additionally excludes synonymous and
  unkeyed variants. The histogram cannot establish the distribution across all
  alleles in a gene or an expected fraction near background risk.
- **Aggregate versus per-variant denominators:** BRCA2's 31 common alleles
  contribute 1,992,988 of 2,006,331 added gnomAD allele counts (99.33%). Only
  1,399/5,962 plotted BRCA2 variants have a positive added anchor; their median
  positive anchor is one. Of 2,907 variants below 0.10, 1,888 have no positive
  anchor. The aggregate two-million figure does not mean most individual
  variants meet a large measured denominator.
- **Prior-driven numerical modes:** 1,855 of the 1,888 low BRCA2 variants with
  no anchor have exactly one literature carrier. BRCA2 missense median prior
  is 0.003300 and median posterior 0.093736; frameshift/nonsense/splice share
  a median prior of 0.255400 and have median posteriors about 0.324. These
  agree with the one-affected-carrier formula `(10 * prior + 1) / 11`.
  The common-allele influence is transmitted through the fitted
  count-weighted feature prior, not only through each variant's own anchor.
- **Unknown is not unaffected:** `fit_protocol.py` adds `gnomad_ac` directly
  to literature carrier denominators and treats the added counts as
  unaffected; it maps missing AC to zero added evidence. Its formula contains
  no explicit population disease incidence, endpoint, age, sex, sampling
  fraction or censoring adjustment. The resulting score is conditional on
  these assumptions and is not validated absolute disease risk.
- **Spread is not calibration:** KCNQ1's reported comparison is 134 variants
  with at least five carriers in both sources, Spearman 0.73397 and median
  absolute difference 0.10393. Shape alone does not validate risk; shared
  literature and a shared construction can prevent this from being an
  independent calibration set. Neither ClinVar color nor lack of a joined
  ClinVar annotation establishes penetrance.

These checks used the saved `analysis/<GENE>_protocol/posterior_variants.csv`
and KCNQ1 comparison JSON in the 2026-09-09 evidence directory, plus read-only
inspection of the sibling BayesianPenetranceEstimator's
`iterations/grant_e2e_20260909/scripts/{build_observations,join_features,fit_protocol}.py`.

## Not adopted from Grok

- **"Most true variants" near baseline and a required large low-risk spike:**
  these are not established for an unspecified universe, consequence class,
  endpoint or population. The appropriate expectation depends on the
  denominator; the table does not settle it.
- **Prior-halving must make a high spike shrink:** the review's proposed
  acceptance rule is mathematically reversed for all-affected observations
  and a high fixed prior. With one affected carrier and prior 0.94, reducing
  strength from 10 to 5 raises the mean from 0.94545 to 0.95. With prior
  0.003 it raises the mean from 0.09364 to 0.16917. Plot actual prior, count
  and anchor sensitivity; do not require a predetermined movement.
- **Proposed pass thresholds:** the suggested 0.10 median sensitivity bound,
  automatic labeling at `N <= 5`, and stable direction after removing all
  unphenotyped population observations are not registered acceptance gates
  and are not justified by the supplied table. They are candidate diagnostic
  choices, not conditions for accepting a biological model.
- **Generic AC cautions:** all five genes are autosomal. Sex-chromosome
  accounting is irrelevant here. Homozygotes and dataset overlap warrant
  checking; common allele AC cannot simply be called unique unaffected
  people. Relatedness and multi-allelic issues are risks to inspect, not
  demonstrated problems in these records.
- **"Leakage" terminology:** count-weighted regression and partial pooling
  can be intentional. The demonstrated concern is a common-allele-heavy,
  unphenotyped target determining rare-variant priors. Whether any method is
  defensible requires a specified estimand and independent calibration.

Keep the standing protocol and original histograms as frozen evidence. The
next useful work is a provenance and variant-universe audit, followed by
source-audited genotype-first/phenotype-measured validation and explicitly
named sensitivity arms. Do not retune a model to manufacture a desired
histogram or silently replace the accepted protocol.
