# Population-inclusive GCK structural run

The completed run uses the full-locus empirical prior and adds observed
population-only missense variants to the structural donor pool. Alpha receives
affected counts; beta receives literature-unaffected plus all gnomAD carriers,
as requested. The primary gnomAD carrier count is allele count minus homozygotes.

The common prior is fitted upstream from **14,553 GCK count-bearing units**,
including 14,181 population-only units and all eligible consequence classes.
Its parameters are alpha = 0.042974, beta = 2.117582, mean = **1.989%**, and
strength = 2.160555. Structural eligibility is applied afterwards. This produces
**634 canonical-WT-matching missense units**: 385 population-only, 166
literature-only, and 83 represented in both sources.

The primary 1V4S biological monomer supports **614 of 634 targets**. The median
variant-excluded density is **18.239%**, compared with about 83.8% in the earlier
literature-only pilot. These runs differ in the prior, donor universe, and
refreshed counts; that change cannot be attributed solely to geometry.

## Fixed method and identity

The engine uses each union `unit_id` as the excluded variant identity. Clinical
aggregate members have already been consolidated in the audited union. Distinct
population genomic alleles remain distinct donors, including alleles at the
same protein position or producing the same protein substitution. Only the
target unit is excluded, and alternative units at the same residue remain.
Every donor contributes its empirical posterior mean, with equal-variant spatial
weights and no extra carrier-count multiplier.

The primary kernel is `2 / (1 + exp(log(3) * d / 3))`, with no distance or
sequence cutoff. The saved primary donor tables contain **278,640 positive-weight
pairs beyond 20 Å**, among 376,382 total donor pairs. Median Kish effective donor
count is 22.88; the median normalized weight contributed beyond 20 Å is 1.04%.

The frozen, numbered biological-monomer geometry comes from the previous pilot's
`geometry/` folder, whose source/mapping evidence is preserved there. Primary
distances use side-chain heavy-atom mass-weighted centers, with C-alpha fallback
for glycine. IDR pairs use `3.8 * sqrt(canonical sequence separation)` only within
the same contiguous segment and chain. Ambiguous/missing geometry and mixed
structured/IDR pairs remain unavailable. Unavailable targets retain missing
density rather than zero.

Twenty prespecified scenarios cover 1V4S, 1V4T, and the AlphaFold monomer;
COM/C-alpha metrics; half-distance 2/3/5 Å; and the earlier documented missing-loop
control and polymer-parameter sensitivity. Both experimental states support 614
targets; AlphaFold supports 628. No bandwidth selection was performed.

## Predictive comparison

The all-common comparison fits and evaluates **610 targets** with available
AlphaMissense and structural/sequence features. The full 634-unit eligible donor
pool remains available. The separate **original 242** comparison restricts both
fitting and evaluation targets to the exact earlier clinical cohort and uses its
archived AlphaMissense values, while retaining the expanded donors and refreshed
empirical posteriors.

For every outer held-out variant, that identity is also removed from every
training row's donor feature; each training row's own identity remains excluded.
The efficient monomer calculation was checked against the engine for three
held-out identities. Full-locus empirical hyperparameters remain fixed by user
choice, so the result is not an independent refit of the entire prior pipeline.

| Targets | Model | MAE vs empirical posterior | MSE vs observed fraction | Mean Beta-binomial negative log score |
|---|---|---:|---:|---:|
| All common, 610 | AlphaMissense | 0.147045 | 0.179127 | 1.000828 |
| All common, 610 | AM + density | 0.146755 | 0.178283 | 0.994985 |
| Original 242 | AlphaMissense | 0.138995 | 0.237735 | 1.547250 |
| Original 242 | AM + density | 0.139427 | 0.237913 | 1.548626 |

The expanded comparison gives a small numerical improvement from adding density;
the controlled original clinical cohort becomes slightly worse. The sequence
control is competitive. These diagnostics do not establish a meaningful added
predictive benefit from the structural feature. The count score describes the
adopted observed-count experiment, not an external population-outcome study.
Origin and Kish-support strata are retained in the metrics table.

AlphaMissense uses an unambiguous exact genomic-member score for 467 units and
the archived clinical-key fallback for 163; four remain missing. Version/value
conflicts do not trigger fallback or arbitrary version selection. The 242-target
control deliberately preserves the earlier archived AM source.

## Prior sensitivity with geometry and donors fixed

These calculations retain the exact primary weights and all 634 donor identities.
They change only empirical prior parameters or the explicitly named gnomAD count
proxy, then update each donor posterior and recompute its weighted density.

| Prior/count arm | Prior mean | Prior strength | Median own posterior | Median density | Mean absolute density change |
|---|---:|---:|---:|---:|---:|
| Primary: full locus, historical MSE | 1.989% | 2.1606 | 1.360% | 18.239% | 0 |
| Coding/splice empirical prior | 23.742% | 2.3201 | 16.591% | 29.838% | 11.430 percentage points |
| Full locus, normalized weighted MSE | 1.989% | 0.1800 | 0.303% | 31.728% | 13.567 percentage points |
| Full locus, allele-count proxy | 1.986% | 2.1564 | 1.357% | 18.251% | 0.0124 percentage points |

The prior's population and variance convention materially affect the density.
Changing carrier count to allele count has a much smaller effect in this donor
set. The primary remains the full-locus historical-MSE carrier-count analysis;
these alternatives are sensitivities and have not been used to retune the
regression or structural kernel.

## Artifacts and verification

- `GCK_structural_input_variants.csv.gz` records counts, empirical posterior
  parameters, union identity/members, origins, and the AM source policy.
- `GCK_primary_variant_density.csv` contains primary density, support, and
  conditional intervals. Scenario tables and normalized donor matrices are
  sharded; all primary donor pairs are retained in chunks of at most 40 targets.
- `GCK_variant_loo_predictions.csv.gz`, `GCK_variant_loo_metrics.csv`, and
  `GCK_original_242_metric_comparison.csv` preserve the fitted comparisons.
- `GCK_fixed_geometry_prior_sensitivity.csv.gz` has one density column per prior
  arm; the companion summary preserves parameters and effect sizes.
- The two PNG figures show density distributions, sources, residue profiles,
  state sensitivity, and the predictive comparisons.
- `run_checks.json`, `artifact_checks.json`, and
  `prior_geometry_sensitivity_checks.json` preserve inputs, hashes, and checks.

The conditional intervals use 8,192 donor Beta draws with one shared draw per
variant across targets, conditioning on fixed counts, empirical hyperparameters,
identity, geometry, and the kernel. They do not include uncertainty in population
ascertainment, shared pedigrees, endpoint labels, or structure selection. The
underlying GCK clinical evidence still pools disease mechanisms/endpoints.

Independent reconstruction from the saved donor rows recovers every supported
density within 3.9e-15 and every conditional variance within 1.1e-16. Maximum
Monte Carlo mean discrepancy is 2.80 simulation standard errors. Both figures
were visually inspected; every artifact is smaller than 1.2 MB. Five focused
tests cover global donor removal, near-dominant donors, same-residue alternatives,
held-out outcome independence, controlled target sets, and AM conflict handling.

Reproduce using the scientific Python environment already available locally:

```bash
/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/.venv/bin/python docs/evidence/population_inclusive_penetrance_20260912/run_structure.py
/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/.venv/bin/python docs/evidence/population_inclusive_penetrance_20260912/prior_geometry_sensitivity.py
/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/.venv/bin/python docs/evidence/population_inclusive_penetrance_20260912/validate_structure.py
/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/.venv/bin/python -m pytest -q docs/evidence/population_inclusive_penetrance_20260912/test_structure_statistics.py
```
