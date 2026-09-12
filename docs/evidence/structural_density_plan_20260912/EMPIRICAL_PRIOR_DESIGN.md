# Empirical Beta prior design audit — 2026-09-12

This is a planning calculation, not a production model or structural fit. gnomAD counts are treated as unaffected as requested. One shared empirical prior is fitted per gene/endpoint from the full available count-bearing dataset; distinct genes/diseases are not pooled. No AlphaMissense or other predictor enters this stage.

## Recovered historical method

The user's 2020 SCN5A R Markdown specifies a saturating reliability weight and weighted intercept-only regression; its residual MSE is inverted into Beta parameters and then updated with each variant's observed counts. The later KCNH2 R Markdown uses a consistent endpoint and the same construction.

Sources (historical R Markdown was retrieved read-only for inspection; full upstream source files are not copied into this evidence folder):

- [2020 SCN5A source](https://github.com/kroncke-lab/Bayes_BrS1_Penetrance/blob/master/SCN5A-BrS1-penetrance-report.Rmd): `historical_scn5a.Rmd` lines 269–273, 286–292, 309–329. The literal source has an apparent endpoint inconsistency: line 280 assigns LQT3 to observed penetrance, whereas line 329 updates using BrS1. This is not carried forward.
- [Later KCNH2 source](https://github.com/kroncke-lab/LQTS-Penetrance-APC-MAVE-Events/blob/main/predict_penetrance_kcnh2-v7.Rmd): `historical_kcnh2.Rmd` lines 157–173 (MSE and moment inversion), 396–400 (counts and weights), 430–452 (endpoint clipping, empirical prior and empirical posterior). Source URLs are branch URLs; `historical_script_hashes.json` records SHA-256 of the inspected bytes. The named historical files refer to upstream filenames and scratch copies, not files required to reproduce this preflight.
- [2022/2023 methodological paper](https://www.gimjournal.org/article/S1098-3600%2822%2901060-7/fulltext) describes the empirical-prior → empirical-posterior → feature-based modeling sequence.

For every count-bearing variant i, define A_i=affected, U_i=literature unaffected + gnomAD, n_i=A_i+U_i, q_i=A_i/n_i, and w_i=1−1/(n_i+0.01). Then

```
mu = sum(w_i * q_i) / sum(w_i)
v_hist = sum(w_i * (q_i - mu)^2) / m
kappa = mu * (1 - mu) / v_hist - 1
alpha_empirical = mu * kappa
beta_empirical = (1 - mu) * kappa
alpha_posterior_empirical_i = alpha_empirical + A_i
beta_posterior_empirical_i = beta_empirical + U_i
posterior_empirical_mean_i = alpha_posterior_empirical_i / (kappa + n_i)
posterior_empirical_variance_i = p_i * (1-p_i) / (kappa + n_i + 1)
```

Here m is the number of count-bearing variants, exactly matching the historical `mean(w * residual^2)`. This particular weight scale and MSE denominator are part of the historical recipe: dividing by sum(w) instead changes the prior strength materially. The prior mean is a reliability-weighted variant mean, **not the pooled carrier fraction**. A variant with a million gnomAD counts approaches weight 1 in the prior fit; it does not count a million times there. Its own empirical posterior still uses all its observed counts.

Use **MSE**, not MAE, to supply the variance. MAE does not uniquely determine variance, and simply substituting MAE produces invalid negative Beta shape parameters for all five full datasets under every evaluated weight convention. Report MAE as a diagnostic. Require 0<mu<1 and 0<v<mu(1−mu); detect invalid boundaries explicitly instead of silently clipping variance or imposing S=10. All evaluated real datasets give valid parameters.

Recommended adaptation: preserve exact q=0 and q=1 observations. The later historical KCNH2 source replaces these with .0005 and .9995, but this is unnecessary for these moment fits. That literal clipping is retained as a sensitivity calculation only. Across the five full datasets it changes mu by at most .000285, strength by at most .012683, and posterior median by at most .000330.

## Full available universe, before prior-model consequence exclusions

The source feature snapshots were `/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/iterations/grant_e2e_20260909/results/{GENE}_protocol/variants_features.csv`. Their relevant count/identity columns are preserved here in `full_universe_source_counts.csv.gz`, including `is_histogram` membership; `source_count_provenance.json` records all external source hashes and the compact snapshot hash. Reproduction reads only this committed snapshot and its provenance, never the external BPE files. The original `fit_protocol.py` lines 171–173 excluded synonymous/unkeyed rows. None of these feature tables contains an unkeyed row. Restoring all count-bearing classes adds synonymous variants; the full prior must not depend on AlphaMissense availability or structural coverage.

| Gene | Full variants | Restored synonymous | Full affected | Full literature unaffected | Full gnomAD unaffected |
| --- | ---: | ---: | ---: | ---: | ---: |
| HNF1A | 234 | 1 | 762 | 115 | 240,189 |
| GCK | 398 | 7 | 932 | 82 | 3,620 |
| LDLR | 880 | 9 | 4,595 | 216 | 576,920 |
| BRCA2 | 6,106 | 144 | 11,958 | 3,731 | 5,484,999 |
| KCNQ1 | 699 | 8 | 8,392 | 2,264 | 77,099 |

All 8,148 histogram rows join one-to-one and exactly reproduce affected, unaffected, gnomAD, and total counts in the wider 8,317-row universe. All full rows have positive literature n. These files therefore contain the complete available **literature-keyed** count dataset; they do not enumerate population-only variants. Adding those is a separate universe expansion, which would require identity and denominator harmonization. Variants with n=0 receive the fitted empirical prior unchanged but do not contribute q_i or train the structural density.

Full-universe historical MSE with raw A/n:

| Gene | mu | Variance | alpha_empirical | beta_empirical | Strength | Median empirical posterior |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| HNF1A | .604917 | .064617 | 1.632431 | 1.066174 | 2.698605 | .711736 |
| GCK | .807345 | .031917 | 3.127017 | .746194 | 3.873211 | .846878 |
| LDLR | .695784 | .051618 | 2.157393 | .943270 | 3.100662 | .769971 |
| BRCA2 | .530730 | .030433 | 3.812698 | 3.371172 | 7.183869 | .588071 |
| KCNQ1 | .645763 | .057064 | 1.942929 | 1.065805 | 3.008735 | .734129 |

These differ sharply from the existing AM priors, particularly BRCA2. That is expected: a broad empirical prior with no predictor no longer assigns AM-dependent near-zero means. The reliability weight also prevents a few common variants from defining the entire prior mean.

## Method sensitivity, not interchangeable definitions

`full_universe_moments.csv` compares all classes, histogram-only, and missense-only universes. Its primary method label is `historical_formula_raw_fraction`; the historical weight/MSE formula is preserved, while exact 0/1 ratios are retained. `historical_formula_endpoint_clipped` additionally reproduces the endpoint clipping. The other methods compare normalized weighted variance, carrier weighting, and equal variant weights. Fit the primary empirical prior on all classes; structural eligibility is determined afterwards. Restricting the prior to only the mapped/missense variants changes the meaning of “full dataset.”

Normalized weighted MSE `sum(w*r²)/sum(w)` produces full-dataset strengths .599/.763/.797/1.034/1.033 (gene order above), whereas historical MSE produces 2.699/3.873/3.101/7.184/3.009. Both are valid Beta moment conversions; the historical method is not scale invariant. Keep the historical denominator explicit and use normalized variance as a labeled sensitivity. The empirical prior is a moment-based regularizer; residual variation in raw A/n includes binomial sampling noise and should not be described as a recovered noise-free biological variance.

## Structural handoff

Use each variant's empirical posterior, excluding the target variant from its own neighbor aggregation while retaining other substitutions at the same residue. Per the user's instruction, the empirical hyperparameters stay fitted to the full gene dataset. This is variant-excluded neighborhood estimation conditional on full-dataset hyperparameters; the global mean still contains each variant's aggregate contribution by design. Its influence is not assumed negligible. Do not reintroduce leave-residue-out.

Do not add the target's counts again to the **empirical posterior**; those counts are already present. If the structural signal later yields a separately calibrated feature-based prior, a final update with focal counts belongs to that distinct stage. Likewise, do not treat repeated biological-assembly copies as independent clinical observations.

## Reproducibility and checks

Run from GVF:

```
.venv/bin/python docs/evidence/structural_density_plan_20260912/full_empirical_universe.py
```

Outputs: `full_universe_totals.csv`, `full_universe_moments.csv`, `full_universe_empirical_posteriors.csv.gz`, and `full_universe_checks.json` (snapshot SHA-256, partitions, uniqueness, positivity, Beta moment roundtrip, exact endpoint handling). An optional `--output-dir` writes results elsewhere for comparison. The compact-snapshot run reproduced every candidate moment, total, and per-variant posterior exactly against the original external-input calculation. `preflight_reproducibility_checks.json` records the second run comparison and initial match to all frozen histogram counts. No production models, source counts, or sibling files changed.
