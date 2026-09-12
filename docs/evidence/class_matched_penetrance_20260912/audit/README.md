# Independent audit of priors by gene and variant type

This audit reads the frozen 180,613-unit population-inclusive union. It selects
canonical-WT-valid **missense** and **nonsense** separately within each gene.
Nonsense includes the equivalent `nonsense` and `stop_gained` source labels.
Frameshift, synonymous, splice, in-frame indel, start/stop loss, noncoding and
unannotated units are excluded from both type-specific fits. The selections
contain 10,512 missense and 788 nonsense units, with no overlap. No frozen
source, counts, prior parameters or density results were edited.

## Class-specific empirical fits

The primary calculations independently reproduce the existing historical
saturating-weight mean and weighted-MSE convention. The prior strength is
estimated; it has not been fixed to make singleton posteriors stay near the mean.

| Gene | Type | Units | Prior mean | Strength | Posterior with A=0, U=1 |
|---|---|---:|---:|---:|---:|
| HNF1A | Missense | 785 | 10.37% | 2.250 | 7.18% |
| HNF1A | Nonsense | 28 | 46.65% | 6.397 | 40.34% |
| GCK | Missense | 634 | 37.04% | 2.919 | 27.59% |
| GCK | Nonsense | 33 | 70.59% | 5.152 | 59.12% |
| LDLR | Missense | 1,376 | 15.38% | 2.634 | 11.15% |
| LDLR | Nonsense | 95 | 56.66% | 3.655 | 44.49% |
| BRCA2 | Missense | 6,656 | 14.86% | 6.107 | 12.77% |
| BRCA2 | Nonsense | 577 | 65.80% | 5.418 | 55.55% |
| KCNQ1 | Missense | 1,061 | 29.48% | 2.165 | 20.16% |
| KCNQ1 | Nonsense | 55 | 52.99% | 4.804 | 43.86% |

For prior mean μ and strength κ, an unaffected singleton has posterior mean
`μκ/(κ+1)`. Its relative reduction from the prior is exactly `1/(κ+1)`:
14.1–31.6% for these missense fits and 13.5–21.5% for nonsense. Thus the
singleton remains near its appropriate class scale, but “slightly below” is
not always a negligible change. In GCK it is 37.04% → 27.59% for missense and
70.59% → 59.12% for nonsense. This follows the estimated strength, without
swapping affected/unaffected roles or forcing a stronger prior.

## Are most variants singletons?

Counts below refer to the actual prior units. A small number of clinical units
represent multiple DNA alleles, as documented in the source union.

| Gene | Type | Population-only units | All units with n=1 | Population-only units with n=1 |
|---|---|---:|---:|---:|
| HNF1A | Missense | 640 / 785 | 348 / 785 (44.3%) | 312 / 640 (48.8%) |
| HNF1A | Nonsense | 7 / 28 | 16 / 28 (57.1%) | 4 / 7 (57.1%) |
| GCK | Missense | 385 / 634 | 332 / 634 (52.4%) | 246 / 385 (63.9%) |
| GCK | Nonsense | 7 / 33 | 16 / 33 (48.5%) | 5 / 7 (71.4%) |
| LDLR | Missense | 1,005 / 1,376 | 595 / 1,376 (43.2%) | 504 / 1,005 (50.1%) |
| LDLR | Nonsense | 33 / 95 | 39 / 95 (41.1%) | 21 / 33 (63.6%) |
| BRCA2 | Missense | 2,971 / 6,656 | 3,791 / 6,656 (57.0%) | 1,454 / 2,971 (48.9%) |
| BRCA2 | Nonsense | 57 / 577 | 293 / 577 (50.8%) | 38 / 57 (66.7%) |
| KCNQ1 | Missense | 666 / 1,061 | 455 / 1,061 (42.9%) | 371 / 666 (55.7%) |
| KCNQ1 | Nonsense | 14 / 55 | 25 / 55 (45.5%) | 11 / 14 (78.6%) |

Across missense units with any population match, exactly one gnomAD carrier is
present in 45.3%, 60.0%, 45.1%, 47.9%, and 47.8% for HNF1A, GCK, LDLR, BRCA2,
and KCNQ1 respectively. That population count differs from total `n` when
clinical observations also exist. Roughly half is a better general description
than assuming all five genes have a large singleton majority.

## The historical weighting matters

`w(n)=1−1/(n+0.01)` gives a singleton weight of **0.00990099**. A two-carrier
unit has 50.75 times that weight. Consequently the many n=1 units contribute
only 0.95–1.81% of the total mean-fit weight in missense, and 0.94–1.97% in
nonsense. The same weight suppresses affected clinical singletons and unaffected
population singletons; it is not a population-only adjustment.

The diagnostics also calculate normalized weighted MSE and ordinary equal-unit
mean/MSE. They are comparisons, not replacement choices. Equal-unit missense
means are HNF1A 10.93%, GCK 31.83%, LDLR 15.62%, BRCA2 41.72%, and KCNQ1
24.00%. BRCA2 changes the most because many affected literature singletons have
little weight in the historical mean. Ordinary unweighted moment strengths are
only 0.138–0.242 for missense; an unaffected singleton then falls much further
below its prior mean. Replacing the weighted moments with unweighted moments
therefore would not automatically produce the desired slight singleton change.

The unweighted variance of observed fractions also contains finite-count
variation; it is not a measurement-error-corrected estimate of latent
penetrance heterogeneity. This audit retains that distinction and does not fit
an alternative hierarchical model.

## Artifacts and reproduction

- `class_count_summary.csv`: exact singleton, origin, carrier and weight shares.
- `class_moment_comparison.csv`: all three formulas, fitted means/variances,
  alpha/beta/strength, singleton means and conditional intervals.
- `count_distributions.csv`: affected, unaffected, total, gnomAD and literature
  unaffected count buckets, by gene/type and source origin.
- `historical_weight_influence.csv`: the count buckets' shares of variant units
  and empirical-mean weight.
- `consequence_inventory.csv`: complete frozen class inventory and explicit
  inclusion/exclusion labels; `frozen_source_hashes.csv` pins every input shard.

Run from the GVF checkout:

```bash
../BayesianPenetranceEstimator/.venv/bin/python docs/evidence/class_matched_penetrance_20260912/audit/audit_class_counts.py
```

Count partitions, population-only A=0, unique unit IDs, class disjointness and
the weighted-variance/strength identity are checked by the recipe. All ten
primary fits and both diagnostic moment alternatives meet the Beta moment
conditions. gnomAD remains assumed unaffected. This separates variant types;
it does not perform new clinical phenotype or allele-identity curation.
