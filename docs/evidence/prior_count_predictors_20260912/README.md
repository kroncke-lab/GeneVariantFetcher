# Priors, affected/unaffected counts, and genomic predictors

Date: 2026-09-12. **For this requested analysis, every `gnomad_added`
observation is assumed unaffected.** This is a descriptive analysis of the
five frozen protocol fits, with no model refit or extraction change.

Large unaffected totals do not imply many low-risk variant rows. Most of those
observations belong to a few variants, while most remaining rows have one or
a few affected observations and no unaffected observations. BRCA2 already has
many near-zero priors; the one-affected-carrier update moves them toward 0.09.
The other genes' high shared priors largely survive their sparse updates.

## Figures

- [Prior, posterior, and affected-fraction distributions; affected and unaffected count distributions](PRIOR_COUNT_DISTRIBUTIONS.png)
- [Prior versus affected count, unaffected count, and affected fraction](PRIOR_VS_COUNTS.png)
- [Prior versus AlphaMissense, GPN-Star-M, and AlphaGenome AVI](PRIOR_VS_PREDICTORS.png)

Each row is a gene. Histograms give every retained variant key equal weight;
probability bins have width 0.04. Count categories retain zero and explicitly
separate literature unaffected from total unaffected. Scatterplot count axes
retain zero, becoming logarithmic above one. Within each figure, corresponding
axes share ranges across genes. Missing predictor scores are omitted, never
converted to zero; no score range is clipped.

## Definitions and prior weight

The exact row set is `analysis/{GENE}_protocol/posterior_variants.csv` in the
archived `grant_e2e_20260909` evidence. It contains retained literature variant
keys with usable counts, not all possible or all population-observed variants.
The totals below are for these final histogram rows, after protocol exclusions.
Counts summed over variants are variant-carrier observations, not unique people.

```text
A = affected
U_literature = literature_n - affected
U_gnomAD = gnomad_added                  [assumed unaffected]
U = U_literature + U_gnomAD
n = A + U
count-only affected fraction = A / n

p0 = insilico_prior_mean
S = 10
prior equivalent affected count = 10 * p0
prior equivalent unaffected count = 10 * (1 - p0)
analytic posterior mean = (A + 10*p0) / (A + U + 10)
nominal prior weight in final update = 10 / (n + 10)
```

The prior's equivalent counts are model weights, not observed people. The prior
was fitted using these same constructed counts and the features truncation,
AlphaMissense, and AlphaMissense missingness. Thus the displayed update weight
does not measure a fraction of independent evidence. Archived posterior means
were sampled; their maximum absolute difference from the analytic formula is
0.00736. Figures preserve the archived means.

| Gene | Variants | Prior median [IQR] | Affected | Literature U | gnomAD U | Variants with U=0 | Prior weight >50% |
|---|---:|---:|---:|---:|---:|---:|---:|
| HNF1A | 233 | .829 [.189, .878] | 761 | 115 | 179,610 | 67.4% | 85.0% |
| GCK | 391 | .928 [.824, .944] | 923 | 82 | 1,067 | 81.1% | 94.4% |
| LDLR | 871 | .931 [.839, .939] | 4,564 | 216 | 17,079 | 69.6% | 84.7% |
| BRCA2 | 5,962 | .0096 [.0029, .2554] | 10,678 | 2,550 | 2,006,331 | 73.4% | 95.4% |
| KCNQ1 | 691 | .735 [.206, .797] | 8,202 | 2,174 | 7,225 | 56.6% | 78.7% |

Median affected counts are 1, 1, 1, 1, and 2 in this gene order; median unaffected
count is zero in every gene. Median variant-level `A/n` is therefore **1.000
in all five genes**. This differs sharply from pooled `sum(A)/sum(n)`:
0.00422, 0.445, 0.209, 0.00529, and 0.466.

The top five variants account for **99.7%, 93.1%, 90.3%, 96.4%, and 82.2%** of
gnomAD observations, respectively. BRCA2 V2466A alone contains 72.7% of its
gnomAD observations. These totals are highly concentrated, not distributed
across all rows. `counts_summary.json` records the largest rows and exact counts.

## Why near-zero priors do not stay near zero

For one affected carrier and no unaffected carriers, the analytic update is
`(1 + 10*p0)/11`. As the prior approaches zero this approaches **0.0909**.
For example, a prior of 0.003 becomes 0.0936. BRCA2 has 3,748 such rows.
Consequently, 50.2% of BRCA2 priors are below 0.01, but only 2.3% of its saved
posteriors are below 0.01. Its median moves from 0.0096 to 0.1037.

The same mechanism moves a BRCA2 truncating prior of 0.2554 to approximately
0.3231 with one affected carrier, explaining the second broad mode without
requiring an intermediate-risk population for every such variant.

For HNF1A, GCK, and LDLR, nontruncating rows missing AlphaMissense receive common
priors of **0.947, 0.962, and 0.956** (18, 47, and 152 rows). Missing a score is
therefore associated with very high fitted risk in those arms. The corresponding
BRCA2 and KCNQ1 constants are 0.643 and 0.206. Truncating priors form additional
horizontal bands; see `feature_strata.csv`. Seven BRCA2 truncating rows also
carry AlphaMissense values, an existing classification inconsistency retained
here and documented by the previous audit.

Neither the prior model nor this update imposes a disease-prevalence reference
level. A distribution of all variants would also require adding variants that
never appeared in the selected literature set. That missing population-only
universe is distinct from the zero-to-0.09 update effect.

## Prior versus AlphaMissense, GPN-Star, and AlphaGenome

"GMT-star" is interpreted as **GPN-Star**, matching the available local predictor
contracts. The primary GPN field is `gpn_star_m447_llr_calibrated`, the fixed M
model used by Variant Browser. Lower/more negative signed LLR means greater
predicted effect. Its negative correlations below therefore indicate agreement
in direction with the prior. The M model was not selected by the correlation.
The pinned score artifact is revision
`5c799b2ec6aa089f0caa8294ae72adb4510f81ae`.
[Primary GPN-Star score source](https://huggingface.co/datasets/songlab/gpn-star-scores).

The AlphaGenome scalar is the Atlas **AVI_SCORE**, stored as `alphagenome_avi`.
Higher means greater predicted impact; it can be negative or exceed one and is
not a penetrance probability. AVI combines AlphaGenome and AlphaMissense
information, so it is not independent of AlphaMissense.
[Primary AlphaGenome Atlas description](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/).

The following Spearman correlations use the **same single-allele missense rows
for all three predictors within each gene**, avoiding different coverage in the
comparison. The scatter figure separately shows all available rows per panel;
its annotations therefore have different n and correlations.

| Gene | Same rows n | Prior vs AlphaMissense | Prior vs GPN-Star-M | Prior vs AlphaGenome AVI |
|---|---:|---:|---:|---:|
| HNF1A | 139 | +1.000 | −.794 | +.783 |
| GCK | 244 | +1.000 | −.535 | +.755 |
| LDLR | 426 | +1.000 | −.690 | +.748 |
| BRCA2 | 3,167 | +1.000 | −.666 | +.681 |
| KCNQ1 | 391 | +1.000 | −.720 | +.741 |

AlphaMissense's perfect rank agreement is mechanical: within a fixed
nontruncating, score-present stratum the prior is a monotonic fitted
transformation of AlphaMissense. It does not validate penetrance. The five
different curves show that the same AlphaMissense score can imply very different
priors across genes because the fitted count relationships differ.

GPN-Star and AVI broadly track the prior, but also vary among variants assigned
the same prior. This makes them candidates for testing whether the current
feature model discards useful distinctions. Agreement with the existing prior
alone does not establish that either would improve outcome prediction.

## Join and provenance

The read-only warehouse extraction reuses the frozen feature table's
`vf_variant_ids`, checks ID existence and gene membership, and preserves the
allele identity details. It does not remap ambiguous protein keys. Primary
GPN/AVI scalars require exactly one mapped allele and no version/value conflict;
multi-allele means/minima/maxima are retained only for later sensitivity checks.
AVI mirrors were checked against authoritative AVI_SCORE context records;
there were no context conflicts for queried IDs. These checks do not repair
upstream protein-to-allele mapping errors.

GPN/AVI coverage in the final row set is 163/233, 274/391, 508/871, 3,619/5,962,
and 466/691. AlphaMissense in every primary analysis is the **archived model
input**, not a newly fetched version. This matters because current BRCA2
warehouse AlphaMissense sources sometimes disagree. GPN source versions,
warehouse fetch timestamps, input hashes, conflict records, and AVI's absent
upstream version identifier are retained in `predictor_sources.json`, the
compressed allele table and `predictor_version_conflicts.csv.gz`. No missing
version is interpreted as a known current
release.

## Recommended next experiment

Keep the requested gnomAD-as-unaffected assumption fixed. First quantify and
review the high-prior, missing-AlphaMissense strata and shared truncation bands;
validate ambiguous/misclassified allele keys before adding their new scores.
Then compare an independently trained or cross-fitted prior using GPN-Star and
AVI with the current AlphaMissense-only feature arm, using the same held-out
outcome rows and grouped cohort/family splits. Check calibration and predictive
scores, not correlation with the old prior. AVI should not be treated as
independent corroboration of AlphaMissense.

Run a prespecified sensitivity to prior strength and the influence of the few
common variants in prior fitting. If the desired estimand is the distribution
across population-observed variants, add an explicit population-only stratum
with its gnomAD unaffected counts and no literature affected observations,
preserving the distinction between absent literature evidence and measured
clinical outcomes. This tests the missing-universe explanation without
manufacturing a desired histogram shape or choosing a prevalence target after
seeing the plot.

## Reproduce and inspect

From the GVF root, using a Python environment with numpy, pandas, scipy and
matplotlib:

```bash
python docs/evidence/prior_count_predictors_20260912/compute_counts.py
MPLCONFIGDIR=tmp/matplotlib python docs/evidence/prior_count_predictors_20260912/build_figures.py
```

These two commands require only committed evidence. The optional
`collect_predictors.py` refresh needs the local variantFeatures warehouse and
the frozen sibling feature tables; it performs indexed, read-only queries.
Current warehouse results may change, so committed snapshots are the inputs for
reproducing this report. Gzip outputs use deterministic `mtime=0` compression.

- `gene_summary.csv`, `count_distributions.csv`, `probability_distributions.csv`:
  aggregate distribution tables.
- `variant_counts.csv.gz`: every retained key's observed counts, prior equivalent
  counts, weights, analytic mean, and archived posterior mean.
- `feature_strata.csv`, `alphamissense_comparison.csv`, `evidence_strength.csv`:
  feature strata and prior/count influence summaries.
- `predictor_scores.csv.gz`, `predictor_alleles.csv.gz`: frozen predictor join,
  ambiguity/status flags and detailed genomic alleles.
- `predictor_correlations.csv`, `figure_check.json`, `counts_summary.json`:
  numerical checks and statistics underlying the figures.

Validation checks unique row identities, nonnegative integer counts,
`n=A+U`, stored versus recomputed affected fractions, histogram membership,
S=10, allele/gene membership and score conflicts. All three figures were
rendered and visually inspected. This evidence-only analysis changes no
extraction or model defaults.
