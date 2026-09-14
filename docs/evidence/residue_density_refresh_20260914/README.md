# Corrected missense penetrance density by residue — 2026-09-14

The five profiles have been rebuilt after additional source and endpoint
corrections. The corrected inputs support these **descriptive neighborhood
features**; they do not establish calibrated disease probabilities. In
particular, a high neighborhood score does not establish that every substitution
at that residue causes the named disease.

[Five-gene overview PNG](plots/ALL_GENES_RESIDUE_DENSITY.png) ·
[Vector PDF](plots/ALL_GENES_RESIDUE_DENSITY.pdf)

| Gene / intended endpoint | Detailed figure | Residue values | Supported residues / protein length | Median plotted residue score |
|---|---|---|---:|---:|
| HNF1A / MODY3 | [PNG](plots/HNF1A_RESIDUE_DENSITY.png) · [PDF](plots/HNF1A_RESIDUE_DENSITY.pdf) | [CSV](tables/HNF1A_residue_density.csv) | 426 / 631 | 6.34% |
| GCK / GCK-MODY, mild hyperglycemia | [PNG](plots/GCK_RESIDUE_DENSITY.png) · [PDF](plots/GCK_RESIDUE_DENSITY.pdf) | [CSV](tables/GCK_residue_density.csv) | 339 / 465 | 34.74% |
| LDLR / familial hypercholesterolemia | [PNG](plots/LDLR_RESIDUE_DENSITY.png) · [PDF](plots/LDLR_RESIDUE_DENSITY.pdf) | [CSV](tables/LDLR_residue_density.csv) | 637 / 860 | 15.81% |
| BRCA2 / hereditary breast/ovarian cancer susceptibility | [PNG](plots/BRCA2_RESIDUE_DENSITY.png) · [PDF](plots/BRCA2_RESIDUE_DENSITY.pdf) | [CSV](tables/BRCA2_residue_density.csv) | 2,481 / 3,418 | 3.47% |
| KCNQ1 / LQTS1 | [PNG](plots/KCNQ1_RESIDUE_DENSITY.png) · [PDF](plots/KCNQ1_RESIDUE_DENSITY.pdf) | [CSV](tables/KCNQ1_residue_density.csv) | 508 / 676 | 23.98% |

Each detailed figure contains the primary posterior neighborhood, a log-scale
comparison with count diagnostics, each variant's **own** count posterior, and
affected/unaffected counts. The blue line is the mean of supported variant-only
leave-one-out scores at a residue. The pale band is the between-variant minimum
and maximum, not an uncertainty interval. Positions without observed eligible
missense variants or usable donor geometry are gaps, not zeros or interpolated
predictions. Geometry tracks distinguish 3D, polymer, mixed and unestimated
positions. Exact-zero diagnostic values have explicit log-axis markers below
all positive values.

## What was fixed

These are exact observation-level changes, not blanket removal of publications
or variants. Valid evidence elsewhere for the same allele is retained. Unknown
phenotypes never become unaffected observations.

- **[HNF1A source audit](source/HNF1A/README.md):** remove reviewed type 2
  diabetes, somatic liver-tumor and unresolved mixed-diabetes contributions from
  the MODY3 primary. Correct measured glycemic partitions among germline
  relatives, retaining genuine MODY evidence from those papers.
- **[GCK source audit](source/GCK/README.md):** separate activating/hypoglycemia,
  homozygous neonatal diabetes and mosaic/unknown-phenotype observations from
  heterozygous GCK-MODY. Restore the documented mild-hyperglycemia parents and
  genuine MODY relatives. Diagnosis-only diabetes observations from PMID
  36208030 are preserved in a labeled sensitivity. Nineteen canonical-WT
  mismatches remain excluded pending identity evidence.
- **[LDLR source audit](source/LDLR/README.md):** quarantine 83 all-type
  myocardial-infarction case assignments from PMID 29802317 because FH status
  was not established. Correct 23 clinically affected, mutation-positive
  relatives in the Greek FH table from unaffected to affected. Retain the
  independently checked clinical FH case counts from PMID 31491741.
- **[BRCA2 source audit](source/BRCA2/README.md):** quarantine the remaining
  227 canonical missense affected assignments previously flagged for review:
  205 from a pan-tumor table without available germline data, and 22 from a
  variant compilation without an endpoint-specific affected/unaffected split.
  Across all types these decisions quarantine 508 affected assignments. They
  are not all proven somatic variants; the required germline/phenotype evidence
  is unavailable. The five previously verified breast-cancer missense carriers
  from PMID 40664060 remain.
- **[KCNQ1 source audit](source/KCNQ1/README.md):** a 492-carrier founder total
  cannot supply a binary LQTS phenotype partition. Cardiac-event/event-free
  splits among LQTS patients are retained only in a labeled sensitivity, not
  treated as LQTS present/absent. Remove a whole-cohort total assigned to one
  modifier allele and two duplicate emissions. Retain the independently
  reproduced clinical/ECG-assessed LQTS case table.

The source reviews are bounded. Remaining papers, repeated founder families,
multivariant genotypes, age, treatment and cohort ownership require further
curation. The [initial metadata screen](source/other_genes/README.md) documents
its scope and why LDLR source ranking was narrowed to actual canonical missense
keys. Neither original source databases nor previous frozen evidence was edited.

## Priors and observed counts

The table below describes the **final canonical missense union**, after clinical
identity checks and population allele reconciliation. These are observation
counts, not deduplicated people across publications. All gnomAD carriers are
assumed unaffected as requested. Population-only variants are included.

| Gene | Missense units | Affected | Literature unaffected | gnomAD unaffected | α empirical | β empirical | Prior mean | Own posterior if A=0, U=1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HNF1A | 781 | 400 | 42 | 1,721,270 | 0.228733 | 1.996270 | 10.28% | 7.09% |
| GCK | 614 | 603 | 32 | 4,369 | 1.026827 | 1.956145 | 34.42% | 25.78% |
| LDLR | 1,361 | 1,346 | 95 | 108,024 | 0.371949 | 2.138991 | 14.81% | 10.59% |
| BRCA2 | 4,444 | 1,858 | 1,482 | 1,516,774 | 0.163563 | 4.122739 | 3.82% | 3.09% |
| KCNQ1 | 1,061 | 5,845 | 1,594 | 22,640 | 0.638317 | 1.527879 | 29.47% | 20.16% |

The singleton column is an example, not a universal lower bound; variants with
many unaffected observations have lower posteriors. Missense and nonsense
priors are fitted separately; nonsense donors never enter these plots. All ten
fits, previous/current comparisons and the historical MSE/MAE diagnostics are
in [prior_summary.csv](tables/prior_summary.csv). No synonymous, frameshift,
splice or noncoding variants enter either requested type's prior.

Clinical source totals and final union totals need not coincide. For KCNQ1,
three incompatible shared-allele keys exclude 8 A/1 U, while two canonical-cDNA
keys resolve into missense and add 5 A/0 U: net −3 A/−1 U from the source table.
BRCA2 releases four population carriers from former nonsense protein aggregates
into their DNA-confirmed frameshift units. All remain in the complete union;
this is consequence reconciliation, not loss of population observations. The
[independent validation](validation/README.md) preserves both traces.

## Sanity checks and interpretation

**The broad BRCA2 baseline is mainly prior shrinkage, not a long-distance tail.**
Its median supported *variant* neighborhood falls from 4.57% in the preceding
catalog-corrected run to 3.36% now; the median plotted *residue mean* is 3.47%.
There are 3,954 zero-affected units, including 2,012 unaffected singletons.
An unaffected singleton still has a 3.09% own posterior under the empirical
prior. Equal-variant averaging repeats those small positive donor posteriors;
1.5 million population observations do not become 1.5 million neighborhood
votes. The median direct prior-component share of the BRCA2 primary score is
70.7%.

At h=3 Å, BRCA2's median normalized weight beyond 20 Å is **0.474%**; the
maximum is 4.70%. Narrowing h to 2 Å lowers its median variant score to 3.04%
and sharpens peaks, but the minimum remains 0.161%. Increasing h to 5 Å raises
the median to 3.69%. This does not support attributing the residual few-percent
floor primarily to distant donors. Some individual peaks remain bandwidth
sensitive. [All five bandwidth comparisons](tables/bandwidth_sensitivity.csv)

**Zero-case diagnostics now visibly reach zero or below 0.1%.** Forty BRCA2
targets at 28 residues have no affected donor in any eligible context: their
raw-fraction neighborhood is exactly zero, whereas their posterior neighborhood
is 0.161–3.094%. The one-prior pooled-count diagnostic reaches 0.00126% among
supported variant targets. These are different estimands, not alternate
calibrations of one probability. The primary is deliberately not clipped to
zero or forced below 0.1%. A residue with no own cases can also have genuinely
affected spatial or same-segment polymer donors.

**GCK does not imply that every residue causes MODY.** Removing the reviewed
opposite or unresolved endpoints lowers the missense prior from 37.04% to
34.42% and the median variant neighborhood from 37.41% to 34.86%. Among 614
missense units, 393 have zero affected observations, including 247 unaffected
singletons. Those singletons each retain a 25.78% posterior, and the median
prior-component share of the neighborhood is 61.7%. Thus a broad elevated line
is mathematically expected under this particular prior/posterior construction.
It is not evidence that every substitution has that MODY risk. The own-posterior
and count panels make the difference inspectable.

**Large source-count corrections need not produce large density changes.**
KCNQ1's median variant neighborhood barely changes, 24.615% to 24.611%, despite
the corrected founder and event counts. Its prior fits variant fractions with
saturating weights, and its neighborhood averages variants rather than people.
An already count-rich variant can change little, and a residue without usable
geometry remains unestimated even if its own counts are large. Numerical
stability does not excuse the original endpoint error or establish calibration.

The [sanity summary](tables/sanity_summary.csv) records each gene's zero-case
donors, singleton posterior, prior retention/component share and distant weight.
The [GCK diabetes-proxy sensitivity](tables/GCK_diabetes_proxy_sensitivity.json)
and [KCNQ1 cardiac-event sensitivity](tables/KCNQ1_cardiac_events_sensitivity.json)
remain separately labeled. They mix the extra endpoint into the primary inputs;
neither is a pure independently phenotyped cohort. An HNF1A broader-diabetes
ledger is preserved as source evidence but was not used for a plotted fit.

## Declared calculation

For each gene and variant type, n=A+U, y=A/n and w=1−1/(n+0.01), using observed
units with n≥1. The historical empirical mean is μ=Σwy/Σw and variance is
v=Σw(y−μ)²/M, where M is the number of variant units, **not** Σw. Then
κ=μ(1−μ)/v−1, α=μκ and β=(1−μ)κ. The posterior is
(α+A)/(α+β+A+U). **Alpha adds affected; beta adds literature and gnomAD
unaffected.** The full-dataset prior stays fixed during variant leave-one-out.

The primary target score is ΣW times donor posterior means. The excluded
identity includes aliases and every structural copy; other variants at the
same residue stay eligible. Raw K(d)=2/[1+exp(log(3)d/h)] has half weight at
h=3 Å and a strictly positive tail. We also run h=2 and h=5 Å. Structured
distances use mass-weighted side-chain heavy-atom centers (glycine CA).
Disorder uses 3.8√|residue separation| only within the same canonical IDR
segment and chain. Missing structure alone does not make a residue disordered.

HNF1A uses the canonically numbered 8PI8 dimer; GCK the 1V4S biological monomer;
LDLR the existing full monomer AlphaFold model; KCNQ1 the partial 9U7F tetramer.
BRCA2 combines local AlphaFold/experimental frames with candidate IDR polymer
contexts. Unaligned frames never supply cross-frame distances. Each donor has
at most its nearest eligible copy per context; supported contexts are averaged
equally. BRCA2 is not a solved full biological assembly, LDLR interdomain
placement retains the prior PAE caveat, and partial structure/IDR assumptions
limit biological interpretation. Geometry is unchanged from the previous
canonical audits and pinned in the run input manifests.

The gold diagnostic averages A/n under the same normalized W. For the purple
diagnostic, each context first sums **unnormalized** K-weighted A and U, adds
one prior, then averages its supported context scores equally. This avoids
repeating one prior per donor but changes the estimand and its dependence on
absolute kernel mass; it is not silently substituted for the primary. The
saved tables also include the corresponding raw pooled-count fraction.

## Verification, reviews and reproduction

[Independent numerical and figure validation](validation/README.md) checks the
full population union, exact observation corrections, all ten priors, every
primary weight shard, selected independent geometry/context calculations,
target exclusion, same-residue donors, posterior decomposition and plotted
residue aggregation. Figure layout checks cover positive values within log
limits and visible labels. [Local six-figure review](reviews/local_plot_review.md)
checks the actual rendered PNGs and source tables.

**Agy, Grok and Claude CLIs all completed fresh consultations.** Their scope was
the generic method and synthetic examples; the actual records, fitted gene
counts and plots were audited locally. No app fallback was needed. Several
reviewer arithmetic suggestions were wrong and were rejected after explicit
checks. Exact prompts, raw replies, usage and the
[review disposition](reviews/REVIEW_DISPOSITION.md) are preserved. These
consultations are not independent clinical validation.

Run from the repository root with the BPE scientific Python environment:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
../BayesianPenetranceEstimator/.venv/bin/python docs/evidence/residue_density_refresh_20260914/run_refresh.py
../BayesianPenetranceEstimator/.venv/bin/python docs/evidence/residue_density_refresh_20260914/run_refresh.py --genes GCK --sensitivity diabetes_proxy --halves 3
../BayesianPenetranceEstimator/.venv/bin/python docs/evidence/residue_density_refresh_20260914/run_refresh.py --genes KCNQ1 --sensitivity cardiac_events --halves 3
../BayesianPenetranceEstimator/.venv/bin/python docs/evidence/residue_density_refresh_20260914/plot_residues.py
../BayesianPenetranceEstimator/.venv/bin/python docs/evidence/residue_density_refresh_20260914/summarize_sanity.py
```

The source subfolders contain their own deterministic ledger builders. Run
their existing source checks before rebuilding counts. Primary raw weights are
in ignored `results/` and can be regenerated; hashes are recorded in each
`analysis/<GENE>/density_checks.json`. Compact evidence, code and all figures
are pinned by [manifest.json](manifest.json); [freeze.py](freeze.py) verifies
working, staged and committed bytes. This analysis changes no production parser
or source DB, and does not claim an extraction benchmark improvement.

The next scientific step is endpoint/person/family curation followed by a fixed,
independent comparison of prior/posterior, raw-count and structure/sequence
features. Choose bandwidth or a different prior based on that comparison,
not on making the plotted line resemble an expected disease rate. The active
forward checklist remains [TASKS.md](../../../TASKS.md).
