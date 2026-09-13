# Fixed-weight cross-gene decomposition

This audit reconstructs the saved neighborhood features without changing
counts, priors, variant identities or weights. HNF1A, LDLR and KCNQ1 use the
frozen missense structural extension; GCK uses the frozen class-specific
monomer replay; BRCA2 uses the corrected local-3D/canonical-polymer primary.

For shared gene-by-missense prior mean `mu` and strength `S`, every donor has
posterior mean `(mu*S + A)/(S+n)`. The exact decomposition is:

```
prior_retention_i = sum_j W_ij * S/(S+n_j)
prior_component_i = mu * prior_retention_i
counts_component_i = sum_j W_ij * A_j/(S+n_j)
density_i = prior_component_i + counts_component_i
raw_count_neighborhood_i = sum_j W_ij * A_j/n_j
```

`prior_share_of_density` is the direct algebraic prior component divided by
density. It is not the fraction of a causal effect or the amount the prior
raises a feature. Shrinkage can move a neighborhood down or up: corrected
BRCA2 has median raw-count neighborhood 41.43% and median posterior-derived
density 17.98%, whereas GCK has 34.91% versus 37.41%. Medians of separate
components need not sum to the median density.

| Gene | Supported / missense units | Median density | Median direct prior share | A=0 control examples meeting rule |
|---|---:|---:|---:|---:|
| HNF1A | 726 / 785 | 6.35% | 77.3% | 5 |
| GCK | 614 / 634 | 37.41% | 60.7% | 3 |
| LDLR | 1,206 / 1,376 | 16.25% | 48.5% | 28 |
| BRCA2 | 6,326 / 6,656 | 17.98% | 60.2% | 68 |
| KCNQ1 | 975 / 1,061 | 24.62% | 56.4% | 8 |

The diagnostic control rule is **A=0, U>=10, own posterior<=1%, density>=10%**.
These thresholds select clear examples for interpretation; they are not new
model gates or biological benign labels. Population carriers retain the user's
unaffected assumption. Complete diagnostics carry both the rule flag and its
population-only subset. The examples file contains up to five unique units per
gene, ordered by density minus own posterior.

Examples make the distinction concrete. GCK D217N has A=0 and U=227, own
posterior 0.47%, but neighbor density 33.00%. Its density consists of 22.14
percentage points from the shared prior and 10.86 points from neighbor
affected counts. HNF1A E275A has A=0 and U=66, own posterior 0.34%, but
neighbor density 41.76%; most of that feature comes from other variants'
affected counts. Neither is contradictory: a target's own counts are
deliberately excluded from its neighborhood feature. Neither density is
therefore the target's measured penetrance.

Corrected BRCA2 covers 1,937 variants through local 3D and 4,389 through
polymer, leaving 330 unavailable. The other genes retain their original
source-specific coverage. `allgene_density_source_coverage.csv` preserves the
structured/polymer/mixed/missing breakdown, and unsupported rows retain
missing decomposition values throughout.

## Validation and limits

Every saved weight is independently reapplied to the unchanged donor counts.
GCK's wide matrices and the older genes' long matrices are reconstructed
directly. Corrected BRCA2 sparse batches are checked against
`analysis/BRCA2/cache_manifest.json`, including source hashes, exact variant
order, every weight shard hash and the full positive-pair count. Duplicate
pairs, self-donors, negative/nonfinite weights, incomplete identities and
unsupported targets carrying weight fail the audit. All five reconstructed
densities agree within 7.1e-15. Complete per-variant diagnostics are sharded
below 1.2 MB.

This validates arithmetic and provenance, not disease calibration. The
underlying internal LOO uses posterior labels from the selected count
collection and fixed gene-by-class empirical hyperparameters. Conditional
intervals omit source, endpoint and modeling uncertainty. The interpretation
must retain the distinction among own empirical posterior, posterior-derived
neighborhood feature, and independently validated disease risk.

Reproduce after any formatter changes to the audit script so its hash remains
current:

```
../BayesianPenetranceEstimator/.venv/bin/python docs/evidence/structural_sanity_20260913/audit_all_genes.py --require-corrected-brca2
```

`allgene_audit_receipt.json` pins all inspected source files and records the
exact definitions. `allgene_summary.csv`, `strongest_control_examples.csv`
and `variant_diagnostics/<GENE>.part*.csv.gz` provide the source-backed outputs.
