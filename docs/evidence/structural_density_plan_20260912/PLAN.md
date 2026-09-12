# Empirical-posterior structural disease density: execution plan

Date: 2026-09-12. Status: method recovery and numerical preflight; production
structural execution has not started. This plan implements the user's four
choices and supersedes the earlier proposed raw-count donor / whole-residue
holdout experiment for this analysis lane.

## Fixed choices

1. Leave out only the target variant. Other variants at the same residue remain
   informative donors. Remove the target from every equivalent chain copy.
2. Fit one empirical Beta prior to the full available eligible count dataset
   for each gene–disease endpoint. Update it with each variant's counts. Use
   those empirical posterior means as density donors. AlphaMissense does not
   enter this step. gnomAD counts remain assumed unaffected.
3. Begin with the historical compact sine kernel, midpoint 3 Å, and run a small
   prespecified distance sensitivity. Use an explicitly documented metric
   based on the historical centroid concept; the exact old atom selection is
   not recovered, so do not claim exact geometry reproduction.
4. Use polymer distances in disordered sections and correctly numbered
   biological units. A true biological monomer is an appropriate first pilot.

## 1. Recover the empirical prior and posterior

The user's original SCN5A implementation and later KCNH2 implementation already
contain the requested empirical-prior workflow. The later KCNH2 code is the
clear endpoint-consistent reference:

- [KCNH2 R analysis](https://github.com/kroncke-lab/LQTS-Penetrance-APC-MAVE-Events/blob/main/predict_penetrance_kcnh2-v7.Rmd)
- [SCN5A R analysis](https://github.com/kroncke-lab/Bayes_BrS1_Penetrance/blob/master/SCN5A-BrS1-penetrance-report.Rmd)

For one gene–disease dataset with M eligible variants:

```text
A_j = observed affected count
U_j = literature unaffected + gnomAD count
n_j = A_j + U_j > 0
y_j = A_j / n_j
w_j = 1 - 1 / (n_j + 0.01)             [historical saturating weight]

mu = sum(w_j * y_j) / sum(w_j)         [weighted intercept-only fit]
v = sum(w_j * (y_j - mu)^2) / M        [historical weighted MSE]
kappa = mu * (1 - mu) / v - 1
alpha_empirical = mu * kappa
beta_empirical = (1 - mu) * kappa

alpha_posterior_empirical,j = alpha_empirical + A_j
beta_posterior_empirical,j = beta_empirical + U_j
posterior_empirical,j = Beta(alpha_posterior_empirical,j,
                             beta_posterior_empirical,j)
p_emp,j = (alpha_empirical + A_j) / (kappa + n_j)
var_emp,j = p_emp,j * (1-p_emp,j) / (kappa + n_j + 1)
```

The empirical alpha/beta are common within a gene–disease dataset; the
posterior alpha/beta differ by variant. The prior strength is estimated from
the moments, not fixed at the previous S=10.

**MSE supplies the variance-like quantity; MAE is reported separately.** MAE is
not itself a variance. The historical MSE divides by M, not by sum(w). We will
retain that deliberate convention in the primary arm and report normalized
weighted MSE as a sensitivity. The historical convention can make a stronger
prior than normalized weighted variance; it is not an unbiased estimator of
latent between-variant variance after subtracting binomial sampling noise.

Use exact observed fractions, including 0 and 1, for the primary moment fit.
The old endpoint replacements (0.0005 and 0.9995) are a reproduction sensitivity,
not changes to the observed counts. Missing outcomes do not become zero-risk
or unaffected observations. Require 0 < mu < 1 and 0 < v < mu(1-mu), and record
an explicit failure if those Beta moment conditions fail; do not hide an invalid
fit by clipping alpha/beta. All five numerical preflights satisfy the conditions.

### Full dataset means before density eligibility filtering

Use all count-bearing consequence classes when estimating the common empirical
prior, including synonymous rows omitted from the earlier histogram. Apply
structural donor eligibility afterwards. Do not estimate the population prior
from just missense, mapped, or high-confidence structural rows.

The archived feature/count tables contain 8,317 such rows. They reproduce every
count in the 8,148-row histogram subset and restore 169 synonymous rows. This is
the full available count-bearing dataset, still selected through the literature;
it does not manufacture additional population-only variant rows.

| Gene | Full variants | mu | alpha_empirical | beta_empirical | kappa |
|---|---:|---:|---:|---:|---:|
| HNF1A | 234 | .604917 | 1.63243 | 1.06617 | 2.69860 |
| GCK | 398 | .807345 | 3.12702 | .74619 | 3.87321 |
| LDLR | 880 | .695784 | 2.15739 | .94327 | 3.10066 |
| BRCA2 | 6,106 | .530730 | 3.81270 | 3.37117 | 7.18387 |
| KCNQ1 | 699 | .645763 | 1.94293 | 1.06581 | 3.00873 |

These are preflight values from the recovered weighted-MSE method using exact
fractions. In particular, this is a different prior from the old gene-specific
AlphaMissense curves; it should not be expected to reproduce their histogram.

## 2. Variant-level leave-one-out density

Use the historical equal-variant density:

```text
D_i = sum over j != i [K(distance(i,j)) * p_emp,j]
      / sum over j != i [K(distance(i,j))]
```

There is no extra carrier-count multiplier in the primary spatial average.
Carrier evidence already informs each donor's empirical posterior; the original
`func_dist_seq.R` then gives variants equal spatial voting weight. Carrier- and
precision-weighted alternatives may be reported as sensitivities, not silently
substituted for this definition. This also avoids giving a common allele millions
of votes simply because its population count is large.

Identity is a canonical variant identifier, not a residue number. Thus, when
predicting R100W, R100Q is retained with same-residue distance zero; every R100W
chain copy and duplicate record is excluded. Duplicate representations of one
biological variant must be consolidated before LOO. Distinct genomic/splice
alleles that map to the same protein substitution require an explicit identity
record; do not accidentally merge or split them to change donor counts.

Retain missense donors with the same endpoint and compatible mechanism.
Truncating, start-loss, splice, synonymous and other nonlocal mechanisms still
receive empirical posteriors but remain separate from this local structural
kernel. Unknown consequences must be reviewed or explicitly excluded rather
than admitted by a compatibility default. Document known gain-/loss-of-function
or inheritance strata before averaging biologically opposed effects.

The full-dataset empirical hyperparameters stay fixed, as requested. This
removes the target's direct donor contribution, while retaining its indirect
influence on the common empirical prior. Label the result **variant-excluded
density with full-dataset empirical hyperparameters**, not wholly independent
leave-one-out refitting of the entire pipeline. An optional hyperparameter-LOO
sensitivity must retain other same-residue substitutions too.

For biological-unit copies, calculate the spatial context for each homologous
target position with chain identities preserved. Collapse equivalent copies
of the same donor variant to one contribution per target context (primary:
maximum kernel weight, equivalently nearest eligible structural copy), then
average supported equivalent target contexts and report their spread and
supported/total context counts. Do not replace an unsupported context with
zero or imply full assembly coverage; all contexts unsupported means missing
density. This incorporates
interfaces without multiplying clinical observations by oligomer copy count.
Record this copy-collapse rule as a new explicit implementation choice.

## 3. Historical compact kernel and sensitivity

The recovered [original function](https://github.com/kroncke-lab/Bayes_BrS1_Penetrance/blob/master/func_dist_seq.R)
uses a compact sine transition, with midpoint a:

```text
K(d; a) = 1                           if d <= a - pi
          0.5 - 0.5*sin((d-a)/2)      if a-pi < d <= a+pi
          0                           if d > a + pi
```

Use a=3 Å as primary: K(3)=0.5, K(0) is approximately 0.999, and support ends at
approximately 6.14 Å. The original code uses 3.14 for pi; mathematical pi is the
documented numerical cleanup. Density normalization makes the negligible
zero-distance scale difference unimportant. Prespecified half-distance
sensitivity: a in {2, 3, 5} Å; report K(0) as well because a<pi does not have a
full unit-weight plateau. Also test the historical sigmoid as a labeled shape
comparison, normalized/tuned to the same relative half-distance.

The original source describes residue-centroid distances. Current PPA density
uses C-alpha coordinates, so add an explicit coordinate/metric contract.
The old centroid-producing script's precise atom selection has not been
recovered. Use PPA's supported side-chain heavy-atom mass-weighted center
(glycine C-alpha fallback) as an explicit new primary metric, not a claim of
bit-for-bit reproduction of the old distance file. Propagate this choice
through dense, streamed and fragment paths, which currently use C-alpha.
Three angstroms must not silently mean minimum heavy-atom distance in one run
and C-alpha distance in another.

## 4. Disorder and the biological unit

Recover the original polymer rule:

```text
d_polymer(i,j) = 3.8 * sqrt(abs(canonical_residue_i - canonical_residue_j)) Å
```

Use canonical sequence separation, not row index or author numbering. Primary
polymer pairs must both be in the same contiguous, same-chain disordered
segment. Mixed ordered/disordered pairs have unavailable geometry rather than
invented coordinates. Do not extend the polymer across a folded domain, between
different disordered segments, or across chains. A separately labeled boundary
model can be tested later; it is not an unrecorded fallback.
Test the existing PPA
5.5*N^0.55 parameterization separately. Do not use arbitrary predicted IDR
coordinates to create long-range contacts, even when only the donor is
disordered. PPA's current target-only routing needs repair.

No cross-chain polymer distance exists merely from matching residue numbers.
Missing experimental coordinates may represent disorder, unresolved folded
regions, or a truncated construct; distinguish those cases using sequence,
construct mapping and independent disorder/structure-confidence evidence.
Materialize missing canonical residues as flagged rows. High-confidence but
unresolved folded residues remain structurally unavailable until mapped to
an appropriate structure, rather than being relabeled disordered.

For predicted coordinates, begin with pLDDT >=70 eligible for the 3D path;
pLDDT <50 is candidate disorder requiring a recorded contiguous-segment call;
50–70 stays flagged ambiguous and contributes no primary 3D edge. Experimental
B-factors are not pLDDT. Preserve independent disorder annotations where
available and record the evidence for every segment assignment. More permissive
predicted-coordinate eligibility is a labeled sensitivity, not a silent change.

The requested short kernel has a concrete consequence: with a=3 Å and
3.8*sqrt(N), polymer donors with N>=3 lie beyond the cutoff. Apart from another
variant at the same residue, only positions one or two residues away can
contribute. Keep that primary setting and explicitly report sparse/no-donor
regions; test wider distance settings rather than inventing contacts. With no
eligible donors, density is missing and the downstream model can fall back to
the gene empirical prior with a missing-density flag.

### First biological unit and expansion order

Start with **GCK**, whose biological unit is a monomer. Recover a correctly
numbered experimental monomer (candidates 1V4S / 1V4T) mapped to canonical
P35557-1, compare conformational states, and use the full-length AlphaFold
monomer for confidence/missing-region context. Confirm isoform and WT residues
against every donor before calculating distances. Both experimental candidates
are annotated as isoform 2, so their author numbering is not assumed canonical.
Verify with [UniProt P35557](https://www.uniprot.org/uniprotkb/P35557/entry) and
[1V4S](https://www.rcsb.org/structure/1V4S) / [1V4T](https://www.rcsb.org/structure/1V4T).
A biological monomer is a
complete unit for this pilot, not an oligomeric protein analyzed as a shortcut.

Then add LDLR after domain/linker geometry review, and KCNQ1 once a suitable
numbered functional assembly is available. The local KCNQ1 3BJ4 tetramer has
only about 5.5% expected canonical coverage and is insufficient for the full
gene. HNF1A needs its relevant dimeric context. BRCA2's fragments and extensive
disorder require separate frame/support handling; no cross-fragment distances
are invented. The availability gate is per structural region and does not
discard all empirical posteriors for a partially mapped protein.

## Implementation sequence and acceptance checks

1. **Freeze count and identity inputs.** Save the full count table, transcript,
   allele/protein keys, outcome definition, source hashes and prior-eligibility
   flags. Confirm count partition, positive denominators, uniqueness, and exact
   agreement with retained historical counts. Save restored synonymous rows.
2. **Implement empirical moment/posterior module.** Reproduce the formulas and
   preflight table; save alpha/beta, strengths, posterior means/variances and
   intervals for every variant. Test known Beta moment roundtrips, endpoint
   observations, invalid-moment failures and no AlphaMissense dependency.
3. **Acquire and validate the GCK unit.** Save assembly IDs, accession/isoform,
   residue numbering/insertion-code/chain map, WT checks, coordinate metric and
   structure-state provenance. Compare the selected experimental states before
   interpreting structural differences.
4. **Extend PPA's variant-specific density API.** Preserve target/donor IDs
   through preprocessing and chain expansion; add target-only variant exclusion,
   the compact kernel, equal-variant donor weights, copy collapse, pairwise
   disorder routing and missing canonical rows. Keep the current historical
   PPA output path reproducible under a separately named legacy mode.
5. **Validate geometry and exclusion.** Tests must show R100W is excluded on
   every equivalent chain while R100Q is retained; duplicate observations/copy
   count do not add carriers; K(3)=0.5 and the cutoff is exact; numbering shifts
   and WT mismatches fail; N uses canonical separation; IDR donor coordinates
   cannot create tertiary contacts; unrelated chains/fragments do not gain
   polymer contacts; zero-donor rows remain missing with an explicit fallback.
6. **Build the pilot map and sensitivities.** Save variant LOO densities,
   empirical posteriors, donor lists, contributing variant/residue counts,
   weight sums, largest donor share, structural/polymer source, and kernel
   settings. Draw from donor Beta posteriors to propagate conditional donor
   uncertainty, reusing each donor draw across all structural copies; label
   the intervals conditional on fixed empirical hyperparameters and geometry.
7. **Compare predictive features on the same variants.** Compare empirical
   baseline, density, AlphaMissense and their combination; include a simple
   sequence-neighbor control. Preserve variant-only LOO, including same-residue
   neighbors. Separate an initial fixed-setting comparison from bandwidth
   selection; tuning and regression fits must not use the held-out variant's
   target value. In an outer variant-LOO regression for target i, remove i from
   every training row's donor features as well as from D_i; each training
   row's own variant also remains excluded from its local feature. Retain
   other substitutions at those residues and keep the requested shared
   empirical hyperparameters fixed. Report MAE/MSE and probabilistic/count-based scores with
   support strata, and label the shared empirical-hyperparameter dependence.
8. **Expand after the pilot's artifact and geometry checks pass.** Use the
   same versioned count/prior/kernel/identity contracts for the next genes;
   do not tune toward a desired histogram shape or publish an improvement
   before the comparison has been performed.

Expected work: a bounded GCK pilot first, followed by validation and biological
unit completion across genes. The distance engine already exists; the largest
new work items are variant-specific exclusion, empirical-prior integration,
pairwise disorder handling, and validated assembly/centroid mapping.

Reproduce the preflight from this folder with:

```bash
python full_empirical_universe.py
python kernel_checks.py
```

Both scripts use committed local snapshots/constants and submit no model or
structure jobs. `full_universe_source_counts.csv.gz` plus
`source_count_provenance.json` preserves the full eligible count input;
`full_universe_moments.csv` and `full_universe_empirical_posteriors.csv.gz`
contain the preflight estimates. `kernel_sensitivity.csv` records the kernel
anchors and polymer support. Companion method-recovery reports distinguish
historical source from chosen adaptations. Numerical checks passed and the
empirical outputs reproduced byte for byte from the local snapshot.

Grok consultation and the disposition of its suggestions are retained under
`reviews/`; the numerical preflight and historical source references accompany
this plan. This document is the detailed execution specification; `TASKS.md`
remains the active forward checklist.
