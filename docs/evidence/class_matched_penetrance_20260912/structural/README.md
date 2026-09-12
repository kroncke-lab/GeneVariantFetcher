# GCK structural replay with a missense-only prior

The completed replay fits the shared empirical prior only within canonical,
WT-matching **missense** variants. Nonsense receives a separate prior upstream
and contributes no structural donors. Synonymous, noncoding, frameshift, and
other variant classes do not inform this missense prior.

The GCK missense prior uses 634 units, including 385 population-only units:

| Quantity | Value |
|---|---:|
| Affected | 670 |
| Literature unaffected | 72 |
| gnomAD unaffected carriers | 4,369 |
| Empirical alpha | 1.081159 |
| Empirical beta | 1.837596 |
| Prior mean | 37.0418% |
| Prior strength | 2.918755 |
| Posterior after one unaffected carrier and no affected carriers | 27.5894% |

Alpha adds affected counts. Beta adds literature-unaffected and all gnomAD
carriers, as requested. The one-unaffected-carrier result follows directly from
`alpha / (alpha + beta + 1)`; it is a finite-evidence posterior under this
missense-class prior, not the observed fraction of zero affected carriers.

## Structural results

The replay preserves the exact earlier 634 missense identities, counts,
AlphaMissense values, biological-monomer geometry, and spatial settings.
**Every saved normalized donor weight is exactly unchanged.** Only the shared
prior and resulting empirical posteriors change.

The primary 1V4S monomer supports **614 targets**, and the median variant-excluded
density is **37.4072%**. The broader-class prior gave 18.2394% with the same
weights. The average density increase from restricting the prior to missense is
18.7891 percentage points. The per-variant comparison is preserved in
`GCK_prior_scope_comparison.csv.gz`.

Twenty scenarios were rerun: the 1V4S and 1V4T biological monomers and the
AlphaFold monomer, COM/C-alpha metrics, half-distances 2/3/5 Å, and the documented
loop/polymer sensitivities. Both experimental structures support 614 targets;
AlphaFold supports 628. The primary remains side-chain heavy-atom mass-weighted
COM, glycine C-alpha fallback, and the positive-tail sigmoid with half weight
at 3 Å. There is no 20 Å cutoff: 278,640 saved donor pairs beyond 20 Å have
positive weights.

Only the target variant unit is excluded. Other variants at the same residue,
including distinct population DNA alleles, remain. Clinical aggregate members
were consolidated upstream. Donors use empirical posterior means and equal
spatial variant weights, with no additional carrier-count multiplier. The
previous pairwise geometry and same-segment IDR polymer rules are unchanged.

## Variant-only outer-LOO comparison

The all-common comparison retains 610 fitting/evaluation targets. The separate
original-242 comparison retains the exact earlier clinical fitting/evaluation
set and its archived AlphaMissense values. Both use the same expanded 634-unit
donor pool. Every held-out target is removed globally from all training donor
features, while each training row also excludes itself. Shared missense prior
hyperparameters remain fixed by design.

| Targets | Model | MAE vs empirical posterior | MSE vs observed fraction | Mean Beta-binomial negative log score |
|---|---|---:|---:|---:|
| All common, 610 | AlphaMissense | 0.124636 | 0.158462 | 0.981904 |
| All common, 610 | AM + density | 0.123400 | 0.157109 | 0.974993 |
| Original 242 | AlphaMissense | 0.111226 | 0.145272 | 1.310909 |
| Original 242 | AM + density | 0.111203 | 0.145228 | 1.311517 |

Adding structure yields modest numerical gains in the expanded comparison.
The original clinical cohort shows essentially no change and mixed score
directions. A sequence-neighborhood control remains competitive. The replay
does not establish a robust incremental predictive benefit from structure.
The underlying clinical dataset still pools GCK disease endpoints/mechanisms;
these are internal diagnostics with shared empirical hyperparameters.

## Reproduction and verification

The adapter imports the frozen population-inclusive runner and supplies new
class-specific inputs and labels. It does not modify or overwrite earlier
source code or evidence. Input provenance includes the new posterior shards,
prior table, adapter, reused runner/statistics sources, original identity/AM
template, and frozen geometry.

The primary donor tables retain 376,382 pairs in chunks of at most 40 targets.
Independent reconstruction recovers all 614 supported densities within 6.8e-15
and conditional variances within 1.0e-16. All donor weight matrices match the
previous run exactly. Three engine-level global exclusion checks and the five
existing LOO/identity/AM tests pass. The 8,192 conditional Beta draws reuse each
donor draw across targets; maximum Monte Carlo mean discrepancy is 2.99 simulation
standard errors. These intervals condition on the selected counts, identity,
prior, geometry, and kernel. Both generated figures were visually inspected,
and every output file remains below 1.2 MB.

```bash
/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/.venv/bin/python docs/evidence/class_matched_penetrance_20260912/run_structure.py
/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/.venv/bin/python docs/evidence/class_matched_penetrance_20260912/validate_structure.py
```

`run_checks.json` records the class-specific run contract and source/output
hashes. `artifact_checks.json` records independent reconstruction and exact
weight preservation. Scenario and donor matrices are sharded; primary outputs,
predictive metrics, and the two PNG figures remain directly inspectable.
