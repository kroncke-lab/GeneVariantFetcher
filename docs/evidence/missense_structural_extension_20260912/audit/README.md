# Independent runner audit

The HNF1A and KCNQ1 runner outputs pass independent reconstruction. This audit
does **not** call the PPA density engine: it derives distances and weights directly
from the saved canonical coordinates, using the agreed sigmoid and polymer
rules. It then reconstructs selected outer-LOO regression inputs. The existing
fractional-logistic optimizer is reused after those independent feature calculations.

| Gene | Missense units | Cartesian only | Polymer only | Mixed contexts | Unavailable |
| --- | --- | --- | --- | --- | --- |
| HNF1A | 785 | 183 | 536 | 7 | 59 |
| KCNQ1 | 1,061 | 586 | 389 | 0 | 86 |

Every saved positive donor weight and primary density agrees with the independent
calculation to numerical precision. Donor copies are collapsed to the closest
eligible copy within each target context, each supported context is normalized
separately, and the contexts receive equal weight. This prevents tetramer copies
from becoming four observations. All target self weights are zero, while other
variant units at the same residue remain available. Unsupported targets retain
NaN densities, rather than zero or a substituted prior. Genomic members are
disjoint across the frozen units, and all donors use the same fixed gene-specific
missense prior and canonical WT reference.

The audit also checks all saved unfitted sequence-control predictions and 45
selected fitted predictions across structured, polymer, and mixed target examples.
For each held-out target, its donor identity is independently removed from **all**
training neighborhoods. The structural feature is renormalized inside each
context before contexts are averaged; the sequence control is recomputed after
the same global exclusion. Training labels, AlphaMissense features, and target
indices remain correctly aligned. This verifies the primary merge order and the
outer-LOO handoff in the runner, separately from the engine's own tests.

HNF1A's seven mixed variants occur at canonical residues **181 and 200**. One
chain context provides experimental Cartesian geometry while the other provides
candidate-IDR polymer geometry. Their context-density spreads are approximately
0.056–0.077. The published mixed density is the equal average of those supported
contexts, not a pure experimental 3D result. The full HNF1A supported count is
therefore dominated by polymer-only estimates. Source-stratified scores are
necessary when discussing a benefit from actual 3D proximity.

Two interpretation points remain relevant:

- **Identity grain is inherited from the prior freeze.** Seven HNF1A and eleven
  KCNQ1 protein substitutions each have two distinct population DNA-allele units.
  Distinct DNA alleles encoding the same protein change may support one another.
  Clinical protein aggregates remain single units. This is the documented
  variant-unit LOO experiment, not a leave-one-protein-substitution-out experiment.
- **Intercept-only LOO Spearman is not a ranking benchmark.** Its training mean
  varies inversely with the held-out observation, giving mechanically negative
  correlations (approximately −0.52 for HNF1A and −0.71 for KCNQ1). This should be
  suppressed or explicitly labeled in a user-facing ranking comparison. Its
  held-out error and count scores remain meaningful baseline diagnostics.

On the common AlphaMissense targets with Cartesian-only support, adding density
to AlphaMissense has posterior MSE 0.034954 for HNF1A versus 0.035406 for adding
the sequence control; KCNQ1 gives 0.043503 versus 0.044050. These are small internal
differences in this ascertained count experiment, with no independent clinical
validation or uncertainty claim about the difference. The full-sample priors
remain fixed as requested; that shared-prior contribution is not removed by LOO.

The runner's subsequent LDLR correction is appropriate in principle: if global
removal of the held-out donor exhausts a training target's neighborhood, drop that
training target from **all** model fits in that fold, preserve the held-out
evaluation, and record the lost-support identities/counts. Do not replace the
missing feature with zero. HNF1A and KCNQ1 had no such training-support losses.

Reproduce from the GVF root:

```sh
/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/.venv/bin/python docs/evidence/missense_structural_extension_20260912/audit/independent_runner_audit.py
```

`independent_runner_checks.json` records checks, inspected source hashes, mixed
target details, and all selected prediction differences. Geometry files remain
unchanged by this audit.
