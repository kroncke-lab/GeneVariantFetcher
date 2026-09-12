# Grok consultation: disposition

Subsequent user clarification supersedes the compact-sine choice discussed in
these reviews: the active plan uses the reviewed normalized sigmoid
`2/(1+exp(log(3)*d/3))`, with a positive tail and no hard 20 Å or sequence cutoff.
Earlier compact-sine recommendations remain as review history. Empirical
posterior donors, variant-only exclusion and the polymer distance rule remain
the same.

The user explicitly requested Grok's help with all four methodological steps.
Two tool-disabled, verbatim reviews were requested using the local Grok CLI:
the first discussed the user's instructions and generic math; the second added
the recovered public historical formulas. No internal per-variant data or
private source code was sent. Prompts, complete returned text, raw results and
usage are retained here.

Both reviews completed on `grok-4.6-build`; combined reported cost was $0.06636.
The extracted Markdown has repository-standard trailing whitespace cleanup;
the raw JSON responses preserve the exact returned text.

The execution plan is controlled by the user's instructions and recovered
method, not by an external reviewer's preferred replacement model.

## Adopted

- Variant-key exclusion across every structural copy, with other substitutions
  at the same residue retained.
- One full gene–disease empirical prior, then count-updated Beta posteriors as
  donors; AlphaMissense is not part of this empirical step.
- MSE for moment inversion, MAE as a separate diagnostic, explicit positivity
  checks for Beta parameters, and no silent invalid-moment repair.
- Equal-variant spatial weighting, consistent with the original function.
- Explicit coordinate metric and kernel version, short-distance sensitivity,
  same-chain/same-segment polymer routing, and donor as well as target quality
  checks.
- Collapse copies of a donor variant before aggregation; do not interpret
  structural multiplicity as extra independent carrier observations.
- Save support and missing-density flags; propagate donor uncertainty while
  labeling fixed hyperparameters and fixed geometry.
- Label the primary run as direct variant exclusion with shared full-dataset
  empirical hyperparameters. Do not call it independent external validation.

## Primary choices retained instead of initial Grok recommendations

- Grok initially preferred unweighted moments with a plug-in binomial-noise
  subtraction. After recovering the historical source, the plan retains
  `w=1-1/(n+.01)` and `sum(w*residual^2)/M`. Normalized weighted MSE is a
  sensitivity. The simple proposed noise correction is especially unreliable
  for singleton and boundary observations, so it is not the primary method.
- Grok initially preferred minimum heavy-atom distance and a rescaled cosine.
  The plan instead starts from the recovered sine with midpoint 3 Angstrom
  and a documented side-chain center-of-mass metric. The old centroid producer
  was not recovered, so this metric is an explicit choice rather than a claim
  of exact historical geometry reproduction.
- Grok initially suggested a polymer prefactor near 5.5 Angstrom. The user’s
  original `3.8*sqrt(N)` is primary; existing PPA `5.5*N^0.55` is a sensitivity.
- No leave-whole-residue-out requirement is introduced. A hyperparameter-LOO
  comparison, if performed, is optional and still retains other substitutions
  at the residue.

## Mathematical corrections to the initial review

- Shared hyperparameter influence is not guaranteed to be small or O(1/M) in
  every finite dataset. Measure it rather than assuming a bound.
- Side-chain centroid distance need not always be smaller than C-alpha
  distance. There is no universal strict ordering of those two metrics.
- A Beta distribution with concentration 2 is uniform only when its mean is
  0.5, so it is not a generic fallback preserving an arbitrary mean.
- For fixed normalized weights and independent Beta donors, the variance of
  their linear weighted average is exactly the weighted sum of donor
  variances. Hyperparameter, shared-cohort and geometry uncertainties remain
  additional sources; a Gaussian interval is not automatically appropriate
  near the boundaries. Use bounded Monte Carlo quantiles for the conditional
  density interval.
- A sum of kernel weights is support, not a calibrated effective number of
  clinical carriers. Report distinct donors, weight concentration and the
  Kish-style donor count separately, without claiming independent outcomes.

The historical-method follow-up checks these corrections and the retained
primary definitions. The final plan records unresolved experimental geometry
as execution work; it does not block starting the empirical preflight.

## Follow-up review: additional corrections

The follow-up accepts the recovered primary formulas, variant-only exclusion,
equal-variant posterior donors, compact sine and original polymer scale. Its
remaining imprecisions do not enter the plan:

- Its opening phrase about gnomAD is ambiguous; the actual plan and both
  prompts explicitly treat gnomAD counts as unaffected.
- Its requested mandatory concentration floor is not adopted. All preflight
  moment fits are valid; future invalid moments fail explicitly. No arbitrary
  floor or fallback is inserted into the requested empirical fit.
- At v=mu(1-mu), kappa=0, not negative. The derivative of moment concentration
  with respect to v is -mu(1-mu)/v^2; the singularity is at v=0, not at the
  Bernoulli upper bound. The plan uses the exact valid-domain check.
- Historical normalization and binomial sampling noise exert competing effects
  on inferred concentration, so the direction of total latent-variance bias
  cannot be guaranteed from the noise term alone.
- Both structural and polymer donors are weighted by their distance kernel;
  "equal variant weights" means no additional n or precision multiplier.
  The follow-up's instruction to average polymer donors equally without K is
  not followed.
- A deterministic average of target-copy contexts is not a mixture distribution.
  Shared donor draws preserve the induced dependence; between-context spread
  is reported separately as structural context sensitivity.
- The effective number of donor variants is not simply sum(w); support mass
  and a Kish-style concentration measure are reported separately.

Independent local review also added exclusion of the outer held-out variant
from every training density feature and explicit reporting of unsupported
target-copy contexts. These changes retain the user's variant-only rule and
the fixed full-dataset empirical prior.
