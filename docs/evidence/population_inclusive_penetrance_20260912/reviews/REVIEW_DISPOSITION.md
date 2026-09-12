# Updated Grok consultation and local disposition

The completed conceptual consultation used the installed Grok CLI, requested
`grok-4.6` with low reasoning for bounded turnaround, and returned
`grok-4.6-build`. The exact corrected payload is `grok_prompt.txt`, response is
`grok_corrected_raw.json`, and readable response is `grok_review.md`. The
completed call reported $0.01232976. Grok had no tools, source access or
independent data access; the payload contained the abstract mathematical
procedure. Local source and saved outputs were inspected separately.

An earlier xhigh call used a mistaken prompt description of inverse-uncertainty
Gaussian weighting. Source inspection established distance-only sigmoid
weighting, so that call was interrupted (exit 130) and superseded. Its exact
payload is preserved as `grok_initial_prompt.txt`; it produced no usable
review or final usage accounting. No inference from that mistaken description
was adopted. The completed corrected call explicitly describes the actual
distance-only sigmoid kernel with fixed weights.

## Accepted findings

1. The historical moment formula and posterior orientation are correctly
   implemented: `alpha + affected`, `beta + literature_unaffected +
   population_carriers`. Repeating the unaffected-carrier objection would
   contradict the user's fixed assumption and was outside review scope.
2. Full-span priors target the mixture of observed gene-region units. They
   are not missense-specific priors. GCK full-span mean is 0.01989008,
   prior-footprint mean 0.11169618 and coding/splice mean 0.23741848. Keeping
   these sensitivities visible is necessary to interpret the lower full-span
   prior and downstream missense density.
3. Variant-only exclusion removes the target donor globally from training
   density features, but shared gene hyperparameters remain fixed. This is
   a conditional analysis. Distinct same-position variants remain eligible
   by design; this answers a variant-level question, not a residue-level
   generalization question.
4. The regression target is the empirical posterior mean, so the reported
   prediction metrics measure recovery of that statistic. They do not by
   themselves establish outcome calibration. Beta-draw density intervals are
   conditional on fixed counts, identities, hyperparameters and geometry.
5. A clinical protein aggregate can contain several genomic alleles, while
   unmatched population units retain genomic-allele grain. The supplied
   protein-aggregate sensitivity measures this mixed-grain effect; exact
   allele membership and quarantined clinical counts remain important.

## Reviewer errors and qualifications

- **Reject Grok's variance-direction claim.** Because `sum(w) < M`, dividing
  the same weighted sum of squared deviations by `M` produces a *smaller*
  variance than dividing by `sum(w)`. The historical prior is therefore
  *more concentrated*. In saved GCK results, historical variance 0.00616805
  gives strength 2.16056, whereas normalized variance 0.01652023 gives
  strength 0.18004. The historical implementation is intentional and
  correctly reproduced; Grok's opposite algebraic interpretation is wrong.
- Same-position biological correlation is not a target-label leak merely
  because correlated variants are retained. The question expressly keeps
  other variants at the same residue. No new residue-holdout requirement or
  pipeline change is introduced by this review.
- The wording about clinical counts being "rolled onto" genomic alleles
  should not imply duplicated numerators. The code instead creates one
  clinical aggregate unit and records its population allele members once.
- Suggested new exclusions, priors or thresholds are not acceptance gates.
  The request preserves the historical estimator and named sensitivities.

## Independent bounded computational check

`loo_hyperparameter_spotcheck.json` contains an exact leave-one-unit-out
moment calculation using sufficient statistics from the 14,553 GCK full-locus
units. Across all 634 eligible missense units, dropping one unit changes the
shared prior mean by at most **0.00016423** (0.01642 percentage points), alpha
by at most 0.00031700 and beta by at most 0.00553678. This bounds the direct
hyperparameter change for this dataset; it does not rerun regressions or
establish independent outcome validation.

Sources inspected locally: `rebuild_union.py`, `analysis/union_coverage.csv`,
`analysis/empirical_prior_comparison.csv`, GCK empirical posterior shards,
`run_structure.py`, `structure_statistics.py`, `structural/run_checks.json`,
and the sibling ProteinProximityAnalysis module
`src/alphafold_rin/empirical_density.py`. No additional blocking computational
error was demonstrated in this bounded review. The scope, conditional
uncertainty and mixed-grain qualifications remain material.
