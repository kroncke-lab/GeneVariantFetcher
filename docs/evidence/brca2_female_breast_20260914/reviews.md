# Method consultation and disposition

2026-09-14. Agy, Grok and Claude all returned successfully from their CLIs.
They reviewed the same generic mathematical brief and synthetic examples;
they did **not** inspect the clinical records, implementation or resulting plot.
Primary-source adjudication and executable checks remain the evidence for this
run. Raw local receipts live in the ignored
`results/brca2_female_breast_20260914/reviews/` directory. Their hashes are in
this evidence folder's manifest; only final critiques are summarized here.

| CLI | Model reported/configured | Exit | Main useful finding | Disposition |
| --- | --- | --- | --- | --- |
| Agy | `gemini-3.1-pro-high`, high effort | 0 | Track A and U attrition independently; use XX-specific homozygotes; exclude n=0 before moments | Implemented in observation, membership and exclusion ledgers and validation |
| Grok | `grok-4.6` (reported `grok-4.6-build`), high effort | 0 | Sex, endpoint, germline and identity must agree before adding counts; preserve variant types and frozen geometry | Source-reviewed subset, separate missense/nonsense unions, exact XX counts, shared density engine |
| Claude | CLI configured default (reported `claude-opus-5[1m]`) | 0 | Filtering can remove controls disproportionately; moment variance and binomial noise are unresolved; avoid a second pipeline | Separate source-corrected all-sex comparison, explicit subset limitations, existing prior/density/plot functions reused |

All three warned that population carriers treated as unaffected do not establish
absence of future disease. This run follows the user's explicit population-U
assumption and labels the result descriptive. Clinical XX/female equivalence is
not inferred: XX is the population proxy; clinical sex comes from primary-source
cohort or patient evidence. Unknown clinical sex, another cancer endpoint, and
tumor-only observations are excluded, never converted to U.

The historical empirical prior is retained to isolate this source/universe
correction. Its weight `1 - 1/(n+0.01)` nearly discards singletons, and its
M-normalized weighted variance is not a beta-binomial marginal-likelihood fit.
Changing only its variance normalization would not remove binomial sampling
noise. The next estimator comparison remains in [TASKS.md](../../../TASKS.md),
with hyperparameters refitted inside validation folds before calibration claims.

Suggestions not adopted:

- Presuming clinical sex, coding uncertain outcomes as A or U, or halving
  population totals would defeat the source restriction. Explicit all-women
  cohort evidence is accepted without requiring a sex tag on each row.
- A hard polymer cutoff and leave-residue-out validation contradict the user's
  specified positive tail and variant-only LOO. Other variants at the same
  residue remain eligible. Whole-alias/copy self-exclusion is checked.
- Both gold and blue use the **same normalized W**. The review brief's notation
  prompted a concern that blue was unnormalized; executable reconstruction
  checks this directly. Raw-kernel count pooling is a separate third feature.
- A>0 with U=0 is not inherently a failed population lookup: a clinical allele
  may be absent from the complete observed XX inventory. Such sparse cases are
  reported rather than forced to low risk. Likewise, missense mean below
  nonsense mean is a diagnostic expectation, not an arithmetic acceptance rule.
- Scores need not remain identical when a target's own counts are unchanged:
  changed donors and a refitted full-data prior can change its neighborhood.
  Unchanged geometry/code hashes, weight identities and count conservation are
  the applicable checks. Unsupported positions retain missing scores.
- Some review arithmetic was imprecise: at n=0 the weight is -99, whereas A/n
  is undefined; Grok's displayed effective-sample-size formula was malformed.
  Actual computations, rather than those expressions, determine acceptance.

Initial sandbox CLI attempts failed on local socket/auth/network access. The
authorized retries above completed; no GUI fallback was needed.
