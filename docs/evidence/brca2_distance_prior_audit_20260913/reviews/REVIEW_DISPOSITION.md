# Fresh CLI mathematical reviews — 2026-09-13

Both new consultations completed. These were **conceptual reviews of a synthetic brief**, not reviews of the real BRCA2 records, source audit, measured distance weights, prior fit or implementation. The exact identical as-sent text is in `shared_prompt.txt`, `agy_prompt.txt` and `grok_prompt.txt`. Raw responses, readable reviews, binary hashes, timing and usage receipts are preserved alongside this document.

| CLI | Requested routing | Reported model | Elapsed | Total tokens | Reported cost |
|---|---|---|---:|---:|---:|
| agy | gemini-3.1-pro-high, high effort | Not echoed | 50.920 s | 20,037 | Not provided |
| grok | grok-4.6, high effort | grok-4.6-build | 76.596 s | 17,588 | $0.01709384 |

The successful headless invocation pattern and unchanged binary identities were reused; these are new responses, not recycled earlier reviews. No source-access tools were authorized to either reviewer.

## Accepted mathematical conclusions

- A normalized average is bounded by its donor posteriors. If every donor has the same sparse-count posterior, distance reweighting cannot change that value. All-zero-affected eligible segments can isolate the shared-prior contribution without requiring remote affected donors. Different unaffected counts can still make the weights matter.
- The synthetic Beta(1,9) example is correct: each A=0/U=1 posterior is 1/11; their normalized average remains 9.09%; pooling 100 unaffected observations once gives 1/110, or 0.909%. An own posterior at most 0.1% needs U≥990 under that synthetic prior. These are illustrations, not BRCA2 estimates.
- Repeated per-variant regularization, averaging raw variant fractions and pooling kernel-weighted counts once answer different questions. Raw-kernel scaling changes the pooled estimator's effective exposure; it cancels from the normalized posterior average. Neither operation automatically produces calibrated disease risk.
- Tiny absolute weights can still normalize to full influence in a sparse context. Report raw kernel mass, nearest donor, normalized weight radii, effective donor count, affected/unaffected evidence and context disagreement. A steep raw tail alone does not settle normalized influence.
- Count ownership, clinical endpoints and source type must be audited before retuning a kernel or prior. Catalogue membership and pathogenicity classes cannot substitute for observed affected carrier counts. This is a general principle from the reviews; they were not given the actual source finding.

## Qualifications and reviewer overstatements

Agy's phrase “purely an artifact” is too strong: finite-count Bayesian shrinkage is a defined model operation, and these mathematics alone do not establish its clinical validity or invalidity. The current score is a posterior-derived feature; calling any of the three quantities “regional risk” requires a separate observation model and validation. A deterministic weighted average itself does not require independent donor variants; dependence matters for uncertainty and evidence ownership.

Raw A/n is defined at n=1 and is a valid descriptive fraction; it is not automatically a reliable individual disease-risk estimate. Averaging context scores does not itself add duplicate counts. Treating repeated context estimates as independent evidence, or summing repeated carrier contributions when pooling, can overstate support. Deduplicate repeated observations of the same person while preserving genuinely distinct relatives; family dependence is a separate issue.

Both reviews suggest support diagnostics and possible alternatives. Abstention, another baseline, count pooling or a changed prior remain explicitly labeled proposals, not replacements for the user-requested primary. Validation must use independently held-out observations/endpoints and preserve person/family/source ownership, rather than merely reproducing the same posterior labels used to build the feature. A narrower kernel or changed estimator must not be chosen to force a 0.1% target.

## Local verification versus external opinion

`local_math_checks.json` records independently computed arithmetic for the actual fixed prior and a hash of the inspected kernel implementation, plus the synthetic arithmetic checks. Those real-data parameter values and implementation receipt were not sent to either CLI. The separate local distance and source audits establish any claims about actual BRCA2 evidence; this review does not independently corroborate those observations.

## Approval scope

Automatic approval review rejected transmitting the detailed internal BRCA2 brief because the user's CLI authorization did not specifically cover that sensitive payload. Neither CLI executed that request. `unsent_detailed_prompt.txt` preserves the rejected draft; `approval_scope.json` records the exact rejection and action. The approved safer alternative replaced it with general mathematics and synthetic examples, excluding actual internal priors, measurements, counts, identifiers and source-provenance findings. The detailed payload was not retried through another channel.
