# Fresh Agy, Grok and Claude reviews — 2026-09-14

All three CLIs completed new, approved consultations on the exact same **generic
methods brief with synthetic examples**. None reviewed the real fitted priors,
gene-level counts, source records, implementation or plots. Actual scientific and
source validation remains local. No app fallback was needed and no approval block
was bypassed.

| CLI | Requested route | Models reported | Elapsed | Reported cost |
|---|---|---|---:|---:|
| Agy | Existing gemini-3.1-pro-high / high | Not echoed | 45.792 s | Not provided |
| Grok | Existing grok-4.6 / high | grok-4.6-build | 48.639 s | $0.01495524 |
| Claude 2.1.259 | Configured default; no model/effort override | claude-opus-5[1m] plus auxiliary claude-haiku-4-5-20251001 | 80.671 s | $0.180874 |

Claude used safe mode, no tools, strict MCP, plan permissions and no session
persistence. All three exited 0 and their prompt bytes match. Agy's 78-byte
stderr log warns that plan mode has no effect while slash-command expansion is
disabled. Grok's 299-byte stderr log reports a Git-discovery warning about the
unsupported `extensions.relativeworktrees` extension. Both returned complete
review responses; these technical warnings did not fail the consultations.
Claude's stderr log is empty. Claude reported zero web requests, zero subagents and no permission
denials. `consultation_summary.json` preserves exact usage, including cache fields
and Claude's auxiliary call. The reported Grok+Claude sum is $0.19582924; Agy cost
is unknown, so this is not a complete total. Raw outputs and exact as-sent prompts
remain alongside the readable reviews. The launcher and CLI hashes are recorded.

## Conclusions useful for the plot audit

The reviews support treating the primary as a target-excluded neighborhood feature
of posterior means, showing the own posterior separately, and retaining the user’s
fixed gene-by-missense prior, population-as-unaffected assumption, positive kernel
tail and same-residue alternatives. Source typing must precede interpretation:
catalogue membership/classification is not a counted patient, and germline,
somatic and opposite clinical endpoints cannot be treated as interchangeable.
Independent functional measurements test molecular effects; independent clinical
outcomes with a specified time horizon are needed for risk calibration.

Sparse-variant shrinkage, unit weighting and kernel support are different effects.
A constant set of donor labels remains constant under every normalized weight
matrix. The synthetic Beta(1,9) singleton value is 1/11; pooling 100 unaffected
observations once gives 1/110. The primary average, raw-fraction average and
one-prior pooled-count diagnostic have different estimands. Keep those labels and
show raw support, nearest donor, weight radii and effective donor number rather
than inferring locality from the kernel formula alone.

## Corrections: reviewer output is not an authority

Independent synthetic arithmetic is saved in `synthetic_review_checks.json`.
Several proposed checks were wrong and must not become acceptance criteria:

| Reviewer claim | Disposition and correct check |
|---|---|
| Agy: fit a prior to all A=100/U=0 units and require alpha much larger than beta | This is a zero-variance/boundary fit, not a valid orientation test. Use nondegenerate asymmetric A/U pairs, then swap A/U and verify alpha/beta swap; separately exercise the declared degenerate-fit handling. |
| Agy: two identical donor labels reveal whether the target was included | They cannot: averaging identical labels masks self-inclusion. Use distinct labels, e.g. 0.1 and 0.9; with only the other unit eligible the scores must be 0.9 and 0.1. |
| Agy: a massive affected count should change which residues receive highest geometric weights | Counts do not set this distance-only kernel. Perturbing counts with fixed geometry/eligibility must leave W unchanged. Test coordinate mapping with known coordinates and canonical amino-acid identity. |
| Agy: internal score–posterior correlation is guaranteed | False. The asymmetric two-unit example gives perfect negative correlation. Internal evaluation still does not establish independent clinical calibration. |
| Agy: high raw mass proves genuine regional consensus; weight radius above 20 Å proves an artifact | Neither follows. High support can still repeat the same prior or biased labels, and 20 Å is not a biological validity threshold. Report support and sensitivity without adding a cutoff. |
| Claude: w(1) is about 0.5 | It is 0.0099009901; w(2) is about 0.5024876. n=0 is invalid for A/n and also gives a negative fitting weight. |
| Grok/Claude: high dispersion can give negative kappa under the valid historical convention; Claude uses a half-variance threshold | With integer n>=1, 0<w<1, y in [0,1] and the stated weighted mean, v <= (sum(w)/M) mu(1-mu). A finite interior fit therefore has positive kappa. kappa<=0 would require v>=mu(1-mu), not half that value. Zero variance, boundary means, invalid inputs and numerical precision remain real checks. |
| Claude: doubling every within-variant count must increase kappa | Not under this convention. Two opposing singleton units have kappa=100; doubling their counts gives kappa≈0.990099. Cloning the entire unit list is a different operation and leaves the moments unchanged. |
| Claude: prior retention is the share of the score supplied by the prior | Let R=sum(W*kappa/(kappa+n)). The direct prior component is mu*R; its share of D is mu*R/D. For all A=0/U=1 with Beta(1,9), R=10/11 but the prior share of D is 1. |
| Claude: d=0 gives normalized W=1 | Only if it is the sole eligible positive-weight donor. Generally d=0 gives raw K=1. Same-residue contributions remain eligible by request. |

Further qualifications: zeros are legitimate masks for unavailable donor pairs;
the unsupported *final score* must stay missing. The polymer formula is a stated
sequence-distance surrogate, not a measured contact or a guaranteed mean over a
physical ensemble. A comparison of the kernel at that surrogate with an ensemble
mean would be a new modeling sensitivity. Opposite endpoints may be mixed or
misassigned; their pooling does not necessarily numerically “cancel.” Functional
assays or structural predictors are not automatically independent of prior source
selection. Catalogue ACMG criteria should be audited if classifications enter an
analysis, but dropping such criteria is not a substitute for validating actual
patient observations.

Do not adopt Claude’s “never by variant” instruction as a replacement for the
requested primary. Fixed-prior variant exclusion remains the descriptive rule;
independent family/cohort validation is a separate evaluation question. Likewise,
no reviewer-defined cutoff, prior clipping, risk target, abstention threshold or
new estimator is installed by this consultation. A properly specified conditional
variance calculation can be investigated, but between-variant spread alone is not
a posterior interval and no uncomputed clinical uncertainty band is justified.

## Corrected local plot checklist

These are proposed checks for the actual outputs, **not claims that those outputs
have passed**:

1. **Count and endpoint ownership:** every A/U contribution has the correct gene,
   allele class, person/cohort ownership and endpoint. Absent or unassessed records
   remain distinct from observed unaffected counts. No catalogue row supplies a
   synthetic patient. Keep the adopted population assumption explicit.
2. **Prior arithmetic:** validate n>=1 and nonnegative counts; reproduce the exact
   historical MSE divisor M. Test nondegenerate alpha/beta reversal and posterior
   monotonicity. Record the actual policy for zero-variance/boundary fits.
3. **Weights and exclusion:** supported rows sum to one, target-identity weight is
   zero, distinct same-residue units remain eligible, and positive rescaling of
   raw K leaves the primary unchanged. Perturbing a target’s counts with its prior
   held fixed must not change its own target-excluded score.
4. **Support and missingness:** show nearest donor, raw kernel mass, r50/r90,
   effective donor number and population/count evidence where available. Keep
   weak but nonzero support distinct from no eligible support. All-absent support
   yields a missing final score, never a zero disease probability.
5. **Decomposition:** reconstruct D=mu*R+sum(W*A/(kappa+n)); label R, the prior
   component and component/D distinctly. A constant donor-label control remains
   constant; broad shrinkage alone cannot establish that the kernel is too wide.
6. **Geometry:** canonical amino-acid identities and residue indices round-trip;
   coordinate metric and chain/frame identity are explicit; no unaligned frame
   distances or cross-IDR-segment polymer paths occur. Repeated structural copies
   do not multiply carrier exposure.
7. **Residue aggregation:** preserve the variant-unit table behind each residue,
   declare the displayed summary and weighting, and expose multiplicity. Different
   substitutions at one position can differ in mechanism and own posterior.
8. **Rendering:** verify plotted values against the saved table, distinguish own
   posterior/primary/diagnostics, mark geometry sources, and break unsupported gaps.
   Do not label between-variant spread as a confidence or credible interval.
9. **Interpretation:** describe an independently uncalibrated feature, not the
   probability that every substitution causes disease. Any later bandwidth or
   estimator choice needs endpoint-matched evaluation and an honest held-out unit.
