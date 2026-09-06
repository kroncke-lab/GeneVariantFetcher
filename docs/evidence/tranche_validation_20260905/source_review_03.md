# Second-tranche source and influence audit

This post-lock audit explains the measured result. It changes no gold value,
locked prediction, threshold or runtime. Paths and hashes are retained in
`source_review_03_receipt.json`; exact paper count changes are in
`paper_count_changes_03.json`.

## SCN5A 20129283: reference agreement is not yet person-count adjudication

Table 2's H558R row contains Number = 408. Both protocols retain identical
row/header evidence from the deterministic table parser. Baseline carrier and
unaffected fields are NULL; candidate supplies 408/408, matching gold. Affected
remains NULL. The primary CLEANED input hashes are identical, and the reviewed
FULL_CONTEXT files also match. This was a previously-unscored paper.

The full context adds essential semantics omitted from the narrow row excerpt:
Table 2's caption describes 2,600 reference alleles (line 118), while methods
(line 28) describe 1,300 healthy volunteers and credit prior reports. The header
says only Number (line 120); the H558R row is line 136. We have not established
that 408 means distinct people rather than allele observations, or resolved how
previously reported controls should be attributed under the extraction contract.
Do not interpret gold agreement as independent clinical validation. Retain the
original gold and flag unit/attribution adjudication separately.

The whole paper has unchanged identity scores: 339 TP / 7 FP / 78 FN. Carrier
error falls 592 → 184 and unaffected error 408 → 0; affected error stays 436.
The H558R row alone explains both 408-unit reductions. The existing structural
outlier exception is a plausible retention mechanism, but an identical-response
replay was not performed here, so this audit does not isolate a causal code
change from fresh-run variation.

The row explains 87.4% of tranche-03 A/U improvement and 79.8% pooled. The
post-hoc count-row exclusion yields 4.01% and 3.43% residual reductions,
respectively; the previously-unscored pooled residual is 5.90%. Carrier gains
are dominated by this paper and tranche-02's SCN5A 27566755. Removing both entire
papers leaves only a 0.37% pooled carrier error reduction. See
`outlier_sensitivity.py` and its JSON; these are influence checks, not new gates.

Other descriptive A/U movements include SCN5A 20102920 affected error 35 → 0,
KCNQ1 18713323 affected 73 → 53 despite supply falling 5 → 3, and SCN5A
28294644 affected 13 → 3. Regressions include SCN5A 29506689 affected 1 → 12
with supply 1 → 0, and KCNQ1 24372464 unaffected 15 → 22 with supply 2 → 0.
These are count-output comparisons, not separately adjudicated source truths.

## RYR2 31970460: one actual-input difference, no score difference

Initial source snapshots match on every attempt. Actual primary text hashes
match for all 120 attempts in tranche 02 and 119/120 in tranche 03. The one
exception is this paper: baseline uses 4,824-byte abstract JSON and candidate
uses 74,456-byte CLEANED body text. Historical validation rejected the staged
body as an abstract/reference shell at 18:57:26; candidate accepted it. The
baseline driver confirms abstract-only expert extraction, while the candidate
uses the staged body. Both expert calls returned zero variants; both final
paper projections score 0 TP / 1 FP / 1 FN with all counts omitted.

This is protocol source-selection behavior on identical initial availability,
not additional acquisition, gold-conditioned source replacement or cross-arm
corpus contamination. Reporting 240/240 identical actual text inputs would be
incorrect. The primary text/representation hashes and both driver-log hashes
are in the receipt. The prospective fixed-availability design permits normal
protocol validation; the result must retain this exception.

## Operational retry sensitivity

The baseline BRCA1 metadata type failure was retried once under the unchanged
runtime before gold access. The first complete retry has no TP or FP and no
counts, so the predeclared empty-failure sensitivity is numerically identical
to the main comparison. Preserve the failed attempt and its $0.09305 cost;
do not present the arm as uninterrupted. No production metadata fix was made
between arms.
