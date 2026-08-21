# Grok 4.6 fail-closed follow-up — 2026-08-21

Grok CLI model `grok-4.6` reviewed the collaborator-readiness, production-eval,
source-identity, and Variant Browser publication surfaces at `xhigh` reasoning
effort. The user requested “4.6 max”; `xhigh` is the highest effort value
accepted by the installed CLI. Tools, web search, and subagents were disabled.

This was an adversarial implementation review, not a scientific adjudication.
The CLI session was interrupted after it had produced actionable findings; it
did not return a clean final release verdict. Accordingly, this evidence must
not be represented as a Grok GO.

## Findings

1. The first Scholar URL guard could accept a target DOI as a prefix of another
   DOI, or a PMID whose digits were concatenated into a longer identifier.
2. The readiness audit could pass empty or unbound extraction payloads, source
   totals that did not match the locked cohort, empty traces, and raw ambiguous
   VariantFeatures identities.
3. Production evaluation and publication helpers could select a database by
   modification time instead of proving the one completed
   `RUN_STATUS.active_db`.
4. A validation-only resume could obscure the stronger canonical extraction
   trace, and a nominal trace could be present without any verified model call.
5. Variant Browser import could stage unmatched or residue-suspect identities
   merely because they were live in the raw GVF database.
6. Missing-PMID overrides and stale staging behavior could make a partial import
   look complete.

## Changes made in response

- DOI and PMID URL identity are matched as exact normalized tokens, with
  regression tests for prefix and digit-concatenation attacks.
- Collaborator readiness now verifies the exact manifest, valid extraction
  payloads, DB/source totals, active database, nonempty write-time-verified
  traces, species scope, trust tiers, and VariantFeatures identity classes.
- Production evaluation and publish helpers require exactly one completed
  active run; modification-time and backup-file selection were removed.
- Canonical extraction traces survive validation-only resumes and are rejected
  unless write-time verified with at least one model call and no errors.
- Variant Browser trusted import publishes matched identities plus only the
  defensible `novel_in_range` and `cdna_only_unmatched` classes. Other ambiguous
  classes remain in the raw audit trail but are held out of staging.
- Publication refuses missing-PMID overrides and records trusted-projection hold
  counts in the import run metadata.

## Verification state

- GVF offline unit suite: 2,164 passed.
- Variant Browser suite: 585 passed.
- Repository-wide Ruff and whitespace checks passed.
- The raw collaborator audit intentionally fails the current BRCA1, BRCA2, and
  BMPR2 databases on 346, 129, and 56 ambiguous live identities respectively.
  That is an honest pre-publication result, not a regression to waive.
- A fresh 119-attempt cardiac production evaluation and a live trusted staging
  import remain required before promotion.
