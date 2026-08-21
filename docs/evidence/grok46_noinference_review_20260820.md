# Grok 4.6 review — no-inference and publication-gate changes

> **Superseded scope note (2026-08-21):** this review covered the no-inference
> patch then presented to Grok. A later repository-wide adversarial pass found
> additional DOI/PMID token, readiness, active-DB, trace, and publication
> bypasses outside that review packet. Those findings are recorded in
> `grok46_fail_closed_followup_20260821.md`; this document is not an unqualified
> release GO for the current collaborator workflow.

**Date:** 2026-08-20
**Reviewer:** Grok 4.6 at `xhigh`, the highest reasoning effort exposed by the
local Grok CLI (`max` is not a valid CLI effort name).
**Scope:** extraction contract, penetrance/phenotype handling, collaborator
gates, and the actual repository diff. No files were delegated to the reviewer
for editing.

## Proposal review

Grok ranked these as P0 requirements before a paid gold run:

1. Remove prompt permission to calculate penetrance from counts.
2. Stop the aggregator from deriving a percentage.
3. Stop claim verification from manufacturing affected/unaffected values from
   diagnoses, enrollment, or arithmetic completion.

It ranked these as P1:

1. Let collaborator/publication runs require VariantFeatures enrichment plus
   false-positive quarantine.
2. Require a pinned PMID manifest for publication and refuse publication after
   any failed stage.
3. Do not publish until the Variant Browser consumes trusted fields and
   quarantines.

## Diff review

After the changes and tests were in place, Grok reviewed the actual diff in a
read-only session. Its verdict was **GO for the paid cardiac acceptance run**:

- no P0 or P1 blocker remained;
- it found no gold leakage or paper-specific tuning;
- nulling unsupported phenotype partitions was judged an intentional safety
  policy rather than a disguised benchmark optimization;
- it cautioned that the resulting lock must not be described as a new
  penetrance/affected baseline without reporting the abstention-driven coverage
  change.

This is advisory evidence, not a substitute for tests, a blinded benchmark, or
source adjudication. The locked benchmark and collaborator audits are recorded
separately.
