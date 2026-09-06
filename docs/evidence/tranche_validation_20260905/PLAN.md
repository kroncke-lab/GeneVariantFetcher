# Prospective tranche validation, 2026-09-05

User request: test another tranche or two, then update the forecast. API costs
only count against the existing $100 campaign; $4.69165 was used in the targeted
tests. Routine runs use Azure. No Anthropic consultation is planned.

## Design fixed before scoring

- Start with registered continuation tranche 02: 120 gene-paper attempts,
  110 distinct articles. Run the registered frozen baseline before current
  candidate, in consumption order. Use tranche 03 as a second validation if
  the first pair leaves the forecast uncertain, including a null result.
- Baseline follows the existing nine-file `506a949c` protocol definition used
  in continuation tranche 01. The remaining runtime infrastructure is common
  to both arms. The new source-snapshot staging branch is also common; it
  prevents live harvesting and cached prediction reuse. Record every file hash
  and archive the complete runtime before each arm.
- Candidate is the implemented main code at `d2996c74`, with no reader edits
  based on the new tranche's gold. Proposed clinical-table continuation and
  patient-ID joins are not implemented and are not being tested here.
- Baseline preparation freezes source text and local assets. Candidate uses
  the same snapshot. Source changes during preprocessing are allowed only as
  ordinary protocol behavior and are recorded by actual-input rebinding.
  Source bytes are frozen; PubMed bibliographic metadata still loads normally.
  This is a reading comparison, not an isolated test of live acquisition gain.
- Read no new gold values before each arm's production completion and lock.
  Keep failed/missing-source attempts in the registered denominator. No cached
  predictions, prior databases, corpus sync or review publication.
- Use the registry's PMID-cluster bootstrap and unchanged acceptance thresholds.
  Also report carrier, affected and unaffected supply, conditional error, and
  end-to-end error, including omissions. A null comparison is not a passing
  discovery; a second tranche after a failure is further calibration.
- Prior exposure was audited from already-locked selection membership before
  scoring. Tranche 02 has 18 previously scored PMIDs and 97 previously unscored
  attempts; tranche 03 has 20 and 99 respectively. Retain the full registered
  scores and separately report the predeclared previously-unscored subset.
  Neither entire tranche is described as a wholly unseen holdout.

## Budget and forecast decision

The registry estimates $24.97 for tranche 02's pair and $25.37 for tranche 03's
pair; together $50.34, or $62.93 with its 25% headroom. Reserve $65 for the two
pairs within the remaining $95.31. Check trace-derived API cost after every
arm, and do not start a further arm unless it fits with headroom. Provider
model and token totals are recorded; public-price estimates are not invoices.

The previous 72–76% recall forecast was an engineering estimate for the hard
opened continuation-01 cohort after remaining work, not an expected absolute
score on differently composed papers. Re-estimate transferable improvement
from within-tranche baseline/candidate deltas. Separate the known 20031634
acquisition gain from general reading gain, and keep missing-supplement upside
conditional. Do not call a provider-noise-sized change a repeatable improvement.
