# Pre-run clinical-code review disposition

Claude, Grok and Agy reviewed the unrun script with tools disabled. No paper
scores or gold values were in the review packet. Raw receipts are adjacent.

Implemented: verify the baseline prediction/selection lock and the frozen source
hashes; preserve existing rationales; append accepted evidence only after all
acceptance checks so reverted fills leave no orphan quote; attach unvalidated
quotes to the raw diagnostic lane; represent failed-call usage as unknown rather
than zero; reject a native lock without the required exact telemetry. The
existing validator checks individual fields, so the overlay additionally holds
all new fills on a row when its person totals conflict and cannot be resolved
without choosing which value is wrong. That is an explicit conservative hold,
not a claim that every held quote is false. Six synthetic tests cover fills,
contradictions, ungrounded/derived values, conflicting proposals, evidence
rollback and a nonzero original row index.

Claims checked against code and not adopted:

- `validate_paper_response` does not emit "phantom omissions" for fields absent
  from a single-field candidate: `raw is None` continues without rejection.
- `variant_id` is an opaque identifier preserved in `RecoveredCount`, not an
  index into the validator's one-element target list. The new nonzero-index
  test proves the intended projection behavior.
- The validator does not consult `paper_derived` to relax source binding.
- Request caps and full returned usage, including Astra's reasoning counter,
  are captured by the common trace wrapper. They are not absent from traces.
- `prepare` uses gold only for PMID eligibility; its answer values and variant
  identities are not passed to the reader or overlay. The setup wording now
  states this boundary directly.
- The raw **scoring lane** contains only unique, nonconflicting proposed fills
  on retained identities. Full proposals, including conflicts and unmatched
  identities, remain in the response trace and audit; an ambiguous field is
  not forced into a scalar score by selecting a winner.
- Intermediate progress files are permitted; scoring waits for the immutable
  standard lock. Baseline production traces are bound as well as new native
  reader/decision traces.

Three concurrent clinical calls are permitted, each capped at 32,000 output
tokens with a 1,200-second request timeout. No retry is dispatched automatically
by this script; a timeout retains the campaign's uncertainty reservation.
