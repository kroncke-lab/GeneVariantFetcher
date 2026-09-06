# Failed-usage contract review disposition

Claude, Grok and Agy reviewed the change before scores. Valid concerns were
addressed before any unknown-usage artifact could lock:

- The exception now joins each failed-call ID, path and SHA to the original
  write-time index, in addition to the normal manifest verification at lock.
  It enumerates same-paper indexed API calls and requires the full set of
  unknown failures. Missing index entries, missing references, duplicate IDs,
  path escape, altered bytes, wrong papers and successful calls are rejected.
- Known usage from successful sibling calls is retained and checked against
  the indexed provider records. Both paper and run totals remain null when
  incomplete; the run must retain the complete known subset. Aggregate tables
  label their known-subset totals when failed calls have unknown usage.
- The binder checks that all final arms lack score reports, verifies an actual
  failed-call budget reservation, preserves baseline usage separately and
  asserts that every primary/comparison scientific prediction is unchanged.
  The receipt distinguishes filesystem verification from operator attestation
  about prior gold exposure. It does not claim investigator blinding.

Several review claims were conditional on missing code and were not reproduced:

- `command_lock` is the only lifecycle caller of `validate_predictions`.
  `score` verifies the immutable lock instead. A direct schema-only validation
  without a trace root intentionally cannot authorize missing native usage.
- The existing manifest builder checks every recorded trace against the
  write-time index. The added explicit join makes the failure exception locally
  reviewable too; a file/ref checksum alone was not the entire lock contract.
- Any failed API stage can make total paper usage unknown, even when a later
  clinical call succeeds. Requiring the known subset prevents that failure
  from erasing successful usage. Restricting the exception to one stage would
  hide other real failures.
- Reasoning helper changes are a separate, already-tested part of this model
  experiment. Responses helpers return the nested shape; older Grok no-ops and
  new-model effort transmission have dedicated tests. Hypothetical flat/None
  helper responses are not current behavior.
- The binder is deliberately for this completed native clinical overlay, not
  legacy/manual inputs missing required keys. Null/malformed index contexts
  now fail validation through the indexed helper rather than crashing lookup.
- The report says incomplete legacy telemetry **may** include traced failures;
  it does not assert every legacy run had a timeout. Zero is never presented as
  a failed call's price. Campaign costs come from the separate budget ledger.

No source, scientific prediction, reference value, count-acceptance rule or
matching arithmetic changed. Failure papers remain in the fixed denominator.
