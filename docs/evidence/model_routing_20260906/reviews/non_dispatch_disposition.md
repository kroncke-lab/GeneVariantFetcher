# Pre-dispatch metadata review disposition

Claude, Grok and Agy reviewed the same experiment-only binder, tests and frozen
budget hook before any new scores. Their outputs remain in this directory.

Accepted Claude's useful hardenings: verify indexed records before filtering by
paper, rejecting disagreement between actual and indexed context; wait for both
API-writing processes to exit; require no live reservations; bind and recheck
the budget digest in the metadata receipt. A regression test covers a dispatched
call whose index context falsely names another paper. This does not change any
scientific prediction, model request or acceptance gate.

Grok's error-prefix concern does not apply to the actual driver: it explicitly
stores `type(exc).__name__ + ": " + str(exc)[:600]`. Keep that exact prefix and
its exact agreement with the recorded refusal instead of weakening it to a
substring match. Keep the `budget_not_dispatched` status/note visible. Failure
to establish non-dispatch aborts binding; it never manufactures a charge or
converts a genuine timeout to zero.

Agy's three proposed relaxations were rejected. The incremental clinical driver
calls only Astra; Grok belongs to the already-completed baseline, so accepting
arbitrary model names would widen the exception unnecessarily. Skipping malformed
index lines could hide a real API call and conflicts with native lock integrity.
The frozen clinical driver makes one guarded call per triggered paper and does
not retry a pre-dispatch refusal; duplicate refusal events therefore remain an
ambiguity. SDK transport retries happen inside an existing reservation and must
never qualify as zero calls.

The binder's valid zero case requires all three: a hash-bound exact refusal
already referenced by that paper, no indexed API call for the paper, and no
reservation for that paper in the completed clinical arm. Actual timeout usage
stays unknown with its reservation retained. The regression run includes both
these cases and the native failed-call telemetry tests.
