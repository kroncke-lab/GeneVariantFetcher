# Does Astra justify routine count extraction?

**Decision: retain the existing default; routine Astra use is not justified by
these tests.** Astra can read difficult source material and sometimes completes
when Grok 4.6 stalls, but the apparent new counts in this panel do not add usable
coverage over the saved functioning Grok 4.3 protocol. The completed cheaper-reader
remediation test recovers every apparent Astra-only number. On seven matching
papers with known usage, Astra's returned-token proxy costs **5.42 times as much**
under the different request packages, with no additional
accepted count. Its demonstrated advantage here is faster completion under the
tested configurations.

This is an opened-source diagnostic, not a production rollout or a new recall
headline. No extraction defaults, literal validators, databases, or unopened
tranches changed.

## What was tested

Eight already-opened papers provide 43 fixed questions per arm: 24 positive
numeric references and 19 valid null references. The packets mix clinical tables,
relatives outside the numbered cohort, diagnostic timepoints, aggregate prose,
an abstract, a 179-person roster, and two sources with variants but no applicable
clinical counts. There are no zero-reference questions. The investigator supplies
target variants, population and endpoint; autonomous discovery of those elements
is outside this test.

Each packet receives independent Grok 4.6 low, repeated Grok 4.6 low, and Astra low
requests. All three receive the same source, prompt and strict provider JSON
schema, a 4,096 completion-token cap, a 180-second timeout and no retries.
Within-packet order is seeded; “first” and “repeat” identify arms rather than
necessarily temporal order. The repeated Grok request receives no prior answer.

Sources, references and requests were frozen before dispatch. All 24 outputs
were locked before clinical inspection or scoring. A source-only independent
review identified one interpretive reference, RYR2 19398417 affected=2; the
43-question primary remains immutable, with a separate 42-question sensitivity
excluding that field. An explained null on it is not a false positive.

After three settled calls, metadata-only observation of slow service led to a
recorded scheduling amendment from serial dispatch to three isolated workers.
No in-flight call was canceled and no scientific input changed. Scheduling and
shared service conditions limit latency comparisons.

## Locked primary result

| Outcome | Grok first | Grok repeat | Astra low |
| --- | ---: | ---: | ---: |
| Successful papers / 8 | 2 | 3 | 7 |
| Exact numeric references / 24 | 4 | 13 | 22 |
| Correct nulls / 19 | 6 | 6 | 15 |
| Wrong nonnull answers returned | 0 | 0 | 0 |
| Failed/missing question slots | 33 | 24 | 6 |
| Mechanically source-bound exact numeric answers | 4 | 12 | 22 |
| Correct numeric claims accepted by existing literal validator | 0 | 0 | 1 |

Failures never earn null credit. Source binding means allowed source IDs and a
verbatim quotation; it does not establish clinical meaning, complete cohorts or
correct timepoints. One Grok repeat quotation splices noncontiguous rows, while
its numeric answer remains correct. Astra's interpretive affected=2 answer
passes mechanical checks but its explanation does not independently prove the
clinical split. These are separate quality dimensions.

The nine Astra-only numeric candidates occur on three papers where **both Grok
calls timed out**, not where Grok completed with wrong numbers:

| Source question | Astra values, carrier / affected / unaffected | Increment over saved Grok 4.3 protocol |
| --- | --- | --- |
| RYR2 25435091, initial diagnostic test | 8 / 7 / 1 | All three already present. Only carrier=8 passes the current literal validator in this diagnostic. |
| RYR2 19398417, abstract | 4 / 2 / 2 | All three already present; affected=2 retains an interpretive caveat. |
| RYR2 25814417, living baseline, previous symptoms | 179 / 45 / 133 | Different cohort/endpoint from existing 185 / 97 / 62; not three recovered errors. |

The four documented-relative carrier counts absent from the saved baseline are
recovered by Astra **and both successful Grok arms**. They are not an Astra
advantage, and a whole-paper merge still requires the population contract.
The original 12-paper full-pipeline experiment likewise failed to justify Astra
as primary: its accepted independent overlay added one carrier field and no
A/U fields at about $6.41 returned incremental API proxy. In the separate
93-person transcription experiment, cheaper Grok and Astra derived the same
counts. Together these findings support a decision about current routine use,
not a claim that Astra can never help.

## Cheaper-reader remediation comparison

This separately frozen post-hoc arm covers all eight packets; all eight
responses were locked before inspection and scoring. It keeps source, questions, references and the local semantic
schema, expresses the contract inline, uses provider `json_object`, raises the
cap to 8,192 and uses the ordinary 1,200-second timeout. These settings change as
a package: this is neither an isolated schema experiment nor an unchanged
production run (production does not explicitly request `json_object`). Earlier
answers and reference values are not included in requests. Unknown charges stay
reserved under the same budget. After two settled serial calls, a second
metadata-only scheduling amendment resumes the remaining six with two isolated
workers and at least 22 seconds between subprocess launches. No in-flight call
was canceled. Azure key retrieval precedes HTTP dispatch, so that launch
spacing is not a guaranteed exact request interval.

| Outcome | Grok remediated package | Astra original package |
| --- | ---: | ---: |
| Successful papers / 8 | 8 | 7 |
| Exact positive references / 24 | 24 | 22 |
| Correct nulls / 19 | 19 | 15 |
| Wrong nonnull / failed query slots | 0 / 0 | 0 / 6 |
| Exact answers in 42-query sensitivity | 42 | 36 |
| Mechanically supported exact answers in sensitivity | 42 | 36 |
| Accepted correct literal counts | 1 | 1 |

Grok recovers all nine previously Astra-only numeric candidates. Eight have
mechanically valid evidence; the ninth is the already-disputed affected=2
inference, for which Grok also stitches a nonverbatim quote. This is excluded
from the 42-query sensitivity, not silently repaired. The only accepted count
is identical in both models, and already exists in the retained Grok 4.3
baseline. An Astra missing-only overlay on remediated Grok adds **zero raw
numeric answers and zero accepted counts**.

The complete eight-paper Grok package costs **$0.097632** returned proxy
(**$0.109575** with input premium) and has median elapsed time **305 seconds**.
On the **seven same papers where both packages completed with known usage**,
Grok costs **$0.088746 total** versus Astra **$0.480560 total**, with observed medians
**324 seconds versus 13 seconds**. Astra is 5.42 times the returned-proxy cost.
This conditional subset excludes Astra's failure; the complete panel above
retains it. An all-eight known-use price ratio would be misleading, so the
scorer suppresses it and retains Astra's unknown reservation separately.
These are packet-package measurements, not whole-paper pricing or a controlled
model-speed ablation. The eight papers are selected and correlated within
paper; 43 questions are not 43 independent trials.

[Remediation scores](remediated_scores.json), [lock](remediation_locked.json),
[baseline comparison](baseline_comparison.md), and
[independent semantic review](remediation_semantic_review.md) confirm 23 numeric
answers and 19 abstentions against the frozen sources after excluding the
interpretive field. The
[transport findings](transport_findings.md) preserve the distinctions.
[Azure deployment metrics](azure_transport_metrics.md) corroborate long remote
completion times and show no 429 series in the observed window, but do not
isolate queueing, reasoning or schema overhead. Deployment aggregates do not
replace our receipt counts or resolve individual timeout billing.

## Practical implication

Keep Astra out of the automatic regular protocol. Human-reviewed Astra low
after a properly configured cheaper call is an **untested option**, not a
fallback policy with demonstrated incremental yield. All eight remediated
cheaper calls completed, so this follow-up supplies no failed case on which
to validate that policy. Medium effort had no count advantage in the
previous paired roster test. No claimed corpus uplift or 75%-of-human ceiling
follows from any of these selected diagnostics.

The next acceptance gate is a source-validated derived-count contract shared
by readers: explicit person/cohort/variant ownership, endpoint and timepoint,
unknown cells, duplicate relatives and index cases, and source coordinates.
Use deterministic table handling and arithmetic when the representation permits
it. SCN5A 20031634 still exposes a table-representation gap; source-correct
sums cannot safely be admitted by simply relaxing the literal-quote validator.
A more expensive reader does not itself resolve that integration problem.

[Claude](reviews/claude_results.md), [Grok](reviews/grok_results.md) and
[Agy](reviews/agy_results.md) CLI reviewers all support the current no-routine-inclusion
decision. Their useful challenges concern the asymmetric request packages,
conditional price comparison and unproven fallback policy. Agreement among
reviewers is not additional empirical evidence. The
[review disposition](reviews/results_review_notes.md) records their challenges
and the rejected factual misreadings; raw receipts are preserved in that directory.

## Costs, verification and reproduction

The final audit verifies **32 API requests**, with 20 known-usage responses
and 12 unknown-usage timeouts. New returned API proxy is **$0.627676**;
input premiums and full unknown reservations bring new conservative accounting
to **$1.634202**. The original envelope now stands at **$149.4599652 / $150**,
with **$0.5400348** unallocated; the Azure invoice remains unreconciled.
[Budget summary](budget_summary.json) and [ledger](budget.json).

The prior $150 envelope carried $147.8257632 of
conservative accounting. This campaign reserves at most another $2.17, uses
Azure extraction only, retains full unknown-call reservations, and treats
returned token prices as proxies rather than invoices. CLI consultation costs
are excluded under the user's instruction. No new Anthropic extraction API was
used.

The full offline unit suite passed: **2,969 tests**, including 18 new diagnostic
scoring regressions. Wheel build and CLI help passed. The canonical phenotype
figure was rebuilt and visually checked; its historical and opened cohorts stay
separate, and these component questions are not pooled into it.

Artifacts: [frozen plan](plan.json), [reference review](reference_review.md),
[source adjudication](source_adjudication.md), [primary lock](outputs_locked.json),
[primary scores](source_scores.json), [semantic review](semantic_review.md),
[remediation plan](remediation_plan.json), and [verification](verification.json).
The source packets and response receipts are under `packets/` and `responses/`.
A stale settings description calling Sol the strongest available deployment was
corrected; its model default is unchanged. The unused health-probe entrypoint
was removed from the imported experimental client only after all requests
settled. The audit proves its executed functions unchanged and records the
narrow cleanup in `client_import_only_cleanup.json`; pre-format code snapshots preserve
exact execution text.

Launch scripts are historical paid-test records with overwrite/budget guards;
do not rerun them to reproduce a score. Offline scorers replay the frozen
receipts and retain failures. Full frozen-source hash verification additionally
requires the original local `validation_runs/model_routing_20260906` snapshot.

```bash
.venv/bin/python docs/evidence/astra_value_20260906/score.py --check
.venv/bin/python docs/evidence/astra_value_20260906/score_remediation.py --check
.venv/bin/python docs/evidence/astra_value_20260906/audit.py
```
