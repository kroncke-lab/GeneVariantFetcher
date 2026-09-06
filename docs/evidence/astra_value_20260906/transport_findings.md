# Transport findings from the initial test and early remediation receipts

This is a **metadata-only snapshot of 26 completed receipt files**: all 24
initial requests and the first two completed remediation requests
(`scn200_table__grok_remediated` and `ryr304_relatives__grok_remediated`).
Remediation choices, messages, and clinical answers were not inspected. Pending
or subsequently completed remediation requests are outside this snapshot.
No API calls, endpoint configuration reads, or credential reads were performed.

**The observable failure is a client read deadline being reached before a
complete response is returned to the caller. The receipts do not identify why Azure/Grok took that
long.** There is no captured quota error, authentication error, or HTTP failure
response. Longer-deadline jobs demonstrably can return after three minutes, but
the changed remediation package does not isolate the cause of the earlier
timeouts.

## Captured outcomes

| Group | Receipts | HTTP 200 returned | Client ReadTimeout | Returned elapsed time |
| --- | ---: | ---: | ---: | --- |
| Initial Grok 4.6, both arms | 16 | 5 | 11 | 18.38–122.64 seconds |
| Initial Astra | 8 | 7 | 1 | 8.52–23.94 seconds |
| Early Grok remediation | 2 | 2 | 0 | 205.67 and 368.92 seconds |

All initial timeouts report `ReadTimeout` with `The read operation timed out`,
at approximately 180.15–180.27 seconds. None of these failed receipts contains
an HTTP status, response headers, provider request ID, or token usage. They
cannot be classified as HTTP 429, 5xx, output-cap exhaustion, empty model
answers, or model reasoning failures from the available metadata. Their
unreported usage remains unknown; reaching a local deadline is not evidence
that the provider stopped processing or incurred zero cost.

All receipts record zero retries. The initial requests used low reasoning
effort, a 4,096 completion-token cap, strict JSON schema, and a 180-second read
deadline. The first remediation requests used low effort, an 8,192 cap, JSON
object mode, a 1,200-second deadline, and revised output instructions. These
multiple changes prevent a clean schema-versus-cap-versus-deadline attribution.
Later concurrency amendments do not change what can be inferred from these two
early completed jobs.

## Rate limits and request identity

Successful initial Grok responses report 48–49 remaining requests out of 50 per
minute and 44,735–48,827 remaining tokens out of 50,000 per minute. Successful
initial Astra responses report 3,499 of 3,500 requests remaining and
3,493,495–3,498,704 of 3,500,000 tokens remaining. All successful receipts report
the abuse-penalty flag as false. The two early remediation Grok responses show
49 requests remaining and 46,536 or 47,294 tokens remaining.

These successful-response snapshots do not show exhausted advertised quota.
They do not prove that the failed requests had the same admission state or
exclude upstream routing, scheduling, concurrency, or temporary service
conditions. There is no defensible basis here to call this an authentication
problem or to assert that increasing quota alone will fix it.

The successful receipts retain both provider and Azure APIM request IDs for
support correlation. Examples:

| Receipt | Provider request ID | Azure APIM request ID |
| --- | --- | --- |
| Initial SCN5A 20031634 Grok repeat | `3eec2afb-8831-46bd-b78b-78a065c0b958` | `167a95cd-53c1-48e0-b404-297a7baf821c` |
| Remediated SCN5A 20031634 | `9c3d4af0-c3e5-4b79-aaa3-a9539739dc6e` | `addef35b-eb17-4b8b-aeb1-4a2ae3b56547` |
| Remediated RYR2 30403697 | `a7537dbb-3793-4c7b-a8d2-47f022634e53` | `d33f5b1a-e181-4cf9-b94e-b3d0924a311b` |

The client deadline failures provide no response request IDs to correlate in
the same way. Local receipt names and dispatch accounting can still identify
the attempted calls, but they do not supply missing provider-side timings.

## Output caps and elapsed time

The five successful initial Grok responses use 521–1,243 tokens when counted as
`total_tokens - prompt_tokens`; this includes separately reported reasoning
tokens. Their reported reasoning token counts are 175–404. The two early
remediation responses use 1,334 and 762 such tokens, including 245 and 270
reported reasoning tokens respectively. These successful requests are well
below both the old and new configured caps.

Therefore the returned requests do not supply evidence of cap exhaustion. The
failed requests have no usage or finish metadata, so cap exhaustion cannot be
either established or excluded for those attempts. Small returned token counts
paired with long elapsed times also do not tell us whether the delay occurred
before generation, between tokens, during provider processing, or elsewhere.

The remediated SCN5A request took 205.67 seconds, and the remediated relative
request took 368.92 seconds. Those actual executions exceeded the earlier
180-second deadline. This establishes that a three-minute cutoff can discard
responses that eventually arrive under the new package. It does not establish
that simply extending the old failed calls would have produced the same
responses, nor that strict schema was the culprit. The same two packets had
also returned successfully under strict schema in at least one initial Grok
attempt, so strict schema is not categorically unsupported.

## What the latency metadata can and cannot distinguish

Successful Astra responses include a provider `latency_checkpoint` object with
fields such as `pre_inference_ms`, `service_ttft_ms`, and `service_ttlt_ms`.
For these seven responses, reported `pre_inference_ms` ranges from 69 to 228
milliseconds, `service_ttft_ms` from 940 to 1,830 milliseconds, and
`service_ttlt_ms` from 7,746 to 22,336 milliseconds. These are useful provider
measurements for the successful Astra calls only; they do not describe Astra's
timeout or Grok's requests.

The successful Grok receipts do not include corresponding latency checkpoints.
None of the failed receipts includes provider latency measurements. The client
used a non-streaming request, so these artifacts also lack independently
observed first-token/inter-token timings. Consequently, the receipts cannot
partition Grok's time into queueing versus reasoning/generation. Token usage is
not a wall-clock trace, and timestamps alone do not resolve that missing
instrumentation.

## Practical interpretation

The supported explanation is: **Grok 4.6 has highly variable end-to-end response
time on this Azure path, and many requests exceeded the selected deadline.**
Astra returned more promptly on most initial requests, but it also had one
deadline failure. The underlying service reason remains unproven. Preserve the
unknown usage reserves, report operational failures separately from clinical
accuracy, and describe remediation as a combined operating-package test until
a controlled one-variable test or provider telemetry isolates the cause.
