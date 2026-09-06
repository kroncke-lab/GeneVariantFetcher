# Azure transport metric snapshot

**Deployment-filtered Azure metrics corroborate substantial remote completion
latency for Grok 4.6. They do not show a 429 throttle event in this snapshot or
identify whether the delay comes from queuing, reasoning, decoding, or output
format handling.** This adds provider-side evidence to the local read timeouts;
it does not establish the fate or billing of any individual timed-out request.

The read-only diagnostic queried `magen-api-2-resource` in resource group
`MAGen`, region `eastus2`, for **2026-09-06 19:30:20–20:12:03 UTC**, with
one-minute aggregation. All reported request and timing series were filtered
by `ModelDeploymentName`; none is an unfiltered account-wide aggregate. These
are nevertheless deployment aggregates, not a request-ID join to this
experiment. Concurrent callers and metric ingestion lag are not controlled.
The query arguments, returned data, and selected metric definitions are saved
in [azure_transport_metrics_selected.json](azure_transport_metrics_selected.json).
No settings, keys, inference calls, or clinical response content were accessed.

## Grok 4.6

`ModelRequests` returned **11 status-499 events and 7 status-200 events**. No
429 or 5xx series appeared in the queried result. The totals align with the
11 initial client read timeouts and the 5 initial plus 2 early remediation
responses recorded in [transport_findings.md](transport_findings.md). This is
aggregate consistency, not proof that a particular 499 event belongs to a
particular local timeout. Absence of a 429 series here is not a universal
guarantee against throttling or queueing.

`TimeToLastByte` returned the following nonzero minute buckets:

| UTC minute | Average, seconds | Maximum, seconds |
| --- | ---: | ---: |
| 19:33 | 122.365 | 122.365 |
| 19:39 | 118.346 | 119.932 |
| 19:51 | 51.324 | 51.324 |
| 19:52 | 18.175 | 18.175 |
| 19:57 | 205.002 | 205.002 |
| 20:01 | 368.358 | 368.358 |

The magnitudes closely track the successful local receipt durations, including
the 205.673-second and 368.923-second early remediation returns. This supports
the inference that the long waits involve remote service completion, rather
than being explained solely by local application overhead. The 19:39 bucket
contains two status-200 events; the others contain one each. The metric is an
estimate for non-streaming calls, and zero minimum/padding points were not
treated as actual instantaneous responses.

`TimeToResponse` had nonzero bucket averages of approximately 107–452 ms, but
the captured definition describes a gateway-level approximation, recommends
it for PTU deployments, and explicitly qualifies non-streaming measurements.
It is insufficient to claim that model generation began within half a second
or that queuing was absent. `NormalizedTimeToFirstToken` had no nonzero values,
so it supplies no usable first-token signal. `NormalizedTimeBetweenTokens`
also carries a non-streaming estimate qualification; it cannot independently
partition reasoning, schema handling, or provider scheduling time.

## Astra and comparability

For `gpt-6-astra`, `AzureOpenAIRequests` returned **29 status-200 events and
1 status-499 event**; the status-400 series contained only zeros. No nonzero
429 or 5xx events were returned. The 29 successful events exceed the 7
successful initial Astra receipts, so this deployment aggregate cannot be
attributed wholly to the captured experiment. It must not replace the local
receipt denominator or be used to inflate Astra's measured success rate.

The queried generic model timing metrics and `AzureOpenAITimeToResponse`
returned no Astra timing series. Thus this snapshot provides no equivalent
Azure aggregate latency comparison for Astra. Its captured per-response
latency checkpoints remain separate evidence, as documented in
[transport_findings.md](transport_findings.md).

## Interpretation limits

The evidence favors a practical diagnosis of **responses taking longer than
the initial three-minute client read deadline**, with observed provider-side
completion times reaching over six minutes. It does not isolate a root cause:
the remediation package simultaneously changed the deadline, token cap,
output mode, and instructions, and subsequent concurrency changes add another
factor. Small successful outputs do not establish what happened inside failed
requests. No timeout usage or cancellation confirmation is recovered by these
aggregate metrics.

The diagnostic deliberately avoided the legacy `Latency`, `ClientErrors`, and
`ServerErrors` metrics, whose captured definitions say not to use them for
Azure OpenAI. Microsoft's monitoring guidance likewise directs latency
investigations to the dedicated response and token timing metrics.
[Microsoft monitoring guidance](https://learn.microsoft.com/en-us/azure/foundry-classic/openai/how-to/monitor-openai)
