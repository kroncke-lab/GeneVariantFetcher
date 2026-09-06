# Cost basis and provider pricing check

The experiment's ledger uses its predeclared uncached-rate proxies, retains
failed/canceled-call reservations and includes a $25 uncertainty margin. It is
not an Azure invoice. CLI consultations are outside the authorized API envelope.

On 2026-09-06, the [OpenAI Astra model page](https://developers.openai.com/api/docs/models/gpt-6-astra)
listed $10 input and $50 output per million tokens, with $1 cached reads and
$12.50 cache writes. Its long-context rule above 272k input tokens uses 2x input
and 1.5x output; the ledger reserves a more conservative 2x output multiplier.
Our calls are below that threshold. The [Grok 4.6 model page](https://docs.x.ai/developers/models/grok-4.6)
listed $2 input, $6 output and $0.50 cached input per million tokens.
These first-party prices do not establish the user's Azure invoice rates.

The original returned-usage proxy does not apply cache discounts or cache-write
premiums. It therefore must not be called an exact bill or an upper bound for
every individual call. The separate uncertainty margin also covers possible
cache-write uplifts, alongside early smoke and uncaptured pilot retry usage.
Keep any cache-aware repricing as a separate diagnostic rather than silently
rewriting the experiment's original cost receipts.

Output usage includes reasoning. For Grok's separately reported reasoning,
the ledger uses the larger of reported completion/output tokens and
`total_tokens - input_tokens`. Failed or operator-canceled calls without usage
remain unknown, with reservations retained. Successful calls after observed
SDK retries retain the earlier-attempt reservation too.
