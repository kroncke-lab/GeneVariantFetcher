# GPT-5.6 Luna Tier-2 relevance shadow

## Decision

**Do not promote Luna to the default Tier-2 route from this result alone.** Luna
passed this compatibility shadow against the historical Sol run and has a strong
estimated price advantage, but the cohort is diagnostic rather than manually
labeled. At maximum deployment-supported effort it also used materially more
output/reasoning tokens and was slightly slower. The next useful experiment is
the same locked cohort at `low` and `none`, followed by blind human adjudication
of disagreements and a sample of agreements.

No production route, model alias, prompt, threshold, guard, or default changed.

## Contract

- Source run: `results/BMPR2/20260807_163246`
- Historical route: `azure_ai/gpt-5.6-sol`, reasoning effort omitted
- Shadow route: `azure_ai/gpt-5.6-luna`
- Requested effort: `max`
- Deployment-effective effort: `xhigh`
- Cohort: 150 papers, 50 per diagnostic group
- Committed PMID/group lock: `cohort_pmids.tsv`
- Cohort SHA-256:
  `f052f3c87637e9527d8c484a6a452bb46bd1f0e58f214169c1a630142efb534d`
- Concurrency: 12 workers
- Luna wall time: 182.7 seconds
- Inputs, prompt/JSON contract, confidence threshold, disease context, and
  deterministic fail-open policy were held constant.

The current Azure deployment rejects literal `max` for this model and reports
`none`, `low`, `medium`, `high`, and `xhigh` as its supported values. GVF's
existing compatibility mapping therefore sent `xhigh`. This record uses
"requested max / effective xhigh" throughout; it does not claim that literal
`max` ran.

## Results

| Measure | Historical Sol | Luna requested max / effective xhigh | Change |
|---|---:|---:|---:|
| Calls / provider failures | 150 / 0 | 150 / 0 | — |
| Final routed-decision agreement | reference | 150/150 (100%) | — |
| Raw model-decision agreement | reference | 145/150 (96.7%) | — |
| Input tokens | 195,141 | 195,141 | 0% |
| Output tokens | 14,247 | 25,944 | +82.1% |
| Reasoning tokens | 3,909 | 14,808 | +278.8% |
| Total tokens | 209,388 | 221,085 | +5.6% |
| Median provider-call latency | 1.861s | 1.951s | +4.8% |
| p95 provider-call latency | 4.029s | 4.514s | +12.0% |
| Summed provider-call time | 335.2s | 350.3s | +4.5% |

Group outcomes:

| Diagnostic group | Expected final decision | Luna expected rate | Raw PASS / FAIL | Guard overrides |
|---|---|---:|---:|---:|
| Productive positive | PASS | 50/50 | 49 / 1 | 1 |
| High-confidence negative | FAIL | 50/50 | 0 / 50 | 0 |
| Fail-open boundary | PASS | 50/50 | 4 / 46 | 46 |

The five raw disagreements did not alter final routing. Luna passed four
historical fail-open boundary papers directly. It rejected one productive paper,
PMID 23592887, but the unchanged deterministic BMPR2 signal guard retained it.
That protected miss is the main reason not to promote from this non-gold shadow
alone.

## Public-list-price projection

OpenAI's public API list prices on the run date were $5 input / $30 output per
million tokens for Sol and $1 / $6 for Luna. Applying those rates to exact token
telemetry estimates $1.403 for the 150 historical Sol calls and $0.351 for the
Luna shadow: **75.0% lower**, despite Luna's extra output. Scaling the same mix
linearly to the source run's 974 Tier-2 calls gives about $9.11 versus $2.28.

This is not an Azure invoice. Contract prices, caching, future rate changes, and
workload mix can change the realized saving. The exact tokens and latency above
are the durable measurements; dollars are a labeled public-price projection.

Official model/pricing source:
<https://developers.openai.com/api/docs/models/gpt-5.6-luna>

## Limitations and next gate

- The positive label means the existing pipeline extracted at least one BMPR2
  variant; it is not a manual relevance judgment.
- Negative labels inherit historical Sol decisions, so agreement can only show
  compatibility with the incumbent—not independent correctness.
- Fail-open boundary papers deliberately test deterministic protection, not raw
  classifier accuracy.
- The fixed cardiac 48 plus eight BRCA2 diagnostic papers cannot evaluate this
  stage because explicit PMID manifests bypass Tier 1 and Tier 2.
- Run `low` and `none` on the identical cohort and require zero provider errors,
  no loss in final productive-positive retention, and manual review of every
  disagreement before changing `TIER2_MODEL`.

Raw output lives under the ignored validation path
`validation_runs/20260810_bmpr2_tier2_luna_max/`; use `run.sh` to replay into a
fresh `_replay` directory. The local historical source run must be present.
