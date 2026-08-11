# Luna xhigh count-semantics result — 56 papers

> **Historical mixed-provenance result.** The active cohort was narrowed on
> 2026-08-11 to the cardiac 48 plus Nate's two lead-approved BRCA2 papers. See
> `../20260811_collaborator_gold_50/`; do not describe this 56-paper result as
> collaborator-gold performance.

## Outcome

Carrier MAE fell from **0.8148 to 0.0794** on the exact A1 56-paper output, a
90.3% reduction in absolute error (308 → 30). The locked predictions still
contribute all 378 observed carrier counts. Count recall is effectively flat,
moving from 378/1111 (34.02%) to 378/1110 (34.05%) only because one duplicate
R176W gold row was consolidated.

| Slice | Carrier observations | Absolute error | MAE |
|---|---:|---:|---:|
| Cardiac 48 | 326 | 16 | 0.0491 |
| BRCA2 8 | 52 | 14 | 0.2692 |
| Combined 56 | 378 | 30 | 0.0794 |

All three count fields are reported so carrier improvements cannot hide a
partition regression:

| Field | Before absolute error / MAE | After absolute error / MAE |
|---|---:|---:|
| Carriers | 308 / 0.8148 | 30 / 0.0794 |
| Affected | 240 / 0.7869 | 241 / 0.7902 |
| Unaffected | 51 / 0.1802 | 37 / 0.1307 |

The one-point affected regression is expected: the source-supported R420Q
mother is symptomatic, while the locked prediction counted only the proband.
This is evidence that adjudication did not simply follow predictions.

The dominant apparent carrier errors were count-scope defects in the answer
key, not extraction failures. The scorer also ignored an existing S1103Y v2
adjudication. Locked prediction digests did not change.

## Luna routing and cost signal

Seven compact Luna xhigh verification/debug calls consumed 27,682 total tokens
in 80.9 seconds. The 64k-budget broad recovery pilot had already consumed
153,010 tokens across six attempted calls while grounding zero additions in its
completed batches. Compact ambiguity cards used about 82% fewer tokens in this
diagnostic comparison and produced useful decisions; the two workloads are not
an apples-to-apples model-efficiency benchmark.

No Anthropic model participated in the original Luna adjudication run. A later
independent trust-chain and blind-source review used Claude Fable 5 Max together
with Grok 4.5 High and AGY Gemini 3.1 Pro High.

## Reproducible local artifacts

- Cardiac report: `validation_runs/20260810_failure_routing_a1_56/fixed48_all_layers/report.json`
- BRCA2 report: `validation_runs/20260810_failure_routing_a1_56/brca2_all_layers/report.json`
- Luna traces: `validation_runs/20260810_failure_routing_a1_56/luna_count_semantics_shadow_traces*`
- Source decisions: `benchmarks/count_semantics_eval/ADJUDICATIONS_20260810.md`
- Independent review: `benchmarks/count_semantics_eval/MULTIMODEL_REVIEW_20260810.md`
