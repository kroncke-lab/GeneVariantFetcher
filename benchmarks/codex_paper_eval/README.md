# Extraction-blinded paper evaluation

This harness compares paper-reading protocols on a fixed, gold-value-blinded
cardiac paper set. `prepare` may use gold only to determine PMID eligibility and
whether all three count fields have assertions. Extraction is finalized and
locked before `score` reads any gold value or row count.

## Current fixed sets

`highcarrier_48_papers_20260723.tsv` is the remembered 48-paper test: exactly 12
papers each for KCNH2, KCNQ1, SCN5A, and RYR2. It is cardiac-only; BRCA2 is not
part of this lock. BRCA2 currently lives in the broader curated diagnostic
fixture, whose overrides are not equivalent to the manual cardiac gold.

The most recent completed production artifact is
`runs/20260726_fixed48_production/`. A hermetic metric-only rescore (no LLM
calls and no local corpus lookup) is:

```bash
GVF_DISABLE_LOCAL_DATA=1 .venv/bin/python \
  benchmarks/codex_paper_eval/run_eval.py score \
  --run-dir benchmarks/codex_paper_eval/runs/20260726_fixed48_production
```

The 2026-08-08 rescore reproduced 789 TP / 985 FP / 212 FN: 78.8% recall,
44.5% raw precision, and 56.9% F1. Count values were supplied for 32.8% of
carrier assertions, 29.3% of affected assertions, and 23.2% of unaffected
assertions. This checks the current matcher/scorer against locked predictions;
it is distinct from a fresh extraction, which spends model tokens and receives
a new run identity, manifest, and traces.

## 2026-08-08 fresh snapshot replay

The current production path was also rerun from the exact preserved 48-paper
source selection with local-data discovery disabled. All four gene runs
completed with no stage failure or warning, and their trace indexes verified at
write time. The external evaluation projections are locked separately:

| projection | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| all layers | 831 | 1060 | 170 | 43.9% | 83.0% | 57.5% |
| paper-derived only | 710 | 592 | 291 | 54.5% | 70.9% | 61.7% |

Against the 2026-07-26 production lock, the all-layer replay gained 42 TP and
75 FP (recall +4.2 points, precision -0.5, F1 +0.6); the paper-only replay
gained 40 TP and 4 FP (recall +4.0, precision +1.2, F1 +2.4). Almost all of the
gain is RYR2 PMID 19926015, which the old markup circuit breaker discarded. The
fresh paper-only result for that PMID is 37 TP / 8 FP / 3 FN; adding external
linkage layers makes it 38 / 78 / 2, exposing a separate ClinVar-to-paper
precision problem.

Exact trace telemetry was 279 calls, 978,972 input tokens, 314,432 output
tokens, and 1,372,842 total tokens. The four production runs took about 17
minutes wall-clock in parallel (812–1014 seconds each); summed provider-call
duration was 3,670.5 seconds. The durable run record and reproducibility caveats
are in `runs/20260808_fixed48_snapshot_replay/`.

The native `prepare` → `extract` → `lock` → `score` harness does not read gold
values until after the lock. The archived `gvf-run` production projection is a
slightly different protocol: its recovery driver automatically writes
read-only per-layer scorecards when a registered gold CSV exists, before the
external evaluation lock. Those scores do not feed back into the DB, and the
ClinVar/PubTator layers use only PMIDs already observed in the DB unless the
explicit `--allow-gold-pmid-enrichment` switch is supplied (it was not). Thus
the production predictions are gold-independent, but a future projection
should suppress the intermediate scorecards as well when strict
lock-before-any-gold-read blinding is required.

## 2026-08-10 failure-routing A1 pilot

`runs/20260810_failure_routing_a1_56/` contains the reproducible scripts and
decision record for an experimental reason-class verifier gate on the same 48
cardiac sources plus the eight BRCA2 diagnostic entries. It reduced fixed-48
claim-verifier tokens 27.6% and total tokens 7.5%; variant metrics stayed within
the predeclared gates. Aggregate carrier MAE failed its gate, however, so
`ENABLE_TIER3_REASON_CLASS_ROUTING` remains default off pending a paired
same-primary-output ablation or independent locked replicate. The BRCA2
`gold_overrides` are curator/derived diagnostics, not manual headline gold.

## Fully traced run

```bash
.venv/bin/python benchmarks/codex_paper_eval/run_eval.py prepare \
  --seed 2026072301 \
  --paper-manifest benchmarks/codex_paper_eval/highcarrier_48_papers_20260723.tsv \
  --run-id my_traced_run

.venv/bin/python benchmarks/codex_paper_eval/run_eval.py extract \
  --run-dir benchmarks/codex_paper_eval/runs/my_traced_run \
  --model gpt-5.6-sol

.venv/bin/python benchmarks/codex_paper_eval/run_eval.py lock \
  --run-dir benchmarks/codex_paper_eval/runs/my_traced_run

.venv/bin/python benchmarks/codex_paper_eval/run_eval.py score \
  --run-dir benchmarks/codex_paper_eval/runs/my_traced_run
```

Every paper gets exact route/extraction call traces and explicit route/final
decision events under `llm_traces/<GENE>/<PMID>/`. The trace manifest is inside
the pre-gold lock and is revalidated during scoring. The self-contained
`llm_trace_report.html` is generated automatically and locked alongside it;
open that file directly in a browser to review the run paper by paper. See
[`../../docs/LLM_TRACING.md`](../../docs/LLM_TRACING.md) for the trace contract
and adjudication workflow.

Runs made before trace schema v2 cannot be treated as exact request/response
audits. Their final predictions and rationales remain useful, but rerunning the
fixed manifest is required to produce authentic raw call traces.
