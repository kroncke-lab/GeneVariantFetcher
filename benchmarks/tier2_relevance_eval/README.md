# Tier-2 relevance shadow evaluation

This harness replays a deterministic, stratified slice of a historical
discovery run through `InternFilter` with one alternate model. It does not
mutate the source run or change production routing.

The committed `cohort_pmids.tsv` locks three diagnostic groups:

- 50 `productive_positive` papers from the pinned review manifest; every paper
  has at least one extracted BMPR2 variant.
- 50 `high_confidence_negative` papers that historical GPT-5.6 Sol rejected at
  confidence 0.95 or higher.
- 50 `fail_open_boundary` papers whose raw historical rejection was changed to
  PASS by the deterministic target-gene/cohort guard.

These labels are not manual relevance gold. They are useful for checking
compatibility, regressions against the existing route, deterministic guard
behavior, exact token use, and latency. A production promotion still needs
manual adjudication or an independently labeled relevance set.

Run from the repository root with a fresh output directory:

```bash
.venv/bin/python benchmarks/tier2_relevance_eval/run_shadow.py \
  --positive-count 50 \
  --negative-count 50 \
  --boundary-count 50 \
  --workers 12 \
  --outdir validation_runs/bmpr2_tier2_luna_shadow
```

The evaluator writes `cohort.json`, `results.jsonl`, `summary.json`, and exact
LLM traces. It refuses a non-empty output directory so separate runs cannot be
silently mixed. By default it replays the committed PMID/group manifest; the
historical source run supplies the original title, abstract, decision, and exact
Sol trace telemetry.

The CLI default requests `max`. On the Azure deployment tested on 2026-08-10,
GVF's provider compatibility layer maps that request to `xhigh`, the highest
value the deployment accepts. A direct literal-`max` smoke was rejected. Always
use the `luna.reasoning_efforts` trace summary as the deployment-effective value.

The first locked run and decision are in
`runs/20260810_bmpr2_luna_max/SUMMARY.md`.
