# Extraction-blinded paper evaluation

This harness compares paper-reading protocols on a fixed, gold-value-blinded
cardiac paper set. `prepare` may use gold only to determine PMID eligibility and
whether all three count fields have assertions. Extraction is finalized and
locked before `score` reads any gold value or row count.

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

## Standard comparison set: 56 papers (48 cardiac + 8 BRCA2)

**As of 2026-08-10 the standard set for comparing paper-processing strategies
is the 56-paper manifest** `highcarrier48_plus_brca2_20260810.tsv`: the fixed
cardiac 48 plus the 8 BRCA2 gold papers from the curated extraction benchmark
(10398279, 15365993, 18489799, 21356067, 22655046, 25802882, 26833046,
26848529). A protocol change is evaluated by running it on this set (or on the
arm it affects) and comparing against the production baselines below. Use the
same `prepare`/`extract`/`lock`/`score` flow by passing the manifest as
`--paper-manifest`; the BRCA2-only arm is `brca2_8_papers_20260810.tsv`.
Both the original 48-paper manifest and these are frozen — do not edit them;
runs scored on a frozen manifest stay comparable across time.

Current-production baselines (the numbers a change must beat):

- Cardiac arm: `runs/20260726_fixed48_production` (full gvf-run pipeline,
  all layers) and `runs/20260726_fixed48_production_paperonly` (paper-derived
  layers only). Single-model references: `runs/20260724_fixed48_sol`,
  `runs/20260724_fixed48_grok`.
- BRCA2 arm: `runs/20260810_brca2_8_production` (same pipeline, same
  `--pmid-file --no-source-recovery` calibrated replay) and its
  `_paperonly` sibling.

Gold provenance differs by arm, and the harness records it per gene in
`selection.json` and `report.json` under `gold_sources`:

- Cardiac genes score against the manual, human-curated gold standard
  (`gene_variant_fetcher_gold_standard/normalized/`).
- BRCA2 has no manual gold; it scores against the curator-adjudicated answer
  key `benchmarks/curated_extraction_eval/gold_overrides/BRCA2_recall_input.csv`
  (Azure lead-approved adjudications). That is adjudicated, not manual gold —
  mirror `docs/RECALL_STATUS.md` scope rules: report BRCA2 results separately
  and never fold them into cardiac headline metrics.

Known open adjudication follow-ups on the BRCA2 key (26833046 count semantics,
26848529 subset) are tracked in the curated benchmark; treat BRCA2 count MAE
as provisional until they close. Seeded random sampling (`--per-gene` without a
manifest) remains cardiac-only: BRCA2 has too few gold papers to sample and
widening the pool would change historical seeded draws.
