# benchmarks/

Runnable evaluation sets for the GeneVariantFetcher (GVF) extraction pipeline.

## Start here: the three active tiers

[`evaluation_tiers/`](evaluation_tiers/README.md) is the only active cohort
index. Protocol rollout proceeds in order:

1. `gold_50`: 50 scored gene–paper attempts (48 cardiac + Nate's two BRCA2).
2. `cardiac_120`: 120 cardiac reviewer attempts (98 unique PMIDs).
3. `reviewer_396`: 396 attempts / 357 unique PMIDs in the eight established
   private reviewer workspaces.

No other directory is an active rollout cohort.

## Specialized harnesses and historical evidence

- [`codex_paper_eval/`](codex_paper_eval/README.md): extraction-blinded,
  hash-locked execution/scoring harness used by the first tier.
- [`curated_extraction_eval/`](curated_extraction_eval/README.md): the
  strategy-diverse regression fixture. It is useful for diagnostics but is not
  a fourth active rollout tier.
- [`count_semantics_eval/`](count_semantics_eval/README.md): frozen
  count-scope/MAE audit artifacts, including the historical 56-paper study.
- `cold_start_eval/` and `tier2_relevance_eval/`: specialized component tests,
  not paper rollout cohorts.

Dated runs remain immutable scientific evidence. New experiments belong in a
dated run directory and do not become active merely because they exist.
`scripts/run_recall_suite.py` plus `docs/RECALL_STATUS.md` remain authoritative
for the full four-gene headline metrics.
