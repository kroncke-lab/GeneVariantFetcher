# Production evaluation gold-blinding failure — 2026-08-21

## Incident

The first current-code `gold_120` production attempt was stopped before lock and
declared invalid. Its extraction launcher did not supply gold explicitly, but
`gvf-run` automatically discovered
`gene_variant_fetcher_gold_standard/normalized/KCNH2_recall_input.csv` after
paper extraction and passed it to `run_all_layers.py`. That subprocess opened
gold and produced `score_0_baseline`, `score_1_clinvar`, and
`score_2_pubtator` before the operator noticed and terminated the four gene
processes.

The recovery inputs were still DB-PMID scoped and the score values were not
intended to feed later layers. That does not rescue the evaluation: the harness
contract says no gold values are opened before predictions are finalized and
locked. The attempt is therefore unusable as acceptance evidence.

## Containment

- All four processes were terminated before any `RUN_STATUS.json` was written.
- No prediction lock or score report was created.
- The invalid attempt had reached 28 KCNH2, 18 KCNQ1, 18 RYR2, and 16 SCN5A
  provider extraction starts. Those paid calls are not reused in the accepted
  evaluation.
- The run directory was moved to the macOS Trash as recoverable incident
  evidence and is not part of the repository benchmark history.

## Fix

- `gvf-run --gold-free-run` disables automatic gold discovery for recovery,
  source diagnostics, and reporting, and refuses the gold-dependent v12 merge.
- The production-evaluation scaffold adds that flag to every gene command and
  its preflight refuses a launcher missing it.
- Regression coverage proves a gold-free run never calls `_find_gold` and passes
  `gold=None` into recovery layers.
- A completely new scaffold and extraction are required; the interrupted
  outputs cannot be resumed or rebound into predictions.
