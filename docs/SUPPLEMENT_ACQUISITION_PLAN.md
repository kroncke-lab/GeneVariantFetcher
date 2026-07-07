# Supplement Acquisition Historical Summary

Generated 2026-06-05 from the supplement-acquisition effort over the four
cardiac gold genes. This is a historical summary, not the active next-run plan.
Use [`TASKS.md`](../TASKS.md) for the forward checklist,
[`RECALL_STATUS.md`](RECALL_STATUS.md) for current measured status, and
[`RECALL_REFRESH_RUNBOOK.md`](RECALL_REFRESH_RUNBOOK.md) for the idempotent
refresh workflow.

## What Landed

- Wired supplement folding into `refresh_run_db` so on-disk supplements become
  visible during refresh.
- Added recursive supplement folding for nested extracted archives.
- Added `scripts/recall_audit/supplement_fold_gap.py` for fold-gap QC.
- Added `ElsevierAPIClient.download_supplements` and
  `scripts/fetch_elsevier_supplements.py` for Elsevier `mmc` recovery.
- Landed promoted canonical DB improvements on 2026-06-05 after gated scoring.

## Measured 2026-06-05 Effect

| Metric | Canonical | With Supplements | Delta |
| --- | --- | --- | --- |
| unique_variants | 82.2% (2473) | 83.8% (2523) | +50, +1.6pp |
| variant_rows | 75.5% (5160) | 78.3% (5350) | +190, +2.8pp |
| patients | 77.8% | 80.6% | +534 |
| affected | 76.6% | 78.7% | +264 |
| carriers MAE | 0.910 | 0.882 | better |

The main recall lever was not Karger/Sage Cloudflare bypassing. The useful gap
was Elsevier body-only API recovery dropping `mmc` supplement tables, plus stale
refresh paths that did not fold already-downloaded supplements back into
`FULL_CONTEXT`.

## Current Operational Path

Use `scripts/refresh_recall.py` for normal recall refreshes. It composes
publisher supplement fetches, corpus-to-harvest bridging, supplement folding,
acceptance-gated re-extraction, scoring, and optional canonical DB landing.

Use `scripts/fetch_elsevier_supplements.py` when Elsevier supplements need a
direct refresh. Springer/Wiley supplement download remains a future enhancement
that needs live verification before it should be documented as operational.

Durable historical details are preserved in [`RECALL_HISTORY.md`](RECALL_HISTORY.md),
especially the 2026-06-05 entry.
