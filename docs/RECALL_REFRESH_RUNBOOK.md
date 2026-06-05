# Idempotent Recall-Refresh Runbook

How to re-run the recall pipeline over time — when **more papers are published**
or **API permissions expand** (a new Springer/Wiley key, off-VPN access) — and
have it pick up *only the new work*. Every step here is idempotent: re-running
on an unchanged corpus is a no-op.

History/benchmarks: [`RECALL_HISTORY.md`](RECALL_HISTORY.md). Current numbers:
[`RECALL_STATUS.md`](RECALL_STATUS.md). Supplement deep-dive:
[`SUPPLEMENT_ACQUISITION_PLAN.md`](SUPPLEMENT_ACQUISITION_PLAN.md). Cold-start a
brand-new gene: [`NEW_GENE_RUNBOOK.md`](NEW_GENE_RUNBOOK.md).

## One command (the common case)

```bash
# measure what new supplements/papers would add (no canonical mutation):
python scripts/refresh_recall.py --gene SCN5A --run-dir <run-dir-with-extractions-and-SCN5A.db>
# measure AND promote into the canonical DB if it improves (backs up first):
python scripts/refresh_recall.py --gene SCN5A --run-dir <run-dir> --land
```

`refresh_recall.py` composes the idempotent pieces: fetch publisher supplements →
bridge the corpus to a flat harvest dir → `refresh_run_db` (fold on-disk
supplements + acceptance-gated re-extract of any paper whose folded source now
yields more variants) → score vs gold → optionally surgically land into the
canonical DB (preserving ClinVar/PubTator/figure layer rows, reverting on any
regression). If nothing new is found, it reports a no-op and exits.

## The two growth scenarios

### A. New papers were published
1. Re-run discovery + harvest for the gene (adds new papers; reuses cached source
   for old ones — `gvf-run` is corpus-as-cache):
   ```bash
   gvf gvf-run <GENE> --email you@example.com --output results/
   ```
   New full text/supplements/figures fold into `corpus/<GENE>/<PMID>/` via the
   `--corpus-sync` step (default on); already-cached usable papers are skipped.
2. `python scripts/refresh_recall.py --gene <GENE> --run-dir <run> --land`
   re-extracts only the new/changed papers and lands the gains.

### B. A new/upgraded API key or access (e.g. Springer key, Wiley off-VPN)
1. Add the key to `.env` (see [`API_KEYS.md`](API_KEYS.md)).
2. Re-fetch supplements — previously-blocked papers are now reachable, already-
   fetched ones are skipped:
   ```bash
   python scripts/fetch_elsevier_supplements.py --gene <GENE>     # Elsevier mmc (today)
   # Springer/Wiley: add SpringerAPIClient.download_supplements (see SUPPLEMENT_ACQUISITION_PLAN.md
   # "T6" — reactivates the moment the Springer key is restored)
   ```
3. `python scripts/refresh_recall.py --gene <GENE> --run-dir <run> --land`.

## Why each step is idempotent (safe to re-run)

| Step | Tool | Idempotency guarantee |
|---|---|---|
| Source acquisition | `gvf gvf-run` (corpus-as-cache) | Reuses usable cached source per PMID; a `stub` falls through to a fresh fetch so a new key re-attempts it. |
| Corpus build | `scripts/build_source_corpus.py --apply` | Never-downgrade: adds new papers, upgrades only compromised categories; a clean re-run is a no-op. |
| Supplement fetch | `scripts/fetch_elsevier_supplements.py` | Skips papers that already have supplements on disk and files already downloaded; converges to 0-new. |
| Supplement fold | `harvesting/supplement_fold.py` (wired into `refresh_run_db`) | Sentinel-delimited, backs up once to `.pre_fold_bak`, rebuilds the folded block each pass — never double-appends. |
| Re-extraction | `scripts/refresh_run_db.py` | Selects candidates by deterministic-variant lift; a paper already extracted at its current variant count is not re-selected. Regression + explosion acceptance gates. |
| Land into canonical | `scripts/refresh_recall.py --land` | Backs up the canonical DB first; deletes only extraction-origin rows (preserves layer rows); reverts if the landed DB scores below canonical. |

## Manual building blocks (if you want the steps individually)

```bash
# 1. fetch supplements the full-text API dropped (idempotent)
python scripts/fetch_elsevier_supplements.py --gene SCN5A
# 2. QC: how many on-disk supplements aren't yet folded into FULL_CONTEXT
python scripts/recall_audit/supplement_fold_gap.py --genes SCN5A
# 3. bridge nested corpus -> flat harvest dir
python scripts/corpus_to_harvest.py --gene SCN5A --out /tmp/harvest/SCN5A
# 4. refresh (folds supplements + acceptance-gated re-extract); writes a NEW db
python scripts/refresh_run_db.py --gene SCN5A --run-dir <run> \
    --harvest-dir /tmp/harvest/SCN5A --stage-extractions \
    --output-db /tmp/SCN5A.refreshed.db --skip-recovery
# 5. score
python scripts/run_recall_suite.py --score --genes SCN5A \
    --db SCN5A=/tmp/SCN5A.refreshed.db --outdir /tmp/score
```

## Notes
- The canonical DBs live under `results/` / `validation_runs/` and are **gitignored
  (local)**; the refresh updates them in place (with backups). Commit only the
  code + docs, never the DBs.
- Do **not** wholesale-rebuild a canonical DB with bare `migrate(extractions)` —
  it drops the ClinVar/PubTator/figure layer rows the canonical DB accumulated.
  Use `refresh_recall.py --land` (surgical, layer-preserving) instead.
- After landing, append the new benchmark to [`RECALL_HISTORY.md`](RECALL_HISTORY.md)
  and update the current numbers in [`RECALL_STATUS.md`](RECALL_STATUS.md).
