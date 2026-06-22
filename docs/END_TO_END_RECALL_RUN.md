# End-to-End Recall Runbook

This is the portable way to run the current GVF recall workflow from another
computer without committing local papers, databases, logs, or credentials.

Current status and metrics live in `docs/RECALL_STATUS.md`. For existing runs
with newly recovered source artifacts, use `scripts/refresh_run_db.py`; it
supersedes the older post-insttoken re-extraction experiment for normal work.

Scope: use this when moving or recreating the recall workflow on another
machine. Use `QUICKSTART.md` for ordinary local setup, `NEW_GENE_RUNBOOK.md` for
a new no-gold gene, and `RECALL_REFRESH_RUNBOOK.md` for ongoing refreshes of an
existing scored run.

## Repository Setup

Follow `docs/QUICKSTART.md` for clone, `.venv`, browser, and first-run setup.
Follow `docs/API_KEYS.md` for the `.env` credential list. For recall refresh
work, publisher credentials and institutional full-text access materially affect
coverage.

## Fresh Cold-Start Run

Prefer the `gvf <command>` console script. When it is not on `PATH`, use
`.venv/bin/python -m cli <command>` as the fallback.

For one gene:

```bash
gvf gvf-run KCNH2 \
  --email "$NCBI_EMAIL" \
  --output validation_runs/fresh_e2e \
  --max-pmids 10000
```

For all four genes from scratch:

```bash
for gene in KCNH2 KCNQ1 RYR2 SCN5A; do
  gvf gvf-run "$gene" \
    --email "$NCBI_EMAIL" \
    --output validation_runs/fresh_e2e \
    --max-pmids 10000
done

.venv/bin/python scripts/run_recall_suite.py --score \
  --genes KCNH2,KCNQ1,RYR2,SCN5A \
  --results-dir validation_runs/fresh_e2e \
  --outdir recall_metrics/fresh_e2e
```

This redownloads papers and regenerates databases. Use the next section if you
want to reuse the full-text corpus already present on this machine.

## Reusing Local Full-Text Artifacts

The current full-text corpus is intentionally not tracked by Git. To reuse
papers already downloaded on this machine, copy the ignored `corpus/` directory
and any specific ignored run directories you plan to refresh or score.

Example transfer:

```bash
rsync -a --info=progress2 \
  corpus/ \
  other-host:/path/to/GeneVariantFetcher/corpus/

rsync -a --info=progress2 \
  validation_runs/<run-to-reuse>/ \
  other-host:/path/to/GeneVariantFetcher/validation_runs/<run-to-reuse>/
```

Then refresh the copied run from its source artifacts into canonical extraction
JSON, rebuild the DB, and run DB-observed recovery layers:

```bash
.venv/bin/python scripts/refresh_run_db.py \
  --gene <GENE> \
  --run-dir validation_runs/<run-to-reuse>/results/<GENE>/<timestamp> \
  --replace-db
```

The legacy `scripts/run_insttoken_reextract_experiment.py` driver is retained
for historical comparisons, but it should not be the default path for consuming
newly recovered sources.

`scripts/run_full_recall_experiment.sh` is a wrapper around
`run_insttoken_reextract_experiment.py` and is part of the historical
post-insttoken re-extraction path; **it is not a fresh-run launcher**. Do not
use it to drive a new `gvf-run`. To run a fresh multi-gene job, invoke
`gvf gvf-run <GENE>` per gene (or script a shell loop over `gvf gvf-run`).

## Full-Text Acquisition Only

To audit or extend the full-text corpus without running extraction, scoring, or
DB migration:

```bash
.venv/bin/python scripts/fulltext_acquisition_pass.py --target gold
```

Copy already-downloaded usable full text from prior local runs into the current
canonical recall directories:

```bash
.venv/bin/python scripts/fulltext_acquisition_pass.py \
  --target gold \
  --consolidate
```

Use the current `.env` publisher credentials, including `ELSEVIER_INSTTOKEN`, to
try downloading the remaining missing or weak full-text files:

```bash
.venv/bin/python scripts/fulltext_acquisition_pass.py \
  --target gold \
  --consolidate \
  --harvest \
  --max-pmids 50
```

For no-gold acquisition coverage, prefer `--target tier2-pass`; it reads each
run's `pmid_status/filter_progress.jsonl` and targets papers that passed Tier 2
triage. For all discovered PMIDs rather than just gold-standard PMIDs, use
`--target discovered`. Either mode can mean thousands of network attempts, so
combine it with `--max-pmids` for bounded passes.

## Scoring and Comparison

The re-extract driver runs scoring automatically. Manual scoring is:

```bash
.venv/bin/python scripts/run_recall_suite.py --score \
  --genes KCNH2,KCNQ1,RYR2,SCN5A \
  --results-dir validation_runs/fresh_e2e \
  --outdir recall_metrics/manual_score
```

Rows-mode MAE:

```bash
.venv/bin/python scripts/recall_mae.py rows \
  --summary recall_metrics/manual_score/summary.json \
  --outdir recall_metrics/mae/manual_rows
```

Run-vs-run recall delta:

```bash
.venv/bin/python scripts/recall_mae.py runs \
  --before recall_metrics/post_insttoken_20260521/summary.json \
  --after recall_metrics/manual_score/summary.json \
  --outdir recall_metrics/mae/manual_runs
```

## What Must Stay Out of Git

These are intentionally ignored and should be copied separately only when a
machine needs the existing paper corpus:

- `.env`
- `validation_runs/`
- `recall_metrics/`
- `results/`
- `*.db`
- `*.log`
- downloaded full-text markdown, PDFs, and supplements

The committed repo should contain the code, gold-standard CSV inputs, run
scripts, and documentation needed to regenerate those local artifacts.
