# End-to-End Recall Runbook

This is the portable way to run the current GVF recall workflow from another
computer without committing local papers, databases, logs, or credentials.

## Repository Setup

```bash
git clone https://github.com/kroncke-lab/GeneVariantFetcher.git
cd GeneVariantFetcher

python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[browser,dev]"
python -m playwright install chromium

cp .env.example .env
```

Edit `.env` with:

- `NCBI_EMAIL`
- at least one LLM key, usually `ANTHROPIC_API_KEY` for the current recall setup
- `ELSEVIER_API_KEY` plus `ELSEVIER_INSTTOKEN` for institutional ScienceDirect full text
- optional `NCBI_API_KEY`, `SPRINGER_API_KEY`, and `WILEY_API_KEY`

## Fresh Cold-Start Run

For one gene:

```bash
.venv/bin/python -m cli gvf-run KCNH2 \
  --email "$NCBI_EMAIL" \
  --output validation_runs/fresh_e2e \
  --max-pmids 10000
```

For all four genes from scratch:

```bash
for gene in KCNH2 KCNQ1 RYR2 SCN5A; do
  .venv/bin/python -m cli gvf-run "$gene" \
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

The current full-text corpus is intentionally not tracked by Git. To reuse the
papers already downloaded on this machine, copy these ignored repo-relative
directories to the same paths on the other computer:

```text
validation_runs/turnkey_e2e_20260518_213934/results/KCNH2/20260518_213938
validation_runs/20260517_203904/results/KCNQ1/20260517_204424
validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938
validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938
```

Example transfer:

```bash
rsync -a --info=progress2 \
  validation_runs/turnkey_e2e_20260518_213934 \
  other-host:/path/to/GeneVariantFetcher/validation_runs/

rsync -a --info=progress2 \
  validation_runs/20260517_203904 \
  other-host:/path/to/GeneVariantFetcher/validation_runs/
```

Then run the post-insttoken re-extraction driver:

```bash
bash scripts/run_full_recall_experiment.sh --dry-run
bash scripts/run_full_recall_experiment.sh --foreground \
  --genes KCNH2,KCNQ1,RYR2,SCN5A \
  --max-pmids 10000
```

If `screen` is installed and you want a detached run:

```bash
bash scripts/run_full_recall_experiment.sh \
  --genes KCNH2,KCNQ1,RYR2,SCN5A \
  --max-pmids 10000

bash scripts/run_full_recall_experiment.sh --status
```

Use `--kill-existing` only when you intentionally want to stop already-running
GVF jobs before launching a replacement. To call the Python driver directly,
use `scripts/run_insttoken_reextract_experiment.py`; add `--full-reextract`
when you want to back up and regenerate every extraction JSON from the full text
that is already on disk.

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

For all discovered PMIDs rather than just gold-standard PMIDs, use
`--target discovered`. That can mean thousands of network attempts, so combine
it with `--max-pmids` for bounded passes.

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
