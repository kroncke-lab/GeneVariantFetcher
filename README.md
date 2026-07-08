# GeneVariantFetcher (GVF)

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-2.1.0-orange.svg)](CHANGELOG.md)

Automated extraction of human genetic variant carriers from biomedical
literature into a normalized SQLite database.

## What GVF Does

Given a gene, GVF discovers relevant papers, downloads full text and
supplements, extracts variants and carrier/count evidence with an LLM-backed
pipeline, and writes a queryable SQLite database for downstream curation and
analysis.

The current production entry point is:

```bash
gvf gvf-run <GENE> --email you@example.com --output ./results [--disease "<phenotype>"]
```

Source recovery, including paywall and supplement acquisition, runs by default.
Pass `--no-source-recovery` only for a fast PMC/free-text-only pass or a
controlled measurement run.

## Current Status

Live recall metrics and blockers live in
[`docs/RECALL_STATUS.md`](docs/RECALL_STATUS.md). The active forward checklist
lives in [`TASKS.md`](TASKS.md). Do not copy recall numbers into new docs; they
drift.

For fast regression checks while changing prompts, extraction logic, guardrails,
or matching, use
[`benchmarks/curated_extraction_eval/`](benchmarks/curated_extraction_eval/README.md).
It is the small additive gold-paper set; confirm headline changes with the full
recall scorer before claiming a metric.

## Start Here

| Need | Start here |
|------|------------|
| Install locally and run once | [`docs/QUICKSTART.md`](docs/QUICKSTART.md) |
| Configure credentials | [`docs/API_KEYS.md`](docs/API_KEYS.md) |
| Run a new gene-disease pair with no gold standard | [`docs/NEW_GENE_RUNBOOK.md`](docs/NEW_GENE_RUNBOOK.md) |
| Re-run recall after new papers, credentials, or recovery code | [`docs/RECALL_REFRESH_RUNBOOK.md`](docs/RECALL_REFRESH_RUNBOOK.md) |
| Run the portable recall workflow on another machine | [`docs/END_TO_END_RECALL_RUN.md`](docs/END_TO_END_RECALL_RUN.md) |
| Understand pipeline internals | [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) |
| Inspect output schema | [`docs/OUTPUT_FORMAT.md`](docs/OUTPUT_FORMAT.md) |
| Publish/adjudicate with Variant_Browser | [`docs/VARIANT_BROWSER_INTEGRATION.md`](docs/VARIANT_BROWSER_INTEGRATION.md) |
| Understand which support scripts to use | [`scripts/README.md`](scripts/README.md) |

## Minimal First Run

Follow [`docs/QUICKSTART.md`](docs/QUICKSTART.md) for the canonical install and
`.env` setup. Then:

```bash
source .venv/bin/activate
gvf gvf-run KCNH2 --email you@example.com --output ./results --disease "Long QT Syndrome"
```

The primary output is:

```text
./results/KCNH2/<timestamp>/KCNH2.db
```

## Main Commands

| Command | Purpose |
|---------|---------|
| `gvf gvf-run` | Turnkey run: doctor checks, extraction, source QC, source recovery, recovery layers, scoring/report handoff, and corpus sync |
| `gvf extract` | Lower-level extraction command for debugging individual stages |
| `gvf scout` | Data Scout on existing downloaded papers |
| `gvf audit-paywalls` | Summarize blocked full-text and supplement acquisition |
| `gvf dashboard` | Static HTML coverage, provenance, and adjudication dashboard |

Prefer `gvf gvf-run` for normal work. Use lower-level commands when debugging
or replaying a specific stage.

## Pipeline Shape

GVF's default workflow is:

1. Discover PMIDs from PubMind, PubMed, and Europe PMC.
2. Fetch abstracts and apply Tier 1/Tier 2 relevance filtering.
3. Download full text and supplements from PMC, publisher APIs, Unpaywall, and
   authenticated recovery routes when configured.
4. Reuse or update the local source corpus under `corpus/`.
5. Extract variants, carrier counts, phenotypes, provenance, and evidence.
6. Migrate extraction JSON to SQLite.
7. Run DB-observed recovery layers and optional recall scoring.

Technical details live in [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).

## Output

A run directory contains the SQLite database, workflow logs, per-paper
extraction JSON, downloaded source, source-QC artifacts, and reports:

```text
results/<GENE>/<timestamp>/
├── <GENE>.db
├── <GENE>_workflow_summary.json
├── abstract_json/
├── extractions/
├── pmc_fulltext/
└── source_qc/
```

See [`docs/OUTPUT_FORMAT.md`](docs/OUTPUT_FORMAT.md) and
[`docs/SQLITE_MIGRATION_GUIDE.md`](docs/SQLITE_MIGRATION_GUIDE.md) for schema
details.

## Development

Use the project `.venv`:

```bash
source .venv/bin/activate
.venv/bin/python -m pytest tests/ -q
```

Live network and institutional-access checks are opt-in:

```bash
GVF_TEST_OUTPUT_DIR=/tmp/gvf_tests .venv/bin/python -m pytest -m requires_network tests/integration -q
```

Project dependencies are defined in [`pyproject.toml`](pyproject.toml). Install
with `pip install -e ".[browser,dev]"`; there is no separate requirements file.

## Repository Hygiene

Do not commit `.env`, local `results/`, SQLite DBs, generated `recall_metrics/`,
`validation_runs/`, the fetched `corpus/`, or agent scratch files. These are
intentionally gitignored.

Tracked runtime data lives in [`data/`](data/README.md): variant aliases and
reference protein sequences used by validation code.

## Citation

```bibtex
@software{genevariantfetcher,
  title = {GeneVariantFetcher: Automated Extraction of Genetic Variant Carriers from Biomedical Literature},
  author = {Kronck, Brett M. and Roden, Dan M.},
  year = {2026},
  version = {2.1.0},
  url = {https://github.com/kroncke-lab/GeneVariantFetcher}
}
```

MIT License; see [`LICENSE`](LICENSE).
