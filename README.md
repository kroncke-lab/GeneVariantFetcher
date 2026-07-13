# GeneVariantFetcher (GVF)

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-2.1.0-orange.svg)](CHANGELOG.md)

Builds an auditable evidence database from gene-focused biomedical literature,
designed to run **autonomously at scale** — hundreds of genes and hundreds of
thousands of papers — with automated quality gates, not human review, as the
primary control.

## What GVF Does

Given a gene, GVF discovers relevant papers, downloads full text and
supplements, extracts variants and carrier/count evidence with an LLM-backed
pipeline, and writes a queryable SQLite database with per-fact provenance so
every extracted count is auditable back to its source table, row, and quote.

The design target is unattended operation across many genes: automated source
QC, acceptance gates, count guards, non-regression checks, and per-fact
provenance decide what to trust, so runs scale without a human in the per-paper
loop. Human adjudication is an **exception-only escape hatch** for the rare
marginal case — routed through
[Variant_Browser](docs/VARIANT_BROWSER_INTEGRATION.md), off the per-paper and
per-run critical path, never a required step. A per-fact **trust gate** (v1) sorts
every extracted fact into a **trusted** or **quarantine** tier using gold-free
structural checks (`pipeline/trust_gate.py`, default-on in `gvf-run`;
`scripts/trust_report.py` inspects the tiers). Making the trusted tier the default
that scoring and downstream tools consume — and calibrating the gate per gene
class — is in progress.

The current production entry point is:

```bash
gvf gvf-run <GENE> --email brett.kroncke@gmail.com --output ./results [--disease "<phenotype>"]
```

Source recovery, including paywall and supplement acquisition, runs by default.
Pass `--no-source-recovery` only for a fast PMC/free-text-only pass or a
controlled measurement run.

## Scope and Intended Use

- **Operating model:** unattended, high-throughput extraction across many genes
  and hundreds of thousands of papers. Automated quality gates are the primary
  control; human adjudication is an exception-only escape hatch for the rare
  marginal case, not a per-record or per-run step.
- **Research use only.** GVF is not a patient-facing or clinical-decision
  system. Do not treat extracted variants, carrier counts, or classifications as
  clinical-grade.
- **Validation surface:** metrics are strongest on the four cardiac gold genes
  (KCNH2, KCNQ1, SCN5A, RYR2), which are in-distribution from repeated tuning.
  Generalization to other gene classes — e.g. BRCA1/BRCA2, whose truncating and
  case-control evidence differs structurally from cardiac missense — and
  arbitrary-gene cold start are exercised separately by
  [`benchmarks/cold_start_eval/`](benchmarks/cold_start_eval/README.md); report
  cardiac-gold, generalization, and cold-start numbers distinctly.

## Source Access and Browser Cookies

Default source recovery loads cookies from your local Chrome profile
(`scripts/fetch_paywalled.py`) so authenticated/institutional publisher access
carries into the fetch. This reads your browser session — only run it on a
machine and account you control, and only where your institution's access terms
permit programmatic full-text/supplement retrieval. To disable cookie loading,
pass `--no-cookies` to `scripts/fetch_paywalled.py`; for a fast run with no
paywall/cookie access at all, use `gvf gvf-run --no-source-recovery`.

## Current Status

Live recall metrics and blockers live in
[`docs/RECALL_STATUS.md`](docs/RECALL_STATUS.md), rendered as a published status
dashboard at
<https://kroncke-lab.github.io/GeneVariantFetcher/dashboard/>. The active forward
checklist lives in [`TASKS.md`](TASKS.md). Do not copy recall numbers into new
docs; they drift.

**Recall, precision, and MAE are reported only for the four cardiac genes
(KCNH2, KCNQ1, SCN5A, RYR2)** — the only gene-disease pairs with a fully
human-curated, manually derived gold standard. Other genes (APOE, BRCA1, BRCA2,
MYBPC3) are review targets scored against curator/derived `gold_overrides`, not
counted in headline metrics. See `docs/RECALL_STATUS.md` for the scope rule.

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
gvf gvf-run KCNH2 --email brett.kroncke@gmail.com --output ./results --disease "Long QT Syndrome"
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

1. Discover PMIDs from PubMind and PubMed (and Europe PMC when `USE_EUROPEPMC=1`).
2. Fetch abstracts and apply Tier 1/Tier 2 relevance filtering.
3. Download full text and supplements from PMC, publisher APIs, Unpaywall, and
   authenticated recovery routes when configured.
4. Reuse or update the local source corpus under `corpus/`.
5. Extract variants, carrier counts, phenotypes, provenance, and evidence.
6. Migrate extraction JSON to SQLite.
7. Run DB-observed recovery layers and the default-on per-fact trust gate.
8. After trust gating, run the default-on `azure_ai/gpt-5.6-sol`/`xhigh`
   final per-paper sniff test (Step 3.8); it persists soft review results and
   never mutates extracted counts.
9. Produce optional recall scoring and report handoff artifacts.

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
.venv/bin/python -m pytest tests/unit -q
```

Live network and institutional-access checks are opt-in:

```bash
GVF_TEST_OUTPUT_DIR=/tmp/gvf_tests .venv/bin/python -m pytest -m requires_network tests/integration -q
```

Project dependencies are defined in [`pyproject.toml`](pyproject.toml). Install
with `pip install -e ".[browser,dev]"`. A pinned snapshot of a known-good
environment is committed as [`requirements.lock`](requirements.lock)
(`uv pip install -r requirements.lock`); each run records that lock's hash, the
git SHA, the prompt/extractor hash, and the resolved model routing in its
`run_manifest.json` under `provenance`.

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
