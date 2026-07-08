# GeneVariantFetcher Handoff

## Active Checkout

The authoritative working directory is:

`/Users/kronckbm/GitRepos/GeneVariantFetcher`

Use that path for current GVF work unless Brett explicitly says otherwise. Do
not treat `.claude/worktrees/`, `.codex/worktrees/`, old `Projects/`, or remote
`/mnt/temp4/` copies as current; they are side worktrees or historical scratch
checkouts.

GVF extracts genetic variants, carrier counts, and phenotype data from
biomedical literature for the Kroncke Lab variant interpretation pipeline.

## Read First

- `TASKS.md` - active forward checklist and next-run plan. If a plan elsewhere
  conflicts with this file, `TASKS.md` wins.
- `docs/RECALL_STATUS.md` - live recall metrics, scored baseline artifacts,
  current blocker shape, and high-yield missing-PMID context. Do not copy recall
  numbers into this handoff.
- `docs/RECALL_HISTORY.md` - append-only benchmark and change history.
- `docs/RECALL_REFRESH_RUNBOOK.md` - idempotent re-run path when source access,
  papers, or recovery logic changes.
- `docs/NEW_GENE_RUNBOOK.md` - operational flow for a new gene-disease pair
  without a gold standard.
- `docs/QUICKSTART.md` - canonical local setup and first-run commands.
- `docs/API_KEYS.md` - required and optional credentials, including publisher
  access notes.
- `docs/ARCHITECTURE.md` - pipeline architecture, module map, model/provider
  settings, and reasoning-effort knobs.
- `docs/VARIANT_BROWSER_INTEGRATION.md` - publish/adjudication round trip with
  the sibling Variant_Browser app.
- `benchmarks/curated_extraction_eval/README.md` - fast curated extraction
  benchmark for prompt, harness, guardrail, and matcher changes.

## Turnkey Commands

```bash
gvf gvf-run <GENE> --email brett.kroncke@gmail.com --output ./results [--disease "<phenotype>"]
```

`gvf-run` runs doctor checks, extraction, source QC, source recovery, DB-observed
recovery layers, scoring/report handoff, and corpus sync. Source recovery
(paywall plus supplement acquisition) runs by default; pass
`--no-source-recovery` for a fast PMC/free-text-only pass or for calibrated
`--pmid-file` measurement runs.

Use the project virtualenv:

```bash
source .venv/bin/activate
```

Default local tests avoid live network calls (bare `pytest` resolves to this
same offline unit suite via `pytest.ini`):

```bash
.venv/bin/python -m pytest tests/unit -q
```

Live network/institutional checks are opt-in:

```bash
GVF_TEST_OUTPUT_DIR=/tmp/gvf_tests .venv/bin/python -m pytest -m requires_network tests/integration -q
```

## Operating Shape

- Corpus cache: `corpus/<GENE>/<PMID>/`, indexed by `corpus/INDEX.json` and
  `corpus/INDEX.csv`, managed by `scripts/build_source_corpus.py`, and
  gitignored. `gvf-run` reuses usable cached source and folds new fetches back by
  default through corpus sync.
- Existing-run refresh: use `scripts/refresh_run_db.py`; do not patch SQLite rows
  directly to consume recovered source.
- Full recall rerun: use `scripts/refresh_recall.py`, then update
  `docs/RECALL_STATUS.md` for current numbers and append durable historical
  context to `docs/RECALL_HISTORY.md`.
- Scoring: `scripts/run_recall_suite.py` and `cli/compare_variants.py`.
- Recovery layers: `scripts/recall_recovery/run_all_layers.py`.
- Paywall recovery entry point: `scripts/fetch_paywalled.py`; authenticated
  browser strategies live under `harvesting/browser_html/`.
- Elsevier supplement recovery: `scripts/fetch_elsevier_supplements.py`.
- Review DB publish/adjudication: `gvf-run --publish-review` and
  `scripts/ingest_review_adjudications.py`; the full contract lives in
  `docs/VARIANT_BROWSER_INTEGRATION.md`.

## Files To Know

- `cli/gvf_run.py` - turnkey orchestration.
- `cli/extract.py` - lower-level extraction command.
- `cli/compare_variants.py` - gold-standard matcher and recall summary.
- `cli/dashboard.py` - static HTML coverage, provenance, and adjudication views.
- `pipeline/steps.py` - reusable pipeline step implementations.
- `pipeline/filters.py` - Tier 1/Tier 2 paper relevance filtering.
- `pipeline/table_router.py` and `pipeline/extraction.py` - clinical table and
  variant extraction.
- `pipeline/prompts.py` - extraction prompt templates.
- `harvesting/orchestrator.py` - full-text and supplement download coordination.
- `harvesting/elsevier_api.py`, `harvesting/springer_api.py`,
  `harvesting/wiley_api.py` - publisher integrations.
- `harvesting/figure_text_extractor.py`, `harvesting/figure_variant_reader.py`,
  and `pipeline/pedigree_extractor.py` - figure/pedigree extraction paths.
- `harvesting/migrate_to_sqlite.py` - extraction JSON to SQLite migration.
- `config/settings.py` - environment-backed settings and provider-aware model
  resolution.
- `utils/llm_utils.py` - LLM calls, rate limits, token limits, and
  reasoning-effort plumbing.
- `utils/variant_normalizer.py` and `utils/variant_scanner.py` - notation cleanup
  and regex pre-scan.
- `gene_variant_fetcher_gold_standard/` - normalized recall inputs and source
  exports.
- `benchmarks/curated_extraction_eval/` - small additive gold-paper benchmark;
  expand this when newly discovered failures need cheap regression coverage.

## House Rules

- Keep `CLAUDE.md`, `AGENTS.md`, and `CODEX.md` as pointers. Do not put live
  metrics, dated plans, setup blocks, or API key inventories here.
- Use the project `.venv` when available.
- Do not commit `.env`, local `results/`, SQLite DBs, generated
  `recall_metrics/`, or agent scratch files (`HEARTBEAT.md`, `IDENTITY.md`,
  `SOUL.md`, `TOOLS.md`, `USER.md`, `.openclaw/`).
- Side folders/worktrees are fine for experiments, but useful work should be
  merged into the active checkout and pushed to `origin/main` before an agent
  calls it current.
- Pre-commit hooks may reformat staged files. Never use `--no-verify`; re-stage
  after hooks run.
