# GeneVariantFetcher Handoff

## Active Checkout

The authoritative working directory is:

`/Users/kronckbm/GitRepos/GeneVariantFetcher`

Use that path and the `main` branch for current GVF work unless Brett explicitly
says otherwise. The 2026-08-12 handoff intentionally has no side worktrees or
local feature branches. Old `Projects/` or remote `/mnt/temp4/` copies are
historical scratch checkouts, not current sources.

GVF extracts genetic variants, carrier counts, and phenotype data from
biomedical literature for the Kroncke Lab variant interpretation pipeline.

## Read First

- `TASKS.md` - active forward checklist and next-run plan. If a plan elsewhere
  conflicts with this file, `TASKS.md` wins.
- `docs/README.md` - documentation authority map: current operating references,
  append-only history, and dated evidence.
- `docs/RECALL_STATUS.md` - live recall metrics, scored baseline artifacts,
  current blocker shape, and high-yield missing-PMID context. Do not copy recall
  numbers into this handoff.
- `docs/RECALL_HISTORY.md` - append-only benchmark and change history.
- `docs/PROTOCOL_CHANGELOG.md` - per-iteration ledger of protocol changes (one
  row per protocol-affecting PR). Append here every iteration; keep it in sync
  with RECALL_HISTORY (metrics) and ARCHITECTURE (current shape).
- `docs/PROTOCOL_COST_EVAL.md` - dated cost/quality measurements. It is evidence,
  not a statement of current defaults; `TASKS.md` owns the next acceptance gate.
- `docs/RECALL_REFRESH_RUNBOOK.md` - idempotent re-run path when source access,
  papers, or recovery logic changes.
- `docs/NEW_GENE_RUNBOOK.md` - operational flow for a new gene-disease pair
  without a gold standard.
- `docs/QUICKSTART.md` - canonical local setup and first-run commands.
- `docs/API_KEYS.md` - required and optional credentials, including publisher
  access notes.
- `docs/ARCHITECTURE.md` - pipeline architecture, module map, model/provider
  settings, and reasoning-effort knobs.
- `docs/EXTRACTION_CONTRACT.md` - meta prompt: what extraction must capture and
  must refuse, the reference prompt, and the map from each rule to the trust-gate
  reason code that backs it (or the note that nothing does).
- `docs/LLM_TRACING.md` - per-call prompt/response and decision trace contract,
  benchmark locking, and the curator adjudication workflow.
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

- Local corpus storage on Brett's current workstation: the repo path `corpus/`
  is an absolute symlink to
  `/Volumes/Ezekers/ResearchData/GeneVariantFetcher/corpus`. Keep code and
  commands on the stable `corpus/` interface (or use `GVF_CORPUS_DIR` only for
  a deliberate override). Before a corpus-reading or corpus-writing job, run
  `test -L corpus && test -d corpus`; if it fails, mount the APFS volume named
  `Ezekers` at `/Volumes/Ezekers`. Do not replace a broken symlink with a new
  local directory or rename the volume. The symlink is intentionally local-only
  and untracked, so a fresh checkout needs it recreated after the external
  target has been verified. `.gitignore` ignores it with an anchored, slashless
  `/corpus` rule — a trailing-slash `corpus/` pattern matches directories only,
  so a recreated symlink would show as untracked in a fresh clone.
- Corpus write guard: `utils/local_storage.py` enforces the above in code. A
  *missing* `corpus` link does not fail on its own — `mkdir(parents=True)` would
  create a real local `corpus/` and the run would build a second corpus on the
  internal disk that no later run reads, re-fetching paywalled source already
  cached on the volume. `require_external_storage()` refuses that write, so
  `scripts/build_source_corpus.py` exits 2 with the mount instructions instead.
  Set `GVF_ALLOW_LOCAL_CORPUS=1` to opt into plain local storage on a machine
  with no external volume. This is the write-side companion to
  `GVF_DISABLE_LOCAL_DATA`, which suppresses read-side discovery in the offline
  suite.
- `gvf-run` doctor (Step 1) sorts `corpus/` into four states and blocks on two of
  them. *dangling* — the volume was expected and is not mounted (drive removed,
  wrong machine) — blocks, so the run stops before Step 2 instead of re-fetching
  the whole corpus over the network. *local* — a real directory where the link
  should be — also blocks, because writing there builds a second corpus on the
  internal disk; set `GVF_ALLOW_LOCAL_CORPUS=1` to opt in on a machine that
  legitimately has no external volume (CI, a collaborator's laptop). *linked*
  and *absent* do not block: an absent link is a legitimate fresh checkout, and
  the write-side guard in `require_external_storage` catches it later with mount
  instructions. `GVF_CORPUS_DIR` bypasses the check entirely, and `--skip doctor`
  overrides it. The logic lives in `_apply_local_storage_checks`, split out of
  `doctor` so it is testable without doctor's network probes.
- When `corpus/` is unreachable but the run proceeds (`--skip doctor`), corpus
  reuse is disabled rather than fatal: `pipeline/steps.py::_resolve_corpus_dir`
  logs a `corpus reuse DISABLED` warning naming the volume, and `gvf-run`'s
  corpus sync records a stage warning. Such a run still completes, but re-fetches
  source it already has.
- A `VARIANTFEATURES_DB` that is set but not a readable file (typically the
  volume holding it is unmounted) now logs a warning in
  `utils/gene_metadata.py` before falling back to built-in gene metadata. The
  fallback is intentional; doing it silently was the problem.
- Any helper that *guesses* a local data path must go through
  `local_data_discovery_disabled()` — see `tests/unit/test_local_data_guard.py`
  for the three-branch pattern (guard on, guard off, explicit wins). This is not
  cosmetic: `scripts/backfill_paper_metadata.py::run_backfill` defaulted to
  `REPO / "corpus"` and `read_text()`-ed every `*_artifacts.json` under it, so
  once the corpus moved to an external volume the offline suite spent minutes
  walking it over USB — `tests/unit/test_gvf_run_pipeline_wiring.py` hung
  outright. With the guard applied, `pytest tests/unit` runs in ~43s.
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
- Linked-supplement recovery: `scripts/fetch_linked_supplements.py` — fetches
  supplements a paper's markup advertised but that never landed on disk. Runs
  by default inside `gvf-run` source recovery; use the script directly for a
  corpus-wide backlog sweep (`--dry-run` first).
- Review DB publish/adjudication: `gvf-run --publish-review` and
  `scripts/ingest_review_adjudications.py`; the full contract lives in
  `docs/VARIANT_BROWSER_INTEGRATION.md`.
- Additive count recovery: `pipeline/count_recovery.py` (gvf-run Step 3.55,
  default OFF via `COUNT_RECOVERY_ENABLED`) with `scripts/recover_counts.py` as
  the standalone CLI. Always `--dry-run` first; it refuses to run under
  `--skip trust-gate`.
- LLM traces: written for the whole `gvf-run` lifetime by `utils/llm_trace.py`;
  the browser report is built by `utils/llm_trace_html.py`
  (`scripts/build_llm_trace_html.py` is a thin CLI). Contract, run isolation,
  integrity levels, and the report size policy are in `docs/LLM_TRACING.md`.

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

- Keep `AGENTS.md` and `CODEX.md` as pointers to this canonical handoff. Keep
  live metrics and dated plans out of `CLAUDE.md`; link to their authorities.
- Use the project `.venv` when available.
- Do not commit `.env`, local `results/`, SQLite DBs, generated
  `recall_metrics/`, or agent scratch files (`HEARTBEAT.md`, `IDENTITY.md`,
  `SOUL.md`, `TOOLS.md`, `USER.md`, `.openclaw/`).
- Maintain one checkout on `main`. Do not create local feature branches or
  experimental worktrees unless Brett explicitly changes this handoff policy.
- Pre-commit hooks may reformat staged files. Never use `--no-verify`; re-stage
  after hooks run.
