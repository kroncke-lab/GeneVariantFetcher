# GeneVariantFetcher Handoff

## Active Main Checkout

The active main working directory is:

`/Users/kronckbm/GitRepos/GeneVariantFetcher`

Use that path for current GVF work unless Brett explicitly says otherwise.
Do not treat `.claude/worktrees/`, `.codex/worktrees/`, old `Projects/`, or
remote `/mnt/temp4/` copies as authoritative. They are historical side
worktrees or scratch checkouts and can lag behind `main`.

Coordination rule: side folders/worktrees are fine for experiments, but useful
work should be merged into this checkout and pushed to `origin/main` before any
agent calls it the current implementation. If `main` is dirty, either commit and
push the coherent change set or state exactly what remains uncommitted.

GVF extracts genetic variants, carrier counts, and phenotype data from
biomedical literature for the Kroncke Lab variant interpretation pipeline. The
grant target is 90% unique-variant recall, with submission in October 2026.

## Current Measured State

Current live metrics, the scored baseline artifact, highest-yield remaining
PMIDs, and the next-run plan are consolidated in:

`docs/RECALL_STATUS.md`

Treat that file as authoritative. Older KCNH2-only and 2026-05-18 closeout
numbers are historical debugging baselines. Gold-PMID-conditioned enrichment
and KCNH2 v12 manual recovery are diagnostic/manual recovery paths, not
cold-start capability.

## What Changed In This Branch

- Source recovery stack: authenticated Playwright recovery in
  `scripts/fetch_paywalled.py`, publisher strategies under
  `harvesting/browser_html/`, quality gating, per-domain pacing, and PMC
  fallback.
- Clinical table extraction: `pipeline/table_router.py` and
  `pipeline/extraction.py` preserve clinical mutation-list tables and infer one
  carrier per clinical row when there is patient/proband/family context but no
  explicit count column.
- Recall runner: `scripts/run_recall_suite.py` scores normalized per-gene
  recall inputs from `gene_variant_fetcher_gold_standard/normalized/`.
- Turnkey driver: `gvf gvf-run <GENE> --email <email> --output <dir>
  [--disease "<phenotype>"]` runs doctor -> extraction -> source QC ->
  source recovery (the no-gold acquisition loop) -> DB-observed recovery layers
  -> report. As of 2026-06-03 source recovery is **on by default**: it fetches
  the source QC queue, summarizes selected-vs-successful acquisition, and
  acceptance-gated staged-refreshes accepted sources before scoring/reporting.
  Pass `--no-source-recovery` for a fast PMC/free-text-only pass or for
  calibrated `--pmid-file` measurement runs (staged refresh mutates the DB; it
  is backed up first). `--disease` scopes discovery + the Tier-2 filter to a
  gene-disease pair. `--with-v12` is KCNH2-only and opt-in.
- Matcher: `cli/compare_variants.py` has positional-digit guards, greedy
  one-to-one assignment, cDNA/protein bridging, and 2026-05-15 fixes for
  frameshift spellings such as `fsTer`, `fs/185`, `fs+*49`, and malformed
  `AlaX14`.
- Gold-standard builders: `scripts/build_gold_standard_from_varbrowser.py` and
  `scripts/build_ryr2_gold_standard_from_xlsx.py`.
- Reasoning effort + model registry (2026-05-29): per-stage
  `TIER2/TIER3/TABLE_ROUTER/VISION_REASONING_EFFORT` env knobs (default off), and
  `utils/llm_utils.py` now registers `gpt-5.5`, `grok-4.3`, and a generic `gpt-5`
  token-limit fallback. See "Model Reasoning Effort" below.

## Manual Acquisition Queue

The manual acquisition queue is the set of passed-filter PMIDs that do not have
usable assembled full text after the automatic harvesters and fallback routes.
They are not "assigned to a human" by code; they are papers GVF could not fetch
with the currently available credentials and network position.

Current high-yield missing PMIDs are listed in
`docs/RECALL_STATUS.md`. Keep PMID rankings there so this
handoff file does not drift.

## Current Recovery Path

**As of 2026-05-21, the Elsevier insttoken unblock has landed.** Vanderbilt's
institutional `X-ELS-Insttoken` is installed in `.env` and 242 of 246
previously-paywalled Elsevier articles across KCNH2/SCN5A/RYR2/KCNE1/KCNQ1 now
return full text via `harvesting/elsevier_api.py`. The unlocked bodies were
saved as `{PMID}_FULL_CONTEXT.md` into each run's existing `pmc_fulltext/`
directory; pre-token stub files were preserved as `*.pre_insttoken_bak`.

Do not patch SQLite rows directly to consume newly recovered sources. Use
`scripts/refresh_run_db.py` for existing runs: it selects stale or under-counted
PMIDs from source artifacts, rewrites canonical extraction JSON, rebuilds a
fresh DB from the full extraction directory, and then runs DB-observed recovery
layers. Gold standards are optional; without gold, the same recovery layers run
and scoring is skipped.

As of 2026-06-01, `gvf-run` writes `source_qc/` artifacts for genes without a
gold standard: a source acquisition worklist, a fetch queue, a refresh source
override CSV, and coverage summary. Use
`scripts/recall_audit/summarize_acquisition_outcome.py` after
`fetch_paywalled.py` to distinguish PMIDs selected for acquisition from PMIDs
that actually landed usable full text; with a gold/report denominator, it also
reports both quantities as PMID recall. The source QC default denominator is
the harvest/extraction surface, not raw PubMed discovery; use
`--include-discovery-pmids` only for intentionally broad diagnostics.
End-to-end no-gold source recovery now runs by default in `gvf-run` (it
composes source QC, paywall fetch, outcome summarization, and acceptance-gated
staged refresh). Pass `--no-source-recovery` to skip it.

The same path has been exercised on KCNH2, RYR2, and SCN5A with staged
extraction replay. **All current recall figures live in
`docs/RECALL_STATUS.md`** (the single source of truth) — including
selected-vs-actually-downloaded PMID recall, refresh-successful PMID recall,
and acceptance-gate details. Do not restate per-gene numbers here; they drift.
Note that the published baselines are cardiac-gene-specific, depend on the
Vanderbilt Elsevier insttoken, and were scored with figures skipped — an
arbitrary new gene-disease pair starts with no gold standard and no recall
measurement.

Residual-gap diagnoses and the latest per-gene numbers are tracked in
`docs/RECALL_STATUS.md` (Next Run Plan + session log), not here. Durable tooling
that came out of those sessions: a Word `.doc`/`.docx` supplement-conversion
fallback (macOS `textutil` recovers table rows `antiword` misses),
scanner gene-context filtering, a stricter protein-notation artifact guard, a
range-deletion parser, and `scripts/refresh_run_db.py --replay-model` for
bounded replays when the default Tier 3 model is slow or experimental. The
acquisition outcome summary reports selected-for-fetch,
usable-fulltext-downloaded, source-refresh-attempted, and
source-refresh-successful PMID recall separately. The recurring residual
pattern is supplement acquisition (Cloudflare/redirect-blocked `mmc1.docx`
files), not generic underextraction.

Residual non-Elsevier paywalls (Karger Cloudflare, Sage CF fingerprint)
remain in `TASKS.md` Blocked. Wiley TDM was verified working off-VPN on
2026-05-26 — the earlier "revoked key" claim was stale; `harvesting/wiley_api.py`
plus the current `WILEY_API_KEY` returns full-text PDFs for `10.1002/humu.*`.
Some Wiley Online Library DOIs still return TDM 403 or Cloudflare in no-cookie
browser sessions, so treat Wiley as partially covered rather than solved.

## Run On Another Computer Or VPN

From a clean checkout:

```bash
python3.11 -m venv .venv   # Python 3.11+ (pyproject requires-python >=3.11)
source .venv/bin/activate
pip install -e ".[browser,dev]"
python -m playwright install chromium
```

Create `.env` with at least:

```bash
# One LLM provider key:
ANTHROPIC_API_KEY=...
# or OPENAI_API_KEY=...
# or AZURE_AI_API_KEY=...
NCBI_EMAIL=...
NCBI_API_KEY=...
ELSEVIER_API_KEY=...
SPRINGER_API_KEY=...
WILEY_API_KEY=...
# Optional but highest value:
ELSEVIER_INSTTOKEN=...
```

Default local tests avoid live network calls:

```bash
.venv/bin/python -m pytest tests/ -q
```

Live network/institutional checks are marked separately:

```bash
GVF_TEST_OUTPUT_DIR=/tmp/gvf_tests .venv/bin/python -m pytest -m requires_network tests/integration -q
```

Recover queued or paywalled papers after VPN/login access is available:

```bash
# <RUN> is a run dir produced by gvf-run, e.g. results/KCNH2/<TIMESTAMP>
.venv/bin/python scripts/fetch_paywalled.py \
  --input <RUN>/pmc_fulltext/paywalled_missing.csv \
  --output <RUN>/pmc_fulltext \
  --no-headless
```

Then rerun extraction/migration and score:

```bash
.venv/bin/python scripts/run_recall_suite.py --score --genes KCNH2 \
  --db KCNH2=results/KCNH2/20260506_102238/end_to_end_20260515_manual_recovery/KCNH2_v12_manual_recovery_20260515.db \
  --outdir recall_metrics/kcnh2_after_vpn_recovery
```

## Model Reasoning Effort (optional, off by default)

OpenAI-style reasoning models (gpt-5 / o-series) expose a reasoning-effort knob.
GVF can set it per pipeline stage via env vars; all default to unset = the
provider default, so behavior is unchanged until you opt a stage in:

- `TIER2_REASONING_EFFORT` — Tier 2 relevance (`pipeline/filters.py`). High
  volume; `minimal` is a cost/latency saver, not a quality lever.
- `TIER3_REASONING_EFFORT` — Tier 3 extraction (`pipeline/extraction.py`) and the
  claim adjudicator, which inherits the Tier 3 value. This is where
  unique-variant recall and carrier-count accuracy are decided.
- `TABLE_ROUTER_REASONING_EFFORT` — table classification
  (`pipeline/table_router.py`).
- `VISION_REASONING_EFFORT` — figure/pedigree models, both the chat-completions
  and the Azure Responses-API paths (`harvesting/figure_text_extractor.py`,
  `harvesting/figure_variant_reader.py`, `pipeline/pedigree_extractor.py`).

Values: `minimal | low | medium | high` (validated in `config/settings.py`). The
knob is plumbed through `utils/llm_utils.build_reasoning_effort_kwargs` (chat
completions) and `build_responses_reasoning_param` (Responses API). Both are
no-ops for models without an effort knob (Anthropic, Grok, Kimi) — point a stage
at a gpt-5/o-series model for the setting to take effect. Anthropic extended
`thinking` is intentionally not wired (it requires `temperature=1`).

Treat effort as a second-order lever: source acquisition and extraction logic
move recall more. Change one stage at a time and re-score with
`scripts/run_recall_suite.py` before keeping a non-default value.

New model ids are registered in `utils/llm_utils.MODEL_TOKEN_LIMITS`: `gpt-5.5`,
`grok-4.3`, and a generic `gpt-5` fallback so future gpt-5.x ids get a full
output-token budget instead of the 4k default.

## Active Work

Current active work is tracked in
`docs/RECALL_STATUS.md`. Do not duplicate metric gaps here.
The short version: source acquisition improved after the Elsevier insttoken
unblock, and the refresh/rebuild path now exists for converting recovered
source artifacts into DB rows without gold-standard-dependent logic. Beyond
that: fix general count/affected-status extraction issues, address the
remaining Wiley/Karger/Sage paywalls, and keep no-gold QC prominent for new
genes.

## Files To Know

- `cli/compare_variants.py` - gold-standard matcher and recall summary.
- `scripts/run_recall_suite.py` - multi-gene recall scoring.
- `scripts/refresh_run_db.py` - safe source -> extraction JSON -> DB refresh
  path for existing runs; gold is optional.
- `scripts/recall_recovery/run_all_layers.py` - DB-observed ClinVar, PubTator,
  and figure recovery layers.
- `scripts/fetch_paywalled.py` - canonical paywall recovery entry point.
- `harvesting/browser_html/` - authenticated browser recovery strategies.
- `harvesting/browser_html/content_quality.py` - paywall/stub quality gate.
- `pipeline/table_router.py` and `pipeline/extraction.py` - clinical table
  extraction.
- `gene_variant_fetcher_gold_standard/` - normalized recall inputs and source
  exports.
- `corpus/<GENE>/<PMID>/` - consolidated, deduplicated home AND cross-run cache
  for ALL fetched source (full text + supplements + figures), one best usable
  copy per paper. Discover via `corpus/INDEX.json` / `corpus/INDEX.csv`
  (per-paper `full_text_status` ok|stub, `source_sha256`, figure/supplement
  counts). As of 2026-06-04 this replaces the old scattered per-run
  `pmc_fulltext/` trees (gitignored; ~11 GB / 6,382 papers). Run DBs +
  `extractions/` remain in their `results/` run dirs.
- `scripts/build_source_corpus.py` - builds/updates the corpus. Incremental +
  idempotent + never-downgrade: `--apply` folds new source in (adds new papers,
  upgrades only compromised categories), a clean re-run is a no-op. `--verify`
  checks `corpus/MANIFEST.sha256`. To move the corpus: `tar czf corpus.tgz
  corpus/` then `--verify` on the far end.
- Corpus-as-cache: `gvf-run` reuses usable cached source per PMID (skips
  re-fetch; a `stub`/compromised entry falls through to a fresh fetch so a new
  publisher key re-attempts it) and folds new fetches back via the
  `--corpus-sync` step (default on). Override the cache dir with
  `GVF_CORPUS_DIR`. Wired in `pipeline/steps.py` (`_consolidate_from_corpus`)
  and `cli/gvf_run.py` (`step_corpus_sync`).
- `cli/dashboard.py` / `gvf dashboard` - static HTML status/coverage/missingness/
  provenance dashboard from `corpus/INDEX.csv` + the scored DB(s). Per gene:
  source-coverage + extraction funnel + a provenance-completeness audit +
  "what's left". Per paper: an ADJUDICATION view with PubMed/DOI/PMC links and
  the exact on-disk full text rendered, where clicking an extracted record
  jumps to + highlights the `evidence_sentence`/`source_location` it came from
  (both 100%-populated). Also emits a per-gene `<GENE>_variants.html`
  variant-centric view (one row per unique variant → papers + affected/unaffected
  patient counts + clickthrough) for the eventual variant website. Read-only.
- Provenance fields (2026-06-04): the schema/migrator now capture
  `papers.first_author`, `variant_papers.count_provenance` (the "why" behind each
  carrier/affected count — was emitted by the prompt but previously dropped), and
  `individual_records.ethnicity`/`geographic_origin` (extracted by the Tier-3
  prompt, on by default). New columns auto-ALTER onto existing DBs.
- `scripts/backfill_paper_metadata.py` - fills `papers.{first_author,journal,
  publication_date,doi,pmc_id}` into existing DBs from local caches
  (`abstract_json/*.json` authors/journal/year + corpus `{pmid}_artifacts.json`
  doi/pmcid) — NO LLM, NO network, idempotent, COALESCE-safe. Run it after
  migrate/refresh; on the four cardiac genes it lifts first_author/journal/year
  from 0% to ~99%. ethnicity/"why" only populate on the NEXT extraction run.

## House Rules

- Use the project `.venv` when available.
- Do not commit agent scratch files: `HEARTBEAT.md`, `IDENTITY.md`, `SOUL.md`,
  `TOOLS.md`, `USER.md`, or `.openclaw/`.
- Do not commit `.env`, local `results/`, SQLite DBs, or generated
  `recall_metrics/`.
- Pre-commit hooks may reformat staged files. Never use `--no-verify`; re-stage
  after hooks run.
