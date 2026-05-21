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
grant target is 90% unique-variant recall by June 2026.

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
- Turnkey driver: `gvf gvf-run <GENE> --email <email> --output <dir>` runs
  doctor -> extraction -> DB-observed recovery layers -> report. It backs up the
  DB before recovery-layer mutation. `--with-v12` is KCNH2-only and opt-in.
- Matcher: `cli/compare_variants.py` has positional-digit guards, greedy
  one-to-one assignment, cDNA/protein bridging, and 2026-05-15 fixes for
  frameshift spellings such as `fsTer`, `fs/185`, `fs+*49`, and malformed
  `AlaX14`.
- Gold-standard builders: `scripts/build_gold_standard_from_varbrowser.py` and
  `scripts/build_ryr2_gold_standard_from_xlsx.py`.

## Manual Acquisition Queue

The manual acquisition queue is the set of passed-filter PMIDs that do not have
usable assembled full text after the automatic harvesters and fallback routes.
They are not "assigned to a human" by code; they are papers GVF could not fetch
with the currently available credentials and network position.

Current high-yield missing PMIDs are listed in
`docs/RECALL_STATUS.md`. Keep PMID rankings there so this
handoff file does not drift.

## Highest ROI Blocker

**As of 2026-05-21, the Elsevier insttoken unblock has landed.** Vanderbilt's
institutional `X-ELS-Insttoken` is installed in `.env` and 242 of 246
previously-paywalled Elsevier articles across KCNH2/SCN5A/RYR2/KCNE1/KCNQ1 now
return full text via `harvesting/elsevier_api.py`. The unlocked bodies were
saved as `{PMID}_FULL_CONTEXT.md` into each run's existing `pmc_fulltext/`
directory; pre-token stub files were preserved as `*.pre_insttoken_bak`.

The new ROI blocker is **re-extraction**: the saved full-text files do not
affect PMID recall until the SQLite DBs are rebuilt against them. After
re-extraction, expect PMID recall to lift substantially (currently
KCNH2 87.8%, KCNQ1 80.3%, RYR2 68.0%, SCN5A 70.9%, aggregate 75.4%; target 90%).

Residual non-Elsevier paywalls (Wiley revoked key, Karger Cloudflare,
Sage CF fingerprint) remain in `TASKS.md` Blocked.

## Run On Another Computer Or VPN

From a clean checkout:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pip install -r gui/requirements.txt
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
.venv/bin/python scripts/fetch_paywalled.py \
  --input results/KCNH2/20260506_102238/pmc_fulltext/paywalled_missing.csv \
  --output results/KCNH2/20260506_102238/pmc_fulltext \
  --no-headless
```

Then rerun extraction/migration and score:

```bash
.venv/bin/python scripts/run_recall_suite.py --score --genes KCNH2 \
  --db KCNH2=results/KCNH2/20260506_102238/end_to_end_20260515_manual_recovery/KCNH2_v12_manual_recovery_20260515.db \
  --outdir recall_metrics/kcnh2_after_vpn_recovery
```

## Active Work

Current active work is tracked in
`docs/RECALL_STATUS.md`. Do not duplicate metric gaps here.
The short version (as of 2026-05-21): Elsevier source acquisition is unblocked;
the next-highest leverage step is re-extracting KCNH2/KCNQ1/RYR2/SCN5A against
the consolidated `pmc_fulltext/` directories so the unlocked bodies become DB
rows. Beyond that: fix general count/affected-status extraction issues, address
the remaining Wiley/Karger/Sage paywalls, and keep no-gold QC prominent for new
genes.

## Files To Know

- `cli/compare_variants.py` - gold-standard matcher and recall summary.
- `scripts/run_recall_suite.py` - multi-gene recall scoring.
- `scripts/fetch_paywalled.py` - canonical paywall recovery entry point.
- `harvesting/browser_html/` - authenticated browser recovery strategies.
- `harvesting/browser_html/content_quality.py` - paywall/stub quality gate.
- `pipeline/table_router.py` and `pipeline/extraction.py` - clinical table
  extraction.
- `gene_variant_fetcher_gold_standard/` - normalized recall inputs and source
  exports.

## House Rules

- Use the project `.venv` when available.
- Do not commit agent scratch files: `HEARTBEAT.md`, `IDENTITY.md`, `SOUL.md`,
  `TOOLS.md`, `USER.md`, or `.openclaw/`.
- Do not commit `.env`, local `results/`, SQLite DBs, or generated
  `recall_metrics/`.
- Pre-commit hooks may reformat staged files. Never use `--no-verify`; re-stage
  after hooks run.
