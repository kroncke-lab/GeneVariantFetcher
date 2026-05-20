# Codex / AI Agent Onboarding

## Active Main Checkout

The active main working directory is:

`/Users/kronckbm/GitRepos/GeneVariantFetcher`

Use that path for current GVF work unless Brett explicitly says otherwise.
Do not treat `.claude/worktrees/`, `.codex/worktrees/`, old `Projects/`, or
remote `/mnt/temp4/` copies as authoritative. They are historical side
worktrees or scratch checkouts and can lag behind `main`.

You're working on **GeneVariantFetcher (GVF)** — a pipeline that extracts
genetic variants from biomedical literature for the Kroncke Lab variant
interpretation toolkit. Goal: **90% unique-variant recall by June 2026**
(R01 grant submission). Current metrics, blockers, and next actions live in
`docs/RECALL_STATUS.md`; treat that file as the source of
truth and do not duplicate live recall tables here.

## Read these first (in this order)

1. **`docs/RECALL_STATUS.md`** — current measured baseline,
   blockers, stale/bloat audit, and next run plan. **Start here.**
2. **`CLAUDE.md`** — project architecture, paywall recovery context, and
   handoff notes. Historical metrics inside it defer to the current-status doc.
3. **`TASKS.md`** — what's done, what's active, what's blocked, what's
   backlog. The "Active Tasks" section is what's open right now.
4. `README.md`, `CONTRIBUTING.md`, `CHANGELOG.md` — standard project docs.

## Hot files you'll touch most often

- `cli/compare_variants.py` — gold-standard matcher. Has positional-digit
  guard, cDNA↔protein bridging, greedy 1-to-1 assignment.
- `harvesting/browser_html/` — Tier 3.5 authenticated Playwright recovery
  for paywalled publishers. `content_quality.py` is the quality gate that
  rejects abstract-only stubs.
- `scripts/fetch_paywalled.py` — paywall-recovery CLI. Cookie-based
  Playwright + per-domain pacing + PMC fallback. Canonical entry point
  for recovery work.
- `pipeline/` — Tier 1/2/3 extraction.

## House rules

- **Pre-commit hooks WILL reformat staged files** (ruff, ruff-format,
  trailing-whitespace). Expect to re-stage and recommit. **Never** use
  `--no-verify`; investigate hook output and fix root causes.
- The agent-scratch files (`HEARTBEAT.md`, `IDENTITY.md`, `SOUL.md`,
  `TOOLS.md`, `USER.md`, `.openclaw/`, `AGENTS.md`) are
  `.gitignore`d. Don't commit them.
- Worktrees under `.claude/worktrees/` are parallel investigation
  branches. **Most predate current `main`** — their `cli/compare_variants.py`
  and `harvesting/browser_html/strategies/*` are usually older than what's
  in main. Diff carefully before porting.
- The repo has a `.venv` at project root. `python -m pytest tests/`
  from project root works.
- This is a research codebase, not a product. Make small, well-tested
  changes. Don't add framework-level abstractions.

## Quick smoke tests after a change

```bash
# Unit tests for the matcher and the quality gate
.venv/bin/python -m pytest tests/unit/test_compare_variants.py -q
.venv/bin/python -m pytest tests/unit/ -q

# Smoke-test imports of the recovery pipeline
.venv/bin/python -c "from scripts.fetch_paywalled import europepmc_lookup_pmcid; \
    from harvesting.browser_html import validate_article_content; \
    print('ok')"
```

## When in doubt

- Brett (the user, lab PI) prefers minimal diffs and honest reports
  ("brutally honest" is a direct quote). Don't sugar-coat recovery rates.
- If you can't recover a paper, say so plainly and document why in
  `CLAUDE.md` / `TASKS.md` under blockers.
- The multi-gene measurement loop is working, but older recovery metrics used
  gold-PMID-conditioned enrichment layers. Current code defaults ClinVar and
  PubTator enrichment to DB-observed PMIDs only; use gold-PMID enrichment only
  as an explicitly labeled diagnostic.
- The current bottleneck is source/table coverage: missing supplements,
  paywalled publisher tables, browser challenges, and incomplete extraction of
  high-loss PMIDs such as KCNH2 PMID 15840476.
