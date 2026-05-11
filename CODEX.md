# Codex / AI Agent Onboarding

You're working on **GeneVariantFetcher (GVF)** — a pipeline that extracts
genetic variants from biomedical literature for the Kroncke Lab variant
interpretation toolkit. Goal: **90% unique-variant recall by June 2026**
(R01 grant submission). Current measured baseline: **59.1%** (2025-12-11).

## Read these first (in this order)

1. **`CLAUDE.md`** — comprehensive project status, architecture, paywall
   recovery pipeline, the Elsevier INSTTOKEN unblock path, publisher
   coverage matrix, next steps. Has an "Notes for AI Agents Picking This
   Up" section at the bottom. **Start here.**
2. **`TASKS.md`** — what's done, what's active, what's blocked, what's
   backlog. The "Active Tasks" section is what's open right now.
3. `README.md`, `CONTRIBUTING.md`, `CHANGELOG.md` — standard project docs.

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
- Re-running KCNH2 extraction end-to-end (and regenerating
  `comparison_results/`) is the single most important unblocking action.
  Almost any other task is downstream of that measurement loop.
