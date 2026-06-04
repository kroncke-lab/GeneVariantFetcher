# AGENTS.md

This file exists so tools that look for `AGENTS.md` find the right entry point.
It intentionally carries **no independent content** — it mirrors the pointers in
`CLAUDE.md`, which is the single canonical agent handoff for this repo. Keep them
in sync by only ever pointing here, never duplicating status or metrics.

## Read first

- **`CLAUDE.md`** — canonical handoff: active checkout, architecture, recovery
  pipeline, house rules, and where everything lives. Read this before doing
  anything substantive.
- **`docs/RECALL_STATUS.md`** — the single source of truth for live recall
  metrics, the scored baseline, highest-yield missing PMIDs, and the next-run
  plan. No other doc (including this one) should restate a recall number.
- **`docs/NEW_GENE_RUNBOOK.md`** — the no-gold operational flow for a new
  gene-disease pair.
- **`docs/API_KEYS.md`** — which keys are required vs. optional.

## Turnkey command

```bash
gvf gvf-run <GENE> --email you@example.com --output ./results [--disease "<phenotype>"]
```

Source recovery (paywall + supplement acquisition) runs **by default**; pass
`--no-source-recovery` for a fast PMC/free-text-only pass.

## House rules (summary — full list in `CLAUDE.md`)

- Use the project `.venv`.
- Do not commit `.env`, local `results/`, SQLite DBs, generated
  `recall_metrics/`, or agent scratch files (`HEARTBEAT.md`, `IDENTITY.md`,
  `SOUL.md`, `TOOLS.md`, `USER.md`, `.openclaw/`).
- Never use `--no-verify`; re-stage after pre-commit hooks run.
- Side worktrees/checkouts are for experiments; merge useful work into the
  active checkout and push to `origin/main` before calling it current.
