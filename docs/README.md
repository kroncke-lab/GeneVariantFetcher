# GVF Documentation Map

Last reviewed: 2026-08-12.

Use this page to distinguish current operating instructions from dated evidence.
When two documents disagree, use the authority order below rather than choosing
the newest-looking filename.

## Current authorities

| Question | Authority |
| --- | --- |
| What should be done next? | [`../TASKS.md`](../TASKS.md) |
| What are the current measured metrics and caveats? | [`RECALL_STATUS.md`](RECALL_STATUS.md) |
| How does the current pipeline work? | [`ARCHITECTURE.md`](ARCHITECTURE.md) and [`EXTRACTION_CONTRACT.md`](EXTRACTION_CONTRACT.md) |
| How do I install and run it? | [`QUICKSTART.md`](QUICKSTART.md) |
| How do I refresh an existing run? | [`RECALL_REFRESH_RUNBOOK.md`](RECALL_REFRESH_RUNBOOK.md) |
| How do I run a new gene or move to another machine? | [`NEW_GENE_RUNBOOK.md`](NEW_GENE_RUNBOOK.md) and [`END_TO_END_RECALL_RUN.md`](END_TO_END_RECALL_RUN.md) |
| What files and schema are produced? | [`OUTPUT_FORMAT.md`](OUTPUT_FORMAT.md); executable schema in `harvesting/migrate_to_sqlite.py` |
| How does review-gold exchange work? | [`VARIANT_BROWSER_INTEGRATION.md`](VARIANT_BROWSER_INTEGRATION.md) |
| Which evaluation cohorts are active? | [`../benchmarks/evaluation_tiers/README.md`](../benchmarks/evaluation_tiers/README.md) |

`CLAUDE.md` is the canonical agent handoff. `AGENTS.md` and `CODEX.md` are
pointer files only. The handoff may summarize stable operating constraints, but
it must defer tasks and metrics to the authorities above.

## Append-only history

- [`RECALL_HISTORY.md`](RECALL_HISTORY.md) records benchmark and metric changes.
- [`PROTOCOL_CHANGELOG.md`](PROTOCOL_CHANGELOG.md) records protocol-affecting
  changes.
- [`../CHANGELOG.md`](../CHANGELOG.md) records release-level changes.

Append new entries. Do not rewrite older measurements to make them resemble the
current protocol.

## Dated or historical evidence

`PROTOCOL_COST_EVAL.md`, `RECALL_SOURCE_STRATIFIED.md`,
`SUPPLEMENT_ACQUISITION_PLAN.md`, `REPO_AUDIT_2026-08-12.md`, files under
`docs/reviews/`, and dated benchmark run directories describe the protocol and
data available when they were produced. Their measurements remain useful, but
their commands, model routes, checklists, and status statements are not the
current operating plan unless a current authority links to them for a specific
step.

The published dashboard is also historical until `RECALL_STATUS.md` says it was
regenerated from an accepted rescore.

## Maintenance rule

Keep one live checklist and one live metrics file. New current-facing Markdown
must link to the relevant authority instead of copying task lists or recall
tables. Label frozen cohorts, dated experiments, and superseded instructions at
the top of their own document.
