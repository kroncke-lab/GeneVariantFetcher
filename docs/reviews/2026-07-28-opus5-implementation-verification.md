# Claude Opus 5 implementation verification — 2026-07-28

Scope: implementation and independent verification of the findings in
`2026-07-28-opus5-code-review.md`. Base HEAD `037a609`; the existing
dirty worktree was preserved. Nothing was staged, committed, pushed, or sent to
a live provider by the verification pass.

## Model provenance

- Review: Claude Code `2.1.218`, explicit `claude-opus-5`, effort `max`,
  canonical model `claude-opus-5`, session
  `b94a088d-1ae2-42c1-a537-a8510d915822`.
- Separate implementation/correction task: explicit `claude-opus-5`, effort
  `max`, canonical model `claude-opus-5`, session
  `cf7ee108-6b34-4ffb-9520-a71671445b85`.
- The implementation task received the full review Markdown as context. Its
  focused correction pass accepted all eight independent objections rather
  than claiming route completeness.

## Outcome

The original installed-package load defect is fixed for the three material
paths identified by the audit:

1. `utils/llm_trace_html.py` now contains the packaged HTML builder; the
   `scripts/build_llm_trace_html.py` file is only a CLI wrapper.
2. `pipeline/count_recovery.py` now contains the packaged recovery/resolver
   implementation; `scripts/recover_counts.py` is only a CLI wrapper.
3. `pipeline/adjudication_contract.py` now contains the overlay contract used by
   `cli/compare_variants.py`, removing its three unguarded `scripts.*` imports.

The exact load error was `pyproject.toml` shipping only `pipeline*`,
`harvesting*`, `config*`, `utils*`, `gene_literature*`, and `cli*` while those
runtime modules imported the excluded `scripts` package. A checkout/editable
install masked the error. The new isolated installed-layout test copies only
the shipped packages and proves that report generation, count recovery imports,
and adjudication-overlay loading work with no `scripts/` directory.

The trace feature now records raw requests/responses, scoped decision events,
monotonic attempts, accepted/discarded/failed links, multi-call accepted-link
arrays, provider-visible reasoning summaries and token telemetry (without
claiming hidden chain-of-thought), write-time integrity levels, and per-run
storage isolation. The report is self-contained, bounded, sharded for large
runs, searchable/filterable, responsive, dark-mode capable, and null/shard-safe
for jump and record links.

## Independent verification

- Focused trace/packaging/trust subset: **248 passed in 117.00s**.
- Full offline unit suite: **1644 passed, 1 skipped in 152.66s**.
- `ruff check .`: clean.
- `git diff --check`: clean.
- Generated a two-paper synthetic report with accepted calls, decision
  rationales, long prompts, source selection, reasoning telemetry, and a failed
  provider call.
- Extracted JavaScript from that generated report: `node --check` passed.
- Parsed the static document skeleton: no duplicate IDs.
- Direct visual rendering was not performed: the available browser rejected the
  local `file://` report under its URL security policy, and that restriction was
  not bypassed. DOM/interaction contract tests cover search, filtering,
  accepted-call links, omitted/sharded records, relative source links, theme
  persistence guards, clipboard fallback, and script-safe data embedding.

## Remaining punch list

### P0 — operational security

- The ignored `.env.swp` remains on disk and has credential-shaped content
  documented by the prior audit. Close the owning editor, remove the swap file,
  rotate every represented credential, and exclude swap files from snapshots
  and archives. This task deliberately did not read, delete, or expose it.

### P1 — measurement before enablement

- Count recovery remains default OFF. Recreate clean database copies from the
  locked baseline and remeasure hardened v2 through the real validation/lock
  harness. Enable it only if recall improves without unacceptable MAE or
  attribution regression, with explicit attention to RYR2.

### P1 — route/configuration completion

- Standalone `cli/extract.py`, `cli/discover.py`,
  `scripts/targeted_land.py`, `scripts/extract_figure_variants.py`, and
  `scripts/recall_audit/*` do not configure tracing by default.
- `scripts/smoke_azure_models.py` has a raw call but no scoped normalized
  decision event.
- Remaining guarded shipped-package `scripts.*` imports do not hard-crash, but
  installed-package metadata backfill, dashboard trust summaries, and
  institutional preflight can degrade. Move those helpers into shipped modules.

### P2 — proof depth and human QA

- Add route-specific accepted-link tests for the six single-call routes that
  currently rely on shared derivation tests: clinical triage,
  extraction-priority triage, claim verification, final check,
  source-grounded summary, and synonym relevance.
- Perform desktop/mobile interaction and visual QA on a real completed run in
  an approved browser surface. The implementation is structurally and
  syntactically verified, but this session could not render a local file.

These open items are also recorded in the canonical forward checklist,
`TASKS.md`.
