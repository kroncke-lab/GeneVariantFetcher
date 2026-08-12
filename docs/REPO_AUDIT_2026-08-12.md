# Repository audit — 2026-08-12

## Decision

Current `main` was healthy before the audit: the offline unit, integration, and
negative-case suites passed. The two open rescue branches were both stale and
conflicted, and neither was safe to merge wholesale. This audit selectively
ported the changes that were independently reviewable and testable, retired one
unsafe count mutation, and preserved negative experiments as historical
evidence.

The review used seven independent reasoning passes:

- Codex GPT-5.6 at maximum effort, including separate correctness, cleanup, and
  recall/benchmark reviews;
- Grok 4.5 at `high`, the highest effort accepted by the installed CLI (`max`
  was rejected);
- Gemini 3.1 Pro High through AGY; and
- Claude Opus 4.6 Thinking through AGY. The native Claude CLI was not signed in.

The reviewers converged on the branch disposition and the highest-risk findings
below.

## Integrated improvements

1. **Preserve raw scientific observations.** Retired
   `refuse_all_unaffected`, which deleted raw affected/unaffected counts based on
   a broad quarantine state. Quarantine is not evidence that a genuine negative
   cohort observation is false. Step 3.45 now performs only context-gated,
   null-only figure-count adoption; the trust projection handles suspicion.
2. **Keep explicit zero distinct from missing.** The comparison scorer now
   carries count-presence markers through Excel and SQLite aggregation. An
   asserted zero participates in error scoring; a null or unparseable value does
   not silently become zero.
3. **Make source-reachability reporting honest.** ALL GOLD remains the primary
   turnkey metric. SOURCE-REACHABLE is a secondary reader diagnostic, PDF text
   and all image locations are inventoried, and missing predictions are charged
   in the coverage-aware loss. Calibration is labeled as conditional on rows
   already matched by a paper layer, not a false-exclusion bound.
4. **Repair benchmark claims.** Native schema-v2 predictions may truthfully
   report exact zero tokens, while missing legacy telemetry renders as
   unavailable. Reports no longer invent zero duration, locked traces, or a
   strict no-prior-gold-read claim for external schema-v1 projections. The
   locked 2026-07-26 run has an append-only erratum rather than rewritten
   artifacts.
5. **Ship runtime data.** Wheels now contain the 4,000-plus KCNH2 alias mapping
   and reference FASTA/metadata assets. Installed-layout tests exercise both;
   package discovery is restricted to real packages so builds do not scan the
   external corpus.
6. **Widen provenance.** The extraction digest and routing snapshot now cover
   settings, orchestration, recovery/repair, guards, trust/final checks,
   migration, figure readers, and other behavior-changing modules.
7. **Prevent cross-gene source misfiling.** An explicit `--assume-gene` wins
   over corpus PMID reuse, and implicit PMID reuse occurs only when the existing
   gene association is unambiguous.
8. **Remove/document stale behavior.** Removed the advertised no-op manifest
   cleanup command, aligned Ruff 0.15.20 across local and CI configuration,
   moved GitHub's official checkout/setup actions to their Node-24 releases,
   corrected the parked final-check default and current orchestrator docs, and
   redirected gene-synonym documentation to the runtime registry.
9. **Keep rescoring hermetic when requested.** `run_recall_suite.py` can skip
   source-backed disagreement artifacts while preserving aggregate scoring,
   avoiding accidental corpus reads during a metric-only replay.

## Rescue branch disposition

### `rescue/count-semantics-wip-20260810` / PR #181

The unique commit is not a count-semantics implementation; it is a Tier-2
relevance shadow harness. The harness, pinned cohort, test, and locked summary
were salvaged as a diagnostic. Its labels are derived from historical pipeline
decisions rather than manual relevance gold, so the result does not justify a
production model switch. The remaining stale documentation should not be
merged.

### `rescue/main-wip-audit-and-reason-routing` / PR #179

Benchmark-report honesty, expanded provenance, and the two durable run summaries
were salvaged. Reason-class routing was rejected: its own locked replay reduced
verifier cost but worsened carrier MAE from 0.723 to 0.902, failing the declared
gate. Its broad, older settings/extraction/docs changes should not be merged.

## Intentionally unmerged or still open

- Figure-count recovery cannot yet enrich a variant-paper link that the text
  extractor already created, and adopted counts lack the stronger structured
  role/locator and pending-quarantine provenance used by count recovery. This is
  the highest-priority count-path follow-up.
- Claude's dirty `unassessed_count` worktree is preserved as WIP. Its schema,
  prompt, migration, and arithmetic changes pass focused tests, but aggregation,
  adjudicator schemas, final-check/trust projections, carrier-guard handling,
  provenance/rule versioning, and unassessed-only migration are incomplete.
  Merging it now could silently drop the new partition or trust an impossible
  one.
- Installed `gvf` still guards optional imports from unshipped `scripts.*` for
  metadata backfill, dashboard trust readers, institutional paywall preflight,
  and EZproxy self-heal. Move the reusable code into packaged modules before
  claiming installed-feature parity.
- Source acquisition remains the dominant recall opportunity. Current error
  inventory is much larger for missing/stub, abstract-only, and table-body
  source gaps than for count semantics. Do not trade source recovery for more
  speculative count mutation.
- Re-score the active evaluation tiers before publishing new headline count
  metrics; older reports may have treated explicit gold zero as missing.

The active forward work is recorded in `TASKS.md`; current measurements and
their caveats remain in `docs/RECALL_STATUS.md`.

## Validation

- Offline unit suite: **1,835 passed, 1 skipped**.
- Offline integration suite: **52 passed**.
- Curated negative/precision guards: **3 passed**.
- Ruff 0.15.20 lint and format check: clean across **407 files**.
- Cold-start harness: LDLR dry-run passed without touching live services.
- Wheel build: succeeded; the archive contains the KCNH2 alias dictionary and
  reference-sequence assets, and installed-layout runtime checks passed.
