# GVF Tasks

## Current Focus
- Improve unique variant recall from 59% → 90%

## Completed
- [x] Springer API integration (full-text retrieval working)
- [x] Wire `gene_literature/supplements/` UnifiedSupplementFetcher into orchestrator
- [x] Fix all test failures and collection errors (18 → 0)
- [x] Canonical form comparison for variant matching
- [x] Remove 50-paper download cap (was silently dropping 350+ papers)
- [x] Fail-open Tier 1 filter (min_keywords 2→1, no-abstract papers pass through)
- [x] Add source-completeness report (Step 3.5: `source_completeness.json`)
- [x] Add zero-variant QA flagging and re-extraction (GVF_QA_MODEL env var)
- [x] Improve cohort/screening table carrier count extraction prompts
- [x] Fix pytest environment (281 unit tests passing)

## Active Tasks
- [ ] **Re-run KCNH2 extraction end-to-end** with patched pipeline; regenerate `comparison_results/`
- [ ] Browser fetch paywalled papers from top-10 missing PMIDs
- [ ] Add Wiley supplement handler
- [ ] Integrate browser fallback for Karger (Cloudflare blocked)
- [ ] Submit Karger TDM request
- [ ] Set weekly recall cadence (Friday re-run + compare)

## Backlog
- [ ] Expand validation to SCN5A, KCNQ1
- [ ] Create golden test sets for non-KCNH2 genes
- [ ] Re-run extraction with regex disabled to measure impact
- [ ] Consolidate duplicate documentation (ARCHITECTURE.md vs GVF_architecture.md)

See CLAUDE.md for current status and blockers.
