# GVF Tasks

## Current Focus
- Improve unique variant recall from 59% → 90%

## Completed
- [x] Springer API integration (full-text retrieval working)
- [x] Wire `gene_literature/supplements/` UnifiedSupplementFetcher into orchestrator
- [x] Fix all test failures and collection errors (18 → 0)
- [x] Canonical form comparison for variant matching

## Active Tasks
- [ ] Browser fetch paywalled papers from top-10 missing PMIDs
- [ ] Add Wiley supplement handler
- [ ] Integrate browser fallback for Karger (Cloudflare blocked)
- [ ] Submit Karger TDM request

## Backlog
- [ ] Expand validation to SCN5A, KCNQ1
- [ ] Create golden test sets for non-KCNH2 genes
- [ ] Re-run extraction with regex disabled to measure impact
- [ ] Consolidate duplicate documentation (ARCHITECTURE.md vs GVF_architecture.md)

See CLAUDE.md for current status and blockers.
