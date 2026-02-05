# GeneVariantFetcher - Tasks

## Goal
90% variant extraction recall

## Current Status (2026-02-04 11:30)

### Recall Analysis Complete
| Metric | Value |
|--------|-------|
| Overall recall | **6%** (60/991) |
| Baseline PMIDs | 262 |
| Processed PMIDs | 249 |
| **Overlapping PMIDs** | **37 (14.1%)** |
| Missing PMIDs | 224 |

### Key Finding
**The bottleneck is PAPER COVERAGE, not extraction quality.**
- 224 baseline PMIDs (832 variant entries) were never downloaded
- For overlapping 37 PMIDs, recall is ~36% (extraction works)

### Files Created
- `/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db` - SQLite with 2233 variants
- `/mnt/temp2/kronckbm/gvf_output/missing_baseline_pmids.txt` - 224 PMIDs to download
- `/mnt/temp2/kronckbm/gvf_output/RECALL_ANALYSIS_20260204.md` - Full analysis

**Note**: GVF extracts MORE variants than Excel baseline for papers it processes, but needs broader paper coverage.

---

## Completed Tasks ✅

### Publisher Handlers
- [x] Audit supplement_scraper.py - documented all gaps
- [x] Springer/BMC handler implemented (2026-02-01)
- [x] Oxford Academic handler implemented (2026-02-01)
- [x] Wiley handler implemented (2026-02-01)
- [x] Wiley API integration (2026-02-03)
- [x] Elsevier API integration (2026-02-03)
- [x] Supplement reference parser created
- [x] Post-download gap validation integrated (2026-02-04)

### Testing & Validation
- [x] Full extraction run with synonyms (749 papers downloaded)
- [x] Golden test set created (13 papers)
- [x] Extraction JSON files generated (249 papers)

---

## Remaining Tasks

### HIGH PRIORITY - Paper Coverage Gap (critical path to 90%)
- [ ] **Download 224 missing baseline PMIDs** - See `missing_baseline_pmids.txt`
  - Top 10 papers cover 381 variants (38% of gap)
  - PMID 15840476 alone has 89 variants!
- [ ] **Re-run extraction** on newly downloaded papers
- [ ] **Re-enable table regex with gene validation** - Currently DISABLED

### Medium Priority
- [ ] **Karger browser fallback** - Blocked by Cloudflare, need TDM request
- [ ] **Springer API** - Brett needs to register with Vanderbilt credentials
- [x] **Aggregate extractions to SQLite** - DONE: `/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db`
- [x] **Run formal recall comparison** - DONE: 6% recall, 14% PMID coverage

### Medium Priority
- [ ] Merge `feature/save-raw-assets` branch to main
- [ ] Update documentation with new handlers
- [ ] Browser fetch remaining paywalled papers (10 in golden test)

### Low Priority
- [ ] Extend to additional genes (SCN5A, KCNQ1, etc.)
- [ ] Create gene synonyms lookup file

---

## Publisher Coverage

| Publisher | Handler | API | Status |
|-----------|---------|-----|--------|
| Nature | ✅ | - | Working |
| Elsevier | ✅ | ✅ | Working |
| Springer/BMC | ✅ | ❌ | Handler working, API needs registration |
| Oxford Academic | ✅ | - | Working |
| Wiley | ✅ | ✅ | Working |
| Karger | ❌ | - | Blocked by Cloudflare |

---

## API Access Status (Vanderbilt)
- **Elsevier**: ✅ API key configured and tested
- **Wiley**: ✅ API key configured and tested  
- **Springer**: ❌ Brett needs to register
- **Karger**: ❌ No institutional agreement, TDM request drafted

---

*Last updated: 2026-02-04 08:30*
