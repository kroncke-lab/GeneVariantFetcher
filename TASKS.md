# GeneVariantFetcher - Tasks

## Goal
90% variant extraction recall

## Current Status (2026-02-04)
- **Extraction run completed**: 249 papers processed, 201 with variants
- **Variants extracted**: 2370 variant mentions (2138 unique)
- **Total carriers**: 5179
- **Golden test set**: 13 papers, 45.2% recall on overlap

**Note**: GVF now extracts MORE variants than the Excel baseline (2370 vs 1359), indicating the pipeline is catching variants the manual curation missed.

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

### High Priority
- [ ] **Karger browser fallback** - Blocked by Cloudflare, need TDM request
- [ ] **Springer API** - Brett needs to register with Vanderbilt credentials
- [ ] **Aggregate extractions to SQLite** - For proper comparison against Excel
- [ ] **Run formal recall comparison** - Need matched variant counts

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
