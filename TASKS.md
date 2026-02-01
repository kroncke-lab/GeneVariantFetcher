# GeneVariantFetcher - Tasks

## Goal
90% variant extraction recall (currently 5.4% - 73/1359 variants matched)

## Root Cause
Supplement scraper missing handlers for most publishers → 94.6% of variants missed

---

## Active Tasks

### Publisher Handlers
- [x] Audit supplement_scraper.py - documented all gaps
- [x] Springer/BMC handler implemented
- [x] Oxford Academic handler implemented
- [ ] Wiley handler
- [ ] Karger browser fallback (blocked by Cloudflare)
- [ ] Test all handlers with real extraction run
- [ ] Measure recall improvement

### Testing & Validation
- [ ] Re-run full extraction pipeline with new handlers
- [ ] Compare to Excel baseline (1359 variants)
- [ ] Target: 90%+ recall (1223+ variants matched)
- [ ] Document coverage by publisher

### Integration
- [ ] Merge `feature/save-raw-assets` branch
- [ ] Deploy updated pipeline
- [ ] Update documentation

---

## Publisher Coverage

| Publisher | Handler | Status |
|-----------|---------|--------|
| Nature | ✅ | Working |
| Elsevier | ✅ | Working |
| Springer/BMC | ✅ | Implemented, needs testing |
| Oxford Academic | ✅ | Implemented, needs testing |
| Karger | ❌ | Blocked by Cloudflare |
| Wiley | ❌ | Not started |

---

*Last updated: 2026-02-01 17:56*
