# Supplement Scraper Audit

**Date:** 2026-01-31
**Auditor:** Boswell

## Current Publisher Coverage

| Publisher | Domain(s) | Handler | Status |
|-----------|-----------|---------|--------|
| Nature | nature.com | `scrape_nature_supplements()` | ✅ Working |
| Elsevier | sciencedirect.com, elsevier.com | `scrape_elsevier_supplements()` | ✅ Working |
| PMC | ncbi.nlm.nih.gov/pmc | Generic + URL variants | ✅ Working |
| Wiley | wiley.com | Fulltext only | ⚠️ No supplement handler |
| Karger | karger.com | None | ❌ **MISSING** |
| Oxford Academic | academic.oup.com | None | ❌ Missing |
| Springer | springer.com, link.springer.com | None | ❌ Missing |
| AHA Journals | ahajournals.org | None | ❌ Missing |
| PNAS | pnas.org | None | ❌ Missing |
| BMJ | bmj.com | None | ❌ Missing |
| Frontiers | frontiersin.org | None | ❌ Missing |

## Critical Gap: Karger

**Impact:** High - Karger publishes major cardiology journals including *Cardiology*.

**Example DOI:** 10.1159/000440608 (PMID 26496715)
- This paper has 54 KCNH2 variants in "online suppl. table 1"
- GVF downloaded main paper PDF but NOT the supplementary table
- Supplement URL pattern: `karger.com/doi/suppl/10.1159/XXXXXX`

**Current behavior:**
1. DOI resolver follows redirect to article page
2. Generic scraper searches page for supplement links
3. Karger hides supplements behind a separate URL (not linked on main page)
4. Result: Supplements not found

## Generic Scraper Analysis

**File:** `harvesting/supplement_scraper.py` → `scrape_generic_supplements()`

**What it finds:**
- Links with keywords: "supplement", "supporting", "appendix", "additional file"
- Links ending in file extensions: .pdf, .docx, .xlsx, .csv, .zip, etc.

**What it misses:**
- Publisher-specific supplement URLs not on the page
- Supplements behind JavaScript loaders
- Supplements requiring navigation to separate pages
- References in text like "online suppl. table 1" (doesn't construct URLs)

## Recommended Fixes (Priority Order)

### P0: Add Karger handler
```python
def scrape_karger_supplements(self, html: str, base_url: str) -> List[Dict]:
    """
    Karger supplement URL pattern:
    - Article: karger.com/crd/article/133/2/73/...
    - Supplements: karger.com/doi/suppl/10.1159/XXXXXX
    """
    # Extract DOI from page
    # Construct supplement URL
    # Fetch supplement page
    # Extract file links
```

### P1: Add Wiley supplement handler
Currently has fulltext extraction but no supplement detection.

### P2: Add Oxford Academic handler
Common for genetics papers.

### P3: Add AHA Journals handler
Critical for cardiac genetics (Circulation, JAHA, etc.)

## URL Patterns to Implement

| Publisher | Article URL | Supplement URL |
|-----------|-------------|----------------|
| Karger | `karger.com/crd/article/...` | `karger.com/doi/suppl/10.1159/XXXXXX` |
| Wiley | `onlinelibrary.wiley.com/doi/...` | `onlinelibrary.wiley.com/action/downloadSupplement?doi=...` |
| Oxford | `academic.oup.com/...` | Same page, `#supplementary-data` section |
| AHA | `ahajournals.org/doi/...` | `ahajournals.org/doi/suppl/...` |

## Code Locations

- Main scraper: `harvesting/supplement_scraper.py`
- DOI resolution: `harvesting/doi_resolver.py`
- Orchestration: `harvesting/orchestrator.py`
- Download logic: `harvesting/orchestrator.py` → `_process_supplements()`

## Test Cases Needed

| PMID | Publisher | Expected Supplements | Current Result |
|------|-----------|---------------------|----------------|
| 26496715 | Karger | 1+ (suppl table with 54 variants) | 0 |
| TBD | Wiley | TBD | TBD |
| TBD | Oxford | TBD | TBD |

---

*Audit complete. Karger handler is the critical missing piece for KCNH2 recall improvement.*
