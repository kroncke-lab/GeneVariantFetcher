# GeneVariantFetcher Previous Progress Report
**Generated:** 2026-02-05 16:29 CST

## Executive Summary

Multiple GVF runs have been completed for KCNH2, with the most comprehensive being the **Jan 28, 2026 run** which discovered 1,277 PMIDs and extracted 2,152 variants. A recall analysis on Feb 4 revealed **~6% recall** against the baseline Excel database, with the primary bottleneck being paper coverage (224 missing PMIDs from baseline).

---

## Output Directories Found

### 1. Primary GVF Output: `/mnt/temp2/kronckbm/gvf_output/`
**Status:** Active, contains multiple analysis runs

| Subdirectory | Description | Date |
|--------------|-------------|------|
| `KCNH2/` | Main run directory | Feb 4, 2026 |
| `KCNH2/20260202_173749/` | Timestamped run | Feb 2-3, 2026 |
| `KCNH2/browser_downloads/` | 72 PDFs from browser fetch | Feb 3, 2026 |
| `KCNH2/bulk_download/` | 269 PMC fulltext papers | Feb 3, 2026 |
| `KCNH2/golden_test_harvest/` | Golden test papers | Feb 3, 2026 |
| `KCNH2_clean/` | **161 extraction JSONs** (cleaned) | Feb 4, 2026 |
| `KCNH2_comparison_20260204/` | Comparison with baseline Excel | Feb 4, 2026 |
| `KCNH2_fresh.db` | Fresh SQLite database (786KB) | Feb 4, 2026 |
| `KCNH2_missing/` | 4 paywalled papers downloaded | Feb 4, 2026 |
| `KCNH2_top10_missing/` | Top 10 missing high-value PMIDs | Feb 5, 2026 |

**Key File:** `KCNH2_fresh.db` - Most recent consolidated database

### 2. GVF Repo Output: `/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/gvf_output/`
**Status:** Older comprehensive run (Jan 28-29, 2026)

| Subdirectory | Description |
|--------------|-------------|
| `KCNH2/20260128_210249/` | Complete workflow run |
| `browser_downloads/` | 474 PDFs |
| `abstract_only_papers.csv` | 85KB - papers without fulltext |
| `retry_timeout_pmids.csv` | PMIDs that timed out |

### 3. GVF New Output: `/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/output/`
**Status:** Fresh run in progress (Feb 5, 2026)

| Path | Description |
|------|-------------|
| `KCNH2/KCNH2/20260205_161223/` | Current active run |

---

## Run Statistics

### Run 1: Jan 28, 2026 (Most Comprehensive)
**Location:** `/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/gvf_output/KCNH2/20260128_210249/`
**Database:** `KCNH2.db` (958 KB)

| Metric | Value |
|--------|-------|
| PMIDs Discovered | 1,277 |
| PMIDs Passed Filters | 687 |
| Papers Downloaded | 415 |
| Papers Extracted | 590 |
| From Fulltext | 319 |
| From Abstract Only | 271 |
| Extraction Failures | 97 |
| **Total Variants Found** | **2,152** |
| Variants with Penetrance | 1,655 |
| Total Carriers | 11,958 |
| Affected Carriers | 8,215 |
| Success Rate | 46.2% |

**Extraction Files:** 592 JSON extractions in `extractions/` folder

### Run 2: Feb 5, 2026 (Current/Smaller Test)
**Location:** `/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/output/KCNH2/KCNH2/20260205_161223/`
**Database:** `KCNH2.db` (159 KB)

| Metric | Value |
|--------|-------|
| PMIDs Discovered | 100 |
| PMIDs Passed Filters | 55 |
| Papers Downloaded | 28 |
| Papers Extracted | 50 |
| **Total Variants Found** | **89** |
| Variants with Penetrance | 75 |
| Success Rate | 50.0% |

---

## Recall Analysis (Feb 4, 2026)

**File:** `/mnt/temp2/kronckbm/gvf_output/RECALL_ANALYSIS_20260204.md`

### Current Performance
- **Recall:** ~6% (60 matches / 991 baseline entries)
- **Target:** 90%
- **Primary Bottleneck:** Paper coverage, NOT extraction quality

### Coverage Gap
| Metric | Count |
|--------|-------|
| Excel Baseline PMIDs | 262 |
| SQLite Extracted PMIDs | 249 |
| **Overlapping PMIDs** | **37 (14.1%)** |
| **Missing PMIDs** | **224** |

### Within-Overlap Performance
For the 37 PMIDs both processed:
- Excel entries: 168
- SQLite variants: 213  
- Matched: 60
- **Within-overlap recall: ~36%**

### Top 10 Missing High-Value Papers
| PMID | Variant Count |
|------|---------------|
| 15840476 | 89 |
| 10973849 | 60 |
| 26496715 | 54 |
| 11854117 | 44 |
| 14661677 | 29 |
| 19038855 | 28 |
| 24667783 | 23 |
| 29650123 | 22 |
| 16922724 | 20 |
| 23631430 | 12 |

**These 10 papers alone = 381 variants (38% of missing entries)**

---

## Downloaded Papers Status

### Browser Downloads
| Location | Count |
|----------|-------|
| `/mnt/temp2/kronckbm/gvf_output/KCNH2/browser_downloads/` | 72 |
| `/mnt/temp2/kronckbm/gvf_output/KCNH2/bulk_download/` | 269 |
| `/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/gvf_output/browser_downloads/` | 474 |

### Download Failures (from Jan 29 analysis)
| Category | Count | Description |
|----------|-------|-------------|
| CAPTCHA Blocked | 7 | Cloudflare detection |
| Crash/Timeout | 34 | Retryable |
| Paywall | 108 | Vanderbilt lacks access |
| **Total Failed** | **149** | |

---

## Partial Results to Resume

### 1. Cleaned Extractions Ready for Analysis
**Path:** `/mnt/temp2/kronckbm/gvf_output/KCNH2_clean/`
**Count:** 161 JSON extraction files
**Status:** ✅ Ready to merge/analyze

### 2. Top 10 Missing Papers (Downloaded)
**Path:** `/mnt/temp2/kronckbm/gvf_output/KCNH2_top10_missing/`
**Papers Downloaded:**
- 19038855 (28 variants expected)
- 24606995 
- 24667783 (23 variants expected)
- 29622001
**Status:** ✅ Downloaded, awaiting extraction

### 3. Missing Baseline PMIDs List
**Path:** `/mnt/temp2/kronckbm/gvf_output/missing_baseline_pmids.txt`
**Count:** 224 PMIDs
**Status:** 📋 Ready to download

### 4. Comparison Results
**Path:** `/mnt/temp2/kronckbm/gvf_output/KCNH2_comparison_20260204/`
**Files:**
- `discrepancies.csv` - 90KB detailed discrepancies
- `missing_in_sqlite.csv` - 62KB (variants GVF didn't find)
- `missing_in_excel.csv` - 22KB (variants GVF found not in baseline)
- `report.md` - Human-readable summary

---

## Recommended Next Steps

### Phase 1: Paper Coverage (High Impact)
1. Download the 224 missing baseline PMIDs
2. Prioritize top 10 high-value papers (already downloaded in `KCNH2_top10_missing/`)
3. Expected gain: +832 potential variant entries

### Phase 2: Extraction Quality
1. Re-enable table regex extraction with gene validation
2. Improve variant normalization
3. Expected gain: 70-75% recall

### Phase 3: Validation
1. Run iterative comparisons after each batch
2. Analyze false negatives in detail
3. Target: 90%+ recall

---

## Key Database Locations

| Description | Path |
|-------------|------|
| Most Recent SQLite | `/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db` |
| Jan 28 Complete Run | `/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/gvf_output/KCNH2/20260128_210249/KCNH2.db` |
| Feb 5 Test Run | `/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/output/KCNH2/KCNH2/20260205_161223/KCNH2.db` |

---

*Report generated by GVF subagent search task*
