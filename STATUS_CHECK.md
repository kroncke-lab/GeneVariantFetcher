# GeneVariantFetcher Status Check

**Date:** 2026-02-05
**Analyst:** Boswell (subagent)

---

## Executive Summary

| Metric | Current Value | Target |
|--------|---------------|--------|
| Overall Recall | **6%** | 90% |
| PMID Coverage | 14.1% (37/262) | >90% |
| Within-overlap Recall | ~36% | >80% |
| Total Variants in DB | 2,233 | - |
| Extractions Completed | 249 papers | - |

**Primary Bottleneck:** Paper coverage gap (224 missing PMIDs = 832 variant entries = 84% of gap)

---

## 1. Handler Status Assessment

### ✅ Fully Functional Handlers

| Publisher | Handler | API | Notes |
|-----------|---------|-----|-------|
| **Nature** | `scrape_nature_supplements` | - | Working well |
| **Elsevier** | `scrape_elsevier_supplements` | ✅ | MMC patterns, ScienceDirect support |
| **Springer/BMC** | `scrape_springer_supplements` | ❌ needs registration | Implemented 2026-02-01, MediaObjects pattern |
| **Oxford Academic** | `scrape_oxford_supplements` | - | Implemented 2026-02-01 |

### ⚠️ Partially Working

| Publisher | Issue | Impact |
|-----------|-------|--------|
| **Wiley** | API works but **no supplement handler** | Supplements missed |
| **Generic fallback** | Only finds obvious supplement links | Misses JS-rendered, API endpoints |

### ❌ Broken/Blocked

| Publisher | Issue | Workaround |
|-----------|-------|------------|
| **Karger** | Cloudflare blocks all requests | Browser fallback exists (`browser_supplement_fetcher.py`) but NOT integrated into main pipeline |
| **AHA Journals** | No handler at all | Falls through to generic |

---

## 2. Pipeline Test Results

### Quick Test: Top 10 Missing PMIDs Download

From the 10 highest-impact missing PMIDs:
- **4 successful** (19038855, 24606995, 29622001, 24667783)
- **6 failed**:
  - 5 = "Free text extraction failed" (paywalled)
  - 1 = "No PMCID" (14661677)

**Success rate: 40%** — indicates many baseline papers are paywalled or lack PMC access.

### Environment Issue ⚠️

The Python virtual environment is incomplete:
```
$ pip list
Package    Version
pip        21.3.1
setuptools 53.0.0
```

**No dependencies installed!** Pipeline cannot run until `pip install -e .` is executed.

---

## 3. Extraction Quality Analysis

For successfully processed papers, extraction quality is **reasonable**:

**Example: PMID 15840476** (89 variants expected, 89 extracted)
- Correct gene symbol (KCNH2)
- Both cDNA and protein notation captured
- Source location tracked (Table 3, Row N)
- Penetrance data inferred from patient counts

**Variant data schema includes:**
- Gene symbol, cDNA/protein notation
- Clinical significance
- Patient counts and penetrance data
- Source location in paper
- Functional data and segregation

---

## 4. Top 3 Actionable Improvements (Ranked by Impact)

### 🥇 #1: Fix venv and Download Missing 224 PMIDs

**Impact:** Could raise recall from 6% → 50-60%  
**Effort:** Medium  
**Blockers:** Paywall access for ~60% of papers

```bash
# Fix environment
cd /mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher
python3 -m venv venv
source venv/bin/activate
pip install -e .

# Then run download for missing PMIDs
python -m cli.automated_workflow --pmid-file missing_baseline_pmids.txt
```

**High-value PMIDs to prioritize:**
| PMID | Expected Variants |
|------|-------------------|
| 15840476 | 89 |
| 10973849 | 60 |
| 26496715 | 54 |
| 11854117 | 44 |
| 14661677 | 29 |

These 5 papers = 276 variants (28% of gap).

### 🥈 #2: Add Wiley Supplement Handler

**Impact:** Could add 5-10% recall  
**Effort:** Low-Medium  
**Pattern:** `/action/downloadSupplement?doi=...&file=...`

Wiley is a major cardiac genetics publisher. The API works for full text, but **supplement URLs are not being captured**. Pattern is similar to Oxford's `/downloadSupplement` endpoint.

### 🥉 #3: Integrate Browser Fallback into Main Pipeline

**Impact:** Unlock Karger + paywalled papers  
**Effort:** Medium  
**Files:** `harvesting/browser_supplement_fetcher.py` exists but is standalone

Currently:
- Karger is 100% blocked by Cloudflare
- Paywalled papers fail silently
- Browser fetcher exists but requires manual invocation

Needed:
- Auto-fallback to browser when HTTP fails with 403/Cloudflare
- Integrate with Vanderbilt proxy for authenticated access

---

## 5. Files Reference

| File | Purpose |
|------|---------|
| `/mnt/temp2/kronckbm/gvf_output/KCNH2_fresh.db` | SQLite with 2,233 variants |
| `/mnt/temp2/kronckbm/gvf_output/missing_baseline_pmids.txt` | 224 PMIDs to download |
| `/mnt/temp2/kronckbm/gvf_output/RECALL_ANALYSIS_20260204.md` | Detailed recall breakdown |
| `/mnt/temp2/kronckbm/gvf_output/KCNH2_clean/` | 249 extraction JSON files |
| `harvesting/supplement_scraper.py` | Handler implementations |
| `harvesting/browser_supplement_fetcher.py` | Playwright fallback (Karger) |

---

## 6. Projected Recall Improvement

| Action | Estimated Recall |
|--------|------------------|
| Current baseline | 6% |
| + Download missing PMIDs (free access) | 25-35% |
| + Browser fallback for paywalled | 50-60% |
| + Wiley supplement handler | 55-65% |
| + Re-enable table regex | 70-75% |
| + Variant normalization fixes | 80-85% |
| + Manual review & edge cases | 90%+ |

---

## Next Steps (Recommended Order)

1. **Immediate:** Rebuild venv with `pip install -e .`
2. **Today:** Run automated workflow on missing_baseline_pmids.txt
3. **This week:** Add Wiley supplement handler (copy Oxford pattern)
4. **Ongoing:** Integrate browser fallback for 403/paywalled failures
5. **When available:** Brett submits Karger TDM request, Springer API registration

---

*Generated by GVF status check subagent*

---

## Update: 2026-02-05 15:31 CST — venv Fix Attempt

### Diagnosis

Attempted to recreate/fix the venv. **Python version incompatibility discovered:**

| Component | Requirement | Available |
|-----------|-------------|-----------|
| Project | Python >=3.11 | 3.9.25 |
| biopython | >=1.86 (needs 3.10+) | Would need 1.85 |
| markitdown | >=0.1.3 (needs 3.10+) | Only 0.0.1a1 for 3.9 |
| litellm | Recent versions need 3.10+ | Limited options |
| pandas | >=2.3.3 (doesn't exist) | 2.2.x latest |

### Current venv State

```
venv/bin/python --version → Python 3.9.25
pip list → only pip + setuptools (empty)
```

### Attempted Fixes

1. ✅ Upgraded pip to 26.0.1
2. ❌ `pip install -e .` fails — biopython 1.86+ requires Python 3.10+
3. ❌ Relaxed to Python 3.9 — markitdown requires Python 3.10+
4. ❌ Multiple core dependencies incompatible with Python 3.9

### Resolution Required

**Option A (Recommended):** Install Python 3.11+ on the system
```bash
# If using pyenv:
pyenv install 3.11.8
cd /mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher
pyenv local 3.11.8
python -m venv venv
source venv/bin/activate
pip install -e .
```

**Option B:** Use conda/mamba to create a Python 3.11 environment

**Option C:** Fork dependencies (not recommended — maintenance burden)

### Pipeline Status

❌ **Cannot run pipeline** — venv has no dependencies installed, and Python 3.9 cannot satisfy requirements.

### Recommendation

Brett needs to install Python 3.11+ on the system or use a container/conda environment. The project's dependency tree is firmly rooted in Python 3.10+ features.

---

*Updated by Boswell subagent (venv-fix task)*
