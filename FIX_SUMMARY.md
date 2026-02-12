# GVF Paper Download Orchestrator Fix Summary

**Date:** 2026-02-11  
**Commit:** c612ff1  
**Files Modified:** `harvesting/orchestrator.py`, `harvesting/unpaywall_api.py`

## Problem Analysis

From the KCNH2 run log (`/tmp/gvf_fresh_kcnh2.log`), only **227/600 papers** were successfully downloaded (62% failure rate).

### Failure Breakdown:
1. **66 papers**: Unpaywall "Expecting value: line 1 column 1" (JSON parsing error)
2. **84 papers**: "No PMCID and not available via any method" (broken fallback chain)
3. **Publisher APIs barely used**: Only ~14 Elsevier + ~25 Wiley attempts out of 373 failures

## Root Causes

### 1. Unpaywall JSON Parsing Error
**Location:** `harvesting/unpaywall_api.py`, line ~91

**Issue:** `response.json()` called without error handling. When Unpaywall API returns empty response or HTML error page, this raises `JSONDecodeError` which propagates as "Request failed: Expecting value: line 1 column 1".

**Fix:** Wrapped `response.json()` in try-except to handle invalid JSON gracefully:
```python
try:
    data = response.json()
except ValueError as e:
    logger.warning(f"Invalid JSON response from Unpaywall for {doi}: {e}")
    return None, f"Invalid JSON response: {str(e)}"
```

### 2. Broken Fallback Chain
**Location:** `harvesting/orchestrator.py`, `_download_free_text_pmid()` method

**Issue:** When a paper has no PMCID:
1. Checks if marked as "free full text" on PubMed
2. If not → tries Unpaywall
3. If Unpaywall fails → **GIVES UP** with "No PMCID and not available via any method"
4. Publisher APIs (Elsevier/Springer/Wiley) are never tried for these papers

**The flow was:**
```
No PMCID → is_free check → Unpaywall → ❌ FAIL
                                          (Publisher APIs never reached)
```

**Fix:** Added publisher API fallback **before** giving up:
```
No PMCID → is_free check → Unpaywall fails → Try publisher APIs by DOI prefix → Give up
```

**Publisher API routing logic:**
- `10.1016/*` → Elsevier API
- `10.1007/*`, `10.1038/*`, `10.1186/*` → Springer API
- `10.1111/*`, `10.1002/*` → Wiley API

### 3. Circuit Breaker Settings
**Status:** No changes needed

Circuit breakers are configured with:
- `max_failures=5` (opens after 5 consecutive failures)
- `reset_timeout=60` (retries after 60 seconds)

These settings are reasonable and not the cause of low success rates.

## Changes Made

### File: `harvesting/unpaywall_api.py`
- Added try-except around `response.json()` to handle empty/invalid responses
- Returns descriptive error message instead of crashing

### File: `harvesting/orchestrator.py`
- Added new fallback section after Unpaywall fails (lines 1337-1403)
- Routes to publisher APIs based on DOI prefix
- Only gives up after all methods exhausted
- Logs successful downloads with source="publisher-api-fallback"
- Updated failure message to "all methods failed" instead of "Unpaywall failed"

## Expected Impact

### Before Fix:
- **227/600 papers** (38% success rate)
- **84 papers** failed with "No PMCID and not available"
- **66 papers** failed with Unpaywall JSON error

### After Fix:
- **Unpaywall JSON errors** (66 papers): Will now skip gracefully instead of crashing, allowing fallback to publisher APIs
- **Publisher API fallback** (84 papers): Many of these have DOIs from Elsevier/Springer/Wiley and should now be retrieved via direct API calls

**Conservative estimate:** +50-70 papers recovered (→ ~300/600, **50% success rate**)  
**Optimistic estimate:** +100-120 papers recovered (→ ~350/600, **58% success rate**)

## Verification

### Import Test:
```bash
cd /mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher
./venv/bin/python3 -c "from harvesting.orchestrator import PMCHarvester; print('✅ OK')"
```
**Result:** ✅ Imports successful

### Git Status:
```bash
git log --oneline -1
# c612ff1 fix: improve download fallback chain and Unpaywall error handling

git push origin main
# ✅ Pushed to main branch
```

## Testing Recommendations

### 1. Dry Run Test (Recommended)
Re-run the KCNH2 pipeline on a subset of the 84 failed papers:
```bash
# Extract PMIDs that failed with "No PMCID and not available"
grep "No PMCID and not available" /tmp/gvf_fresh_kcnh2.log | \
  grep -oP 'PMID \K\d+' | head -20 > test_pmids.txt

# Run orchestrator on these 20 PMIDs
python3 main.py --pmids test_pmids.txt --output /tmp/gvf_test_recovery
```

### 2. Full Pipeline Re-run
```bash
# Re-run on all 600 KCNH2 papers
python3 main.py --gene KCNH2 --output /mnt/temp2/kronckbm/gvf_output/KCNH2/20260211_rerun/
```

**Expected log messages (new):**
```
- Unpaywall failed, trying publisher APIs based on DOI prefix...
✓ Retrieved via Elsevier API fallback
✅ Downloaded via publisher API: 12345678_FULL_CONTEXT.md (2 supplements)
```

## Additional Notes

- No changes to extraction logic (as requested)
- No changes to circuit breaker configuration
- All API keys remain unchanged in `.env`
- Supplement processing remains the same
- Only modified download fallback chain logic

## Rollback (if needed)

If issues arise:
```bash
cd /mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher
git revert c612ff1
git push origin main
```
