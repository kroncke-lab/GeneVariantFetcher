# CODEX_TASKS.md - GVF Recall Improvement

## Status: MOSTLY COMPLETE ✅
Last updated: 2026-02-04

### Summary
- Tasks 1-5: ✅ Complete
- Task 6: ⏸️ Karger blocked by Cloudflare (requires TDM request)

---

## Goal
Improve variant extraction recall from 5% to 90% by fixing supplementary table downloads.

## Current Results
- **2370 variants** extracted from 249 papers
- **5179 total carriers** identified
- **201 papers** with variant data
- Exceeds Excel baseline (1359 variant-PMID pairs)

---

## Task Status

### ✅ Task 1: Audit supplement_scraper.py - COMPLETE
Comprehensive docstring added documenting:
- Supported publishers (Nature, Elsevier, Karger, Springer/BMC, Oxford)
- Generic fallback patterns
- Known gaps and limitations

### ✅ Task 2: Create Karger publisher handler - COMPLETE (but blocked)
`scrape_karger_supplements()` method implemented but Karger uses Cloudflare
protection that blocks automated requests. Browser fallback exists at
`harvesting/browser_supplement_fetcher.py`.

### ✅ Task 3: Supplement reference parser - COMPLETE
Created `harvesting/supplement_reference_parser.py` with:
- `parse_supplement_references()` - finds table/figure refs in text
- `extract_supplement_urls_from_text()` - extracts explicit URLs
- `check_supplement_gap()` - validates downloads vs expectations

### ✅ Task 4: URL extraction integration - COMPLETE
Supplement URLs are now extracted from article text and added to download list.

### ✅ Task 5: Post-download validation - COMPLETE (2026-02-04)
`check_supplement_gap()` integrated into `orchestrator.py`:
- Warns when expected supplements weren't downloaded
- Logs variant count mismatches
- Flags papers needing manual review

### ⏸️ Task 6: Karger E2E test - BLOCKED
Cannot test PMID 26496715 until Karger TDM access is obtained.
Draft request prepared; awaiting Brett to submit.

---

## Next Steps
1. Brett submits Karger TDM request
2. Brett registers for Springer API access
3. Run browser fetch for remaining paywalled papers
4. Aggregate JSON extractions to SQLite for comparison
5. Extend pipeline to additional genes

---

## Test Commands

```bash
# Run golden test harvest
cd /mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher
python -m cli.automated_workflow --gene KCNH2 --output /mnt/temp2/kronckbm/gvf_output/KCNH2/test_run

# Check extraction results
ls /mnt/temp2/kronckbm/gvf_output/KCNH2/20260202_173749/extractions/*.json | wc -l
python /tmp/measure_recall_v2.py
```
