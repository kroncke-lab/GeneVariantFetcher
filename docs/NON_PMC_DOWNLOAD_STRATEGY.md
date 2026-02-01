# Non-PMC Paper Download Strategy

## Test Results Summary (2026-01-28)

Tested 29 non-PMC papers from the test_pmids.txt file:

| Method | Success | Rate |
|--------|---------|------|
| Elsevier API | 13/14 | 93% |
| Wiley TDM API | 3/4 | 75% |
| DOI Resolver | 2/11 | 18% |
| **Total** | **18/29** | **62%** |

## What's Working

### Elsevier API
- Excellent coverage for recent papers
- Returns full XML converted to markdown (100k-350k chars typical)
- API key configured in `.env` as `ELSEVIER_API_KEY`

### Wiley TDM API  
- PDFs downloaded and converted via `markitdown`
- Returns 30k-70k chars of markdown
- API key configured in `.env` as `WILEY_API_KEY`
- Old SICI-style DOIs (pre-2001) may not be available

## What's Not Working (Yet)

### AHA Journals (5 papers)
DOIs like `10.1161/...`

**Issue**: Bot detection — ahajournals.org returns 403 Forbidden for automated requests

**Root cause**: AHA's server detects the requests library as a bot and blocks access, even with browser-like User-Agent headers.

**Solutions**:
1. **Browser automation** — Install Playwright to bypass bot detection
2. **Institutional proxy** — Route through Vanderbilt library proxy
3. **Abstract-only fallback** — The pipeline already handles this; variants can be extracted from abstracts

**Note**: Many older AHA papers become freely available after the embargo period, but bot detection still blocks automated access.

### Other Publishers (7 papers)
- Springer (`10.1007/`)
- Karger (`10.1159/`)
- Mayo Clinic (`10.4065/`)
- Taylor & Francis (`10.1080/`)
- Circulation Journal (`10.1253/`)
- Mary Ann Liebert (`10.1089/`)

**Issue**: Same DOI resolver code bug

**Potential solutions**:
1. Fix DOI resolver scraper argument
2. Add publisher-specific API clients (Springer has an API)
3. Fall back to `browser_fetch.py` with Playwright

### Truly Paywalled (manual needed)
Some papers may never be freely accessible. Use `fetch_manager.py` for semi-automated manual download.

## Next Steps

### Priority 1: Fix DOI Resolver
```python
# In test_non_pmc_download.py, the DOI resolver needs:
from harvesting.supplement_scraper import SupplementScraper
scraper = SupplementScraper(session, target_dir)
doi_resolver = DOIResolver(session=session, paywalled_log=paywalled_log)
content, error = resolver.resolve_and_fetch_fulltext(pmid, doi, scraper)
```

### Priority 2: Add AHA Scraper
Many older AHA papers are freely available. Create an AHA-specific scraper.

### Priority 3: Browser Automation (Optional)
Install Playwright for `browser_fetch.py`:
```bash
pip install playwright
playwright install chromium
```

### Priority 4: Manual Fallback
Use `fetch_manager.py` for remaining paywalled papers:
```bash
python cli/fetch_manager.py tests/fixtures/pmids/non_pmc_papers.csv --convert --run-scout --gene KCNH2
```

## Files Created

- `tests/fixtures/pmids/non_pmc_pmids.txt` - List of non-PMC PMIDs
- `tests/fixtures/pmids/non_pmc_papers.csv` - Full details with DOIs and publishers
- `tests/fixtures/pmids/download_test_results.json` - Test run results
- `scripts/check_pmc_status.py` - Check PMC availability for PMIDs
- `scripts/check_publishers.py` - Identify publishers from DOIs
- `scripts/test_non_pmc_download.py` - Test download via APIs

## Supplement Handling

Note: The current test only checks full-text download. Supplements need separate handling:
- Elsevier supplements: Usually linked in article metadata
- Wiley supplements: Need to scrape article landing page
- Use `SupplementScraper` with publisher-specific handlers
