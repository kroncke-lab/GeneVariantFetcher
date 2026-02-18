# Manual Assist Module

Tools for handling papers that can't be automatically downloaded (~10% of corpus).

## Problem

Some papers require human intervention:
- CAPTCHA challenges
- Institutional login required
- Complex paywalls
- Non-standard PDF locations

## Solution

Semi-automated workflow that:
1. Opens browser to paper URL (via institutional proxy when possible)
2. Human solves CAPTCHA / logs in / clicks download
3. Tool watches Downloads folder and organizes the PDF
4. Progress tracked in CSV

## Current Tools

### `semi_manual_download.py`

Interactive download assistant using Playwright.

```bash
# Install dependencies
pip install pandas playwright nest_asyncio
playwright install chromium

# Run
python semi_manual_download.py --input papers_to_download.csv --output downloads/
```

**Input CSV format:**
```csv
pmid,title,publisher,doi
12345678,Paper Title,Wiley,10.1002/example
```

**Supported proxies (Vanderbilt):**
- Wiley
- Elsevier  
- Springer
- Oxford
- Nature
- Cell
- NEJM
- BMJ
- AHA Journals

## Roadmap

### Phase 1: Current (v0.1)
- [x] Basic interactive download with proxy support
- [x] Progress tracking
- [x] Skip/quit commands

### Phase 2: GVF Integration
- [ ] Auto-generate `papers_to_download.csv` from GVF failed downloads
- [ ] Import downloaded PDFs back into GVF output structure
- [ ] Update GVF status tracking

### Phase 3: Smarter Automation
- [ ] CAPTCHA detection (pause only when needed)
- [ ] Session persistence (stay logged in)
- [ ] Batch mode for same-publisher papers
- [ ] Retry logic for transient failures

### Phase 4: Multi-Institution Support
- [ ] Config file for proxy URLs
- [ ] Support for different auth methods (Shibboleth, etc.)
- [ ] VPN integration option

## Integration with GVF Pipeline

```
GVF Pipeline
    │
    ├── Auto-download (~90%)
    │       ↓
    │   gvf_output/pdfs/
    │
    └── Failed downloads (~10%)
            ↓
        papers_to_download.csv
            ↓
        [Human runs semi_manual_download.py]
            ↓
        manual_assist/downloads/
            ↓
        [Import script moves to gvf_output/pdfs/]
```

## Notes

- Requires human presence (not headless)
- Best run in batches of 20-50 papers
- Vanderbilt proxy requires VUNet login
- Some publishers have download rate limits
