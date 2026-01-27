# GeneVariantFetcher Download Pipeline

## Overview

Download papers for the test PMIDs in `tests/test_pmids.txt` (51 PMIDs) for KCNH2 gene. **Designed for Claude Code cowork mode** - fully automated browser-based downloading with VPN access.

---

## Prerequisites

- **Vanderbilt VPN active** - Required for accessing paywalled journals
- **Playwright installed** with chromium browser

---

## Quick Start

```bash
cd /Users/kronckbm/GitRepos/GeneVariantFetcher

# Step 1: Download papers via browser automation (uses VPN + Chrome profile)
./venv/bin/python tests/run_browser_download.py

# Step 2: Process downloaded PDFs (convert to markdown + run data scout)
./venv/bin/python tests/run_post_download.py
```

---

## Step 1: Browser-Based Download

Downloads papers using Playwright browser automation with your Chrome profile (for VPN/cookies):

```bash
./venv/bin/python tests/run_browser_download.py
```

This will:
- Open each paper URL in a browser (uses real Chrome profile)
- Attempt to find and download PDFs automatically
- Auto-retry up to 3 times on browser crashes
- Save PDFs to `tests/test_run_output/KCNH2/downloads/pmc_fulltext/`
- Update CSV with status (browser_done, paywall, captcha_blocked)

### Download Options

```bash
# Download only first 5 papers (for testing)
./venv/bin/python tests/run_browser_download.py --max 5

# Retry only paywalled papers (after fixing VPN)
./venv/bin/python tests/run_browser_download.py --retry-paywall

# Download specific PMIDs
./venv/bin/python tests/run_browser_download.py --pmids 15840476,10973849
```

### Check Download Results

```bash
# Check download status counts
cat tests/test_run_output/KCNH2/downloads/test_pmids_download.csv | cut -d, -f3 | sort | uniq -c

# Count downloaded PDFs
ls tests/test_run_output/KCNH2/downloads/pmc_fulltext/*.pdf 2>/dev/null | wc -l

# Check what failed
grep -E "paywall|captcha" tests/test_run_output/KCNH2/downloads/test_pmids_download.csv
```

---

## Step 2: Post-Download Processing

Convert downloaded PDFs to markdown and run data scout:

```bash
./venv/bin/python tests/run_post_download.py
```

This will:
- Convert all PDFs to markdown format
- Combine main text + supplements into `{PMID}_FULL_CONTEXT.md`
- Run Data Scout to create condensed `{PMID}_DATA_ZONES.md` files

---

## Retry Failed Downloads

If some papers are paywalled or CAPTCHA-blocked:

```bash
# Retry paywalled papers (make sure VPN is active)
./venv/bin/python tests/run_browser_download.py --retry-paywall

# For CAPTCHA-blocked papers, use interactive mode
./venv/bin/python -m cli.browser_fetch \
  tests/test_run_output/KCNH2/downloads/test_pmids_download.csv \
  --retry-captcha --wait-for-captcha --interactive
```

---

## Check Progress

```bash
# Count FULL_CONTEXT.md files created
ls tests/test_run_output/KCNH2/downloads/pmc_fulltext/*_FULL_CONTEXT.md 2>/dev/null | wc -l

# Check which PMIDs are still missing
./venv/bin/python -c "
from pathlib import Path

with open('tests/test_pmids.txt') as f:
    expected = {line.split()[0] for line in f if line.strip() and not line.startswith('#') and line.split()[0].isdigit()}

processed = {p.stem.split('_')[0] for p in Path('tests/test_run_output/KCNH2/downloads/pmc_fulltext').glob('*_FULL_CONTEXT.md')}

missing = expected - processed
print(f'Downloaded: {len(processed)}/{len(expected)} papers')
if missing:
    print(f'Missing {len(missing)} papers:')
    for pmid in sorted(missing):
        print(f'  {pmid}')
"
```

---

## Output Structure

```
tests/test_run_output/KCNH2/downloads/
├── test_pmids_download.csv        # Download status tracking
└── pmc_fulltext/
    ├── {pmid}_Main_Text.pdf       # Downloaded PDFs
    ├── {pmid}_FULL_CONTEXT.md     # Converted markdown
    ├── {pmid}_DATA_ZONES.md       # Condensed high-value content
    └── {pmid}_supplements/        # Supplement files
```

---

## Troubleshooting

1. **Browser crashes during download**
   - Handled automatically with 3x retry logic
   - If still failing: `./venv/bin/playwright install chromium`

2. **VPN not working / Paywall errors**
   - Make sure Vanderbilt VPN is connected
   - Retry with: `./venv/bin/python tests/run_browser_download.py --retry-paywall`

3. **CAPTCHA blocks**
   - Use interactive mode to solve manually

4. **PDF conversion fails**
   - Some PDFs have DRM or are image-only
   - Check the PDF manually and re-download if corrupted

---

## Target

Download success rate with VPN: **80%+** (40+ of 51 papers)

---

*Last updated: 2026-01-27*
