# GeneVariantFetcher Test Pipeline - Revised v3

## Overview

Run the complete GeneVariantFetcher test pipeline on `tests/test_pmids.txt` (51 PMIDs) for the KCNH2 gene. **Designed for Claude Code cowork mode** - fully automated browser-based downloading with VPN access.

**Updated 2026-01-26:**
- New automated runner scripts for Claude Code cowork
- Browser crash auto-recovery with 3x retry logic
- Uses real Chrome profile for VPN/institutional access
- Fuzzy variant matching with notation normalization (HGVS ↔ legacy)

---

## Prerequisites

### Environment Setup
- **Vanderbilt VPN active** - Required for accessing paywalled journals
- **API keys configured** in `.env` (do not modify this file)
- **Playwright installed** with chromium browser

---

## Quick Start (Claude Code Cowork Mode)

**For fully automated execution with Claude Code:**

```bash
cd /Users/kronckbm/GitRepos/GeneVariantFetcher

# Step 1: Download papers via browser automation (uses VPN + Chrome profile)
./venv/bin/python tests/run_browser_download.py

# Step 2: Process downloaded PDFs (convert to markdown + run data scout)
./venv/bin/python tests/run_post_download.py

# Step 3: Run full extraction pipeline
./venv/bin/python tests/run_test_pipeline.py

# Step 4: Compare results (if reference database available)
./venv/bin/python tests/compare_variants.py \
  --sqlite tests/test_run_output/KCNH2/*/KCNH2.db \
  --output tests/test_run_output/KCNH2/comparison/
```

### Runner Script Options

```bash
# Download only first 5 papers (for testing)
./venv/bin/python tests/run_browser_download.py --max 5

# Use Claude API to help find download links
./venv/bin/python tests/run_browser_download.py --use-claude

# Retry only paywalled papers (after fixing VPN)
./venv/bin/python tests/run_browser_download.py --retry-paywall

# Download specific PMIDs
./venv/bin/python tests/run_browser_download.py --pmids 15840476,10973849
```

---

## Step 1: Browser-Based Download

Downloads papers using Playwright browser automation with your Chrome profile (for VPN/cookies):

```bash
cd /Users/kronckbm/GitRepos/GeneVariantFetcher
./venv/bin/python tests/run_browser_download.py
```

This will:
- Open each paper URL in a browser (uses real Chrome profile)
- Attempt to find and download PDFs automatically
- **Auto-retry up to 3 times on browser crashes**
- Save PDFs to `tests/test_run_output/KCNH2/downloads/pmc_fulltext/`
- Update CSV with status (browser_done, paywall, captcha_blocked)

**Output location:** `tests/test_run_output/KCNH2/downloads/pmc_fulltext/`

### Check Download Results

```bash
# Check download status
cat tests/test_run_output/KCNH2/downloads/test_pmids_download.csv | grep -c "browser_done"

# List downloaded PDFs
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
- Organize supplement files into `{PMID}_supplements/` folders

### Retry Failed Downloads

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

## Step 3: Run Extraction Pipeline

Run the full extraction pipeline on downloaded papers:

```bash
./venv/bin/python tests/run_test_pipeline.py
```

This will:
- Use the downloaded `FULL_CONTEXT.md` files
- Run LLM-based variant extraction (Claude/GPT-4o cascade)
- Aggregate results into SQLite database
- Create `KCNH2_penetrance_summary.json`

Alternatively, run extraction on a specific folder:

```bash
./venv/bin/python tests/run_folder_extraction.py \
  --folder tests/test_run_output/KCNH2/downloads/pmc_fulltext \
  --gene KCNH2
```

---

## Step 4: Verify Results

### Check Extraction Results

```bash
# Count downloaded papers
ls tests/test_run_output/KCNH2/downloads/pmc_fulltext/*_FULL_CONTEXT.md 2>/dev/null | wc -l

# Count extractions
ls tests/test_run_output/KCNH2/*/extractions/*.json 2>/dev/null | wc -l

# Check which PMIDs are still missing
./venv/bin/python -c "
from pathlib import Path

with open('tests/test_pmids.txt') as f:
    expected = {line.split()[0] for line in f if line.strip() and not line.startswith('#') and line.split()[0].isdigit()}

processed = {p.stem.split('_')[0] for p in Path('tests/test_run_output/KCNH2').rglob('*_FULL_CONTEXT.md')}

missing = expected - processed
print(f'Downloaded: {len(processed)}/{len(expected)} papers')
print(f'Missing {len(missing)} papers:')
for pmid in sorted(missing)[:10]:
    print(f'  {pmid}')
if len(missing) > 10:
    print(f'  ... and {len(missing)-10} more')
"
```

### Compare Against Reference Database (Optional)

If you have a reference Excel file:

```bash
./venv/bin/python tests/compare_variants.py \
  --sqlite tests/test_run_output/KCNH2/*/KCNH2.db \
  --excel /path/to/KCNH2_reference.xlsx \
  --output tests/test_run_output/KCNH2/comparison/
```

**Output files:**
- `comparison/report.md` - Summary report
- `comparison/exact_matches.csv` - Exact variant matches
- `comparison/fuzzy_matches.csv` - Notation-normalized matches
- `comparison/missing_in_sqlite.csv` - Coverage gaps

---

## Expected Output Structure

```
tests/test_run_output/KCNH2/
├── downloads/                         # Browser download output
│   ├── test_pmids_download.csv        # Download status tracking
│   └── pmc_fulltext/
│       ├── {pmid}_Main_Text.pdf       # Downloaded PDFs
│       ├── {pmid}_FULL_CONTEXT.md     # Converted markdown
│       ├── {pmid}_DATA_ZONES.md       # Condensed high-value content
│       └── {pmid}_supplements/        # Supplement files
├── {timestamp}/                       # Pipeline run output
│   ├── extractions/
│   │   └── KCNH2_PMID_{pmid}.json     # LLM extraction results
│   ├── KCNH2.db                       # SQLite database
│   └── KCNH2_penetrance_summary.json  # Aggregated data
└── comparison/                        # Optional comparison output
    ├── report.md
    └── *.csv
```

---

## Troubleshooting

### Common Issues

1. **Browser crashes during download**
   - Handled automatically with 3x retry logic
   - If still failing: `./venv/bin/playwright install chromium`

2. **VPN not working / Paywall errors**
   - Make sure Vanderbilt VPN is connected
   - Retry with: `./venv/bin/python tests/run_browser_download.py --retry-paywall`

3. **CAPTCHA blocks**
   - Use interactive mode to solve manually:
   ```bash
   ./venv/bin/python -m cli.browser_fetch \
     tests/test_run_output/KCNH2/downloads/test_pmids_download.csv \
     --retry-captcha --wait-for-captcha --interactive
   ```

4. **PDF conversion fails**
   - Some PDFs have DRM or are image-only
   - Check the PDF manually and re-download if corrupted

5. **"Text too short" extraction errors**
   - Paper may be abstract-only (paywall blocked full text)
   - Re-download with VPN active

### Re-running Specific Papers

```bash
# Download specific PMIDs
./venv/bin/python tests/run_browser_download.py --pmids "15840476,10973849"

# Re-process downloads only
./venv/bin/python tests/run_post_download.py
```

---

## Performance Targets

| Metric | Target |
|--------|--------|
| Download success rate (with VPN) | 80%+ |
| Extraction success rate | 60%+ |
| Papers with usable full text | 40/51 |

---

## Notes for Claude Code Cowork

When executing this pipeline:

1. **Run scripts in order**: download → post-process → extract
2. **Check VPN first**: Run `curl -s ifconfig.me` - should show Vanderbilt IP
3. **Browser uses real profile**: Cookies/sessions from Chrome are available
4. **Auto-retry on failures**: Browser crashes are handled automatically
5. **Check CSV status**: `grep -c browser_done test_pmids_download.csv`
6. **Retry paywalled papers**: Use `--retry-paywall` after fixing VPN

### Quick Commands for Claude Code

```bash
# Full automated run
cd /Users/kronckbm/GitRepos/GeneVariantFetcher
./venv/bin/python tests/run_browser_download.py
./venv/bin/python tests/run_post_download.py
./venv/bin/python tests/run_test_pipeline.py

# Check progress
cat tests/test_run_output/KCNH2/downloads/test_pmids_download.csv | cut -d, -f3 | sort | uniq -c
```

---

*Last updated: 2026-01-26*
