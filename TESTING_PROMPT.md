# GeneVariantFetcher Test Pipeline - Revised v2

## Overview

Run the complete GeneVariantFetcher test pipeline on `tests/test_pmids.txt` for the KCNH2 gene. This includes automated downloads with browser crash recovery, enhanced HTML parsing, VPN-assisted manual downloads, figure extraction, data scouting, and LLM-based variant extraction with improved notation normalization.

**Recent Fixes Applied (2026-01-26):**
- Browser crash auto-recovery with 3x retry logic
- Enhanced HTML parsing for AHA, Oxford, Springer, PNAS publishers
- SOCKS proxy support installed
- Fuzzy variant matching with notation normalization (HGVS ↔ legacy)

---

## Prerequisites

### Required Packages
```bash
# Install required dependencies (run once)
pip install "httpx[socks]" --break-system-packages
pip install rapidfuzz --break-system-packages
```

### Environment Setup
- **Vanderbilt VPN active** - Required for accessing paywalled Elsevier, Wiley, and AHA journal articles
- **Claude in Chrome extension installed** - For browser-based downloads of paywalled papers
- **API keys configured** in `.env`:
  - `ANTHROPIC_API_KEY` (required)
  - `OPENAI_API_KEY` (optional fallback)
  - `ELSEVIER_API_KEY` (optional, improves Elsevier access)

---

## Quick Start (Full Pipeline)

For a complete end-to-end test run:

```bash
cd /path/to/GeneVariantFetcher

# Step 1: Run automated pipeline
python tests/run_test_pipeline.py

# Step 2: Check what's missing and download manually (see Step 2 below)

# Step 3: Re-run extraction on new downloads
python cli/run_extraction.py --gene KCNH2 --input tests/test_run_output/KCNH2/*/pmc_fulltext/

# Step 4: Compare results with reference database (uses fuzzy matching by default)
python tests/compare_variants.py \
  --sqlite tests/test_run_output/KCNH2/*/KCNH2.db \
  --excel tests/reference/KCNH2_HeterozygoteDatabase.xlsx \
  --output tests/test_run_output/KCNH2/*/comparison/
```

---

## Step 1: Run Automated Pipeline

Execute the initial automated pipeline which handles PMC open-access papers:

```bash
cd /path/to/GeneVariantFetcher
python tests/run_test_pipeline.py
```

This will:
- Parse PMIDs from `tests/test_pmids.txt`
- Attempt automated downloads from PMC, Elsevier API, Wiley API
- **Auto-retry up to 3 times on browser crashes** (new fix)
- **Enhanced HTML extraction** for more publisher formats (new fix)
- Generate `paywalled_missing.csv` for papers that need manual download
- Create markdown files and run initial data scout on successful downloads

**Output location:** `tests/test_run_output/KCNH2/{timestamp}/pmc_fulltext/`

### Check Pipeline Results

```bash
# Find the latest test run
LATEST=$(ls -td tests/test_run_output/KCNH2/*/ | head -1)
echo "Latest run: $LATEST"

# Check download success rate
echo "=== Download Stats ==="
wc -l ${LATEST}pmc_fulltext/successful_downloads.csv
wc -l ${LATEST}pmc_fulltext/paywalled_missing.csv

# Check extraction success rate
echo "=== Extraction Stats ==="
ls ${LATEST}extractions/*.json 2>/dev/null | wc -l
```

---

## Step 2: Manual Download of Paywalled Papers

### 2.1 Identify Papers Needing Manual Download

```bash
# List papers that failed automated download
cat tests/test_run_output/KCNH2/*/pmc_fulltext/paywalled_missing.csv

# Prioritize by expected variant count (if reference data available)
python tests/compare_variants.py --list-priority-missing \
  --sqlite tests/test_run_output/KCNH2/*/KCNH2.db \
  --excel tests/reference/KCNH2_HeterozygoteDatabase.xlsx
```

### 2.2 Download Papers by Publisher

**Priority order (by typical success rate with Vanderbilt VPN):**

#### A. PMC Direct PDFs (Open Access)
- Pattern: `https://pmc.ncbi.nlm.nih.gov/articles/PMC{id}/pdf/{filename}.pdf`
- Should download automatically; failures usually due to rate limiting

#### B. AHA/Circulation Journals (VPN IP-based access)
- Navigate to: `https://www.ahajournals.org/doi/pdf/{DOI}`
- Example: `https://www.ahajournals.org/doi/pdf/10.1161/01.cir.102.10.1178`
- VPN should auto-grant access - no login required
- **Check for supplemental materials** in "Data Supplement" section

#### C. Wiley Journals (May require institutional login)
- Navigate to: `https://onlinelibrary.wiley.com/doi/pdf/{DOI}`
- If paywall shown, click "Access through your institution" → search "Vanderbilt"

#### D. Elsevier/ScienceDirect (Often blocked by CAPTCHA)
- Navigate to: `https://www.sciencedirect.com/science/article/pii/{PII}/pdfft`
- If CAPTCHA appears, solve manually then proceed

### 2.3 Download Supplemental Materials

**CRITICAL: Supplements often contain the most valuable variant tables!**

For each paper, navigate to the article's main page and look for:
- "Supplementary Materials" / "Supporting Information" / "Data Supplement"

**File naming convention:**
```
PMID_{pmid}_supplement_1.pdf
PMID_{pmid}_supplement_table_S1.xlsx
```

### 2.4 Organize Downloaded Files

```bash
# Create pdfs directory
mkdir -p tests/test_run_output/KCNH2/{timestamp}/pmc_fulltext/pdfs/

# Use organize script (auto-matches recent downloads to PMIDs)
python scripts/organize_downloads.py \
  tests/test_run_output/KCNH2/{timestamp}/pmc_fulltext/pdfs/ \
  --pmids "12345,67890" \
  --since 120
```

---

## Step 3: Process Downloaded PDFs

### 3.1 Convert PDFs to Markdown

```python
from harvesting.format_converters import FormatConverter

converter = FormatConverter()

# Convert main PDF
md_content = converter.pdf_to_markdown("path/to/PMID_12345.pdf")
with open("12345_main.md", "w") as f:
    f.write(md_content)

# Convert supplement PDFs
supp_md = converter.pdf_to_markdown("path/to/PMID_12345_supplement.pdf")
with open("12345_supplement.md", "w") as f:
    f.write(supp_md)
```

### 3.2 Create FULL_CONTEXT.md Files

```python
# Template for FULL_CONTEXT.md
full_context = f"""# PMID: {pmid}
# Title: {title}
# Gene: KCNH2

## Main Text
{main_markdown}

## Supplementary Materials
{supplement_markdown}
"""

with open(f"{pmid}_FULL_CONTEXT.md", "w") as f:
    f.write(full_context)
```

---

## Step 4: Extract Figures (Optional)

For papers with pedigrees, tables, or variant diagrams:

```python
from harvesting.format_converters import FormatConverter
from pathlib import Path

converter = FormatConverter()
figure_dir = Path(f"tests/test_run_output/KCNH2/{timestamp}/pmc_fulltext/{pmid}_figures")
figure_dir.mkdir(exist_ok=True)

converter.extract_pdf_figures(
    pdf_path="path/to/PMID_12345.pdf",
    output_dir=figure_dir
)
```

---

## Step 5: Run Data Scout (LLM Analysis)

```python
from pipeline.data_scout import GeneticDataScout

scout = GeneticDataScout()

with open(f"{pmid}_FULL_CONTEXT.md") as f:
    full_text = f.read()

# Identify high-value data zones
zones = scout.analyze(text=full_text, gene_symbol="KCNH2", pmid=pmid)

# Save results
with open(f"{pmid}_DATA_ZONES.json", "w") as f:
    json.dump(zones.to_dict(), f, indent=2)

# Generate condensed markdown
condensed_md = scout.generate_condensed_markdown(zones, full_text)
with open(f"{pmid}_DATA_ZONES.md", "w") as f:
    f.write(condensed_md)
```

---

## Step 6: Run LLM Extraction

### 6.1 Extract Variant Data

```python
from pipeline.extraction import ExpertExtractor
from utils.models import Paper

extractor = ExpertExtractor(
    models=["claude-sonnet-4-20250514", "gpt-4o"],
    gene_symbol="KCNH2"
)

paper = Paper(
    pmid=pmid,
    title=title,
    full_text=full_context_md,
    data_zones_md=data_zones_md
)

result = extractor.extract(paper)

with open(f"{pmid}_extraction.json", "w") as f:
    json.dump(result.to_dict(), f, indent=2)
```

### 6.2 Aggregate Results

```python
from pipeline.aggregation import ResultAggregator

aggregator = ResultAggregator(output_dir="tests/test_run_output/KCNH2/{timestamp}")
summary = aggregator.aggregate_all()
aggregator.export_to_csv("KCNH2_variants_extracted.csv")
```

---

## Step 7: Verify Results

### 7.1 Compare Against Reference Database

The comparison script now uses **fuzzy matching by default** with automatic notation normalization (HGVS ↔ legacy format).

```bash
# Full comparison with fuzzy matching (default)
python tests/compare_variants.py \
  --sqlite tests/test_run_output/KCNH2/*/KCNH2.db \
  --excel tests/reference/KCNH2_HeterozygoteDatabase.xlsx \
  --output tests/test_run_output/KCNH2/*/comparison/

# Exact matching only (stricter)
python tests/compare_variants.py \
  --sqlite tests/test_run_output/KCNH2/*/KCNH2.db \
  --excel tests/reference/KCNH2_HeterozygoteDatabase.xlsx \
  --variant_match_mode exact \
  --output tests/test_run_output/KCNH2/*/comparison/

# Custom fuzzy threshold (default is 0.80)
python tests/compare_variants.py \
  --sqlite tests/test_run_output/KCNH2/*/KCNH2.db \
  --excel tests/reference/KCNH2_HeterozygoteDatabase.xlsx \
  --fuzzy_threshold 0.85 \
  --output tests/test_run_output/KCNH2/*/comparison/
```

**Output files generated:**
- `comparison/report.md` - Summary comparison report
- `comparison/exact_matches.csv` - Variants matching exactly
- `comparison/fuzzy_matches.csv` - Variants matching via fuzzy/notation normalization
- `comparison/missing_in_sqlite.csv` - Coverage gaps (in reference, not extracted)
- `comparison/missing_in_excel.csv` - New discoveries (extracted, not in reference)
- `comparison/count_mismatches.csv` - Carrier count discrepancies

### 7.2 Understand Notation Normalization

The comparison script automatically converts between formats:

| Extraction (HGVS) | Reference (Legacy) | Status |
|-------------------|-------------------|--------|
| p.Arg328Cys | R328C | ✅ Matched |
| p.Ala490Thr | A490T | ✅ Matched |
| p.Pro926Alafs*14 | P926fsX | ✅ Fuzzy match |
| p.Arg123_Lys130del | R123_K130del | ✅ Matched |

### 7.3 Check for Missing Papers

```bash
# Count successfully processed papers
ls tests/test_run_output/KCNH2/{timestamp}/pmc_fulltext/*_FULL_CONTEXT.md | wc -l

# Check which PMIDs are still missing
python -c "
from pathlib import Path

with open('tests/test_pmids.txt') as f:
    expected = {line.split()[0] for line in f if line.strip() and not line.startswith('#') and line.split()[0].isdigit()}

processed = {p.stem.split('_')[0] for p in Path('tests/test_run_output/KCNH2').rglob('*_FULL_CONTEXT.md')}

missing = expected - processed
print(f'Missing {len(missing)} papers:')
for pmid in sorted(missing):
    print(f'  {pmid}')
"
```

---

## Expected Output Structure

```
tests/test_run_output/KCNH2/{timestamp}/
├── pmc_fulltext/
│   ├── paywalled_missing.csv      # Papers needing manual download
│   ├── successful_downloads.csv   # Successfully downloaded papers
│   ├── pdfs/                      # Manually downloaded PDFs
│   ├── {pmid}_FULL_CONTEXT.md     # Combined main + supplements
│   ├── {pmid}_DATA_ZONES.json     # Scout analysis results
│   ├── {pmid}_DATA_ZONES.md       # Condensed high-value content
│   └── {pmid}_figures/            # Extracted figures
├── extractions/
│   └── {pmid}_extraction.json     # LLM extraction results
├── comparison/
│   ├── report.md                  # Summary report
│   ├── exact_matches.csv          # Perfect matches
│   ├── fuzzy_matches.csv          # Notation-normalized matches
│   ├── missing_in_sqlite.csv      # Coverage gaps
│   └── missing_in_excel.csv       # New discoveries
├── KCNH2.db                       # SQLite database
├── KCNH2_penetrance_summary.json  # Aggregated data
└── pipeline_log.txt               # Full pipeline log
```

---

## Troubleshooting

### Common Issues

1. **Browser crashes during download**
   - Now handled automatically with 3x retry logic
   - If still failing, check Playwright installation: `playwright install chromium`

2. **"Text too short" extraction errors**
   - Enhanced HTML parser now handles more publisher formats
   - If still occurring, check if article is behind paywall (abstract-only)

3. **SOCKS proxy errors**
   - Fixed by installing: `pip install "httpx[socks]"`
   - Verify: `python -c "import httpx; print('OK')"`

4. **Low match rate with reference database**
   - Now uses fuzzy matching by default (0.80 threshold)
   - Check `fuzzy_matches.csv` for notation-based matches
   - Try lowering threshold: `--fuzzy_threshold 0.75`

5. **Elsevier CAPTCHA blocks**
   - User must solve CAPTCHA manually in browser
   - Consider using Elsevier API key for better access

6. **PDF conversion fails**
   - Some PDFs have DRM or are image-only
   - Use OCR fallback: `converter.pdf_to_markdown(path, use_ocr=True)`

### Re-running Failed Papers

```bash
# Re-run pipeline for specific PMIDs only
python cli/run_extraction.py \
  --gene KCNH2 \
  --pmids "12345,67890" \
  --output tests/test_run_output/KCNH2/{timestamp}/
```

### Viewing Fix Details

```bash
# See all fixes applied in this test run
cat tests/test_run_output/KCNH2/*/FIXES_APPLIED.md
```

---

## Performance Targets

| Metric | Baseline | Target | Current |
|--------|----------|--------|---------|
| Download success rate | 62.7% | 80% | Track in results |
| Extraction success rate | 35.3% | 60% | Track in results |
| Variant coverage (vs reference) | ~10% | 80% | Track in comparison |
| Exact + fuzzy match rate | 0.2% | 90% | **4.5%** (improved) |

---

## Notes for Claude (AI Assistant)

When executing this pipeline:

1. **Browser crashes are now auto-recovered** - No need to manually restart
2. **HTML parsing is enhanced** - Most publishers should extract full text now
3. **Always check VPN status first** - Many papers require institutional access
4. **Download supplements proactively** - They often contain the most valuable variant tables
5. **Use fuzzy comparison** - Default mode now handles notation differences automatically
6. **Name files consistently** - Include PMID in all filenames for easy tracking
7. **Verify downloads** - Check PDF file sizes (>10KB usually means real content)
8. **Run data scout before extraction** - It identifies the most valuable sections
9. **Save intermediate results** - Pipeline can be resumed from any step
10. **Check FIXES_APPLIED.md** - Documents all improvements made to the pipeline

---

*Last updated: 2026-01-26*
*See FIXES_APPLIED.md for detailed fix documentation*
