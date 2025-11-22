# GeneVariantFetcher - Command-Line Tools

A collection of Python tools for extracting genetic variant data, patient-level information, and full-text articles from PubMed literature. All tools are designed to run locally on your computer without requiring a web browser.

## Table of Contents

- [What This Tool Does](#what-this-tool-does)
- [Local Setup](#local-setup)
- [Available Tools](#available-tools)
- [Using the PMC Harvester](#using-the-pmc-harvester)
- [Using the Extraction Pipeline](#using-the-extraction-pipeline)
- [Using the Clinical Data Triage Tool](#using-the-clinical-data-triage-tool)
- [Troubleshooting](#troubleshooting)

---

## What This Tool Does

**GeneVariantFetcher** provides command-line tools for:

1. **Full-Text Harvester** - Download complete PubMed Central articles with:
   - Full-text content (not just abstracts)
   - All supplemental materials (Excel, Word, PDFs)
   - Everything consolidated into markdown files for LLM processing

2. **Extraction Pipeline** - Extract structured biomedical data from articles:
   - Individual-level variant and phenotype data
   - Clinical status (affected/unaffected)
   - Age, sex, and demographic information
   - Evidence sentences supporting each extraction

3. **Clinical Data Triage** - Process and triage clinical data:
   - Extract and organize patient-level information
   - Filter and classify clinical variants
   - Generate structured outputs for analysis

---

## Local Setup

### Prerequisites

- Python 3.11 or higher
- pip (Python package manager)
- Git

### Step 1: Clone the Repository

```bash
git clone https://github.com/kroncke-lab/GeneVariantFetcher.git
cd GeneVariantFetcher
```

### Step 2: Set Up Python Environment

**Option A: Using venv (recommended)**

```bash
# Create virtual environment
python -m venv venv

# Activate it
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate
```

**Option B: Using conda**

```bash
conda create -n genevariant python=3.11
conda activate genevariant
```

### Step 3: Install Dependencies

```bash
pip install -e .
```

This installs all required packages listed in `pyproject.toml`.

### Step 4: Set Up OpenAI API Key (for AI extraction features)

The extraction features require an OpenAI API key.

**Option A: Using environment variables (recommended)**

```bash
# On macOS/Linux:
export AI_INTEGRATIONS_OPENAI_API_KEY="your-api-key-here"
export AI_INTEGRATIONS_OPENAI_BASE_URL="https://api.openai.com/v1"

# On Windows (PowerShell):
$env:AI_INTEGRATIONS_OPENAI_API_KEY="your-api-key-here"
$env:AI_INTEGRATIONS_OPENAI_BASE_URL="https://api.openai.com/v1"
```

**Option B: Create a .env file**

Create a file named `.env` in the project root:

```env
AI_INTEGRATIONS_OPENAI_API_KEY=your-api-key-here
AI_INTEGRATIONS_OPENAI_BASE_URL=https://api.openai.com/v1
```

---

## Available Tools

### 1. Automated PMID Discovery (`pubmind_fetcher.py`) 🆕

Automatically discover relevant papers using PubMind - no manual PMID list needed!

**Quick Start:**
```bash
python example_automated_workflow.py BRCA1 --email your@email.com
```

**See detailed documentation:** [PUBMIND_README.md](PUBMIND_README.md)

**Features:**
- Queries PubMind database for variant-focused papers
- Falls back to enhanced PubMed queries if needed
- Integrates with harvester and extraction pipeline
- Complete hands-free workflow from gene to data

### 2. PMC Full-Text Harvester (`harvest_pmc_fulltext.py`)

Download full-text articles and supplemental materials from PubMed Central.

**Quick Start:**
```bash
python harvest_pmc_fulltext.py
```

**See detailed documentation:** [HARVESTER_QUICKSTART.md](HARVESTER_QUICKSTART.md)

### 3. Extraction Pipeline (`pipeline.py`)

Extract structured biomedical data from articles using AI.

**Quick Start:**
```bash
python example_usage.py
```

**See detailed documentation:** [README_PIPELINE.md](README_PIPELINE.md)

### 4. Clinical Data Triage (`clinical_data_triage.py`)

Process and triage clinical data from various sources.

**Quick Start:**
```bash
python example_triage.py
```

**See detailed documentation:** [TRIAGE_README.md](TRIAGE_README.md)

---

## Using the PMC Harvester

The harvesting tool downloads **full-text** PubMed Central articles with all supplemental materials. This provides much richer data than abstracts alone.

### When to Use the Harvester

- You need full article text, not just abstracts
- You want supplemental tables with patient-level data
- You're preparing articles for detailed LLM extraction
- You have a list of PMIDs to process

### Option 1: Quick Start (Edit the Script)

1. **Get your PMIDs**
   - Create a list of PubMed IDs you want to download

2. **Edit the harvesting script**
   - Open `harvest_pmc_fulltext.py`
   - Find line 333 (look for `pmids = [`)
   - Replace with your PMIDs:

   ```python
   pmids = [
       '35443093',
       '33442691',
       '34931732',
   ]
   ```

3. **Set your email** (required by NCBI)
   - Find line 32: `Entrez.email = "your.email@example.com"`
   - Change to your actual email

4. **Run the script**

   ```bash
   python harvest_pmc_fulltext.py
   ```

5. **Check the output**
   - Look in the `pmc_harvest/` folder
   - Each article creates a `{PMID}_FULL_CONTEXT.md` file
   - Also creates:
     - `successful_downloads.csv` - what worked
     - `paywalled_missing.csv` - what didn't work

### Option 2: Using a PMID File

1. **Prepare your PMID file**
   - Create a text file (e.g., `my_pmids.txt`)
   - One PMID per line:
   ```
   35443093
   33442691
   34931732
   ```

2. **Use the example script**
   ```bash
   python example_harvest_from_pubmind.py
   ```

3. **Choose option 1**
   - Enter the path to your PMID file
   - Script will process all PMIDs

### Option 3: Python API (Advanced)

```python
from harvest_pmc_fulltext import PMCHarvester

# Load PMIDs from a file
with open('my_pmids.txt', 'r') as f:
    pmids = [line.strip() for line in f if line.strip()]

# Create harvester
harvester = PMCHarvester(output_dir="my_harvest")

# Harvest (delay=2.0 means 2 seconds between requests)
harvester.harvest(pmids[:50], delay=2.0)  # Start with first 50

print(f"Check the my_harvest/ folder for results")
```

### Understanding Harvester Output

Each successful article creates a markdown file with:

```markdown
# MAIN TEXT

## [Article Title]

### Abstract
Full abstract text...

### Introduction
Introduction section...

### Methods
Methods section...

### Results
Results with all details...

# SUPPLEMENTAL FILE 1: Table_S1.xlsx

#### Sheet: Patient_Cohort

| Patient_ID | Gene | Variant | Age | Phenotype |
|------------|------|---------|-----|-----------|
| P001       | BRCA1 | c.1234G>A | 45 | Breast cancer |
| P002       | BRCA1 | c.5678C>T | 52 | Ovarian cancer |
```

**This is perfect for:**
- Feeding to LLMs for extraction
- Manual review of case details
- Finding patient-level data in supplemental tables

### Harvester Tips

**Rate Limiting:**
- Use `delay=2.0` (2 seconds) for moderate batches (recommended)
- Use `delay=3.0` for large batches (>500 PMIDs)
- Use `delay=1.0` for very small batches (<20 PMIDs)

**Success Rates:**
- Only ~30% of PubMed articles have full text in PMC
- Paywalled articles will be logged in `paywalled_missing.csv`
- This is normal - focus on successful downloads

**Processing Large Lists:**
```python
# For 1000+ PMIDs, process in batches
pmids = load_your_pmids()  # e.g., 2000 PMIDs

batch_size = 100
for i in range(0, len(pmids), batch_size):
    batch = pmids[i:i+batch_size]
    harvester = PMCHarvester(output_dir=f"harvest_batch_{i//batch_size}")
    harvester.harvest(batch, delay=2.0)
```

---

## Using the Extraction Pipeline

See [README_PIPELINE.md](README_PIPELINE.md) for detailed documentation on:
- Extracting individual-level data from articles
- Tiered extraction with cost optimization
- Filtering and post-processing
- Output formats and usage

**Quick Start:**
```bash
python example_usage.py
```

---

## Using the Clinical Data Triage Tool

See [TRIAGE_README.md](TRIAGE_README.md) for detailed documentation on:
- Clinical data processing and triage
- Variant classification
- Quality control and filtering
- Integration with extraction pipeline

**Quick Start:**
```bash
python example_triage.py
```

---

## Troubleshooting

### Harvester Issues

**Problem: "No PMCID found"**
- **Solution**: This is normal. Only ~30% of PubMed articles are in PMC. These will be logged and skipped automatically.

**Problem: "Error fetching supplemental files"**
- **Solution**: Many articles don't have supplements. The script continues and creates markdown from main text.

**Problem: Empty markdown files**
- **Solution**:
  - Check the raw XML at: `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC{PMCID}/fullTextXML`
  - Some PMC articles have parsing issues
  - File a bug report with the PMCID

**Problem: "Email required" error**
- **Solution**: Edit line 32 in `harvest_pmc_fulltext.py` and set your email:
  ```python
  Entrez.email = "your.email@example.com"
  ```

**Problem: Slow processing**
- **Solution**:
  - This is normal - processing is intentionally rate-limited
  - Expect ~1800 articles/hour with default settings
  - Don't reduce delay below 1.0 seconds (respect NCBI servers)

### Installation Issues

**Problem: `pip install -e .` fails**
- **Solution**:
  - Update pip: `pip install --upgrade pip`
  - Install setuptools: `pip install setuptools wheel`
  - Try again

**Problem: Import errors when running**
- **Solution**:
  - Verify virtual environment is activated (you should see `(venv)` in terminal)
  - Reinstall: `pip install -e . --force-reinstall`

### Extraction Issues

**Problem: OpenAI API errors**
- **Solution**:
  - Check that API key is set correctly
  - Verify you have API credits available
  - Check the error message for rate limit issues

**Problem: No data extracted**
- **Solution**:
  - Works best with case reports and detailed clinical studies
  - Limited to abstracts unless using harvested full-text
  - Some articles don't contain individual-level data

---

## Workflow Examples

### Example 1: Fully Automated Workflow (Recommended) 🆕

**Goal**: Get all available variant and patient data for a gene with zero manual work

Simply run:
```bash
python example_automated_workflow.py BRCA1 --email your@email.com
```

This automatically:
1. ✅ Discovers relevant papers from PubMind
2. ✅ Downloads full-text articles from PMC
3. ✅ Extracts variant and patient data
4. ✅ Saves structured JSON results

**Example with custom limits:**
```bash
python example_automated_workflow.py SCN5A \
  --email your@email.com \
  --max-pmids 200 \
  --max-downloads 100
```

### Example 2: Manual Workflow (Step-by-Step)

**Goal**: Get all available data for SCN5A variants with manual control

1. **Discover PMIDs** using PubMind:
   ```python
   from pubmind_fetcher import fetch_pmids_for_gene

   pmids = fetch_pmids_for_gene("SCN5A", email="your@email.com",
                                 output_file="scn5a_pmids.txt")
   ```

2. **Harvest full-text articles:**
   ```bash
   python harvest_pmc_fulltext.py
   ```
   (Edit the script to use your PMID list)

3. **Run extraction pipeline on harvested files:**
   ```bash
   python example_usage.py
   ```

4. **Process and triage the results:**
   ```bash
   python example_triage.py
   ```

### Example 3: Processing a Large PMID List

**Goal**: Process 1000+ articles efficiently

1. Split your PMID list into batches
2. Process each batch:
   ```python
   from harvest_pmc_fulltext import PMCHarvester

   batch_size = 100
   for i in range(0, len(pmids), batch_size):
       batch = pmids[i:i+batch_size]
       harvester = PMCHarvester(output_dir=f"batch_{i//batch_size}")
       harvester.harvest(batch, delay=2.0)
   ```

---

## Additional Resources

- **Full Harvester Documentation**: See `HARVESTER_QUICKSTART.md` and `README_HARVEST.md`
- **Extraction Pipeline Documentation**: See `README_PIPELINE.md`
- **Clinical Triage Documentation**: See `TRIAGE_README.md`
- **NCBI E-utilities**: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **OpenAI API Docs**: https://platform.openai.com/docs

---

## Getting Help

- Check the troubleshooting section above
- Review the detailed guides in the documentation files
- For bugs or feature requests, contact the repository maintainers

---

## Credits

- **NCBI E-utilities**: National Library of Medicine
- **EuropePMC API**: Europe PubMed Central
- **MarkItDown**: Microsoft Research
- **OpenAI**: GPT models for extraction
