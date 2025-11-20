# GeneVariantFetcher - Usage Guide

A powerful tool for searching genetic variants and extracting patient-level data from PubMed literature. This application combines the PubMind-DB variant database with AI-powered text extraction to help researchers analyze gene variants and associated publications.

## Table of Contents

- [What This Tool Does](#what-this-tool-does)
- [Quick Start (Browser - Replit)](#quick-start-browser---replit)
- [Local Setup (PyCharm/VSCode)](#local-setup-pycharmintellijvscode)
- [Using the Web Interface](#using-the-web-interface)
- [Using the Harvesting Tool](#using-the-harvesting-tool)
- [Features Overview](#features-overview)
- [Troubleshooting](#troubleshooting)

---

## What This Tool Does

**GeneVariantFetcher** provides two main capabilities:

1. **Web Interface (Streamlit App)** - Search for genetic variants by gene name and:
   - View variant details with pathogenicity classifications
   - See associated PubMed publications
   - Extract individual patient-level data using AI
   - Download results as CSV/JSON files

2. **Full-Text Harvester** - Download complete PubMed Central articles with:
   - Full-text content (not just abstracts)
   - All supplemental materials (Excel, Word, PDFs)
   - Everything consolidated into markdown files for LLM processing

---

## Quick Start (Browser - Replit)

### Option 1: Running on Replit (Easiest)

This is the **simplest way** to use the tool - no installation needed!

1. **Open the project on Replit**
   - If you're viewing this in Replit, you're already there!
   - Look for the green "Run" button at the top of the page

2. **Click the "Run" button**
   - The application will start automatically
   - Wait 10-30 seconds for the server to start
   - You'll see a web browser preview open on the right side

3. **Start using the app**
   - The web interface will appear automatically
   - Skip to [Using the Web Interface](#using-the-web-interface) below

**That's it!** The browser interface is now ready to use.

---

## Local Setup (PyCharm/IntelliJ/VSCode)

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

The individual-level data extraction feature requires an OpenAI API key.

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

Then modify `app.py` to load the .env file (add at the top):

```python
from dotenv import load_dotenv
load_dotenv()
```

And install python-dotenv:

```bash
pip install python-dotenv
```

### Step 5: Run the Application

**In PyCharm:**

1. Open the project folder in PyCharm
2. Right-click on `app.py`
3. Select "Run 'app'"
4. Or use the terminal:

```bash
streamlit run app.py --server.port 5000
```

**In VSCode:**

1. Open the project folder
2. Open a terminal (View → Terminal)
3. Run:

```bash
streamlit run app.py --server.port 5000
```

**In the terminal:**

```bash
streamlit run app.py --server.port 5000
```

4. **Open your browser**
   - The app will automatically open at `http://localhost:5000`
   - If not, manually navigate to that URL

---

## Using the Web Interface

### Step-by-Step: Searching for Variants

1. **Enter a Gene Symbol**
   - Look at the left sidebar
   - In the text box labeled "Enter Gene Symbol"
   - Type a gene name (e.g., `BRCA1`, `TP53`, `SCN5A`, `KCNQ1`)
   - Gene symbols are case-insensitive

2. **Choose Search Mode (Optional)**
   - **Single Gene**: Search one gene at a time (default)
   - **Multiple Genes**: Enter multiple genes, one per line

3. **Set Filters (Optional)**
   - **Min PMIDs per variant**: Show only variants with at least X publications
   - Leave at 0 to see all variants

4. **Click the "Search" Button**
   - Wait a few seconds for results
   - You'll see a progress indicator

### Understanding the Results

The results are organized into 4 tabs:

#### Tab 1: 📊 Visualizations

- **Left chart**: Distribution of how many PMIDs each variant has
- **Right chart**: Top 15 most-studied variants
- Use these to identify well-researched variants

#### Tab 2: 📋 Variant Details

- **Table view** with columns:
  - `Gene`: The gene symbol you searched
  - `Variant_ID`: PubMind unique identifier (PVID)
  - `HGVS`: DNA-level variant notation (e.g., c.1234G>A)
  - `HGVS_Protein`: Protein-level change (e.g., p.Arg412His)
  - `rsID`: dbSNP identifier (clickable link)
  - `Pathogenicity`: AI-curated classification (Pathogenic/Benign/VUS)
  - `PMID_Count`: Number of PubMed publications mentioning this variant
  - `First 5 PMIDs`: Quick links to publications

- **Features:**
  - Click on PVID links to view full details on PubMind-DB
  - Click on rsID links to view in dbSNP database
  - Expand "View LLM Reasoning" to see AI explanations for pathogenicity
  - Click PMID numbers to open articles on PubMed

- **Download:**
  - Click "Download Variant Data (CSV)" at the bottom
  - Opens in Excel, Google Sheets, etc.

#### Tab 3: 📚 PMID List

- **Full list** of all unique PubMed IDs across all variants
- **Quick links** to the first 50 PMIDs
- **Download:**
  - Click "Download PMIDs (TXT)" to get a plain text file
  - One PMID per line
  - Use this for the Harvesting Tool (see below)

#### Tab 4: 🧬 Individual Extractions

**What this does:**
Analyzes PubMed articles using AI to extract structured patient-level data including:
- Individual patients with specific variants
- Phenotypes (symptoms, conditions)
- Clinical status (affected/unaffected)
- Age, sex, and other demographics

**How to use:**

1. **Choose number of articles**
   - Use the number input (default: 5)
   - Start small (5-10) to test
   - Maximum: 20 articles per run
   - **Note:** Uses OpenAI API credits (~$0.02 per article)

2. **Click "🚀 Extract Individual Data"**
   - Wait for processing (30-60 seconds per article)
   - Progress indicator shows status

3. **Review Results**
   - Table shows extracted individuals
   - Columns: individual_id, PMID, gene, age, sex, affected_status, phenotypes, variants
   - View extraction statistics in expandable section

4. **Download Results**
   - **JSON format**: For programmatic analysis
   - **CSV format**: For Excel/Google Sheets

**Important Notes:**
- Works best with case reports and detailed clinical studies
- Limited to abstracts (full text requires harvesting - see below)
- Requires valid OpenAI API key
- Costs are typically $0.02-0.05 per article

---

## Using the Harvesting Tool

The harvesting tool downloads **full-text** PubMed Central articles with all supplemental materials. This provides much richer data than abstracts alone.

### When to Use the Harvester

- You need full article text, not just abstracts
- You want supplemental tables with patient-level data
- You're preparing articles for detailed LLM extraction
- You have a list of PMIDs from the web interface

### Option 1: Quick Start (Edit the Script)

1. **Get your PMIDs**
   - From the web interface: Tab 3 → Download PMIDs (TXT)
   - Or manually create a list

2. **Edit the harvesting script**
   - Open `harvest_pmc_fulltext.py` in PyCharm/VSCode
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

   **In PyCharm:**
   - Right-click `harvest_pmc_fulltext.py`
   - Click "Run 'harvest_pmc_fulltext'"

   **In terminal:**
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

# Load PMIDs from the web app download
with open('BRCA1_pmids.txt', 'r') as f:
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

## Features Overview

### Web Interface Features

| Feature | Description | Use Case |
|---------|-------------|----------|
| **Single Gene Search** | Search one gene at a time | Quick lookups, focused research |
| **Multi-Gene Search** | Search multiple genes | Comparative studies, gene panels |
| **Pathogenicity AI** | AI-curated variant classifications | Clinical interpretation |
| **Publication Links** | Direct links to PubMed articles | Literature review |
| **Individual Extraction** | AI-powered patient data extraction | Building databases, meta-analysis |
| **Data Export** | Download as CSV/JSON | Further analysis, Excel processing |
| **Visualizations** | Charts and graphs | Presentations, quick insights |

### Harvester Features

| Feature | Description | Benefit |
|---------|-------------|---------|
| **Full-Text Download** | Complete article text | 10-100x more data than abstracts |
| **Supplemental Files** | Excel, Word, PDF supplements | Patient-level data tables |
| **Markdown Conversion** | Everything in one .md file | LLM-ready format |
| **Batch Processing** | Process hundreds of PMIDs | Automated large-scale extraction |
| **Error Handling** | Automatic retry, logging | Robust processing |

---

## Troubleshooting

### Web Interface Issues

**Problem: "API Error" when searching**
- **Solution**: Check your internet connection. The PubMind API might be temporarily unavailable. Try again in a few minutes.

**Problem: "No variants found"**
- **Solution**:
  - Check spelling of gene symbol
  - Try alternative gene names (e.g., "SCN5A" vs "Sodium Channel")
  - Not all genes have data in PubMind-DB

**Problem: Individual extraction fails**
- **Solution**:
  - Check that OpenAI API key is set correctly
  - Verify you have API credits available
  - Check the error message for rate limit issues
  - Try with fewer articles (start with 1-2)

**Problem: App won't start in PyCharm**
- **Solution**:
  - Verify Python 3.11+ is installed: `python --version`
  - Reinstall dependencies: `pip install -e .`
  - Try running from terminal: `streamlit run app.py`
  - Check if port 5000 is already in use

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

**Problem: Port already in use**
- **Solution**:
  - Use a different port: `streamlit run app.py --server.port 8501`
  - Or kill the process using port 5000

---

## Workflow Examples

### Example 1: Clinical Variant Review

**Goal**: Review pathogenicity of BRCA1 variants

1. Open the web interface
2. Search for "BRCA1"
3. Set "Min PMIDs" to 5 (focus on well-studied variants)
4. Click Search
5. Go to Tab 2: Variant Details
6. Review pathogenicity classifications
7. Expand "View LLM Reasoning" for AI explanations
8. Click on rsID links to verify in dbSNP
9. Download CSV for reporting

### Example 2: Building a Patient Database

**Goal**: Extract individual patient data from SCN5A papers

1. Open the web interface
2. Search for "SCN5A"
3. Go to Tab 3: PMID List
4. Download PMIDs (TXT)
5. Go to Tab 4: Individual Extractions
6. Set to 10 articles
7. Click "Extract Individual Data"
8. Wait for processing
9. Download JSON results
10. Use the harvester to get full-text for top articles:
    ```bash
    # Edit harvest_pmc_fulltext.py with PMIDs
    python harvest_pmc_fulltext.py
    ```
11. Feed full-text markdown files to your own LLM for deeper extraction

### Example 3: Multi-Gene Panel Analysis

**Goal**: Compare variants across a cardiac gene panel

1. Open the web interface
2. Select "Multiple Genes" mode
3. Enter genes (one per line):
   ```
   SCN5A
   KCNQ1
   KCNH2
   SCN1B
   ```
4. Click Search
5. Go to Tab 1: Visualizations
6. Compare distribution charts
7. Go to Tab 2: Download CSV
8. Analyze in Excel/Python

---

## Additional Resources

- **Full Harvester Documentation**: See `HARVESTER_QUICKSTART.md`
- **Detailed Harvester Guide**: See `README_HARVEST.md`
- **PubMind-DB Website**: https://pubmind.wglab.org
- **NCBI E-utilities**: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **OpenAI API Docs**: https://platform.openai.com/docs

---

## Getting Help

- Check the troubleshooting section above
- Review the detailed guides in `HARVESTER_QUICKSTART.md`
- For bugs or feature requests, contact the repository maintainers

---

## Credits

- **PubMind-DB**: Wang Genomics Lab
- **NCBI E-utilities**: National Library of Medicine
- **MarkItDown**: Microsoft Research
- **Streamlit**: Streamlit Inc.
