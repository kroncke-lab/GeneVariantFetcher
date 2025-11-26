# GeneVariantFetcher

Command-line tools for extracting genetic variant and patient data from PubMed literature using AI-powered analysis.

## What It Does

GeneVariantFetcher automates the **complete end-to-end workflow** from gene name to SQLite database:

1. **Discovers papers** - Automatically finds relevant papers for a gene using PubMind/PubMed
2. **Downloads full-text** - Retrieves complete articles with all supplemental materials from PubMed Central
3. **Extracts data** - Uses AI to extract patient-level variants, phenotypes, and clinical data
4. **Aggregates penetrance** - Summarizes variant penetrance data across all papers
5. **Stores in SQLite** - Saves normalized data to a compressed SQLite database for easy querying

**Key Feature**: Downloads full-text articles with supplemental Excel/Word files, not just abstracts, providing 10-100x more data for extraction.

## ⚡ Quick Start (One Command!)

```bash
# Complete workflow: Papers → Full-text → Extraction → SQLite Database
python automated_workflow.py BRCA1 --email your@email.com

# This creates:
# - automated_output/BRCA1/{timestamp}/BRCA1.db  ← SQLite database
# - automated_output/BRCA1/{timestamp}/extractions/  ← JSON files
# - automated_output/BRCA1/{timestamp}/pmc_fulltext/  ← Full-text papers
```

That's it! You now have a complete SQLite database with all variant and penetrance data.

### 🎯 What Changed? (Simplified!)

This codebase has been **dramatically simplified**:
- ✅ **One workflow**: `automated_workflow.py` does everything end-to-end
- ✅ **SQLite integrated**: No manual migration step needed
- ✅ **No duplicates**: Removed redundant code (3 duplicate modules deleted)
- ✅ **Clear structure**: Straightforward pipeline without complex tiering

The old complex tiered filtering system has been archived to `pipeline_tiered_old.py` if you need it.

---

## Quick Start

### Prerequisites

```bash
# Python 3.11+ required
python --version

# Clone repository
git clone https://github.com/kroncke-lab/GeneVariantFetcher.git
cd GeneVariantFetcher

# Install dependencies
pip install -e .

# Set OpenAI API key (required for extraction)
export AI_INTEGRATIONS_OPENAI_API_KEY="your-api-key"
export AI_INTEGRATIONS_OPENAI_BASE_URL="https://api.openai.com/v1"
```

### Automated Workflow (Recommended)

Run the **complete end-to-end workflow** with a single command:

```bash
python automated_workflow.py BRCA1 --email your@email.com
```

This will:
1. Find relevant papers from PubMind
2. Download full-text articles from PMC
3. Extract variant and patient data
4. Aggregate penetrance across papers
5. **Migrate to SQLite database** ← NEW!

**Custom limits:**
```bash
python automated_workflow.py SCN5A \
  --email your@email.com \
  --max-pmids 200 \
  --max-downloads 100
```

**Query your database:**
```bash
python query_variants_db.py automated_output/SCN5A/{timestamp}/SCN5A.db --stats
```

---

## Core Workflows

### 1. Full Automated Pipeline

**Best for**: Getting all available data with minimal effort

```bash
python example_automated_workflow.py BRCA1 --email your@email.com
```

**Output**: `automated_output/BRCA1/{timestamp}/extracted_data.json`

### 2. Step-by-Step Manual Control

**Best for**: Custom PMID lists or detailed control

**Step 1: Get PMIDs**
```python
from pubmind_fetcher import fetch_pmids_for_gene

pmids = fetch_pmids_for_gene("SCN5A", email="your@email.com",
                              output_file="scn5a_pmids.txt")
```

**Step 2: Harvest full-text**
```python
from harvesting import PMCHarvester

harvester = PMCHarvester(output_dir="scn5a_papers")
harvester.harvest(pmids, delay=2.0)
```

**Step 3: Extract data**
```bash
python example_usage.py  # Edit script to point to your harvested papers
```

### 3. Using Your Own PMID List

**Best for**: Working with specific papers or custom searches

```bash
# Create pmids.txt with one PMID per line
# Then edit harvest_pmc_fulltext.py line 333 to load your file

python harvest_pmc_fulltext.py
```

---

## Tools Reference

### PMC Full-Text Harvester

Downloads complete PubMed Central articles with all supplemental materials.

**Features:**
- Converts PMID → PMCID automatically
- Downloads full article text (not just abstracts)
- Retrieves ALL supplemental files (Excel, Word, PDFs)
- Converts Excel tables to markdown for AI processing
- Creates unified `{PMID}_FULL_CONTEXT.md` files

**Usage:**
```python
from harvesting import PMCHarvester

harvester = PMCHarvester(output_dir="pmc_harvest")
pmids = ['35443093', '33442691', '34931732']
harvester.harvest(pmids, delay=2.0)  # 2 second delay between requests
```

**Output Structure:**
```
pmc_harvest/
├── 35443093_FULL_CONTEXT.md          # Unified markdown (article + supplements)
├── 35443093_supplements/              # Raw supplemental files
│   ├── Table_S1.xlsx
│   └── Figure_S2.pdf
├── successful_downloads.csv           # Log of successful downloads
└── paywalled_missing.csv              # Log of unavailable articles
```

**Important**: Only ~30% of PubMed articles are available full-text in PMC. Paywalled articles are automatically logged and skipped.

**Configuration:**
```python
# Required: Set your email (NCBI requirement)
# Edit line 32 in harvest_pmc_fulltext.py
Entrez.email = "your.email@example.com"

# Rate limiting (be polite to APIs)
harvester.harvest(pmids, delay=2.0)  # Recommended: 2-3 seconds
```

### Extraction Pipeline

Extracts structured individual-level variant and phenotype data from articles using AI.

**What it extracts:**
- Genetic variants (HGVS notation: c.1234G>A, p.Arg412Gln)
- Patient phenotypes (with HPO codes when possible)
- Clinical status (affected/unaffected)
- Demographics (age, sex, ancestry)
- Evidence sentences supporting each extraction

**Usage:**
```python
from pipeline import BiomedicalExtractionPipeline

pipeline = BiomedicalExtractionPipeline(email="your@email.com")
results = pipeline.run_for_pmids(['35443093', '33442691'])
```

**Output Format (JSON):**
```json
{
  "pmid": "35443093",
  "individuals": [
    {
      "individual_id": "P1",
      "variants": ["c.1234G>A", "p.Arg412Gln"],
      "phenotypes": ["Cardiomyopathy", "Arrhythmia"],
      "clinical_status": "affected",
      "age": "45",
      "sex": "male",
      "evidence": "Patient P1 presented with..."
    }
  ]
}
```

**Cost**: Uses OpenAI GPT models. Typical cost: $0.01-0.10 per paper depending on length.

### Clinical Data Triage

Filters papers to identify those with original clinical data (case reports, cohorts) vs. reviews or basic science.

**Purpose**: Save time and money by processing only papers with extractable patient data.

**Usage:**
```python
from filters import ClinicalDataTriageFilter

triage = ClinicalDataTriageFilter()
result = triage.triage(
    title="Novel BRCA1 mutation in breast cancer patient",
    abstract="We report a 32-year-old woman...",
    gene="BRCA1"
)

print(result['decision'])    # "KEEP" or "DROP"
print(result['confidence'])  # 0.0-1.0
```

**What gets KEPT:**
- Case reports
- Case series
- Clinical cohorts with patient data
- Functional studies that include patient phenotypes

**What gets DROPPED:**
- Review articles
- Meta-analyses
- Animal studies only
- Cell/in vitro studies only
- Computational predictions only

**Cost**: ~$0.0001 per paper (1 cent per 100 papers)

### PubMind Integration

Automatically discovers relevant papers using PubMind, a literature database of variant-disease associations.

**Usage:**
```python
from pubmind_fetcher import fetch_pmids_for_gene

# Fetch PMIDs for a gene
pmids = fetch_pmids_for_gene(
    gene_symbol="BRCA1",
    email="your@email.com",
    max_results=100,
    output_file="brca1_pmids.txt"
)

print(f"Found {len(pmids)} papers")
```

**How it works:**
1. Queries PubMind database (https://pubmind.wglab.org/)
2. Falls back to enhanced PubMed queries if PubMind unavailable
3. Returns deduplicated PMID list

**Fallback**: If PubMind is inaccessible, automatically uses enhanced PubMed queries optimized for clinical variant papers.

---

## Configuration

### Required Setup

**NCBI Email** (required for PubMed API):
```python
# Edit line 32 in harvest_pmc_fulltext.py
Entrez.email = "your.email@example.com"
```

**OpenAI API Key** (required for extraction):
```bash
export AI_INTEGRATIONS_OPENAI_API_KEY="your-api-key"
export AI_INTEGRATIONS_OPENAI_BASE_URL="https://api.openai.com/v1"
```

Or create `.env` file:
```env
AI_INTEGRATIONS_OPENAI_API_KEY=your-api-key
AI_INTEGRATIONS_OPENAI_BASE_URL=https://api.openai.com/v1
```

### Optional Settings

**Rate Limiting:**
```python
# Harvester delay (seconds between requests)
harvester.harvest(pmids, delay=2.0)  # Recommended: 2-3 seconds

# For large batches (1000+ PMIDs)
harvester.harvest(pmids, delay=3.0)  # Be extra polite
```

**LLM Model Selection:**
```python
# Default: gpt-4o-mini (fast, cheap)
pipeline = BiomedicalExtractionPipeline(model="gpt-4o-mini")

# Higher accuracy (more expensive)
pipeline = BiomedicalExtractionPipeline(model="gpt-4o")
```

---

## Troubleshooting

### Harvester Issues

**"No PMCID found"**
- Normal: Only ~30% of PubMed articles are in PMC
- These are automatically logged to `paywalled_missing.csv`
- Not an error, just unavailable full-text

**"Full-text not available"**
- Article has PMCID but XML not accessible
- Logged and skipped automatically

**"Email required"**
- Set your email in `harvest_pmc_fulltext.py` line 32:
  ```python
  Entrez.email = "your.email@example.com"
  ```

**Slow processing**
- Normal: Rate-limited to respect API servers
- Expected: ~1800 articles/hour with 2s delay
- Don't reduce delay below 1 second

### Extraction Issues

**OpenAI API errors**
- Check API key is set: `echo $AI_INTEGRATIONS_OPENAI_API_KEY`
- Verify you have API credits available
- Check for rate limit messages in error output

**No data extracted**
- Works best with case reports and clinical studies
- Some papers don't contain individual-level data
- Try full-text (harvested papers) instead of abstracts only

**Import errors**
- Verify virtual environment is activated: `(venv)` in prompt
- Reinstall: `pip install -e . --force-reinstall`

### Installation Issues

**`pip install -e .` fails**
- Update pip: `pip install --upgrade pip`
- Install build tools: `pip install setuptools wheel`
- Try again

---

## Output Files

### Harvester Output

```
pmc_harvest/
├── {PMID}_FULL_CONTEXT.md      # Article + supplements as markdown
├── {PMID}_supplements/          # Raw supplemental files
├── successful_downloads.csv     # What worked
└── paywalled_missing.csv        # What didn't work
```

### Extraction Output

```
automated_output/{GENE}/{TIMESTAMP}/
├── extracted_data.json          # Structured extraction results
├── pmc_fulltext/                # Harvested full-text papers
└── pipeline_log.txt             # Processing log
```

### PMID Lists

```
{gene}_pmids.txt                 # One PMID per line
```

---

## Performance Tips

**Large PMID lists (1000+)**: Process in batches of 100
```python
batch_size = 100
for i in range(0, len(pmids), batch_size):
    batch = pmids[i:i+batch_size]
    harvester.harvest(batch, delay=2.0)
```

**Cost optimization**: Use triage filter before extraction
```python
# Triage costs ~$0.0001/paper
# Extraction costs ~$0.05/paper
# Filter first to save money
triage = ClinicalDataTriageFilter()
for paper in papers:
    if triage.triage_paper(paper)['decision'] == 'KEEP':
        pipeline.process_paper(paper)  # Only extract clinical papers
```

**Speed vs. politeness**:
- Fast (1s delay): 3600 articles/hour - use for small batches only
- Recommended (2s delay): 1800 articles/hour - good for most use cases
- Polite (3s delay): 1200 articles/hour - use for large batches (1000+)

---

## API Details

### NCBI E-utilities
- **Purpose**: Convert PMID → PMCID, fetch abstracts
- **Rate Limit**: 3 requests/second (no API key), 10 requests/second (with key)
- **Authentication**: Email required, API key optional

### EuropePMC
- **Purpose**: Download full-text XML and supplemental files
- **Endpoints**:
  - Full-text: `https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML`
  - Supplements: `https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles`
- **Rate Limit**: No official limit (be polite with 2s delays)
- **Authentication**: None required

### OpenAI API
- **Purpose**: AI-powered biomedical data extraction
- **Models**: gpt-4o-mini (recommended), gpt-4o (higher accuracy)
- **Cost**: ~$0.15-1.50 per 1M tokens
- **Authentication**: API key required

---

## Example Workflows

### Workflow 1: Quick Gene Analysis

```bash
# Get all data for a gene in one command
python example_automated_workflow.py KCNQ1 --email your@email.com
```

### Workflow 2: Custom PMID List

```bash
# 1. Create pmids.txt (one PMID per line)
echo "35443093" > my_pmids.txt
echo "33442691" >> my_pmids.txt

# 2. Harvest full-text
python harvest_pmc_fulltext.py  # Edit script to load my_pmids.txt

# 3. Extract data
python example_usage.py  # Edit script to process harvested papers
```

### Workflow 3: Large-Scale Processing

```python
from pubmind_fetcher import fetch_pmids_for_gene
from harvesting import PMCHarvester
from pipeline import BiomedicalExtractionPipeline

# 1. Get PMIDs
pmids = fetch_pmids_for_gene("TTR", email="your@email.com", max_results=500)

# 2. Harvest in batches
batch_size = 100
for i in range(0, len(pmids), batch_size):
    batch = pmids[i:i+batch_size]
    harvester = PMCHarvester(output_dir=f"ttr_batch_{i//batch_size}")
    harvester.harvest(batch, delay=2.0)

# 3. Extract from harvested papers
pipeline = BiomedicalExtractionPipeline(email="your@email.com")
# Process harvested markdown files...
```

---

## Why Full-Text Matters

**Abstracts only (typical approach):**
- Summary statistics
- Limited individual patient data
- No supplemental tables

**Full-text + supplements (this tool):**
- Complete case reports with clinical details
- Supplemental Excel tables with individual patient rows
- Methods sections with cohort descriptions
- Figure legends with additional phenotype data
- **10-100x more extractable data**

**Example**: A paper abstract mentions "15 patients with BRCA1 mutations" but supplemental Table S1 contains individual rows with specific variants, ages, phenotypes for all 15 patients. Harvesting captures this.

---

## Development

**Run tests:**
```bash
python test_triage.py  # Test clinical triage
```

**Project structure:**
```
GeneVariantFetcher/
├── harvest_pmc_fulltext.py      # Full-text harvester
├── pipeline.py                  # Extraction pipeline
├── clinical_data_triage.py      # Clinical data triage
├── pubmind_fetcher.py           # PubMind integration
├── filters.py                   # Filtering utilities
├── sourcer.py                   # Paper sourcing utilities
├── example_*.py                 # Example scripts
└── README.md                    # This file
```

---

## Support

- **Troubleshooting**: See sections above
- **NCBI E-utilities docs**: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **EuropePMC API**: https://europepmc.org/RestfulWebService
- **OpenAI API**: https://platform.openai.com/docs
- **PubMind**: https://pubmind.wglab.org/

---

## License

See LICENSE file in repository.

---

## Credits

- **NCBI E-utilities**: National Library of Medicine
- **EuropePMC API**: Europe PubMed Central
- **PubMind**: Wang Genomics Lab
- **OpenAI**: GPT models for extraction
