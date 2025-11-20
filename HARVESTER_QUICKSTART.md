# PMC Full-Text Harvester - Quick Start Guide

## What This Does

Converts a list of PMIDs into comprehensive markdown files containing:
- **Full-text article** (not just abstracts!)
- **ALL supplemental materials** (Excel tables, Word docs, PDFs)
- **Everything unified** into one `{PMID}_FULL_CONTEXT.md` file per article

Perfect for preparing rich datasets for LLM extraction of patient-level data.

---

## Installation

```bash
pip install -r requirements_harvest.txt
```

**Required packages:**
- biopython (NCBI E-utilities)
- requests (HTTP client)
- markitdown (file conversion)
- pandas (Excel tables)
- openpyxl (Excel reading)
- python-docx (Word docs)
- tabulate (markdown tables)

---

## Basic Usage

### Option 1: Edit the Script Directly

Open `harvest_pmc_fulltext.py` and modify line 333:

```python
pmids = [
    '35443093',   # Your PMID here
    '33442691',   # Another PMID
    '34931732',   # And another
]
```

Then run:

```bash
python harvest_pmc_fulltext.py
```

### Option 2: Use the Example Script

```bash
python example_harvest_from_pubmind.py
```

Choose from:
1. Load PMIDs from text file (`example_pmids.txt`)
2. Load PMIDs from PubMind CSV export
3. Manually specify PMIDs
4. Harvest + prepare for LLM extraction

### Option 3: Python API

```python
from harvest_pmc_fulltext import PMCHarvester

harvester = PMCHarvester(output_dir="my_harvest")
pmids = ['12345678', '87654321']
harvester.harvest(pmids, delay=2.0)
```

---

## Integration with PubMind App

### Step 1: Export PMIDs from PubMind

In your Streamlit app, search for a gene (e.g., BRCA2) and click **"Download PMIDs (TXT)"** button.

This creates a file like `BRCA2_pmids.txt` with PMIDs.

### Step 2: Harvest Full-Text Papers

```python
from harvest_pmc_fulltext import PMCHarvester

with open('BRCA2_pmids.txt', 'r') as f:
    pmids = [line.strip() for line in f if line.strip()]

harvester = PMCHarvester(output_dir="harvest_brca2")
harvester.harvest(pmids[:50], delay=2.0)  # Start with first 50
```

### Step 3: Review Output

Check `harvest_brca2/` folder:
- `{PMID}_FULL_CONTEXT.md` files ready for LLM extraction
- `successful_downloads.csv` - which ones worked
- `paywalled_missing.csv` - which ones are unavailable

---

## Output Format

Each `{PMID}_FULL_CONTEXT.md` file contains:

```markdown
# MAIN TEXT

## Article Title Here

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
| P001       | BRCA2 | c.1234G>A | 45 | Breast cancer |
| P002       | BRCA2 | c.5678C>T | 52 | Ovarian cancer |

# SUPPLEMENTAL FILE 2: Extended_Methods.docx

Detailed supplemental methods...
```

This is **exactly** what you need for comprehensive LLM extraction!

---

## Important Configuration

### Set Your Email (Required)

NCBI E-utilities requires an email. Edit line 32 in `harvest_pmc_fulltext.py`:

```python
Entrez.email = "your.email@example.com"  # Change this!
```

### Adjust Rate Limiting

EuropePMC has no official rate limit, but be polite:

```python
harvester.harvest(pmids, delay=3.0)  # 3 seconds between requests
```

---

## What Gets Logged

### Successful Downloads
`successful_downloads.csv`:
```csv
PMID,PMCID,Supplements_Downloaded
35443093,PMC10078423,3
33442691,PMC7805448,0
```

### Failed/Paywalled
`paywalled_missing.csv`:
```csv
PMID,Reason
34931732,No PMCID found
12345678,Full-text not available
```

---

## Common Issues

### "No PMCID found"
- Article is paywalled or not in PubMed Central
- Only ~30% of PubMed articles have full-text in PMC
- These are logged and skipped automatically

### "Full-text not available"
- Article has a PMCID but XML not accessible via EuropePMC
- May only be available as PDF on PMC website
- Logged and skipped

### "Error fetching supplemental files"
- Many articles don't have supplemental files
- This is normal and not an error condition
- Script continues and creates markdown from main text

### Empty Markdown Files
- Rare, but can happen if XML parsing fails
- Check the raw XML at: `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC{ID}/fullTextXML`

---

## Performance Tips

### For Large PMID Lists (1000+)

Process in batches:

```python
pmids = load_your_pmids()  # e.g., 2000 PMIDs

batch_size = 100
for i in range(0, len(pmids), batch_size):
    batch = pmids[i:i+batch_size]
    
    harvester = PMCHarvester(output_dir=f"harvest_batch_{i//batch_size}")
    harvester.harvest(batch, delay=2.0)
    
    print(f"Completed batch {i//batch_size}")
```

### Speed vs. Politeness

- **Fast** (1s delay): ~3600 articles/hour (use for small batches)
- **Moderate** (2s delay): ~1800 articles/hour (recommended)
- **Polite** (3s delay): ~1200 articles/hour (use for large batches)

---

## Why Full-Text > Abstracts?

### Abstracts Only (Current PubMind Extraction):
- Summary statistics
- Limited individual patient data
- No supplemental tables

### Full-Text + Supplements (This Harvester):
- **Complete case reports** with full clinical details
- **Supplemental tables** with individual patient rows
- **Methods sections** with cohort descriptions
- **Figure legends** with additional phenotype data
- **10-100x more data** for LLM extraction

---

## Example Workflow

```bash
# 1. Search BRCA2 in PubMind app
# 2. Download PMIDs → BRCA2_pmids.txt

# 3. Harvest full-text papers
python harvest_pmc_fulltext.py

# 4. Review harvested papers
ls pmc_harvest/
# → 35443093_FULL_CONTEXT.md (56 KB!)
# → 33442691_FULL_CONTEXT.md (657 bytes)
# → successful_downloads.csv
# → paywalled_missing.csv

# 5. Use with your LLM extraction pipeline
# Now you have FULL CONTEXT instead of just abstracts!
```

---

## Next Steps

1. **Test with small sample**: Try 5-10 PMIDs first
2. **Review output quality**: Check the markdown files
3. **Scale up**: Process larger PMID lists in batches
4. **LLM Extraction**: Feed the markdown into OpenAI/GPT-5
5. **Structured Data**: Extract patient-level variant/phenotype data

---

## Support Files

- `harvest_pmc_fulltext.py` - Main harvester script
- `example_harvest_from_pubmind.py` - Integration examples
- `requirements_harvest.txt` - Python dependencies
- `README_HARVEST.md` - Detailed documentation
- `HARVESTER_QUICKSTART.md` - This file

---

## Questions?

Check `README_HARVEST.md` for:
- API details
- Advanced configuration
- Troubleshooting
- Full documentation
