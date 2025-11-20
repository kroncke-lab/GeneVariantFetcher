# PubMed Central Full-Text Harvester

A standalone Python script to download full-text articles and ALL supplemental materials from PubMed Central, converting everything to unified markdown files for LLM processing.

## Features

- **PMID → PMCID Conversion**: Uses NCBI E-utilities to find PubMed Central IDs
- **Full-Text Download**: Retrieves complete article XML from EuropePMC API
- **Supplemental Files**: Downloads ALL attachments (Excel, Word, PDF, etc.)
- **Unified Markdown**: Converts everything to a single `{PMID}_FULL_CONTEXT.md` file
- **Excel Tables**: Converts spreadsheets to markdown tables
- **Error Handling**: Logs paywalled/missing articles to `paywalled_missing.csv`

## Installation

```bash
pip install -r requirements_harvest.txt
```

## Usage

### Basic Usage

Edit `harvest_pmc_fulltext.py` and modify the `pmids` list in the `main()` function:

```python
pmids = [
    '34931732',
    '35443093', 
    '33442691',
]
```

Then run:

```bash
python harvest_pmc_fulltext.py
```

### Programmatic Usage

```python
from harvest_pmc_fulltext import PMCHarvester

harvester = PMCHarvester(output_dir="pmc_harvest")

pmids = ['12345678', '87654321']
harvester.harvest(pmids, delay=2.0)
```

## Output Structure

```
pmc_harvest/
├── 34931732_FULL_CONTEXT.md          # Unified markdown file
├── 34931732_supplements/              # Supplemental files folder
│   ├── Table_S1.xlsx
│   ├── Figure_S2.pdf
│   └── Methods_S3.docx
├── 35443093_FULL_CONTEXT.md
├── successful_downloads.csv           # Log of successful harvests
└── paywalled_missing.csv              # Log of failed/paywalled articles
```

## Unified Markdown Format

Each `{PMID}_FULL_CONTEXT.md` file contains:

```markdown
# MAIN TEXT

## Article Title

### Abstract
[Abstract content...]

### Introduction
[Introduction content...]

### Methods
[Methods content...]

# SUPPLEMENTAL FILE 1: Table_S1.xlsx

#### Sheet: Patient_Data

| Patient_ID | Variant | Phenotype |
|------------|---------|-----------|
| P001       | c.123G>A | Affected  |

# SUPPLEMENTAL FILE 2: Methods_S3.docx

[Supplemental methods text...]
```

## Configuration

### Set Your Email (Required)

NCBI E-utilities requires an email address. Edit line 27 in `harvest_pmc_fulltext.py`:

```python
Entrez.email = "your.email@example.com"  # Change this!
```

### Adjust Rate Limiting

Change the delay between requests (default 2 seconds):

```python
harvester.harvest(pmids, delay=3.0)  # 3 second delay
```

## API Details

### NCBI E-utilities
- **Purpose**: Convert PMID → PMCID
- **Rate Limit**: 3 requests/second (without API key)
- **Authentication**: None required (but email recommended)

### EuropePMC API
- **Full-Text Endpoint**: `https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML`
- **Supplements Endpoint**: `https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles`
- **Rate Limit**: No official limit, but script uses 2-second delays
- **Authentication**: None required

## Limitations

1. **Paywalled Articles**: Articles without PMCIDs or full-text access are logged and skipped
2. **Large Excel Files**: Only first 100 rows shown in markdown (full file still downloaded)
3. **PDF Conversion**: Requires `markitdown` library; otherwise just notes PDF location
4. **Compressed Archives**: ZIP files are downloaded but not extracted

## Dependencies

- **biopython**: NCBI E-utilities access
- **requests**: HTTP requests for APIs
- **pandas**: Excel to markdown conversion
- **openpyxl**: Excel file reading
- **markitdown**: Advanced file conversion (optional but recommended)
- **python-docx**: Word document reading
- **tabulate**: Markdown table formatting

## Troubleshooting

### "No PMCID found"
- Article is likely behind a paywall or not in PubMed Central
- Check `paywalled_missing.csv` for details

### "Full-text not available"
- Article has a PMCID but full-text XML is not accessible
- May be available only as PDF on PMC website

### "markitdown not available"
- Install with: `pip install markitdown`
- Script will fall back to basic conversion without it

### Rate Limiting
- If you get 429 errors, increase the delay parameter
- Consider getting an NCBI API key for higher limits

## Example Output

Running with 3 PMIDs might produce:

```
Processing PMID: 34931732
  ✓ PMCID: PMC8675309
  ✓ Full-text XML retrieved
  Found 5 supplemental files
    Downloading: Table_S1.xlsx
    Downloading: Table_S2.xlsx
    Downloading: Figure_S1.pdf
  ✅ Created: 34931732_FULL_CONTEXT.md (3 supplements)

============================
Harvest complete!
  ✅ Successful: 2
  ❌ Failed: 1
============================
```

## License

This script is for educational and research purposes. Ensure compliance with PubMed Central and EuropePMC terms of service.
