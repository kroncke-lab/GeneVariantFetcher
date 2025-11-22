# PubMind Integration Guide

## Overview

GeneVariantFetcher now integrates with **PubMind** ([https://pubmind.wglab.org/](https://pubmind.wglab.org/)), a literature-derived knowledgebase of variant-disease-pathogenicity associations. PubMind uses LLM-assisted extraction to identify papers containing genetic variant data, making it ideal for automatically discovering papers with patient-level variant information.

## What is PubMind?

PubMind is a large language model (LLM)-assisted framework for Publication Mutation and Information Discovery developed by the Wang Genomics Lab. It:

- **Curates variant-disease associations** from biomedical literature
- **Extracts pathogenicity relationships** using AI
- **Links variants to PMIDs/PMCIDs** for original publications
- **Provides comprehensive coverage** of genetic variants in literature

This makes PubMind perfect for finding papers that contain the detailed variant and patient data you need.

## Features

✅ **Automatic PMID discovery** - No need to manually search PubMed
✅ **Variant-focused results** - Papers are pre-filtered for genetic variant content
✅ **Fallback to PubMed** - Automatically uses enhanced PubMed queries if PubMind is unavailable
✅ **Integration with existing pipeline** - Works seamlessly with the harvester and extraction tools

## Installation

The PubMind integration is included in GeneVariantFetcher. Install dependencies:

```bash
pip install -e .
```

Required packages:
- `requests` - For web scraping
- `beautifulsoup4` - For HTML parsing
- `biopython` - For PubMed fallback queries

## Quick Start

### Option 1: Standalone PMID Fetching

Fetch PMIDs for a gene and save to file:

```python
from pubmind_fetcher import fetch_pmids_for_gene

# Fetch PMIDs for BRCA1
pmids = fetch_pmids_for_gene(
    gene_symbol="BRCA1",
    email="your@email.com",
    max_results=100,
    output_file="BRCA1_pmids.txt"
)

print(f"Found {len(pmids)} papers for BRCA1")
```

### Option 2: Integrated Pipeline

The PubMind fetcher is automatically used by the sourcer:

```python
from sourcer import PaperSourcer

sourcer = PaperSourcer(email="your@email.com")

# This now queries PubMind + PubMed + EuropePMC
pmids = sourcer.fetch_papers(
    gene_symbol="SCN5A",
    use_pubmind=True  # Default: True
)

print(f"Found {len(pmids)} deduplicated PMIDs")
```

### Option 3: Fully Automated Workflow

Run the complete workflow from gene to extracted data:

```bash
python example_automated_workflow.py BRCA1 --email your@email.com
```

This will:
1. Fetch relevant PMIDs from PubMind
2. Download full-text papers from PMC
3. Extract variant and patient data
4. Save structured JSON results

## Usage Examples

### Example 1: Basic Gene Query

```python
from pubmind_fetcher import PubMindFetcher

fetcher = PubMindFetcher(email="your@email.com")

# Fetch PMIDs for a gene
pmids = fetcher.fetch_pmids_for_gene("TP53", max_results=200)

# Save to file
fetcher.save_pmids_to_file(pmids, "TP53_papers.txt")
```

### Example 2: Variant-Specific Query

```python
from pubmind_fetcher import fetch_pmids_for_variant

# Find papers discussing a specific variant
pmids = fetch_pmids_for_variant(
    variant="c.5946delT",
    gene_symbol="BRCA2",
    email="your@email.com"
)

print(f"Found {len(pmids)} papers mentioning c.5946delT in BRCA2")
```

### Example 3: Custom Workflow

```python
from pubmind_fetcher import PubMindFetcher
from harvest_pmc_fulltext import PMCHarvester
from pipeline import BiomedicalExtractionPipeline

# Step 1: Get PMIDs
fetcher = PubMindFetcher(email="your@email.com")
pmids = fetcher.fetch_pmids_for_gene("KCNQ1", max_results=50)

# Step 2: Download papers
harvester = PMCHarvester(output_dir="kcnq1_papers")
harvester.harvest(pmids, delay=2.0)

# Step 3: Extract data
pipeline = BiomedicalExtractionPipeline(email="your@email.com")
results = pipeline.run_for_pmids(pmids)
```

## Command-Line Usage

### Fetch PMIDs for a Gene

```bash
python pubmind_fetcher.py BRCA1
```

This will:
- Prompt for your email
- Fetch up to 50 PMIDs for BRCA1
- Display results
- Offer to save to file

### Run Full Automated Workflow

```bash
# Basic usage
python example_automated_workflow.py SCN5A --email your@email.com

# With custom limits
python example_automated_workflow.py BRCA1 \
  --email your@email.com \
  --max-pmids 200 \
  --max-downloads 100 \
  --output my_results

# Quick test (small dataset)
python example_automated_workflow.py TP53 \
  --email your@email.com \
  --max-pmids 10 \
  --max-downloads 5
```

## How It Works

### 1. PubMind Scraping (Primary Method)

The fetcher attempts to query the PubMind web interface:

```python
# Constructs search URL
https://pubmind.wglab.org/search?gene=BRCA1&type=gene

# Parses HTML to extract PMIDs from:
- Links to PubMed (pubmed.ncbi.nlm.nih.gov/12345678)
- Text containing "PMID: 12345678"
- Table cells with PMIDs
```

### 2. PubMed Fallback (Secondary Method)

If PubMind is unavailable or returns few results, uses enhanced PubMed queries:

```python
# Enhanced query optimized for clinical variant papers
query = """
BRCA1[Gene] AND (
    (mutation OR variant OR variants) AND
    (patient OR case OR cohort OR clinical OR phenotype)
)
"""
```

### 3. Deduplication

All PMIDs from multiple sources are deduplicated:

```python
# Combines results from:
- PubMind
- PubMed
- EuropePMC

# Returns unique, sorted list
```

## Configuration Options

### PubMindFetcher Constructor

```python
PubMindFetcher(
    email="your@email.com",  # Required for PubMed fallback
    use_fallback=True        # Use PubMed if PubMind fails
)
```

### Fetch Methods

```python
# For genes
fetch_pmids_for_gene(
    gene_symbol="BRCA1",
    max_results=1000,    # Maximum PMIDs to return
    delay=1.0            # Delay between requests (seconds)
)

# For variants
fetch_pmids_for_variant(
    variant="c.1234G>A",
    gene_symbol="TP53",  # Optional but recommended
    max_results=500,
    delay=1.0
)
```

## Troubleshooting

### Issue: PubMind returns no results

**Possible causes:**
- Site access is restricted (403 error)
- Site structure has changed
- Gene symbol is not recognized

**Solution:**
The fetcher automatically falls back to PubMed with enhanced queries. You can check logs to see which source was used:

```python
import logging
logging.basicConfig(level=logging.INFO)

# Check logs for:
# "PubMind returned no results, using PubMed fallback"
```

### Issue: Too many PMIDs returned

**Solution:**
Limit results with the `max_results` parameter:

```python
pmids = fetcher.fetch_pmids_for_gene("BRCA1", max_results=50)
```

### Issue: Slow performance

**Solution:**
- Reduce `max_results`
- Increase `delay` between requests
- Use `use_fallback=False` to skip PubMed queries

```python
fetcher = PubMindFetcher(
    email="your@email.com",
    use_fallback=False  # Skip fallback for speed
)
```

### Issue: "Access denied (403)"

**Cause:** PubMind website may require authentication or block automated access

**Solution:** The fallback to PubMed will activate automatically:

```python
# Fallback activates automatically
pmids = fetcher.fetch_pmids_for_gene("SCN5A")
# Uses PubMed enhanced queries instead
```

## Advanced Usage

### Custom Fallback Queries

You can customize the PubMed fallback query:

```python
class CustomPubMindFetcher(PubMindFetcher):
    def _fetch_from_pubmed_gene(self, gene_symbol, max_results):
        # Your custom query here
        query = f'{gene_symbol}[Gene] AND "case report"[Publication Type]'
        return self._query_pubmed(query, max_results)
```

### Combining Multiple Genes

```python
from pubmind_fetcher import PubMindFetcher

fetcher = PubMindFetcher(email="your@email.com")
all_pmids = set()

for gene in ["BRCA1", "BRCA2", "TP53"]:
    pmids = fetcher.fetch_pmids_for_gene(gene, max_results=100)
    all_pmids.update(pmids)

print(f"Total unique PMIDs across all genes: {len(all_pmids)}")
```

### Filtering Results

```python
from pubmind_fetcher import fetch_pmids_for_gene
from sourcer import PaperSourcer

# Get PMIDs
pmids = fetch_pmids_for_gene("KCNQ1", email="your@email.com")

# Filter by checking metadata
sourcer = PaperSourcer(email="your@email.com")
clinical_pmids = []

for pmid in pmids:
    paper = sourcer.fetch_paper_metadata(pmid)
    if paper and paper.abstract:
        # Filter for case reports or clinical studies
        if "case report" in paper.abstract.lower() or "patient" in paper.abstract.lower():
            clinical_pmids.append(pmid)

print(f"Found {len(clinical_pmids)} clinical papers out of {len(pmids)}")
```

## Integration with Existing Tools

### With PMC Harvester

```python
from pubmind_fetcher import fetch_pmids_for_gene
from harvest_pmc_fulltext import PMCHarvester

# Fetch PMIDs
pmids = fetch_pmids_for_gene("SCN5A", email="your@email.com", max_results=100)

# Download papers
harvester = PMCHarvester(output_dir="scn5a_harvest")
harvester.harvest(pmids, delay=2.0)
```

### With Extraction Pipeline

```python
from pubmind_fetcher import fetch_pmids_for_gene
from pipeline import run_pipeline_for_gene

# Option 1: Use PubMind manually
pmids = fetch_pmids_for_gene("BRCA1", email="your@email.com")
# Then process with pipeline...

# Option 2: Pipeline automatically uses PubMind
results = run_pipeline_for_gene(
    gene_symbol="BRCA1",
    email="your@email.com",
    max_papers=50
)
# PubMind is queried automatically via sourcer
```

### With Clinical Triage

```python
from pubmind_fetcher import fetch_pmids_for_gene
from filters import ClinicalDataTriageFilter
from sourcer import PaperSourcer

# Fetch PMIDs
pmids = fetch_pmids_for_gene("TP53", email="your@email.com")

# Triage to find clinical papers
triage = ClinicalDataTriageFilter()
sourcer = PaperSourcer(email="your@email.com")

keep_pmids = []
for pmid in pmids:
    paper = sourcer.fetch_paper_metadata(pmid)
    if paper:
        result = triage.triage_paper(paper, gene="TP53")
        if result['decision'] == 'KEEP':
            keep_pmids.append(pmid)

print(f"Kept {len(keep_pmids)} out of {len(pmids)} papers")
```

## API Reference

### `PubMindFetcher`

Main class for fetching PMIDs from PubMind.

**Constructor:**
```python
PubMindFetcher(
    email: str = "your.email@example.com",
    use_fallback: bool = True
)
```

**Methods:**

#### `fetch_pmids_for_gene(gene_symbol, max_results=1000, delay=1.0)`
Fetch PMIDs for papers about a gene.

**Returns:** `List[str]` - List of PMIDs

#### `fetch_pmids_for_variant(variant, gene_symbol=None, max_results=500, delay=1.0)`
Fetch PMIDs for papers about a specific variant.

**Returns:** `List[str]` - List of PMIDs

#### `save_pmids_to_file(pmids, output_file)`
Save PMIDs to a text file (one per line).

### Convenience Functions

#### `fetch_pmids_for_gene(gene_symbol, email, max_results=1000, output_file=None)`
Quick function to fetch PMIDs for a gene.

#### `fetch_pmids_for_variant(variant, gene_symbol=None, email, max_results=500, output_file=None)`
Quick function to fetch PMIDs for a variant.

## Best Practices

1. **Always provide an email** - Required by NCBI for PubMed fallback
2. **Use reasonable limits** - Start with 50-100 max_results for testing
3. **Respect rate limits** - Use delay >= 1.0 second between requests
4. **Save PMID lists** - Always save to file for reproducibility
5. **Enable fallback** - Keep `use_fallback=True` for reliability
6. **Check logs** - Monitor which source provided results

## Performance Tips

- **For large queries (>500 PMIDs):** Increase delay to 2-3 seconds
- **For quick tests:** Use `max_results=10` and `use_fallback=False`
- **For comprehensive coverage:** Use all sources (PubMind + PubMed + EuropePMC)
- **For variant-specific:** Use `fetch_pmids_for_variant()` instead of gene queries

## Known Limitations

1. **PubMind access** - Site may restrict automated access (falls back to PubMed)
2. **HTML parsing** - Site structure changes may affect scraping (fallback handles this)
3. **Result limits** - PubMind may have internal limits on results returned
4. **No authentication** - Current implementation doesn't support PubMind API authentication

## Future Enhancements

Planned improvements:
- [ ] Official PubMind API support (when available)
- [ ] Better variant notation normalization
- [ ] Caching of PMID lists to avoid re-fetching
- [ ] Integration with variant databases (ClinVar, gnomAD)
- [ ] Batch processing for multiple genes

## Support and Feedback

For issues, questions, or suggestions:
- Check the troubleshooting section above
- Review logs for error messages
- File an issue on the GitHub repository

## References

- **PubMind Database:** https://pubmind.wglab.org/
- **Wang Genomics Lab:** https://github.com/WGLab
- **PubMed E-utilities:** https://www.ncbi.nlm.nih.gov/books/NBK25501/

## License

Part of the GeneVariantFetcher project. See main LICENSE file.

---

**Note:** PubMind is developed and maintained by the Wang Genomics Lab. This integration is an independent implementation and is not officially affiliated with or endorsed by PubMind.
