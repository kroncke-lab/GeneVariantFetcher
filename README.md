# GeneVariantFetcher

Extracts human genetic variant carriers from biomedical literature (including supplements), classifies them as Affected/Unaffected/Ambiguous, and stores results in SQLite.

## What It Does

GeneVariantFetcher automates the complete pipeline from gene name to SQLite database:

1. **Discovery** - Finds relevant papers using PubMind/PubMed
2. **Download** - Retrieves full-text articles with all supplemental materials
3. **Filter** - Identifies papers with clinical data (Tier 1 keywords + Tier 2 LLM)
4. **Extract** - Uses AI to extract variant and patient data
5. **Storage** - Writes normalized data to SQLite database

**Key Feature**: Downloads full-text articles with supplemental Excel/Word files (70-80% of variant data is in supplements), providing 10-100x more data than abstracts alone.

## Quick Start

### Prerequisites

```bash
# Python 3.11+ required
python --version

# Install dependencies
pip install -e .

# Set required environment variables
export AI_INTEGRATIONS_OPENAI_API_KEY="your-api-key"
export NCBI_EMAIL="your@email.com"
export NCBI_API_KEY="your-key"  # Optional, increases rate limits
```

### Run Full Pipeline

```bash
# Complete workflow: Discovery → Download → Filter → Extract → SQLite
python automated_workflow.py BRCA1 --email your@email.com

# With custom limits
python automated_workflow.py SCN5A --email your@email.com --max-pmids 200 --max-downloads 100
```

**Output**: `automated_output/{GENE}/{TIMESTAMP}/{GENE}.db` (SQLite database)

### Query Database

```bash
# View statistics
python query_variants_db.py automated_output/BRCA1/{timestamp}/BRCA1.db --stats

# Export to CSV
python query_variants_db.py automated_output/BRCA1/{timestamp}/BRCA1.db --export output.csv
```

## Pipeline Stages

| Stage | Module | Description |
|-------|--------|-------------|
| Discovery | `gene_literature/pubmind_fetcher.py` | Fetch PMIDs from PubMind/PubMed |
| Download | `harvesting/orchestrator.py` | Full-text + supplements → `{PMID}_FULL_CONTEXT.md` |
| Filter | `gene_literature/clinical_data_triage.py` | Tier 1 (keywords) + Tier 2 (cheap LLM) |
| Extract | `automated_workflow.py` | LLM extraction with model cascade |
| Storage | `migrate_to_sqlite.py` | Write to SQLite database |

## Key Modules

- `automated_workflow.py` - Main production entry point, orchestrates everything
- `harvesting/orchestrator.py` - PMCHarvester class for downloading full-text
- `gene_literature/clinical_data_triage.py` - Filters papers before expensive extraction
- `migrate_to_sqlite.py` - JSON → SQLite migration
- `query_variants_db.py` - Query the database

## Output Structure

```
automated_output/{GENE}/{TIMESTAMP}/
├── {GENE}_pmids.txt           # Discovered PMIDs
├── pmc_fulltext/
│   └── {PMID}_FULL_CONTEXT.md # Article + supplements
├── extractions/
│   └── {PMID}_extraction.json  # Extracted data
└── {GENE}.db                  # CANONICAL DATABASE
```

## Database Schema

- `papers` - Paper metadata (PMID, title, DOI)
- `individuals` - Case-level data (patient_id, variant, affected_status, phenotype)
- `aggregates` - Cohort statistics (carrier counts, penetrance)
- `variants` - Normalized variant info

## Configuration

### Required Environment Variables

```bash
AI_INTEGRATIONS_OPENAI_API_KEY=your-key  # Required for extraction
NCBI_EMAIL=your@email.com                 # Required for PubMed API
NCBI_API_KEY=your-key                     # Optional, increases rate limits
```

### Command-Line Options

```bash
python automated_workflow.py GENE --email EMAIL [OPTIONS]

Options:
  --max-pmids N      Limit number of PMIDs to discover (default: no limit)
  --max-downloads N  Limit number of papers to download (default: no limit)
```

## Re-running Extraction

To re-run extraction on existing downloads:

```bash
python rerun_extraction.py automated_output/GENE/{timestamp}/
```

## Important Notes

1. **SQLite is canonical** - All data must be written to `{GENE}.db`
2. **Supplements matter** - 70-80% of variant data is in Excel/Word supplements
3. **Production entry point** - Always use `automated_workflow.py`, never `example_*.py` files

## Troubleshooting

**"No PMCID found"**
- Normal: Only ~30% of PubMed articles are in PMC
- These are automatically logged and skipped

**OpenAI API errors**
- Check API key: `echo $AI_INTEGRATIONS_OPENAI_API_KEY`
- Verify you have API credits available

**Slow processing**
- Normal: Rate-limited to respect API servers
- Expected: ~1800 articles/hour with 2s delay

## Development

**Code Style:**
- Python 3.11+, PEP 8, type hints, Google-style docstrings
- Logging via `logger = logging.getLogger(__name__)`
- Imperative commit messages: "Add feature" not "Added feature"

## License

See LICENSE file in repository.
