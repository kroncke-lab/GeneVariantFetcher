# GeneVariantFetcher

A unified workflow for extracting human genetic variant carriers from biomedical literature.

## Mission

**Goal:** Extract all human carriers of genetic variants from biomedical literature (including full text and supplements).

**Core Criteria:** Only humans with clinical phenotyping that allows classification as:
- **Affected** - Manifests the disease phenotype
- **Unaffected** - Carrier without disease manifestation
- **Ambiguous** - Insufficient clinical data for classification

**Output:** A single harmonized SQLite database of variant carriers per gene.

---

## What It Does

GeneVariantFetcher automates the complete pipeline from gene name to SQLite database:

1. **Discovery** - Finds variant-specific papers via PubMind
2. **Download** - Retrieves full-text articles with all supplemental materials (Excel, Word)
3. **Extract** - Uses LLM cascade to extract variant and patient data (parallel processing)
4. **Aggregate** - Validates penetrance data and calculates cross-paper statistics
5. **Storage** - Writes normalized data to SQLite database

**Key Insight:** 70-80% of variant data is buried in supplemental Excel/Word files. This tool extracts it all, providing 10-100x more data than abstracts alone.

---

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
# Complete workflow: Discovery → Download → Extract → Aggregate → SQLite
python automated_workflow.py BRCA1 --email your@email.com

# With limits
python automated_workflow.py SCN5A --email your@email.com --max-pmids 200 --max-downloads 100
```

**Output:** `automated_output/{GENE}/{TIMESTAMP}/{GENE}.db` (SQLite database)

---

## Pipeline Architecture

```
Discovery → Download → Extract → Aggregate → SQLite
```

| Stage | Module | Description |
|-------|--------|-------------|
| Discovery | `gene_literature/pubmind_fetcher.py` | Fetch PMIDs from PubMind (variant-specific literature) |
| Download | `harvesting/orchestrator.py` | Full-text + supplements → `{PMID}_FULL_CONTEXT.md` |
| Extract | `pipeline/extraction.py` | LLM extraction with model cascade |
| Aggregate | `pipeline/aggregation.py` | Penetrance validation and statistics |
| Storage | `harvesting/migrate_to_sqlite.py` | JSON → normalized SQLite database |

---

## Output Structure

```
automated_output/{GENE}/{TIMESTAMP}/
├── {GENE}_pmids.txt                # Discovered PMIDs
├── {GENE}_workflow_summary.json    # Execution statistics
├── {GENE}_penetrance_summary.json  # Aggregated penetrance data
├── pmc_fulltext/
│   ├── {PMID}_FULL_CONTEXT.md      # Article + supplements (unified markdown)
│   ├── successful_downloads.csv    # Download log
│   └── paywalled_missing.csv       # Unavailable papers
├── extractions/
│   └── {GENE}_PMID_{pmid}.json     # Per-paper extraction
└── {GENE}.db                       # CANONICAL DATABASE
```

---

## Database Schema

### Core Tables
- **`papers`** - Paper metadata (pmid, title, journal, doi, pmc_id)
- **`variants`** - Normalized variant info (gene_symbol, cdna_notation, protein_notation)
- **`variant_papers`** - Many-to-many linking variants ↔ papers with source quotes

### Data Tables
- **`individual_records`** - Person-level data (patient_id, age, sex, **affected_status**, phenotype_details)
- **`penetrance_data`** - Cohort statistics (affected_count, unaffected_count, uncertain_count)
- **`age_dependent_penetrance`** - Age-stratified penetrance data
- **`functional_data`** - In vitro assay results
- **`phenotypes`** - Detailed phenotype descriptions

---

## Configuration

### Required Environment Variables

```bash
AI_INTEGRATIONS_OPENAI_API_KEY=your-key  # Or OPENAI_API_KEY
NCBI_EMAIL=your@email.com
```

### Optional Environment Variables

```bash
NCBI_API_KEY=your-key                    # Increases rate limits

# Model configuration
TIER3_MODELS=gpt-4o-mini,gpt-4o          # Cascades if too few variants found
TIER3_THRESHOLD=1                         # Minimum variants before cascade
```

### Command-Line Options

```bash
python automated_workflow.py GENE --email EMAIL [OPTIONS]

Options:
  --max-pmids N      Limit number of PMIDs to discover (default: no limit)
  --max-downloads N  Limit number of papers to download (default: no limit)
```

---

## Key Features

- **PubMind-First Strategy** - Variant-specific literature reduces noise vs. broad PubMed searches
- **Supplement Extraction** - Captures data from Excel/Word files that most tools miss
- **Model Cascade** - Starts with cheap LLM, escalates to expensive model if needed
- **Parallel Processing** - ThreadPoolExecutor for 3-5x speedup on extractions
- **Penetrance Validation** - Cross-paper aggregation and consistency checks
- **Normalized Database** - SQLite with proper schema, foreign keys, and indexes

---

## Troubleshooting

**"No PMCID found"**
- Normal: Only ~30% of PubMed articles are in PMC
- Paywalled papers are logged to `paywalled_missing.csv`

**OpenAI API errors**
- Check API key: `echo $AI_INTEGRATIONS_OPENAI_API_KEY`
- Verify you have API credits available

**Slow processing**
- Normal: Rate-limited to respect API servers
- Parallel extraction enabled by default (up to 8 workers)

---

## Development

See [CLAUDE.md](CLAUDE.md) for detailed architecture and module documentation.

**Code Style:**
- Python 3.11+, PEP 8, type hints, Google-style docstrings
- Logging via `logger = logging.getLogger(__name__)`
- Imperative commit messages: "Add feature" not "Added feature"

---

## License

See LICENSE file in repository.
