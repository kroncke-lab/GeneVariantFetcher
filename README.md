# GeneVariantFetcher

Extracts human genetic variant carriers from biomedical literature (including supplements), classifies them as Affected/Unaffected/Ambiguous, and stores results in SQLite.

**Key insight**: 70-80% of variant data is in Excel/Word supplements, not article text. This tool downloads and processes both.

## Quick Start

```bash
# Install
pip install -e .

# Set environment variables
export AI_INTEGRATIONS_OPENAI_API_KEY="your-key"
export NCBI_EMAIL="your@email.com"

# Run complete pipeline (one command)
python automated_workflow.py BRCA1 --email your@email.com
```

That's it. You now have a SQLite database with all extracted variant and patient data.

## Pipeline

```
Discovery → Download → Filter → Extract → SQLite
```

| Stage | Module | Description |
|-------|--------|-------------|
| Discovery | `gene_literature/pubmind_fetcher.py` | Fetch PMIDs from PubMind/PubMed |
| Download | `harvesting/orchestrator.py` | Full-text + supplements → `{PMID}_FULL_CONTEXT.md` |
| Filter | `clinical_data_triage.py` | Tier 1 (keywords) + Tier 2 (cheap LLM) |
| Extract | `automated_workflow.py` | LLM extraction with model cascade |
| Storage | `migrate_to_sqlite.py` | Write to SQLite database |

## Output Structure

```
automated_output/{GENE}/{TIMESTAMP}/
├── {GENE}_pmids.txt           # Discovered PMIDs
├── pmc_fulltext/
│   └── {PMID}_FULL_CONTEXT.md # Article + supplements as markdown
├── extractions/
│   └── {PMID}_extraction.json # Extracted variant/patient data
└── {GENE}.db                  # CANONICAL DATABASE
```

## Database Schema

| Table | Contents |
|-------|----------|
| `papers` | Paper metadata (PMID, title, DOI) |
| `individuals` | Case-level data (patient_id, variant, affected_status, phenotype) |
| `aggregates` | Cohort statistics (carrier counts, penetrance) |
| `variants` | Normalized variant info |

## CLI Reference

```bash
# Run full pipeline
python automated_workflow.py BRCA1 --email you@email.com --max-pmids 100

# Query database
python query_variants_db.py GENE.db --stats
python query_variants_db.py GENE.db --export output.csv

# Re-run extraction on existing downloads
python rerun_extraction.py automated_output/GENE/{timestamp}/
```

## Environment Variables

```bash
AI_INTEGRATIONS_OPENAI_API_KEY=your-key  # Required - OpenAI API key
NCBI_EMAIL=your@email.com                 # Required - for PubMed API
NCBI_API_KEY=your-key                     # Optional - increases rate limits
```

## Requirements

- Python 3.11+
- OpenAI API key (for LLM extraction)
- NCBI email (for PubMed API compliance)

## Key Modules

| Module | Purpose |
|--------|---------|
| `automated_workflow.py` | Main entry point, orchestrates everything |
| `harvesting/orchestrator.py` | PMCHarvester for downloading full-text |
| `clinical_data_triage.py` | Filters papers before expensive extraction |
| `migrate_to_sqlite.py` | JSON → SQLite migration |
| `query_variants_db.py` | Query the database |

## Notes

- Only ~30% of PubMed articles are available full-text in PMC. Paywalled articles are logged and skipped.
- The `example_*.py` files are demos only—use `automated_workflow.py` for production.
- SQLite database is the canonical output format.
