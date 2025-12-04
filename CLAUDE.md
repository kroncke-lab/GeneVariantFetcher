# CLAUDE.md

**GeneVariantFetcher** extracts human genetic variant carriers from biomedical literature (including supplements), classifies them as Affected/Unaffected/Ambiguous, and stores results in SQLite.

## Critical Rules

1. **Entry point:** `python automated_workflow.py GENE --email EMAIL`
2. **SQLite is canonical** - all data must be written to `{GENE}.db`
3. **Never use `example_*.py` in production** - they're demos only
4. **Supplements matter** - 70-80% of variant data is in Excel/Word supplements

## Pipeline Stages

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

## Key Modules

- `automated_workflow.py` - Main production entry point, orchestrates everything
- `harvesting/orchestrator.py` - PMCHarvester class for downloading full-text
- `clinical_data_triage.py` - Filters papers before expensive extraction
- `migrate_to_sqlite.py` - JSON → SQLite migration
- `query_variants_db.py` - Query the database

## Output Structure

```
automated_output/{GENE}/{TIMESTAMP}/
├── {GENE}_pmids.txt           # Discovered PMIDs
├── pmc_fulltext/
│   └── {PMID}_FULL_CONTEXT.md # Article + supplements
├── extractions/
│   └── {PMID}_extraction.json # Extracted data
└── {GENE}.db                  # CANONICAL DATABASE
```

## Database Tables

- `papers` - Paper metadata (PMID, title, DOI)
- `individuals` - Case-level data (patient_id, variant, affected_status, phenotype)
- `aggregates` - Cohort statistics (carrier counts, penetrance)
- `variants` - Normalized variant info

## CLI Quick Reference

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
AI_INTEGRATIONS_OPENAI_API_KEY=your-key  # Required
NCBI_EMAIL=your@email.com                 # Required
NCBI_API_KEY=your-key                     # Optional, increases rate limits
```

## Code Style

- Python 3.11+, PEP 8, type hints, Google-style docstrings
- Logging via `logger = logging.getLogger(__name__)`
- Imperative commit messages: "Add feature" not "Added feature"
