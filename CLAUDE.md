# CLAUDE.md

## Mission

Extract human genetic variant carriers from biomedical literature (full text + supplements) into a normalized SQLite database. Each individual is classified as **Affected**, **Unaffected**, or **Ambiguous**.

## Quick Reference

```bash
python automated_workflow.py GENE --email EMAIL --output OUTPUT_DIR
```

## Pipeline (6 Steps)

1. **Discovery** - Fetch PMIDs from PubMind (`gene_literature/pubmind_fetcher.py`)
2. **Abstracts** - Fetch metadata for discovered PMIDs (`harvesting/abstracts.py`)
3. **Download** - Full-text + supplements → `{PMID}_FULL_CONTEXT.md` (`harvesting/orchestrator.py`)
4. **Extract** - LLM extraction with model cascade (`pipeline/extraction.py`)
5. **Aggregate** - Penetrance validation (`pipeline/aggregation.py`)
6. **SQLite** - JSON → normalized database (`harvesting/migrate_to_sqlite.py`)

## Project Structure

```
automated_workflow.py          # Entry point
config/settings.py             # Pydantic configuration
gene_literature/               # Paper discovery (PubMind)
harvesting/                    # Download & SQLite migration
pipeline/                      # Extraction & aggregation
utils/                         # Shared utilities
```

## Output

```
{OUTPUT_DIR}/{GENE}/{TIMESTAMP}/
├── {GENE}_pmids.txt           # Discovered PMIDs
├── abstract_json/             # Paper metadata
├── pmc_fulltext/              # Full-text markdown + logs
├── extractions/               # Per-paper JSON extractions
├── {GENE}_penetrance_summary.json
└── {GENE}.db                  # Canonical SQLite database
```

## Environment Variables

```bash
AI_INTEGRATIONS_OPENAI_API_KEY=...  # Required (or OPENAI_API_KEY)
NCBI_EMAIL=...                       # Required
NCBI_API_KEY=...                     # Optional (higher rate limits)
```

## Code Style

- Python 3.11+, PEP 8, type hints
- Logging: `logger = logging.getLogger(__name__)`
- Commits: imperative mood ("Add feature" not "Added feature")
