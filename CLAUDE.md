# CLAUDE.md

## Mission

Extract human genetic variant carriers from biomedical literature (full text + supplements) into a normalized SQLite database. Each individual is classified as **Affected**, **Unaffected**, or **Ambiguous**.

## Quick Reference

```bash
# GUI (recommended - supports checkpointing/resume)
python main.py

# CLI mode
python main.py --cli GENE --email EMAIL --output OUTPUT_DIR
```

## Entry Points

| Entry | Mode | Features |
|-------|------|----------|
| `python main.py` | GUI (default) | Web interface, checkpointing, resume |
| `python main.py --cli GENE` | CLI | Full pipeline, no checkpointing |
| `python automated_workflow.py GENE` | CLI (direct) | Same as --cli mode |

## Pipeline (9 Steps)

The GUI/CLI run identical logic. GUI adds checkpointing for resume after interruption.

1. **Discover Synonyms** (optional) - Auto-find gene aliases via NCBI (`gene_literature/synonym_finder.py`)
2. **Fetch PMIDs** - Query PubMind, PubMed, Europe PMC (`gene_literature/discovery.py`)
3. **Fetch Abstracts** - Get metadata for all PMIDs (`harvesting/abstracts.py`)
4. **Filter Papers** - Tier 1 keyword + Tier 2 LLM triage (`pipeline/filters.py`)
5. **Download Full-Text** - PMC articles + supplements → `{PMID}_FULL_CONTEXT.md` (`harvesting/orchestrator.py`)
6. **Scout Data** - Identify high-value zones → `{PMID}_DATA_ZONES.md` (`pipeline/data_scout.py`)
7. **Extract Variants** - LLM extraction with model cascade (`pipeline/extraction.py`)
8. **Aggregate** - Penetrance validation (`pipeline/aggregation.py`)
9. **Migrate to SQLite** - JSON → normalized database (`harvesting/migrate_to_sqlite.py`)

### Filter Tiers

| Tier | Component | Purpose | Model |
|------|-----------|---------|-------|
| **1** | `KeywordFilter` | Fast regex filter | None |
| **2** | `InternFilter` / `ClinicalDataTriageFilter` | LLM relevance check | gpt-4o-mini |
| **3** | `ExpertExtractor` | Full extraction with cascade | gpt-4o-mini → gpt-4o |

## Project Structure

```
main.py                        # Primary entry point (GUI/CLI)
automated_workflow.py          # CLI workflow implementation
config/settings.py             # Pydantic configuration
gene_literature/               # Paper discovery (PubMind, PubMed, Europe PMC)
harvesting/                    # Download, abstracts, SQLite migration
pipeline/                      # Filters, extraction, aggregation, data scout
gui/                           # FastAPI web interface + checkpointing
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

## Known Unused/Stale Files

These files exist but are not used in the main pipeline:

| File | Status | Notes |
|------|--------|-------|
| `pipeline/abstract_extraction.py` | **Unused** | Abstract-only extraction is now integrated into `ExpertExtractor.extract()` |
| `pipeline/cli.py` | **Stub** | Only implements `version` command; `gvf` console script is non-functional |
| `gene_literature/collect_literature.py` | **Standalone** | Separate CLI tool, not part of main pipeline |
| `run_gui.py` | **Redundant** | Duplicates `main.py` GUI startup logic |

## Key Implementation Details

- **Abstract-only fallback**: Papers that fail download are still extracted using their abstracts (handled in `ExpertExtractor`)
- **DATA_ZONES preference**: Extraction prefers `DATA_ZONES.md` over `FULL_CONTEXT.md` to reduce token usage
- **Checkpoint system**: GUI jobs save state after each step to `~/.gvf_jobs/{job_id}/checkpoint.json`
- **Folder jobs**: GUI supports processing existing folders (skip discovery/download, start at extraction)
