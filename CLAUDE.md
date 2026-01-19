# CLAUDE.md

## CRITICAL: Do Not Modify

**NEVER read, write, edit, or overwrite the `.env` file.** It contains API keys and credentials that are difficult to replace. If you need to reference environment variables, only document their names - never touch the actual file.

---

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
| `python automated_workflow.py GENE --email EMAIL --output OUTPUT_DIR` | CLI (direct) | Full CLI with advanced options |

## Pipeline Steps

The GUI/CLI run identical logic. GUI adds checkpointing for resume after interruption.

1. **Discover Synonyms** (optional) - Auto-find gene aliases via NCBI (`gene_literature/synonym_finder.py`)
2. **Fetch PMIDs** - Query PubMind, PubMed, Europe PMC (`gene_literature/discovery.py`)
3. **Fetch Abstracts** - Get metadata for all PMIDs (`harvesting/abstracts.py`)
4. **Filter Papers** - Tier 1 keyword + Tier 2 LLM triage (`pipeline/filters.py`)
5. **Download Full-Text** - PMC articles + supplements → `{PMID}_FULL_CONTEXT.md` (`harvesting/orchestrator.py`)
   - **Automatically runs Data Scout** during download → creates `{PMID}_DATA_ZONES.md` (`pipeline/data_scout.py`)
6. **Extract Variants** - LLM extraction with model cascade (`pipeline/extraction.py`)
   - Prefers `DATA_ZONES.md` if available, falls back to `FULL_CONTEXT.md` or abstract
   - Handles abstract-only extraction for papers that failed download
7. **Aggregate** - Penetrance validation (`pipeline/aggregation.py`)
8. **Migrate to SQLite** - JSON → normalized database (`harvesting/migrate_to_sqlite.py`)

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
├── {GENE}_pmids.txt                    # Discovered PMIDs
├── {GENE}_workflow_summary.json        # Overall workflow statistics
├── {GENE}_penetrance_summary.json      # Aggregated penetrance statistics
├── {GENE}_workflow.log                 # Workflow execution log
├── abstract_json/                      # Paper metadata (JSON per PMID)
├── pmid_status/                        # Filter decisions and failures
│   ├── filtered_out.csv
│   └── extraction_failures.csv
├── pmc_fulltext/                       # Full-text markdown + logs
│   ├── PMID_*_FULL_CONTEXT.md          # Full paper content
│   ├── PMID_*_DATA_ZONES.md            # Condensed sections (auto-created)
│   ├── successful_downloads.csv
│   └── paywalled_missing.csv           # Use with fetch_manager.py
├── extractions/                        # Per-paper JSON extractions
│   └── {GENE}_PMID_*.json
└── {GENE}.db                           # Canonical SQLite database
```

## Environment Variables

All keys are loaded from `.env` file:

```bash
OPENAI_API_KEY=...                   # Required
NCBI_EMAIL=...                       # Required
NCBI_API_KEY=...                     # Optional (higher rate limits)
ANTHROPIC_API_KEY=...                # Optional
GEMINI_API_KEY=...                   # Optional
ELSEVIER_API_KEY=...                 # Optional (publisher access)
WILEY_API_KEY=...                    # Optional (publisher access)
```

## Code Style

- Python 3.11+, PEP 8, type hints
- Logging: `logger = logging.getLogger(__name__)`
- Commits: imperative mood ("Add feature" not "Added feature")

## Standalone Tools

These files are standalone utilities, not part of the automated pipeline:

| File | Purpose | Usage |
|------|---------|-------|
| `fetch_manager.py` | Semi-manual download helper for paywalled papers | `python fetch_manager.py paywalled_missing.csv --convert --run-scout --gene GENE` |
| `gene_literature/collect_literature.py` | Standalone literature collection CLI | Direct invocation for one-off literature searches |

## Known Unused/Stale Files

These files exist but are not used in the main pipeline:

| File | Status | Notes |
|------|--------|-------|
| `pipeline/abstract_extraction.py` | **Unused** | Abstract-only extraction is now integrated into `ExpertExtractor.extract()` |
| `pipeline/cli.py` | **Stub** | Only implements `version` command; `gvf` console script is non-functional |
| `run_gui.py` | **Redundant** | Duplicates `main.py` GUI startup logic |

## Key Implementation Details

- **Abstract-only fallback**: Papers that fail download are still extracted using their abstracts (handled in `ExpertExtractor`)
- **DATA_ZONES preference**: Extraction prefers `DATA_ZONES.md` over `FULL_CONTEXT.md` to reduce token usage (created automatically during download)
- **Scout Data integration**: Data Scout runs automatically during `PMCHarvester.harvest()` - no separate step needed
- **Checkpoint system**: GUI jobs save state after each step to `~/.gvf_jobs/{job_id}/checkpoint.json`
- **Folder jobs**: GUI supports processing existing folders (skip discovery/download, start at extraction)
- **Two summary files**: Both `{GENE}_workflow_summary.json` (overall stats) and `{GENE}_penetrance_summary.json` (variant aggregations) are created
- **Functional study detection**: Extraction prompts distinguish functional studies (massively parallel assays) from clinical studies to avoid misinterpreting assay replicates as patient carrier counts