# CLAUDE.md

## Mission

**Goal:** Build a unified workflow that extracts all human carriers of genetic variants from biomedical literature (including full text and supplements).

**Core Criteria:** Only humans with clinical phenotyping that allows classification as **Affected**, **Unaffected**, or **Ambiguous** based on known gene-disease associations.

**Output:** A single harmonized SQLite database of variant carriers.

---

## Critical Rules

1. **Entry point:** `python automated_workflow.py GENE --email EMAIL`
2. **SQLite is canonical** - all data must be written to `{GENE}.db`
3. **Supplements are essential** - 70-80% of variant data is in Excel/Word supplements
4. **Classification matters** - every individual must be classified as Affected/Unaffected/Ambiguous

---

## Pipeline Architecture

```
Discovery → Download → Extract → Aggregate → SQLite
```

| Stage | Module | Description |
|-------|--------|-------------|
| Discovery | `gene_literature/pubmind_fetcher.py` | Fetch PMIDs from PubMind (variant-specific literature) |
| Download | `harvesting/orchestrator.py` | Full-text + supplements → `{PMID}_FULL_CONTEXT.md` |
| Extract | `pipeline/extraction.py` | LLM extraction with model cascade (parallel processing) |
| Aggregate | `pipeline/aggregation.py` | Penetrance validation and cross-paper statistics |
| Storage | `harvesting/migrate_to_sqlite.py` | JSON → normalized SQLite database |

---

## Key Modules

### Core Pipeline
- `automated_workflow.py` - Main entry point, orchestrates all 6 stages
- `harvesting/orchestrator.py` - PMCHarvester class for full-text + supplement downloads
- `pipeline/extraction.py` - ExpertExtractor for structured variant extraction
- `pipeline/aggregation.py` - DataAggregator for penetrance validation
- `harvesting/migrate_to_sqlite.py` - JSON → SQLite migration

### Supporting Modules
- `gene_literature/pubmind_fetcher.py` - Paper discovery via PubMind
- `gene_literature/clinical_data_triage.py` - Standalone paper filtering tool
- `pipeline/filters.py` - Three-tier filtering (Keyword → LLM → Expert)
- `config/settings.py` - Pydantic-based configuration management

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
- **`penetrance_data`** - Cohort statistics (affected_count, unaffected_count, uncertain_count, penetrance_percentage)
- **`age_dependent_penetrance`** - Age-stratified penetrance data
- **`functional_data`** - In vitro assay results
- **`phenotypes`** - Detailed phenotype descriptions

---

## CLI Reference

```bash
# Run full pipeline
python automated_workflow.py BRCA1 --email you@email.com

# With limits
python automated_workflow.py SCN5A --email you@email.com --max-pmids 200 --max-downloads 100
```

---

## Environment Variables

```bash
# Required
AI_INTEGRATIONS_OPENAI_API_KEY=your-key  # Or OPENAI_API_KEY
NCBI_EMAIL=your@email.com

# Optional
NCBI_API_KEY=your-key                    # Increases rate limits

# Model configuration (via config/settings.py)
TIER3_MODELS=gpt-4o-mini,gpt-4o          # Cascades if too few variants found
TIER3_THRESHOLD=1                         # Minimum variants before cascade
```

---

## Code Style

- Python 3.11+, PEP 8, type hints, Google-style docstrings
- Logging via `logger = logging.getLogger(__name__)`
- Imperative commit messages: "Add feature" not "Added feature"
