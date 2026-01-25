# GeneVariantFetcher

Extract human genetic variant carriers from biomedical literature into a SQLite database.

## What It Do

Given a gene name, this tool runs a pipeline with these major steps:

1. **Discover Synonyms** (optional) – Auto-find gene aliases via NCBI
2. **Fetch PMIDs** – Query PubMind, PubMed, and Europe PMC
3. **Fetch Abstracts** – Get metadata for all discovered papers
4. **Filter Papers** – Tier 1 keyword filter + Tier 2 LLM relevance triage
5. **Download Full-Text** – Fetch PMC articles and supplements (Excel, Word, PDFs)
   - Automatically creates condensed `DATA_ZONES.md` files during download to reduce token usage
6. **Extract Variants** – LLM extraction with model cascade (uses `DATA_ZONES.md` if available, falls back to full text or abstracts)
7. **Aggregate** – Validate and merge penetrance data
8. **Migrate to SQLite** – Create normalized database

**Note:** The Data Scout runs automatically during the Download Full-Text step (step 5), creating condensed `DATA_ZONES.md` files. Papers that can't be downloaded are still processed using abstract-only extraction.

**Key Insight:** 70-80% of variant data is in supplemental files. This tool extracts it all.

## Quick Start

```bash
# Install
pip install -e .
pip install -r gui/requirements.txt

# Set required environment variables (or use .env file)
export OPENAI_API_KEY="your-key"
export NCBI_EMAIL="your@email.com"

# Launch the GUI (recommended)
python main.py
```

On first launch, go to the **Settings** tab to configure your API keys.

### CLI Mode

```bash
python main.py --cli BRCA1 --email your@email.com --output ./results
```

For advanced CLI options (limits, triage, synonyms), use `automated_workflow.py` directly.

## Prerequisites

- Python 3.11+
- System dependencies for PDF processing:

```bash
# Ubuntu/Debian
sudo apt-get install poppler-utils

# macOS
brew install poppler
```

## Pipeline Architecture

### Entry Points

| Command | Mode | Features |
|---------|------|----------|
| `python main.py` | **GUI** (default) | Web interface, checkpointing, pause/resume |
| `python main.py --cli GENE` | CLI | Full pipeline, no checkpointing |
| `python automated_workflow.py GENE` | CLI (direct) | Advanced options (`--auto-synonyms`, `--max-downloads`, etc.) |

The GUI is recommended for most use cases. It provides real-time progress, job management, and automatic resume after interruption.

### Filter Tiers

| Tier | Component | Purpose | Model |
|------|-----------|---------|-------|
| **1** | KeywordFilter | Fast heuristic filter | None (regex) |
| **2** | InternFilter | LLM relevance classification | gpt-4o-mini |
| **3** | ExpertExtractor | Full variant extraction | gpt-4o-mini → gpt-4o |

**Model Cascade:** If Tier 3's first model finds fewer variants than threshold, it retries with the next model.

### Genetic Data Scout

During the download step, the pipeline automatically identifies relevant "data zones" in each paper and creates condensed `DATA_ZONES.md` files. The extraction step prefers these condensed files over full text, reducing token usage by 60-80%. This happens automatically - no separate step needed.

### Abstract-Only Fallback

Papers that cannot be downloaded (paywalled, missing from PMC, or failed download) are still processed using their abstracts. The `ExpertExtractor` handles both fulltext and abstract-only extraction transparently. These papers are tracked in `paywalled_missing.csv` and can be manually downloaded using `fetch_manager.py` if needed.

## Environment Variables

Create a `.env` file in the project root:

```bash
OPENAI_API_KEY=...       # Required - for LLM extraction
NCBI_EMAIL=...           # Required - for PubMed API access
NCBI_API_KEY=...         # Optional (higher rate limits)
ANTHROPIC_API_KEY=...    # Optional (for browser_fetch --use-claude)
GEMINI_API_KEY=...       # Optional (alternative LLM provider)
ELSEVIER_API_KEY=...     # Optional (publisher API access)
WILEY_API_KEY=...        # Optional (publisher API access)
```

See `config/settings.py` for all configuration options.

## CLI Options

```bash
python main.py --cli GENE --email EMAIL --output DIR
# for advanced options:
python automated_workflow.py GENE --email EMAIL --output DIR [OPTIONS]
```

Options below apply to `automated_workflow.py`.

| Option | Description |
|--------|-------------|
| `--max-pmids N` | Limit PMIDs to discover (default: 100) |
| `--max-downloads N` | Limit papers to download (default: 50) |
| `--tier-threshold N` | Model cascade threshold (default: 1) |
| `--clinical-triage` | Use ClinicalDataTriageFilter for Tier 2 |
| `--auto-synonyms` | Auto-discover gene synonyms from NCBI |
| `--synonym SYN` | Manual synonym (repeatable) |
| `--verbose` | Enable verbose logging |

### Examples

```bash
# Basic run
python main.py --cli TTR --email user@email.com --output ./results

# With synonym discovery
python main.py --cli TTR --email user@email.com --output ./results --auto-synonyms

# Skip to stronger model
python main.py --cli SCN5A --email user@email.com --output ./results --tier-threshold 0
```

## Output Structure

```
{OUTPUT_DIR}/{GENE}/{TIMESTAMP}/
├── {GENE}.db                        # SQLite database (final output)
├── {GENE}_pmids.txt                 # Discovered PMIDs
├── {GENE}_workflow_summary.json     # Overall workflow statistics
├── {GENE}_penetrance_summary.json   # Aggregated penetrance statistics
├── {GENE}_workflow.log              # Workflow execution log
├── abstract_json/                   # Paper metadata (JSON files per PMID)
├── pmid_status/                     # Filter decisions and failures
│   ├── filtered_out.csv
│   └── extraction_failures.csv
├── pmc_fulltext/                    # Full-text + supplements
│   ├── PMID_*_FULL_CONTEXT.md       # Full paper content (main + supplements)
│   ├── PMID_*_DATA_ZONES.md         # Condensed high-value sections (created during download)
│   ├── PMID_*_supplements/          # Original supplement files
│   ├── successful_downloads.csv     # Successfully downloaded papers
│   └── paywalled_missing.csv        # Papers that couldn't be downloaded
└── extractions/                     # Per-paper JSON extractions
    └── {GENE}_PMID_*.json           # Structured variant data per paper
```

## Database Schema

| Table | Key Columns |
|-------|-------------|
| `papers` | pmid (PK), title, journal, doi, pmc_id, gene_symbol, extraction_timestamp |
| `variants` | variant_id (PK), gene_symbol, cdna_notation, protein_notation, genomic_position, clinical_significance |
| `individual_records` | record_id (PK), variant_id (FK), pmid (FK), age_at_onset, sex, affected_status (affected/unaffected/uncertain), phenotype_details |
| `penetrance_data` | penetrance_id (PK), variant_id (FK), pmid (FK), total_carriers_observed, affected_count, unaffected_count, penetrance_percentage |
| `age_dependent_penetrance` | age_penetrance_id (PK), penetrance_id (FK), age_range, penetrance_percentage, carriers_in_range |
| `functional_data` | functional_id (PK), variant_id (FK), pmid (FK), summary, assays (JSON) |
| `phenotypes` | phenotype_id (PK), variant_id (FK), pmid (FK), patient_count, phenotype_description |
| `variant_papers` | variant_id + pmid (composite PK), source_location, key_quotes (JSON) |

See [docs/SQLITE_MIGRATION_GUIDE.md](docs/SQLITE_MIGRATION_GUIDE.md) for full schema details.

## Additional Tools

### Export to CSV

```bash
python scripts/extract_ttr_to_csv.py --db path/to/GENE.db --output variants.csv
```

### Compare Curated vs Automated

```bash
python tests/compare_variants.py --excel curated.xlsx --sqlite GENE.db
```

### Fetching Paywalled Papers

Two options for downloading papers that couldn't be auto-fetched:

**Option 1: Browser Automation (recommended)**
```bash
# Install browser automation (one time)
pip install playwright && playwright install chromium

# Interactive mode: browser opens, you log in/download, files auto-organized
python browser_fetch.py results/GENE/TIMESTAMP/pmc_fulltext/paywalled_missing.csv --interactive

# With CAPTCHA support: waits up to 5 min for manual CAPTCHA completion
python browser_fetch.py paywalled_missing.csv --interactive --wait-for-captcha

# Retry only previously failed papers
python browser_fetch.py paywalled_missing.csv --retry-failures

# Fully automated (works for open-access papers)
python browser_fetch.py paywalled_missing.csv --use-claude
```

**Option 2: Manual with fetch_manager**
```bash
python fetch_manager.py results/GENE/TIMESTAMP/pmc_fulltext/paywalled_missing.csv
```

**Post-download processing:**
```bash
# Convert PDFs to markdown and run scout
python fetch_manager.py paywalled_missing.csv --convert --run-scout --gene GENE
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "No papers downloaded" | Check `paywalled_missing.csv`, use `fetch_manager.py` |
| "Extraction failed" | Check `extraction_failures.csv`, increase `TIER3_MAX_TOKENS` |
| "All papers filtered" | Lower `TIER2_CONFIDENCE_THRESHOLD` or disable filters |
| Rate limiting | Set `NCBI_API_KEY` for higher limits |
| GUI won't start | Install: `pip install -r gui/requirements.txt` |

## GUI Features

### Checkpointing & Resume

Jobs automatically save state after each pipeline step to `~/.gvf_jobs/`. If interrupted, jobs can be resumed from the last completed step via the GUI's job management interface.

### Folder Jobs

The GUI supports "folder jobs" for processing existing paper collections:
- Skip discovery and download steps
- Start directly at scouting or extraction
- Useful for re-processing or adding new extractions

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidelines and [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for design details.

```bash
# Run tests
pytest tests/

# Example scripts
python docs/examples/example_harvest_from_pubmind.py
python docs/examples/example_triage.py
```

## License

Wild Wild West
