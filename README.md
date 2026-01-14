# GeneVariantFetcher

Extract human genetic variant carriers from biomedical literature into a SQLite database.

## What It Does

Given a gene name, this tool runs a 9-step pipeline:

1. **Discover Synonyms** (optional) – Auto-find gene aliases via NCBI
2. **Fetch PMIDs** – Query PubMind, PubMed, and Europe PMC
3. **Fetch Abstracts** – Get metadata for all discovered papers
4. **Filter Papers** – Tier 1 keyword filter + Tier 2 LLM relevance triage
5. **Download Full-Text** – Fetch PMC articles and supplements (Excel, Word, PDFs)
6. **Scout Data** – Identify high-value data zones to reduce token usage
7. **Extract Variants** – LLM extraction with model cascade
8. **Aggregate** – Validate and merge penetrance data
9. **Migrate to SQLite** – Create normalized database

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

The GUI is recommended for most use cases. It provides real-time progress, job management, and automatic resume after interruption.

### Filter Tiers

| Tier | Component | Purpose | Model |
|------|-----------|---------|-------|
| **1** | KeywordFilter | Fast heuristic filter | None (regex) |
| **2** | InternFilter | LLM relevance classification | gpt-4o-mini |
| **3** | ExpertExtractor | Full variant extraction | gpt-4o-mini → gpt-4o |

**Model Cascade:** If Tier 3's first model finds fewer variants than threshold, it retries with the next model.

### Genetic Data Scout

Before extraction, the pipeline identifies relevant "data zones" in each paper and creates condensed `DATA_ZONES.md` files. This reduces token usage by 60-80%.

### Abstract-Only Fallback

Papers that cannot be downloaded (paywalled, missing from PMC) are still processed using their abstracts. The `ExpertExtractor` handles both fulltext and abstract-only extraction transparently.

## Environment Variables

Create a `.env` file in the project root:

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| **API Keys** |
| `OPENAI_API_KEY` | Yes | - | OpenAI API key (primary LLM provider) |
| `NCBI_EMAIL` | Yes | - | Email for NCBI E-utilities |
| `NCBI_API_KEY` | No | - | NCBI API key (10 req/s vs 3 req/s) |
| **Paper Sources** |
| `USE_PUBMIND` | No | `true` | Use PubMind as primary source |
| `USE_PUBMED` | No | `true` | Use PubMed keyword search |
| `USE_EUROPEPMC` | No | `false` | Use Europe PMC |
| `PUBMIND_ONLY` | No | `false` | Disable all sources except PubMind |
| `MAX_PAPERS_PER_SOURCE` | No | `100` | Max papers per source |
| **Tier Configuration** |
| `ENABLE_TIER1` | No | `true` | Enable keyword filtering |
| `ENABLE_TIER2` | No | `true` | Enable LLM triage |
| `TIER1_MIN_KEYWORDS` | No | `2` | Min keyword matches for Tier 1 |
| `TIER2_MODEL` | No | `gpt-4o-mini` | Model for Tier 2 |
| `TIER2_CONFIDENCE_THRESHOLD` | No | `0.5` | Min confidence for Tier 2 |
| `TIER3_MODELS` | No | `gpt-4o-mini,gpt-4o` | Model cascade for Tier 3 |
| `TIER3_THRESHOLD` | No | `1` | Retry threshold for cascade |
| **Data Scout** |
| `SCOUT_ENABLED` | No | `true` | Generate DATA_ZONES.md files |
| `SCOUT_MIN_RELEVANCE` | No | `0.3` | Min zone relevance (0.0-1.0) |
| `SCOUT_MAX_ZONES` | No | `30` | Max zones per paper |

## CLI Options

```bash
python main.py --cli GENE --email EMAIL --output DIR [OPTIONS]
# or directly:
python automated_workflow.py GENE --email EMAIL --output DIR [OPTIONS]
```

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
├── {GENE}_penetrance_summary.json   # Aggregated statistics
├── abstract_json/                   # Paper metadata
├── pmid_status/                     # Filter decisions
│   ├── filtered_out.csv
│   └── extraction_failures.csv
├── pmc_fulltext/                    # Full-text + supplements
│   ├── PMID_*_FULL_CONTEXT.md
│   ├── PMID_*_DATA_ZONES.md
│   └── PMID_*_supplements/
└── extractions/                     # Per-paper JSON
```

## Database Schema

| Table | Description |
|-------|-------------|
| `papers` | Paper metadata (pmid, title, doi) |
| `variants` | Variant info (gene, cDNA, protein notation) |
| `individual_records` | Patient data (age, sex, affected_status) |
| `penetrance_data` | Cohort statistics |

See [docs/SQLITE_MIGRATION_GUIDE.md](docs/SQLITE_MIGRATION_GUIDE.md) for details.

## Additional Tools

### Export to CSV

```bash
python tests/extract_ttr_to_csv.py --db path/to/GENE.db --output variants.csv
```

### Compare Curated vs Automated

```bash
python tests/compare_variants.py --excel curated.xlsx --sqlite GENE.db
```

### Semi-Manual Fetch for Paywalled Papers

```bash
python fetch_manager.py results/GENE/TIMESTAMP/pmc_fulltext/paywalled_missing.csv
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

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for design details.

```bash
# Run tests
pytest tests/

# Example scripts
python tests/example_harvest_from_pubmind.py
python tests/example_triage.py
```

## License

See LICENSE file for details.
