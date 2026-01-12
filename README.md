# GeneVariantFetcher

Extract human genetic variant carriers from biomedical literature into a SQLite database.

## What It Does

Given a gene name, this tool:
1. **Discovers** variant-related papers via PubMind, PubMed, and Europe PMC
2. **Filters** papers through a 3-tier relevance pipeline (keywords → LLM triage → extraction)
3. **Downloads** full-text articles and supplemental materials (Excel, Word, PDFs)
4. **Scouts** data zones in papers to focus extraction on relevant sections
5. **Extracts** patient-level variant data using an LLM model cascade
6. **Aggregates** penetrance statistics across papers
7. **Writes** normalized data to SQLite

**Key Insight:** 70-80% of variant data is in supplemental files. This tool extracts it all.

## Prerequisites

- Python 3.11+
- System dependencies for PDF/document processing:

```bash
# Ubuntu/Debian
sudo apt-get install poppler-utils

# macOS
brew install poppler
```

## Quick Start

```bash
# Install
pip install -e .

# Set environment variables
export OPENAI_API_KEY="your-key"
export ANTHROPIC_API_KEY="your-key"
export NCBI_EMAIL="your@email.com"

# Run full pipeline
python automated_workflow.py BRCA1 --email your@email.com --output ./results
```

---

## Pipeline Architecture

The pipeline uses a 3-tier filtering system to efficiently process papers:

| Tier | Component | Purpose | Model | Config |
|------|-----------|---------|-------|--------|
| **Tier 1** | KeywordFilter | Fast heuristic filter on abstracts | None (regex) | `ENABLE_TIER1`, `TIER1_MIN_KEYWORDS` |
| **Tier 2** | InternFilter | LLM-based relevance classification | gpt-4o-mini | `ENABLE_TIER2`, `TIER2_CONFIDENCE_THRESHOLD` |
| **Tier 3** | ExpertExtractor | Full variant/patient extraction | Model cascade | `TIER3_MODELS`, `TIER3_THRESHOLD` |

**Model Cascade:** If Tier 3's first model finds fewer variants than `TIER3_THRESHOLD`, it automatically retries with the next model in the cascade (e.g., `gpt-4o-mini` → `gpt-4o`).

Papers dropped at any tier are logged to `pmid_status/filtered_out.csv` with reasons.

### Genetic Data Scout

Before extraction, the pipeline identifies relevant "data zones" in each paper—tables, patient cohort sections, variant lists—and creates condensed `DATA_ZONES.md` files. This reduces token usage by 60-80% and improves extraction accuracy.

| Setting | Default | Description |
|---------|---------|-------------|
| `SCOUT_ENABLED` | `true` | Enable/disable zone detection |
| `SCOUT_MIN_RELEVANCE` | `0.3` | Minimum relevance score (0.0-1.0) for zones |
| `SCOUT_MAX_ZONES` | `30` | Maximum zones per paper |

---

## Use Cases

### 1. End-to-End Pipeline (Primary Use Case)

Run the complete automated workflow from gene name to SQLite database.

**Script:** `automated_workflow.py`

```bash
python automated_workflow.py GENE --email EMAIL --output DIR [OPTIONS]
```

**Options:**
| Option | Description |
|--------|-------------|
| `--max-pmids N` | Limit PMIDs to discover (default: 100) |
| `--max-downloads N` | Limit papers to download (default: 50) |
| `--tier-threshold N` | Model cascade threshold (default: 1) |
| `--clinical-triage` | Use ClinicalDataTriageFilter for Tier 2 |
| `--auto-synonyms` | Auto-discover gene synonyms from NCBI |
| `--synonym SYN` | Manual synonym (repeatable) |
| `--verbose` | Enable verbose logging |

**Examples:**
```bash
# Basic run
python automated_workflow.py TTR --email user@email.com --output ./results

# With automatic synonym discovery
python automated_workflow.py TTR --email user@email.com --output ./results --auto-synonyms

# With manual synonyms
python automated_workflow.py TTR --email user@email.com --output ./results \
  --synonym PALB --synonym prealbumin

# PubMind-only (skip PubMed keyword search)
PUBMIND_ONLY=true python automated_workflow.py KCNH2 --email user@email.com --output ./results

# Skip to stronger model (disable cascade)
python automated_workflow.py SCN5A --email user@email.com --output ./results --tier-threshold 0
```

---

### 2. Browser-Based GUI

Web interface for running pipelines with background job support.

**Script:** `run_gui.py`

```bash
# Install GUI dependencies
pip install -r gui/requirements.txt

# Start server
python run_gui.py                    # localhost:8000
python run_gui.py --port 8080        # Custom port
python run_gui.py --host 0.0.0.0     # Allow external connections
```

**Features:**
- Background job execution (survives browser close)
- Checkpoint/resume for interrupted jobs
- Real-time progress updates via WebSocket
- Job history and status tracking

---

### 3. Export Database to CSV

Export variant data from SQLite to CSV for analysis in Excel, R, or Python.

**Script:** `extract_ttr_to_csv.py`

```bash
python extract_ttr_to_csv.py --db path/to/GENE.db [OPTIONS]
```

**Options:**
| Option | Description |
|--------|-------------|
| `--db PATH` | Path to SQLite database (required) |
| `--output FILE` | Output filename (default: `ttr_variants.csv`) |
| `--no-individual-records` | Skip individual records CSV |

**Example:**
```bash
python extract_ttr_to_csv.py --db results/TTR/20250101_120000/TTR.db --output ttr_analysis.csv
```

**Output:**
- `ttr_analysis.csv` - Variant-level aggregations with affected/unaffected counts
- `ttr_analysis_individual_records.csv` - Individual patient records with onset ages

---

### 4. Compare Curated vs Automated Results

Validate automated extractions against manually curated Excel sheets.

**Script:** `compare_variants.py`

```bash
python compare_variants.py --excel CURATED.xlsx --sqlite GENE.db [OPTIONS]
```

**Options:**
| Option | Description |
|--------|-------------|
| `--excel FILE` | Curated Excel file (required) |
| `--sqlite DB` | SQLite database (required) |
| `--variant-match-mode` | `exact` or `fuzzy` (default: fuzzy) |
| `--fuzzy-threshold` | Similarity threshold (default: 0.85) |
| `--output REPORT` | Output report path |

**Example:**
```bash
python compare_variants.py \
  --excel manual_curation.xlsx \
  --sqlite results/KCNH2/20250101_120000/KCNH2.db \
  --variant-match-mode fuzzy \
  --fuzzy-threshold 0.85
```

---

### 5. Semi-Manual Fetch for Paywalled Papers

Browser-assisted workflow to recover papers the automated harvester couldn't access.

**Script:** `fetch_manager.py`

```bash
python fetch_manager.py <csv_file> [OPTIONS]
```

**How it works:**
1. Opens paper URL in your browser (DOI preferred)
2. Waits for you to manually download PDF/supplements
3. Detects new files in your Downloads folder
4. Renames and moves files to target directory
5. Tracks progress in CSV

**Options:**
| Option | Description |
|--------|-------------|
| `--doi-column` | Custom DOI column name |
| `--target-dir` | Custom target directory |
| `--downloads-dir` | Custom Downloads folder path |

**Example:**
```bash
# Use the paywalled_missing.csv generated by the pipeline
python fetch_manager.py results/TTR/20250101_120000/pmc_fulltext/paywalled_missing.csv
```

---

### 6. SQLite Migration (Standalone)

Convert JSON extraction files to SQLite database. Runs automatically as part of the pipeline, but can be used standalone for re-processing.

**Script:** `harvesting/migrate_to_sqlite.py`

```bash
python harvesting/migrate_to_sqlite.py --data-dir /path/to/output/GENE/TIMESTAMP [OPTIONS]
```

**Options:**
| Option | Description |
|--------|-------------|
| `--data-dir PATH` | Directory with extraction files (required) |
| `--db DBNAME` | Custom database name |
| `--cleanup` | Delete empty directories after migration |
| `--delete-pmc-after-archive` | Archive and delete full-text directory |
| `--dry-run` | Preview changes without applying |

**Example:**
```bash
python harvesting/migrate_to_sqlite.py \
  --data-dir results/TTR/20250101_120000 \
  --cleanup
```

---

### 7. Library Usage (Component-Level)

Import individual components for custom workflows.

#### Full-Text Harvesting Only

```python
from harvesting import PMCHarvester

harvester = PMCHarvester(output_dir="my_papers")
harvester.harvest(['35443093', '33442691'], delay=2.0)
```

#### Clinical Data Triage Only

```python
from pipeline.filters import ClinicalDataTriageFilter

triage = ClinicalDataTriageFilter()
result = triage.triage(
    title="Novel SCN5A Mutation...",
    abstract="We report a 28-year-old...",
    gene="SCN5A",
    pmid="12345678"
)
print(result.decision)  # "EXTRACT" or "SKIP"
```

#### PMID Discovery Only

```python
from gene_literature.discovery import discover_pmids_for_gene

result = discover_pmids_for_gene(
    gene_symbol="BRCA1",
    email="user@email.com",
    max_results=100,
    synonyms=["BRCA1", "PSCP"]
)
print(result.pubmind_pmids)
print(result.pubmed_pmids)
```

---

## Output Structure

```
{OUTPUT_DIR}/{GENE}/{TIMESTAMP}/
├── {GENE}.db                        # SQLite database (final output)
├── {GENE}_pmids.txt                 # Merged/deduplicated PMIDs
├── {GENE}_pmids_pubmind.txt         # PubMind-only PMIDs
├── {GENE}_pmids_pubmed.txt          # PubMed-only PMIDs
├── {GENE}_penetrance_summary.json   # Aggregated penetrance statistics
├── {GENE}_workflow_summary.json     # Pipeline run statistics
├── abstract_json/                   # Paper metadata
├── pmid_status/                     # PMID tracking and diagnostics
│   ├── filtered_out.csv             # PMIDs dropped by Tier 1/2 filters
│   └── extraction_failures.csv      # PMIDs that failed Tier 3 extraction
├── pmc_fulltext/                    # Full-text + supplements
│   ├── PMID_35443093_FULL_CONTEXT.md
│   ├── PMID_35443093_DATA_ZONES.md  # Condensed version (if Scout enabled)
│   ├── PMID_35443093_supplements/
│   │   ├── Table_S1.xlsx
│   │   └── Figure_S2.pdf
│   ├── successful_downloads.csv     # Successfully downloaded papers
│   └── paywalled_missing.csv        # Papers needing manual download
└── extractions/                     # Per-paper JSON extractions
    ├── GENE_PMID_35443093.json
    └── GENE_PMID_33442691.json
```

## Database Schema

| Table | Description |
|-------|-------------|
| **papers** | Paper metadata (pmid, title, doi) |
| **variants** | Variant info (gene, cDNA, protein notation) |
| **individual_records** | Patient data (age, sex, affected_status, phenotype) |
| **penetrance_data** | Cohort statistics |

---

## Environment Variables

### Required

| Variable | Description |
|----------|-------------|
| `OPENAI_API_KEY` | OpenAI API key (or `AI_INTEGRATIONS_OPENAI_API_KEY`) |
| `ANTHROPIC_API_KEY` | Anthropic API key |
| `NCBI_EMAIL` | Email for NCBI E-utilities |

### Optional - API Keys

| Variable | Default | Description |
|----------|---------|-------------|
| `NCBI_API_KEY` | - | Increases NCBI rate limits (10 req/s vs 3 req/s) |

### Optional - Paper Sources

| Variable | Default | Description |
|----------|---------|-------------|
| `USE_PUBMIND` | `true` | Use PubMind as primary source |
| `USE_PUBMED` | `true` | Use PubMed keyword search |
| `USE_EUROPEPMC` | `false` | Use Europe PMC as additional source |
| `PUBMIND_ONLY` | `false` | Disable all sources except PubMind |
| `MAX_PAPERS_PER_SOURCE` | `100` | Max papers to fetch per source |

### Optional - Tier Configuration

| Variable | Default | Description |
|----------|---------|-------------|
| `ENABLE_TIER1` | `true` | Enable keyword filtering |
| `ENABLE_TIER2` | `true` | Enable LLM triage |
| `ENABLE_TIER3` | `true` | Enable expert extraction |
| `TIER1_MIN_KEYWORDS` | `2` | Minimum keyword matches to pass Tier 1 |
| `TIER2_MODEL` | `gpt-4o-mini` | Model for Tier 2 classification |
| `TIER2_CONFIDENCE_THRESHOLD` | `0.5` | Minimum confidence to pass Tier 2 |
| `TIER3_MODELS` | `gpt-4o-mini,gpt-4o` | Comma-separated model cascade |
| `TIER3_THRESHOLD` | `1` | Retry next model if variants < threshold |
| `TIER3_TEMPERATURE` | `0.0` | Temperature for extraction |
| `TIER3_MAX_TOKENS` | `16000` | Max tokens for extraction response |

### Optional - Data Scout

| Variable | Default | Description |
|----------|---------|-------------|
| `SCOUT_ENABLED` | `true` | Generate DATA_ZONES.md files |
| `SCOUT_MIN_RELEVANCE` | `0.3` | Minimum zone relevance (0.0-1.0) |
| `SCOUT_MAX_ZONES` | `30` | Maximum zones per paper |
| `SCOUT_USE_CONDENSED` | `true` | Prefer DATA_ZONES.md for extraction |

---

## Troubleshooting

### "No papers were successfully downloaded"

- Many papers are paywalled; check `pmc_fulltext/paywalled_missing.csv`
- Use `fetch_manager.py` for semi-manual recovery of important papers
- Some PMIDs may not have PMC full-text available

### "Extraction failed for PMID X"

- Check `pmid_status/extraction_failures.csv` for error details
- Large papers may exceed token limits—try increasing `TIER3_MAX_TOKENS`
- Malformed papers may fail parsing—check `pmc_fulltext/` logs

### "All papers filtered out before download"

- Tier 1/2 filters may be too strict for your gene
- Disable filters temporarily: `ENABLE_TIER1=false ENABLE_TIER2=false python automated_workflow.py ...`
- Lower confidence threshold: `TIER2_CONFIDENCE_THRESHOLD=0.3`
- Check `pmid_status/filtered_out.csv` to understand filter decisions

### Rate limiting errors

- Set `NCBI_API_KEY` for higher NCBI rate limits (10 req/s vs 3 req/s)
- Default delays are 2 seconds between downloads
- OpenAI rate limits may require spreading large jobs over time

### "Missing required settings" error

- Ensure all required environment variables are set (not placeholder values)
- Check your `.env` file or export variables directly
- Required: `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `NCBI_EMAIL`

---

## Development

See [CLAUDE.md](CLAUDE.md) for architecture details.

```bash
# Run tests
pytest tests/

# Example scripts
python tests/example_harvest_from_pubmind.py
python tests/example_triage.py
```
