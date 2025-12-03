# CLAUDE.md - AI Assistant Guide for GeneVariantFetcher

This document provides comprehensive guidance for AI assistants working on the GeneVariantFetcher codebase. It explains the project structure, development workflows, conventions, and best practices.

## Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Architecture & Data Flow](#architecture--data-flow)
- [Key Modules & Components](#key-modules--components)
- [Development Workflows](#development-workflows)
- [Configuration System](#configuration-system)
- [Testing Strategy](#testing-strategy)
- [Code Conventions](#code-conventions)
- [Common Tasks](#common-tasks)
- [Important Files Reference](#important-files-reference)
- [Troubleshooting Guidelines](#troubleshooting-guidelines)

---

## Project Overview

**GeneVariantFetcher** is a tiered biomedical extraction pipeline that automates the complete workflow from gene name to SQLite database containing genetic variant data, patient-level information, and penetrance statistics.

### What It Does

1. **Discovers papers** - Finds relevant papers using PubMind/PubMed
2. **Downloads full-text** - Retrieves complete articles with supplemental materials
3. **Extracts data** - Uses AI to extract patient-level variants and phenotypes
4. **Aggregates penetrance** - Summarizes variant penetrance across papers
5. **Stores in SQLite** - Saves normalized data for easy querying

### Key Features

- **PubMind-First Strategy**: Prioritizes variant-focused literature over broad searches
- **Tiered Filtering**: Progressive cost reduction (85% savings vs. direct extraction)
- **Full-Text + Supplements**: Processes Excel/Word files, not just abstracts
- **Configuration-Driven**: All behavior controlled via environment variables
- **SQLite Integration**: Built-in migration, no separate tools needed

### Technology Stack

- **Python**: 3.11+ required
- **LLM APIs**: OpenAI (GPT-4o, GPT-4o-mini) or Anthropic (Claude)
- **Data Sources**: PubMind, PubMed/NCBI, EuropePMC
- **Key Libraries**: litellm, biopython, pydantic, tenacity, markitdown

---

## Repository Structure

```
GeneVariantFetcher/
├── 📄 Root-Level Entry Points
│   ├── automated_workflow.py ......... MAIN ENTRY: Complete end-to-end pipeline
│   ├── main.py ....................... Simple help message
│   ├── config.py ..................... Configuration singleton with validation
│   ├── models.py ..................... Pydantic data models (Paper, FilterResult, etc.)
│   │
├── 🧠 Core Pipeline Modules
│   ├── pipeline/
│   │   ├── sourcing.py ............... Stage 1: Paper discovery (PubMind PRIMARY)
│   │   ├── filters.py ................ Stage 2: Tiered filtering
│   │   │                               - Tier 1: KeywordFilter (fast, free)
│   │   │                               - Tier 2: InternFilter (gpt-4o-mini)
│   │   ├── extraction.py ............. Stage 3: Expert extraction (gpt-4o/claude)
│   │   ├── aggregation.py ............ Stage 4: Penetrance aggregation
│   │   └── cli.py .................... Typer-based CLI commands
│   │
├── 📥 Harvesting Modules
│   ├── harvesting/
│   │   ├── orchestrator.py ........... Main harvesting coordinator (PMCHarvester)
│   │   ├── pmc_api.py ................ NCBI PMC API client
│   │   ├── doi_resolver.py ........... DOI resolution with fallback
│   │   ├── supplement_scraper.py ..... Web scraping (Nature, Elsevier, generic)
│   │   └── format_converters.py ...... PDF/DOCX/XLSX to Markdown
│   │
├── 🛠️ Utility Modules
│   ├── utils/
│   │   ├── llm_utils.py .............. BaseLLMCaller for consistent LLM calls
│   │   ├── pubmed_utils.py ........... PubMed/Entrez API utilities
│   │   ├── retry_utils.py ............ Retry decorators with exponential backoff
│   │   └── html_utils.py ............. HTML parsing utilities
│   │
├── ⚙️ Configuration Package
│   ├── config/
│   │   ├── settings.py ............... Pydantic BaseSettings
│   │   └── test_settings.py .......... Configuration tests
│   │
├── 📋 Supporting Scripts
│   ├── pubmind_fetcher.py ............ PubMind API client
│   ├── migrate_to_sqlite.py .......... JSON → SQLite migration
│   ├── query_variants_db.py .......... Database query utility
│   ├── rerun_extraction.py ........... Re-run extraction helper
│   ├── example_*.py .................. Usage examples (6 files)
│   │
├── 🧪 Test Suite
│   ├── test_*.py ..................... Unit and integration tests (8 files)
│   │
├── 📚 Documentation
│   ├── README.md ..................... Quick start and overview
│   ├── ARCHITECTURE.md ............... System design and tiered pipeline
│   ├── SIMPLE_USAGE_GUIDE.md ......... Step-by-step workflows
│   ├── SQLITE_MIGRATION_GUIDE.md ..... Database schema
│   ├── IMPLEMENTATION_NOTES.md ....... DOI resolution details
│   ├── CLAUDE.md ..................... This file
│   │
├── 📊 Output Directory
│   └── automated_output/
│       └── {GENE}/{TIMESTAMP}/
│           ├── {GENE}_pmids.txt ...... Discovered PMIDs
│           ├── pmc_fulltext/ ......... Downloaded articles
│           ├── extractions/ .......... JSON extraction results
│           └── {GENE}.db ............. Final SQLite database
│
└── ⚙️ Configuration Files
    ├── pyproject.toml ................ Project metadata and dependencies
    ├── .env.example .................. Configuration template
    └── .gitignore .................... Git ignore patterns
```

---

## Architecture & Data Flow

### Sequential Pipeline Stages

```
Gene Symbol
    ↓
┌──────────────────────────────────────────────┐
│ STAGE 1: SOURCING                            │
│ Module: pipeline/sourcing.py                 │
│ Class: PaperSourcer                          │
│ - Query PubMind API (PRIMARY)                │
│ - Optional: PubMed/EuropePMC (fallback)      │
│ - Return deduplicated PMIDs                  │
└──────────────────────────────────────────────┘
    ↓ (Set of PMIDs)
┌──────────────────────────────────────────────┐
│ STAGE 2: TIERED FILTERING                    │
│ Module: pipeline/filters.py                  │
│ Tier 1: KeywordFilter                        │
│   - Fast keyword matching (~1000 papers/sec) │
│   - No API calls (FREE)                      │
│   - Filters ~40-60% of papers                │
│ Tier 2: InternFilter                         │
│   - LLM classification (gpt-4o-mini)         │
│   - Cost: ~$0.0001/paper                     │
│   - Identifies clinical data                 │
│   - Filters ~70% of remaining papers         │
└──────────────────────────────────────────────┘
    ↓ (Filtered PMIDs)
┌──────────────────────────────────────────────┐
│ STAGE 2b: FULL-TEXT HARVESTING               │
│ Module: harvesting/orchestrator.py           │
│ Class: PMCHarvester                          │
│ - PMID → PMCID conversion                    │
│ - Download full-text XML from PMC            │
│ - Download supplemental files                │
│ - Convert PDF/DOCX/XLSX to Markdown          │
│ - Create unified {PMID}_FULL_CONTEXT.md      │
└──────────────────────────────────────────────┘
    ↓ (Markdown files)
┌──────────────────────────────────────────────┐
│ STAGE 3: EXPERT EXTRACTION                   │
│ Module: pipeline/extraction.py               │
│ Class: ExpertExtractor                       │
│ - Uses gpt-4o or claude-3-opus               │
│ - Extracts variants (HGVS notation)          │
│ - Extracts patient phenotypes                │
│ - Extracts penetrance data                   │
│ - Cost: ~$0.05/paper                         │
└──────────────────────────────────────────────┘
    ↓ (JSON extraction files)
┌──────────────────────────────────────────────┐
│ STAGE 4: DATA AGGREGATION                    │
│ Module: pipeline/aggregation.py              │
│ Class: DataAggregator                        │
│ - Validate extracted data                    │
│ - Combine data across papers                 │
│ - Calculate penetrance statistics            │
│ - Generate aggregated JSON                   │
└──────────────────────────────────────────────┘
    ↓ (Aggregated JSON)
┌──────────────────────────────────────────────┐
│ STAGE 5: SQLITE MIGRATION                    │
│ Script: migrate_to_sqlite.py                 │
│ - Parse JSON files                           │
│ - Create normalized schema                   │
│ - Insert into SQLite database                │
│ - Final output: {GENE}.db                    │
└──────────────────────────────────────────────┘
```

### Cost Optimization Through Tiering

```
Example: Processing 1000 papers

WITHOUT FILTERING:
1000 papers × $0.10/paper = $100.00

WITH TIERED FILTERING:
├─ Tier 1 (KeywordFilter): 1000 papers × $0.00 = $0.00
│  └─ 500 papers FAIL (50%), 500 papers PASS
├─ Tier 2 (InternFilter): 500 papers × $0.0001 = $0.05
│  └─ 350 papers FAIL (70%), 150 papers PASS
└─ Tier 3 (ExpertExtractor): 150 papers × $0.10 = $15.00

TOTAL COST: $15.05 (85% savings)
```

---

## Key Modules & Components

### 1. Paper Sourcing (`pipeline/sourcing.py`)

**Purpose:** Discover relevant papers for a gene

**Main Class:** `PaperSourcer`

**Key Methods:**
- `fetch_papers(gene_symbol, max_results_per_source, use_pubmind, use_pubmed, use_europepmc)` - Query configured sources and return deduplicated PMIDs
- `fetch_paper_metadata(pmid)` - Get paper title, abstract, authors, DOI

**Configuration:**
- `USE_PUBMIND=true` - Use PubMind (recommended)
- `PUBMIND_ONLY=true` - Ignore other sources (default)
- `USE_PUBMED=false` - Disable PubMed (default)
- `USE_EUROPEPMC=false` - Disable EuropePMC (default)
- `MAX_PAPERS_PER_SOURCE=100` - Limit papers per source

**PubMind-First Strategy:**
- PubMind provides variant-focused literature (higher relevance)
- Reduces noise compared to broad PubMed searches
- Better cost efficiency (fewer irrelevant papers)

---

### 2. Tiered Filtering (`pipeline/filters.py`)

**Purpose:** Progressively filter papers to reduce expensive LLM calls

#### Tier 1: KeywordFilter

**Purpose:** Fast elimination of obviously irrelevant papers

**Implementation:**
- Keyword/regex matching on abstracts
- No API calls (completely free)
- Speed: ~1000 papers/second

**Configuration:**
- `ENABLE_TIER1=true` - Enable keyword filtering
- `TIER1_MIN_KEYWORDS=2` - Minimum matches to pass
- `TIER1_USE_LLM=false` - Use keywords (not LLM)

**Default Keywords:** variant, mutation, clinical, patient, case, cohort, carrier, affected, unaffected, penetrance, phenotype, genotype, pathogenic, benign, etc. (50+ terms)

**Returns:** `FilterResult` with `decision=PASS/FAIL`

#### Tier 2: InternFilter (LLM Classification)

**Purpose:** Identify papers with original clinical data

**Implementation:**
- Uses cheap LLM (gpt-4o-mini by default)
- Classifies papers as clinical vs. non-clinical
- Confidence threshold filtering

**Configuration:**
- `ENABLE_TIER2=true` - Enable LLM classification
- `TIER2_MODEL=gpt-4o-mini` - Model choice
- `TIER2_TEMPERATURE=0.1` - Low temperature for consistency
- `TIER2_MAX_TOKENS=150` - Limit response size
- `TIER2_CONFIDENCE_THRESHOLD=0.5` - Minimum confidence to pass

**Cost:** ~$0.0001 per paper

**What Gets KEPT:**
- Case reports
- Case series
- Clinical cohorts with patient data
- Functional studies with phenotypes

**What Gets DROPPED:**
- Review articles
- Meta-analyses
- Animal studies only
- Pure computational predictions

**Returns:** `FilterResult` with `confidence` score (0.0-1.0)

---

### 3. Full-Text Harvesting (`harvesting/orchestrator.py`)

**Purpose:** Download complete articles with supplemental materials

**Main Class:** `PMCHarvester`

**Key Methods:**
- `harvest(pmids, delay=2.0)` - Download full-text for list of PMIDs
- `_get_pmc_id(pmid)` - Convert PMID to PMCID
- `_download_full_text(pmcid)` - Get article XML from PMC
- `_download_supplements(pmcid)` - Get supplemental files
- `_convert_to_markdown(file_path)` - Convert to Markdown

**Workflow:**
1. PMID → PMCID conversion (via NCBI API)
2. Download full-text XML from PMC
3. Download supplemental files (Excel, Word, PDF)
4. Convert all files to Markdown
5. Create unified `{PMID}_FULL_CONTEXT.md` file

**Fallback Strategy:**
- **Step 1:** Try PMC Open Access API
- **Step 2:** Try DOI resolution (doi.org)
- **Step 3:** Try web scraping (Nature, Elsevier)
- **Step 4:** Construct direct publisher URL

**Supplemental Scraping:**
- `harvesting/supplement_scraper.py` contains domain-specific scrapers
- Nature: Looks for `<div id="supplementary-information">`
- Elsevier: Three-tier strategy (CSS classes, JSON data, regex patterns)
- Generic: Keyword-based detection

**Rate Limiting:**
- Default delay: 2 seconds between requests
- Respects NCBI guidelines (3 req/sec without API key, 10 req/sec with key)

**Output Structure:**
```
automated_output/{GENE}/{TIMESTAMP}/pmc_fulltext/
├── {PMID}_FULL_CONTEXT.md          # Unified article + supplements
├── {PMID}_supplements/              # Raw supplemental files
│   ├── Table_S1.xlsx
│   └── Figure_S2.pdf
├── successful_downloads.csv         # Success log
└── paywalled_missing.csv           # Unavailable articles log
```

**Important:** Only ~30% of PubMed articles are available full-text in PMC

---

### 4. Expert Extraction (`pipeline/extraction.py`)

**Purpose:** Extract structured variant and phenotype data using AI

**Main Class:** `ExpertExtractor` (extends `BaseLLMCaller`)

**Key Methods:**
- `extract(paper)` - Extract data from a Paper object
- `_build_extraction_prompt(paper)` - Create extraction prompt
- `_parse_extraction_response(response)` - Parse LLM JSON output

**Configuration:**
- `ENABLE_TIER3=true` - Enable extraction
- `TIER3_MODEL=gpt-4o` - Model choice (or claude-3-opus)
- `TIER3_TEMPERATURE=0.0` - Deterministic extraction
- `TIER3_MAX_TOKENS=8000` - Allow detailed responses

**Model Options:**
- `gpt-4o` - Balanced cost/performance (recommended)
- `gpt-4o-mini` - Cheaper, less accurate
- `claude-3-opus-20240229` - Higher accuracy
- `claude-3-sonnet-20240229` - Good balance for Claude

**Extracted Data:**
- **Variants:** HGVS notation (cDNA, protein, genomic)
- **Clinical Significance:** Pathogenic, benign, VUS
- **Patient Demographics:** Age, sex, ancestry
- **Phenotypes:** Clinical features with HPO codes
- **Penetrance Data:** Affected vs. unaffected carriers
- **Functional Studies:** In vitro/in vivo results
- **Segregation:** Family co-segregation data
- **Evidence:** Sentences supporting each extraction

**Output Format:**
```json
{
  "pmid": "35443093",
  "individuals": [
    {
      "individual_id": "P1",
      "variants": ["c.1234G>A", "p.Arg412Gln"],
      "phenotypes": ["Cardiomyopathy", "Arrhythmia"],
      "clinical_status": "affected",
      "age": "45",
      "sex": "male",
      "evidence": "Patient P1 presented with..."
    }
  ],
  "penetrance_data": {
    "variant": "c.1234G>A",
    "affected_carriers": 15,
    "unaffected_carriers": 3,
    "penetrance": 0.833
  }
}
```

**Cost:** ~$0.05 per paper (varies with model and paper length)

**Text Truncation:**
- Abstracts: 2000 characters
- Full text: 30,000 characters
- Reduces cost without losing key information

---

### 5. Data Aggregation (`pipeline/aggregation.py`)

**Purpose:** Validate and aggregate extracted data across papers

**Main Class:** `DataAggregator`

**Key Methods:**
- `aggregate(extraction_files)` - Combine data from multiple JSON files
- `validate_penetrance_data(data)` - Check data quality
- `calculate_statistics(variant_data)` - Compute penetrance stats

**Operations:**
- Validates extracted data completeness
- Deduplicates variants across papers
- Aggregates penetrance statistics (affected/unaffected ratios)
- Calculates confidence intervals
- Generates penetrance reports

**Output:** `penetrance_aggregated.json` with cross-paper summaries

---

### 6. SQLite Migration (`migrate_to_sqlite.py`)

**Purpose:** Convert JSON extractions to normalized SQLite database

**Main Function:** `migrate_to_sqlite(json_dir, output_db_path)`

**Database Schema:**
- `papers` table - Paper metadata
- `variants` table - Genetic variants
- `phenotypes` table - Clinical phenotypes
- `individuals` table - Patient records
- `penetrance` table - Penetrance statistics
- Junction tables for many-to-many relationships

**Usage:**
```bash
python migrate_to_sqlite.py automated_output/BRCA1/{timestamp}/extractions/ BRCA1.db
```

**Features:**
- Automatic schema creation
- Referential integrity constraints
- Indexed columns for fast queries
- Compressed storage

---

### 7. Shared Utilities

#### `utils/llm_utils.py` - LLM Integration

**Main Class:** `BaseLLMCaller`

**Purpose:** Consistent LLM calling across all pipeline stages

**Key Methods:**
- `call_llm_json(prompt)` - Call LLM and parse JSON response
- `_retry_llm_call(prompt)` - Automatic retry with exponential backoff
- `parse_llm_json_response(response)` - Clean and parse JSON

**Features:**
- Unified API across OpenAI and Anthropic
- Automatic retry on failures (3 attempts)
- JSON parsing with error recovery
- Token usage tracking

**Usage Pattern:**
```python
from utils.llm_utils import BaseLLMCaller

class MyExtractor(BaseLLMCaller):
    def __init__(self):
        super().__init__(model="gpt-4o", temperature=0.0)

    def extract(self, text):
        prompt = f"Extract data from: {text}"
        return self.call_llm_json(prompt)
```

#### `utils/pubmed_utils.py` - PubMed API

**Key Functions:**
- `query_pubmed_with_entrez(gene_symbol, max_results)` - Standard PubMed query
- `query_pubmed_for_gene(gene_symbol)` - Gene-specific query builder
- `fetch_paper_metadata(pmid)` - Get paper metadata
- `fetch_paper_abstract(pmid)` - Get abstract text
- `get_doi_from_pmid(pmid)` - Resolve DOI
- `query_europepmc(gene_symbol, max_results)` - Europe PMC query
- `batch_fetch_metadata(pmids)` - Efficient bulk retrieval

**Important:** All functions require `NCBI_EMAIL` to be set

#### `utils/retry_utils.py` - Retry Logic

**Pre-configured Decorators:**
- `@standard_retry` - 3 attempts, exponential backoff (1s, 2s, 4s)
- `@api_retry` - Optimized for API calls
- `@llm_retry` - Longer backoff for rate limits (2s, 4s, 8s)
- `@scraping_retry` - Fewer retries (2 attempts)

**Usage:**
```python
from utils.retry_utils import api_retry

@api_retry
def fetch_data():
    return requests.get("https://api.example.com")
```

---

## Development Workflows

### Setting Up Development Environment

```bash
# Clone repository
git clone https://github.com/kroncke-lab/GeneVariantFetcher.git
cd GeneVariantFetcher

# Create virtual environment (Python 3.11+ required)
python3.11 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -e .

# Copy environment template
cp .env.example .env

# Edit .env with your API keys
nano .env  # Or use your preferred editor
```

### Required Environment Variables

**Minimum Configuration:**
```bash
# .env file
OPENAI_API_KEY=your-openai-api-key-here
NCBI_EMAIL=your_email@example.com
```

**Recommended Configuration:**
```bash
# API Keys
OPENAI_API_KEY=your-openai-api-key-here
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=your-ncbi-api-key-here  # Optional, increases rate limits

# Model Selection
TIER2_MODEL=gpt-4o-mini
TIER3_MODEL=gpt-4o

# Paper Sourcing (PubMind-first)
USE_PUBMIND=true
PUBMIND_ONLY=true
USE_PUBMED=false
USE_EUROPEPMC=false

# Tier Configuration
ENABLE_TIER1=true
ENABLE_TIER2=true
ENABLE_TIER3=true
TIER1_MIN_KEYWORDS=2
TIER2_CONFIDENCE_THRESHOLD=0.5
```

### Running the Main Workflow

**Complete End-to-End Pipeline:**
```bash
python automated_workflow.py BRCA1 --email your@email.com
```

**With Custom Limits:**
```bash
python automated_workflow.py SCN5A \
  --email your@email.com \
  --max-pmids 200 \
  --max-downloads 100
```

**Step-by-Step Manual Control:**
```python
# Step 1: Source papers
from pipeline.sourcing import PaperSourcer
sourcer = PaperSourcer()
pmids = sourcer.fetch_papers("BRCA1", max_results_per_source=100)

# Step 2: Harvest full-text
from harvesting import PMCHarvester
harvester = PMCHarvester(output_dir="output/pmc_fulltext")
harvester.harvest(pmids, delay=2.0)

# Step 3: Extract data
from pipeline.extraction import ExpertExtractor
extractor = ExpertExtractor(model="gpt-4o")
# Process harvested markdown files...

# Step 4: Aggregate
from pipeline.aggregation import DataAggregator
aggregator = DataAggregator()
# Combine extraction results...

# Step 5: Migrate to SQLite
python migrate_to_sqlite.py output/extractions/ BRCA1.db
```

### Branching Strategy

**Branch Naming Convention:**
- `main` - Stable production code
- `claude/{description}-{session-id}` - Claude-generated branches
- `codex/{description}` - Codex-generated branches

**Example:** `claude/add-claude-documentation-014ge91B1oC8dhxv88gQPv9f`

**Important:** When pushing, branch names must:
- Start with `claude/`
- End with matching session ID
- Otherwise push will fail with 403 HTTP code

### Commit Message Style

**Observed Patterns:**
```
Add method to sort non-null PMIDs and refactor affected status checks
Fix PMID parsing in automated workflow extractions
Simplify workflow: Remove duplicates, integrate SQLite, streamline to single pipeline
Add `rerun_extraction.py` to support re-extracting variant data from existing workflow outputs.
Create abstraction layer for extractors
```

**Conventions:**
- Start with imperative verb (Add, Fix, Simplify, Create, Update, Remove)
- Be specific and descriptive
- Use present tense
- Don't use periods at the end
- Reference functions/files with backticks when helpful

**Good Examples:**
```
Add retry logic to DOI resolution with exponential backoff
Fix keyword filtering to handle multi-line abstracts
Update ExpertExtractor to support Claude models
Refactor harvesting module into separate submodules
```

**Bad Examples:**
```
Fixed stuff
Updated code
Changes
WIP
```

### Making Changes

**When Adding New Features:**
1. Read relevant documentation (ARCHITECTURE.md, README.md)
2. Understand the module's purpose and dependencies
3. Follow existing code patterns (see Code Conventions)
4. Add tests if applicable
5. Update docstrings
6. Test thoroughly before committing

**When Fixing Bugs:**
1. Reproduce the issue
2. Identify root cause
3. Check if similar patterns exist elsewhere
4. Fix all occurrences
5. Add test case if missing
6. Verify fix doesn't break other functionality

**When Refactoring:**
1. Understand current implementation completely
2. Ensure changes maintain backward compatibility
3. Update all callers if API changes
4. Run all tests
5. Update documentation

---

## Configuration System

### Dual Configuration Architecture

The project uses two configuration systems:

#### 1. Root-Level Config (`config.py`)

**Purpose:** Simple, legacy configuration singleton

**Usage:**
```python
from config import config

# Access settings
openai_key = config.get("OPENAI_API_KEY")
model = config.get("TIER2_MODEL", "gpt-4o-mini")

# Validate required settings
config.validate()
```

**Key Methods:**
- `config.get(key, default=None)` - Get setting with optional default
- `config.validate()` - Ensure required keys are set
- `config.as_dict()` - Export all settings as dictionary

#### 2. Package-Level Config (`config/settings.py`)

**Purpose:** Advanced Pydantic-based configuration with validation

**Usage:**
```python
from config.settings import get_settings

settings = get_settings()  # Cached singleton

# Type-safe access
print(settings.tier2_model)  # gpt-4o-mini
print(settings.pubmind_only)  # True
print(settings.tier2_confidence_threshold)  # 0.5
```

**Key Features:**
- Type validation via Pydantic
- Environment variable loading
- Default values
- Cached singleton pattern

**Settings Class:**
```python
class Settings(BaseSettings):
    # API Keys
    openai_api_key: str | None
    anthropic_api_key: str | None
    ncbi_email: str | None
    ncbi_api_key: str | None

    # Paper Sourcing
    use_pubmind: bool = True
    pubmind_only: bool = True
    use_pubmed: bool = False
    use_europepmc: bool = False
    max_papers_per_source: int = 100

    # Tiered Classification
    enable_tier1: bool = True
    enable_tier2: bool = True
    enable_tier3: bool = True

    tier1_min_keywords: int = 2
    tier1_use_llm: bool = False

    tier2_model: str = "gpt-4o-mini"
    tier2_temperature: float = 0.1
    tier2_max_tokens: int = 150
    tier2_confidence_threshold: float = 0.5

    tier3_model: str = "gpt-4o"
    tier3_temperature: float = 0.0
    tier3_max_tokens: int = 8000

    # Logging
    log_level: str = "INFO"

    # Retry Settings
    max_retries: int = 3
    retry_min_wait: int = 2
    retry_max_wait: int = 10
```

### Configuration Hierarchy

1. **Environment Variables** (highest priority)
2. **`.env` file** (loaded automatically)
3. **Default values** (in Settings class)

### Adding New Configuration Options

**Step 1:** Add to `.env.example`
```bash
# New feature configuration
NEW_FEATURE_ENABLED=true
NEW_FEATURE_THRESHOLD=0.8
```

**Step 2:** Add to `config/settings.py`
```python
class Settings(BaseSettings):
    # Existing settings...

    # New feature settings
    new_feature_enabled: bool = True
    new_feature_threshold: float = 0.8
```

**Step 3:** Document in comments
```python
class Settings(BaseSettings):
    """Application settings loaded from environment variables.

    All settings can be overridden via environment variables.
    For example, NEW_FEATURE_ENABLED=false will disable the feature.
    """
```

**Step 4:** Use in code
```python
from config.settings import get_settings

settings = get_settings()
if settings.new_feature_enabled:
    # Use the new feature
    threshold = settings.new_feature_threshold
```

---

## Testing Strategy

### Test Organization

**Test Files (8 total):**
- `test_scraper.py` - Supplement scraper unit tests
- `test_pubmed_utils.py` - PubMed API utility tests
- `test_pubmed_email_propagation.py` - Email configuration tests
- `test_with_pmids_file.py` - Full workflow integration test
- `test_doi_resolution.py` - DOI resolution tests
- `test_triage.py` - Clinical triage filter tests
- `test_sqlite_migration.py` - Database migration tests
- `config/test_settings.py` - Configuration loading tests

### Running Tests

**Run All Tests:**
```bash
pytest
```

**Run Specific Test:**
```bash
pytest test_triage.py
pytest test_pubmed_utils.py -v  # Verbose output
```

**Run with Coverage:**
```bash
pytest --cov=pipeline --cov=harvesting --cov=utils
```

### Test Patterns

**Unit Test Example:**
```python
def test_keyword_filter_passes_relevant_paper():
    """Test that papers with clinical keywords pass Tier 1 filter."""
    filter = KeywordFilter(min_keywords=2)
    paper = Paper(
        pmid="12345678",
        title="Novel BRCA1 variant in breast cancer patient",
        abstract="We report a case of a 32-year-old woman with breast cancer..."
    )
    result = filter.filter(paper)
    assert result.decision == FilterDecision.PASS
    assert result.tier == FilterTier.TIER_1
```

**Integration Test Example:**
```python
def test_full_workflow_with_pmid_file():
    """Test complete workflow from PMID file to SQLite database."""
    # Setup
    pmid_file = "test_pmids.txt"
    output_dir = "test_output"

    # Run workflow
    result = run_automated_workflow(
        gene_symbol="BRCA1",
        pmid_file=pmid_file,
        output_dir=output_dir
    )

    # Verify
    assert result.success
    assert os.path.exists(f"{output_dir}/BRCA1.db")
```

### Testing Best Practices

**DO:**
- Write tests for new features
- Test edge cases and error conditions
- Use descriptive test names
- Mock external API calls when possible
- Clean up test artifacts in teardown

**DON'T:**
- Commit failing tests
- Test implementation details
- Make tests depend on external services
- Leave test data in repository

---

## Code Conventions

### Python Style

**Follow PEP 8:**
- 4 spaces for indentation (no tabs)
- Maximum line length: 100 characters (flexible for readability)
- Two blank lines between top-level definitions
- One blank line between methods

**Type Hints:**
```python
def fetch_papers(
    gene_symbol: str,
    max_results: int = 100,
    use_pubmind: bool = True
) -> list[str]:
    """Fetch PMIDs for a gene symbol.

    Args:
        gene_symbol: Gene name (e.g., "BRCA1")
        max_results: Maximum papers to return
        use_pubmind: Whether to query PubMind

    Returns:
        List of PMIDs as strings
    """
    ...
```

**Use Modern Python Features:**
- Type hints with `|` for unions (Python 3.10+)
- f-strings for formatting
- Pathlib for file paths
- Context managers for resources
- List/dict comprehensions where readable

### Docstring Format

**Module-Level:**
```python
"""Full-text harvesting orchestrator.

This module coordinates downloading complete articles with supplemental
materials from PubMed Central and publisher websites.

Classes:
    PMCHarvester: Main harvesting class

Functions:
    convert_pmid_to_pmcid: Convert PMID to PMCID via NCBI API
"""
```

**Class-Level:**
```python
class PMCHarvester:
    """Downloads full-text articles and supplements from PubMed Central.

    This class handles the complete workflow of:
    1. Converting PMIDs to PMCIDs
    2. Downloading full-text XML
    3. Downloading supplemental files
    4. Converting to Markdown

    Attributes:
        output_dir: Directory for downloaded files
        session: HTTP session with browser headers

    Example:
        >>> harvester = PMCHarvester(output_dir="pmc_harvest")
        >>> harvester.harvest(["35443093", "33442691"], delay=2.0)
    """
```

**Function-Level:**
```python
def harvest(self, pmids: list[str], delay: float = 2.0) -> dict[str, bool]:
    """Download full-text articles for a list of PMIDs.

    Args:
        pmids: List of PubMed IDs to harvest
        delay: Seconds to wait between requests (default: 2.0)

    Returns:
        Dictionary mapping PMID to success status

    Raises:
        ValueError: If pmids list is empty
        ConnectionError: If API is unreachable

    Note:
        Only ~30% of PubMed articles are available full-text in PMC.
        Unavailable articles are logged to paywalled_missing.csv.
    """
```

### Error Handling

**Use Specific Exceptions:**
```python
# Good
try:
    response = requests.get(url, timeout=10)
    response.raise_for_status()
except requests.HTTPError as e:
    logger.error(f"HTTP error fetching {url}: {e}")
    raise
except requests.Timeout:
    logger.warning(f"Timeout fetching {url}, will retry")
    raise

# Bad
try:
    response = requests.get(url)
except Exception:
    pass
```

**Use Retry Decorators:**
```python
from utils.retry_utils import api_retry

@api_retry
def fetch_from_api(url: str) -> dict:
    """Fetch data from API with automatic retry."""
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return response.json()
```

**Fail Gracefully:**
```python
def download_supplement(url: str) -> bytes | None:
    """Download supplemental file, return None on failure."""
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.content
    except requests.RequestException as e:
        logger.warning(f"Failed to download supplement {url}: {e}")
        return None
```

### Logging

**Use Standard Logging Module:**
```python
import logging

logger = logging.getLogger(__name__)

# Log levels
logger.debug("Detailed information for debugging")
logger.info("General informational messages")
logger.warning("Warning messages (non-critical issues)")
logger.error("Error messages (recoverable failures)")
logger.critical("Critical messages (unrecoverable failures)")
```

**Log Format:**
```python
# Good: Structured, informative
logger.info(f"Fetched {len(pmids)} PMIDs for gene {gene_symbol}")
logger.warning(f"No PMCID found for PMID {pmid}, skipping")
logger.error(f"Failed to parse extraction result for {pmid}: {error}")

# Bad: Vague, unhelpful
logger.info("Success")
logger.error("Error occurred")
```

### Data Models

**Use Pydantic for Validation:**
```python
from pydantic import BaseModel, Field, validator

class Paper(BaseModel):
    """Represents a scientific paper."""

    pmid: str = Field(..., description="PubMed ID")
    title: str = Field(..., description="Paper title")
    abstract: str = Field(default="", description="Abstract text")
    full_text: str = Field(default="", description="Full-text content")
    doi: str | None = Field(default=None, description="Digital Object Identifier")

    @validator('pmid')
    def validate_pmid(cls, v):
        """Ensure PMID is numeric."""
        if not v.isdigit():
            raise ValueError(f"PMID must be numeric, got: {v}")
        return v

    class Config:
        """Pydantic configuration."""
        validate_assignment = True
        extra = "forbid"
```

### LLM Integration Pattern

**Always Extend BaseLLMCaller:**
```python
from utils.llm_utils import BaseLLMCaller

class MyExtractor(BaseLLMCaller):
    """Custom extractor using LLM."""

    def __init__(self, model: str = "gpt-4o", temperature: float = 0.0):
        """Initialize extractor.

        Args:
            model: LLM model identifier
            temperature: Sampling temperature (0.0-1.0)
        """
        super().__init__(model=model, temperature=temperature)

    def extract(self, text: str) -> dict:
        """Extract data from text.

        Args:
            text: Input text to process

        Returns:
            Extracted data as dictionary
        """
        prompt = self._build_prompt(text)
        result = self.call_llm_json(prompt)
        return result

    def _build_prompt(self, text: str) -> str:
        """Build extraction prompt."""
        return f"""Extract genetic variants from this text:

{text}

Return JSON with this format:
{{
    "variants": ["c.1234G>A", "p.Arg412Gln"],
    "phenotypes": ["breast cancer", "early onset"]
}}"""
```

---

## Common Tasks

### Task 1: Add a New Filter Tier

**Scenario:** Add a Tier 2b filter for detecting cohort studies

**Steps:**

1. **Create filter class in `pipeline/filters.py`:**
```python
class CohortFilter(BaseLLMCaller):
    """Identifies papers with cohort study design."""

    def __init__(self, model: str = "gpt-4o-mini", temperature: float = 0.1):
        super().__init__(model=model, temperature=temperature)

    def filter(self, paper: Paper) -> FilterResult:
        """Filter paper based on cohort study design."""
        prompt = f"""Is this a cohort study?

Title: {paper.title}
Abstract: {paper.abstract}

Return JSON: {{"is_cohort": true/false, "confidence": 0.0-1.0, "reason": "..."}}
"""
        result = self.call_llm_json(prompt)

        decision = FilterDecision.PASS if result["is_cohort"] else FilterDecision.FAIL
        return FilterResult(
            decision=decision,
            tier=FilterTier.TIER_2,  # Or create TIER_2B
            reason=result["reason"],
            confidence=result["confidence"]
        )
```

2. **Add to pipeline workflow:**
```python
# In automated_workflow.py or pipeline orchestrator
cohort_filter = CohortFilter(model=settings.tier2_model)
result = cohort_filter.filter(paper)
if result.decision == FilterDecision.FAIL:
    logger.info(f"Paper {paper.pmid} filtered by cohort filter")
    continue
```

3. **Add configuration:**
```python
# In config/settings.py
class Settings(BaseSettings):
    # ...
    enable_cohort_filter: bool = True
    cohort_filter_threshold: float = 0.7
```

4. **Add tests:**
```python
# In test_filters.py
def test_cohort_filter_detects_cohort_study():
    filter = CohortFilter()
    paper = Paper(
        pmid="12345678",
        title="10-year follow-up of BRCA1 carriers",
        abstract="We followed 500 BRCA1 carriers for 10 years..."
    )
    result = filter.filter(paper)
    assert result.decision == FilterDecision.PASS
```

### Task 2: Add Support for a New Publisher

**Scenario:** Add scraping support for Wiley journals

**Steps:**

1. **Add scraper to `harvesting/supplement_scraper.py`:**
```python
def _scrape_wiley_supplements(self, url: str, soup: BeautifulSoup) -> list[str]:
    """Scrape supplemental files from Wiley journal pages.

    Args:
        url: Wiley article URL
        soup: BeautifulSoup object of the page

    Returns:
        List of supplement download URLs
    """
    supplements = []

    # Strategy 1: Look for supporting information section
    supp_section = soup.find("section", {"class": "article-section__supporting"})
    if supp_section:
        links = supp_section.find_all("a", href=True)
        for link in links:
            href = link["href"]
            if any(ext in href.lower() for ext in [".pdf", ".xlsx", ".docx"]):
                full_url = urljoin(url, href)
                supplements.append(full_url)

    # Strategy 2: Look for data availability section
    if not supplements:
        data_section = soup.find("div", {"id": "support-information"})
        if data_section:
            links = data_section.find_all("a", href=True)
            for link in links:
                # Similar processing...

    return supplements
```

2. **Add to domain router:**
```python
def _scrape_supplements_by_domain(self, url: str, html: str) -> list[str]:
    """Route to domain-specific scraper."""
    soup = BeautifulSoup(html, "html.parser")

    if "nature.com" in url:
        return self._scrape_nature_supplements(url, soup)
    elif "sciencedirect.com" in url or "gimjournal.org" in url:
        return self._scrape_elsevier_supplements(url, soup)
    elif "wiley.com" in url or "onlinelibrary.wiley.com" in url:
        return self._scrape_wiley_supplements(url, soup)  # NEW
    else:
        return self._scrape_generic_supplements(url, soup)
```

3. **Add tests:**
```python
def test_wiley_supplement_scraping():
    scraper = SupplementScraper()
    url = "https://onlinelibrary.wiley.com/doi/10.1002/example"
    # Mock HTML response with Wiley structure
    html = """<section class="article-section__supporting">
        <a href="/doi/suppl/10.1002/example/suppl_file/table_s1.xlsx">Table S1</a>
    </section>"""

    supplements = scraper._scrape_supplements_by_domain(url, html)
    assert len(supplements) > 0
    assert "table_s1.xlsx" in supplements[0]
```

### Task 3: Add a New Extraction Field

**Scenario:** Extract family history information

**Steps:**

1. **Update extraction prompt in `pipeline/extraction.py`:**
```python
def _build_extraction_prompt(self, paper: Paper) -> str:
    """Build extraction prompt with family history field."""
    return f"""Extract genetic variant data from this paper:

{paper.full_text}

Return JSON with this structure:
{{
    "individuals": [
        {{
            "individual_id": "P1",
            "variants": ["c.1234G>A"],
            "phenotypes": ["breast cancer"],
            "clinical_status": "affected",
            "age": "45",
            "sex": "female",
            "family_history": {{  # NEW FIELD
                "has_family_history": true,
                "affected_relatives": 3,
                "relationship": ["mother", "sister", "aunt"],
                "details": "Mother diagnosed at age 42..."
            }}
        }}
    ]
}}
"""
```

2. **Update data model in `models.py`:**
```python
class FamilyHistory(BaseModel):
    """Family history information."""
    has_family_history: bool
    affected_relatives: int | None = None
    relationship: list[str] = Field(default_factory=list)
    details: str = ""

class Individual(BaseModel):
    """Individual patient record."""
    individual_id: str
    variants: list[str]
    phenotypes: list[str]
    clinical_status: str
    age: str | None = None
    sex: str | None = None
    family_history: FamilyHistory | None = None  # NEW FIELD
```

3. **Update database schema in `migrate_to_sqlite.py`:**
```python
CREATE TABLE IF NOT EXISTS family_history (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    individual_id INTEGER NOT NULL,
    has_family_history BOOLEAN NOT NULL,
    affected_relatives INTEGER,
    relationship TEXT,  -- JSON array
    details TEXT,
    FOREIGN KEY (individual_id) REFERENCES individuals(id)
)
```

4. **Update migration logic:**
```python
def insert_family_history(conn, individual_id, family_history):
    """Insert family history record."""
    if not family_history:
        return

    cursor = conn.cursor()
    cursor.execute("""
        INSERT INTO family_history
        (individual_id, has_family_history, affected_relatives, relationship, details)
        VALUES (?, ?, ?, ?, ?)
    """, (
        individual_id,
        family_history.has_family_history,
        family_history.affected_relatives,
        json.dumps(family_history.relationship),
        family_history.details
    ))
```

### Task 4: Debug Failed Extractions

**Scenario:** ExpertExtractor failing to parse some papers

**Debugging Steps:**

1. **Check logs:**
```bash
# Look for error messages in output
grep "ERROR" automated_output/BRCA1/{timestamp}/pipeline.log

# Check which PMIDs failed
grep "Failed to extract" automated_output/BRCA1/{timestamp}/pipeline.log
```

2. **Inspect failing paper:**
```python
from harvesting import PMCHarvester

# Read the problematic paper
with open("automated_output/BRCA1/{timestamp}/pmc_fulltext/12345678_FULL_CONTEXT.md") as f:
    content = f.read()

print(f"Paper length: {len(content)} characters")
print(f"First 500 chars: {content[:500]}")
```

3. **Test extraction manually:**
```python
from pipeline.extraction import ExpertExtractor
from models import Paper

extractor = ExpertExtractor(model="gpt-4o")
paper = Paper(
    pmid="12345678",
    title="...",
    full_text=content
)

try:
    result = extractor.extract(paper)
    print(result)
except Exception as e:
    print(f"Extraction failed: {e}")
    import traceback
    traceback.print_exc()
```

4. **Common issues and fixes:**

**Issue: Paper too long (>30,000 characters)**
```python
# Truncate in extraction.py
full_text = paper.full_text[:30000]
```

**Issue: LLM returns invalid JSON**
```python
# Add more robust JSON parsing in llm_utils.py
def parse_llm_json_response(response: str) -> dict:
    # Remove markdown code blocks
    response = re.sub(r"```json\s*", "", response)
    response = re.sub(r"```\s*$", "", response)

    # Try parsing
    try:
        return json.loads(response)
    except json.JSONDecodeError:
        # Try to extract JSON from text
        match = re.search(r"\{.*\}", response, re.DOTALL)
        if match:
            return json.loads(match.group(0))
        raise
```

**Issue: Rate limit exceeded**
```python
# Increase retry delay in retry_utils.py
@retry(
    wait=wait_exponential(min=4, max=30),  # Longer waits
    stop=stop_after_attempt(5),  # More attempts
    retry=retry_if_exception_type((openai.RateLimitError,))
)
```

5. **Re-run failed extractions:**
```bash
python rerun_extraction.py \
  automated_output/BRCA1/{timestamp}/ \
  --failed-only
```

---

## Important Files Reference

### Must-Read Documentation

| File | Purpose | When to Read |
|------|---------|--------------|
| `README.md` | Quick start, workflows, features | First time setup |
| `ARCHITECTURE.md` | System design, tiered pipeline | Understanding architecture |
| `CLAUDE.md` | This file | Working with codebase |
| `.env.example` | Configuration options | Setting up environment |
| `SIMPLE_USAGE_GUIDE.md` | Step-by-step examples | Learning workflows |

### Key Source Files

| File | Lines | Purpose | When to Modify |
|------|-------|---------|----------------|
| `automated_workflow.py` | 460+ | Main entry point | Adding workflow steps |
| `pipeline/sourcing.py` | 300+ | Paper discovery | Adding data sources |
| `pipeline/filters.py` | 400+ | Tiered filtering | Adding filter logic |
| `pipeline/extraction.py` | 300+ | Expert extraction | Changing extraction logic |
| `harvesting/orchestrator.py` | 500+ | Full-text harvesting | Modifying download process |
| `migrate_to_sqlite.py` | 386 | SQLite migration | Changing database schema |
| `config/settings.py` | 150+ | Configuration | Adding config options |
| `utils/llm_utils.py` | 200+ | LLM integration | Changing LLM behavior |
| `models.py` | 200+ | Data models | Adding new fields |

### Configuration Files

| File | Purpose |
|------|---------|
| `pyproject.toml` | Project metadata, dependencies, packaging |
| `.env.example` | Environment variable template |
| `.env` | Actual environment variables (not in git) |
| `.gitignore` | Files to exclude from version control |

### Output Files

| File/Directory | Purpose |
|----------------|---------|
| `automated_output/{GENE}/{TIMESTAMP}/` | Workflow output directory |
| `{GENE}_pmids.txt` | Discovered PMIDs |
| `pmc_fulltext/` | Downloaded articles |
| `extractions/` | JSON extraction results |
| `{GENE}.db` | Final SQLite database |
| `successful_downloads.csv` | Success log |
| `paywalled_missing.csv` | Unavailable articles log |
| `pipeline.log` | Execution log |

---

## Troubleshooting Guidelines

### Common Issues

#### Issue 1: Import Errors

**Symptom:**
```
ModuleNotFoundError: No module named 'pipeline'
```

**Cause:** Package not installed or wrong directory

**Fix:**
```bash
# Ensure you're in repository root
cd GeneVariantFetcher

# Install in editable mode
pip install -e .

# Verify installation
python -c "import pipeline; print(pipeline.__file__)"
```

#### Issue 2: Missing API Keys

**Symptom:**
```
ValueError: OPENAI_API_KEY not found in environment
```

**Cause:** Environment variables not set

**Fix:**
```bash
# Check .env file exists
ls -la .env

# Check key is set
grep OPENAI_API_KEY .env

# If missing, copy from template
cp .env.example .env
nano .env  # Add your actual key
```

#### Issue 3: PubMed Rate Limiting

**Symptom:**
```
HTTP 429 Too Many Requests from NCBI
```

**Cause:** Too many requests without API key

**Fix:**
```bash
# Get NCBI API key from https://www.ncbi.nlm.nih.gov/account/settings/
# Add to .env
NCBI_API_KEY=your-key-here

# This increases limit from 3 req/sec to 10 req/sec
```

#### Issue 4: LLM JSON Parsing Errors

**Symptom:**
```
JSONDecodeError: Expecting value: line 1 column 1
```

**Cause:** LLM returned invalid JSON

**Fix:**
- Check `llm_utils.py` has robust JSON parsing
- Increase `TIER3_MAX_TOKENS` to allow complete responses
- Try different model (e.g., gpt-4o instead of gpt-4o-mini)
- Add explicit JSON format instructions in prompt

#### Issue 5: Empty Extraction Results

**Symptom:** Extraction JSON has empty or minimal data

**Possible Causes:**
- Paper doesn't contain individual-level data (only summary stats)
- Paper is review/meta-analysis (no original data)
- Full-text not available (using abstract only)
- Extraction prompt not specific enough

**Debugging:**
```python
# Check what was actually sent to LLM
from pipeline.extraction import ExpertExtractor

extractor = ExpertExtractor(model="gpt-4o")
paper = ...  # Load paper

# Add debug logging to see full prompt
prompt = extractor._build_extraction_prompt(paper)
print(prompt)  # Check if enough context included
```

**Fixes:**
- Ensure using full-text (not just abstracts)
- Check paper passed Tier 2 filter (should have clinical data)
- Adjust extraction prompt to be more specific
- Use more powerful model (gpt-4o or claude-3-opus)

#### Issue 6: Harvester Fails with 403 Forbidden

**Symptom:**
```
HTTP 403 Forbidden when downloading from publisher
```

**Cause:** Anti-bot protection from publisher websites

**Fix:**
- Normal for cloud/datacenter IPs
- Run from residential network if possible
- Check `paywalled_missing.csv` for logging
- These papers may need manual retrieval
- PMC Open Access should still work (not affected)

#### Issue 7: SQLite Migration Errors

**Symptom:**
```
sqlite3.IntegrityError: FOREIGN KEY constraint failed
```

**Cause:** Data inconsistency or schema mismatch

**Fix:**
```python
# Check extraction JSON structure
import json
with open("extractions/12345678_extraction.json") as f:
    data = json.load(f)
    print(json.dumps(data, indent=2))

# Validate against models
from models import ExtractionResult
result = ExtractionResult(**data)  # Will raise if invalid

# Re-run migration with fresh database
rm output.db
python migrate_to_sqlite.py extractions/ output.db
```

### Getting Help

**Before asking for help:**
1. Check this documentation
2. Read error messages carefully
3. Check logs (`pipeline.log`)
4. Search existing issues on GitHub
5. Try reproducing with minimal example

**When reporting issues:**
1. Provide complete error message
2. Include relevant configuration (`.env` settings)
3. Share logs if possible
4. Describe what you tried
5. Include code snippets if applicable

---

## Best Practices Summary

### DO:

- ✅ Read existing documentation before starting
- ✅ Follow PEP 8 style guidelines
- ✅ Use type hints consistently
- ✅ Write descriptive docstrings
- ✅ Add tests for new features
- ✅ Use retry decorators for API calls
- ✅ Log important operations
- ✅ Validate inputs with Pydantic
- ✅ Handle errors gracefully
- ✅ Use configuration via environment variables
- ✅ Commit with descriptive messages
- ✅ Test thoroughly before pushing

### DON'T:

- ❌ Commit API keys or secrets
- ❌ Push test output or generated files
- ❌ Make breaking changes without discussion
- ❌ Skip error handling
- ❌ Use bare `except:` clauses
- ❌ Hard-code configuration values
- ❌ Write code without docstrings
- ❌ Skip testing edge cases
- ❌ Leave debugging print statements
- ❌ Commit code that doesn't run

### Security Considerations:

- Always use `.env` for sensitive data
- Never commit `.env` file
- Use `.env.example` for templates
- Validate all user inputs
- Use parameterized SQL queries
- Be careful with web scraping (rate limits)
- Handle API keys securely

---

## Quick Reference

### Environment Variables Cheat Sheet

```bash
# Required
OPENAI_API_KEY=sk-...
NCBI_EMAIL=your@email.com

# Model Selection
TIER2_MODEL=gpt-4o-mini      # Cheap filter
TIER3_MODEL=gpt-4o           # Smart extractor

# Paper Sources (PubMind-first recommended)
USE_PUBMIND=true
PUBMIND_ONLY=true            # Ignore PubMed/EuropePMC
USE_PUBMED=false
USE_EUROPEPMC=false

# Tier Configuration
ENABLE_TIER1=true            # Keyword filter
ENABLE_TIER2=true            # LLM classification
ENABLE_TIER3=true            # Expert extraction
TIER1_MIN_KEYWORDS=2
TIER2_CONFIDENCE_THRESHOLD=0.5
```

### CLI Commands Cheat Sheet

```bash
# Complete workflow (recommended)
python automated_workflow.py BRCA1 --email your@email.com

# Custom limits
python automated_workflow.py BRCA1 --max-pmids 200 --max-downloads 100

# Query database
python query_variants_db.py output/BRCA1.db --stats

# Re-run extraction
python rerun_extraction.py output/BRCA1/{timestamp}/

# Run tests
pytest
pytest test_triage.py -v
pytest --cov=pipeline
```

### Common Import Patterns

```python
# Pipeline modules
from pipeline.sourcing import PaperSourcer
from pipeline.filters import KeywordFilter, InternFilter
from pipeline.extraction import ExpertExtractor
from pipeline.aggregation import DataAggregator

# Harvesting
from harvesting import PMCHarvester

# Utilities
from utils.llm_utils import BaseLLMCaller
from utils.pubmed_utils import query_pubmed_for_gene, fetch_paper_metadata
from utils.retry_utils import api_retry, llm_retry

# Configuration
from config.settings import get_settings
from config import config

# Models
from models import Paper, FilterResult, ExtractionResult, FilterDecision
```

---

## Conclusion

This guide provides comprehensive information for AI assistants working on GeneVariantFetcher. For specific implementation details, consult the source code files referenced throughout this document.

**Key Takeaways:**
1. **Tiered pipeline** reduces costs by 85% through progressive filtering
2. **PubMind-first** strategy prioritizes high-quality variant-focused literature
3. **Configuration-driven** design allows easy customization via environment variables
4. **Full-text + supplements** provides 10-100x more data than abstracts alone
5. **Modular architecture** makes components independently testable and reusable

**For Updates:**
- This documentation should be updated when major architectural changes occur
- Keep code examples synchronized with actual implementation
- Update version numbers and dependency information as project evolves

**Contact:**
- For questions or issues, open a GitHub issue
- For security concerns, contact maintainers directly

---

*Last Updated: 2025-12-03*
*Document Version: 1.0*
*Project Version: 0.1.0*


NOW ALSO ADDED:# CLAUDE.md - AI Assistant Guide for GenePhenExtract

**Last Updated:** 2025-12-02
**Project:** Gene Literature Collector
**Version:** 0.2.0

## Project Overview

GenePhenExtract is a Python-based pipeline for gathering structured publication metadata about genes of interest from PubMed. The tool retrieves matching articles, extracts useful metadata, evaluates whether articles likely report patient-level information, and provides downloadable URLs for literature access.

### Key Features

- Automatic synonym discovery from NCBI Gene database
- LLM-based relevance filtering using Claude AI (Anthropic)
- LLM-based synonym relevance checking
- Patient-level evidence detection using keyword heuristics
- Multiple output formats (JSON, CSV, SQLite, URLs)
- Automated file renaming and organization for downloaded PDFs
- PMC XML availability detection

### Core Use Case

Researchers studying specific genes (e.g., BRCA1, SCN5A, TP53) need to collect relevant literature. This tool automates PubMed searches, filters irrelevant papers (e.g., distinguishing "TTR" as Transthyretin gene vs "time to reimbursement"), and organizes the results for further analysis.

---

## Repository Structure

```
GenePhenExtract/
├── src/gene_literature/         # Main package source
│   ├── __init__.py              # Package exports
│   ├── collector.py             # High-level orchestration and query building
│   ├── pubmed_client.py         # PubMed E-utilities API client
│   ├── synonym_finder.py        # NCBI Gene database synonym discovery
│   ├── synonym_relevance_checker.py  # LLM-based synonym assessment
│   ├── relevance_checker.py     # LLM-based paper relevance filtering
│   └── writer.py                # Output formatting (JSON, CSV, SQLite, URLs)
├── tests/                       # Unit tests
│   └── test_collector.py        # Main test suite
├── collect_literature.py        # CLI entry point for literature collection
├── rename_downloads.py          # CLI tool for organizing downloaded files
├── pyproject.toml               # Package configuration and dependencies
├── README.md                    # User-facing documentation
└── CLAUDE.md                    # This file - AI assistant guide

```

### Key Modules

#### `collector.py`
- **Purpose:** High-level orchestration of literature collection workflow
- **Key Functions:**
  - `build_gene_query()`: Constructs PubMed search queries from gene + synonyms
  - `LiteratureCollector.collect()`: Main entry point - searches PubMed, fetches metadata, filters results
- **Dependencies:** `PubMedClient`, optional `RelevanceChecker`

#### `pubmed_client.py`
- **Purpose:** Direct interaction with NCBI PubMed E-utilities API
- **Key Classes:**
  - `PubMedClient`: HTTP client with retry logic, rate limiting
  - `ArticleMetadata`: Dataclass holding all article metadata fields
- **Rate Limiting:** Respects NCBI guidelines (~3 requests/second)
- **Batching:** Processes PMIDs in batches of 200 to avoid URI length limits (HTTP 414)
- **Retry Logic:** Exponential backoff on failures (up to 3 retries)

#### `synonym_finder.py`
- **Purpose:** Automatic gene synonym discovery from NCBI Gene database
- **Key Functions:**
  - `SynonymFinder.find_gene_synonyms()`: Queries NCBI Gene for aliases
  - `interactive_synonym_selection()`: CLI prompt for synonym selection
- **Integration:** Can optionally use `SynonymRelevanceChecker` for LLM-based filtering

#### `relevance_checker.py`
- **Purpose:** LLM-based filtering of irrelevant papers
- **Model Used:** `claude-3-5-haiku-20241022` (cost-effective)
- **Key Pattern:** Uses Anthropic API to assess if papers are genuinely about the gene vs. abbreviation collisions
- **Fallback:** Gracefully degrades if API key missing or anthropic package unavailable

#### `writer.py`
- **Purpose:** Export collected metadata in multiple formats
- **Formats Supported:**
  - `json`: Pretty-printed JSON with all fields
  - `csv`: Flat CSV with standard fields
  - `sqlite`: SQLite database with `articles` table
  - `urls`: Formatted text file with downloadable URLs

---

## Development Workflows

### Setting Up the Environment

```bash
# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # or .venv\Scripts\activate on Windows

# Install package in editable mode
pip install -e .

# Install with optional dependencies
pip install -e ".[relevance]"  # Adds anthropic for LLM features
pip install -e ".[test]"       # Adds pytest for testing
pip install -e ".[dev]"        # Includes test dependencies
```

### Running Tests

```bash
# Run all tests
pytest

# Run with verbose output
pytest -v

# Run specific test file
pytest tests/test_collector.py

# Run with coverage (if installed)
pytest --cov=gene_literature
```

**Test Patterns:**
- Uses `DummyPubMedClient` for mocking PubMed API calls
- Provides sample XML payloads (`SAMPLE_XML`, `SAMPLE_XML_PMC_FORMAT`)
- Tests both old (`IdType="pmcid"`) and new (`IdType="pmc"`) PubMed formats
- Uses `tmp_path` fixture for file I/O tests

### Git Workflow

**Branch Naming Convention:**
- Feature branches: `claude/feature-description-<session-id>`
- Main branch: `main` (or default branch - check with `git branch -a`)

**Commit Message Patterns:**
Based on recent commits, follow these conventions:
- Use imperative mood: "Add feature" not "Added feature"
- Be descriptive but concise
- Reference purpose/context when helpful

Examples from history:
- "Add LLM-based relevance filtering to remove irrelevant papers"
- "Fix PMCID and XML detection to support multiple PubMed formats"
- "Fix HTTP 414 error by batching PubMed metadata requests"

**Never:**
- Push to `main` directly
- Force push to shared branches
- Commit sensitive API keys or credentials

---

## Code Conventions and Patterns

### Python Style

**Type Hints:**
- Use throughout the codebase
- Import from `__future__ import annotations` for forward references
- Optional types clearly marked with `Optional[T]`

**Example:**
```python
from __future__ import annotations

from typing import List, Optional, Sequence

def build_gene_query(gene: str, synonyms: Optional[Sequence[str]] = None) -> str:
    """Build a simple PubMed query using the provided gene and synonyms."""
    ...
```

**Logging:**
- Use module-level logger: `logger = logging.getLogger(__name__)`
- Log levels:
  - `DEBUG`: Low-level details (query construction, XML parsing)
  - `INFO`: Major workflow steps (searching PubMed, fetching metadata)
  - `WARNING`: Degraded functionality (missing API key, fallback behavior)
  - `ERROR`: Recoverable errors (failed API calls with retry)
  - `CRITICAL`: Not used in this codebase

**Example:**
```python
logger.info("Searching PubMed with query: %s", query)
logger.debug("Constructed PubMed query: %s", query)
logger.warning("Relevance filtering requested but no RelevanceChecker provided")
```

**Dataclasses:**
- Use `@dataclass` for data-holding classes
- Include `to_dict()` method for serialization
- Example: `ArticleMetadata`, `RelevanceScore`, `GeneSynonym`

**Error Handling:**
- Custom exceptions inherit from appropriate base (e.g., `PubMedError(RuntimeError)`)
- Graceful degradation for optional features (LLM checking)
- Retry logic with exponential backoff for network calls

**Docstrings:**
- Use triple-quoted strings for all public functions/classes
- Google-style format preferred
- Include Args, Returns sections when applicable

### API Integration Patterns

**PubMed E-utilities:**
- Base URL: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils`
- User-Agent header required: `GeneLiteratureCollector/0.1 (+https://github.com/openai)`
- Rate limiting: ~3 requests/second (0.34s delay between batches)
- Authentication: Optional API key increases rate limits

**Anthropic Claude:**
- Model: `claude-3-5-haiku-20241022` for cost-effectiveness
- API key from environment variable `ANTHROPIC_API_KEY` or parameter
- Short max_tokens (200) for structured responses
- JSON response parsing with regex fallback

### Optional Dependencies Pattern

The project uses optional dependency groups defined in `pyproject.toml`:

```toml
[project.optional-dependencies]
relevance = [
    "anthropic>=0.39.0",
]
test = [
    "pytest>=7.0.0",
]
```

**Pattern for Optional Imports:**
```python
try:
    from .relevance_checker import RelevanceChecker
    RELEVANCE_CHECKER_AVAILABLE = True
except ImportError:
    RELEVANCE_CHECKER_AVAILABLE = False
    logger.debug("RelevanceChecker not available (anthropic package not installed)")
```

This allows the package to work without `anthropic` installed, with graceful degradation.

---

## Common Tasks for AI Assistants

### 1. Adding a New Output Format

**Steps:**
1. Add format choice to `collect_literature.py` argument parser
2. Implement `_write_<format>()` function in `writer.py`
3. Add format handling in `write_metadata()` dispatcher
4. Update tests in `test_collector.py`
5. Document in README.md

**Key Files:**
- `src/gene_literature/writer.py`
- `collect_literature.py`
- `tests/test_collector.py`
- `README.md`

### 2. Adding New Metadata Fields

**Steps:**
1. Add field to `ArticleMetadata` dataclass in `pubmed_client.py`
2. Extract field in `PubMedClient.fetch_metadata()` XML parsing
3. Update CSV fieldnames in `writer.py:_write_csv()`
4. Update SQLite schema in `writer.py:_write_sqlite()`
5. Update test assertions in `test_collector.py`
6. Document in README.md

**Key Patterns:**
- Use `Optional[T]` for fields that may be missing
- Add to end of dataclass to maintain backward compatibility
- Update `to_dict()` automatically via `@dataclass` + `asdict()`

### 3. Modifying PubMed Query Logic

**Location:** `src/gene_literature/collector.py:build_gene_query()`

**Current Pattern:**
```python
terms = [gene, *(synonyms or [])]
quoted = [f'"{term}"[Title/Abstract]' for term in sanitized]
query = " OR ".join(quoted)
```

**Considerations:**
- PubMed search field tags: `[Title/Abstract]`, `[MeSH Terms]`, etc.
- Boolean operators: `OR`, `AND`, `NOT`
- Quote terms for exact matching
- Test with actual PubMed API to verify syntax

### 4. Improving Relevance Filtering

**Location:** `src/gene_literature/relevance_checker.py:check_relevance()`

**Current Approach:**
- Sends title + abstract to Claude Haiku
- Uses structured prompt with JSON response format
- Parses JSON with regex fallback

**Improvement Ideas:**
- Batch processing (currently sequential)
- Caching results to avoid re-checking same papers
- Fine-tuning prompt based on specific gene types
- Adding few-shot examples to prompt

**Cost Optimization:**
- Current: ~$0.001-0.003 per paper
- Uses Haiku model (cheapest)
- max_tokens=200 keeps responses short
- Consider batching to reduce API overhead

### 5. Adding New CLI Options

**Steps:**
1. Add argument in `collect_literature.py:parse_args()`
2. Pass argument through to relevant function
3. Update help text
4. Document in README.md examples section

**Naming Convention:**
- Use `--kebab-case` for CLI flags
- Use `snake_case` for Python variable names
- Boolean flags: use `action="store_true"`

### 6. Handling XML Format Changes

**Context:** PubMed periodically updates XML schema

**Key Patterns:**
- Check both uppercase and lowercase attribute names: `IdType` vs `idtype`
- Support multiple identifier formats: `"pmc"` and `"pmcid"`
- Use robust XPath selectors: `.//ArticleIdList/ArticleId`
- Add defensive `None` checks

**Example from codebase:**
```python
id_type = article_id.get("IdType") or article_id.get("idtype") or ""
id_type_lower = id_type.lower()
if id_type_lower in ("pmc", "pmcid"):
    # Handle both formats
```

---

## Dependencies and Requirements

### Core Dependencies
- **Python:** >= 3.9
- **Standard Library Only:** No dependencies for basic functionality
  - `urllib` for HTTP requests
  - `xml.etree.ElementTree` for XML parsing
  - `sqlite3` for database export
  - `json`, `csv` for data serialization

### Optional Dependencies
- **anthropic** >= 0.39.0: Required for LLM-based relevance filtering and synonym checking
- **pytest** >= 7.0.0: Required for running tests

### External APIs
- **NCBI E-utilities:** PubMed search and metadata retrieval
  - Free tier: ~3 requests/second
  - With API key: ~10 requests/second
  - Email recommended for compliance
- **Anthropic Claude API:** LLM-based filtering
  - Requires API key
  - Uses Haiku model for cost optimization

---

## Testing Strategy

### Test Coverage

**Current Test File:** `tests/test_collector.py`

**Patterns:**
- Mock PubMed API with `DummyPubMedClient` class
- Use sample XML payloads for predictable parsing
- Test both success paths and error handling
- Use `pytest.mark.parametrize` for multiple test cases

**What to Test:**
- Query construction with various inputs
- XML parsing for different PubMed formats
- Output format generation
- Error handling and edge cases

**What NOT to Test:**
- Actual PubMed API calls (use mocks)
- Actual Claude API calls (use mocks or skip if API key missing)
- Network-dependent behavior

### Adding Tests

**Pattern:**
```python
def test_feature_name(tmp_path: Path):
    # Arrange
    client = DummyPubMedClient(["12345"], SAMPLE_XML)
    collector = LiteratureCollector(client)

    # Act
    results = collector.collect("GENE")

    # Assert
    assert len(results) == 1
    assert results[0].pmid == "12345"
```

---

## Important Implementation Details

### Rate Limiting and Batching

**Why:** NCBI requires respectful API usage

**Implementation:**
- 0.34 second delay between batches (≈3 requests/second)
- Batch size: 200 PMIDs per request (avoids HTTP 414 URI Too Long)
- Located in: `pubmed_client.py:fetch_metadata()`

**Pattern:**
```python
for batch_num in range(total_batches):
    # ... fetch batch ...
    if batch_num < total_batches - 1:
        time.sleep(0.34)  # Respect rate limits
```

### XML Parsing Robustness

**Challenge:** PubMed XML schema has evolved over time

**Solutions:**
- Check multiple attribute name variations (uppercase/lowercase)
- Accept multiple identifier types (`"pmc"` and `"pmcid"`)
- Gracefully handle missing fields with `Optional` types
- Use defensive `None` checks throughout

**Pattern:**
```python
def _find_text(root: ET.Element, selector: str) -> Optional[str]:
    element = root.find(selector)
    if element is None:
        return None
    text = element.text or ""
    text = text.strip()
    return text or None
```

### Patient-Level Evidence Detection

**Method:** Simple keyword heuristics

**Keywords:**
```python
PATIENT_KEYWORDS = {
    "patient", "patients", "case", "case report", "cases",
    "cohort", "subjects", "clinical"
}
```

**Limitation:** Basic pattern matching, may have false positives/negatives

**Future Improvement Opportunity:** Could use LLM-based detection similar to relevance checking

### Relevance Filtering Prompt Strategy

**Key Insight:** The prompt is tuned to err on the side of inclusion

**Pattern:**
- Explicitly lists what to EXCLUDE (non-genetic topics, wrong abbreviation usage)
- States "when uncertain, treat it as relevant"
- Asks for JSON response format for structured parsing
- Includes confidence score to allow threshold tuning

**Customization Point:** The prompt in `relevance_checker.py` can be tuned for specific research needs

---

## Environment Variables

- `ANTHROPIC_API_KEY`: Required for LLM-based relevance filtering and synonym checking
- `NCBI_API_KEY`: Optional, increases PubMed API rate limits (can also pass via `--api-key`)

---

## Debugging Tips

### Enable Debug Logging

```bash
python collect_literature.py GENE --log-level DEBUG
```

This shows:
- Exact PubMed query constructed
- Number of PMIDs returned
- XML parsing details
- Relevance check results for each paper

### Common Issues

**Issue:** HTTP 414 URI Too Long
**Solution:** Already fixed via batching in `pubmed_client.py`
**Commit:** `27d26d1 Fix HTTP 414 error by batching PubMed metadata requests`

**Issue:** PMCID not detected
**Solution:** Support both `IdType="pmc"` and `IdType="pmcid"`
**Commit:** `9c16760 Fix PMCID and XML detection to support multiple PubMed formats`

**Issue:** Relevance filtering not working
**Check:**
1. Is `anthropic` package installed? (`pip install anthropic`)
2. Is `ANTHROPIC_API_KEY` set?
3. Is `--filter-irrelevant` flag used?

**Issue:** Too many/few papers filtered
**Solution:** Adjust `--min-relevance-score` threshold (default: 0.7)

---

## Key Conventions Summary

### When Writing Code
- Use type hints consistently
- Import with `from __future__ import annotations`
- Add docstrings to public functions
- Use module-level loggers
- Handle optional dependencies gracefully
- Prefer dataclasses for data structures

### When Modifying Features
- Read existing code first to understand patterns
- Maintain backward compatibility for data formats
- Update all affected files (module, CLI, tests, docs)
- Test with actual PubMed API when changing query logic
- Consider rate limiting and API costs

### When Adding Tests
- Use `DummyPubMedClient` for mocking
- Use `tmp_path` for file I/O
- Test edge cases and error handling
- Don't rely on network/external APIs

### When Updating Documentation
- Update README.md for user-facing changes
- Update CLAUDE.md for architectural changes
- Include examples in CLI help text
- Document cost implications for LLM features

---

## Recent Development History

Based on recent commits:

1. **Relevance Filtering** (PR #22, #21): Adjusted LLM prompt for better filtering criteria
2. **LLM Synonym Checking** (PR #20): Added LLM-based relevance checking for auto-discovered synonyms
3. **Initial Relevance Filtering** (PR #19): Implemented LLM-based paper filtering
4. **Synonym Discovery** (PR #18): Added automatic NCBI Gene synonym discovery
5. **PMCID Detection Fix** (PR #17): Fixed XML parsing for multiple PubMed formats
6. **Batching Fix** (PR #16): Fixed HTTP 414 errors via batching
7. **URL Extraction** (PR #15): Added URL output format and file renaming tool

---

## Contact and Resources

- **Repository:** https://github.com/kroncke-lab/GenePhenExtract (assumed based on branch structure)
- **PubMed E-utilities Documentation:** https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **Anthropic API Documentation:** https://docs.anthropic.com/
- **Package Name:** gene-literature-collector

---

## AI Assistant Quick Reference

### When asked to add a feature:
1. Read relevant module files first
2. Understand existing patterns
3. Maintain consistency with current code style
4. Update tests, docs, and CLI together
5. Consider backward compatibility

### When asked to fix a bug:
1. Check if similar issues were fixed before (see git log)
2. Add test case that reproduces the bug
3. Fix the bug
4. Verify test passes
5. Check for similar patterns elsewhere in codebase

### When asked about the project:
1. Refer to README.md for user-facing information
2. Refer to this CLAUDE.md for architecture and conventions
3. Check module docstrings for implementation details
4. Look at tests for usage examples

---

**End of CLAUDE.md**
