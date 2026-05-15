# GeneVariantFetcher (GVF)

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-2.1.0-orange.svg)](CHANGELOG.md)

**Automated extraction of human genetic variant carriers from biomedical literature into a normalized SQLite database.**

---

## Overview

GeneVariantFetcher (GVF) is an end-to-end pipeline that discovers, downloads, and extracts genetic variant data from published literature. Given a gene name, it queries multiple databases (PubMed, PubMind, Europe PMC), downloads full-text articles and supplemental materials, uses LLM-powered extraction to identify variant carriers, and outputs a normalized SQLite database ready for downstream analysis. Designed for genetics researchers studying penetrance, variant pathogenicity, and genotype-phenotype correlations, GVF automates what would otherwise take weeks of manual curation — with the key insight that **70-80% of variant data lives in supplemental files** (Excel tables, PDFs, Word documents), not just main text.

---

## Key Features

- 🔍 **Multi-source literature discovery** — PubMind, PubMed, and Europe PMC integration
- 📄 **Full-text + supplement harvesting** — PMC, Elsevier, Springer, Wiley APIs with automatic format conversion
- 🧠 **Intelligent Data Scout** — Condenses papers to high-value "data zones," reducing LLM token usage by 60-80%
- 🤖 **Tiered LLM extraction** — Keyword filter → LLM triage → Expert extraction with model cascade
- 📊 **Rich table extraction** — Excel, Word, and PDF tables parsed with comprehensive header recognition
- 🔄 **Fuzzy variant matching** — Normalizes all frameshift/nonsense/deletion conventions for deduplication
- 🔤 **Unicode normalization** — Handles arrow variants (→, ➔, ⟶), concatenated gene+variant patterns (HERGG604S)
- 💾 **SQLite output** — Normalized relational database with penetrance, phenotypes, and functional data
- ⏸️ **Checkpoint/resume** — Jobs survive interruption and resume from last completed step
- 🖥️ **GUI + CLI** — Web interface for interactive use, CLI for automation and scripting

---

## Performance

Current measured KCNH2 performance against the normalized 2026 gold input
(`gene_variant_fetcher_gold_standard/normalized/KCNH2_recall_input.csv`) and the
manual-recovery v12 database scored on 2026-05-15:

| Metric | Recall |
|--------|--------|
| **PMIDs** | 184/262 (70.2%) |
| **Variant rows** | 542/991 (54.7%) |
| **Unique variants** | 323/530 (60.9%) |
| **Patients/carriers** | 1758/2674 (65.7%) |
| **Affected individuals** | 1095/1635 (67.0%) |

The older 59.1% number from the 2025-12-11 KCNH2 baseline is historical. The
current branch adds multi-gene recall inputs, authenticated paywall recovery,
clinical mutation-list table handling, and matcher improvements; source access
is now the main limiter for KCNH2 recall.

> **Note:** Performance varies by gene and literature availability. Cardiac genes (8 validated) have optimized configurations with known protein lengths; other genes work generically but may have lower recall.

---

## Quick Start

```bash
# 1. Clone and install
git clone https://github.com/your-org/GeneVariantFetcher.git
cd GeneVariantFetcher
python -m venv .venv && source .venv/bin/activate
pip install -e .

# 2. Configure API keys (create .env file)
cat > .env << 'EOF'
OPENAI_API_KEY=sk-...
NCBI_EMAIL=your@email.com
EOF

# 3. Run your first extraction
gvf extract BRCA1 --email your@email.com --output ./results

# Or launch the GUI
python main.py
```

Your results will be in `./results/BRCA1/{timestamp}/BRCA1.db`

---

## Installation

### Requirements

- **Python 3.11+** (tested on 3.11.8)
- **poppler-utils** (for PDF processing)

### System Dependencies

```bash
# Ubuntu/Debian
sudo apt-get install poppler-utils

# macOS
brew install poppler

# RHEL/CentOS
sudo dnf install poppler-utils
```

### Python Setup

```bash
# Create virtual environment
python3.11 -m venv .venv
source .venv/bin/activate  # Linux/macOS
# or: .venv\Scripts\activate  # Windows

# Install GVF
pip install -e .

# For GUI support
pip install -r gui/requirements.txt

# For browser automation (optional, for paywalled papers)
pip install playwright && playwright install chromium
```

### Verify Installation

```bash
gvf --help
# Should show: extract, scout, audit-paywalls commands
```

---

## Usage

### CLI Commands

#### Basic Extraction

```bash
gvf extract GENE --email EMAIL --output DIR
```

**Required arguments:**
- `GENE` — Gene symbol (e.g., KCNH2, BRCA1, SCN5A)
- `--email` — Your email for NCBI API access
- `--output` — Directory for results

**Examples:**

```bash
# Basic run
gvf extract TTN --email researcher@university.edu --output ./results

# With automatic synonym discovery
gvf extract KCNH2 --email you@email.com --output ./results --auto-synonyms

# Limit scope for quick testing
gvf extract SCN5A --email you@email.com --output ./results \
    --max-pmids 50 --max-downloads 25

# Use clinical triage filter (stricter relevance)
gvf extract RYR2 --email you@email.com --output ./results --clinical-triage

# Skip to stronger model immediately (slower, more thorough)
gvf extract CACNA1C --email you@email.com --output ./results --tier-threshold 0

# Run Data Scout first for better token efficiency
gvf extract KCNQ1 --email you@email.com --output ./results --scout-first
```

#### Key Flags

| Flag | Description | Default |
|------|-------------|---------|
| `--max-pmids N` | Maximum PMIDs to discover | 100 |
| `--max-downloads N` | Maximum papers to download | 50 |
| `--tier-threshold N` | Model cascade threshold (0=always use strong model) | 1 |
| `--clinical-triage` | Use stricter clinical relevance filter | off |
| `--auto-synonyms` | Auto-discover gene aliases from NCBI | off |
| `--synonym SYN` | Add manual synonym (repeatable) | — |
| `--scout-first` | Run Data Scout before extraction | off |
| `--verbose` | Enable verbose logging | off |

#### Data Scout (Standalone)

Run the Data Scout on existing downloaded papers:

```bash
gvf scout ./results/KCNH2/20260210_143000/pmc_fulltext --gene KCNH2
```

### GUI Mode

```bash
python main.py
# Opens http://localhost:8000 in your browser
```

The GUI provides:
- Real-time progress tracking
- Job management (pause, resume, cancel)
- Settings configuration
- Folder job support (process existing paper collections)

---

## Pipeline Architecture

```
INPUT: Gene Symbol (e.g., "KCNH2")
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  STAGE 1: SYNONYM DISCOVERY (optional)                              │
│    • Query NCBI Gene database for aliases                           │
│    • Output: ["HERG", "LQT2", "Kv11.1", ...]                        │
└─────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  STAGE 2: PMID DISCOVERY                                            │
│    • PubMind API (gene-variant focused)                             │
│    • PubMed E-Utilities (broad search)                              │
│    • Europe PMC (optional, additional coverage)                     │
│    • Output: {gene}_pmids.txt (merged, deduplicated)               │
└─────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  STAGE 3: ABSTRACT FETCH                                            │
│    • PubMed E-Utilities efetch                                      │
│    • Output: abstract_json/{PMID}.json                             │
└─────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  STAGE 4: PAPER FILTERING                                           │
│    • Tier 1: KeywordFilter (regex, ~65 clinical terms)             │
│    • Tier 2: InternFilter (LLM relevance, gpt-4o-mini)             │
│    • Output: pmid_status/filtered_out.csv                          │
└─────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  STAGE 5: FULL-TEXT + SUPPLEMENT DOWNLOAD                            │
│    Full-text sources (priority order):                              │
│      1. PMC Open Access (BioC XML)                                  │
│      2. Elsevier API (ScienceDirect)                               │
│      3. Springer API (SpringerLink)                                │
│      4. Wiley API                                                   │
│      5. Unpaywall (OA PDF links)                                   │
│    Supplement sources (UnifiedSupplementFetcher):                   │
│      Tier 1: Europe PMC + NCBI OA (free)                           │
│      Tier 2: Elsevier supplement API                               │
│      Tier 3: DOI-based web scraping (publisher-specific)           │
│    Conversion: XML/PDF/Excel/Word → Markdown                        │
│    Output: pmc_fulltext/{PMID}_FULL_CONTEXT.md                     │
└─────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  STAGE 6: DATA SCOUT                                                │
│    • Identifies tables, variant lists, patient data                │
│    • Creates condensed files for efficient LLM processing          │
│    • Output: pmc_fulltext/{PMID}_DATA_ZONES.md                     │
└─────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  STAGE 7: VARIANT EXTRACTION (Tier 3)                               │
│    • Input: DATA_ZONES.md > FULL_CONTEXT.md > abstract             │
│    • Model cascade: gpt-4o-mini → gpt-4o (if low yield)            │
│    • Pre-scan: Regex scanner on FULL_CONTEXT.md (not condensed)    │
│      - Supports concatenated patterns (e.g., HERGG604S, KCNH2A561V)│
│      - Unicode arrow normalization (→, ➔, ⟶ → standard arrow)     │
│    • Output: extractions/{gene}_PMID_{pmid}.json                   │
└─────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  STAGE 8: AGGREGATION & SQLITE MIGRATION                            │
│    • Normalize variant names (HGVS standardization)                │
│    • Fuzzy matching for deduplication                              │
│    • Aggregate penetrance statistics                               │
│    • Output: {gene}.db (SQLite database)                           │
└─────────────────────────────────────────────────────────────────────┘
         │
         ▼
OUTPUT: SQLite database + JSON summaries + workflow logs
```

---

## API Keys

### Required

| API | Purpose | How to Get |
|-----|---------|-----------|
| **OpenAI** | LLM extraction (Tier 2 & 3) | [platform.openai.com](https://platform.openai.com/api-keys) |
| **NCBI Email** | PubMed access (required for E-utilities) | Your institutional email |

### Optional (Improves Coverage)

| API | Purpose | How to Get |
|-----|---------|-----------|
| **NCBI API Key** | 10x higher rate limits (3→10 req/sec) | [ncbi.nlm.nih.gov/account](https://www.ncbi.nlm.nih.gov/account/) |
| **Elsevier** | ScienceDirect full-text | [dev.elsevier.com](https://dev.elsevier.com/) |
| **Springer** | SpringerLink content | [dev.springernature.com](https://dev.springernature.com/) |
| **Wiley** | Wiley Online Library | [onlinelibrary.wiley.com/library-info/resources/text-and-datamining](https://onlinelibrary.wiley.com/library-info/resources/text-and-datamining) |

### Free (No Key Required)

- **PMC Open Access** — Free BioC XML for all open access papers
- **Unpaywall** — Uses your NCBI email for OA link resolution

---

## Output Format

### Directory Structure

```
{output}/{GENE}/{TIMESTAMP}/
├── {GENE}.db                           # ← Primary output: SQLite database
├── {GENE}_pmids.txt                    # Discovered PMIDs
├── {GENE}_workflow_summary.json        # Run statistics
├── {GENE}_penetrance_summary.json      # Aggregated variant data
├── {GENE}_workflow.log                 # Execution log
├── run_manifest.json                   # Full execution metadata
│
├── abstract_json/                      # Paper metadata
│   └── {PMID}.json
│
├── pmid_status/                        # Filter decisions
│   ├── filtered_out.csv
│   └── extraction_failures.csv
│
├── pmc_fulltext/                       # Downloaded papers
│   ├── {PMID}_FULL_CONTEXT.md          # Full paper (main + supplements)
│   ├── {PMID}_DATA_ZONES.md            # Condensed high-value sections
│   ├── {PMID}_supplements/             # Original supplement files
│   ├── successful_downloads.csv
│   └── paywalled_missing.csv           # Papers needing manual fetch
│
└── extractions/                        # Per-paper extractions
    └── {GENE}_PMID_{PMID}.json
```

### SQLite Schema

| Table | Description | Key Columns |
|-------|-------------|-------------|
| `papers` | Paper metadata | pmid (PK), title, journal, doi, pmc_id |
| `variants` | Unique variants | variant_id (PK), cdna_notation, protein_notation, clinical_significance |
| `individual_records` | Per-patient data | record_id (PK), variant_id (FK), age_at_onset, sex, affected_status |
| `penetrance_data` | Carrier statistics | total_carriers, affected_count, unaffected_count, penetrance_percentage |
| `age_dependent_penetrance` | Age-stratified penetrance | age_range, penetrance_percentage, carriers_in_range |
| `functional_data` | Functional studies | summary, assays (JSON) |
| `phenotypes` | Clinical phenotypes | patient_count, phenotype_description |
| `variant_papers` | Variant↔paper links | source_location, key_quotes (JSON) |

See [docs/SQLITE_MIGRATION_GUIDE.md](docs/SQLITE_MIGRATION_GUIDE.md) for full schema details.

---

## Configuration

### Environment Variables (.env)

Create a `.env` file in the project root:

```bash
# Required
OPENAI_API_KEY=sk-...              # For LLM extraction
NCBI_EMAIL=your@email.com          # For PubMed API

# Recommended
NCBI_API_KEY=...                   # 10x rate limits

# Optional (publisher access)
ELSEVIER_API_KEY=...               # ScienceDirect
SPRINGER_API_KEY=...               # SpringerLink
WILEY_API_KEY=...                  # Wiley Online

# Alternative LLM providers
ANTHROPIC_API_KEY=...              # For Claude
GEMINI_API_KEY=...                 # For Gemini
```

### Pipeline Settings (config/settings.py)

Key tunable parameters:

```python
# Filtering
TIER1_MIN_KEYWORD_MATCHES = 2      # Keywords required to pass Tier 1
TIER2_CONFIDENCE_THRESHOLD = 0.5   # LLM confidence for Tier 2

# Extraction
TIER3_MODELS = ["gpt-4o-mini", "gpt-4o"]  # Model cascade
TIER3_MAX_TOKENS = 16000           # Max tokens for extraction
TIER3_TEMPERATURE = 0.1            # LLM temperature

# Downloads
DOWNLOAD_DELAY = 2.0               # Seconds between API requests
MAX_RETRIES = 3                    # Retry attempts for failed downloads
```

---

## Supported Genes

### Validated (Cardiac Ion Channels)

These 8 genes have known protein lengths in `utils/variant_normalizer.py` for position validation:

| Gene | Associated Condition | Protein Length |
|------|---------------------|---------------|
| KCNH2 | Long QT Syndrome Type 2 | 1159 aa |
| KCNQ1 | Long QT Syndrome Type 1 | 676 aa |
| SCN5A | Brugada Syndrome, Long QT Type 3 | 2016 aa |
| KCNE1 | Long QT Syndrome Type 5 | 129 aa |
| KCNE2 | Long QT Syndrome Type 6 | 123 aa |
| KCNJ2 | Andersen-Tawil Syndrome | 427 aa |
| CACNA1C | Timothy Syndrome | 2221 aa |
| RYR2 | CPVT | 4967 aa |

### Generic Support

Any gene symbol works with the pipeline. However:
- Synonym discovery depends on NCBI Gene database coverage
- Variant normalization uses generic HGVS rules (no gene-specific aliases)
- No pre-tuned keyword filters for non-cardiac domains

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "No papers downloaded" | Check `paywalled_missing.csv`; add publisher API keys or use `scripts/fetch_paywalled.py` |
| "Extraction failed" | Check `extraction_failures.csv`; increase `TIER3_MAX_TOKENS` |
| "All papers filtered out" | Lower `TIER2_CONFIDENCE_THRESHOLD` or use `--tier-threshold 0` |
| Rate limiting errors | Set `NCBI_API_KEY` for higher limits; increase `DOWNLOAD_DELAY` |
| GUI won't start | Install: `pip install -r gui/requirements.txt` |
| PDF conversion fails | Install poppler: `apt install poppler-utils` |

### Fetching Paywalled Papers

For papers that couldn't be auto-downloaded:

```bash
# Canonical recovery path: authenticated browser HTML plus quality gate plus PMC fallback
python scripts/fetch_paywalled.py \
  --input results/GENE/TIMESTAMP/pmc_fulltext/paywalled_missing.csv \
  --output-dir results/GENE/TIMESTAMP/pmc_fulltext

# Manual fallback for papers that remain hard-blocked
python -m cli.fetch_manager results/GENE/TIMESTAMP/pmc_fulltext/paywalled_missing.csv
```

Validated recovery artifacts should be copied into the canonical run's
`pmc_fulltext/` directory and then re-extracted before running recall scoring.

---

## Contributing

We welcome contributions! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Adding Support for New Genes

1. **Add protein length** in `utils/variant_normalizer.py`:
   ```python
   PROTEIN_LENGTHS = {
       ...
       "NEW_GENE": 1234,  # Protein length in amino acids
   }
   ```

2. **Add variant aliases** (optional) in `utils/variant_normalizer.py`

3. **Create golden test set** in `golden_test_set/` for validation

### Extending Extractors

- **New publisher API**: Add module in `harvesting/`, register in `orchestrator.py`
- **New extraction field**: Modify prompts in `pipeline/prompts.py`, update schema

### Development Setup

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/

# Run live network/institutional tests explicitly
GVF_TEST_OUTPUT_DIR=/tmp/gvf_tests pytest -m requires_network tests/integration

# Run with coverage
pytest --cov=. tests/

# Format code
ruff format .
```

---

## Citation

If you use GeneVariantFetcher in your research, please cite:

```bibtex
@software{genevariantfetcher,
  title = {GeneVariantFetcher: Automated Extraction of Genetic Variant Carriers from Biomedical Literature},
  author = {Kronck, Brett M. and Roden, Dan M.},
  year = {2026},
  version = {2.1.0},
  url = {https://github.com/your-org/GeneVariantFetcher}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Acknowledgments

- Built at Vanderbilt University Medical Center
- Supported by NIH/NHLBI grants [grant numbers]
- Uses [LiteLLM](https://github.com/BerriAI/litellm) for LLM provider abstraction
- PDF processing via [pdfplumber](https://github.com/jsvine/pdfplumber) and [PyMuPDF](https://github.com/pymupdf/PyMuPDF)
