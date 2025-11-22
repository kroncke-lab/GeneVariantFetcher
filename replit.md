# Overview

This is a collection of Python command-line tools for extracting genetic variant data, patient-level information, and full-text articles from PubMed literature. The core functionality revolves around harvesting full-text papers from PubMed Central and extracting structured biomedical data using AI-powered analysis.

The tools provide comprehensive functionality for researchers to:
- Download full-text articles with supplemental materials from PubMed Central
- Extract individual-level variant and phenotype data from literature
- Process and triage clinical data
- Generate structured outputs for analysis

All tools are designed to run locally on your computer without requiring a web browser.

# User Preferences

Preferred communication style: Simple, everyday language.

# System Architecture

## Application Framework
- **Technology**: Python 3.11+ with command-line scripts
- **Rationale**: Simple, straightforward command-line tools that can be run locally without web dependencies
- **Design Pattern**: Modular scripts with clear separation of concerns

## Core Components

### 1. PMC Full-Text Harvester (`harvest_pmc_fulltext.py`)
**Purpose**: Download complete PubMed Central articles with all supplemental materials.

**APIs Used**:
- **NCBI E-utilities**: PMID → PMCID conversion via `Bio.Entrez.elink()`
- **EuropePMC fullTextXML**: `https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML`
- **EuropePMC supplementaryFiles**: `https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles`

**Key Features**:
1. PMID → PMCID Conversion
2. Full-Text Download (complete article XML)
3. Supplemental Files Download (Excel, Word, PDF, ZIP)
4. Unified Markdown Output (everything in one `{PMID}_FULL_CONTEXT.md` file)
5. Excel Table Conversion (to markdown tables for LLM processing)
6. Error Handling (logs paywalled/unavailable articles)

**Output Structure**:
- `{PMID}_FULL_CONTEXT.md` - Unified markdown with full-text + supplements
- `{PMID}_supplements/` - Downloaded supplemental files
- `successful_downloads.csv` - Log of successful harvests
- `paywalled_missing.csv` - Log of unavailable articles

### 2. Extraction Pipeline (`pipeline.py`, `extractor.py`)
**Purpose**: Extract structured individual-level variant and phenotype data from articles using AI.

**Features**:
- Individual-level data extraction (patients/subjects with variant data)
- Genetic variants in HGVS notation (cDNA and protein)
- Phenotypes with HPO codes when possible
- Clinical status (affected/unaffected)
- Demographics (age, sex, ancestry)
- Evidence sentences supporting each extraction

**APIs Used**:
- **OpenAI GPT-5**: Biomedical text extraction
- **NCBI E-utilities**: PubMed abstract fetching

### 3. Clinical Data Triage (`clinical_data_triage.py`)
**Purpose**: Process and triage clinical data from various sources.

**Features**:
- Variant classification and filtering
- Quality control and validation
- Structured output generation
- Integration with extraction pipeline

## Data Processing
- **Format**: JSON, CSV, and markdown outputs
- **Data Transformation**: Pandas DataFrames for tabular processing
- **Export Formats**: CSV, TXT, and JSON
- **Error Handling**: Comprehensive logging and error tracking

## Error Handling Approach
- Network request exceptions
- JSON/XML parsing errors
- HTTP status errors via `raise_for_status()`
- Timeout protection
- Rate limit detection and automatic retry

# External Dependencies

## Third-Party APIs
- **NCBI E-utilities** (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/)
  - Used to convert PMIDs to PMCIDs and fetch PubMed abstracts
  - No authentication required (public API)
  - Returns XML format

- **EuropePMC API** (https://www.ebi.ac.uk/europepmc/)
  - Used to fetch full-text articles and supplemental files
  - No authentication required (public API)
  - Returns XML format with comprehensive article data

- **OpenAI GPT-5** (via OpenAI API)
  - Used for biomedical text extraction from PubMed articles
  - Requires API key (set via environment variables)
  - Analyzes article text to extract individual-level variant and phenotype data

## Python Libraries
- **biopython**: NCBI Entrez API access
- **requests**: HTTP client for API calls
- **pandas**: Data manipulation and tabular processing
- **markitdown**: Markdown conversion for supplemental files
- **openai**: OpenAI Python SDK for LLM calls
- **tenacity**: Retry logic with exponential backoff
- **python-docx**: Word document processing
- **openpyxl**: Excel file processing
- **tabulate**: Table formatting

## Data Sources
- **PubMed**: Primary source for literature and abstracts
- **PubMed Central (PMC)**: Full-text articles and supplemental materials
- **dbSNP**: Referenced through rsID for variant information
- **HPO**: Human Phenotype Ontology for standardized phenotype coding

# Recent Changes (November 22, 2025)

## Removal of Streamlit Web Interface
**Motivation**: User requested removal of browser-based interface to run everything locally via command-line.

**Changes**:
1. Removed `app.py` (Streamlit web application)
2. Removed Streamlit and Plotly dependencies from `pyproject.toml`
3. Removed `.streamlit/` configuration directory
4. Updated `.replit` to remove web-based run commands
5. Rewrote `README.md` to focus on command-line usage
6. Updated documentation to emphasize local execution

**Impact**:
- Repository now focuses entirely on command-line tools
- No web server or browser required
- Simpler dependency tree
- All functionality available via Python scripts

## Standalone PMC Full-Text Harvester (November 20, 2025)

**Purpose**: Created standalone Python script to harvest full-text papers and ALL supplemental materials from PubMed Central for comprehensive LLM-based data extraction.

**Script**: `harvest_pmc_fulltext.py`

**Key Features**:
1. **PMID → PMCID Conversion**: Uses NCBI E-utilities to identify PubMed Central IDs
2. **Full-Text Download**: Retrieves complete article XML from EuropePMC API
3. **Supplemental Files**: Downloads ALL attachments (Excel, Word, PDF, ZIP)
4. **Unified Markdown**: Converts everything to a single `{PMID}_FULL_CONTEXT.md` file per article
5. **Excel Table Conversion**: Converts spreadsheets to markdown tables for LLM processing
6. **Error Handling**: Logs paywalled/unavailable articles to `paywalled_missing.csv`

**Usage**:
```python
from harvest_pmc_fulltext import PMCHarvester

harvester = PMCHarvester(output_dir="pmc_harvest")
pmids = ['35443093', '33442691', '34931732']
harvester.harvest(pmids, delay=2.0)
```

# Available Tools

## 1. harvest_pmc_fulltext.py
Download full-text articles and supplemental materials from PubMed Central.

**Usage**: `python harvest_pmc_fulltext.py`

**Documentation**: See `HARVESTER_QUICKSTART.md` and `README_HARVEST.md`

## 2. example_harvest_from_pubmind.py
Example script showing how to harvest articles from PMID lists.

**Usage**: `python example_harvest_from_pubmind.py`

## 3. pipeline.py / example_usage.py
Extract structured biomedical data from articles using AI.

**Usage**: `python example_usage.py`

**Documentation**: See `README_PIPELINE.md`

## 4. clinical_data_triage.py / example_triage.py
Process and triage clinical data from various sources.

**Usage**: `python example_triage.py`

**Documentation**: See `TRIAGE_README.md`
