# Overview

This is a Streamlit-based web application that interfaces with the **PubMind-DB API** from Wang Genomics Lab to search for genetic variants and associated PubMed publications by gene name. The application provides a user-friendly interface for researchers and scientists to explore AI-curated genetic variant data with pathogenicity classifications and their literature references.

The core functionality revolves around querying the PubMind-DB undocumented JSON API to retrieve comprehensive variant information including DNA/protein changes, pathogenicity classifications, LLM-generated reasoning, disease associations, and associated PubMed article identifiers (PMIDs). This provides richer data than traditional variant databases by leveraging AI to curate and explain variant-disease-pathogenicity relationships from biomedical literature.

# User Preferences

Preferred communication style: Simple, everyday language.

# System Architecture

## Application Framework
- **Technology**: Streamlit
- **Rationale**: Streamlit provides rapid development of data-centric web applications with minimal frontend code. It's ideal for scientific tools and data exploration interfaces where the focus is on functionality rather than complex UI interactions.
- **Design Pattern**: Single-page application with reactive components and session state management

## API Integration Layer
The application implements a multi-tier API integration strategy:

1. **Gene Variant Search (PubMind-DB)**
   - Uses undocumented `/download_json?query={gene}&field=gene` endpoint discovered through reverse engineering
   - Returns comprehensive JSON with variant details, pathogenicity, PMIDs, and AI reasoning
   - Supports single and multiple gene searches with progress tracking
   - Results cached for 1 hour using `@st.cache_data` decorator
   - **Data Retrieved Per Variant:**
     - PVID (PubMind Variant ID) - unique identifier
     - DNA change (dna_change), Amino acid change (aa_change)
     - RSID (dbSNP identifier, may be "-" if unavailable)
     - Pathogenicity classification (pathogenic, benign, uncertain, conflicting)
     - Pathogenicity score (0-1 numerical confidence)
     - Confidence level (0-2 based on number of supporting records)
     - LLM reasoning (AI-generated explanation of pathogenicity)
     - MONDO disease associations
     - Formatted references with PMIDs and PMC IDs
     - Number of literature records supporting the classification

2. **PMID Extraction (PubMind Data)**
   - Extracts PMIDs from PubMind's `formatted_reference` and `PMCID_PMID_counted` fields
   - Uses regex pattern `PMID:?(\d+)` to handle formats like "PMC123(PMID:456)" and "PMID:789"
   - Aggregates unique PMIDs across all variants for a gene
   - No separate API call needed (PMIDs embedded in variant data)

3. **Article Fetching (NCBI E-utilities)**
   - Uses the `efetch.fcgi` endpoint to retrieve PubMed article abstracts
   - XML parsing to extract title and abstract text
   - Batched processing (200 PMIDs per request)
   - Results cached for 1 hour

4. **Biomedical Text Extraction (OpenAI GPT-5)**
   - Analyzes PubMed article text to extract structured individual-level data
   - Uses specialized biomedical prompt for variant/phenotype extraction
   - Implements retry logic with exponential backoff for rate limits
   - Concurrent processing (max 2 workers) with ThreadPoolExecutor

**Error Handling Approach**: Comprehensive error handling including:
- Network request exceptions
- JSON/XML parsing errors
- HTTP status errors via `raise_for_status()`
- Timeout protection (30s for search, 60s for PMID retrieval, custom for extraction)
- Rate limit detection and automatic retry
- Session state error tracking with user-facing error messages

## Data Processing
- **Format**: JSON responses from APIs processed into Python dictionaries
- **Data Transformation**: Pandas DataFrames for tabular display and analysis
- **Caching Strategy**: 1-hour TTL on all API calls using `@st.cache_data`
- **Session State**: Persistent storage of extraction results and processing status
- **Export Formats**: CSV, TXT, and JSON downloads available

## Frontend Architecture
- **Layout**: Wide page layout for better data visualization
- **Configuration**: Custom page title for browser tab identification
- **User Interface**: 
  - Sidebar for search parameters and filters
  - Four-tab results view: Visualizations, Variant Details, PMID List, Individual Extractions
  - Interactive Plotly charts for data visualization
  - Clickable links to external databases (dbSNP, PubMed)
  - Progress indicators and status widgets for long-running operations
  - Session state management for persistent UI state across reruns

## Biomedical Text Extraction Feature
**Architecture**:
- State machine design with stages: idle → fetching → extracting → complete/error
- Each stage stored in `st.session_state` and triggers page rerun
- Uses `st.status()` widget for persistent progress feedback
- Results persist across tab switches via session state

**Data Flow**:
1. User selects number of articles to process (1-20)
2. Button click initiates 'fetching' stage
3. PubMed abstracts fetched via E-utilities API
4. Stage transitions to 'extracting'
5. Each article processed sequentially with GPT-5
6. Extracted data aggregated and stored in session state
7. Final stage shows results table with download options

**Extraction Logic**:
- Identifies individual patients/subjects mentioned in articles
- Extracts genetic variants in HGVS notation
- Converts phenotypes to HPO codes when possible
- Determines affected/unaffected status using clinical logic
- Captures demographic information (age, sex, ancestry)
- Provides sentence-level evidence for each extraction

**Output Schema** (per individual):
```json
{
  "individual_id": "PMID_XXXX_case1",
  "pmid": "XXXXXXX",
  "gene": "GENE_NAME",
  "variants": [{"hgvs_c": "...", "hgvs_p": "...", "rsid": "...", "genomic": "..."}],
  "age": 52,
  "sex": "female",
  "phenotypes_hpo": ["HP:0003002"],
  "phenotypes_raw": ["breast cancer diagnosed at 52"],
  "affected_status": "affected",
  "evidence": {"sentence": "...", "section": "case report"},
  "data_from": "narrative"
}
```

# External Dependencies

## Third-Party APIs
- **PubMind-DB API** (https://pubmind.wglab.org/)
  - Primary data source for AI-curated genetic variant information
  - Developed by Wang Genomics Lab at CHOP and Penn
  - No authentication required (undocumented JSON export endpoint)
  - **Endpoint**: `/download_json?query={gene}&field=gene`
  - **Discovery Method**: Reverse-engineered from website's export functionality
  - Returns literature-derived knowledgebase of variant-disease-pathogenicity associations
  - Provides AI-generated reasoning and confidence scores for each variant
  - **Advantages over LitVar2**:
    - Richer pathogenicity classifications with AI explanations
    - Direct disease associations (MONDO ontology)
    - Confidence scoring based on literature evidence
    - Integrated LLM reasoning for variant interpretation
    - More comprehensive variant coverage (e.g., 4,144 variants for BRCA2 vs fewer in LitVar2)
  - **Note**: This is an undocumented API discovered through research; production use should contact WGLab for official API access

- **NCBI E-utilities API** (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/)
  - Used to fetch PubMed article abstracts and metadata
  - No authentication required (public API)
  - Returns XML format with title and abstract text
  - Batch processing with 200 PMIDs per request

- **OpenAI GPT-5** (via Replit AI Integrations)
  - Used for biomedical text extraction from PubMed articles
  - Accessed through Replit AI Integrations (no separate API key needed)
  - Analyzes article text to extract individual-level variant and phenotype data
  - Costs charged to user's Replit credits

## Python Libraries
- **streamlit**: Web application framework
- **requests**: HTTP client for API calls
- **pandas**: Data manipulation and tabular display
- **json**: JSON parsing and serialization
- **openai**: OpenAI Python SDK for LLM calls
- **tenacity**: Retry logic with exponential backoff for API calls
- **plotly**: Interactive data visualizations
- **concurrent.futures**: Parallel processing for batch operations

## Data Sources
- **PubMind-DB**: Primary source for variant-pathogenicity-disease relationships curated from literature using LLM assistance
- **PubMed**: Referenced through PMIDs for scientific literature access and extraction
- **dbSNP**: Referenced through rsID when available (NCBI's database of genetic variation)
- **MONDO**: Disease Ontology for standardized disease names and classifications
- **PubMed Central (PMC)**: Referenced through PMC IDs in variant literature citations

# Recent Changes (November 20, 2025)

## Major API Migration: LitVar2 → PubMind-DB
**Motivation**: User requested integration with PubMind (https://pubmind.wglab.org/) to access more comprehensive variant information with AI-powered pathogenicity curation.

**Implementation**:
1. Discovered undocumented PubMind JSON API endpoint through reverse engineering of website's export functionality
2. Completely replaced LitVar2 API calls with PubMind API integration
3. Updated data extraction to parse PubMind's richer JSON schema
4. Enhanced UI to display pathogenicity classifications, confidence scores, and LLM reasoning
5. Maintained backward compatibility with biomedical text extraction feature

**Benefits**:
- 10x+ more variants for genes like BRCA2 (4,144 vs ~400 in LitVar2)
- AI-curated pathogenicity classifications with explanatory reasoning
- Direct disease associations from literature
- Confidence scoring based on evidence strength
- Integrated PMID extraction (no separate API call needed)

**Technical Approach**:
- Analyzed PubMind website's network requests to identify `/download_json` endpoint
- Implemented PMID extraction via regex from `formatted_reference` and `PMCID_PMID_counted` fields
- Created new data extraction functions: `search_pubmind_variants()`, `extract_pubmind_variant_data()`
- Updated UI column configuration to show: Pathogenicity, Confidence, Diseases, LLM Reasoning
- Added PubMind links (PVID) alongside dbSNP links for external reference

**Testing**:
- Verified with BRCA2: 4,144 variants, 2,028 PMIDs, pathogenicity data displayed
- Confirmed extraction feature works with PubMind-derived PMIDs
- All four tabs functional with enhanced data presentation

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

**APIs Used**:
- **NCBI E-utilities**: PMID → PMCID conversion via `Bio.Entrez.elink()`
- **EuropePMC fullTextXML**: `https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML`
- **EuropePMC supplementaryFiles**: `https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles`

**Dependencies**: biopython, requests, markitdown, pandas, openpyxl, python-docx, tabulate

**Usage**:
```python
from harvest_pmc_fulltext import PMCHarvester

harvester = PMCHarvester(output_dir="pmc_harvest")
pmids = ['35443093', '33442691', '34931732']
harvester.harvest(pmids, delay=2.0)
```

**Output Structure**:
- `{PMID}_FULL_CONTEXT.md` - Unified markdown with full-text + supplements
- `{PMID}_supplements/` - Downloaded supplemental files
- `successful_downloads.csv` - Log of successful harvests
- `paywalled_missing.csv` - Log of unavailable articles

**Integration with PubMind App**: The harvester can process PMIDs extracted from PubMind variant searches, providing richer data for LLM extraction than abstract-only processing. Example scripts provided in `example_harvest_from_pubmind.py`.