# Overview

This is a Streamlit-based web application that interfaces with the NCBI LitVar2 API to search for genetic variants and associated PubMed publications by gene name. The application provides a user-friendly interface for researchers and scientists to explore genetic variant data and their literature references.

The core functionality revolves around querying the LitVar2 API endpoints to retrieve variant information (rsIDs) for specific genes and then fetching associated PubMed article identifiers (PMIDs) for those variants.

# User Preferences

Preferred communication style: Simple, everyday language.

# System Architecture

## Application Framework
- **Technology**: Streamlit
- **Rationale**: Streamlit provides rapid development of data-centric web applications with minimal frontend code. It's ideal for scientific tools and data exploration interfaces where the focus is on functionality rather than complex UI interactions.
- **Design Pattern**: Single-page application with reactive components and session state management

## API Integration Layer
The application implements a multi-tier API integration strategy:

1. **Gene Variant Search (NCBI LitVar2)**
   - Uses the `/entity/search/{gene_name}` endpoint to retrieve variants associated with a given gene
   - Supports single and multiple gene searches with progress tracking
   - Results cached for 1 hour using `@st.cache_data` decorator

2. **PMID Retrieval (NCBI LitVar2)**
   - Uses the `/public/rsids2pmids` endpoint to fetch PubMed article identifiers for specific rsIDs
   - Batched requests (100 rsIDs per batch) to avoid URL length limits
   - Per-batch error handling for resilience
   - Results cached for 1 hour

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
- **NCBI LitVar2 API** (https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1)
  - Primary data source for genetic variant information
  - No authentication required (public API)
  - Rate limiting considerations should be implemented for production use
  - Two main endpoints utilized:
    - Entity search: Gene name to variant mapping
    - rsIDs to PMIDs: Variant to literature mapping

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
- **PubMed**: Referenced through PMIDs for scientific literature
- **dbSNP**: Implicit dependency through rsID references (NCBI's database of genetic variation)