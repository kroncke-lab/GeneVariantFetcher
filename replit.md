# Overview

This is a Streamlit-based web application that interfaces with the NCBI LitVar2 API to search for genetic variants and associated PubMed publications by gene name. The application provides a user-friendly interface for researchers and scientists to explore genetic variant data and their literature references.

The core functionality revolves around querying the LitVar2 API endpoints to retrieve variant information (rsIDs) for specific genes and then fetching associated PubMed article identifiers (PMIDs) for those variants.

# User Preferences

Preferred communication style: Simple, everyday language.

# System Architecture

## Application Framework
- **Technology**: Streamlit
- **Rationale**: Streamlit provides rapid development of data-centric web applications with minimal frontend code. It's ideal for scientific tools and data exploration interfaces where the focus is on functionality rather than complex UI interactions.
- **Design Pattern**: Single-page application with reactive components

## API Integration Layer
The application implements a two-tier API integration strategy with the NCBI LitVar2 service:

1. **Gene Variant Search**: Uses the `/entity/search/{gene_name}` endpoint to retrieve variants associated with a given gene
2. **PMID Retrieval**: Uses the `/public/rsids2pmids` endpoint to fetch PubMed article identifiers for specific rsIDs (reference SNP identifiers)

**Error Handling Approach**: Implements try-catch blocks with specific handling for:
- Network request exceptions
- JSON parsing errors
- HTTP status errors via `raise_for_status()`
- Timeout protection (30s for search, 60s for PMID retrieval)

## Data Processing
- **Format**: JSON responses from API are processed into Python dictionaries
- **Data Transformation**: Pandas DataFrames are imported (suggesting tabular display of results)
- **CSV Export**: StringIO import indicates capability for data export functionality

## Frontend Architecture
- **Layout**: Wide page layout for better data visualization
- **Configuration**: Custom page title for browser tab identification
- **User Interface**: Markdown-based instructions and documentation

# External Dependencies

## Third-Party APIs
- **NCBI LitVar2 API** (https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1)
  - Primary data source for genetic variant information
  - No authentication required (public API)
  - Rate limiting considerations should be implemented for production use
  - Two main endpoints utilized:
    - Entity search: Gene name to variant mapping
    - rsIDs to PMIDs: Variant to literature mapping

## Python Libraries
- **streamlit**: Web application framework
- **requests**: HTTP client for API calls
- **pandas**: Data manipulation and tabular display
- **json**: JSON parsing and serialization

## Data Sources
- **PubMed**: Referenced through PMIDs for scientific literature
- **dbSNP**: Implicit dependency through rsID references (NCBI's database of genetic variation)