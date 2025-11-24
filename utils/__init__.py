"""
Shared utility modules for GeneVariantFetcher.

This package provides common functionality used across multiple modules:
- llm_utils: LLM calling and response parsing
- html_utils: HTML parsing and PMID extraction
- retry_utils: Retry configuration and decorators
- pubmed_utils: Unified PubMed/NCBI API access
"""

from .llm_utils import BaseLLMCaller, parse_llm_json_response
from .html_utils import extract_pmids_from_html, create_scraping_session
from .retry_utils import get_standard_retry_decorator, standard_retry
from .pubmed_utils import query_pubmed_with_entrez, fetch_paper_metadata

__all__ = [
    "BaseLLMCaller",
    "parse_llm_json_response",
    "extract_pmids_from_html",
    "create_scraping_session",
    "get_standard_retry_decorator",
    "standard_retry",
    "query_pubmed_with_entrez",
    "fetch_paper_metadata",
]
