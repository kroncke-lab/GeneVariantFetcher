"""
Shared utility modules for GeneVariantFetcher.

This package provides common functionality used across multiple modules:
- llm_utils: LLM calling and response parsing
- html_utils: HTML parsing and PMID extraction
- retry_utils: Retry configuration and decorators
- pubmed_utils: Unified PubMed/NCBI API access
- pmid_utils: PMID extraction and validation
- logging_utils: Centralized logging configuration
"""

from .html_utils import create_scraping_session, extract_pmids_from_html
from .llm_utils import BaseLLMCaller, parse_llm_json_response
from .logging_utils import get_logger, setup_logging
from .pmid_utils import (
    extract_pmid_from_filename,
    extract_pmids_from_text,
    is_valid_pmid,
)
from .pubmed_utils import fetch_paper_metadata, query_pubmed_with_entrez
from .retry_utils import get_standard_retry_decorator, standard_retry

__all__ = [
    "BaseLLMCaller",
    "parse_llm_json_response",
    "extract_pmids_from_html",
    "create_scraping_session",
    "get_standard_retry_decorator",
    "standard_retry",
    "query_pubmed_with_entrez",
    "fetch_paper_metadata",
    "extract_pmid_from_filename",
    "is_valid_pmid",
    "extract_pmids_from_text",
    "setup_logging",
    "get_logger",
]
