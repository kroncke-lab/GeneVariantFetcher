"""Gene-focused literature discovery utilities.

This package provides tools for discovering relevant papers about genes
from PubMind and PubMed, plus clinical data triage filtering.

Public API:
    PubMindFetcher - Fetches PMIDs from PubMind for genes/variants
    fetch_pmids_for_gene - Convenience function for gene-based discovery
    fetch_pmids_for_variant - Convenience function for variant-based discovery
"""

from .pubmind_fetcher import (
    PubMindFetcher,
    fetch_pmids_for_gene,
    fetch_pmids_for_variant,
)

__all__ = [
    "PubMindFetcher",
    "fetch_pmids_for_gene",
    "fetch_pmids_for_variant",
]
