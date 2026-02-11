"""
Gene Variant Fetcher Pipeline Package

This package contains the core pipeline components for extracting genetic variant data
from scientific literature through a tiered filtering and extraction process.

Pipeline stages:
    1. Filters (KeywordFilter, InternFilter) - Filter papers by relevance
    2. Extraction (ExpertExtractor) - Extract variant data from papers
    3. Aggregation (DataAggregator) - Aggregate and validate results
"""

from .aggregation import DataAggregator
from .extraction import ExpertExtractor
from .filters import ClinicalDataTriageFilter, InternFilter, KeywordFilter

__all__ = [
    "KeywordFilter",
    "InternFilter",
    "ClinicalDataTriageFilter",
    "ExpertExtractor",
    "DataAggregator",
]
