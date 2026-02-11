"""Unified supplement retrieval system.

Tiered fetcher architecture:
    Tier 1 (free): PMC via NCBI eUtils / Europe PMC
    Tier 2 (API key): Elsevier, Wiley, Springer
    Tier 3 (scraping): Existing harvesting/supplement_scraper.py handlers

Public API:
    UnifiedSupplementFetcher - Tries all tiers, deduplicates results
    PMCSupplementFetcher    - Free PMC supplement access
    ElsevierSupplementFetcher - Elsevier API-based access
    SupplementFetcher       - Abstract base class
    SupplementFile          - Dataclass for supplement metadata
"""

from .base import SupplementFetcher, SupplementFile
from .elsevier_fetcher import ElsevierSupplementFetcher
from .pmc_fetcher import PMCSupplementFetcher
from .unified import UnifiedSupplementFetcher

__all__ = [
    "SupplementFetcher",
    "SupplementFile",
    "PMCSupplementFetcher",
    "ElsevierSupplementFetcher",
    "UnifiedSupplementFetcher",
]
