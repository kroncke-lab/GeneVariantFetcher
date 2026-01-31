"""Centralized constants for the Gene Variant Fetcher pipeline.

This module provides a single source of truth for keyword lists, thresholds,
and other magic values used across the codebase.
"""

from typing import List

# =============================================================================
# Clinical Keywords for Tier 1 Keyword Filter
# =============================================================================
# Used by pipeline/filters.py KeywordFilter for fast initial screening

FILTER_CLINICAL_KEYWORDS: List[str] = [
    # Variant/mutation terms
    "variant",
    "mutation",
    "polymorphism",
    "SNP",
    "deletion",
    "insertion",
    "substitution",
    "missense",
    "nonsense",
    "frameshift",
    "splice",
    "indel",
    "copy number",
    "CNV",
    # Clinical terms
    "patient",
    "patients",
    "clinical",
    "disease",
    "syndrome",
    "phenotype",
    "diagnosis",
    "treatment",
    "therapy",
    "outcome",
    "prognosis",
    "symptom",
    "symptoms",
    "manifestation",
    "pathogenic",
    "benign",
    # Study types
    "case report",
    "case series",
    "cohort",
    "clinical trial",
    "study",
    "analysis",
    "association",
    "genotype",
    # Medical/genetic terms
    "pathology",
    "molecular",
    "genetic",
    "genomic",
    "exome",
    "sequencing",
    "gene",
    "chromosome",
    "allele",
    "heterozygous",
    "homozygous",
    "carrier",
    "inheritance",
    "familial",
    "sporadic",
]

# =============================================================================
# Clinical Keywords for Data Scout Relevance Scoring
# =============================================================================
# Used by pipeline/data_scout.py for identifying high-value text zones

SCOUT_CLINICAL_KEYWORDS: List[str] = [
    "patient",
    "proband",
    "carrier",
    "affected",
    "unaffected",
    "phenotype",
    "family",
    "pedigree",
    "segregation",
    "penetrance",
    "symptomatic",
    "asymptomatic",
    "diagnosis",
    "onset",
    "presentation",
    "heterozygous",
    "homozygous",
    "compound heterozygous",
    "de novo",
]

# =============================================================================
# Extraction Thresholds
# =============================================================================

# Minimum size (chars) for DATA_ZONES.md to be used instead of FULL_CONTEXT.md
MIN_CONDENSED_SIZE: int = 500

# Minimum size (chars) for extraction input to be considered usable
# Below this threshold, LLM extraction is skipped (circuit breaker)
MIN_EXTRACTION_INPUT_SIZE: int = 500

# Minimum ratio of alphanumeric content to total content
# Below this indicates garbage/placeholder text
MIN_ALPHANUMERIC_RATIO: float = 0.3

# =============================================================================
# HTTP Request Constants
# =============================================================================

# Standard browser User-Agent for web requests
DEFAULT_USER_AGENT: str = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/120.0.0.0 Safari/537.36"
)

# macOS Chrome User-Agent (alternative)
MACOS_USER_AGENT: str = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/131.0.0.0 Safari/537.36"
)

# Browser-like headers for HTTP requests
BROWSER_HEADERS: dict[str, str] = {
    "User-Agent": DEFAULT_USER_AGENT,
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.5",
    "Accept-Encoding": "gzip, deflate, br",
    "Connection": "keep-alive",
    "Upgrade-Insecure-Requests": "1",
}

# =============================================================================
# Retry Configuration
# =============================================================================

# Default retry settings for HTTP requests
DEFAULT_MAX_RETRIES: int = 3
DEFAULT_BACKOFF_SECONDS: float = 1.0
