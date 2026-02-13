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

# Max characters to send to the LLM in a single extraction prompt
TEXT_TRUNCATION_MAX_CHARS: int = 60_000

# Lines of context to include around gene mentions during focused truncation
GENE_CONTEXT_WINDOW_LINES: int = 80

# Minimum number of pipe-delimited columns for a line to count as a table row
TABLE_MIN_COLUMNS: int = 3

# Table row count above which we try the deterministic parser
LARGE_TABLE_ROW_THRESHOLD: int = 100

# Deterministic parser minimum variant yield to bypass LLM
DETERMINISTIC_PARSER_MIN_VARIANTS: int = 50

# Minimum variant gap between extracted and expected to trigger continuation
CONTINUATION_VARIANT_GAP: int = 5

# Scanner merge minimum confidence threshold
SCANNER_MERGE_MIN_CONFIDENCE: float = 0.6

# Maximum number of scanner hints to include in the LLM prompt
SCANNER_MAX_HINTS: int = 50

# Table row hint threshold for raising the adaptive model-cascade bar
ADAPTIVE_TABLE_THRESHOLD: int = 50

# =============================================================================
# HTTP Timeouts (seconds)
# =============================================================================

HTTP_TIMEOUT_DEFAULT: int = 30
HTTP_TIMEOUT_CONVERSION: int = 120
HTTP_TIMEOUT_SHORT: int = 10

# =============================================================================
# Concurrency
# =============================================================================

DEFAULT_MAX_WORKERS: int = 8

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
