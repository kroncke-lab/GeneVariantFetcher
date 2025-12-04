"""
Pydantic data models for the Tiered Biomedical Extraction Pipeline.
"""

from typing import Optional, List, Dict, Any
from pydantic import BaseModel, Field
from enum import Enum


class FilterDecision(str, Enum):
    """Filter decision outcomes."""
    PASS = "pass"
    FAIL = "fail"


class FilterTier(str, Enum):
    """Pipeline filter tiers."""
    TIER_1_KEYWORD = "Tier 1: Keyword Filter"
    TIER_2_INTERN = "Tier 2: Intern Filter (LLM)"
    TIER_3_EXTRACTOR = "Tier 3: Expert Extractor"


class Paper(BaseModel):
    """Represents a scientific paper."""
    pmid: str
    title: Optional[str] = None
    abstract: Optional[str] = None
    full_text: Optional[str] = None
    authors: Optional[List[str]] = None
    journal: Optional[str] = None
    publication_date: Optional[str] = None
    doi: Optional[str] = None
    pmc_id: Optional[str] = None

    # Metadata tracking
    source: Optional[str] = Field(default=None, description="Source API (PubMed/EuropePMC)")
    gene_symbol: Optional[str] = Field(default=None, description="Query gene symbol")


class FilterResult(BaseModel):
    """Result from a filter evaluation."""
    decision: FilterDecision
    tier: FilterTier
    reason: str
    pmid: str
    confidence: Optional[float] = Field(default=None, ge=0.0, le=1.0)
    metadata: Optional[Dict[str, Any]] = Field(default_factory=dict)


class ExtractionResult(BaseModel):
    """Result from the expert extraction."""
    pmid: str
    success: bool
    extracted_data: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    model_used: Optional[str] = None
    tokens_used: Optional[int] = None


class PipelineResult(BaseModel):
    """Complete pipeline result for a single paper."""
    pmid: str
    passed_all_filters: bool
    final_tier_reached: FilterTier
    filter_results: List[FilterResult] = Field(default_factory=list)
    extraction_result: Optional[ExtractionResult] = None
    total_cost_estimate: Optional[float] = Field(default=None, description="Estimated cost in USD")
