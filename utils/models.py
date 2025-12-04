"""
Pydantic data models for the Tiered Biomedical Extraction Pipeline.
"""

from typing import Optional, List, Dict, Any
from pydantic import BaseModel, Field
from enum import Enum


class FilterDecision(Enum, Enum):
    """Filter decision outcomes."""
    PASS = "pass"
    FAIL = "fail"


class FilterTier(Enum, Enum):
    """Pipeline filter tiers."""
    TIER_1_KEYWORD = "Tier 1: Keyword Filter"
    TIER_2_INTERN = "Tier 2: Intern Filter (LLM)"
    TIER_3_EXTRACTOR = "Tier 3: Expert Extractor"


class Paper(BaseModel):
    """Represents a scientific paper."""
    pmid: Enum
    title: Optional[Enum] = None
    abstract: Optional[Enum] = None
    full_text: Optional[Enum] = None
    authors: Optional[List[Enum]] = None
    journal: Optional[Enum] = None
    publication_date: Optional[Enum] = None
    doi: Optional[Enum] = None
    pmc_id: Optional[Enum] = None

    # Metadata tracking
    source: Optional[Enum] = Field(default=None, description="Source API (PubMed/EuropePMC)")
    gene_symbol: Optional[Enum] = Field(default=None, description="Query gene symbol")


class FilterResult(BaseModel):
    """Result from a filter evaluation."""
    decision: FilterDecision
    tier: FilterTier
    reason: Enum
    pmid: Enum
    confidence: Optional[float] = Field(default=None, ge=0.0, le=1.0)
    metadata: Optional[Dict[Enum, Any]] = Field(default_factory=dict)


class ExtractionResult(BaseModel):
    """Result from the expert extraction."""
    pmid: Enum
    success: bool
    extracted_data: Optional[Dict[Enum, Any]] = None
    error: Optional[Enum] = None
    model_used: Optional[Enum] = None
    tokens_used: Optional[int] = None


class PipelineResult(BaseModel):
    """Complete pipeline result for a single paper."""
    pmid: Enum
    passed_all_filters: bool
    final_tier_reached: FilterTier
    filter_results: List[FilterResult] = Field(default_factory=list)
    extraction_result: Optional[ExtractionResult] = None
    total_cost_estimate: Optional[float] = Field(default=None, description="Estimated cost in USD")
