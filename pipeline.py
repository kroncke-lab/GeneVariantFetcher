"""
Pipeline Orchestrator for the Tiered Biomedical Extraction Pipeline.

Coordinates the flow: Sourcer → Keyword Filter → Intern Filter → Expert Extractor
CRITICAL: Discards papers at each tier to save processing costs and time.
"""

import logging
from typing import List, Optional, Dict, Any
from dataclasses import dataclass
from datetime import datetime

from models import (
    Paper,
    PipelineResult,
    FilterResult,
    FilterDecision,
    FilterTier,
    ExtractionResult
)
from pipeline.sourcing import PaperSourcer
from pipeline.filters import KeywordFilter, InternFilter
from pipeline.extraction import ExpertExtractor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class PipelineStats:
    """Statistics for a pipeline run."""
    total_papers_sourced: int = 0
    passed_tier1_keyword: int = 0
    passed_tier2_intern: int = 0
    completed_tier3_extraction: int = 0
    failed_tier1: int = 0
    failed_tier2: int = 0
    failed_tier3: int = 0
    total_extraction_successes: int = 0
    total_variants_extracted: int = 0
    estimated_cost_savings: float = 0.0  # In USD
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None


class BiomedicalExtractionPipeline:
    """
    The main pipeline orchestrator that connects all components.

    Flow:
    1. Sourcer: Query APIs for PMIDs
    2. Tier 1 (Keyword Filter): Fast keyword check on abstracts
       → DROP papers that fail
    3. Tier 2 (Intern Filter): LLM classification with gpt-4o-mini
       → DROP papers that fail
    4. Tier 3 (Expert Extractor): Deep extraction with advanced LLM
    """

    def __init__(
        self,
        email: str = "your_email@example.com",
        keyword_filter: Optional[KeywordFilter] = None,
        intern_filter: Optional[InternFilter] = None,
        expert_extractor: Optional[ExpertExtractor] = None,
        enable_tier1: bool = True,
        enable_tier2: bool = True,
        enable_tier3: bool = True
    ):
        """
        Initialize the pipeline with all components.

        Args:
            email: Email for PubMed API.
            keyword_filter: Custom keyword filter (or uses default).
            intern_filter: Custom intern filter (or uses default).
            expert_extractor: Custom expert extractor (or uses default).
            enable_tier1: Enable Tier 1 keyword filtering.
            enable_tier2: Enable Tier 2 intern filtering.
            enable_tier3: Enable Tier 3 extraction.
        """
        self.sourcer = PaperSourcer(email=email)

        # Initialize filters with defaults if not provided
        self.keyword_filter = keyword_filter or KeywordFilter(min_keyword_matches=2)
        self.intern_filter = intern_filter or InternFilter(model="gpt-4o-mini")
        self.expert_extractor = expert_extractor or ExpertExtractor(model="gpt-4o")

        # Tier toggles
        self.enable_tier1 = enable_tier1
        self.enable_tier2 = enable_tier2
        self.enable_tier3 = enable_tier3

        self.stats = PipelineStats()

    def _log_paper_dropped(self, pmid: str, tier: FilterTier, reason: str):
        """
        Log when a paper is dropped at a specific tier.

        Args:
            pmid: Paper PMID.
            tier: Filter tier where it was dropped.
            reason: Reason for dropping.
        """
        logger.info(f"🚫 DROPPED - PMID {pmid} at {tier.value}: {reason}")

    def _log_paper_passed(self, pmid: str, tier: FilterTier):
        """
        Log when a paper passes a tier.

        Args:
            pmid: Paper PMID.
            tier: Filter tier that was passed.
        """
        logger.info(f"✓ PASSED - PMID {pmid} through {tier.value}")

    def _estimate_cost_savings(self, papers_dropped_tier2: int, papers_dropped_tier1: int) -> float:
        """
        Estimate cost savings from early filtering.

        Assumptions:
        - Tier 2 (Intern): ~$0.0001 per paper (gpt-4o-mini, ~500 tokens)
        - Tier 3 (Extractor): ~$0.05 per paper (gpt-4o, ~5000 tokens)

        Args:
            papers_dropped_tier2: Papers dropped after Tier 2.
            papers_dropped_tier1: Papers dropped after Tier 1.

        Returns:
            Estimated savings in USD.
        """
        cost_per_tier2 = 0.0001
        cost_per_tier3 = 0.05

        # Papers dropped at Tier 1 saved both Tier 2 and Tier 3 costs
        savings_tier1 = papers_dropped_tier1 * (cost_per_tier2 + cost_per_tier3)

        # Papers dropped at Tier 2 saved Tier 3 costs
        savings_tier2 = papers_dropped_tier2 * cost_per_tier3

        return savings_tier1 + savings_tier2

    def process_paper(self, paper: Paper) -> PipelineResult:
        """
        Process a single paper through all pipeline tiers.

        Args:
            paper: Paper object with metadata and abstract.

        Returns:
            PipelineResult with filter results and extraction data.
        """
        result = PipelineResult(
            pmid=paper.pmid,
            passed_all_filters=False,
            final_tier_reached=FilterTier.TIER_1_KEYWORD,
            filter_results=[],
            extraction_result=None
        )

        # ============================================
        # TIER 1: Keyword Filter (Fast & Cheap)
        # ============================================
        if self.enable_tier1:
            tier1_result = self.keyword_filter.filter(paper)
            result.filter_results.append(tier1_result)

            if tier1_result.decision == FilterDecision.FAIL:
                self._log_paper_dropped(paper.pmid, FilterTier.TIER_1_KEYWORD, tier1_result.reason)
                self.stats.failed_tier1 += 1
                return result  # 🚫 DISCARD - Save money by not processing further

            self._log_paper_passed(paper.pmid, FilterTier.TIER_1_KEYWORD)
            self.stats.passed_tier1_keyword += 1

        # ============================================
        # TIER 2: Intern Filter (LLM Classification)
        # ============================================
        result.final_tier_reached = FilterTier.TIER_2_INTERN

        if self.enable_tier2:
            tier2_result = self.intern_filter.filter(paper)
            result.filter_results.append(tier2_result)

            if tier2_result.decision == FilterDecision.FAIL:
                self._log_paper_dropped(paper.pmid, FilterTier.TIER_2_INTERN, tier2_result.reason)
                self.stats.failed_tier2 += 1
                return result  # 🚫 DISCARD - Not worth expensive extraction

            self._log_paper_passed(paper.pmid, FilterTier.TIER_2_INTERN)
            self.stats.passed_tier2_intern += 1

        # ============================================
        # TIER 3: Expert Extraction (Heavy Lifting)
        # ============================================
        result.final_tier_reached = FilterTier.TIER_3_EXTRACTOR

        if self.enable_tier3:
            extraction_result = self.expert_extractor.extract(paper)
            result.extraction_result = extraction_result
            self.stats.completed_tier3_extraction += 1

            if extraction_result.success:
                self.stats.total_extraction_successes += 1

                # Count variants extracted
                variants_count = (
                    extraction_result.extracted_data
                    .get("extraction_metadata", {})
                    .get("total_variants_found", 0)
                )
                self.stats.total_variants_extracted += variants_count

                logger.info(
                    f"✓ EXTRACTED - PMID {paper.pmid}: {variants_count} variants found"
                )
            else:
                self.stats.failed_tier3 += 1
                logger.warning(
                    f"⚠ EXTRACTION FAILED - PMID {paper.pmid}: {extraction_result.error}"
                )

        # Made it through all tiers!
        result.passed_all_filters = True
        return result

    def run(
        self,
        gene_symbol: str,
        max_papers: Optional[int] = None,
        fetch_full_text: bool = False
    ) -> Dict[str, Any]:
        """
        Run the complete pipeline for a gene symbol.

        Args:
            gene_symbol: Gene symbol to search for.
            max_papers: Maximum number of papers to process (None = all).
            fetch_full_text: Whether to fetch full text (if False, uses abstracts).

        Returns:
            Dictionary with results and statistics.
        """
        self.stats = PipelineStats(start_time=datetime.now())

        logger.info(f"{'='*80}")
        logger.info(f"🚀 STARTING PIPELINE FOR GENE: {gene_symbol}")
        logger.info(f"{'='*80}")

        # ============================================
        # STEP 1: Source Papers
        # ============================================
        logger.info(f"📚 STEP 1: Sourcing papers from PubMind (primary), PubMed, and EuropePMC...")
        pmids = self.sourcer.fetch_papers(gene_symbol, use_pubmind=True)
        self.stats.total_papers_sourced = len(pmids)

        logger.info(f"Found {len(pmids)} unique PMIDs for {gene_symbol}")

        if max_papers:
            pmids = pmids[:max_papers]
            logger.info(f"Limiting to first {max_papers} papers")

        # ============================================
        # STEP 2: Fetch Metadata
        # ============================================
        logger.info(f"📥 STEP 2: Fetching metadata for {len(pmids)} papers...")
        papers: List[Paper] = []

        for pmid in pmids:
            paper = self.sourcer.fetch_paper_metadata(pmid)
            if paper:
                paper.gene_symbol = gene_symbol
                papers.append(paper)
            else:
                logger.warning(f"Could not fetch metadata for PMID {pmid}")

        logger.info(f"Successfully fetched metadata for {len(papers)} papers")

        # ============================================
        # STEP 3: Process Through Tiered Pipeline
        # ============================================
        logger.info(f"⚙️ STEP 3: Processing papers through tiered pipeline...")

        results: List[PipelineResult] = []

        for i, paper in enumerate(papers, 1):
            logger.info(f"\n--- Processing paper {i}/{len(papers)} - PMID {paper.pmid} ---")

            result = self.process_paper(paper)
            results.append(result)

        # ============================================
        # STEP 4: Calculate Statistics
        # ============================================
        self.stats.end_time = datetime.now()
        self.stats.estimated_cost_savings = self._estimate_cost_savings(
            self.stats.failed_tier2,
            self.stats.failed_tier1
        )

        # ============================================
        # STEP 5: Generate Report
        # ============================================
        logger.info(f"\n{'='*80}")
        logger.info(f"📊 PIPELINE SUMMARY FOR {gene_symbol}")
        logger.info(f"{'='*80}")
        logger.info(f"Total papers sourced: {self.stats.total_papers_sourced}")
        logger.info(f"Passed Tier 1 (Keyword): {self.stats.passed_tier1_keyword}")
        logger.info(f"Passed Tier 2 (Intern): {self.stats.passed_tier2_intern}")
        logger.info(f"Completed Tier 3 (Extraction): {self.stats.completed_tier3_extraction}")
        logger.info(f"")
        logger.info(f"❌ Dropped at Tier 1: {self.stats.failed_tier1}")
        logger.info(f"❌ Dropped at Tier 2: {self.stats.failed_tier2}")
        logger.info(f"❌ Failed at Tier 3: {self.stats.failed_tier3}")
        logger.info(f"")
        logger.info(f"✓ Successful extractions: {self.stats.total_extraction_successes}")
        logger.info(f"🧬 Total variants extracted: {self.stats.total_variants_extracted}")
        logger.info(f"💰 Estimated cost savings: ${self.stats.estimated_cost_savings:.2f}")

        if self.stats.start_time and self.stats.end_time:
            duration = (self.stats.end_time - self.stats.start_time).total_seconds()
            logger.info(f"⏱️ Total processing time: {duration:.1f} seconds")

        logger.info(f"{'='*80}\n")

        return {
            "gene_symbol": gene_symbol,
            "stats": self.stats,
            "results": results,
            "successful_extractions": [
                r.extraction_result for r in results
                if r.extraction_result and r.extraction_result.success
            ]
        }


def run_pipeline_for_gene(
    gene_symbol: str,
    email: str = "your_email@example.com",
    max_papers: Optional[int] = None
) -> Dict[str, Any]:
    """
    Convenience function to run the pipeline for a gene.

    Args:
        gene_symbol: Gene symbol to search for.
        email: Email for PubMed API.
        max_papers: Maximum papers to process.

    Returns:
        Pipeline results and statistics.
    """
    pipeline = BiomedicalExtractionPipeline(email=email)
    return pipeline.run(gene_symbol, max_papers=max_papers)
