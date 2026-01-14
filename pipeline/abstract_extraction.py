"""
Abstract-based carrier extraction module.

Extracts affected/unaffected carrier counts from paper abstracts using LLM.
This allows prioritizing papers for manual follow-up when full-text is unavailable.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional

from utils.llm_utils import BaseLLMCaller
from utils.models import Paper
from config.settings import get_settings

logger = logging.getLogger(__name__)


@dataclass
class AbstractCarrierResult:
    """Result of abstract-based carrier extraction."""

    pmid: str
    success: bool

    # Carrier counts extracted from abstract
    total_carriers: Optional[int] = None
    affected_count: Optional[int] = None
    unaffected_count: Optional[int] = None

    # Confidence and metadata
    confidence: float = 0.0
    extraction_notes: str = ""

    # Variant information if available
    variants_mentioned: list = field(default_factory=list)

    # Error information
    error: Optional[str] = None

    @property
    def has_carrier_counts(self) -> bool:
        """Returns True if any carrier counts were successfully extracted."""
        return (
            self.total_carriers is not None
            or self.affected_count is not None
            or self.unaffected_count is not None
        )

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "pmid": self.pmid,
            "success": self.success,
            "has_carrier_counts": self.has_carrier_counts,
            "total_carriers": self.total_carriers,
            "affected_count": self.affected_count,
            "unaffected_count": self.unaffected_count,
            "confidence": self.confidence,
            "extraction_notes": self.extraction_notes,
            "variants_mentioned": self.variants_mentioned,
            "error": self.error,
        }


# LLM prompt for abstract carrier extraction
ABSTRACT_CARRIER_EXTRACTION_PROMPT = """You are a medical geneticist extracting carrier count information from a paper abstract.

TARGET GENE: {gene_symbol}

Paper Title: {title}

Abstract:
{abstract}

TASK: Determine if this abstract contains extractable carrier counts for {gene_symbol} variants.

Look for:
1. Total number of carriers/patients with variants in {gene_symbol}
2. Number of AFFECTED carriers (those with disease/phenotype)
3. Number of UNAFFECTED carriers (healthy/asymptomatic carriers)
4. Any specific variant notations mentioned (e.g., c.1234G>A, p.Arg412His)

IMPORTANT DISTINCTIONS:
- "Affected" = has the disease phenotype or clinical manifestations
- "Unaffected" = carrier without disease phenotype (healthy carrier, asymptomatic)
- If only patient counts are given without distinguishing affected/unaffected, report as total_carriers only
- For disease-associated studies (case reports, clinical cohorts), patients are typically AFFECTED carriers
- For population/screening studies, distinguish affected vs unaffected if stated

Return a JSON object:
{{
    "extraction_possible": true or false,
    "total_carriers": integer or null (total number of carriers mentioned),
    "affected_count": integer or null (carriers with disease),
    "unaffected_count": integer or null (carriers without disease),
    "variants_mentioned": ["list of specific variant notations if any"],
    "confidence": 0.0-1.0 (how confident are you in these counts),
    "notes": "Brief explanation of what was found or why extraction wasn't possible"
}}

CRITICAL RULES:
- Only extract counts that are EXPLICITLY stated in the abstract
- Do NOT infer or calculate counts that aren't directly stated
- If the abstract mentions variants but no specific carrier counts, set extraction_possible to false
- If numbers are ambiguous (could be patients, families, or variants), note the ambiguity
- Set confidence lower if counts are unclear or potentially overlapping
"""


class AbstractCarrierExtractor(BaseLLMCaller):
    """
    Extracts carrier counts from paper abstracts using LLM.

    This is a lightweight extraction step that runs after filtering but before
    full-text download, to identify papers where carrier information can be
    obtained from the abstract alone.
    """

    def __init__(
        self,
        model: Optional[str] = None,
        temperature: Optional[float] = None,
        max_tokens: Optional[int] = None,
    ):
        """
        Initialize the abstract carrier extractor.

        Args:
            model: LLM model to use. Defaults to tier2_model from settings.
            temperature: Model temperature. Defaults to 0.2 for consistent extraction.
            max_tokens: Maximum tokens. Defaults to 500 (abstracts need less).
        """
        settings = get_settings()

        model = model or settings.tier2_model or "gpt-4o-mini"
        temperature = temperature if temperature is not None else 0.2
        max_tokens = max_tokens if max_tokens is not None else 500

        super().__init__(model=model, temperature=temperature, max_tokens=max_tokens)

        logger.debug(f"AbstractCarrierExtractor initialized with model={model}")

    def extract(
        self,
        paper: Paper,
        gene_symbol: Optional[str] = None,
    ) -> AbstractCarrierResult:
        """
        Extract carrier counts from a paper's abstract.

        Args:
            paper: Paper object with title and abstract.
            gene_symbol: Target gene symbol. Uses paper.gene_symbol if not provided.

        Returns:
            AbstractCarrierResult with extraction results.
        """
        gene = gene_symbol or paper.gene_symbol or "the gene"

        if not paper.abstract:
            logger.warning(f"PMID {paper.pmid}: No abstract available for carrier extraction")
            return AbstractCarrierResult(
                pmid=paper.pmid,
                success=False,
                error="No abstract available",
            )

        if not paper.title:
            logger.warning(f"PMID {paper.pmid}: No title available for carrier extraction")
            # Continue anyway with empty title

        # Construct prompt
        prompt = ABSTRACT_CARRIER_EXTRACTION_PROMPT.format(
            gene_symbol=gene,
            title=paper.title or "Unknown Title",
            abstract=paper.abstract[:3000],  # Truncate very long abstracts
        )

        try:
            logger.debug(f"PMID {paper.pmid}: Extracting carrier counts from abstract")

            result_data = self.call_llm_json(prompt)

            extraction_possible = result_data.get("extraction_possible", False)

            if not extraction_possible:
                return AbstractCarrierResult(
                    pmid=paper.pmid,
                    success=True,  # Extraction ran successfully, just no counts found
                    confidence=result_data.get("confidence", 0.0),
                    extraction_notes=result_data.get("notes", "No carrier counts found in abstract"),
                    variants_mentioned=result_data.get("variants_mentioned", []),
                )

            # Extract counts
            total_carriers = result_data.get("total_carriers")
            affected_count = result_data.get("affected_count")
            unaffected_count = result_data.get("unaffected_count")

            # Validate: convert non-integer values to None
            if total_carriers is not None and not isinstance(total_carriers, int):
                try:
                    total_carriers = int(total_carriers)
                except (ValueError, TypeError):
                    total_carriers = None

            if affected_count is not None and not isinstance(affected_count, int):
                try:
                    affected_count = int(affected_count)
                except (ValueError, TypeError):
                    affected_count = None

            if unaffected_count is not None and not isinstance(unaffected_count, int):
                try:
                    unaffected_count = int(unaffected_count)
                except (ValueError, TypeError):
                    unaffected_count = None

            result = AbstractCarrierResult(
                pmid=paper.pmid,
                success=True,
                total_carriers=total_carriers,
                affected_count=affected_count,
                unaffected_count=unaffected_count,
                confidence=float(result_data.get("confidence", 0.5)),
                extraction_notes=result_data.get("notes", ""),
                variants_mentioned=result_data.get("variants_mentioned", []),
            )

            logger.info(
                f"PMID {paper.pmid}: Abstract extraction - "
                f"total={total_carriers}, affected={affected_count}, "
                f"unaffected={unaffected_count} (confidence: {result.confidence:.2f})"
            )

            return result

        except Exception as e:
            logger.error(f"PMID {paper.pmid}: Abstract carrier extraction failed: {e}")
            return AbstractCarrierResult(
                pmid=paper.pmid,
                success=False,
                error=str(e),
            )

    def extract_from_text(
        self,
        pmid: str,
        title: str,
        abstract: str,
        gene_symbol: str,
    ) -> AbstractCarrierResult:
        """
        Extract carrier counts from raw text (convenience method).

        Args:
            pmid: PubMed ID.
            title: Paper title.
            abstract: Paper abstract.
            gene_symbol: Target gene symbol.

        Returns:
            AbstractCarrierResult with extraction results.
        """
        paper = Paper(
            pmid=pmid,
            title=title,
            abstract=abstract,
            gene_symbol=gene_symbol,
        )
        return self.extract(paper, gene_symbol)


def batch_extract_from_abstracts(
    papers: list[Paper],
    gene_symbol: str,
    max_workers: int = 4,
) -> dict[str, AbstractCarrierResult]:
    """
    Extract carrier counts from multiple papers in parallel.

    Args:
        papers: List of Paper objects with abstracts.
        gene_symbol: Target gene symbol.
        max_workers: Maximum parallel workers.

    Returns:
        Dictionary mapping PMID to AbstractCarrierResult.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    extractor = AbstractCarrierExtractor()
    results = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pmid = {
            executor.submit(extractor.extract, paper, gene_symbol): paper.pmid
            for paper in papers
        }

        for future in as_completed(future_to_pmid):
            pmid = future_to_pmid[future]
            try:
                result = future.result()
                results[pmid] = result
            except Exception as e:
                logger.error(f"PMID {pmid}: Batch extraction failed: {e}")
                results[pmid] = AbstractCarrierResult(
                    pmid=pmid,
                    success=False,
                    error=str(e),
                )

    return results
