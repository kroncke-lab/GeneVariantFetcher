"""Unified LLM-based relevance checking for papers and synonyms.

This module consolidates the common LLM calling logic used for relevance
assessments, reducing code duplication between paper and synonym checking.
"""

from __future__ import annotations

import json
import logging
import os
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# =============================================================================
# Result Dataclasses
# =============================================================================


@dataclass
class RelevanceScore:
    """Relevance assessment for a paper."""

    is_relevant: bool
    confidence: float  # 0.0 to 1.0
    reasoning: str
    pmid: str


@dataclass
class SynonymRelevance:
    """Relevance assessment for a gene synonym."""

    synonym: str
    is_relevant: bool
    confidence: float  # 0.0 to 1.0
    reasoning: str


# =============================================================================
# Base LLM Relevance Checker
# =============================================================================


class BaseLLMRelevanceChecker(ABC):
    """Base class for LLM-based relevance checking.

    Handles common initialization, API setup, and response parsing.
    Subclasses implement specific prompts and result handling.
    """

    DEFAULT_MODEL = "claude-3-5-haiku-20241022"
    DEFAULT_MAX_TOKENS = 200

    def __init__(
        self,
        api_key: Optional[str] = None,
        model: Optional[str] = None,
    ) -> None:
        """Initialize the relevance checker.

        Args:
            api_key: Anthropic API key (or set ANTHROPIC_API_KEY env var)
            model: Claude model to use (haiku is cost-effective)
        """
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        self.model = model or self.DEFAULT_MODEL
        self._client = None

        if not self.api_key:
            logger.warning(
                "No Anthropic API key provided. Set ANTHROPIC_API_KEY "
                "environment variable or pass api_key parameter."
            )

    def _get_client(self):
        """Lazily initialize and return the Anthropic client."""
        if self._client is not None:
            return self._client

        try:
            import anthropic

            self._client = anthropic.Anthropic(api_key=self.api_key)
            return self._client
        except ImportError:
            logger.error(
                "anthropic package not installed. Install with: pip install anthropic"
            )
            return None

    def _call_llm(self, prompt: str) -> Optional[Dict[str, Any]]:
        """Make an LLM call and parse the JSON response.

        Args:
            prompt: The prompt to send to the LLM

        Returns:
            Parsed JSON dict or None if failed
        """
        if not self.api_key:
            return None

        client = self._get_client()
        if client is None:
            return None

        try:
            message = client.messages.create(
                model=self.model,
                max_tokens=self.DEFAULT_MAX_TOKENS,
                messages=[{"role": "user", "content": prompt}],
            )

            response_text = message.content[0].text.strip()

            # Extract JSON from response
            json_match = re.search(r"\{[^}]+\}", response_text, re.DOTALL)
            if json_match:
                return json.loads(json_match.group())
            else:
                logger.warning(f"Could not parse LLM response: {response_text[:100]}")
                return None

        except Exception as e:
            logger.error(f"LLM call failed: {e}")
            return None

    @abstractmethod
    def _build_prompt(self, **kwargs) -> str:
        """Build the prompt for the LLM. Implemented by subclasses."""
        pass

    @abstractmethod
    def _create_fallback_result(self, reason: str, **kwargs) -> Any:
        """Create a fallback result when LLM is unavailable. Implemented by subclasses."""
        pass


# =============================================================================
# Paper Relevance Checker
# =============================================================================


class RelevanceChecker(BaseLLMRelevanceChecker):
    """Check if papers are actually about the gene using LLM assessment."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        model: str = "claude-3-5-haiku-20241022",
        batch_size: int = 10,
    ) -> None:
        """Initialize the relevance checker.

        Args:
            api_key: Anthropic API key (or set ANTHROPIC_API_KEY env var)
            model: Claude model to use (haiku is cost-effective)
            batch_size: Number of papers to process in parallel (reserved for future use)
        """
        super().__init__(api_key=api_key, model=model)
        self.batch_size = batch_size

    def _build_prompt(
        self, gene_name: str, title: str, abstract: Optional[str], **kwargs
    ) -> str:
        """Build the paper relevance prompt."""
        text_to_check = f"Title: {title}\n"
        if abstract:
            text_to_check += f"\nAbstract: {abstract}"

        return f"""You are analyzing whether a scientific paper is actually about the gene "{gene_name}" and whether it likely contains variant- or individual-level genetic information useful for extraction.

Sometimes gene symbols match common abbreviations that have nothing to do with the gene. For example:
- "TTR" could mean the gene Transthyretin OR "time to reimbursement"
- "CAT" could mean the gene Catalase OR "computed axial tomography"

Analyze this paper and determine if it's genuinely about the gene {gene_name}. Only mark NOT RELEVANT when the abstract clearly shows the study is not about genetic variation, such as:
- Discussing the gene only in passing without analyzing genetic variants
- Focusing on expression, protein function, pathways, or other non-variant biology
- Describing group-level mutations with no hint of variant-level detail
- Covering non-genetic topics (clinical workflows, imaging, ML predictions without genotypes, epidemiology without variants)
- Using non-human organisms without human-disease genetic variation
- Mechanistic/functional biology with no specific mutations, alleles, genotypes, or patient cases
- Risk or outcome studies that lack individual genotypes or variant types

Do NOT exclude a paper when mutations/variants are mentioned generally, the genetic focus is unclear but plausible, or the abstract might hide useful variant details in tables/figures/supplements. When uncertain, treat it as relevant.

{text_to_check}

Respond with ONLY a JSON object with this exact format:
{{
  "is_relevant": true or false,
  "confidence": 0.0 to 1.0,
  "reasoning": "brief explanation"
}}

Lean toward relevance unless the abstract unmistakably indicates it is not about variant-level genetic information."""

    def _create_fallback_result(
        self, reason: str, pmid: str = "", **kwargs
    ) -> RelevanceScore:
        """Create a fallback result assuming relevance."""
        return RelevanceScore(
            is_relevant=True,
            confidence=0.5,
            reasoning=reason,
            pmid=pmid,
        )

    def check_relevance(
        self,
        gene_name: str,
        title: str,
        abstract: Optional[str],
        pmid: str,
    ) -> RelevanceScore:
        """Check if a single paper is relevant to the gene.

        Args:
            gene_name: The gene being searched for
            title: Paper title
            abstract: Paper abstract (if available)
            pmid: PubMed ID

        Returns:
            RelevanceScore with assessment
        """
        if not self.api_key:
            return self._create_fallback_result(
                "No API key provided; cannot assess relevance", pmid=pmid
            )

        if self._get_client() is None:
            return self._create_fallback_result(
                "anthropic package not installed", pmid=pmid
            )

        prompt = self._build_prompt(gene_name=gene_name, title=title, abstract=abstract)
        result = self._call_llm(prompt)

        if result:
            return RelevanceScore(
                is_relevant=result.get("is_relevant", False),
                confidence=float(result.get("confidence", 0.5)),
                reasoning=result.get("reasoning", ""),
                pmid=pmid,
            )
        else:
            return self._create_fallback_result(
                "Could not parse LLM response", pmid=pmid
            )

    def check_batch(
        self,
        gene_name: str,
        papers: List[Tuple[str, Optional[str], str]],
    ) -> List[RelevanceScore]:
        """Check relevance for multiple papers.

        Args:
            gene_name: The gene being searched for
            papers: List of (title, abstract, pmid) tuples

        Returns:
            List of RelevanceScore objects
        """
        return [
            self.check_relevance(gene_name, title, abstract, pmid)
            for title, abstract, pmid in papers
        ]


# =============================================================================
# Synonym Relevance Checker
# =============================================================================


class SynonymRelevanceChecker(BaseLLMRelevanceChecker):
    """Check if synonyms are actually relevant gene names using LLM assessment."""

    def _build_prompt(self, gene_name: str, synonym: str, source: str, **kwargs) -> str:
        """Build the synonym relevance prompt."""
        return f"""You are evaluating whether a gene synonym is actually relevant and useful for literature searches.

Primary gene name: {gene_name}
Synonym candidate: {synonym}
Source: {source}

Your task is to determine if this synonym is:
1. A legitimate alternate name for the gene (not just a misspelling or overly generic term)
2. Likely to appear in scientific literature when referring to this gene
3. Sufficiently specific to avoid excessive false positives

Consider:
- Official gene symbols and common aliases are usually highly relevant
- Very generic terms (like "protein", "factor", "component") are NOT relevant
- Overly verbose descriptions are less useful than concise terms
- Terms that could easily refer to other things are less relevant
- Historical names that are rarely used anymore may have lower relevance

Examples of RELEVANT synonyms:
- BRCA1 → BRCC1 (common historical alias)
- TP53 → P53 (widely used alternative)
- SCN5A → LQT3 (disease-associated designation)

Examples of NOT RELEVANT synonyms:
- BRCA1 → "breast cancer 1, early onset" (too verbose)
- CAT → "catalase" (could match too many unrelated terms)
- GENE → "protein" (too generic)

Respond with ONLY a JSON object with this exact format:
{{
  "is_relevant": true or false,
  "confidence": 0.0 to 1.0,
  "reasoning": "brief explanation (max 50 words)"
}}

Be reasonably selective - the goal is to expand searches with useful terms while avoiding noise."""

    def _create_fallback_result(
        self, reason: str, synonym: str = "", **kwargs
    ) -> SynonymRelevance:
        """Create a fallback result assuming relevance."""
        return SynonymRelevance(
            synonym=synonym,
            is_relevant=True,
            confidence=0.5,
            reasoning=reason,
        )

    def check_synonym_relevance(
        self,
        gene_name: str,
        synonym: str,
        source: str,
    ) -> SynonymRelevance:
        """Check if a synonym is relevant for the gene.

        Args:
            gene_name: The primary gene name
            synonym: The synonym candidate to evaluate
            source: Source of the synonym (e.g., 'official_symbol', 'alias')

        Returns:
            SynonymRelevance with assessment
        """
        if not self.api_key:
            return self._create_fallback_result(
                "No API key provided; cannot assess relevance", synonym=synonym
            )

        if self._get_client() is None:
            return self._create_fallback_result(
                "anthropic package not installed", synonym=synonym
            )

        prompt = self._build_prompt(gene_name=gene_name, synonym=synonym, source=source)
        result = self._call_llm(prompt)

        if result:
            return SynonymRelevance(
                synonym=synonym,
                is_relevant=result.get("is_relevant", False),
                confidence=float(result.get("confidence", 0.5)),
                reasoning=result.get("reasoning", ""),
            )
        else:
            return self._create_fallback_result(
                "Could not parse LLM response", synonym=synonym
            )

    def check_synonyms_batch(
        self,
        gene_name: str,
        synonyms: List[Tuple[str, str]],
    ) -> List[SynonymRelevance]:
        """Check relevance for multiple synonyms.

        Args:
            gene_name: The primary gene name
            synonyms: List of (synonym, source) tuples

        Returns:
            List of SynonymRelevance objects
        """
        return [
            self.check_synonym_relevance(gene_name, synonym, source)
            for synonym, source in synonyms
        ]
