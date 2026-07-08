"""
Automatic synonym discovery for genes using NCBI Gene database.

This module provides functionality to:
1. Query NCBI Gene database for official gene synonyms and aliases
2. Use LLM to assess relevance of synonyms
3. Present synonyms to users for interactive selection
4. Return selected synonyms for use in PubMed searches
"""

import logging
import time
from dataclasses import dataclass
from typing import List, Optional, Set

import requests
from requests.adapters import HTTPAdapter, Retry

logger = logging.getLogger(__name__)

try:
    from gene_literature.llm_relevance import SynonymRelevanceChecker

    RELEVANCE_CHECKER_AVAILABLE = True
except ImportError:
    RELEVANCE_CHECKER_AVAILABLE = False


class SynonymFinderError(Exception):
    """Raised when synonym lookup fails."""


@dataclass
class GeneSynonym:
    """Represents a gene synonym with metadata."""

    term: str
    source: str  # e.g., "official_symbol", "alias", "other_designations"
    gene_id: Optional[int] = None
    is_relevant: Optional[bool] = None  # LLM-assessed relevance
    relevance_confidence: Optional[float] = None  # 0.0 to 1.0
    relevance_reasoning: Optional[str] = None


class SynonymFinder:
    """
    Finds gene synonyms using NCBI Gene database.

    Uses the NCBI E-utilities API to search for genes and retrieve their
    official symbols, aliases, and other designations.
    """

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(
        self,
        email: Optional[str] = None,
        api_key: Optional[str] = None,
        retry_attempts: int = 3,
        anthropic_api_key: Optional[str] = None,
    ):
        """
        Initialize the SynonymFinder.

        Args:
            email: Email address for NCBI API (recommended)
            api_key: NCBI API key for higher rate limits (optional)
            retry_attempts: Number of retry attempts for failed requests
            anthropic_api_key: Anthropic API key for LLM-based relevance checking (optional)
        """
        self.email = email
        self.api_key = api_key
        self.retry_attempts = retry_attempts
        self.anthropic_api_key = anthropic_api_key

        # Configure session with retry logic
        self.session = requests.Session()
        retry_strategy = Retry(
            total=retry_attempts,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    def find_gene_synonyms(
        self,
        gene: str,
        include_other_designations: bool = True,
    ) -> List[GeneSynonym]:
        """
        Find synonyms for a given gene.

        Args:
            gene: Gene name or symbol to search for
            include_other_designations: Include other gene designations

        Returns:
            List of GeneSynonym objects

        Raises:
            SynonymFinderError: If the API request fails
        """
        logger.info("Searching for synonyms of '%s'", gene)

        # Step 1: Search for the gene to get Gene ID
        gene_id = self._search_gene(gene)
        if gene_id is None:
            logger.warning("No gene found for '%s'", gene)
            return []

        logger.info("Found Gene ID: %s", gene_id)

        # Step 2: Fetch gene summary to get synonyms
        synonyms = self._fetch_gene_summary(gene_id, include_other_designations)

        logger.info("Found %d synonyms for '%s'", len(synonyms), gene)

        # Step 3: Check relevance with LLM if API key provided
        if self.anthropic_api_key and RELEVANCE_CHECKER_AVAILABLE:
            logger.info("Checking synonym relevance with LLM...")
            synonyms = self._check_synonyms_relevance(gene, synonyms)
        elif self.anthropic_api_key and not RELEVANCE_CHECKER_AVAILABLE:
            logger.warning(
                "Anthropic API key provided but anthropic package not installed. "
                "Install with: pip install anthropic"
            )

        return synonyms

    def _check_synonyms_relevance(
        self,
        gene: str,
        synonyms: List[GeneSynonym],
    ) -> List[GeneSynonym]:
        """
        Check relevance of synonyms using LLM.

        Args:
            gene: Primary gene name
            synonyms: List of synonym candidates

        Returns:
            List of GeneSynonym objects with relevance information added
        """
        checker = SynonymRelevanceChecker(api_key=self.anthropic_api_key)

        # Check each synonym
        for synonym in synonyms:
            relevance = checker.check_synonym_relevance(
                gene_name=gene,
                synonym=synonym.term,
                source=synonym.source,
            )

            # Update synonym with relevance info
            synonym.is_relevant = relevance.is_relevant
            synonym.relevance_confidence = relevance.confidence
            synonym.relevance_reasoning = relevance.reasoning

            logger.debug(
                "Synonym '%s': %s (confidence: %.2f) - %s",
                synonym.term,
                "RELEVANT" if relevance.is_relevant else "NOT RELEVANT",
                relevance.confidence,
                relevance.reasoning,
            )

        # Sort synonyms: relevant first, then by confidence
        synonyms_sorted = sorted(
            synonyms,
            key=lambda s: (
                not s.is_relevant if s.is_relevant is not None else False,
                -(s.relevance_confidence or 0.0),
            ),
        )

        relevant_count = sum(1 for s in synonyms if s.is_relevant)
        logger.info(
            "LLM assessed %d/%d synonyms as relevant",
            relevant_count,
            len(synonyms),
        )

        return synonyms_sorted

    def _search_gene(self, gene: str) -> Optional[int]:
        """
        Search for a gene in NCBI Gene database.

        Args:
            gene: Gene name or symbol

        Returns:
            Gene ID if found, None otherwise
        """
        params = {
            "db": "gene",
            "term": f"{gene}[Gene Name] AND human[Organism]",
            "retmode": "json",
            "retmax": 1,  # Only need the top result
        }

        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key

        url = f"{self.BASE_URL}/esearch.fcgi"

        try:
            response = self._request(url, params)
            data = response.json()

            id_list = data.get("esearchresult", {}).get("idlist", [])
            if id_list:
                return int(id_list[0])
            return None

        except Exception as e:
            logger.error("Failed to search for gene '%s': %s", gene, e)
            raise SynonymFinderError(f"Gene search failed: {e}") from e

    def _fetch_gene_summary(
        self,
        gene_id: int,
        include_other_designations: bool,
    ) -> List[GeneSynonym]:
        """
        Fetch gene summary including synonyms.

        Args:
            gene_id: NCBI Gene ID
            include_other_designations: Include other gene designations

        Returns:
            List of GeneSynonym objects
        """
        params = {
            "db": "gene",
            "id": str(gene_id),
            "retmode": "json",
        }

        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key

        url = f"{self.BASE_URL}/esummary.fcgi"

        try:
            response = self._request(url, params)
            data = response.json()

            result = data.get("result", {}).get(str(gene_id), {})
            synonyms = []

            # Official symbol
            official_symbol = result.get("name")
            if official_symbol:
                synonyms.append(
                    GeneSynonym(
                        term=official_symbol,
                        source="official_symbol",
                        gene_id=gene_id,
                    )
                )

            # Aliases (common synonyms)
            aliases = result.get("otheraliases", "").split(", ")
            for alias in aliases:
                alias = alias.strip()
                if alias:
                    synonyms.append(
                        GeneSynonym(
                            term=alias,
                            source="alias",
                            gene_id=gene_id,
                        )
                    )

            # Other designations (if requested)
            if include_other_designations:
                other_designations = result.get("otherdesignations", "").split("|")
                for designation in other_designations:
                    designation = designation.strip()
                    if designation:
                        synonyms.append(
                            GeneSynonym(
                                term=designation,
                                source="other_designation",
                                gene_id=gene_id,
                            )
                        )

            return synonyms

        except Exception as e:
            logger.error("Failed to fetch gene summary for ID %s: %s", gene_id, e)
            raise SynonymFinderError(f"Gene summary fetch failed: {e}") from e

    def _request(self, url: str, params: dict) -> requests.Response:
        """
        Make an HTTP request with retry logic.

        Args:
            url: URL to request
            params: Query parameters

        Returns:
            Response object

        Raises:
            SynonymFinderError: If request fails after retries
        """
        for attempt in range(self.retry_attempts):
            try:
                # Add delay to respect NCBI rate limits (3 requests/sec without API key)
                if attempt > 0:
                    delay = 2**attempt  # Exponential backoff
                    logger.debug("Retrying request after %ds delay...", delay)
                    time.sleep(delay)
                else:
                    # Small delay between requests
                    time.sleep(0.34)

                response = self.session.get(url, params=params, timeout=30)
                response.raise_for_status()
                return response

            except requests.RequestException as e:
                logger.warning("Request attempt %d failed: %s", attempt + 1, e)
                if attempt == self.retry_attempts - 1:
                    raise SynonymFinderError(
                        f"Request failed after {self.retry_attempts} attempts"
                    ) from e

        # Should not reach here, but just in case
        raise SynonymFinderError("Request failed")


def automatic_synonym_selection(
    gene: str,
    synonyms: List[GeneSynonym],
    include_official: bool = True,
    include_aliases: bool = True,
    include_other_designations: bool = False,
    only_relevant: bool = True,
    min_confidence: float = 0.7,
) -> List[str]:
    """
    Automatically select synonyms for batch/non-interactive mode.

    Args:
        gene: Original gene name
        synonyms: List of found synonyms
        include_official: Include official gene symbols
        include_aliases: Include gene aliases
        include_other_designations: Include verbose other designations
        only_relevant: Only include LLM-assessed relevant synonyms (if available)
        min_confidence: Minimum confidence threshold for relevance

    Returns:
        List of selected synonym terms
    """
    if not synonyms:
        logger.info("No synonyms found for '%s'", gene)
        return []

    selected: Set[str] = set()
    has_relevance_info = any(s.is_relevant is not None for s in synonyms)

    for syn in synonyms:
        # Skip if relevance filtering is enabled and synonym is not relevant
        if only_relevant and has_relevance_info:
            if syn.is_relevant is False:
                continue
            if (
                syn.relevance_confidence is not None
                and syn.relevance_confidence < min_confidence
            ):
                continue

        # Filter by source type
        if syn.source == "official_symbol" and include_official:
            selected.add(syn.term)
        elif syn.source == "alias" and include_aliases:
            selected.add(syn.term)
        elif syn.source == "other_designation" and include_other_designations:
            selected.add(syn.term)

    # Remove the primary gene if it's in the list (avoid duplication)
    selected.discard(gene)
    selected.discard(gene.upper())
    selected.discard(gene.lower())

    result = sorted(selected)
    logger.info(
        "Automatically selected %d synonyms for '%s': %s",
        len(result),
        gene,
        ", ".join(result) if result else "(none)",
    )
    return result
