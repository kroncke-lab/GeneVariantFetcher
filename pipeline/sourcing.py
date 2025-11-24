"""
Paper sourcing module - queries PubMind (primary), PubMed, and EuropePMC for papers.
PubMind is prioritized as it provides more relevant variant-level data.
"""

import logging
from typing import List, Set, Optional

from models import Paper
from utils.pubmed_utils import query_pubmed_for_gene, query_europepmc
from pipeline.utils.retry_utils import standard_retry
from config.settings import get_settings

logger = logging.getLogger(__name__)

# Import PubMindFetcher if available
try:
    from pubmind_fetcher import PubMindFetcher
    PUBMIND_AVAILABLE = True
except ImportError:
    PUBMIND_AVAILABLE = False
    logger.warning("PubMindFetcher not available. Install dependencies or check pubmind_fetcher.py")


class PaperSourcer:
    """
    Queries PubMind (primary source), PubMed, and EuropePMC APIs to find papers related to a gene symbol.
    PubMind is prioritized as it provides more relevant variant-level data and reduces false positives.
    Returns deduplicated list of PMIDs.
    """

    def __init__(self, email: Optional[str] = None):
        """
        Initialize the sourcer.

        Args:
            email: Email for PubMed API (required by NCBI). If None, uses config settings.
        """
        settings = get_settings()
        self.email = email or settings.ncbi_email

    def _query_pubmed(self, gene_symbol: str, max_results: int = 100) -> Set[str]:
        """
        Query PubMed (via Entrez) for papers mentioning the gene symbol.

        Args:
            gene_symbol: Gene symbol to search for.
            max_results: Maximum number of results to return.

        Returns:
            Set of PMIDs.
        """
        # Use shared utility function
        return query_pubmed_for_gene(
            gene_symbol,
            max_results=max_results,
            email=self.email,
        )

    def _query_europepmc(self, gene_symbol: str, max_results: int = 100) -> Set[str]:
        """
        Query EuropePMC for papers mentioning the gene symbol.

        Args:
            gene_symbol: Gene symbol to search for.
            max_results: Maximum number of results to return.

        Returns:
            Set of PMIDs.
        """
        # Use shared utility function
        return query_europepmc(gene_symbol, max_results=max_results)

    def fetch_papers(
        self,
        gene_symbol: str,
        max_results_per_source: int = 100,
        use_pubmed: bool = True,
        use_europepmc: bool = True,
        use_pubmind: bool = True,
        pubmind_query: Optional[str] = None,
    ) -> List[str]:
        """
        Fetch deduplicated list of PMIDs from multiple sources.

        Args:
            gene_symbol: Gene symbol to search for.
            max_results_per_source: Maximum results per API source.
            use_pubmed: Query PubMed API.
            use_europepmc: Query EuropePMC API.
            use_pubmind: Query PubMind (via PubMindFetcher) for variants.
            pubmind_query: Optional override for the PubMind search term (defaults to gene_symbol).

        Returns:
            Deduplicated sorted list of PMIDs.
        """
        all_pmids: Set[str] = set()

        # Query PubMind first (most relevant for variant-level data)
        if use_pubmind and PUBMIND_AVAILABLE:
            try:
                logger.info(f"Querying PubMind for {gene_symbol}...")
                pubmind_fetcher = PubMindFetcher(email=self.email)

                # Use custom query if provided, otherwise use gene_symbol
                search_term = pubmind_query or gene_symbol
                pubmind_pmids = pubmind_fetcher.fetch_pmids_for_gene(
                    search_term,
                    max_results=max_results_per_source
                )
                all_pmids.update(pubmind_pmids)
                logger.info(f"PubMind contributed {len(pubmind_pmids)} PMIDs")
            except Exception as e:
                logger.warning(f"PubMind query failed: {e}, continuing with other sources")

        if use_pubmed:
            pubmed_pmids = self._query_pubmed(gene_symbol, max_results_per_source)
            all_pmids.update(pubmed_pmids)

        if use_europepmc:
            europepmc_pmids = self._query_europepmc(gene_symbol, max_results_per_source)
            all_pmids.update(europepmc_pmids)

        deduplicated_pmids = sorted(list(all_pmids))
        logger.info(
            f"Total deduplicated PMIDs for {gene_symbol}: {len(deduplicated_pmids)}"
        )

        return deduplicated_pmids

    @standard_retry
    def fetch_paper_metadata(self, pmid: str) -> Optional[Paper]:
        """
        Fetch metadata for a single paper from PubMed.

        Args:
            pmid: PubMed ID.

        Returns:
            Paper object with metadata, or None if fetch failed.
        """
        from utils.pubmed_utils import fetch_paper_metadata

        logger.debug(f"Fetching metadata for PMID: {pmid}")

        try:
            # Use shared utility to fetch metadata
            metadata = fetch_paper_metadata(pmid, email=self.email)

            if not metadata:
                logger.warning(f"No metadata found for PMID: {pmid}")
                return None

            # Parse metadata into Paper object
            # ESummary returns different structure than the old API
            paper = Paper(
                pmid=pmid,
                title=metadata.get("Title", ""),
                authors=[author.get("Name", "") for author in metadata.get("AuthorList", [])],
                journal=metadata.get("FullJournalName", ""),
                publication_date=metadata.get("PubDate", ""),
                doi=metadata.get("DOI"),
                pmc_id=next((id_obj["Value"] for id_obj in metadata.get("ArticleIds", [])
                            if id_obj.get("IdType") == "pmc"), None),
                source="PubMed"
            )

            return paper

        except Exception as e:
            logger.error(f"Failed to fetch metadata for PMID {pmid}: {e}")
            return None


def query_papers_for_gene(gene_symbol: str, email: Optional[str] = None) -> List[str]:
    """
    Convenience function to query papers for a gene symbol.

    Args:
        gene_symbol: Gene symbol to search for.
        email: Email for PubMed API.

    Returns:
        Deduplicated list of PMIDs.
    """
    sourcer = PaperSourcer(email=email)
    return sourcer.fetch_papers(gene_symbol)