"""
Paper sourcing module - queries PubMed, EuropePMC, and PubMind for papers.
"""

import logging
import requests
from typing import List, Set, Optional
from tenacity import retry, stop_after_attempt, wait_exponential
from models import Paper

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
    Queries PubMed and EuropePMC APIs to find papers related to a gene symbol.
    Returns deduplicated list of PMIDs.
    """

    def __init__(self, email: str = "your_email@example.com"):
        """
        Initialize the sourcer.

        Args:
            email: Email for PubMed API (required by NCBI).
        """
        self.email = email
        self.pubmed_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.europepmc_base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest"

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
    def _query_pubmed(self, gene_symbol: str, max_results: int = 100) -> Set[str]:
        """
        Query PubMed (via Entrez) for papers mentioning the gene symbol.

        Args:
            gene_symbol: Gene symbol to search for.
            max_results: Maximum number of results to return.

        Returns:
            Set of PMIDs.
        """
        logger.info(f"Querying PubMed for gene symbol: {gene_symbol}")

        # Step 1: Search for papers
        search_url = f"{self.pubmed_base_url}/esearch.fcgi"
        search_params = {
            "db": "pubmed",
            "term": f"{gene_symbol}[Gene Symbol] OR {gene_symbol}[Title/Abstract]",
            "retmax": max_results,
            "retmode": "json",
            "email": self.email,
            "tool": "BiomedicalExtractionPipeline"
        }

        try:
            response = requests.get(search_url, params=search_params, timeout=30)
            response.raise_for_status()
            data = response.json()

            pmids = data.get("esearchresult", {}).get("idlist", [])
            logger.info(f"PubMed returned {len(pmids)} PMIDs for {gene_symbol}")

            return set(pmids)

        except requests.RequestException as e:
            logger.error(f"PubMed query failed: {e}")
            return set()

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
    def _query_europepmc(self, gene_symbol: str, max_results: int = 100) -> Set[str]:
        """
        Query EuropePMC for papers mentioning the gene symbol.

        Args:
            gene_symbol: Gene symbol to search for.
            max_results: Maximum number of results to return.

        Returns:
            Set of PMIDs.
        """
        logger.info(f"Querying EuropePMC for gene symbol: {gene_symbol}")

        search_url = f"{self.europepmc_base_url}/search"
        search_params = {
            "query": f"(GENE:{gene_symbol}) OR (ABSTRACT:{gene_symbol})",
            "format": "json",
            "pageSize": max_results,
            "resultType": "core"
        }

        try:
            response = requests.get(search_url, params=search_params, timeout=30)
            response.raise_for_status()
            data = response.json()

            # Extract PMIDs from results
            pmids = set()
            for result in data.get("resultList", {}).get("result", []):
                pmid = result.get("pmid")
                if pmid:
                    pmids.add(str(pmid))

            logger.info(f"EuropePMC returned {len(pmids)} PMIDs for {gene_symbol}")
            return pmids

        except requests.RequestException as e:
            logger.error(f"EuropePMC query failed: {e}")
            return set()

    def fetch_papers(
        self,
        gene_symbol: str,
        max_results_per_source: int = 100,
        use_pubmed: bool = True,
        use_europepmc: bool = True,
        use_pubmind: bool = True
    ) -> List[str]:
        """
        Fetch deduplicated list of PMIDs from multiple sources.

        Args:
            gene_symbol: Gene symbol to search for.
            max_results_per_source: Maximum results per API source.
            use_pubmed: Query PubMed API.
            use_europepmc: Query EuropePMC API.
            use_pubmind: Query PubMind database (if available).

        Returns:
            Deduplicated sorted list of PMIDs.
        """
        all_pmids: Set[str] = set()

        # Query PubMind first (most relevant for variant-level data)
        if use_pubmind and PUBMIND_AVAILABLE:
            try:
                logger.info(f"Querying PubMind for {gene_symbol}...")
                pubmind_fetcher = PubMindFetcher(email=self.email, use_fallback=False)
                pubmind_pmids = pubmind_fetcher.fetch_pmids_for_gene(
                    gene_symbol,
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

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
    def fetch_paper_metadata(self, pmid: str) -> Optional[Paper]:
        """
        Fetch metadata for a single paper from PubMed.

        Args:
            pmid: PubMed ID.

        Returns:
            Paper object with metadata, or None if fetch failed.
        """
        logger.debug(f"Fetching metadata for PMID: {pmid}")

        fetch_url = f"{self.pubmed_base_url}/esummary.fcgi"
        params = {
            "db": "pubmed",
            "id": pmid,
            "retmode": "json",
            "email": self.email
        }

        try:
            response = requests.get(fetch_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            result = data.get("result", {}).get(pmid, {})

            if not result:
                logger.warning(f"No metadata found for PMID: {pmid}")
                return None

            # Parse metadata
            paper = Paper(
                pmid=pmid,
                title=result.get("title"),
                authors=[author.get("name") for author in result.get("authors", [])],
                journal=result.get("fulljournalname"),
                publication_date=result.get("pubdate"),
                doi=next((id_dict["value"] for id_dict in result.get("articleids", [])
                         if id_dict.get("idtype") == "doi"), None),
                pmc_id=next((id_dict["value"] for id_dict in result.get("articleids", [])
                            if id_dict.get("idtype") == "pmc"), None),
                source="PubMed"
            )

            return paper

        except requests.RequestException as e:
            logger.error(f"Failed to fetch metadata for PMID {pmid}: {e}")
            return None


def query_papers_for_gene(gene_symbol: str, email: str = "your_email@example.com") -> List[str]:
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
