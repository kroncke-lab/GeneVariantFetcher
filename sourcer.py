"""
Paper sourcing module - queries PubMed and EuropePMC APIs for papers.
"""

import logging
import re
from typing import List, Set, Optional, Iterable, Any

import requests
from bs4 import BeautifulSoup
from tenacity import retry, stop_after_attempt, wait_exponential
from models import Paper

logger = logging.getLogger(__name__)


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
        use_pubmind: bool = False,
        pubmind_query: Optional[str] = None,
    ) -> List[str]:
        """
        Fetch deduplicated list of PMIDs from multiple sources.

        Args:
            gene_symbol: Gene symbol to search for.
            max_results_per_source: Maximum results per API source.
            use_pubmed: Query PubMed API.
            use_europepmc: Query EuropePMC API.
            use_pubmind: Query PubMind (via API or HTML scrape) for variants.
            pubmind_query: Optional override for the PubMind search term.

        Returns:
            Deduplicated sorted list of PMIDs.
        """
        all_pmids: Set[str] = set()

        if use_pubmed:
            pubmed_pmids = self._query_pubmed(gene_symbol, max_results_per_source)
            all_pmids.update(pubmed_pmids)

        if use_europepmc:
            europepmc_pmids = self._query_europepmc(gene_symbol, max_results_per_source)
            all_pmids.update(europepmc_pmids)

        if use_pubmind:
            pubmind_sourcer = PubMindSourcer()
            pubmind_pmids = pubmind_sourcer.search(
                pubmind_query or gene_symbol, max_results=max_results_per_source
            )
            all_pmids.update(pubmind_pmids)

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


class PubMindSourcer:
    """Utility to source PMIDs from PubMind (https://pubmind.wglab.org/).

    The client first attempts the documented JSON API endpoint and gracefully
    falls back to scraping the HTML search results when the API is unavailable
    (common in restricted network environments).
    """

    api_search_url = "https://pubmind.wglab.org/api/search"
    html_search_url = "https://pubmind.wglab.org/search"

    def __init__(self, session: Optional[requests.Session] = None):
        self.session = session or requests.Session()

    def search(self, query: str, max_results: int = 200) -> List[str]:
        """Return a list of PMIDs relevant to the query.

        Args:
            query: Gene symbol, variant, or free text for PubMind search.
            max_results: Maximum records to inspect per source.
        """

        pmids: Set[str] = set()

        # Try JSON API first
        try:
            response = self.session.get(
                self.api_search_url,
                params={"query": query, "size": max_results},
                timeout=30,
            )
            response.raise_for_status()

            if "application/json" in response.headers.get("content-type", ""):
                data = response.json()
                pmids.update(self._extract_pmids_from_json(data))
        except Exception as exc:  # noqa: BLE001
            logger.warning("PubMind API lookup failed, falling back to HTML scrape: %s", exc)

        # Fallback to HTML search page scraping
        if not pmids:
            try:
                html_response = self.session.get(
                    self.html_search_url,
                    params={"keyword": query},
                    timeout=30,
                )
                html_response.raise_for_status()
                pmids.update(self._extract_pmids_from_html(html_response.text))
            except requests.RequestException as exc:
                logger.error("PubMind HTML scrape failed: %s", exc)

        return sorted(pmids)

    def _extract_pmids_from_json(self, payload: Any) -> Set[str]:
        pmids: Set[str] = set()

        def walk(node: Any):
            if isinstance(node, dict):
                for key, value in node.items():
                    if key.lower() == "pmid" and isinstance(value, (str, int)):
                        pmids.add(str(value))
                    else:
                        walk(value)
            elif isinstance(node, list):
                for item in node:
                    walk(item)

        walk(payload)
        return pmids

    def _extract_pmids_from_html(self, html: str) -> Set[str]:
        pmids: Set[str] = set()
        soup = BeautifulSoup(html, "html.parser")

        def _iter_texts(nodes: Iterable) -> Iterable[str]:
            for node in nodes:
                text = node.get_text(" ", strip=True)
                if text:
                    yield text
                href = node.get("href")
                if href:
                    yield href

        for text in _iter_texts(soup.find_all("a")):
            pmids.update(re.findall(r"PMID[:\s]*([0-9]{4,10})", text, flags=re.IGNORECASE))
            pmids.update(re.findall(r"(?<!\d)([0-9]{7,10})(?!\d)", text))

        # Some pages display PMIDs in table cells or plain text
        if not pmids:
            page_text = soup.get_text(" ", strip=True)
            pmids.update(re.findall(r"PMID[:\s]*([0-9]{4,10})", page_text, flags=re.IGNORECASE))

        return pmids
