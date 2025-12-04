"""
PubMed/NCBI utilities for querying and fetching paper metadata.

This module provides unified access to PubMed and NCBI Entrez APIs,
consolidating different query approaches into a single consistent interface.
"""

import logging
from typing import Set, List, Dict, Any, Optional
from Bio import Entrez
from Bio.Entrez.Parser import ValidationError
import requests

from config.settings import get_settings
from .retry_utils import api_retry

logger = logging.getLogger(__name__)

# Exceptions that should trigger a retry
_RETRYABLE_EXCEPTIONS = (
    requests.exceptions.RequestException,
    requests.exceptions.Timeout,
    ConnectionError,
    TimeoutError,
)


def _set_entrez_email(custom_email: Optional[str] = None) -> str:
    """Set ``Entrez.email`` to the provided value or the configured default."""

    email = custom_email or get_settings().ncbi_email
    Entrez.email = email
    return email


@api_retry
def query_pubmed_with_entrez(
    query: str,
    max_results: int = 100,
    sort: str = "relevance",
    email: Optional[str] = None,
) -> List[str]:
    """
    Query PubMed using Bio.Entrez and return PMIDs.

    This is the recommended method for querying PubMed as it uses the
    official Biopython Entrez interface.

    Args:
        query: PubMed search query (e.g., "BRCA1[Gene Symbol]")
        max_results: Maximum number of results to return (default: 100)
        sort: Sort order - "relevance", "pub_date", or "first_author" (default: "relevance")

    Returns:
        List of PMID strings

    Example:
        >>> pmids = query_pubmed_with_entrez("BRCA1[Gene Symbol]", max_results=50)
        >>> print(f"Found {len(pmids)} papers")
    """
    logger.info(f"Querying PubMed with query: {query}")

    try:
        _set_entrez_email(email)
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            sort=sort
        )
        record = Entrez.read(handle)
        handle.close()

        pmids = record.get("IdList", [])
        logger.info(f"PubMed query returned {len(pmids)} PMIDs")

        return pmids

    except _RETRYABLE_EXCEPTIONS as exc:
        logger.warning(
            "Transient error during PubMed query for '%s'; raising for retry", query, exc_info=exc
        )
        raise
    except (ValidationError, ValueError) as exc:
        logger.error("PubMed query returned malformed data for '%s': %s", query, exc)
        return []


def query_pubmed_for_gene(
    gene_symbol: str,
    max_results: int = 100,
    include_abstract: bool = True,
    email: Optional[str] = None,
) -> Set[str]:
    """
    Query PubMed for papers mentioning a specific gene.

    This function constructs a gene-specific query that searches both
    the Gene Symbol field and Title/Abstract.

    Args:
        gene_symbol: Gene symbol to search for (e.g., "BRCA1")
        max_results: Maximum number of results to return (default: 100)
        include_abstract: Whether to search in abstracts (default: True)

    Returns:
        Set of PMID strings

    Example:
        >>> pmids = query_pubmed_for_gene("BRCA1", max_results=50)
        >>> print(pmids)
    """
    logger.info(f"Querying PubMed for gene symbol: {gene_symbol}")

    # Build query
    if include_abstract:
        query = f"{gene_symbol}[Gene Symbol] OR {gene_symbol}[Title/Abstract]"
    else:
        query = f"{gene_symbol}[Gene Symbol]"

    pmids = query_pubmed_with_entrez(query, max_results=max_results, email=email)
    return set(pmids)


@api_retry
def fetch_paper_metadata(pmid: str, email: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """
    Fetch metadata for a single paper from PubMed.

    Uses the Entrez ESummary API to retrieve comprehensive metadata
    including title, authors, journal, publication date, etc.

    Args:
        pmid: PubMed ID (as string)

    Returns:
        Dictionary containing paper metadata, or None if fetch fails

    Example:
        >>> metadata = fetch_paper_metadata("12345678")
        >>> print(metadata.get("Title"))
    """
    logger.debug(f"Fetching metadata for PMID: {pmid}")

    try:
        _set_entrez_email(email)
        handle = Entrez.esummary(db="pubmed", id=pmid)
        record = Entrez.read(handle)
        handle.close()

        if not record:
            logger.warning(f"No metadata found for PMID: {pmid}")
            return None

        # ESummary returns a list with one document
        metadata = record[0] if isinstance(record, list) else record

        logger.debug(f"Successfully fetched metadata for PMID: {pmid}")
        return metadata

    except _RETRYABLE_EXCEPTIONS as exc:
        logger.warning(
            "Transient error fetching metadata for PMID %s; raising for retry", pmid, exc_info=exc
        )
        raise
    except (ValidationError, ValueError, IndexError, KeyError) as exc:
        logger.error("Malformed metadata response for PMID %s: %s", pmid, exc)
        return None


@api_retry
def fetch_paper_abstract(pmid: str, email: Optional[str] = None) -> Optional[str]:
    """
    Fetch the abstract text for a paper from PubMed.

    Uses the Entrez EFetch API to retrieve the full abstract text.

    Args:
        pmid: PubMed ID (as string)

    Returns:
        Abstract text, or None if not available

    Example:
        >>> abstract = fetch_paper_abstract("12345678")
        >>> print(abstract[:100])
    """
    logger.debug(f"Fetching abstract for PMID: {pmid}")

    try:
        _set_entrez_email(email)
        handle = Entrez.efetch(
            db="pubmed",
            id=pmid,
            rettype="abstract",
            retmode="text"
        )
        abstract = handle.read()
        handle.close()

        if abstract:
            logger.debug(f"Successfully fetched abstract for PMID: {pmid}")
            return abstract
        else:
            logger.warning(f"No abstract found for PMID: {pmid}")
            return None

    except _RETRYABLE_EXCEPTIONS as exc:
        logger.warning(
            "Transient error fetching abstract for PMID %s; raising for retry", pmid, exc_info=exc
        )
        raise
    except (ValidationError, ValueError) as exc:
        logger.error("Malformed abstract response for PMID %s: %s", pmid, exc)
        return None


@api_retry
def get_doi_from_pmid(pmid: str, email: Optional[str] = None) -> Optional[str]:
    """
    Get the DOI for a paper given its PMID.

    Args:
        pmid: PubMed ID (as string)

    Returns:
        DOI string, or None if not found

    Example:
        >>> doi = get_doi_from_pmid("12345678")
        >>> print(doi)
        '10.1038/nature12345'
    """
    logger.debug(f"Fetching DOI for PMID: {pmid}")

    try:
        _set_entrez_email(email)
        # Fetch the full record in XML format
        handle = Entrez.efetch(
            db="pubmed",
            id=pmid,
            rettype="xml",
            retmode="xml"
        )
        records = Entrez.read(handle)
        handle.close()

        # Navigate the XML structure to find DOI
        if records and "PubmedArticle" in records:
            for article in records["PubmedArticle"]:
                article_ids = article.get("PubmedData", {}).get("ArticleIdList", [])
                for article_id in article_ids:
                    if hasattr(article_id, 'attributes') and \
                       article_id.attributes.get("IdType") == "doi":
                        doi = str(article_id)
                        logger.debug(f"Found DOI {doi} for PMID {pmid}")
                        return doi

        logger.warning(f"No DOI found for PMID: {pmid}")
        return None

    except _RETRYABLE_EXCEPTIONS as exc:
        logger.warning(
            "Transient error fetching DOI for PMID %s; raising for retry", pmid, exc_info=exc
        )
        raise
    except (ValidationError, ValueError, KeyError, IndexError) as exc:
        logger.error("Malformed DOI response for PMID %s: %s", pmid, exc)
        return None


@api_retry
def query_europepmc(
    gene_symbol: str,
    max_results: int = 100
) -> Set[str]:
    """
    Query Europe PMC for papers mentioning a gene symbol.

    Europe PMC is an alternative to PubMed with broader coverage of
    European biomedical literature.

    Args:
        gene_symbol: Gene symbol to search for
        max_results: Maximum number of results to return (default: 100)

    Returns:
        Set of PMID strings

    Example:
        >>> pmids = query_europepmc("BRCA1", max_results=50)
        >>> print(pmids)
    """
    logger.info(f"Querying Europe PMC for gene symbol: {gene_symbol}")

    search_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": f"{gene_symbol}",
        "format": "json",
        "pageSize": max_results,
        "resultType": "core"
    }

    try:
        response = requests.get(search_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        pmids = set()
        results = data.get("resultList", {}).get("result", [])

        for result in results:
            # Europe PMC results may have PMID or other IDs
            if "pmid" in result and result["pmid"]:
                pmids.add(str(result["pmid"]))

        logger.info(f"Europe PMC returned {len(pmids)} PMIDs for {gene_symbol}")
        return pmids

    except _RETRYABLE_EXCEPTIONS as exc:
        logger.warning(
            "Transient error querying Europe PMC for %s; raising for retry", gene_symbol, exc_info=exc
        )
        raise
    except ValueError as exc:
        logger.error("Europe PMC returned malformed data for %s: %s", gene_symbol, exc)
        return set()


def batch_fetch_metadata(
    pmids: List[str],
    batch_size: int = 200,
    email: Optional[str] = None
) -> Dict[str, Dict[str, Any]]:
    """
    Fetch metadata for multiple papers in batches.

    This is more efficient than fetching metadata one paper at a time.

    Args:
        pmids: List of PMIDs
        batch_size: Number of PMIDs to fetch per batch (default: 200)

    Returns:
        Dictionary mapping PMIDs to their metadata

    Example:
        >>> pmids = ["12345678", "87654321"]
        >>> metadata_dict = batch_fetch_metadata(pmids)
        >>> print(metadata_dict["12345678"]["Title"])
    """
    logger.info(f"Fetching metadata for {len(pmids)} papers in batches of {batch_size}")

    metadata_dict = {}

    # Process in batches
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        batch_ids = ",".join(batch)

        try:
            _set_entrez_email(email)
            handle = Entrez.esummary(db="pubmed", id=batch_ids)
            records = Entrez.read(handle)
            handle.close()

            # Parse the batch results
            for record in records:
                if "Id" in record:
                    pmid = str(record["Id"])
                    metadata_dict[pmid] = record

            logger.debug(f"Fetched metadata for batch {i // batch_size + 1}")

        except Exception as e:
            logger.error(f"Failed to fetch metadata for batch starting at {i}: {e}")
            continue

    logger.info(f"Successfully fetched metadata for {len(metadata_dict)} papers")
    return metadata_dict


def validate_pmid(pmid: str) -> bool:
    """
    Validate that a PMID exists in PubMed.

    Args:
        pmid: PubMed ID to validate

    Returns:
        True if the PMID exists in PubMed

    Example:
        >>> is_valid = validate_pmid("12345678")
        >>> print(is_valid)
    """
    metadata = fetch_paper_metadata(pmid)
    return metadata is not None
