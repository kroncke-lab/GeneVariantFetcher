"""Utilities for combining PubMind, PubMed, and EuropePMC PMID discovery."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence

from config.settings import Settings, get_settings
from .collector import build_gene_query
from .pubmed_client import PubMedClient, PubMedError
from .pubmind_fetcher import PubMindFetcher
from utils.pubmed_utils import query_europepmc

logger = logging.getLogger(__name__)


@dataclass
class PMIDDiscoveryResult:
    """Structured representation of PMID discovery results."""

    pubmind_pmids: List[str]
    pubmed_pmids: List[str]
    europepmc_pmids: List[str]
    combined_pmids: List[str]
    synonyms_used: List[str] = field(default_factory=list)


def build_gene_keyword_queries(
    gene_symbol: str,
    synonyms: Optional[Sequence[str]] = None,
) -> List[str]:
    """Return a set of PubMed queries focused on a gene and variant keywords.

    Args:
        gene_symbol: Primary gene symbol to search for
        synonyms: Optional list of gene synonyms to include in queries

    Returns:
        List of PubMed query strings
    """
    base_query = build_gene_query(gene_symbol, synonyms)
    keyword_clauses: Sequence[str] = (
        "",  # base gene query
        "AND (variant OR mutation OR polymorphism)",
        "AND (pathogenic OR germline OR somatic OR hereditary)",
        "AND (case report OR cohort OR patient)",
    )

    queries: List[str] = []
    for clause in keyword_clauses:
        query = f"{base_query} {clause}".strip()
        if query not in queries:
            queries.append(query)

    synonym_info = f" (+ {len(synonyms)} synonyms)" if synonyms else ""
    logger.debug(
        "Built %d PubMed keyword queries for %s%s",
        len(queries),
        gene_symbol,
        synonym_info,
    )
    return queries


def discover_pmids_for_gene(
    gene_symbol: str,
    *,
    email: str | None = None,
    max_results: int | None = None,
    pubmind_output: Path | None = None,
    pubmed_output: Path | None = None,
    combined_output: Path | None = None,
    api_key: str | None = None,
    settings: Settings | None = None,
    synonyms: Sequence[str] | None = None,
) -> PMIDDiscoveryResult:
    """Fetch PMIDs from configured sources, merge, and persist them.

    Args:
        gene_symbol: Primary gene symbol to search for
        email: NCBI email for API access
        max_results: Maximum number of PMIDs to return per source
        pubmind_output: Path to save PubMind PMIDs
        pubmed_output: Path to save PubMed PMIDs
        combined_output: Path to save combined PMIDs
        api_key: NCBI API key for higher rate limits
        settings: Configuration settings
        synonyms: Optional list of gene synonyms to include in PubMed queries

    Returns:
        PMIDDiscoveryResult containing PMIDs from all sources
    """
    settings = settings or get_settings()
    effective_email = email or settings.ncbi_email
    if not effective_email:
        raise ValueError("NCBI email is required for PMID discovery.")

    # Convert synonyms to list for consistent handling
    synonyms_list = list(synonyms) if synonyms else []
    if synonyms_list:
        logger.info(
            "Using %d gene synonyms: %s", len(synonyms_list), ", ".join(synonyms_list)
        )

    effective_max_results = max_results or settings.max_papers_per_source

    # Resolve which sources are active
    use_pubmind = settings.use_pubmind
    use_pubmed = settings.use_pubmed
    use_europepmc = settings.use_europepmc

    if settings.pubmind_only:
        use_pubmed = False
        use_europepmc = False

    if not any([use_pubmind, use_pubmed, use_europepmc]):
        raise ValueError(
            "No literature sources are enabled. Set PUBMIND_ONLY=true to limit to PubMind or toggle USE_PUBMED/USE_EUROPEPMC to disable specific sources."
        )

    pmid_sets: list[set[str]] = []

    # PubMind discovery
    pubmind_pmids: list[str] = []
    if use_pubmind:
        pubmind_fetcher = PubMindFetcher(email=effective_email)
        pubmind_pmids = pubmind_fetcher.fetch_pmids_for_gene(
            gene_symbol, max_results=effective_max_results
        )
        pmid_sets.append(set(pubmind_pmids))
        if pubmind_output:
            pubmind_fetcher.save_pmids_to_file(pubmind_pmids, pubmind_output)

    # PubMed discovery using keyword-focused queries (with synonyms if provided)
    pubmed_pmids_list: list[str] = []
    if use_pubmed:
        pubmed_client = PubMedClient(api_key=api_key, email=effective_email)
        pubmed_pmids: set[str] = set()
        for query in build_gene_keyword_queries(
            gene_symbol, synonyms=synonyms_list or None
        ):
            try:
                pubmed_pmids.update(
                    pubmed_client.search(query, retmax=effective_max_results)
                )
            except PubMedError as exc:
                logger.warning("PubMed query failed for '%s': %s", query, exc)

        pubmed_pmids_list = sorted(pubmed_pmids)[:effective_max_results]
        pmid_sets.append(set(pubmed_pmids_list))
        if pubmed_output:
            _save_pmids(pubmed_pmids_list, pubmed_output)

    # EuropePMC discovery
    europepmc_pmids: list[str] = []
    if use_europepmc:
        europepmc_pmids = sorted(
            query_europepmc(gene_symbol, max_results=effective_max_results)
        )
        pmid_sets.append(set(europepmc_pmids))

    # Merge and de-duplicate
    combined_pmids = sorted(set().union(*pmid_sets)) if pmid_sets else []
    if len(combined_pmids) > effective_max_results:
        combined_pmids = combined_pmids[:effective_max_results]

    if combined_output:
        _save_pmids(combined_pmids, combined_output)

    return PMIDDiscoveryResult(
        pubmind_pmids=pubmind_pmids,
        pubmed_pmids=pubmed_pmids_list,
        europepmc_pmids=europepmc_pmids,
        combined_pmids=combined_pmids,
        synonyms_used=synonyms_list,
    )


def _save_pmids(pmids: List[str], output_file: Path) -> None:
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as handle:
        for pmid in pmids:
            handle.write(f"{pmid}\n")
    logger.info("Saved %d PMIDs to %s", len(pmids), output_file)
