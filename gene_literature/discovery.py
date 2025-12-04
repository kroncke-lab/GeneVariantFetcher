"""Utilities for combining PubMind and PubMed PMID discovery."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence

from .collector import build_gene_query
from .pubmed_client import PubMedClient, PubMedError
from .pubmind_fetcher import PubMindFetcher

logger = logging.getLogger(__name__)


@dataclass
class PMIDDiscoveryResult:
    """Structured representation of PMID discovery results."""

    pubmind_pmids: List[str]
    pubmed_pmids: List[str]
    combined_pmids: List[str]


def build_gene_keyword_queries(gene_symbol: str) -> List[str]:
    """Return a set of PubMed queries focused on a gene and variant keywords."""

    base_query = build_gene_query(gene_symbol)
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

    logger.debug("Built %d PubMed keyword queries for %s", len(queries), gene_symbol)
    return queries


def discover_pmids_for_gene(
    gene_symbol: str,
    *,
    email: str,
    max_results: int,
    pubmind_output: Path | None = None,
    pubmed_output: Path | None = None,
    combined_output: Path | None = None,
    api_key: str | None = None,
) -> PMIDDiscoveryResult:
    """Fetch PMIDs from PubMind and PubMed, merge, and persist them."""

    # PubMind discovery
    pubmind_fetcher = PubMindFetcher(email=email)
    pubmind_pmids = pubmind_fetcher.fetch_pmids_for_gene(
        gene_symbol, max_results=max_results
    )
    if pubmind_output:
        pubmind_fetcher.save_pmids_to_file(pubmind_pmids, pubmind_output)

    # PubMed discovery using keyword-focused queries
    pubmed_client = PubMedClient(api_key=api_key, email=email)
    pubmed_pmids: set[str] = set()
    for query in build_gene_keyword_queries(gene_symbol):
        try:
            pubmed_pmids.update(
                pubmed_client.search(query, retmax=max_results)
            )
        except PubMedError as exc:
            logger.warning("PubMed query failed for '%s': %s", query, exc)

    pubmed_pmids_list = sorted(pubmed_pmids)[:max_results]
    if pubmed_output:
        _save_pmids(pubmed_pmids_list, pubmed_output)

    # Merge and de-duplicate
    combined_pmids = sorted(set(pubmind_pmids) | set(pubmed_pmids_list))
    if len(combined_pmids) > max_results:
        combined_pmids = combined_pmids[:max_results]

    if combined_output:
        _save_pmids(combined_pmids, combined_output)

    return PMIDDiscoveryResult(
        pubmind_pmids=pubmind_pmids,
        pubmed_pmids=pubmed_pmids_list,
        combined_pmids=combined_pmids,
    )


def _save_pmids(pmids: List[str], output_file: Path) -> None:
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as handle:
        for pmid in pmids:
            handle.write(f"{pmid}\n")
    logger.info("Saved %d PMIDs to %s", len(pmids), output_file)
