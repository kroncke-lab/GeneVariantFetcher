"""High-level orchestration for gene-focused literature collection."""

from __future__ import annotations

import logging
from typing import Optional, Sequence

logger = logging.getLogger(__name__)


def build_gene_query(gene: str, synonyms: Optional[Sequence[str]] = None) -> str:
    """Build a simple PubMed query using the provided gene and synonyms."""

    terms = [gene, *(synonyms or [])]
    sanitized = [term.strip() for term in terms if term and term.strip()]
    quoted = [f'"{term}"[Title/Abstract]' for term in sanitized]
    if not quoted:
        raise ValueError("At least one search term must be provided")
    query = " OR ".join(quoted)
    logger.debug("Constructed PubMed query: %s", query)
    return query
