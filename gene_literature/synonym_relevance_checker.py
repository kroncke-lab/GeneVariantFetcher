"""LLM-based relevance checker for filtering gene synonyms.

DEPRECATED: This module is maintained for backwards compatibility.
New code should import from gene_literature.llm_relevance instead.
"""

from gene_literature.llm_relevance import SynonymRelevanceChecker, SynonymRelevance

__all__ = ["SynonymRelevanceChecker", "SynonymRelevance"]
