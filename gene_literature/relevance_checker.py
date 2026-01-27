"""LLM-based relevance checker for filtering irrelevant papers.

DEPRECATED: This module is maintained for backwards compatibility.
New code should import from gene_literature.llm_relevance instead.
"""

from gene_literature.llm_relevance import RelevanceChecker, RelevanceScore

__all__ = ["RelevanceChecker", "RelevanceScore"]
