"""
HTTP utilities for consistent browser-like requests across the project.

This module provides standardized HTTP session configuration and headers
to avoid duplication across harvesting and fetching modules.
"""

from typing import Dict

import requests

from config.constants import BROWSER_HEADERS


def get_browser_session() -> requests.Session:
    """
    Create a requests session with browser-like headers.

    Returns:
        A configured requests.Session with headers that mimic a real browser,
        reducing the likelihood of being blocked by websites.

    Example:
        >>> session = get_browser_session()
        >>> response = session.get("https://example.com")
    """
    session = requests.Session()
    session.headers.update(BROWSER_HEADERS)
    return session


def get_headers_for_domain(domain: str) -> Dict[str, str]:
    """
    Get appropriate headers for a specific domain.

    Some domains may require specific headers. This function returns
    the standard browser headers with any domain-specific modifications.

    Args:
        domain: The domain being accessed (e.g., "pubmed.ncbi.nlm.nih.gov")

    Returns:
        Dictionary of headers appropriate for the domain
    """
    headers = BROWSER_HEADERS.copy()

    # Add domain-specific headers if needed
    if "ncbi.nlm.nih.gov" in domain:
        headers["Accept"] = (
            "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8"
        )

    return headers
