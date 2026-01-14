"""
HTTP utilities for consistent browser-like requests across the project.

This module provides standardized HTTP session configuration and headers
to avoid duplication across harvesting and fetching modules.
"""

import requests
from typing import Dict

# Browser-like headers to avoid blocks from websites
BROWSER_HEADERS: Dict[str, str] = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.5",
    "Accept-Encoding": "gzip, deflate, br",
    "Connection": "keep-alive",
    "Upgrade-Insecure-Requests": "1",
}


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
        headers["Accept"] = "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8"

    return headers
