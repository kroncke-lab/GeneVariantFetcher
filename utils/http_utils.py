"""
HTTP utilities for consistent browser-like requests across the project.

This module provides standardized HTTP session configuration and headers
to avoid duplication across harvesting and fetching modules.
"""

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
