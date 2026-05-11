"""Tier 3.5 — Browser-based HTML harvester for free-after-embargo papers.

Sits between the requests-based publisher API tier and the manual-login
browser_fetch tier. Drives a real Playwright browser to fetch publisher HTML
that requires JavaScript or stricter fingerprinting than `requests` provides,
then delegates parsing back to the existing `SupplementScraper`.

Usage::

    from harvesting.browser_html import BrowserHTMLFetcher

    fetcher = BrowserHTMLFetcher(scraper, converter, session, settings, output_dir)
    result = fetcher.fetch(pmid="12345", doi="10.1161/...", pub_date=None)
    if result and result.main_markdown:
        ...

The tier can reuse an institutional Chrome profile when
``BROWSER_HTML_USE_PROFILE=true`` and ``BROWSER_HTML_PROFILE_PATH`` points at a
dedicated logged-in profile. Strategies are auto-discovered from the
``strategies/`` package — adding a new publisher is just dropping in a file.
"""

from .base import FetchContext, FetchResult, PublisherStrategy
from .embargo import EmbargoChecker, get_pub_date_from_pmid
from .fetcher import BrowserHTMLFetcher

__all__ = [
    "BrowserHTMLFetcher",
    "PublisherStrategy",
    "FetchResult",
    "FetchContext",
    "EmbargoChecker",
    "get_pub_date_from_pmid",
]
