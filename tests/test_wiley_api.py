#!/usr/bin/env python3
"""
Test script for the Wiley TDM API.

Usage:
    python tests/test_wiley_api.py                    # Run all tests
    python tests/test_wiley_api.py --doi "10.1002/xxx"  # Test specific DOI
    python tests/test_wiley_api.py --verbose          # Show full response
    python tests/test_wiley_api.py --web-only         # Test web scraping only

This tests the Wiley TDM API configuration and functionality.
"""

import argparse
import os
import sys
import logging
import time
from typing import Optional, Tuple
from urllib.parse import quote

import requests

# Add parent directory to path to import the actual module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Try to import the actual WileyAPIClient, fall back to inline version if needed
try:
    from harvesting.wiley_api import WileyAPIClient, WILEY_DOI_PREFIXES
    USING_REAL_MODULE = True
    logger.info("Using actual WileyAPIClient from harvesting.wiley_api")
except ImportError as e:
    logger.warning(f"Could not import WileyAPIClient: {e}")
    logger.warning("Using inline test version (web scraping not available)")
    USING_REAL_MODULE = False

    # Inline fallback for when module import fails
    WILEY_DOI_PREFIXES = (
        "10.1002/",   # Wiley main prefix
        "10.1111/",   # Wiley-Blackwell
        "10.1113/",   # The Physiological Society (Wiley)
        "10.1096/",   # FASEB Journal (Wiley)
        "10.1634/",   # Stem Cells (Wiley)
        "10.1111/j.", # Wiley journal articles
    )

    class WileyAPIClient:
        """Fallback client for testing when module import fails."""

        BASE_URL = "https://api.wiley.com/onlinelibrary/tdm/v1/articles"

        def __init__(self, api_key: Optional[str] = None, session: Optional[requests.Session] = None):
            self.api_key = api_key
            self.session = session or requests.Session()
            self._last_request_time = 0
            self._min_request_interval = 0.5

        @property
        def is_available(self) -> bool:
            return bool(self.api_key and self.api_key.strip())

        @staticmethod
        def is_wiley_doi(doi: str) -> bool:
            if not doi:
                return False
            return doi.lower().startswith(WILEY_DOI_PREFIXES)

        def _rate_limit(self):
            elapsed = time.time() - self._last_request_time
            if elapsed < self._min_request_interval:
                time.sleep(self._min_request_interval - elapsed)
            self._last_request_time = time.time()

        def get_fulltext_by_doi(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
            if not self.is_available:
                return None, "Wiley API key not configured"

            self._rate_limit()

            encoded_doi = quote(doi, safe='/:')
            url = f"{self.BASE_URL}/{encoded_doi}"
            headers = {
                "Wiley-TDM-Client-Token": self.api_key,
                "Accept": "application/xml, text/xml, application/xhtml+xml, text/html",
            }

            try:
                logger.info(f"Fetching full text from Wiley API for DOI: {doi}")
                logger.info(f"Request URL: {url}")
                response = self.session.get(url, headers=headers, timeout=30)
                logger.info(f"Response status: {response.status_code}")

                if response.status_code == 200:
                    return response.text, None
                elif response.status_code == 401:
                    return None, "Invalid or unauthorized API key"
                elif response.status_code == 403:
                    return None, "Access forbidden - API key may lack permissions or article not available"
                elif response.status_code == 404:
                    return None, "Article not found via Wiley API"
                elif response.status_code == 429:
                    return None, "Rate limit exceeded"
                else:
                    return None, f"HTTP {response.status_code}: {response.reason}"

            except requests.exceptions.Timeout:
                return None, "Request timed out"
            except requests.exceptions.RequestException as e:
                return None, f"Request failed: {str(e)}"

        def content_to_markdown(self, content: str) -> Optional[str]:
            if not content:
                return None
            return f"Content received: {len(content)} chars"

        def fetch_fulltext(self, doi: str = None, url: str = None, try_web_scraping: bool = True) -> Tuple[Optional[str], Optional[str]]:
            return None, "Web scraping not available in fallback mode"

        def scrape_fulltext_from_web(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
            return None, "Web scraping not available in fallback mode"

# Known Wiley test DOIs (mix of newer and older articles)
TEST_DOIS = [
    # Recent Human Mutation articles (should work)
    ("10.1002/humu.24231", "Recent Human Mutation article"),
    ("10.1002/humu.24114", "Recent Human Mutation article"),

    # Recent JAHA (Wiley) articles
    ("10.1161/JAHA.121.022754", "Recent JAHA article"),

    # Older SICI-style DOI (may not work - like the user's failing case)
    ("10.1002/(SICI)1098-1004(200005)15:5<483::AID-HUMU18>3.0.CO;2-T", "Old SICI-style DOI"),

    # Another format test
    ("10.1111/cge.13351", "Wiley-Blackwell journal"),
]


def test_api_configuration():
    """Test that the API client is properly configured."""
    print("=" * 60)
    print("Testing Wiley API Configuration")
    print("=" * 60)

    # Check environment variable
    api_key = os.environ.get("WILEY_API_KEY")
    if api_key:
        print(f"✓ WILEY_API_KEY is set (length: {len(api_key)} chars)")
        print(f"  Key prefix: {api_key[:8]}...")
    else:
        print("✗ WILEY_API_KEY is NOT set in environment")

        # Try loading from settings
        try:
            from config.settings import Settings
            settings = Settings()
            if settings.wiley_api_key:
                api_key = settings.wiley_api_key
                print(f"✓ Found key in settings (length: {len(api_key)} chars)")
            else:
                print("✗ No key found in settings either")
        except Exception as e:
            print(f"  Could not load settings: {e}")

    # Initialize client
    client = WileyAPIClient(api_key=api_key)
    print(f"\nClient is_available: {client.is_available}")

    return client


def test_doi(client: WileyAPIClient, doi: str, description: str = "", verbose: bool = False):
    """Test fetching a single DOI via TDM API only."""
    print(f"\n{'─' * 50}")
    print(f"Testing DOI: {doi}")
    if description:
        print(f"Description: {description}")
    print(f"Is Wiley DOI: {client.is_wiley_doi(doi)}")

    if not client.is_available:
        print("⚠ Cannot test - API key not configured")
        return False

    # Test raw API call first
    print("Making API request...")
    content, error = client.get_fulltext_by_doi(doi)

    if error:
        print(f"✗ Error: {error}")
        return False

    if content:
        print(f"✓ Received response ({len(content)} chars)")

        if verbose:
            print(f"\nResponse preview (first 500 chars):")
            print("-" * 40)
            print(content[:500])
            print("-" * 40)

        # Try converting to markdown
        markdown = client.content_to_markdown(content)
        if markdown:
            print(f"✓ Converted to markdown ({len(markdown)} chars)")
            if verbose:
                print(f"\nMarkdown preview (first 500 chars):")
                print("-" * 40)
                print(markdown[:500])
                print("-" * 40)
        else:
            print("✗ Failed to convert to markdown")

        return True

    print("✗ No content received (but no error either)")
    return False


def test_web_scraping(client: WileyAPIClient, doi: str, description: str = "", verbose: bool = False):
    """Test fetching a DOI via web scraping (fallback method)."""
    print(f"\n{'─' * 50}")
    print(f"Testing Web Scraping for DOI: {doi}")
    if description:
        print(f"Description: {description}")

    if not USING_REAL_MODULE:
        print("⚠ Cannot test - real module not available (web scraping requires actual module)")
        return False

    if not hasattr(client, 'scrape_fulltext_from_web'):
        print("⚠ Cannot test - scrape_fulltext_from_web method not available")
        return False

    print("Making web scraping request...")
    markdown, error = client.scrape_fulltext_from_web(doi)

    if error:
        print(f"✗ Error: {error}")
        return False

    if markdown:
        print(f"✓ Successfully scraped content ({len(markdown)} chars)")

        if verbose:
            print(f"\nMarkdown preview (first 1000 chars):")
            print("-" * 40)
            print(markdown[:1000])
            print("-" * 40)

        # Check for key article sections
        has_abstract = '### Abstract' in markdown or '### abstract' in markdown.lower()
        has_body = any(section in markdown.lower() for section in
                      ['### introduction', '### methods', '### results', '### discussion'])

        print(f"  Has Abstract: {'✓' if has_abstract else '✗'}")
        print(f"  Has Body Sections: {'✓' if has_body else '✗'}")

        return True

    print("✗ No content received")
    return False


def test_fetch_fulltext(client: WileyAPIClient, doi: str, description: str = "", verbose: bool = False):
    """Test the full fetch_fulltext method (API + web scraping fallback)."""
    print(f"\n{'─' * 50}")
    print(f"Testing fetch_fulltext (API + fallback) for DOI: {doi}")
    if description:
        print(f"Description: {description}")

    if not USING_REAL_MODULE:
        print("⚠ Cannot test - real module not available")
        return False

    if not hasattr(client, 'fetch_fulltext'):
        print("⚠ Cannot test - fetch_fulltext method not available")
        return False

    print("Making fetch_fulltext request (will try API then web scraping)...")
    markdown, error = client.fetch_fulltext(doi=doi, try_web_scraping=True)

    if error:
        print(f"✗ Error: {error}")
        return False

    if markdown:
        print(f"✓ Successfully fetched content ({len(markdown)} chars)")

        if verbose:
            print(f"\nMarkdown preview (first 1000 chars):")
            print("-" * 40)
            print(markdown[:1000])
            print("-" * 40)

        return True

    print("✗ No content received")
    return False


def main():
    parser = argparse.ArgumentParser(description="Test Wiley TDM API")
    parser.add_argument("--doi", help="Test a specific DOI")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show response content")
    parser.add_argument("--all", "-a", action="store_true", help="Test all sample DOIs")
    parser.add_argument("--key", "-k", help="Wiley API key (or set WILEY_API_KEY env var)")
    parser.add_argument("--web-only", "-w", action="store_true", help="Test web scraping only (no API)")
    parser.add_argument("--full", "-f", action="store_true", help="Test full fetch_fulltext (API + web fallback)")
    args = parser.parse_args()

    # If key provided via arg, set it in env for test_api_configuration
    if args.key:
        os.environ["WILEY_API_KEY"] = args.key

    client = test_api_configuration()

    if args.doi:
        if args.web_only:
            # Test web scraping only
            success = test_web_scraping(client, args.doi, "User-provided DOI", args.verbose)
        elif args.full:
            # Test full fetch with fallback
            success = test_fetch_fulltext(client, args.doi, "User-provided DOI", args.verbose)
        else:
            # Test API only
            success = test_doi(client, args.doi, "User-provided DOI", args.verbose)
        sys.exit(0 if success else 1)

    if args.all or not args.doi:
        # Test sample DOIs
        print("\n" + "=" * 60)
        if args.web_only:
            print("Testing Web Scraping for Sample DOIs")
        elif args.full:
            print("Testing Full Fetch (API + Web Fallback) for Sample DOIs")
        else:
            print("Testing TDM API for Sample DOIs")
        print("=" * 60)

        results = []
        for doi, description in TEST_DOIS:
            if args.web_only:
                success = test_web_scraping(client, doi, description, args.verbose)
            elif args.full:
                success = test_fetch_fulltext(client, doi, description, args.verbose)
            else:
                success = test_doi(client, doi, description, args.verbose)
            results.append((doi, success))

        # Summary
        print("\n" + "=" * 60)
        print("Summary")
        print("=" * 60)

        passed = sum(1 for _, success in results if success)
        failed = len(results) - passed

        for doi, success in results:
            status = "✓" if success else "✗"
            print(f"  {status} {doi[:50]}...")

        print(f"\nPassed: {passed}/{len(results)}")

        if failed > 0:
            if args.web_only:
                print("\nNote: Some DOIs may fail because:")
                print("  - The article is behind a paywall")
                print("  - The article page structure changed")
                print("  - Network issues")
            else:
                print("\nNote: Some DOIs may fail because:")
                print("  - The article is not available via TDM API")
                print("  - It's an older article not indexed in TDM")
                print("  - The DOI format is not supported (SICI-style)")
                print("  - API rate limiting")
                if args.full:
                    print("  - Web scraping fallback also failed (paywall/access issue)")


if __name__ == "__main__":
    main()
