#!/usr/bin/env python3
"""
Test downloading non-PMC papers using available APIs and methods.

This tests the full download pipeline for papers that aren't in PMC,
using Elsevier API, Wiley API, DOI resolution, and web scraping.
"""

import json
import logging
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

# Add parent to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load environment
from dotenv import load_dotenv

load_dotenv(Path(__file__).parent.parent / ".env")

from config.settings import get_settings
from harvesting.doi_resolver import DOIResolver
from harvesting.elsevier_api import ElsevierAPIClient
from harvesting.supplement_scraper import SupplementScraper
from harvesting.wiley_api import WileyAPIClient

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


@dataclass
class DownloadResult:
    """Result of a download attempt."""

    pmid: str
    doi: Optional[str]
    publisher: str
    success: bool
    method: str = "none"
    content_length: int = 0
    has_fulltext: bool = False
    has_supplements: bool = False
    error: Optional[str] = None
    supplements_found: List[str] = field(default_factory=list)


def test_elsevier_download(
    client: ElsevierAPIClient, pmid: str, doi: str
) -> DownloadResult:
    """Test Elsevier API download."""
    result = DownloadResult(pmid=pmid, doi=doi, publisher="Elsevier", success=False)

    if not client.is_available:
        result.error = "Elsevier API key not configured"
        return result

    try:
        # Try fetching full text
        markdown, error = client.get_fulltext_by_doi(doi)

        if markdown:
            result.success = True
            result.method = "elsevier_api"
            result.content_length = len(markdown)
            result.has_fulltext = len(markdown) > 5000  # Reasonable threshold

            # Check for supplements - Elsevier API doesn't provide them directly
            # but we can check if the content mentions supplementary materials
            if "supplement" in markdown.lower() or "appendix" in markdown.lower():
                result.has_supplements = True
        else:
            result.error = error

    except Exception as e:
        result.error = str(e)

    return result


def test_wiley_download(client: WileyAPIClient, pmid: str, doi: str) -> DownloadResult:
    """Test Wiley API download."""
    result = DownloadResult(pmid=pmid, doi=doi, publisher="Wiley", success=False)

    if not client.is_available:
        result.error = "Wiley API key not configured"
        return result

    try:
        # Try API first - get_fulltext_by_doi already returns markdown
        # (it converts PDFs internally)
        content, error = client.get_fulltext_by_doi(doi)

        if content:
            # Content is already markdown if PDF was converted, or XML/HTML
            # Check if it looks like content (not an error message)
            if len(content) > 1000:
                result.success = True
                result.method = "wiley_api"
                result.content_length = len(content)
                result.has_fulltext = len(content) > 5000
                return result
            else:
                # Very short content - might be XML/HTML, try converting
                markdown = client.content_to_markdown(content)
                if markdown and len(markdown) > 1000:
                    result.success = True
                    result.method = "wiley_api"
                    result.content_length = len(markdown)
                    result.has_fulltext = len(markdown) > 5000
                    return result

        # Try web scraping fallback if available
        if hasattr(client, "scrape_fulltext_from_web"):
            markdown, scrape_error = client.scrape_fulltext_from_web(doi)
            if markdown:
                result.success = True
                result.method = "wiley_web_scrape"
                result.content_length = len(markdown)
                result.has_fulltext = len(markdown) > 5000
                return result

        result.error = error

    except Exception as e:
        result.error = str(e)

    return result


def test_doi_resolver(
    resolver: DOIResolver,
    scraper: SupplementScraper,
    pmid: str,
    doi: str,
    publisher: str,
) -> DownloadResult:
    """Test DOI resolution and free text fetch."""
    result = DownloadResult(pmid=pmid, doi=doi, publisher=publisher, success=False)

    try:
        # Try to resolve DOI and fetch content
        # resolve_and_fetch_fulltext returns (markdown, final_url, supplements_list)
        content, final_url, supplements = resolver.resolve_and_fetch_fulltext(
            doi, pmid, scraper
        )

        if content:
            result.success = True
            result.method = "doi_resolver"
            result.content_length = len(content)
            result.has_fulltext = len(content) > 5000
            if supplements:
                result.has_supplements = True
                result.supplements_found = [
                    s.get("name", s.get("url", "unknown")) for s in supplements
                ]
            return result

        result.error = f"No content from DOI resolver (URL: {final_url})"

    except Exception as e:
        result.error = str(e)

    return result


def main():
    settings = get_settings()

    # Initialize clients
    elsevier_client = ElsevierAPIClient(api_key=settings.elsevier_api_key)
    wiley_client = WileyAPIClient(api_key=settings.wiley_api_key)

    # DOI Resolver needs a session and log path
    import requests

    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"
        }
    )
    paywalled_log = (
        Path(__file__).parent.parent / "tests/fixtures/pmids/test_paywalled.csv"
    )
    doi_resolver = DOIResolver(session=session, paywalled_log=paywalled_log)

    # SupplementScraper for full-text extraction
    scraper = SupplementScraper()

    print("=" * 70)
    print("NON-PMC PAPER DOWNLOAD TEST")
    print("=" * 70)
    print(
        f"Elsevier API: {'✓ Configured' if elsevier_client.is_available else '✗ Not configured'}"
    )
    print(
        f"Wiley API: {'✓ Configured' if wiley_client.is_available else '✗ Not configured'}"
    )
    print()

    # Load test data
    csv_file = Path(__file__).parent.parent / "tests/fixtures/pmids/non_pmc_papers.csv"
    if not csv_file.exists():
        print(f"Error: {csv_file} not found. Run check_publishers.py first.")
        sys.exit(1)

    import csv

    papers = []
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        papers = list(reader)

    print(f"Testing {len(papers)} non-PMC papers...\n")

    results = []

    for i, paper in enumerate(papers, 1):
        pmid = paper["PMID"]
        doi = paper["DOI"]
        publisher = paper["Publisher"]

        print(f"[{i}/{len(papers)}] PMID {pmid} ({publisher})")
        print(f"   DOI: {doi}")

        result = None

        # Route to appropriate API
        if "Elsevier" in publisher:
            result = test_elsevier_download(elsevier_client, pmid, doi)
        elif "Wiley" in publisher:
            result = test_wiley_download(wiley_client, pmid, doi)
        else:
            # Try DOI resolver for other publishers (AHA, Springer, etc.)
            result = test_doi_resolver(doi_resolver, scraper, pmid, doi, publisher)

        results.append(result)

        if result.success:
            print(f"   ✓ SUCCESS via {result.method}")
            print(
                f"     Content: {result.content_length:,} chars, Fulltext: {result.has_fulltext}"
            )
        else:
            print(f"   ✗ FAILED: {result.error}")

        # Rate limiting
        time.sleep(0.5)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    successes = [r for r in results if r.success]
    failures = [r for r in results if not r.success]

    print(
        f"\nSuccessful: {len(successes)}/{len(results)} ({100 * len(successes) / len(results):.0f}%)"
    )

    if successes:
        print("\n✓ Successful downloads:")
        by_method = {}
        for r in successes:
            by_method.setdefault(r.method, []).append(r)
        for method, items in by_method.items():
            print(f"   {method}: {len(items)} papers")
            for r in items:
                print(f"      {r.pmid} ({r.publisher}) - {r.content_length:,} chars")

    if failures:
        print(f"\n✗ Failed downloads ({len(failures)}):")
        for r in failures:
            print(f"   {r.pmid} ({r.publisher}): {r.error}")

    # Write results JSON
    output_file = (
        Path(__file__).parent.parent / "tests/fixtures/pmids/download_test_results.json"
    )
    with open(output_file, "w") as f:
        json.dump(
            [
                {
                    "pmid": r.pmid,
                    "doi": r.doi,
                    "publisher": r.publisher,
                    "success": r.success,
                    "method": r.method,
                    "content_length": r.content_length,
                    "has_fulltext": r.has_fulltext,
                    "error": r.error,
                }
                for r in results
            ],
            f,
            indent=2,
        )

    print(f"\nResults written to: {output_file}")


if __name__ == "__main__":
    main()
