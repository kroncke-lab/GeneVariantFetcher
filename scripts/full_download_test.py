#!/usr/bin/env python3
"""
Full download test for all test PMIDs.
Tries: Elsevier API → Wiley API → DOI Resolver → Browser Fetch (Playwright)
"""

import json
import os
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dotenv import load_dotenv

load_dotenv(Path(__file__).parent.parent / ".env")

from config.settings import get_settings
from harvesting.elsevier_api import ElsevierAPIClient
from harvesting.pmc_api import PMCAPIClient
from harvesting.wiley_api import WileyAPIClient

# Playwright imports
try:
    from playwright.sync_api import sync_playwright

    PLAYWRIGHT_AVAILABLE = True
except ImportError:
    PLAYWRIGHT_AVAILABLE = False

import logging

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


@dataclass
class DownloadResult:
    pmid: str
    doi: Optional[str] = None
    pmcid: Optional[str] = None
    publisher: str = "unknown"
    success: bool = False
    method: str = "none"
    content_length: int = 0
    error: Optional[str] = None


def get_doi_for_pmid(pmc_client: PMCAPIClient, pmid: str) -> Optional[str]:
    """Get DOI for a PMID."""
    try:
        return pmc_client.get_doi_from_pmid(pmid)
    except:
        return None


def identify_publisher(doi: str) -> str:
    """Identify publisher from DOI prefix."""
    if not doi:
        return "unknown"
    prefixes = {
        "10.1016/": "Elsevier",
        "10.1053/": "Elsevier",
        "10.1002/": "Wiley",
        "10.1111/": "Wiley",
        "10.1161/": "AHA",
        "10.1093/": "Oxford",
        "10.1038/": "Nature",
        "10.1007/": "Springer",
    }
    for prefix, publisher in prefixes.items():
        if doi.startswith(prefix):
            return publisher
    return "Other"


def try_elsevier(client: ElsevierAPIClient, doi: str) -> tuple[bool, str, int]:
    """Try Elsevier API."""
    if not client.is_available or not doi or not doi.startswith("10.1016"):
        return False, "Not Elsevier DOI", 0
    content, error = client.get_fulltext_by_doi(doi)
    if content and len(content) > 1000:
        return True, content, len(content)
    return False, error or "No content", 0


def try_wiley(client: WileyAPIClient, doi: str) -> tuple[bool, str, int]:
    """Try Wiley API."""
    if not client.is_available or not doi:
        return False, "Wiley API not available", 0
    if not (doi.startswith("10.1002") or doi.startswith("10.1111")):
        return False, "Not Wiley DOI", 0
    content, error = client.get_fulltext_by_doi(doi)
    if content and len(content) > 1000:
        return True, content, len(content)
    return False, error or "No content", 0


def try_browser_fetch(doi: str, pmid: str) -> tuple[bool, str, int]:
    """Try Playwright browser fetch."""
    if not PLAYWRIGHT_AVAILABLE:
        return False, "Playwright not installed", 0

    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(
                headless=True,
                args=[
                    "--disable-blink-features=AutomationControlled",
                    "--disable-dev-shm-usage",
                ],
            )
            context = browser.new_context(
                user_agent="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
                viewport={"width": 1920, "height": 1080},
            )
            page = context.new_page()
            page.add_init_script(
                "Object.defineProperty(navigator, 'webdriver', {get: () => undefined});"
            )

            url = (
                f"https://doi.org/{doi}"
                if doi
                else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            )
            page.goto(url, timeout=30000, wait_until="networkidle")

            # Wait for any challenge
            for _ in range(5):
                if "Just a moment" not in page.title():
                    break
                time.sleep(1)

            body = page.query_selector("body")
            content = body.inner_text() if body else ""
            browser.close()

            # Check if we got real content
            if len(content) > 2000 and (
                "abstract" in content.lower() or "method" in content.lower()
            ):
                return True, content, len(content)
            return False, "No usable content extracted", len(content)

    except Exception as e:
        return False, str(e), 0


def main():
    settings = get_settings()

    # Initialize clients
    elsevier = ElsevierAPIClient(api_key=settings.elsevier_api_key)
    wiley = WileyAPIClient(api_key=settings.wiley_api_key)
    pmc = PMCAPIClient()

    print("=" * 70)
    print("FULL DOWNLOAD TEST - ALL TEST PMIDs")
    print("=" * 70)
    print(f"Elsevier API: {'✓' if elsevier.is_available else '✗'}")
    print(f"Wiley API: {'✓' if wiley.is_available else '✗'}")
    print(f"Playwright: {'✓' if PLAYWRIGHT_AVAILABLE else '✗'}")
    print()

    # Load test PMIDs
    test_file = Path(__file__).parent.parent / "tests/fixtures/pmids/test_pmids.txt"
    pmids = []
    with open(test_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                pmids.append(line.split()[0])

    print(f"Testing {len(pmids)} PMIDs...\n")

    results = []

    for i, pmid in enumerate(pmids, 1):
        result = DownloadResult(pmid=pmid)

        # Get metadata
        result.doi = get_doi_for_pmid(pmc, pmid)
        result.pmcid = pmc.pmid_to_pmcid(pmid)
        result.publisher = identify_publisher(result.doi)

        print(f"[{i}/{len(pmids)}] PMID {pmid} ({result.publisher})")
        if result.doi:
            print(f"    DOI: {result.doi}")
        if result.pmcid:
            print(f"    PMCID: {result.pmcid} (in PMC - standard path works)")
            result.success = True
            result.method = "pmc"
            results.append(result)
            continue

        # Try download methods in order

        # 1. Elsevier API
        if result.publisher == "Elsevier":
            success, content, length = try_elsevier(elsevier, result.doi)
            if success:
                result.success = True
                result.method = "elsevier_api"
                result.content_length = length
                print(f"    ✓ Elsevier API: {length:,} chars")
                results.append(result)
                continue

        # 2. Wiley API
        if result.publisher == "Wiley":
            success, content, length = try_wiley(wiley, result.doi)
            if success:
                result.success = True
                result.method = "wiley_api"
                result.content_length = length
                print(f"    ✓ Wiley API: {length:,} chars")
                results.append(result)
                continue

        # 3. Browser fetch (for everything else, including AHA)
        print("    Trying browser fetch...")
        success, content, length = try_browser_fetch(result.doi, pmid)
        if success:
            result.success = True
            result.method = "browser_fetch"
            result.content_length = length
            print(f"    ✓ Browser fetch: {length:,} chars")
        else:
            result.error = content  # Error message stored in content
            print(f"    ✗ Failed: {content[:80]}")

        results.append(result)
        time.sleep(0.5)  # Rate limiting

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    success_count = sum(1 for r in results if r.success)
    print(
        f"\nTotal: {success_count}/{len(results)} ({100 * success_count / len(results):.0f}%)"
    )

    # By method
    by_method = {}
    for r in results:
        if r.success:
            by_method.setdefault(r.method, []).append(r)

    print("\nBy method:")
    for method, items in sorted(by_method.items(), key=lambda x: -len(x[1])):
        print(f"  {method}: {len(items)}")

    # Failed
    failed = [r for r in results if not r.success]
    if failed:
        print(f"\nFailed ({len(failed)}):")
        for r in failed:
            print(
                f"  {r.pmid} ({r.publisher}): {r.error[:60] if r.error else 'Unknown'}..."
            )

    # Save results
    output_file = (
        Path(__file__).parent.parent / "tests/fixtures/pmids/full_test_results.json"
    )
    with open(output_file, "w") as f:
        json.dump([asdict(r) for r in results], f, indent=2)
    print(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
