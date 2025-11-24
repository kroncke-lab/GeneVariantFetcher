"""
DOI Resolver Module

Handles DOI resolution and routing to appropriate supplement scrapers
based on the publisher domain.
"""

import csv
from typing import List, Dict
from urllib.parse import urlparse
from pathlib import Path
import requests


class DOIResolver:
    """Resolves DOIs and routes to domain-specific supplement scrapers."""

    def __init__(
        self,
        session: requests.Session,
        paywalled_log: Path
    ):
        """
        Initialize DOI resolver.

        Args:
            session: Requests session with configured headers
            paywalled_log: Path to log file for failed DOI resolutions
        """
        self.session = session
        self.paywalled_log = paywalled_log

    def resolve_and_scrape_supplements(
        self,
        doi: str,
        pmid: str,
        scraper
    ) -> List[Dict]:
        """
        Resolves a DOI to its final URL and routes to the appropriate scraper.

        Args:
            doi: Digital Object Identifier
            pmid: PubMed ID (for logging)
            scraper: SupplementScraper instance to route to

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys
        """
        try:
            # Use the session with browser-like headers to resolve the DOI
            # allow_redirects=True follows the redirect chain to the final publisher page
            print(f"  Resolving DOI: https://doi.org/{doi}")
            response = self.session.get(f"https://doi.org/{doi}", allow_redirects=True, timeout=30)
            response.raise_for_status()

            final_url = response.url
            domain = urlparse(final_url).netloc
            print(f"  ✓ DOI resolved to: {final_url}")

        except requests.exceptions.RequestException as e:
            print(f"  ❌ DOI resolution failed for {doi}: {e}")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, f'DOI resolution failed: {doi}', f"https://doi.org/{doi}"])
            return []

        # Route to the specific scraper based on the resolved domain
        if "nature.com" in domain:
            return scraper.scrape_nature_supplements(response.text, final_url)
        elif "gimjournal.org" in domain or "sciencedirect.com" in domain:
            return scraper.scrape_elsevier_supplements(response.text, final_url)
        else:
            print(f"  - No specific scraper for domain: {domain}. Using generic scraper.")
            return scraper.scrape_generic_supplements(response.text, final_url)
