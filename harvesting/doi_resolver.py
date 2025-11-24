"""
DOI Resolver Module

Handles DOI resolution and routing to appropriate supplement scrapers
based on the publisher domain.
"""

import csv
import re
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

        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 403:
                print(f"  ❌ DOI resolution blocked (403 Forbidden) for {doi} - publisher may require authentication")
            else:
                print(f"  ❌ DOI resolution failed for {doi}: HTTP {e.response.status_code} - {e}")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, f'DOI resolution failed: {doi}', f"https://doi.org/{doi}"])
            return []
        except requests.exceptions.RequestException as e:
            print(f"  ❌ DOI resolution failed for {doi}: {e}")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, f'DOI resolution failed: {doi}', f"https://doi.org/{doi}"])
            return []

        # Route to the specific scraper based on the resolved domain
        if "nature.com" in domain:
            return scraper.scrape_nature_supplements(response.text, final_url)
        elif "gimjournal.org" in domain or "sciencedirect.com" in domain or "linkinghub.elsevier.com" in domain or "elsevier.com" in domain:
            # For linkinghub URLs, try to follow one more redirect to get to the actual article page
            if "linkinghub.elsevier.com" in domain:
                try:
                    # Try to extract PII from URL and construct ScienceDirect URL
                    # linkinghub URLs often have format: /retrieve/pii/S1547527109005682
                    pii_match = re.search(r'/pii/([^/?]+)', final_url)
                    if pii_match:
                        pii = pii_match.group(1)
                        sciencedirect_url = f"https://www.sciencedirect.com/science/article/pii/{pii}"
                        print(f"  → Attempting to access ScienceDirect page: {sciencedirect_url}")
                        redirect_response = self.session.get(sciencedirect_url, allow_redirects=True, timeout=30)
                        redirect_response.raise_for_status()
                        if redirect_response.url != final_url:
                            print(f"  → Following redirect to: {redirect_response.url}")
                            final_url = redirect_response.url
                            response = redirect_response
                    else:
                        # Fallback: try following redirect from linkinghub
                        redirect_response = self.session.get(final_url, allow_redirects=True, timeout=30)
                        redirect_response.raise_for_status()
                        if redirect_response.url != final_url:
                            print(f"  → Following redirect to: {redirect_response.url}")
                            final_url = redirect_response.url
                            response = redirect_response
                except Exception as e:
                    print(f"  - Could not follow redirect from linkinghub: {e}")
                    # Continue with linkinghub page anyway
            return scraper.scrape_elsevier_supplements(response.text, final_url)
        else:
            print(f"  - No specific scraper for domain: {domain}. Using generic scraper.")
            return scraper.scrape_generic_supplements(response.text, final_url)
