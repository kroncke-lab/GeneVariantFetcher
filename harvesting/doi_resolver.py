"""
DOI Resolver Module

Handles DOI resolution and routing to appropriate supplement scrapers
based on the publisher domain.
Also handles full-text retrieval from free articles without PMCIDs.
"""

import csv
import re
from typing import List, Dict, Optional, Tuple
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

    def resolve_and_fetch_fulltext(
        self,
        doi: str,
        pmid: str,
        scraper
    ) -> Tuple[Optional[str], Optional[str], List[Dict]]:
        """
        Resolves a DOI to its final URL and fetches full text + supplements.

        This is used for free articles without PMCIDs. It:
        1. Resolves the DOI to the publisher's article page
        2. Extracts the full text content from the page
        3. Finds any supplemental files on the page

        Args:
            doi: Digital Object Identifier
            pmid: PubMed ID (for logging)
            scraper: SupplementScraper instance with fulltext extraction methods

        Returns:
            Tuple of (markdown_content, final_url, supplements_list)
            - markdown_content: Full text as markdown, or None if extraction failed
            - final_url: The resolved URL of the article
            - supplements_list: List of supplement file dictionaries
        """
        try:
            # Use the session with browser-like headers to resolve the DOI
            print(f"  Resolving DOI for full text: https://doi.org/{doi}")
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
                writer.writerow([pmid, f'DOI resolution failed (free text): {doi}', f"https://doi.org/{doi}"])
            return None, None, []
        except requests.exceptions.RequestException as e:
            print(f"  ❌ DOI resolution failed for {doi}: {e}")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, f'DOI resolution failed (free text): {doi}', f"https://doi.org/{doi}"])
            return None, None, []

        # Handle Elsevier linkinghub redirects
        if "linkinghub.elsevier.com" in domain:
            try:
                pii_match = re.search(r'/pii/([^/?]+)', final_url)
                if pii_match:
                    pii = pii_match.group(1)
                    sciencedirect_url = f"https://www.sciencedirect.com/science/article/pii/{pii}"
                    print(f"  → Attempting to access ScienceDirect page: {sciencedirect_url}")
                    redirect_response = self.session.get(sciencedirect_url, allow_redirects=True, timeout=30)
                    redirect_response.raise_for_status()
                    final_url = redirect_response.url
                    response = redirect_response
                    domain = urlparse(final_url).netloc
            except Exception as e:
                print(f"  - Could not follow redirect from linkinghub: {e}")

        # Extract full text using the scraper
        html_content = response.text
        markdown_content, title = scraper.extract_fulltext(html_content, final_url)

        if markdown_content:
            print(f"  ✓ Extracted full text from publisher page ({len(markdown_content)} characters)")
        else:
            print(f"  ❌ Could not extract full text from publisher page")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Full text extraction failed from publisher', final_url])

        # Also get supplements
        if "nature.com" in domain:
            supplements = scraper.scrape_nature_supplements(html_content, final_url)
        elif any(d in domain for d in ["gimjournal.org", "sciencedirect.com", "elsevier.com"]):
            supplements = scraper.scrape_elsevier_supplements(html_content, final_url)
        else:
            supplements = scraper.scrape_generic_supplements(html_content, final_url)

        return markdown_content, final_url, supplements

    def resolve_doi_url(self, doi: str, pmid: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Resolve a DOI to its final URL without fetching full content.

        Args:
            doi: Digital Object Identifier
            pmid: PubMed ID (for logging)

        Returns:
            Tuple of (final_url, html_content) or (None, None) if resolution failed
        """
        try:
            print(f"  Resolving DOI: https://doi.org/{doi}")
            response = self.session.get(f"https://doi.org/{doi}", allow_redirects=True, timeout=30)
            response.raise_for_status()

            final_url = response.url
            print(f"  ✓ DOI resolved to: {final_url}")

            return final_url, response.text

        except requests.exceptions.HTTPError as e:
            print(f"  ❌ DOI resolution failed for {doi}: HTTP {e.response.status_code}")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, f'DOI resolution failed: {doi}', f"https://doi.org/{doi}"])
            return None, None
        except requests.exceptions.RequestException as e:
            print(f"  ❌ DOI resolution failed for {doi}: {e}")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, f'DOI resolution failed: {doi}', f"https://doi.org/{doi}"])
            return None, None
