"""
DOI Resolver Module

Handles DOI resolution and routing to appropriate supplement scrapers
based on the publisher domain.
Also handles full-text retrieval from free articles without PMCIDs.
"""

import csv
import logging
import random
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote, urlparse

import requests

from config.constants import HTTP_TIMEOUT_DEFAULT

logger = logging.getLogger(__name__)


def encode_doi_for_url(doi: str) -> str:
    """
    Encode a DOI for use in URLs.

    DOIs can contain special characters like <, >, (, ), etc. that need to be
    percent-encoded when used in URLs. The safe characters are those that
    don't need encoding in the path component of a URL.

    Args:
        doi: The raw DOI string

    Returns:
        URL-encoded DOI string
    """
    # Encode the DOI, keeping / and : as safe since they're common in DOIs
    # and handled correctly by doi.org
    return quote(doi, safe="/:")  # Keep common DOI characters unencoded


# User-Agent strings to rotate through when encountering 403 errors
# These represent various browsers and research tools
USER_AGENTS = [
    # Chrome on macOS
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36",
    # Chrome on Windows
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36",
    # Firefox on Windows
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:133.0) Gecko/20100101 Firefox/133.0",
    # Safari on macOS
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.2 Safari/605.1.15",
    # Edge on Windows
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36 Edg/131.0.0.0",
    # Research-oriented bot (some publishers whitelist these)
    "Mozilla/5.0 (compatible; Googlebot/2.1; +http://www.google.com/bot.html)",
]


class DOIResolver:
    """Resolves DOIs and routes to domain-specific supplement scrapers."""

    # Number of retry attempts for 403 errors
    MAX_RETRIES = 3
    # Initial backoff delay in seconds
    INITIAL_BACKOFF = 2.0

    def __init__(self, session: requests.Session, paywalled_log: Path):
        """
        Initialize DOI resolver.

        Args:
            session: Requests session with configured headers
            paywalled_log: Path to log file for failed DOI resolutions
        """
        self.session = session
        self.paywalled_log = paywalled_log
        self._user_agent_index = 0

    def _log_paywalled(
        self, pmid: str, reason: str, url: str, classification: str = ""
    ) -> None:
        """
        Log a paper to the paywalled/missing CSV with placeholder columns.

        Args:
            pmid: PubMed ID
            reason: Why the paper couldn't be downloaded
            url: URL attempted
            classification: One of PAYWALLED, CAPTCHA_BLOCKED,
                INSTITUTIONAL_ACCESS, SUPPLEMENT_ONLY, API_LIMIT, or empty.
        """
        with open(self.paywalled_log, "a", newline="") as f:
            writer = csv.writer(f)
            # Write with empty placeholders for carrier columns (matching orchestrator format)
            writer.writerow(
                [
                    pmid,
                    reason,
                    url,
                    classification,
                    "",
                    "",
                    "",  # Abstract_Carriers, Affected_Count, Unaffected_Count
                    "",
                    "",  # Variants_Mentioned, Extraction_Confidence
                    "",
                    "",
                    "",  # More_In_Fulltext_Probability, Priority_Score, Notes
                ]
            )

    def _get_next_user_agent(self) -> str:
        """
        Get the next User-Agent in rotation.

        Returns:
            User-Agent string
        """
        ua = USER_AGENTS[self._user_agent_index % len(USER_AGENTS)]
        self._user_agent_index += 1
        return ua

    def _resolve_with_retry(
        self, url: str, description: str = "DOI"
    ) -> Optional[requests.Response]:
        """
        Resolve a URL with retry logic and User-Agent rotation for 403 errors.

        Args:
            url: URL to resolve
            description: Description for logging (e.g., "DOI", "ScienceDirect page")

        Returns:
            Response object if successful, None if all retries failed
        """
        last_exception = None
        last_status = None

        for attempt in range(self.MAX_RETRIES):
            try:
                # Rotate User-Agent for retries
                if attempt > 0:
                    new_ua = self._get_next_user_agent()
                    self.session.headers["User-Agent"] = new_ua
                    # Add a small random delay to avoid detection patterns
                    backoff = self.INITIAL_BACKOFF * (
                        2 ** (attempt - 1)
                    ) + random.uniform(0.5, 1.5)
                    logger.debug(
                        f"Retry {attempt + 1}/{self.MAX_RETRIES} for {description} with new User-Agent, waiting {backoff:.1f}s"
                    )
                    time.sleep(backoff)

                response = self.session.get(url, allow_redirects=True, timeout=HTTP_TIMEOUT_DEFAULT)
                response.raise_for_status()
                return response

            except requests.exceptions.HTTPError as e:
                last_exception = e
                last_status = e.response.status_code if e.response else None

                # Only retry on 403 errors (access denied)
                if last_status == 403:
                    if attempt < self.MAX_RETRIES - 1:
                        logger.debug(
                            f"Got 403 for {description}, will retry with different User-Agent"
                        )
                        continue
                    else:
                        logger.debug(
                            f"All retries exhausted for {description} (403 Forbidden)"
                        )
                        break
                else:
                    # Don't retry other HTTP errors
                    break

            except requests.exceptions.RequestException as e:
                last_exception = e
                # Network errors - retry with backoff
                if attempt < self.MAX_RETRIES - 1:
                    backoff = self.INITIAL_BACKOFF * (2**attempt)
                    logger.debug(
                        f"Network error for {description}, retrying in {backoff:.1f}s: {e}"
                    )
                    time.sleep(backoff)
                    continue
                break

        # All retries failed
        if last_status == 403:
            raise requests.exceptions.HTTPError(
                f"403 Forbidden after {self.MAX_RETRIES} attempts",
                response=last_exception.response
                if hasattr(last_exception, "response")
                else None,
            )
        elif last_exception:
            raise last_exception
        return None

    def resolve_and_scrape_supplements(
        self, doi: str, pmid: str, scraper
    ) -> List[Dict]:
        """
        Resolves a DOI to its final URL and routes to the appropriate scraper.

        Uses retry logic with User-Agent rotation to handle 403 errors from
        publishers with anti-bot measures.

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
            encoded_doi = encode_doi_for_url(doi)
            response = self._resolve_with_retry(
                f"https://doi.org/{encoded_doi}", f"DOI {doi}"
            )

            if not response:
                print(f"  ❌ DOI resolution failed for {doi}: No response received")
                self._log_paywalled(
                    pmid, f"DOI resolution failed: {doi}", f"https://doi.org/{doi}",
                    classification="PAYWALLED"
                )
                return []

            final_url = response.url
            domain = urlparse(final_url).netloc
            print(f"  ✓ DOI resolved to: {final_url}")

        except requests.exceptions.HTTPError as e:
            status_code = e.response.status_code if e.response else "unknown"
            if status_code == 403:
                print(
                    f"  ❌ DOI resolution blocked (403 Forbidden) for {doi} after {self.MAX_RETRIES} attempts - publisher may require authentication"
                )
            else:
                print(f"  ❌ DOI resolution failed for {doi}: HTTP {status_code}")
            self._log_paywalled(
                pmid,
                f"DOI resolution failed (HTTP {status_code}): {doi}",
                f"https://doi.org/{doi}",
                classification="CAPTCHA_BLOCKED" if status_code == 403 else "PAYWALLED"
            )
            return []
        except requests.exceptions.RequestException as e:
            print(f"  ❌ DOI resolution failed for {doi}: {e}")
            self._log_paywalled(
                pmid, f"DOI resolution failed: {doi}", f"https://doi.org/{doi}",
                classification="PAYWALLED"
            )
            return []

        # Route to the specific scraper based on the resolved domain
        if "nature.com" in domain:
            return scraper.scrape_nature_supplements(response.text, final_url)
        elif (
            "springer.com" in domain
            or "biomedcentral.com" in domain
            or "springeropen.com" in domain
        ):
            # Springer/BMC is a critical publisher - many variants are in their supplements
            return scraper.scrape_springer_supplements(response.text, final_url)
        elif "oup.com" in domain or "academic.oup.com" in domain:
            # Oxford Academic - another major publisher with supplements
            return scraper.scrape_oxford_supplements(response.text, final_url)
        elif "wiley.com" in domain or "onlinelibrary.wiley.com" in domain:
            # Wiley Online Library - major publisher for genetics/medical journals
            return scraper.scrape_wiley_supplements(response.text, final_url)
        elif "karger.com" in domain:
            return scraper.scrape_karger_supplements(response.text, final_url)
        elif (
            "gimjournal.org" in domain
            or "sciencedirect.com" in domain
            or "linkinghub.elsevier.com" in domain
            or "elsevier.com" in domain
        ):
            # For linkinghub URLs, try to follow one more redirect to get to the actual article page
            if "linkinghub.elsevier.com" in domain:
                try:
                    # Try to extract PII from URL and construct ScienceDirect URL
                    # linkinghub URLs often have format: /retrieve/pii/S1547527109005682
                    pii_match = re.search(r"/pii/([^/?]+)", final_url)
                    if pii_match:
                        pii = pii_match.group(1)
                        sciencedirect_url = (
                            f"https://www.sciencedirect.com/science/article/pii/{pii}"
                        )
                        print(
                            f"  → Attempting to access ScienceDirect page: {sciencedirect_url}"
                        )
                        redirect_response = self.session.get(
                            sciencedirect_url, allow_redirects=True, timeout=HTTP_TIMEOUT_DEFAULT
                        )
                        redirect_response.raise_for_status()
                        if redirect_response.url != final_url:
                            print(f"  → Following redirect to: {redirect_response.url}")
                            final_url = redirect_response.url
                            response = redirect_response
                    else:
                        # Fallback: try following redirect from linkinghub
                        redirect_response = self.session.get(
                            final_url, allow_redirects=True, timeout=HTTP_TIMEOUT_DEFAULT
                        )
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
            print(
                f"  - No specific scraper for domain: {domain}. Using generic scraper."
            )
            return scraper.scrape_generic_supplements(response.text, final_url)

    def resolve_and_fetch_fulltext(
        self, doi: str, pmid: str, scraper
    ) -> Tuple[Optional[str], Optional[str], List[Dict]]:
        """
        Resolves a DOI to its final URL and fetches full text + supplements.

        This is used for free articles without PMCIDs. It:
        1. Resolves the DOI to the publisher's article page (with retry on 403)
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
            # Uses retry logic with User-Agent rotation for 403 errors
            print(f"  Resolving DOI for full text: https://doi.org/{doi}")
            encoded_doi = encode_doi_for_url(doi)
            response = self._resolve_with_retry(
                f"https://doi.org/{encoded_doi}", f"DOI {doi} (full text)"
            )

            if not response:
                print(f"  ❌ DOI resolution failed for {doi}: No response received")
                self._log_paywalled(
                    pmid,
                    f"DOI resolution failed (free text): {doi}",
                    f"https://doi.org/{doi}",
                    classification="PAYWALLED"
                )
                return None, None, []

            final_url = response.url
            domain = urlparse(final_url).netloc
            print(f"  ✓ DOI resolved to: {final_url}")

        except requests.exceptions.HTTPError as e:
            status_code = e.response.status_code if e.response else "unknown"
            if status_code == 403:
                print(
                    f"  ❌ DOI resolution blocked (403 Forbidden) for {doi} after {self.MAX_RETRIES} attempts - publisher may require authentication"
                )
            else:
                print(f"  ❌ DOI resolution failed for {doi}: HTTP {status_code}")
            self._log_paywalled(
                pmid,
                f"DOI resolution failed (free text, HTTP {status_code}): {doi}",
                f"https://doi.org/{doi}",
                classification="CAPTCHA_BLOCKED" if status_code == 403 else "PAYWALLED"
            )
            return None, None, []
        except requests.exceptions.RequestException as e:
            print(f"  ❌ DOI resolution failed for {doi}: {e}")
            self._log_paywalled(
                pmid,
                f"DOI resolution failed (free text): {doi}",
                f"https://doi.org/{doi}",
                classification="PAYWALLED"
            )
            return None, None, []

        # Handle Elsevier linkinghub redirects
        if "linkinghub.elsevier.com" in domain:
            try:
                pii_match = re.search(r"/pii/([^/?]+)", final_url)
                if pii_match:
                    pii = pii_match.group(1)
                    sciencedirect_url = (
                        f"https://www.sciencedirect.com/science/article/pii/{pii}"
                    )
                    print(
                        f"  → Attempting to access ScienceDirect page: {sciencedirect_url}"
                    )
                    redirect_response = self.session.get(
                        sciencedirect_url, allow_redirects=True, timeout=HTTP_TIMEOUT_DEFAULT
                    )
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
            print(
                f"  ✓ Extracted full text from publisher page ({len(markdown_content)} characters)"
            )
        else:
            print("  ❌ Could not extract full text from publisher page")
            self._log_paywalled(
                pmid, "Full text extraction failed from publisher", final_url,
                classification="PAYWALLED"
            )

        # Also get supplements
        if "nature.com" in domain:
            supplements = scraper.scrape_nature_supplements(html_content, final_url)
        elif (
            "springer.com" in domain
            or "biomedcentral.com" in domain
            or "springeropen.com" in domain
        ):
            supplements = scraper.scrape_springer_supplements(html_content, final_url)
        elif "oup.com" in domain:
            supplements = scraper.scrape_oxford_supplements(html_content, final_url)
        elif "wiley.com" in domain or "onlinelibrary.wiley.com" in domain:
            supplements = scraper.scrape_wiley_supplements(html_content, final_url)
        elif "karger.com" in domain:
            supplements = scraper.scrape_karger_supplements(html_content, final_url)
        elif any(
            d in domain for d in ["gimjournal.org", "sciencedirect.com", "elsevier.com"]
        ):
            supplements = scraper.scrape_elsevier_supplements(html_content, final_url)
        else:
            supplements = scraper.scrape_generic_supplements(html_content, final_url)

        return markdown_content, final_url, supplements

    def resolve_doi_url(
        self, doi: str, pmid: str
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Resolve a DOI to its final URL without fetching full content.

        Uses retry logic with User-Agent rotation to handle 403 errors.

        Args:
            doi: Digital Object Identifier
            pmid: PubMed ID (for logging)

        Returns:
            Tuple of (final_url, html_content) or (None, None) if resolution failed
        """
        try:
            print(f"  Resolving DOI: https://doi.org/{doi}")
            encoded_doi = encode_doi_for_url(doi)
            response = self._resolve_with_retry(
                f"https://doi.org/{encoded_doi}", f"DOI {doi}"
            )

            if not response:
                print(f"  ❌ DOI resolution failed for {doi}: No response received")
                self._log_paywalled(
                    pmid, f"DOI resolution failed: {doi}", f"https://doi.org/{doi}",
                    classification="PAYWALLED"
                )
                return None, None

            final_url = response.url
            print(f"  ✓ DOI resolved to: {final_url}")

            return final_url, response.text

        except requests.exceptions.HTTPError as e:
            status_code = e.response.status_code if e.response else "unknown"
            print(f"  ❌ DOI resolution failed for {doi}: HTTP {status_code}")
            self._log_paywalled(
                pmid,
                f"DOI resolution failed (HTTP {status_code}): {doi}",
                f"https://doi.org/{doi}",
                classification="CAPTCHA_BLOCKED" if status_code == 403 else "PAYWALLED"
            )
            return None, None
        except requests.exceptions.RequestException as e:
            print(f"  ❌ DOI resolution failed for {doi}: {e}")
            self._log_paywalled(
                pmid, f"DOI resolution failed: {doi}", f"https://doi.org/{doi}",
                classification="PAYWALLED"
            )
            return None, None
