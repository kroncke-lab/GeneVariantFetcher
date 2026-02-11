"""
Unpaywall API Module

Provides access to open access versions of paywalled papers via Unpaywall.
This is a FREE service that finds legal open access versions of papers.

API Documentation: https://unpaywall.org/products/api
No API key required - just provide an email address.

Usage:
    client = UnpaywallClient(email="your@email.com")
    result = client.find_open_access("10.1016/j.example.2024.12345")
    if result:
        pdf_url = result.get('pdf_url')
"""

import logging
import time
from typing import Optional, Dict, Any, Tuple

import requests

logger = logging.getLogger(__name__)


class UnpaywallClient:
    """Client for the Unpaywall API."""

    BASE_URL = "https://api.unpaywall.org/v2"

    def __init__(self, email: str, session: Optional[requests.Session] = None):
        """
        Initialize Unpaywall client.

        Args:
            email: Your email address (required by Unpaywall ToS)
            session: Optional requests session for connection pooling
        """
        self.email = email
        self.session = session or requests.Session()
        self._last_request_time = 0
        self._min_request_interval = 0.1  # Be nice to the free service

    def _rate_limit(self):
        """Enforce rate limiting."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()

    def find_open_access(
        self, doi: str
    ) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
        """
        Find open access version of a paper by DOI.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (result_dict, error_message)
            result_dict contains:
                - pdf_url: Direct link to PDF (if available)
                - is_oa: Whether paper is open access
                - oa_status: Type of OA (gold, green, hybrid, bronze)
                - best_oa_location: Best location details
        """
        if not doi:
            return None, "No DOI provided"

        # Clean DOI
        doi = doi.strip()
        if doi.startswith("https://doi.org/"):
            doi = doi[16:]
        elif doi.startswith("http://doi.org/"):
            doi = doi[15:]
        elif doi.startswith("doi:"):
            doi = doi[4:]

        self._rate_limit()

        url = f"{self.BASE_URL}/{doi}"
        params = {"email": self.email}

        try:
            logger.info(f"Checking Unpaywall for DOI: {doi}")
            response = self.session.get(url, params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()

                # Extract useful info
                result = {
                    "doi": data.get("doi"),
                    "title": data.get("title"),
                    "is_oa": data.get("is_oa", False),
                    "oa_status": data.get("oa_status"),
                    "journal_name": data.get("journal_name"),
                    "publisher": data.get("publisher"),
                    "pdf_url": None,
                    "best_oa_location": None,
                }

                # Find best OA location
                best_location = data.get("best_oa_location")
                if best_location:
                    result["best_oa_location"] = best_location
                    result["pdf_url"] = best_location.get("url_for_pdf")

                    # If no direct PDF, try landing page
                    if not result["pdf_url"]:
                        result["landing_page"] = best_location.get(
                            "url_for_landing_page"
                        )

                # Check all OA locations for PDFs
                if not result["pdf_url"]:
                    for location in data.get("oa_locations", []):
                        pdf_url = location.get("url_for_pdf")
                        if pdf_url:
                            result["pdf_url"] = pdf_url
                            break

                if result["is_oa"]:
                    logger.info(
                        f"Found OA version: {result.get('pdf_url') or result.get('landing_page')}"
                    )
                else:
                    logger.info(f"No OA version found for {doi}")

                return result, None

            elif response.status_code == 404:
                return None, "DOI not found in Unpaywall"
            elif response.status_code == 422:
                return None, "Invalid DOI format"
            else:
                return None, f"HTTP {response.status_code}: {response.reason}"

        except requests.exceptions.Timeout:
            return None, "Request timed out"
        except requests.exceptions.RequestException as e:
            return None, f"Request failed: {str(e)}"
        except Exception as e:
            return None, f"Unexpected error: {str(e)}"

    def download_pdf(
        self, pdf_url: str, output_path: str
    ) -> Tuple[bool, Optional[str]]:
        """
        Download PDF from Unpaywall URL.

        Args:
            pdf_url: URL to PDF
            output_path: Where to save the PDF

        Returns:
            Tuple of (success, error_message)
        """
        if not pdf_url:
            return False, "No PDF URL provided"

        try:
            logger.info(f"Downloading PDF from: {pdf_url}")
            response = self.session.get(
                pdf_url,
                timeout=60,
                headers={
                    "User-Agent": "GeneVariantFetcher/1.0 (mailto:{})".format(
                        self.email
                    )
                },
            )

            if response.status_code == 200:
                content_type = response.headers.get("content-type", "")
                if "pdf" in content_type.lower() or pdf_url.endswith(".pdf"):
                    with open(output_path, "wb") as f:
                        f.write(response.content)
                    logger.info(f"PDF saved to: {output_path}")
                    return True, None
                else:
                    return (
                        False,
                        f"Response was not a PDF (content-type: {content_type})",
                    )
            else:
                return False, f"HTTP {response.status_code}"

        except Exception as e:
            return False, f"Download failed: {str(e)}"


def get_unpaywall_client(email: Optional[str] = None) -> UnpaywallClient:
    """
    Factory function to create Unpaywall client.

    Args:
        email: Email address (if not provided, reads from NCBI_EMAIL env var)

    Returns:
        Configured UnpaywallClient
    """
    import os

    if email is None:
        email = os.environ.get("NCBI_EMAIL", "gvf@example.com")

    return UnpaywallClient(email=email)
