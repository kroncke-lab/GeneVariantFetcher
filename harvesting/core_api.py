"""
CORE API Module - Access to 200M+ open access papers.

CORE (core.ac.uk) is the world's largest aggregator of open access research.
Free API access with registration.

API Documentation: https://core.ac.uk/documentation/api
"""

import logging
import time
from typing import Optional, Dict, Any, Tuple, List

import requests

logger = logging.getLogger(__name__)


class COREAPIClient:
    """Client for the CORE API."""

    BASE_URL = "https://api.core.ac.uk/v3"

    def __init__(
        self,
        api_key: Optional[str] = None,
        session: Optional[requests.Session] = None,
    ):
        """
        Initialize CORE API client.

        Args:
            api_key: CORE API key (get free at https://core.ac.uk/api-keys/register)
            session: Optional requests session
        """
        self.api_key = api_key
        self.session = session or requests.Session()
        self._last_request_time = 0
        self._min_request_interval = 1.0  # Rate limit: 1 req/sec for free tier
        
        if api_key:
            self.session.headers.update({"Authorization": f"Bearer {api_key}"})

    @property
    def is_available(self) -> bool:
        """Check if API key is configured."""
        return bool(self.api_key and self.api_key.strip())

    def _rate_limit(self):
        """Enforce rate limiting."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()

    def search_by_doi(self, doi: str) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
        """
        Search for a paper by DOI.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (result_dict, error_message)
        """
        if not self.is_available:
            return None, "CORE API key not configured"

        self._rate_limit()

        try:
            # CORE uses DOI search
            url = f"{self.BASE_URL}/search/works"
            params = {
                "q": f"doi:{doi}",
                "limit": 1,
            }

            response = self.session.get(url, params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()
                results = data.get("results", [])
                if results:
                    return results[0], None
                return None, "DOI not found in CORE"
            elif response.status_code == 401:
                return None, "Invalid CORE API key"
            elif response.status_code == 429:
                return None, "Rate limit exceeded"
            else:
                return None, f"HTTP {response.status_code}"

        except requests.exceptions.Timeout:
            return None, "Request timed out"
        except Exception as e:
            return None, f"Request failed: {str(e)}"

    def get_fulltext_by_id(self, core_id: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Get full text of a paper by CORE ID.

        Args:
            core_id: CORE internal ID

        Returns:
            Tuple of (fulltext, error_message)
        """
        if not self.is_available:
            return None, "CORE API key not configured"

        self._rate_limit()

        try:
            url = f"{self.BASE_URL}/works/{core_id}"
            response = self.session.get(url, timeout=30)

            if response.status_code == 200:
                data = response.json()
                fulltext = data.get("fullText")
                if fulltext:
                    return fulltext, None
                # Try download URL
                download_url = data.get("downloadUrl")
                if download_url:
                    return None, f"Full text not available, but download URL: {download_url}"
                return None, "No full text available"
            else:
                return None, f"HTTP {response.status_code}"

        except Exception as e:
            return None, f"Request failed: {str(e)}"

    def search_papers(
        self,
        query: str,
        limit: int = 10,
        fulltext_only: bool = True,
    ) -> Tuple[Optional[List[Dict]], Optional[str]]:
        """
        Search for papers by query.

        Args:
            query: Search query (e.g., "KCNH2 variant")
            limit: Maximum results to return
            fulltext_only: Only return papers with full text available

        Returns:
            Tuple of (results_list, error_message)
        """
        if not self.is_available:
            return None, "CORE API key not configured"

        self._rate_limit()

        try:
            url = f"{self.BASE_URL}/search/works"
            params = {
                "q": query,
                "limit": limit,
            }
            
            if fulltext_only:
                params["q"] += " AND _exists_:fullText"

            response = self.session.get(url, params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()
                return data.get("results", []), None
            else:
                return None, f"HTTP {response.status_code}"

        except Exception as e:
            return None, f"Request failed: {str(e)}"

    def find_fulltext_for_doi(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Try to get full text for a DOI via CORE.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (fulltext_content, error_message)
        """
        # First search for the paper
        result, error = self.search_by_doi(doi)
        if error:
            return None, error
        if not result:
            return None, "Paper not found in CORE"

        # Check if full text is directly available
        fulltext = result.get("fullText")
        if fulltext and len(fulltext) > 500:
            logger.info(f"Found full text in CORE for DOI: {doi}")
            return fulltext, None

        # Try to get by CORE ID
        core_id = result.get("id")
        if core_id:
            fulltext, error = self.get_fulltext_by_id(str(core_id))
            if fulltext:
                return fulltext, None

        # Check for download URL
        download_url = result.get("downloadUrl")
        if download_url:
            return None, f"Full text not indexed, but available at: {download_url}"

        return None, "Full text not available in CORE"


def get_core_client(api_key: Optional[str] = None) -> COREAPIClient:
    """
    Factory function to create CORE API client.

    Args:
        api_key: API key (if not provided, reads from CORE_API_KEY env var)

    Returns:
        Configured COREAPIClient
    """
    import os

    if api_key is None:
        api_key = os.environ.get("CORE_API_KEY")

    return COREAPIClient(api_key=api_key)
