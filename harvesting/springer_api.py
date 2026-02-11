"""
Springer Nature API Module

Provides access to Springer Nature (including Nature, BMC, Springer journals)
full-text content via the official OpenAccess API.

This is used as a preferred method for fetching Springer/Nature/BMC articles when
an API key is available, before falling back to web scraping.

API Documentation: https://dev.springernature.com/
Note: Free API keys only provide access to OpenAccess content.
"""

import logging
import re
import time
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import quote
from xml.etree import ElementTree as ET

import requests
from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)


# DOI prefixes known to be Springer Nature
SPRINGER_DOI_PREFIXES = (
    "10.1007/",  # Springer main
    "10.1038/",  # Nature
    "10.1186/",  # BMC (BioMed Central)
    "10.1140/",  # European Physical Journal
    "10.1057/",  # Palgrave Macmillan
    "10.1023/",  # Springer (older)
    "10.1134/",  # Pleiades (Springer)
    "10.1365/",  # Springer Science
    "10.26508/",  # Some Springer journals
)

# Domains that indicate Springer Nature publisher
SPRINGER_DOMAINS = (
    "springer.com",
    "link.springer.com",
    "nature.com",
    "biomedcentral.com",
    "bmcgenomics.biomedcentral.com",
    "springeropen.com",
)


class SpringerAPIClient:
    """Client for the Springer Nature OpenAccess API."""

    # OpenAccess API endpoint (works with free API keys)
    OPENACCESS_URL = "https://api.springernature.com/openaccess/json"
    OPENACCESS_JATS_URL = "https://api.springernature.com/openaccess/jats"

    # Meta API endpoint (requires higher-tier API key)
    META_URL = "https://api.springernature.com/meta/v2/json"

    def __init__(
        self, api_key: Optional[str] = None, session: Optional[requests.Session] = None
    ):
        """
        Initialize the Springer Nature API client.

        Args:
            api_key: Springer Nature API key (from dev.springernature.com)
            session: Optional requests session to use (for connection pooling)
        """
        self.api_key = api_key
        self.session = session or requests.Session()
        self._last_request_time = 0
        self._min_request_interval = 0.5  # Rate limiting: max 2 req/sec

    @property
    def is_available(self) -> bool:
        """Check if the API client is configured with a valid API key."""
        return bool(self.api_key and self.api_key.strip())

    @staticmethod
    def is_springer_doi(doi: str) -> bool:
        """
        Check if a DOI belongs to a Springer Nature publication.

        Args:
            doi: Digital Object Identifier

        Returns:
            True if the DOI is from Springer Nature, False otherwise
        """
        if not doi:
            return False
        return doi.lower().startswith(SPRINGER_DOI_PREFIXES)

    @staticmethod
    def is_springer_url(url: str) -> bool:
        """
        Check if a URL is from a Springer Nature domain.

        Args:
            url: URL to check

        Returns:
            True if the URL is from Springer Nature, False otherwise
        """
        if not url:
            return False
        from urllib.parse import urlparse

        domain = urlparse(url).netloc.lower()
        return any(springer_domain in domain for springer_domain in SPRINGER_DOMAINS)

    def _rate_limit(self):
        """Enforce rate limiting between API requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()

    def search_openaccess(
        self, query: str, max_results: int = 10
    ) -> Tuple[Optional[List[Dict[str, Any]]], Optional[str]]:
        """
        Search for articles in the Springer Nature OpenAccess collection.

        Args:
            query: Search query (e.g., "KCNH2 mutation")
            max_results: Maximum number of results to return

        Returns:
            Tuple of (records_list, error_message)
        """
        if not self.is_available:
            return None, "Springer API key not configured"

        self._rate_limit()

        params = {
            "q": query,
            "api_key": self.api_key,
            "p": max_results,
        }

        try:
            logger.info(f"Searching Springer OpenAccess for: {query}")
            response = self.session.get(self.OPENACCESS_URL, params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()
                records = data.get("records", [])
                return records, None
            elif response.status_code == 401:
                return None, "Invalid or unauthorized API key"
            elif response.status_code == 429:
                return None, "Rate limit exceeded"
            else:
                return None, f"HTTP {response.status_code}: {response.reason}"

        except requests.exceptions.Timeout:
            return None, "Request timed out"
        except requests.exceptions.RequestException as e:
            return None, f"Request failed: {str(e)}"

    def get_fulltext_by_doi(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Fetch full-text content from Springer OpenAccess API using DOI.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (content, error_message)
            - content: Full-text content (JATS XML) if successful, None otherwise
            - error_message: Error description if failed, None otherwise
        """
        if not self.is_available:
            return None, "Springer API key not configured"

        self._rate_limit()

        # Query the OpenAccess JATS endpoint
        params = {
            "q": f"doi:{doi}",
            "api_key": self.api_key,
        }

        try:
            logger.info(f"Fetching full text from Springer OpenAccess for DOI: {doi}")
            response = self.session.get(
                self.OPENACCESS_JATS_URL, params=params, timeout=30
            )

            if response.status_code == 200:
                content = response.text
                # Check if we got actual content
                if len(content) > 500 and "<article" in content.lower():
                    return content, None
                else:
                    return None, "Article not available in OpenAccess"
            elif response.status_code == 401:
                return None, "Invalid or unauthorized API key"
            elif response.status_code == 403:
                return None, "Access forbidden - article may not be OpenAccess"
            elif response.status_code == 404:
                return None, "Article not found"
            elif response.status_code == 429:
                return None, "Rate limit exceeded"
            else:
                return None, f"HTTP {response.status_code}: {response.reason}"

        except requests.exceptions.Timeout:
            return None, "Request timed out"
        except requests.exceptions.RequestException as e:
            return None, f"Request failed: {str(e)}"

    def get_metadata_by_doi(self, doi: str) -> Tuple[Optional[Dict], Optional[str]]:
        """
        Fetch article metadata from Springer OpenAccess API using DOI.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (metadata_dict, error_message)
        """
        if not self.is_available:
            return None, "Springer API key not configured"

        self._rate_limit()

        params = {
            "q": f"doi:{doi}",
            "api_key": self.api_key,
            "p": 1,
        }

        try:
            logger.info(f"Fetching metadata from Springer for DOI: {doi}")
            response = self.session.get(self.OPENACCESS_URL, params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()
                records = data.get("records", [])
                if records:
                    return records[0], None
                return None, "Article not found in OpenAccess"
            elif response.status_code == 401:
                return None, "Invalid or unauthorized API key"
            else:
                return None, f"HTTP {response.status_code}"

        except requests.exceptions.Timeout:
            return None, "Request timed out"
        except requests.exceptions.RequestException as e:
            return None, f"Request failed: {str(e)}"

    def jats_to_markdown(self, jats_content: str) -> Optional[str]:
        """
        Convert JATS XML to markdown format.

        JATS (Journal Article Tag Suite) is the standard XML format for
        scientific articles. This method extracts the article content and
        converts it to markdown suitable for LLM processing.

        Args:
            jats_content: Raw JATS XML from Springer API

        Returns:
            Markdown-formatted article content, or None if parsing fails
        """
        try:
            # Parse the JATS XML
            root = ET.fromstring(jats_content)

            # Define namespace (JATS often uses namespaces)
            ns = {}
            # Try to extract namespace from root
            if root.tag.startswith("{"):
                ns_end = root.tag.find("}")
                ns["jats"] = root.tag[1:ns_end]

            parts = []

            # Extract title
            title_elem = root.find(".//article-title", ns) or root.find(
                ".//title-group/article-title"
            )
            if title_elem is not None and title_elem.text:
                parts.append(f"# {title_elem.text.strip()}\n")

            # Extract abstract
            abstract_elem = root.find(".//abstract", ns) or root.find(".//abstract/p")
            if abstract_elem is not None:
                abstract_text = "".join(abstract_elem.itertext()).strip()
                if abstract_text:
                    parts.append(f"## Abstract\n\n{abstract_text}\n")

            # Extract body sections
            body = root.find(".//body", ns)
            if body is not None:
                for section in body.findall(".//sec", ns) or body.findall(".//sec"):
                    # Section title
                    sec_title = section.find("title")
                    if sec_title is not None and sec_title.text:
                        parts.append(f"\n## {sec_title.text.strip()}\n")

                    # Section paragraphs
                    for para in section.findall("p"):
                        para_text = "".join(para.itertext()).strip()
                        if para_text:
                            parts.append(f"\n{para_text}\n")

            # If body parsing didn't yield much, try to get all text
            if len(parts) < 3:
                all_text = "".join(root.itertext())
                if len(all_text) > 500:
                    # Clean up whitespace
                    all_text = re.sub(r"\s+", " ", all_text).strip()
                    parts.append(f"\n{all_text}\n")

            if parts:
                return "\n".join(parts)
            return None

        except ET.ParseError as e:
            logger.warning(f"Failed to parse JATS XML: {e}")
            # Try with BeautifulSoup as fallback
            try:
                soup = BeautifulSoup(jats_content, "xml")
                text = soup.get_text(separator="\n\n", strip=True)
                if len(text) > 500:
                    return text
            except Exception:
                pass
            return None
        except Exception as e:
            logger.warning(f"Error converting JATS to markdown: {e}")
            return None

    def fetch_article(
        self, doi: str
    ) -> Tuple[Optional[str], Optional[Dict], Optional[str]]:
        """
        High-level method to fetch an article's full text and metadata.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (markdown_content, metadata, error_message)
        """
        # Get full text
        jats_content, error = self.get_fulltext_by_doi(doi)

        if error:
            # Try to at least get metadata
            metadata, meta_error = self.get_metadata_by_doi(doi)
            if metadata:
                # Return metadata even without full text
                return None, metadata, f"Full text unavailable: {error}"
            return None, None, error

        # Convert JATS to markdown
        markdown = self.jats_to_markdown(jats_content)

        # Get metadata
        metadata, _ = self.get_metadata_by_doi(doi)

        return markdown, metadata, None


def get_springer_client(api_key: Optional[str] = None) -> SpringerAPIClient:
    """
    Factory function to create a Springer API client.

    If no API key is provided, attempts to read from environment variable.

    Args:
        api_key: Optional API key (if not provided, reads from SPRINGER_API_KEY env var)

    Returns:
        Configured SpringerAPIClient instance
    """
    import os

    if api_key is None:
        api_key = os.environ.get("SPRINGER_API_KEY")

    return SpringerAPIClient(api_key=api_key)
