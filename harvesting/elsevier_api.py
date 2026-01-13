"""
Elsevier API Module

Provides access to Elsevier/ScienceDirect full-text content via the official API.
This is used as a preferred method for fetching Elsevier articles when an API key
is available, before falling back to web scraping.

API Documentation: https://dev.elsevier.com/documentation/ArticleRetrievalAPI.wadl
"""

import logging
import re
import time
from typing import Optional, Tuple
from urllib.parse import urlparse
from xml.etree import ElementTree as ET

import requests

logger = logging.getLogger(__name__)


# DOI prefixes known to be Elsevier
ELSEVIER_DOI_PREFIXES = (
    "10.1016/",   # ScienceDirect (main)
    "10.1053/",   # Some Elsevier journals
    "10.1067/",   # Mosby (Elsevier)
    "10.1054/",   # Academic Press (Elsevier)
    "10.1006/",   # Academic Press (Elsevier)
    "10.1038/",   # Nature (some under Elsevier agreement)
)

# Domains that indicate Elsevier publisher
ELSEVIER_DOMAINS = (
    "sciencedirect.com",
    "elsevier.com",
    "linkinghub.elsevier.com",
    "gimjournal.org",
    "cell.com",
)


class ElsevierAPIClient:
    """Client for the Elsevier Article Retrieval API."""

    BASE_URL = "https://api.elsevier.com/content/article"

    def __init__(self, api_key: Optional[str] = None, session: Optional[requests.Session] = None):
        """
        Initialize the Elsevier API client.

        Args:
            api_key: Elsevier API key (from dev.elsevier.com)
            session: Optional requests session to use (for connection pooling)
        """
        self.api_key = api_key
        self.session = session or requests.Session()
        self._last_request_time = 0
        self._min_request_interval = 0.2  # Rate limiting: max 5 req/sec

    @property
    def is_available(self) -> bool:
        """Check if the API client is configured with a valid API key."""
        return bool(self.api_key and self.api_key.strip())

    @staticmethod
    def is_elsevier_doi(doi: str) -> bool:
        """
        Check if a DOI belongs to an Elsevier publication.

        Args:
            doi: Digital Object Identifier

        Returns:
            True if the DOI is from Elsevier, False otherwise
        """
        if not doi:
            return False
        return doi.lower().startswith(ELSEVIER_DOI_PREFIXES)

    @staticmethod
    def is_elsevier_url(url: str) -> bool:
        """
        Check if a URL is from an Elsevier domain.

        Args:
            url: URL to check

        Returns:
            True if the URL is from Elsevier, False otherwise
        """
        if not url:
            return False
        domain = urlparse(url).netloc.lower()
        return any(elsevier_domain in domain for elsevier_domain in ELSEVIER_DOMAINS)

    @staticmethod
    def extract_pii_from_url(url: str) -> Optional[str]:
        """
        Extract PII (Publisher Item Identifier) from a ScienceDirect URL.

        Args:
            url: ScienceDirect URL

        Returns:
            PII string if found, None otherwise
        """
        # Pattern: /pii/S0000000000000000
        pii_match = re.search(r'/pii/([A-Z0-9]+)', url, re.IGNORECASE)
        if pii_match:
            return pii_match.group(1)
        return None

    def _rate_limit(self):
        """Enforce rate limiting between API requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()

    def get_fulltext_by_doi(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Fetch full-text content from Elsevier API using DOI.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (xml_content, error_message)
            - xml_content: Full-text XML if successful, None otherwise
            - error_message: Error description if failed, None otherwise
        """
        if not self.is_available:
            return None, "Elsevier API key not configured"

        self._rate_limit()

        url = f"{self.BASE_URL}/doi/{doi}"
        headers = {
            "X-ELS-APIKey": self.api_key,
            "Accept": "application/xml",
        }

        try:
            logger.info(f"Fetching full text from Elsevier API for DOI: {doi}")
            response = self.session.get(url, headers=headers, timeout=30)

            if response.status_code == 200:
                return response.text, None
            elif response.status_code == 401:
                return None, "Invalid or unauthorized API key"
            elif response.status_code == 403:
                return None, "Access forbidden - API key may lack permissions or article not available"
            elif response.status_code == 404:
                return None, "Article not found via Elsevier API"
            elif response.status_code == 429:
                return None, "Rate limit exceeded"
            else:
                return None, f"HTTP {response.status_code}: {response.reason}"

        except requests.exceptions.Timeout:
            return None, "Request timed out"
        except requests.exceptions.RequestException as e:
            return None, f"Request failed: {str(e)}"

    def get_fulltext_by_pii(self, pii: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Fetch full-text content from Elsevier API using PII.

        Args:
            pii: Publisher Item Identifier

        Returns:
            Tuple of (xml_content, error_message)
            - xml_content: Full-text XML if successful, None otherwise
            - error_message: Error description if failed, None otherwise
        """
        if not self.is_available:
            return None, "Elsevier API key not configured"

        self._rate_limit()

        url = f"{self.BASE_URL}/pii/{pii}"
        headers = {
            "X-ELS-APIKey": self.api_key,
            "Accept": "application/xml",
        }

        try:
            logger.info(f"Fetching full text from Elsevier API for PII: {pii}")
            response = self.session.get(url, headers=headers, timeout=30)

            if response.status_code == 200:
                return response.text, None
            elif response.status_code == 401:
                return None, "Invalid or unauthorized API key"
            elif response.status_code == 403:
                return None, "Access forbidden - API key may lack permissions or article not available"
            elif response.status_code == 404:
                return None, "Article not found via Elsevier API"
            elif response.status_code == 429:
                return None, "Rate limit exceeded"
            else:
                return None, f"HTTP {response.status_code}: {response.reason}"

        except requests.exceptions.Timeout:
            return None, "Request timed out"
        except requests.exceptions.RequestException as e:
            return None, f"Request failed: {str(e)}"

    def xml_to_markdown(self, xml_content: str) -> Optional[str]:
        """
        Convert Elsevier full-text XML to markdown format.

        The Elsevier API returns XML in their proprietary format. This method
        extracts the article content and converts it to markdown suitable for
        LLM processing.

        Args:
            xml_content: Raw XML from Elsevier API

        Returns:
            Markdown string if conversion successful, None otherwise
        """
        if not xml_content:
            return None

        try:
            # Remove XML namespace prefixes for easier parsing
            xml_content = re.sub(r'xmlns[^=]*="[^"]*"', '', xml_content)
            xml_content = re.sub(r'<([a-z]+):', r'<', xml_content)
            xml_content = re.sub(r'</([a-z]+):', r'</', xml_content)

            root = ET.fromstring(xml_content)
            markdown = "# MAIN TEXT\n\n"

            # Extract title
            title_elem = root.find('.//dc:title', {'dc': 'http://purl.org/dc/elements/1.1/'})
            if title_elem is None:
                title_elem = root.find('.//title')
            if title_elem is None:
                title_elem = root.find('.//{*}title')
            if title_elem is not None and title_elem.text:
                markdown += f"## {title_elem.text.strip()}\n\n"

            # Extract abstract
            abstract_elem = root.find('.//abstract')
            if abstract_elem is None:
                abstract_elem = root.find('.//{*}abstract')
            if abstract_elem is not None:
                abstract_text = self._extract_text(abstract_elem)
                if abstract_text:
                    markdown += f"### Abstract\n\n{abstract_text}\n\n"

            # Extract body/sections
            body_elem = root.find('.//body')
            if body_elem is None:
                body_elem = root.find('.//{*}body')
            if body_elem is None:
                # Try to find raw-text or originalText
                body_elem = root.find('.//rawtext')
                if body_elem is None:
                    body_elem = root.find('.//{*}rawtext')
                if body_elem is None:
                    body_elem = root.find('.//originalText')
                if body_elem is None:
                    body_elem = root.find('.//{*}originalText')

            if body_elem is not None:
                # Process sections
                sections = body_elem.findall('.//section')
                if not sections:
                    sections = body_elem.findall('.//{*}section')
                if not sections:
                    sections = body_elem.findall('.//ce:section')
                if not sections:
                    sections = body_elem.findall('.//sec')

                if sections:
                    for section in sections:
                        section_md = self._process_section(section)
                        if section_md:
                            markdown += section_md
                else:
                    # No sections found, extract all text from body
                    body_text = self._extract_text(body_elem)
                    if body_text:
                        markdown += f"### Content\n\n{body_text}\n\n"

            # If we got very little content, try to extract from full-text-retrieval-response
            if len(markdown) < 500:
                # Try originalText element (common in Elsevier responses)
                original_text = root.find('.//originalText')
                if original_text is None:
                    original_text = root.find('.//{*}originalText')
                if original_text is not None:
                    full_text = self._extract_text(original_text)
                    if full_text and len(full_text) > len(markdown):
                        markdown = f"# MAIN TEXT\n\n{full_text}\n\n"

            return markdown if len(markdown) > 200 else None

        except ET.ParseError as e:
            logger.error(f"Failed to parse Elsevier XML: {e}")
            return None
        except Exception as e:
            logger.error(f"Error converting Elsevier XML to markdown: {e}")
            return None

    def _process_section(self, section_elem: ET.Element, level: int = 3) -> str:
        """
        Process a section element and convert to markdown.

        Args:
            section_elem: XML element representing a section
            level: Heading level for this section (default: 3)

        Returns:
            Markdown string for this section
        """
        markdown = ""

        # Get section title
        title_elem = section_elem.find('./section-title')
        if title_elem is None:
            title_elem = section_elem.find('./title')
        if title_elem is None:
            title_elem = section_elem.find('./{*}section-title')
        if title_elem is None:
            title_elem = section_elem.find('./{*}title')

        if title_elem is not None and title_elem.text:
            heading_prefix = "#" * min(level, 6)
            markdown += f"{heading_prefix} {title_elem.text.strip()}\n\n"

        # Extract paragraphs
        paragraphs = section_elem.findall('./para')
        if not paragraphs:
            paragraphs = section_elem.findall('./{*}para')
        if not paragraphs:
            paragraphs = section_elem.findall('./p')
        if not paragraphs:
            paragraphs = section_elem.findall('./{*}p')

        for para in paragraphs:
            para_text = self._extract_text(para)
            if para_text:
                markdown += f"{para_text}\n\n"

        # Process nested sections
        nested_sections = section_elem.findall('./section')
        if not nested_sections:
            nested_sections = section_elem.findall('./{*}section')
        if not nested_sections:
            nested_sections = section_elem.findall('./sec')

        for nested in nested_sections:
            markdown += self._process_section(nested, level + 1)

        return markdown

    def _extract_text(self, elem: ET.Element) -> str:
        """
        Extract all text content from an XML element, including nested elements.

        Args:
            elem: XML element to extract text from

        Returns:
            Concatenated text content
        """
        if elem is None:
            return ""

        # Get all text including from child elements
        text_parts = []

        if elem.text:
            text_parts.append(elem.text.strip())

        for child in elem:
            child_text = self._extract_text(child)
            if child_text:
                text_parts.append(child_text)
            if child.tail:
                text_parts.append(child.tail.strip())

        return " ".join(filter(None, text_parts))

    def fetch_fulltext(self, doi: Optional[str] = None, pii: Optional[str] = None,
                       url: Optional[str] = None) -> Tuple[Optional[str], Optional[str]]:
        """
        Fetch full-text content and convert to markdown.

        Tries to fetch using DOI first, then PII if DOI fails.
        If a URL is provided, extracts PII from it as fallback.

        Args:
            doi: Digital Object Identifier
            pii: Publisher Item Identifier
            url: URL (used to extract PII if not provided)

        Returns:
            Tuple of (markdown_content, error_message)
            - markdown_content: Markdown text if successful, None otherwise
            - error_message: Error description if failed, None otherwise
        """
        if not self.is_available:
            return None, "Elsevier API key not configured"

        # Try DOI first
        if doi:
            xml_content, error = self.get_fulltext_by_doi(doi)
            if xml_content:
                markdown = self.xml_to_markdown(xml_content)
                if markdown and len(markdown) > 500:
                    logger.info(f"Successfully fetched Elsevier article via DOI: {doi}")
                    return markdown, None
                elif not error:
                    error = "XML conversion produced insufficient content"

        # Try PII
        if not pii and url:
            pii = self.extract_pii_from_url(url)

        if pii:
            xml_content, pii_error = self.get_fulltext_by_pii(pii)
            if xml_content:
                markdown = self.xml_to_markdown(xml_content)
                if markdown and len(markdown) > 500:
                    logger.info(f"Successfully fetched Elsevier article via PII: {pii}")
                    return markdown, None
                elif not pii_error:
                    pii_error = "XML conversion produced insufficient content"
            # Use PII error if DOI wasn't tried or also failed
            if not doi:
                error = pii_error

        return None, error or "No valid DOI or PII available"
