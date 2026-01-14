"""
Wiley API Module

Provides access to Wiley full-text content via the official TDM (Text and Data Mining) API.
This is used as a preferred method for fetching Wiley articles when an API key
is available, before falling back to web scraping.

API Documentation: https://onlinelibrary.wiley.com/library-info/resources/text-and-datamining
"""

import logging
import re
import time
from typing import Optional, Tuple
from urllib.parse import urlparse, quote
from xml.etree import ElementTree as ET

import requests

logger = logging.getLogger(__name__)


# DOI prefixes known to be Wiley
WILEY_DOI_PREFIXES = (
    "10.1002/",   # Wiley main prefix
    "10.1111/",   # Wiley-Blackwell
    "10.1113/",   # The Physiological Society (Wiley)
    "10.1096/",   # FASEB Journal (Wiley)
    "10.1634/",   # Stem Cells (Wiley)
    "10.1111/j.", # Wiley journal articles
)

# Domains that indicate Wiley publisher
WILEY_DOMAINS = (
    "wiley.com",
    "onlinelibrary.wiley.com",
    "wiley-vch.de",
    "fasebj.org",
)


class WileyAPIClient:
    """Client for the Wiley TDM (Text and Data Mining) API."""

    BASE_URL = "https://api.wiley.com/onlinelibrary/tdm/v1/articles"

    def __init__(self, api_key: Optional[str] = None, session: Optional[requests.Session] = None):
        """
        Initialize the Wiley API client.

        Args:
            api_key: Wiley TDM API key
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
    def is_wiley_doi(doi: str) -> bool:
        """
        Check if a DOI belongs to a Wiley publication.

        Args:
            doi: Digital Object Identifier

        Returns:
            True if the DOI is from Wiley, False otherwise
        """
        if not doi:
            return False
        return doi.lower().startswith(WILEY_DOI_PREFIXES)

    @staticmethod
    def is_wiley_url(url: str) -> bool:
        """
        Check if a URL is from a Wiley domain.

        Args:
            url: URL to check

        Returns:
            True if the URL is from Wiley, False otherwise
        """
        if not url:
            return False
        domain = urlparse(url).netloc.lower()
        return any(wiley_domain in domain for wiley_domain in WILEY_DOMAINS)

    @staticmethod
    def extract_doi_from_url(url: str) -> Optional[str]:
        """
        Extract DOI from a Wiley Online Library URL.

        Args:
            url: Wiley URL

        Returns:
            DOI string if found, None otherwise
        """
        # Pattern: /doi/full/10.1002/xxx or /doi/10.1002/xxx or /doi/abs/10.1002/xxx
        doi_match = re.search(r'/doi/(?:full/|abs/|epdf/|pdf/)?(\d+\.\d+/[^\s?#]+)', url, re.IGNORECASE)
        if doi_match:
            return doi_match.group(1)
        return None

    def _rate_limit(self):
        """Enforce rate limiting between API requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()

    def get_fulltext_by_doi(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Fetch full-text content from Wiley TDM API using DOI.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (content, error_message)
            - content: Full-text XML/HTML if successful, None otherwise
            - error_message: Error description if failed, None otherwise
        """
        if not self.is_available:
            return None, "Wiley API key not configured"

        self._rate_limit()

        # URL-encode the DOI (can contain special characters like <, >, parentheses)
        encoded_doi = quote(doi, safe='/:')
        url = f"{self.BASE_URL}/{encoded_doi}"
        headers = {
            "Wiley-TDM-Client-Token": self.api_key,
            "Accept": "application/xml, text/xml, application/xhtml+xml, text/html",
        }

        try:
            logger.info(f"Fetching full text from Wiley API for DOI: {doi}")
            response = self.session.get(url, headers=headers, timeout=30)

            if response.status_code == 200:
                return response.text, None
            elif response.status_code == 401:
                return None, "Invalid or unauthorized API key"
            elif response.status_code == 403:
                return None, "Access forbidden - API key may lack permissions or article not available"
            elif response.status_code == 404:
                return None, "Article not found via Wiley API"
            elif response.status_code == 429:
                return None, "Rate limit exceeded"
            else:
                return None, f"HTTP {response.status_code}: {response.reason}"

        except requests.exceptions.Timeout:
            return None, "Request timed out"
        except requests.exceptions.RequestException as e:
            return None, f"Request failed: {str(e)}"

    def content_to_markdown(self, content: str) -> Optional[str]:
        """
        Convert Wiley full-text content (XML or HTML) to markdown format.

        The Wiley TDM API can return content in various formats. This method
        extracts the article content and converts it to markdown suitable for
        LLM processing.

        Args:
            content: Raw content from Wiley API (XML or HTML)

        Returns:
            Markdown string if conversion successful, None otherwise
        """
        if not content:
            return None

        # Try XML parsing first
        markdown = self._parse_xml_content(content)
        if markdown and len(markdown) > 500:
            return markdown

        # Fall back to HTML parsing
        markdown = self._parse_html_content(content)
        if markdown and len(markdown) > 200:
            return markdown

        return None

    def _parse_xml_content(self, xml_content: str) -> Optional[str]:
        """
        Parse Wiley XML content to markdown.

        Args:
            xml_content: Raw XML from Wiley API

        Returns:
            Markdown string if parsing successful, None otherwise
        """
        try:
            # Remove XML namespace prefixes for easier parsing
            xml_content = re.sub(r'xmlns[^=]*="[^"]*"', '', xml_content)
            xml_content = re.sub(r'<([a-z]+):', r'<', xml_content)
            xml_content = re.sub(r'</([a-z]+):', r'</', xml_content)

            root = ET.fromstring(xml_content)
            markdown = "# MAIN TEXT\n\n"

            # Extract title
            title_elem = root.find('.//article-title')
            if title_elem is None:
                title_elem = root.find('.//title')
            if title_elem is None:
                title_elem = root.find('.//{*}title')
            if title_elem is not None:
                title_text = self._extract_text(title_elem)
                if title_text:
                    markdown += f"## {title_text.strip()}\n\n"

            # Extract abstract
            abstract_elem = root.find('.//abstract')
            if abstract_elem is None:
                abstract_elem = root.find('.//{*}abstract')
            if abstract_elem is not None:
                abstract_text = self._extract_text(abstract_elem)
                if abstract_text:
                    markdown += f"### Abstract\n\n{abstract_text}\n\n"

            # Extract body sections
            body_elem = root.find('.//body')
            if body_elem is None:
                body_elem = root.find('.//{*}body')

            if body_elem is not None:
                # Process sections
                sections = body_elem.findall('.//sec')
                if not sections:
                    sections = body_elem.findall('.//{*}sec')
                if not sections:
                    sections = body_elem.findall('.//section')
                if not sections:
                    sections = body_elem.findall('.//{*}section')

                if sections:
                    for section in sections:
                        section_md = self._process_section(section)
                        if section_md:
                            markdown += section_md
                else:
                    # No sections found, extract all paragraphs from body
                    paragraphs = body_elem.findall('.//p')
                    if not paragraphs:
                        paragraphs = body_elem.findall('.//{*}p')
                    for para in paragraphs:
                        para_text = self._extract_text(para)
                        if para_text:
                            markdown += f"{para_text}\n\n"

            return markdown if len(markdown) > 200 else None

        except ET.ParseError as e:
            logger.debug(f"Failed to parse Wiley XML: {e}")
            return None
        except Exception as e:
            logger.debug(f"Error converting Wiley XML to markdown: {e}")
            return None

    def _parse_html_content(self, html_content: str) -> Optional[str]:
        """
        Parse Wiley HTML content to markdown.

        Args:
            html_content: Raw HTML from Wiley API

        Returns:
            Markdown string if parsing successful, None otherwise
        """
        try:
            from bs4 import BeautifulSoup
        except ImportError:
            logger.warning("BeautifulSoup not available for HTML parsing")
            return None

        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            markdown = "# MAIN TEXT\n\n"

            # Extract title
            title = soup.find('h1', class_='citation__title')
            if not title:
                title = soup.find('h1')
            if title:
                markdown += f"## {title.get_text(strip=True)}\n\n"

            # Extract abstract
            abstract = soup.find('section', class_='article-section__abstract')
            if not abstract:
                abstract = soup.find('div', class_='abstract')
            if abstract:
                abstract_text = abstract.get_text(separator=' ', strip=True)
                markdown += f"### Abstract\n\n{abstract_text}\n\n"

            # Extract main content sections
            content_sections = soup.find_all('section', class_='article-section__content')
            if not content_sections:
                content_sections = soup.find_all('div', class_='article__body')

            for section in content_sections:
                # Get section heading
                heading = section.find(['h2', 'h3', 'h4'])
                if heading:
                    markdown += f"### {heading.get_text(strip=True)}\n\n"

                # Get paragraphs
                paragraphs = section.find_all('p')
                for para in paragraphs:
                    para_text = para.get_text(separator=' ', strip=True)
                    if para_text and len(para_text) > 20:
                        markdown += f"{para_text}\n\n"

            # Fallback: extract all paragraphs if sections not found
            if len(markdown) < 500:
                all_paragraphs = soup.find_all('p')
                for para in all_paragraphs:
                    para_text = para.get_text(separator=' ', strip=True)
                    if para_text and len(para_text) > 50:
                        markdown += f"{para_text}\n\n"

            return markdown if len(markdown) > 200 else None

        except Exception as e:
            logger.debug(f"Error parsing Wiley HTML: {e}")
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
        title_elem = section_elem.find('./title')
        if title_elem is None:
            title_elem = section_elem.find('./{*}title')

        if title_elem is not None:
            title_text = self._extract_text(title_elem)
            if title_text:
                heading_prefix = "#" * min(level, 6)
                markdown += f"{heading_prefix} {title_text.strip()}\n\n"

        # Extract paragraphs
        paragraphs = section_elem.findall('./p')
        if not paragraphs:
            paragraphs = section_elem.findall('./{*}p')

        for para in paragraphs:
            para_text = self._extract_text(para)
            if para_text:
                markdown += f"{para_text}\n\n"

        # Process nested sections
        nested_sections = section_elem.findall('./sec')
        if not nested_sections:
            nested_sections = section_elem.findall('./{*}sec')
        if not nested_sections:
            nested_sections = section_elem.findall('./section')

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

    def fetch_fulltext(self, doi: Optional[str] = None,
                       url: Optional[str] = None) -> Tuple[Optional[str], Optional[str]]:
        """
        Fetch full-text content and convert to markdown.

        Tries to fetch using DOI. If a URL is provided and no DOI,
        extracts DOI from the URL.

        Args:
            doi: Digital Object Identifier
            url: URL (used to extract DOI if not provided)

        Returns:
            Tuple of (markdown_content, error_message)
            - markdown_content: Markdown text if successful, None otherwise
            - error_message: Error description if failed, None otherwise
        """
        if not self.is_available:
            return None, "Wiley API key not configured"

        # Extract DOI from URL if not provided
        if not doi and url:
            doi = self.extract_doi_from_url(url)

        if not doi:
            return None, "No valid DOI available"

        # Fetch content
        content, error = self.get_fulltext_by_doi(doi)
        if content:
            markdown = self.content_to_markdown(content)
            # Check for sufficient content - abstracts are typically 200-400 words
            # Full articles should be at least 2000 characters (roughly 350+ words)
            # Also check for section headings beyond just Abstract
            MIN_FULLTEXT_LENGTH = 2000
            has_body_sections = markdown and any(
                section in markdown.lower()
                for section in ['### introduction', '### methods', '### results',
                                '### discussion', '### materials', '### content']
            )
            if markdown and (len(markdown) > MIN_FULLTEXT_LENGTH or has_body_sections):
                logger.info(f"Successfully fetched Wiley article via DOI: {doi}")
                return markdown, None
            elif markdown and len(markdown) > 500:
                # Got some content but likely just abstract
                error = "Content conversion produced insufficient content (abstract only)"
            elif not error:
                error = "Content conversion produced insufficient content"

        return None, error or "Failed to fetch content"
