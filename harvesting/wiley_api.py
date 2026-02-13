"""
Wiley API Module

Provides access to Wiley full-text content via the official TDM (Text and Data Mining) API.
This is used as a preferred method for fetching Wiley articles when an API key
is available, before falling back to web scraping.

API Documentation: https://onlinelibrary.wiley.com/library-info/resources/text-and-datamining
"""

import logging
import re
import tempfile
import time
from pathlib import Path
from typing import Optional, Tuple
from urllib.parse import quote, unquote, urlparse
from xml.etree import ElementTree as ET

import requests
from bs4 import BeautifulSoup

from config.constants import HTTP_TIMEOUT_DEFAULT
from utils.http_utils import BROWSER_HEADERS

from .format_converters import FormatConverter

logger = logging.getLogger(__name__)


# DOI prefixes known to be Wiley
WILEY_DOI_PREFIXES = (
    "10.1002/",  # Wiley main prefix
    "10.1111/",  # Wiley-Blackwell
    "10.1113/",  # The Physiological Society (Wiley)
    "10.1096/",  # FASEB Journal (Wiley)
    "10.1634/",  # Stem Cells (Wiley)
    "10.1111/j.",  # Wiley journal articles
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

    def __init__(
        self, api_key: Optional[str] = None, session: Optional[requests.Session] = None
    ):
        """
        Initialize the Wiley API client.

        Args:
            api_key: Wiley TDM API key
            session: Optional requests session to use (for connection pooling)
        """
        self.api_key = api_key
        self.session = session or requests.Session()
        # Ensure browser-like defaults for sessions created outside the orchestrator.
        for key, value in BROWSER_HEADERS.items():
            self.session.headers.setdefault(key, value)
        self._last_request_time = 0
        self._min_request_interval = 0.5  # Rate limiting: max 2 req/sec

    def _get_wiley_headers(self, referer: Optional[str] = None) -> dict:
        headers = dict(self.session.headers)
        if referer:
            headers["Referer"] = referer
        return headers

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
    def encode_doi_for_web(doi: str) -> str:
        """
        Encode a DOI for use in Wiley Online Library URLs.

        For SICI-style DOIs (older format with special characters like <, >, ::, etc.),
        Wiley only requires < and > to be URL-encoded. Other special characters
        like parentheses, colons, and semicolons should remain unencoded to match
        Wiley's actual URL format.

        Example: DOI 10.1002/(SICI)1098-1004(200005)15:5<483::AID-HUMU18>3.0.CO;2-T
        becomes: 10.1002/(SICI)1098-1004(200005)15:5%3C483::AID-HUMU18%3E3.0.CO;2-T

        Args:
            doi: Digital Object Identifier

        Returns:
            URL-encoded DOI string suitable for Wiley web URLs
        """
        # Only encode < and > - these are the only characters that need encoding
        # for Wiley URLs. Parentheses, colons, and semicolons stay unencoded.
        return doi.replace("<", "%3C").replace(">", "%3E")

    @staticmethod
    def extract_doi_from_url(url: str) -> Optional[str]:
        """
        Extract DOI from a Wiley Online Library URL.

        Handles URL-encoded DOIs (e.g., SICI-style DOIs where < and > are
        encoded as %3C and %3E) by decoding them back to raw DOI format.

        Args:
            url: Wiley URL

        Returns:
            DOI string if found, None otherwise (decoded/unescaped)
        """
        # Pattern: /doi/full/10.1002/xxx or /doi/10.1002/xxx or /doi/abs/10.1002/xxx
        doi_match = re.search(
            r"/doi/(?:full/|abs/|epdf/|pdf/)?(\d+\.\d+/[^\s?#]+)", url, re.IGNORECASE
        )
        if doi_match:
            # Decode URL-encoded DOI (e.g., %3C -> <, %3E -> >)
            # This handles SICI-style DOIs that may be pre-encoded in URLs
            raw_doi = unquote(doi_match.group(1))
            return raw_doi
        return None

    def _rate_limit(self):
        """Enforce rate limiting between API requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()

    @staticmethod
    def encode_doi_for_api(doi: str) -> str:
        """
        Encode a DOI for use in Wiley TDM API URLs.

        Based on the Wiley TDM API curl example, DOIs should have minimal encoding.
        The API expects the DOI path segment with only truly unsafe URL characters
        encoded. For SICI-style DOIs, only < and > need encoding.

        Example: DOI 10.1002/(SICI)1098-1004(200005)15:5<483::AID-HUMU18>3.0.CO;2-T
        becomes: 10.1002/(SICI)1098-1004(200005)15:5%3C483::AID-HUMU18%3E3.0.CO;2-T

        Args:
            doi: Digital Object Identifier

        Returns:
            URL-encoded DOI string suitable for Wiley TDM API URLs
        """
        # Match the curl example: minimal encoding - only < and > need escaping
        # Parentheses, colons, semicolons, and slashes stay unencoded
        return doi.replace("<", "%3C").replace(">", "%3E")

    def get_fulltext_by_doi(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Fetch full-text content from Wiley TDM API using DOI.

        Uses the same URL format as the official curl example:
        curl -L -H "Wiley-TDM-Client-Token: TOKEN" \\
             https://api.wiley.com/onlinelibrary/tdm/v1/articles/DOI

        The API may redirect (hence -L in curl) and can return PDF or XML/HTML.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (content, error_message)
            - content: Full-text XML/HTML/PDF-converted if successful, None otherwise
            - error_message: Error description if failed, None otherwise
        """
        if not self.is_available:
            return None, "Wiley API key not configured"

        # TDM API headers - match curl example closely
        # Only the token header is required; redirects are followed automatically
        headers = {
            "Wiley-TDM-Client-Token": self.api_key,
            # Accept both text formats and PDF
            "Accept": "application/pdf, application/xml, text/xml, application/xhtml+xml, text/html, */*",
        }

        def _request_fulltext(target_url: str) -> requests.Response:
            self._rate_limit()
            logger.info(f"Fetching full text from Wiley API for DOI: {doi}")
            logger.debug(f"TDM API URL: {target_url}")
            # Follow redirects (curl -L behavior) - this is critical per Wiley docs
            return self.session.get(
                target_url, headers=headers, timeout=HTTP_TIMEOUT_DEFAULT, allow_redirects=True
            )

        def _handle_response(
            response: requests.Response,
        ) -> Tuple[Optional[bytes], Optional[str], str]:
            """Returns (content_bytes, error, content_type)"""
            if response.status_code == 200:
                content_type = response.headers.get("Content-Type", "").lower()
                return response.content, None, content_type
            if response.status_code == 401:
                return None, "Invalid or unauthorized API key", ""
            if response.status_code == 403:
                return (
                    None,
                    "Access forbidden - API key may lack permissions or article not available",
                    "",
                )
            if response.status_code == 404:
                return None, "Article not found via Wiley API", ""
            if response.status_code == 429:
                return None, "Rate limit exceeded", ""
            return None, f"HTTP {response.status_code}: {response.reason}", ""

        def _process_content(
            content_bytes: bytes, content_type: str
        ) -> Tuple[Optional[str], Optional[str]]:
            """Process response content - handle both PDF and text formats."""
            # Check if response is PDF (by content-type or magic bytes)
            is_pdf = "pdf" in content_type or content_bytes.startswith(b"%PDF")

            if is_pdf:
                # Convert PDF to markdown using FormatConverter
                try:
                    converter = FormatConverter()
                    with tempfile.NamedTemporaryFile(
                        suffix=".pdf", delete=True
                    ) as tmp_file:
                        tmp_file.write(content_bytes)
                        tmp_file.flush()
                        markdown = converter.pdf_to_markdown(Path(tmp_file.name))
                    if markdown and len(markdown) > 500:
                        logger.info(
                            f"Successfully converted TDM API PDF response for DOI: {doi}"
                        )
                        return markdown, None
                    return None, "PDF conversion produced insufficient content"
                except Exception as e:
                    logger.debug(f"Error converting TDM API PDF: {e}")
                    return None, f"PDF conversion failed: {str(e)}"
            else:
                # Text content (XML/HTML)
                try:
                    text_content = content_bytes.decode("utf-8")
                    return text_content, None
                except UnicodeDecodeError:
                    try:
                        text_content = content_bytes.decode("latin-1")
                        return text_content, None
                    except Exception:
                        return None, "Failed to decode response content"

        # Try multiple URL encoding strategies for SICI-style DOIs
        # Strategy 1: Minimal encoding (matches curl example) - only < and > encoded
        encoding_strategies = [
            ("minimal", self.encode_doi_for_api(doi)),
            # Strategy 2: Keep slash and common chars unencoded
            ("standard", quote(doi, safe="/:();-")),
            # Strategy 3: Full encoding as fallback
            ("full", quote(doi, safe="")),
        ]

        # Remove duplicate encodings
        seen = set()
        unique_strategies = []
        for name, encoded in encoding_strategies:
            if encoded not in seen:
                seen.add(encoded)
                unique_strategies.append((name, encoded))

        last_error = None
        try:
            for strategy_name, encoded_doi in unique_strategies:
                url = f"{self.BASE_URL}/{encoded_doi}"
                logger.debug(f"Trying TDM API with {strategy_name} encoding: {url}")

                response = _request_fulltext(url)
                content_bytes, error, content_type = _handle_response(response)

                if content_bytes:
                    text_content, process_error = _process_content(
                        content_bytes, content_type
                    )
                    if text_content:
                        return text_content, None
                    last_error = process_error
                elif response.status_code == 404:
                    # Try next encoding strategy
                    last_error = error
                    continue
                else:
                    # Non-404 error, don't try other encodings
                    return None, error

            return None, last_error or "Article not found via Wiley API"

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
            xml_content = re.sub(r'xmlns[^=]*="[^"]*"', "", xml_content)
            xml_content = re.sub(r"<([a-z]+):", r"<", xml_content)
            xml_content = re.sub(r"</([a-z]+):", r"</", xml_content)

            root = ET.fromstring(xml_content)
            markdown = "# MAIN TEXT\n\n"

            # Extract title
            title_elem = root.find(".//article-title")
            if title_elem is None:
                title_elem = root.find(".//title")
            if title_elem is None:
                title_elem = root.find(".//{*}title")
            if title_elem is not None:
                title_text = self._extract_text(title_elem)
                if title_text:
                    markdown += f"## {title_text.strip()}\n\n"

            # Extract abstract
            abstract_elem = root.find(".//abstract")
            if abstract_elem is None:
                abstract_elem = root.find(".//{*}abstract")
            if abstract_elem is not None:
                abstract_text = self._extract_text(abstract_elem)
                if abstract_text:
                    markdown += f"### Abstract\n\n{abstract_text}\n\n"

            # Extract body sections
            body_elem = root.find(".//body")
            if body_elem is None:
                body_elem = root.find(".//{*}body")

            if body_elem is not None:
                # Process sections
                sections = body_elem.findall(".//sec")
                if not sections:
                    sections = body_elem.findall(".//{*}sec")
                if not sections:
                    sections = body_elem.findall(".//section")
                if not sections:
                    sections = body_elem.findall(".//{*}section")

                if sections:
                    for section in sections:
                        section_md = self._process_section(section)
                        if section_md:
                            markdown += section_md
                else:
                    # No sections found, extract all paragraphs from body
                    paragraphs = body_elem.findall(".//p")
                    if not paragraphs:
                        paragraphs = body_elem.findall(".//{*}p")
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
            soup = BeautifulSoup(html_content, "html.parser")
            markdown = "# MAIN TEXT\n\n"

            # Extract title
            title = soup.find("h1", class_="citation__title")
            if not title:
                title = soup.find("h1")
            if title:
                markdown += f"## {title.get_text(strip=True)}\n\n"

            # Extract abstract
            abstract = soup.find("section", class_="article-section__abstract")
            if not abstract:
                abstract = soup.find("div", class_="abstract")
            if abstract:
                abstract_text = abstract.get_text(separator=" ", strip=True)
                markdown += f"### Abstract\n\n{abstract_text}\n\n"

            # Extract main content sections
            content_sections = soup.find_all(
                "section", class_="article-section__content"
            )
            if not content_sections:
                content_sections = soup.find_all("div", class_="article__body")

            for section in content_sections:
                # Get section heading
                heading = section.find(["h2", "h3", "h4"])
                if heading:
                    markdown += f"### {heading.get_text(strip=True)}\n\n"

                # Get paragraphs
                paragraphs = section.find_all("p")
                for para in paragraphs:
                    para_text = para.get_text(separator=" ", strip=True)
                    if para_text and len(para_text) > 20:
                        markdown += f"{para_text}\n\n"

            # Fallback: extract all paragraphs if sections not found
            if len(markdown) < 500:
                all_paragraphs = soup.find_all("p")
                for para in all_paragraphs:
                    para_text = para.get_text(separator=" ", strip=True)
                    if para_text and len(para_text) > 50:
                        markdown += f"{para_text}\n\n"

            return markdown if len(markdown) > 200 else None

        except Exception as e:
            logger.debug(f"Error parsing Wiley HTML: {e}")
            return None

    def scrape_fulltext_from_web(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Scrape full-text content directly from Wiley Online Library website.

        This is a fallback method for when the TDM API fails (e.g., for older
        SICI-style DOIs that aren't indexed in the TDM system).

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (markdown_content, error_message)
            - markdown_content: Markdown text if successful, None otherwise
            - error_message: Error description if failed, None otherwise
        """
        # URL-encode the DOI for the web URL (same encoding as browser)
        encoded_doi = self.encode_doi_for_web(doi)
        # Wiley HTML full-text URL format
        full_text_url = f"https://onlinelibrary.wiley.com/doi/full/{encoded_doi}"

        logger.info(f"Attempting to scrape Wiley web page for DOI: {doi}")
        logger.debug(f"Wiley web URL: {full_text_url}")

        last_error = None

        try:
            self._rate_limit()
            response = self.session.get(
                full_text_url,
                timeout=HTTP_TIMEOUT_DEFAULT,
                allow_redirects=True,
                headers=self._get_wiley_headers(
                    referer="https://onlinelibrary.wiley.com/"
                ),
            )

            if response.status_code == 403:
                pdf_markdown, pdf_error = self._fetch_pdf_fulltext(doi)
                if pdf_markdown:
                    return pdf_markdown, None
                last_error = pdf_error or "Access forbidden - article may be paywalled"
            elif response.status_code == 404:
                last_error = "Article not found at /doi/full/ URL"
            elif response.status_code != 200:
                last_error = f"HTTP {response.status_code}: {response.reason}"
            else:
                html_content = response.text

                # Check for paywall indicators
                paywall_indicators = [
                    "purchase this article",
                    "get access",
                    "institutional access",
                    "sign in to access",
                    "subscription required",
                ]
                html_lower = html_content.lower()
                if any(indicator in html_lower for indicator in paywall_indicators):
                    # Check if we also have full article content (some pages show both)
                    if "article-section__content" not in html_lower:
                        last_error = "Article is behind paywall"

                if not last_error:
                    # Extract content using BeautifulSoup
                    markdown = self._scrape_wiley_html(html_content)
                    if markdown and len(markdown) > 2000:
                        logger.info(
                            f"Successfully scraped Wiley article via web for DOI: {doi}"
                        )
                        return markdown, None
                    elif markdown and len(markdown) > 500:
                        pdf_markdown, pdf_error = self._fetch_pdf_fulltext(doi)
                        if pdf_markdown:
                            return pdf_markdown, None
                        last_error = (
                            pdf_error
                            or "Web scrape produced insufficient content (abstract only)"
                        )
                    else:
                        pdf_markdown, pdf_error = self._fetch_pdf_fulltext(doi)
                        if pdf_markdown:
                            return pdf_markdown, None
                        last_error = (
                            pdf_error or "Could not extract content from Wiley web page"
                        )

        except requests.exceptions.Timeout:
            last_error = "Web request timed out"
        except requests.exceptions.RequestException as e:
            last_error = f"Web request failed: {str(e)}"

        # If /doi/full/ failed, try the landing page via DOI resolution
        # This works better for old SICI-style DOIs
        if last_error:
            logger.info(f"Trying DOI resolution fallback for: {doi}")
            landing_markdown, landing_error = self._try_doi_landing_page(doi)
            if landing_markdown:
                return landing_markdown, None
            # Combine errors for debugging
            if landing_error:
                return None, f"{last_error}; Landing page: {landing_error}"

        return None, last_error

    def _try_doi_landing_page(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Try to fetch full text by resolving DOI and scraping the landing page.

        For old SICI-style DOIs, the /doi/full/ URL may not work, but resolving
        through doi.org often leads to a working page with full text links.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (markdown_content, error_message)
        """
        # Encode DOI for doi.org URL - keep only slash unencoded to handle
        # SICI-style DOIs with special characters like parentheses, colons, etc.
        encoded_doi = quote(doi, safe="/")
        doi_url = f"https://doi.org/{encoded_doi}"

        try:
            self._rate_limit()
            response = self.session.get(doi_url, timeout=HTTP_TIMEOUT_DEFAULT, allow_redirects=True)

            if (
                response.status_code == 403
                and response.url
                and "wiley.com" in response.url.lower()
            ):
                response = self.session.get(
                    response.url,
                    timeout=HTTP_TIMEOUT_DEFAULT,
                    allow_redirects=True,
                    headers=self._get_wiley_headers(referer="https://doi.org/"),
                )

            if response.status_code != 200:
                return None, f"DOI resolution failed: HTTP {response.status_code}"

            final_url = response.url
            html_content = response.text

            # Check if we landed on a Wiley page
            if "wiley.com" not in final_url.lower():
                return None, f"DOI resolved to non-Wiley site: {final_url}"

            logger.info(f"DOI resolved to: {final_url}")

            # First, try to find and follow a "Full Text" link on the landing page
            soup = BeautifulSoup(html_content, "html.parser")

            # Look for full text links (various patterns used by Wiley)
            fulltext_link = None
            fulltext_patterns = [
                ("a", {"class": re.compile(r"full-?text", re.I)}),
                ("a", {"href": re.compile(r"/doi/full/", re.I)}),
                ("a", {"title": re.compile(r"full\s*text", re.I)}),
                ("a", {"data-test": "full-text-link"}),
            ]

            for tag, attrs in fulltext_patterns:
                link = soup.find(tag, attrs)
                if link and link.get("href"):
                    fulltext_link = link["href"]
                    break

            # Also check for common link text patterns
            if not fulltext_link:
                for a in soup.find_all("a", href=True):
                    link_text = a.get_text(strip=True).lower()
                    if any(
                        pattern in link_text
                        for pattern in ["full text", "full-text", "read full"]
                    ):
                        fulltext_link = a["href"]
                        break

            # If found a full text link, follow it
            if fulltext_link:
                if not fulltext_link.startswith("http"):
                    # Relative URL - construct absolute
                    from urllib.parse import urljoin

                    fulltext_link = urljoin(final_url, fulltext_link)

                logger.info(f"Following full text link: {fulltext_link}")
                self._rate_limit()
                ft_response = self.session.get(
                    fulltext_link,
                    timeout=HTTP_TIMEOUT_DEFAULT,
                    allow_redirects=True,
                    headers=self._get_wiley_headers(referer=final_url),
                )

                if ft_response.status_code == 200:
                    markdown = self._scrape_wiley_html(ft_response.text)
                    if markdown and len(markdown) > 2000:
                        logger.info(
                            "Successfully scraped Wiley article via full text link"
                        )
                        return markdown, None

            # Try to extract from the landing page itself (some old articles have inline text)
            markdown = self._scrape_wiley_html(html_content)
            if markdown and len(markdown) > 2000:
                logger.info("Successfully scraped Wiley article from landing page")
                return markdown, None

            # Try PDF as final fallback
            pdf_markdown, pdf_error = self._fetch_pdf_fulltext(doi)
            if pdf_markdown:
                return pdf_markdown, None

            return None, pdf_error or "Could not extract full text from landing page"

        except requests.exceptions.Timeout:
            return None, "DOI resolution timed out"
        except requests.exceptions.RequestException as e:
            return None, f"DOI resolution failed: {str(e)}"

    def _scrape_wiley_html(self, html_content: str) -> Optional[str]:
        """
        Extract article content from Wiley Online Library HTML page.

        Args:
            html_content: Raw HTML from Wiley website

        Returns:
            Markdown string if extraction successful, None otherwise
        """
        try:
            soup = BeautifulSoup(html_content, "html.parser")
            markdown = "# MAIN TEXT\n\n"

            # Extract title - try multiple patterns for old and new layouts
            title_elem = soup.find("h1", class_="citation__title")
            if not title_elem:
                title_elem = soup.find("h1", class_="article-header__title")
            if not title_elem:
                title_elem = soup.find("h1", id="articleTitle")
            if not title_elem:
                title_elem = soup.find("h1")
            if title_elem:
                title = title_elem.get_text(strip=True)
                markdown += f"## {title}\n\n"

            # Extract abstract - try multiple patterns
            abstract = soup.find("section", class_="article-section__abstract")
            if not abstract:
                abstract = soup.find("div", class_="abstract")
            if not abstract:
                abstract = soup.find("div", id="abstract")
            if not abstract:
                # Old layout: abstract might be in a paragraph with specific id
                abstract = soup.find("p", id="abstract")
            if abstract:
                abstract_text = abstract.get_text(separator=" ", strip=True)
                markdown += f"### Abstract\n\n{abstract_text}\n\n"

            # Extract main content sections
            # Modern Wiley uses article-section__full divs for each section
            content_sections = soup.find_all("section", class_="article-section__full")
            if not content_sections:
                content_sections = soup.find_all(
                    "section", class_="article-section__content"
                )
            if not content_sections:
                content_sections = soup.find_all("div", class_="article__body")

            # Old Wiley layout: try different container patterns
            if not content_sections:
                content_sections = soup.find_all(
                    "div", class_=re.compile(r"article-?body", re.I)
                )
            if not content_sections:
                content_sections = soup.find_all(
                    "div",
                    id=re.compile(r"(fulltext|article-?content|main-?content)", re.I),
                )
            if not content_sections:
                # Very old layouts: look for main article div
                main_article = soup.find("div", class_="mainContent")
                if main_article:
                    content_sections = [main_article]

            for section in content_sections:
                # Get section heading
                heading = section.find(
                    ["h2", "h3", "h4"], class_="article-section__title"
                )
                if not heading:
                    heading = section.find(["h2", "h3", "h4"])
                if heading:
                    heading_text = heading.get_text(strip=True)
                    # Skip duplicate abstract heading
                    if heading_text.lower() != "abstract":
                        markdown += f"### {heading_text}\n\n"

                # Get paragraphs within this section
                paragraphs = section.find_all("p")
                for para in paragraphs:
                    para_text = para.get_text(separator=" ", strip=True)
                    if para_text and len(para_text) > 20:
                        markdown += f"{para_text}\n\n"

            # Fallback 1: if structured extraction didn't find enough content,
            # try extracting all paragraphs from the main content area
            if len(markdown) < 1000:
                main_content = soup.find("div", class_="article__content")
                if not main_content:
                    main_content = soup.find("article")
                if not main_content:
                    main_content = soup.find("div", id="articleBody")
                if main_content:
                    for para in main_content.find_all("p"):
                        para_text = para.get_text(separator=" ", strip=True)
                        if para_text and len(para_text) > 50:
                            markdown += f"{para_text}\n\n"

            # Fallback 2: Try to find content in old-style Wiley pages with
            # nested divs containing article sections
            if len(markdown) < 1000:
                # Look for sections by header text patterns
                for header in soup.find_all(["h2", "h3", "h4"]):
                    header_text = header.get_text(strip=True)
                    header_lower = header_text.lower()
                    if any(
                        kw in header_lower
                        for kw in [
                            "introduction",
                            "methods",
                            "results",
                            "discussion",
                            "materials",
                            "patients",
                            "background",
                            "conclusion",
                        ]
                    ):
                        markdown += f"### {header_text}\n\n"
                        # Get following paragraphs until next header
                        next_elem = header.find_next_sibling()
                        while next_elem and next_elem.name not in ["h2", "h3", "h4"]:
                            if next_elem.name == "p":
                                para_text = next_elem.get_text(
                                    separator=" ", strip=True
                                )
                                if para_text and len(para_text) > 20:
                                    markdown += f"{para_text}\n\n"
                            next_elem = next_elem.find_next_sibling()

            # Ultimate fallback: extract all substantial paragraphs
            if len(markdown) < 500:
                all_paragraphs = soup.find_all("p")
                for para in all_paragraphs:
                    para_text = para.get_text(separator=" ", strip=True)
                    if para_text and len(para_text) > 100:
                        markdown += f"{para_text}\n\n"

            return markdown if len(markdown) > 200 else None

        except Exception as e:
            logger.debug(f"Error scraping Wiley HTML: {e}")
            return None

    def _fetch_pdf_fulltext(self, doi: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Attempt to download a Wiley PDF/EPDF and convert it to markdown.

        Args:
            doi: Digital Object Identifier

        Returns:
            Tuple of (markdown_content, error_message)
        """
        encoded_doi = self.encode_doi_for_web(doi)
        pdf_urls = [
            f"https://onlinelibrary.wiley.com/doi/epdf/{encoded_doi}",
            f"https://onlinelibrary.wiley.com/doi/pdf/{encoded_doi}",
        ]
        converter = FormatConverter()

        for pdf_url in pdf_urls:
            try:
                self._rate_limit()
                response = self.session.get(
                    pdf_url,
                    timeout=HTTP_TIMEOUT_DEFAULT,
                    allow_redirects=True,
                    headers=self._get_wiley_headers(
                        referer="https://onlinelibrary.wiley.com/"
                    ),
                )
                if response.status_code != 200:
                    continue
                content_type = response.headers.get("Content-Type", "").lower()
                if "pdf" not in content_type and not response.content.startswith(
                    b"%PDF"
                ):
                    continue
                with tempfile.NamedTemporaryFile(
                    suffix=".pdf", delete=True
                ) as tmp_file:
                    tmp_file.write(response.content)
                    tmp_file.flush()
                    markdown = converter.pdf_to_markdown(Path(tmp_file.name))
                if markdown and len(markdown) > 500:
                    logger.info(
                        f"Successfully converted Wiley PDF to markdown for DOI: {doi}"
                    )
                    return markdown, None
            except requests.exceptions.Timeout:
                continue
            except requests.exceptions.RequestException:
                continue
            except Exception as e:
                logger.debug(f"Error converting Wiley PDF to markdown: {e}")
                continue

        return None, "Wiley PDF full text not available"

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
        title_elem = section_elem.find("./title")
        if title_elem is None:
            title_elem = section_elem.find("./{*}title")

        if title_elem is not None:
            title_text = self._extract_text(title_elem)
            if title_text:
                heading_prefix = "#" * min(level, 6)
                markdown += f"{heading_prefix} {title_text.strip()}\n\n"

        # Extract paragraphs
        paragraphs = section_elem.findall("./p")
        if not paragraphs:
            paragraphs = section_elem.findall("./{*}p")

        for para in paragraphs:
            para_text = self._extract_text(para)
            if para_text:
                markdown += f"{para_text}\n\n"

        # Process nested sections
        nested_sections = section_elem.findall("./sec")
        if not nested_sections:
            nested_sections = section_elem.findall("./{*}sec")
        if not nested_sections:
            nested_sections = section_elem.findall("./section")

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

    def fetch_fulltext(
        self,
        doi: Optional[str] = None,
        url: Optional[str] = None,
        try_web_scraping: bool = True,
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Fetch full-text content and convert to markdown.

        Tries to fetch using DOI via the TDM API first. If the API fails
        (e.g., for older SICI-style DOIs), falls back to web scraping.
        If a URL is provided and no DOI, extracts DOI from the URL.

        Args:
            doi: Digital Object Identifier
            url: URL (used to extract DOI if not provided)
            try_web_scraping: If True, fall back to web scraping when API fails

        Returns:
            Tuple of (markdown_content, error_message)
            - markdown_content: Markdown text if successful, None otherwise
            - error_message: Error description if failed, None otherwise
        """
        # Extract DOI from URL if not provided
        if not doi and url:
            doi = self.extract_doi_from_url(url)

        if not doi:
            return None, "No valid DOI available"

        api_error = None

        # Try TDM API first (if API key is available)
        if self.is_available:
            content, api_error = self.get_fulltext_by_doi(doi)
            if content:
                markdown = self.content_to_markdown(content)
                # Check for sufficient content - abstracts are typically 200-400 words
                # Full articles should be at least 2000 characters (roughly 350+ words)
                # Also check for section headings beyond just Abstract
                MIN_FULLTEXT_LENGTH = 2000
                has_body_sections = markdown and any(
                    section in markdown.lower()
                    for section in [
                        "### introduction",
                        "### methods",
                        "### results",
                        "### discussion",
                        "### materials",
                        "### content",
                    ]
                )
                if markdown and (
                    len(markdown) > MIN_FULLTEXT_LENGTH or has_body_sections
                ):
                    logger.info(
                        f"Successfully fetched Wiley article via TDM API: {doi}"
                    )
                    return markdown, None
                elif markdown and len(markdown) > 500:
                    # Got some content but likely just abstract
                    api_error = "TDM API returned insufficient content (abstract only)"
                elif not api_error:
                    api_error = (
                        "TDM API content conversion produced insufficient content"
                    )

            # Log API failure for debugging
            if api_error:
                logger.info(f"Wiley TDM API failed for DOI {doi}: {api_error}")

        # Fall back to web scraping if API failed or isn't available
        if try_web_scraping:
            logger.info(f"Falling back to Wiley web scraping for DOI: {doi}")
            web_markdown, web_error = self.scrape_fulltext_from_web(doi)
            if web_markdown:
                logger.info(
                    f"Successfully fetched Wiley article via web scraping: {doi}"
                )
                return web_markdown, None

            # Combine errors for debugging
            if web_error:
                if api_error:
                    return None, f"TDM API: {api_error}; Web scraping: {web_error}"
                return None, f"Web scraping: {web_error}"

        return None, api_error or "Failed to fetch content"
