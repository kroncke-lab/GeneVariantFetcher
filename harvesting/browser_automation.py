"""
Browser Automation Module - Playwright-based download for stubborn papers.

Uses headless browser to bypass JavaScript-heavy pages and some anti-bot measures.
This is the fallback when API-based methods fail.
"""

import asyncio
import logging
import os
import re
from pathlib import Path
from typing import Optional, Tuple, List, TYPE_CHECKING, Any
from dataclasses import dataclass

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from playwright.async_api import Page, Browser

# Check if playwright is available
try:
    from playwright.async_api import async_playwright, Page, Browser

    PLAYWRIGHT_AVAILABLE = True
except ImportError:
    PLAYWRIGHT_AVAILABLE = False
    Page = None  # type: ignore
    Browser = None  # type: ignore
    logger.warning("Playwright not installed. Browser automation unavailable.")


@dataclass
class BrowserDownloadResult:
    """Result of a browser-based download attempt."""

    success: bool
    content: Optional[str] = None  # Markdown content
    pdf_path: Optional[str] = None  # Path to downloaded PDF
    supplements: List[str] = None  # List of supplement file paths
    error: Optional[str] = None
    source_url: Optional[str] = None


class BrowserDownloader:
    """
    Playwright-based downloader for papers that resist API access.

    Usage:
        downloader = BrowserDownloader(output_dir=Path("downloads"))
        result = await downloader.download_paper(doi="10.1234/example")
    """

    # Publisher URL patterns and their handlers
    PUBLISHER_PATTERNS = {
        "elsevier": [
            r"sciencedirect\.com",
            r"cell\.com",
            r"thelancet\.com",
        ],
        "wiley": [
            r"onlinelibrary\.wiley\.com",
            r"ahajournals\.org",
        ],
        "springer": [
            r"link\.springer\.com",
            r"nature\.com",
            r"biomedcentral\.com",
        ],
        "oxford": [
            r"academic\.oup\.com",
        ],
        "karger": [
            r"karger\.com",
        ],
    }

    # Supplement link patterns
    SUPPLEMENT_PATTERNS = [
        r"supplement",
        r"supporting\s*info",
        r"additional\s*file",
        r"extended\s*data",
        r"data\s*availability",
        r"online\s*resource",
        r"appendix",
        r"ESM",
        r"MOESM",
    ]

    def __init__(
        self,
        output_dir: Path,
        headless: bool = True,
        timeout: int = 30000,  # milliseconds
    ):
        """
        Initialize browser downloader.

        Args:
            output_dir: Directory to save downloads
            headless: Run browser in headless mode
            timeout: Page load timeout in milliseconds
        """
        if not PLAYWRIGHT_AVAILABLE:
            raise RuntimeError(
                "Playwright is not installed. Run: pip install playwright && playwright install chromium"
            )

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.headless = headless
        self.timeout = timeout
        self._browser: Optional["Browser"] = None

    async def __aenter__(self):
        """Async context manager entry."""
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        if self._browser:
            await self._browser.close()

    def _detect_publisher(self, url: str) -> Optional[str]:
        """Detect publisher from URL."""
        for publisher, patterns in self.PUBLISHER_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, url, re.IGNORECASE):
                    return publisher
        return None

    async def _setup_browser(self) -> "Browser":
        """Initialize browser with stealth settings."""
        if self._browser:
            return self._browser

        playwright = await async_playwright().start()

        self._browser = await playwright.chromium.launch(
            headless=self.headless,
            args=[
                "--no-sandbox",
                "--disable-blink-features=AutomationControlled",
            ],
        )

        return self._browser

    async def _create_context(self):
        """Create a new browser context with realistic settings."""
        browser = await self._setup_browser()

        context = await browser.new_context(
            viewport={"width": 1280, "height": 900},
            user_agent="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            java_script_enabled=True,
        )

        return context

    async def download_paper(
        self,
        url: Optional[str] = None,
        doi: Optional[str] = None,
        pmid: Optional[str] = None,
    ) -> BrowserDownloadResult:
        """
        Download paper using browser automation.

        Args:
            url: Direct URL to paper
            doi: DOI (will resolve to URL)
            pmid: PubMed ID (will get URL from PubMed)

        Returns:
            BrowserDownloadResult with content/paths
        """
        supplements = []

        # Resolve URL if not provided
        if not url:
            if doi:
                url = f"https://doi.org/{doi}"
            elif pmid:
                url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            else:
                return BrowserDownloadResult(
                    success=False, error="No URL, DOI, or PMID provided"
                )

        try:
            context = await self._create_context()
            page = await context.new_page()

            logger.info(f"Browser navigating to: {url}")

            # Navigate to page
            response = await page.goto(
                url, wait_until="domcontentloaded", timeout=self.timeout
            )

            if not response or response.status >= 400:
                return BrowserDownloadResult(
                    success=False,
                    error=f"HTTP {response.status if response else 'No response'}",
                    source_url=url,
                )

            # Wait for content to load
            await page.wait_for_timeout(3000)

            # Get final URL (after redirects)
            final_url = page.url
            publisher = self._detect_publisher(final_url)

            logger.info(
                f"Resolved to: {final_url} (publisher: {publisher or 'unknown'})"
            )

            # Try to get full text content
            content = await self._extract_content(page, publisher)

            # Try to find and download supplements
            supplements = await self._find_supplements(page, pmid or doi or "paper")

            await context.close()

            if content and len(content) > 500:
                return BrowserDownloadResult(
                    success=True,
                    content=content,
                    supplements=supplements,
                    source_url=final_url,
                )
            else:
                return BrowserDownloadResult(
                    success=False,
                    error="Content too short or not found",
                    supplements=supplements,
                    source_url=final_url,
                )

        except Exception as e:
            logger.error(f"Browser download failed: {e}")
            return BrowserDownloadResult(success=False, error=str(e), source_url=url)

    async def _extract_content(
        self, page: "Page", publisher: Optional[str]
    ) -> Optional[str]:
        """Extract article content from page."""

        # Publisher-specific selectors
        selectors = {
            "elsevier": [
                "#body",
                ".article-body",
                ".Body",
            ],
            "wiley": [
                ".article-body-main",
                "#article__content",
                ".article__body",
            ],
            "springer": [
                "#main-content",
                ".c-article-body",
                "article",
            ],
            "oxford": [
                ".article-body",
                "#ContentColumn",
            ],
            "karger": [
                ".article-body",
                "#articleBody",
            ],
        }

        # Get publisher-specific selectors first, then generic
        try_selectors = selectors.get(publisher, []) + [
            "article",
            ".article",
            ".article-content",
            "#content",
            "main",
            ".main-content",
        ]

        for selector in try_selectors:
            try:
                element = await page.query_selector(selector)
                if element:
                    text = await element.inner_text()
                    if text and len(text) > 500:
                        # Basic markdown formatting
                        content = f"# Article Content\n\nSource: {page.url}\n\n{text}"
                        return content
            except Exception:
                continue

        # Fallback: get all body text
        try:
            body_text = await page.inner_text("body")
            if body_text and len(body_text) > 1000:
                return (
                    f"# Article Content (fallback)\n\nSource: {page.url}\n\n{body_text}"
                )
        except Exception:
            pass

        return None

    async def _find_supplements(self, page: Page, identifier: str) -> List[str]:
        """Find and download supplement files."""
        supplements = []

        # Look for supplement links
        links = await page.query_selector_all("a[href]")

        supplement_urls = []
        for link in links:
            try:
                href = await link.get_attribute("href")
                text = await link.inner_text()

                if not href:
                    continue

                # Check if this looks like a supplement
                combined = f"{href} {text}".lower()
                for pattern in self.SUPPLEMENT_PATTERNS:
                    if re.search(pattern, combined, re.IGNORECASE):
                        # Check for downloadable file extensions
                        if re.search(
                            r"\.(pdf|docx?|xlsx?|csv|zip|rar)(\?|$)",
                            href,
                            re.IGNORECASE,
                        ):
                            supplement_urls.append(href)
                            break
            except Exception:
                continue

        # Download each supplement
        for i, url in enumerate(supplement_urls[:10]):  # Limit to 10 supplements
            try:
                # Make absolute URL if needed
                if url.startswith("/"):
                    base = page.url.split("/")[0:3]
                    url = "/".join(base) + url
                elif not url.startswith("http"):
                    url = page.url.rsplit("/", 1)[0] + "/" + url

                # Download file
                ext = re.search(r"\.(\w+)(\?|$)", url)
                ext = ext.group(1) if ext else "pdf"

                output_path = self.output_dir / f"{identifier}_supp_{i + 1}.{ext}"

                # Use page context to download (maintains cookies)
                response = await page.request.get(url)
                if response.ok:
                    content = await response.body()
                    with open(output_path, "wb") as f:
                        f.write(content)
                    supplements.append(str(output_path))
                    logger.info(f"Downloaded supplement: {output_path.name}")
            except Exception as e:
                logger.warning(f"Failed to download supplement {url}: {e}")

        return supplements


async def download_with_browser(
    url: Optional[str] = None,
    doi: Optional[str] = None,
    pmid: Optional[str] = None,
    output_dir: Path = Path("downloads"),
    headless: bool = True,
) -> BrowserDownloadResult:
    """
    Convenience function for one-off browser downloads.

    Usage:
        result = await download_with_browser(doi="10.1234/example")
    """
    async with BrowserDownloader(
        output_dir=output_dir, headless=headless
    ) as downloader:
        return await downloader.download_paper(url=url, doi=doi, pmid=pmid)


def download_paper_sync(
    url: Optional[str] = None,
    doi: Optional[str] = None,
    pmid: Optional[str] = None,
    output_dir: Path = Path("downloads"),
    headless: bool = True,
) -> BrowserDownloadResult:
    """
    Synchronous wrapper for browser download.

    Usage:
        result = download_paper_sync(doi="10.1234/example")
    """
    return asyncio.run(download_with_browser(url, doi, pmid, output_dir, headless))
