#!/usr/bin/env python3
"""
Browser-Based Paper Fetcher with LLM Assistance

Automates downloading paywalled papers using browser automation (Playwright)
with optional LLM assistance to find download links on complex publisher sites.

This avoids unreliable publisher APIs by using actual browser navigation,
similar to how a human would manually download papers.

Requirements:
    pip install playwright anthropic
    playwright install chromium

Usage:
    # Basic usage - opens browser and attempts to find/download PDFs
    python browser_fetch.py paywalled_missing.csv

    # With Claude assistance for finding download links
    python browser_fetch.py paywalled_missing.csv --use-claude

    # Headless mode (faster, no visible browser)
    python browser_fetch.py paywalled_missing.csv --headless

    # Process specific PMIDs
    python browser_fetch.py paywalled_missing.csv --pmids 12345678,87654321
"""

import json
import logging
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import requests

# Import supplement scraper from harvesting module
try:
    from harvesting.supplement_scraper import SupplementScraper

    SCRAPER_AVAILABLE = True
except ImportError:
    SCRAPER_AVAILABLE = False
    SupplementScraper = None


# Setup logging with immediate flush for real-time output
class FlushingStreamHandler(logging.StreamHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
    handlers=[FlushingStreamHandler()],
)
logger = logging.getLogger(__name__)

# Try to import Playwright
try:
    from playwright.sync_api import (
        Browser,
        Page,
        sync_playwright,
    )
    from playwright.sync_api import (
        TimeoutError as PlaywrightTimeout,
    )

    PLAYWRIGHT_AVAILABLE = True
except ImportError:
    PLAYWRIGHT_AVAILABLE = False
    # Define dummy types for type hints when Playwright not installed
    Page = Any
    Browser = Any
    PlaywrightTimeout = Exception
    sync_playwright = None
    logger.warning(
        "Playwright not installed. Run: pip install playwright && playwright install chromium"
    )

# Try to import Anthropic for Claude assistance
try:
    import anthropic

    ANTHROPIC_AVAILABLE = True
except ImportError:
    ANTHROPIC_AVAILABLE = False


@dataclass
class DownloadResult:
    """Result of a download attempt."""

    pmid: str
    success: bool
    files_downloaded: List[str]
    error: Optional[str] = None
    method: str = "unknown"  # "direct", "claude_assisted", "manual"
    failure_type: Optional[str] = (
        None  # "captcha_blocked", "paywall", "timeout", "crash", None for success
    )


class PageAnalyzer:
    """Analyzes web pages to find PDF download links using Claude."""

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        self.client = None
        if self.api_key and ANTHROPIC_AVAILABLE:
            self.client = anthropic.Anthropic(api_key=self.api_key)

    def find_pdf_links(self, page_content: str, page_url: str) -> List[Dict[str, str]]:
        """
        Use Claude to analyze page content and find PDF download links.

        Args:
            page_content: HTML or text content of the page
            page_url: URL of the page being analyzed

        Returns:
            List of dicts with 'url', 'type' (main/supplement), and 'confidence'
        """
        if not self.client:
            return []

        # Truncate content if too long
        max_content = 50000
        if len(page_content) > max_content:
            page_content = page_content[:max_content] + "\n... [truncated]"

        prompt = f"""Analyze this webpage to find PDF download links for a scientific paper.

Page URL: {page_url}

Page content:
{page_content}

Find all PDF download links. Look for:
1. Main article PDF (full text)
2. Supplementary materials (tables, figures, data files)
3. Both direct download links and links that lead to download pages

Return a JSON array of objects with:
- "url": the download URL (make absolute if relative)
- "type": "main" or "supplement"
- "description": brief description of what the file contains
- "confidence": 0.0-1.0 how confident you are this is a valid download link

Only include links that are likely to work. Skip broken links, login walls, or purchase pages.

Return ONLY the JSON array, no other text."""

        try:
            response = self.client.messages.create(
                model="claude-3-5-haiku-20241022",
                max_tokens=2000,
                messages=[{"role": "user", "content": prompt}],
            )

            # Parse JSON from response
            text = response.content[0].text.strip()
            # Find JSON array in response
            match = re.search(r"\[[\s\S]*\]", text)
            if match:
                return json.loads(match.group())
            return []

        except Exception as e:
            logger.error(f"Claude analysis failed: {e}")
            return []


class BrowserFetcher:
    """Automated browser-based paper fetcher."""

    # Common PDF link patterns on publisher sites
    PDF_PATTERNS = [
        r'href=["\']([^"\']*\.pdf[^"\']*)["\']',
        r'href=["\']([^"\']*download[^"\']*pdf[^"\']*)["\']',
        r'href=["\']([^"\']*fulltext[^"\']*pdf[^"\']*)["\']',
        r'href=["\']([^"\']*article[^"\']*pdf[^"\']*)["\']',
        r'data-pdf-url=["\']([^"\']+)["\']',
    ]

    # Publisher-specific selectors for PDF buttons
    PUBLISHER_SELECTORS = {
        "sciencedirect.com": [
            "a.pdf-download-btn-link",
            'a[data-aa-name="PDF link"]',
            "#pdfLink",
        ],
        "nature.com": [
            "a[data-article-pdf]",
            "a.c-pdf-download__link",
        ],
        "springer.com": [
            "a.c-pdf-download__link",
            'a[data-track-action="Download PDF"]',
        ],
        "wiley.com": [
            "a.pdf-download",
            'a[title="PDF"]',
        ],
        "cell.com": [
            "a.pdf-download",
            'a[data-action="download-pdf"]',
        ],
        "pnas.org": [
            "a.article-dl-pdf-link",
        ],
        "ahajournals.org": [
            "a.article__ctrl--pdf",
        ],
        "oup.com": [  # Oxford Academic
            "a.article-pdf-download",
            'a[data-track="pdf-download"]',
        ],
        "pubmed.ncbi.nlm.nih.gov": [
            "a.link-item.pmc",  # Link to PMC
            "a.id-link",  # DOI link
        ],
    }

    def __init__(
        self,
        downloads_dir: Path,
        target_dir: Path,
        headless: bool = False,
        use_claude: bool = False,
        timeout: int = 60000,
        captcha_wait_time: int = 30,
        fallback_to_manual: bool = True,
        use_profile: bool = False,
        profile_path: Optional[Path] = None,
        slow_mo: int = 0,
    ):
        self.downloads_dir = Path(downloads_dir)
        self.target_dir = Path(target_dir)
        self.headless = headless
        self.timeout = timeout
        self.use_claude = use_claude
        self.captcha_wait_time = captcha_wait_time
        self.fallback_to_manual = fallback_to_manual
        self.use_profile = use_profile
        self.profile_path = profile_path
        self.slow_mo = slow_mo
        self._cloudflare_loop_detected = False  # Flag for loop detection

        self.analyzer = PageAnalyzer() if use_claude else None
        self.browser: Optional[Browser] = None
        self.page: Optional[Page] = None

    def __enter__(self):
        if not PLAYWRIGHT_AVAILABLE:
            raise RuntimeError("Playwright not installed")

        self.playwright = sync_playwright().start()

        # Common launch args for anti-detection
        launch_args = [
            "--disable-blink-features=AutomationControlled",
            "--disable-dev-shm-usage",
            "--no-first-run",
            "--no-default-browser-check",
            "--disable-infobars",
            "--window-size=1920,1080",
            "--start-maximized",
        ]

        if self.use_profile:
            # Use persistent context with real Chrome profile
            # This makes the browser look much more legitimate
            if self.profile_path:
                user_data_dir = str(self.profile_path)
            else:
                # Default Chrome profile location on macOS
                user_data_dir = str(
                    Path.home() / "Library/Application Support/Google/Chrome"
                )

            logger.info(f"  [Profile] Using Chrome profile: {user_data_dir}")

            # Note: launch_persistent_context combines browser + context
            self.context = self.playwright.chromium.launch_persistent_context(
                user_data_dir,
                headless=self.headless,
                downloads_path=str(self.downloads_dir),
                slow_mo=self.slow_mo,
                args=launch_args,
                accept_downloads=True,
                viewport={"width": 1920, "height": 1080},
                locale="en-US",
                timezone_id="America/New_York",
                color_scheme="light",
                channel="chrome",  # Use installed Chrome, not Chromium
            )
            self.browser = None  # Not used with persistent context
            self.page = (
                self.context.pages[0] if self.context.pages else self.context.new_page()
            )
        else:
            # Standard launch with fresh browser
            self.browser = self.playwright.chromium.launch(
                headless=self.headless,
                downloads_path=str(self.downloads_dir),
                slow_mo=self.slow_mo,
                args=launch_args,
            )

            # Create context with more realistic browser fingerprint
            self.context = self.browser.new_context(
                accept_downloads=True,
                user_agent=(
                    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/122.0.0.0 Safari/537.36"
                ),
                viewport={"width": 1920, "height": 1080},
                locale="en-US",
                timezone_id="America/New_York",
                permissions=["geolocation"],
                color_scheme="light",
            )
            self.page = self.context.new_page()

        # Apply stealth patches to mask automation fingerprints
        self._apply_stealth_patches(self.page)

        return self

    def _apply_stealth_patches(self, page: Page):
        """Apply JavaScript patches to mask Playwright automation fingerprints."""
        # These scripts run before page content loads
        stealth_scripts = """
        // Remove webdriver property
        Object.defineProperty(navigator, 'webdriver', {
            get: () => undefined
        });

        // Mock plugins array (real browsers have plugins)
        Object.defineProperty(navigator, 'plugins', {
            get: () => [
                { name: 'Chrome PDF Plugin', filename: 'internal-pdf-viewer' },
                { name: 'Chrome PDF Viewer', filename: 'mhjfbmdgcfjbbpaeojofohoefgiehjai' },
                { name: 'Native Client', filename: 'internal-nacl-plugin' }
            ]
        });

        // Mock languages
        Object.defineProperty(navigator, 'languages', {
            get: () => ['en-US', 'en']
        });

        // Hide automation flags in chrome object
        window.chrome = {
            runtime: {},
            loadTimes: function() {},
            csi: function() {},
            app: {}
        };

        // Mock permissions API
        const originalQuery = window.navigator.permissions.query;
        window.navigator.permissions.query = (parameters) => (
            parameters.name === 'notifications' ?
                Promise.resolve({ state: Notification.permission }) :
                originalQuery(parameters)
        );

        // Add realistic screen properties
        Object.defineProperty(screen, 'availWidth', { get: () => 1920 });
        Object.defineProperty(screen, 'availHeight', { get: () => 1080 });
        Object.defineProperty(screen, 'colorDepth', { get: () => 24 });
        """

        page.add_init_script(stealth_scripts)
        logger.debug("  [Stealth] Applied anti-detection patches")

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.browser:
            self.browser.close()
        if hasattr(self, "playwright"):
            self.playwright.stop()

    def _restart_browser(self):
        """
        Restart the browser after a crash.

        Closes existing browser/context and creates a new one.
        """
        logger.info("  [Recovery] Restarting browser...")

        # Close existing browser if possible
        try:
            if self.page:
                self.page.close()
        except Exception:
            pass

        try:
            if hasattr(self, "context") and self.context:
                self.context.close()
        except Exception:
            pass

        try:
            if self.browser:
                self.browser.close()
        except Exception:
            pass

        # Common launch args for anti-detection
        launch_args = [
            "--disable-blink-features=AutomationControlled",
            "--disable-dev-shm-usage",
            "--no-first-run",
            "--no-default-browser-check",
            "--disable-infobars",
            "--window-size=1920,1080",
            "--start-maximized",
        ]

        # Re-create browser
        self.browser = self.playwright.chromium.launch(
            headless=self.headless,
            downloads_path=str(self.downloads_dir),
            slow_mo=self.slow_mo,
            args=launch_args,
        )

        # Create new context
        self.context = self.browser.new_context(
            accept_downloads=True,
            user_agent=(
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/122.0.0.0 Safari/537.36"
            ),
            viewport={"width": 1920, "height": 1080},
            locale="en-US",
            timezone_id="America/New_York",
            permissions=["geolocation"],
            color_scheme="light",
        )

        # Create new page
        self.page = self.context.new_page()

        # Apply stealth patches
        self._apply_stealth_patches(self.page)

        logger.info("  [Recovery] Browser restarted successfully")

    def _get_existing_files(self) -> set:
        """Get set of existing files in downloads directory."""
        return {f.name for f in self.downloads_dir.iterdir() if f.is_file()}

    def _handle_cookie_consent(self, page: Page) -> bool:
        """
        Detect and dismiss cookie consent dialogs.

        Returns True if a dialog was found and dismissed.
        """
        # Common cookie consent button selectors
        cookie_selectors = [
            'button:has-text("Accept All Cookies")',
            'button:has-text("Accept All")',
            'button:has-text("Accept Cookies")',
            'button:has-text("I Accept")',
            'button:has-text("I Agree")',
            'button:has-text("OK")',
            'button:has-text("Got it")',
            'button:has-text("Allow All")',
            'button:has-text("Allow all cookies")',
            "#onetrust-accept-btn-handler",  # OneTrust
            ".cc-accept",  # Cookie Consent
            ".cookie-accept",
            '[data-testid="cookie-accept"]',
            '[aria-label="Accept cookies"]',
        ]

        for selector in cookie_selectors:
            try:
                button = page.query_selector(selector)
                if button and button.is_visible():
                    logger.info(
                        f"  [Cookies] Found consent dialog, clicking: {selector}"
                    )
                    button.click()
                    time.sleep(1)  # Wait for dialog to close
                    return True
            except Exception:
                continue

        return False

    def _extract_and_save_figures(
        self, page: Page, pmid: str, article_html: str
    ) -> List[Dict]:
        """
        Extract figure images from the page and save them locally.

        Returns list of dicts with figure info (local_path, caption, original_url).
        """
        from urllib.parse import urljoin

        figures_dir = self.target_dir / f"{pmid}_figures"
        figures = []

        try:
            # Find figure images using various selectors
            figure_selectors = [
                "figure img",
                ".figure img",
                ".article-figure img",
                "img.figure",
                "img[class*='figure']",
                "img[src*='figure']",
                "img[src*='fig']",
                ".highres-img",
                "picture img",
            ]

            seen_urls = set()
            img_elements = []

            for selector in figure_selectors:
                try:
                    elements = page.query_selector_all(selector)
                    img_elements.extend(elements)
                except Exception:
                    continue

            if not img_elements:
                logger.info("  [Figures] No figure images found")
                return figures

            # Create figures directory
            figures_dir.mkdir(exist_ok=True)
            logger.info(
                f"  [Figures] Found {len(img_elements)} potential figure images"
            )

            fig_count = 0
            captions_data = []

            for img in img_elements:
                try:
                    src = img.get_attribute("src") or img.get_attribute("data-src")
                    if not src or src in seen_urls:
                        continue
                    seen_urls.add(src)

                    # Skip small icons, logos, etc.
                    if any(
                        skip in src.lower()
                        for skip in ["icon", "logo", "button", "arrow", "spinner"]
                    ):
                        continue

                    # Make absolute URL
                    full_url = urljoin(page.url, src)

                    # Get caption (look for nearby figcaption or alt text)
                    caption = img.get_attribute("alt") or ""
                    try:
                        parent = img.evaluate_handle("el => el.closest('figure')")
                        if parent:
                            figcaption = parent.query_selector("figcaption")
                            if figcaption:
                                caption = figcaption.text_content().strip()
                    except Exception:
                        pass

                    # Determine file extension
                    ext = Path(src.split("?")[0]).suffix.lower()
                    if ext not in [".jpg", ".jpeg", ".png", ".gif", ".webp", ".svg"]:
                        ext = ".jpg"  # Default

                    fig_count += 1
                    filename = f"figure_{fig_count:03d}{ext}"
                    local_path = figures_dir / filename

                    # Download the image
                    try:
                        response = requests.get(
                            full_url,
                            timeout=30,
                            headers={
                                "User-Agent": "Mozilla/5.0 (compatible; research-bot)"
                            },
                        )
                        if response.status_code == 200 and len(response.content) > 1000:
                            local_path.write_bytes(response.content)
                            figures.append(
                                {
                                    "local_path": str(local_path),
                                    "filename": filename,
                                    "caption": caption,
                                    "original_url": full_url,
                                }
                            )
                            captions_data.append(
                                {
                                    "filename": filename,
                                    "caption": caption,
                                    "original_url": full_url,
                                }
                            )
                            logger.info(f"  [Figures] Saved {filename}")
                    except Exception as e:
                        logger.warning(
                            f"  [Figures] Failed to download {full_url}: {e}"
                        )

                except Exception as e:
                    logger.debug(f"  [Figures] Error processing image: {e}")
                    continue

            # Save captions to JSON
            if captions_data:
                captions_path = figures_dir / "captions.json"
                captions_path.write_text(json.dumps(captions_data, indent=2))
                logger.info(
                    f"  [Figures] Saved {len(captions_data)} figures to {figures_dir}"
                )

        except Exception as e:
            logger.error(f"  [Figures] Error extracting figures: {e}")

        return figures

    def _scrape_article_html(self, page: Page, pmid: str) -> Optional[Path]:
        """
        Scrape article content from HTML when PDF is not available.

        For older papers or sites without downloadable PDFs, extract the
        article text from the HTML and save as markdown. Also saves raw HTML
        and extracts figure images.

        Returns path to saved markdown file, or None if extraction failed.
        """
        logger.info("  [HTMLScrape] Attempting to extract article content from HTML...")

        try:
            # Get the full page HTML first (for raw saving)
            full_page_html = page.content()

            # Save raw HTML
            raw_html_path = self.target_dir / f"{pmid}_raw.html"
            raw_html_path.write_text(full_page_html, encoding="utf-8")
            logger.info(f"  [HTMLScrape] Saved raw HTML to {pmid}_raw.html")

            # Common selectors for article content
            content_selectors = [
                "article.article",
                "article",
                ".article-body",
                ".article__body",
                ".article-content",
                ".fulltext",
                ".full-text",
                "#article-content",
                "main.content",
                ".hlFld-Fulltext",  # AHA journals
            ]

            article_html = None
            for selector in content_selectors:
                try:
                    element = page.query_selector(selector)
                    if element:
                        article_html = element.inner_html()
                        if len(article_html) > 1000:  # Meaningful content
                            logger.info(
                                f"  [HTMLScrape] Found content with selector: {selector}"
                            )
                            break
                except Exception:
                    continue

            if not article_html or len(article_html) < 1000:
                # Fallback: get the whole page
                article_html = full_page_html
                logger.info("  [HTMLScrape] Using full page content")

            # Extract and save figures
            figures = self._extract_and_save_figures(page, pmid, article_html)

            # Extract title
            title = ""
            try:
                title_el = page.query_selector(
                    "h1.article-title, h1.citation__title, h1"
                )
                if title_el:
                    title = title_el.text_content().strip()
            except Exception:
                pass

            # Extract abstract
            abstract = ""
            abstract_selectors = [
                ".abstractSection",
                ".abstract",
                "#abstract",
                '[role="doc-abstract"]',
            ]
            for selector in abstract_selectors:
                try:
                    abs_el = page.query_selector(selector)
                    if abs_el:
                        abstract = abs_el.text_content().strip()
                        break
                except Exception:
                    continue

            # Convert HTML to simple text/markdown
            # Remove script and style tags
            import re

            clean_html = re.sub(
                r"<script[^>]*>.*?</script>",
                "",
                article_html,
                flags=re.DOTALL | re.IGNORECASE,
            )
            clean_html = re.sub(
                r"<style[^>]*>.*?</style>",
                "",
                clean_html,
                flags=re.DOTALL | re.IGNORECASE,
            )
            clean_html = re.sub(
                r"<nav[^>]*>.*?</nav>", "", clean_html, flags=re.DOTALL | re.IGNORECASE
            )
            clean_html = re.sub(
                r"<header[^>]*>.*?</header>",
                "",
                clean_html,
                flags=re.DOTALL | re.IGNORECASE,
            )
            clean_html = re.sub(
                r"<footer[^>]*>.*?</footer>",
                "",
                clean_html,
                flags=re.DOTALL | re.IGNORECASE,
            )

            # Convert common HTML to markdown-ish
            clean_html = re.sub(r"<h1[^>]*>", "\n# ", clean_html)
            clean_html = re.sub(r"</h1>", "\n", clean_html)
            clean_html = re.sub(r"<h2[^>]*>", "\n## ", clean_html)
            clean_html = re.sub(r"</h2>", "\n", clean_html)
            clean_html = re.sub(r"<h3[^>]*>", "\n### ", clean_html)
            clean_html = re.sub(r"</h3>", "\n", clean_html)
            clean_html = re.sub(r"<p[^>]*>", "\n", clean_html)
            clean_html = re.sub(r"</p>", "\n", clean_html)
            clean_html = re.sub(r"<br\s*/?>", "\n", clean_html)
            clean_html = re.sub(r"<li[^>]*>", "\n- ", clean_html)

            # Remove remaining HTML tags
            clean_text = re.sub(r"<[^>]+>", "", clean_html)

            # Clean up whitespace
            clean_text = re.sub(r"\n\s*\n\s*\n+", "\n\n", clean_text)
            clean_text = re.sub(r"[ \t]+", " ", clean_text)
            clean_text = clean_text.strip()

            if len(clean_text) < 500:
                logger.info(
                    f"  [HTMLScrape] Content too short ({len(clean_text)} chars), skipping"
                )
                return None

            # Build markdown content
            markdown_content = f"# {title}\n\n" if title else ""
            markdown_content += f"**Source URL:** {page.url}\n\n"
            markdown_content += f"**PMID:** {pmid}\n\n"
            if abstract:
                markdown_content += f"## Abstract\n\n{abstract}\n\n"
            markdown_content += "## Full Text\n\n"
            markdown_content += clean_text

            # Save to file
            filename = f"PMID_{pmid}_FULL_CONTEXT.md"
            dest_path = self.target_dir / filename
            dest_path.write_text(markdown_content, encoding="utf-8")

            logger.info(
                f"  [HTMLScrape] SUCCESS: Saved {len(markdown_content)} chars to {filename}"
            )
            return dest_path

        except Exception as e:
            logger.error(f"  [HTMLScrape] ERROR: {e}")
            return None

    def _wait_for_new_file(
        self, before_files: set, timeout: int = 60
    ) -> Optional[Path]:
        """Wait for a new file to appear in downloads directory."""
        start = time.time()
        while time.time() - start < timeout:
            current_files = {
                f.name for f in self.downloads_dir.iterdir() if f.is_file()
            }
            new_files = current_files - before_files

            # Filter out partial downloads
            complete_files = [
                f for f in new_files if not f.endswith((".crdownload", ".part", ".tmp"))
            ]

            if complete_files:
                # Return the newest file
                newest = max(
                    complete_files,
                    key=lambda f: (self.downloads_dir / f).stat().st_mtime,
                )
                return self.downloads_dir / newest

            time.sleep(0.5)

        return None

    def _manual_download_fallback(
        self, pmid: str, timeout: int = 300
    ) -> Optional[DownloadResult]:
        """
        Fallback to manual download when Cloudflare blocks automation.

        Waits for user to manually download files in the browser window,
        monitors the downloads folder, and organizes the files.

        Args:
            pmid: PubMed ID for file naming
            timeout: Maximum time to wait for manual download (default 5 minutes)

        Returns:
            DownloadResult if files were downloaded, None if timed out/skipped
        """
        if self.headless:
            logger.warning(
                "  [ManualFallback] Cannot use manual fallback in headless mode"
            )
            return None

        logger.info("")
        logger.info("=" * 60)
        logger.info("  MANUAL DOWNLOAD REQUIRED")
        logger.info("=" * 60)
        logger.info(f"  PMID: {pmid}")
        logger.info("  The browser window is open - please download the PDF manually:")
        logger.info("    1. Navigate to the PDF download link")
        logger.info("    2. Click download (files go to your Downloads folder)")
        logger.info("    3. Files will be auto-detected and organized")
        logger.info(f"  Waiting up to {timeout}s for downloads...")
        logger.info("=" * 60)

        # Take snapshot of current files
        before_files = self._get_existing_files()
        start = time.time()
        last_check_count = 0
        downloaded_files = []

        while time.time() - start < timeout:
            elapsed = int(time.time() - start)

            # Check for new files
            current_files = {
                f.name for f in self.downloads_dir.iterdir() if f.is_file()
            }
            new_files = current_files - before_files

            # Filter out partial downloads
            complete_files = [
                f
                for f in new_files
                if not f.endswith((".crdownload", ".part", ".tmp", ".download"))
                and not f.startswith(".")
            ]

            if len(complete_files) > last_check_count:
                # New file detected
                for filename in complete_files:
                    if filename not in [Path(f).name for f in downloaded_files]:
                        src = self.downloads_dir / filename
                        ext = src.suffix.lower()

                        # Determine new filename
                        if ext == ".pdf" and not any(
                            "Main_Text" in f for f in downloaded_files
                        ):
                            new_name = f"{pmid}_Main_Text.pdf"
                        elif ext == ".pdf":
                            idx = len([f for f in downloaded_files if ".pdf" in f])
                            new_name = f"{pmid}_Supplement_{idx}.pdf"
                        elif ext in (".xlsx", ".xls", ".csv"):
                            idx = len(
                                [
                                    f
                                    for f in downloaded_files
                                    if any(e in f for e in [".xlsx", ".xls", ".csv"])
                                ]
                            )
                            new_name = f"{pmid}_Supp_Data_{idx}{ext}"
                        else:
                            idx = len(downloaded_files)
                            new_name = f"{pmid}_file_{idx}{ext}"

                        dest = self.target_dir / new_name
                        try:
                            src.rename(dest)
                            logger.info(
                                f"  [ManualFallback] Downloaded: {filename} -> {new_name}"
                            )
                            downloaded_files.append(str(dest))
                        except Exception as e:
                            logger.error(
                                f"  [ManualFallback] Failed to move {filename}: {e}"
                            )

                last_check_count = len(complete_files)

            # Progress indicator every 30 seconds
            if elapsed > 0 and elapsed % 30 == 0:
                if downloaded_files:
                    logger.info(
                        f"  [ManualFallback] {elapsed}s elapsed, {len(downloaded_files)} file(s) downloaded"
                    )
                else:
                    logger.info(
                        f"  [ManualFallback] {elapsed}s/{timeout}s - waiting for downloads..."
                    )

            # If we have files and no new files for 10 seconds, consider done
            if downloaded_files and elapsed > 10:
                # Check if files are still being added
                time.sleep(5)
                current_files2 = {
                    f.name for f in self.downloads_dir.iterdir() if f.is_file()
                }
                new_files2 = current_files2 - before_files
                complete_files2 = [
                    f
                    for f in new_files2
                    if not f.endswith((".crdownload", ".part", ".tmp", ".download"))
                    and not f.startswith(".")
                ]
                if len(complete_files2) == len(complete_files):
                    # No new files in 5 seconds, assume done
                    logger.info(
                        f"  [ManualFallback] Download complete! {len(downloaded_files)} file(s)"
                    )
                    break

            time.sleep(1)

        if downloaded_files:
            return DownloadResult(
                pmid=pmid,
                success=True,
                files_downloaded=downloaded_files,
                method="manual_fallback",
            )
        else:
            logger.warning(f"  [ManualFallback] No files downloaded after {timeout}s")
            return DownloadResult(
                pmid=pmid,
                success=False,
                files_downloaded=[],
                error="Manual download timed out - no files detected",
                method="manual_fallback_timeout",
                failure_type="timeout",
            )

    def _save_embedded_pdf(self, page: Page, pmid: str) -> Optional[Path]:
        """
        Save a PDF that is rendered inline in the browser.

        When a PDF is loaded directly in the browser (content-type: application/pdf),
        we need to fetch it via the response or use JS to save it.
        """
        try:
            url = page.url
            logger.info(f"  [EmbeddedPDF] Attempting to save from: {url[:60]}...")

            # Method 1: Use requests to fetch the URL now that browser has established session
            # Get cookies from browser to use with requests
            import requests as req

            cookies = page.context.cookies()
            cookie_dict = {c["name"]: c["value"] for c in cookies}

            headers = {
                "User-Agent": (
                    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/120.0.0.0 Safari/537.36"
                ),
            }

            response = req.get(url, headers=headers, cookies=cookie_dict, timeout=60)
            logger.info(
                f"  [EmbeddedPDF] Response: HTTP {response.status_code}, {len(response.content)} bytes"
            )

            if response.status_code == 200 and response.content[:4] == b"%PDF":
                filename = url.split("/")[-1]
                if not filename.endswith(".pdf"):
                    filename = f"{pmid}_Main_Text.pdf"
                else:
                    filename = f"{pmid}_{filename}"

                dest_path = self.target_dir / filename
                dest_path.write_bytes(response.content)
                logger.info(f"  [EmbeddedPDF] SUCCESS: {dest_path.name}")
                return dest_path

            logger.info("  [EmbeddedPDF] Response is not a valid PDF")

        except Exception as e:
            logger.info(f"  [EmbeddedPDF] ERROR: {e}")

        return None

    def _wait_for_cloudflare(
        self,
        page: Page,
        max_wait: Optional[int] = None,
        has_early_content: bool = False,
    ) -> bool:
        """
        Detect and wait for Cloudflare challenge to be completed.

        Returns True if page is ready (no Cloudflare or challenge passed).
        Returns False if still blocked after max_wait seconds or in infinite loop.

        Args:
            page: The Playwright page object
            max_wait: Maximum seconds to wait for challenge to clear
            has_early_content: If True, we already captured content before hitting
                              Cloudflare, so we can exit faster when stuck
        """
        if max_wait is None:
            max_wait = self.captcha_wait_time
        start = time.time()
        check_count = 0
        clear_count = 0  # Track consecutive clear checks for stability
        initial_url = page.url
        cloudflare_detected = False
        challenge_solve_count = 0  # Track how many times we've "solved" a challenge
        last_challenge_time = None  # Track when we last saw a challenge
        urls_seen = set()  # Track URL changes to detect redirect loops
        stuck_on_challenge_start = None  # Track when we first saw current challenge
        STUCK_THRESHOLD = (
            45  # Exit after 45 seconds stuck on same challenge if we have content
        )

        while time.time() - start < max_wait:
            check_count += 1
            # Check for common Cloudflare indicators
            try:
                page_content = page.content().lower()
                title = page.title().lower() if page.title() else ""
                current_url = page.url
                urls_seen.add(current_url)
            except Exception as e:
                logger.info(f"  [Cloudflare] Error getting page content: {e}")
                clear_count = 0  # Reset stability counter on error
                time.sleep(2)
                continue

            cloudflare_indicators = [
                "checking your browser" in page_content,
                "just a moment" in title,
                "cf-browser-verification" in page_content,
                "challenge-running" in page_content,
                "turnstile" in page_content,  # Cloudflare Turnstile widget
                "challenge-form" in page_content,  # Additional indicator
                # Additional CAPTCHA providers and patterns
                "g-recaptcha" in page_content,  # Google reCAPTCHA
                "recaptcha" in page_content and "challenge" in page_content,
                "h-captcha" in page_content,  # hCaptcha
                "hcaptcha" in page_content,
                "verify you are human" in page_content,
                "verify you're human" in page_content,
                "human verification" in page_content,
                "bot detection" in page_content,
                "security check" in title,
                "access denied" in title and "verify" in page_content,
                "please verify" in page_content,
                "complete the captcha" in page_content,
                "solve the captcha" in page_content,
            ]

            # Check for actual content indicators (means we're past Cloudflare)
            content_indicators = [
                "<article" in page_content,
                "<main" in page_content,
                "abstract" in page_content and len(page_content) > 5000,
                "pdf" in page_content and "download" in page_content,
                "full text" in page_content,
            ]

            if any(cloudflare_indicators):
                now = time.time()

                # Track how long we've been stuck on the same challenge
                if stuck_on_challenge_start is None:
                    stuck_on_challenge_start = now
                elif has_early_content:
                    stuck_time = now - stuck_on_challenge_start
                    if stuck_time > STUCK_THRESHOLD:
                        logger.warning(
                            f"  [Cloudflare] Stuck on challenge for {int(stuck_time)}s with no progress."
                        )
                        logger.info(
                            "  [Cloudflare] Early content capture exists - exiting to use that instead."
                        )
                        self._cloudflare_loop_detected = True
                        return False

                # Log which type of challenge was detected
                if "g-recaptcha" in page_content or (
                    "recaptcha" in page_content and "challenge" in page_content
                ):
                    logger.info("  [Cloudflare] Detected: Google reCAPTCHA")
                elif "h-captcha" in page_content or "hcaptcha" in page_content:
                    logger.info("  [Cloudflare] Detected: hCaptcha")
                elif "turnstile" in page_content:
                    logger.info("  [Cloudflare] Detected: Cloudflare Turnstile")
                elif (
                    "cf-browser-verification" in page_content
                    or "challenge-running" in page_content
                ):
                    logger.info("  [Cloudflare] Detected: Cloudflare challenge")
                elif any(
                    x in page_content
                    for x in [
                        "verify you are human",
                        "verify you're human",
                        "human verification",
                    ]
                ):
                    logger.info("  [Cloudflare] Detected: Human verification challenge")
                else:
                    logger.info("  [Cloudflare] Detected: Generic security challenge")
                # If we were previously clear and now see challenge again, we're in a loop
                if clear_count > 0:
                    challenge_solve_count += 1
                    # Detect loop after just 2 cycles (reduced from 3)
                    if challenge_solve_count >= 2:
                        elapsed = int(now - start)
                        logger.warning(
                            f"  [Cloudflare] INFINITE LOOP DETECTED! "
                            f"({challenge_solve_count} cycles, {elapsed}s elapsed)"
                        )
                        logger.warning(
                            "  [Cloudflare] The site keeps re-challenging after CAPTCHA solve."
                        )
                        logger.warning("  [Cloudflare] This usually means:")
                        logger.warning(
                            "  [Cloudflare]   1. The site blocks automated browsers"
                        )
                        logger.warning(
                            "  [Cloudflare]   2. Session cookies aren't persisting"
                        )
                        logger.warning(
                            "  [Cloudflare]   3. Will attempt manual download fallback..."
                        )
                        self._cloudflare_loop_detected = True
                        return False
                    logger.warning(
                        f"  [Cloudflare] Challenge REAPPEARED after solving! (cycle {challenge_solve_count}/2)"
                    )
                    logger.info(
                        "  [Cloudflare] If this keeps happening, the site may block automation"
                    )

                cloudflare_detected = True
                clear_count = 0  # Reset stability counter
                last_challenge_time = now
                elapsed = int(now - start)
                logger.info(
                    f"  [Cloudflare] Challenge detected, waiting for solve... ({elapsed}s/{max_wait}s)"
                )
                time.sleep(2)
                continue

            # No Cloudflare indicators detected
            clear_count += 1
            stuck_on_challenge_start = None  # Reset stuck timer when clear

            # If we detected Cloudflare earlier and now see content, we passed
            if cloudflare_detected and any(content_indicators):
                # Wait longer (5s instead of 3s) to catch reload loops
                logger.info(
                    "  [Cloudflare] Content detected, verifying stability (5s)..."
                )
                time.sleep(5)

                # Re-check for Cloudflare after the delay
                try:
                    recheck_content = page.content().lower()
                    recheck_title = page.title().lower() if page.title() else ""
                    recheck_indicators = [
                        "checking your browser" in recheck_content,
                        "just a moment" in recheck_title,
                        "cf-browser-verification" in recheck_content,
                        "challenge-running" in recheck_content,
                        "turnstile" in recheck_content,
                        "challenge-form" in recheck_content,
                        "g-recaptcha" in recheck_content,
                        "recaptcha" in recheck_content
                        and "challenge" in recheck_content,
                        "h-captcha" in recheck_content,
                        "hcaptcha" in recheck_content,
                        "verify you are human" in recheck_content,
                        "verify you're human" in recheck_content,
                        "human verification" in recheck_content,
                        "security check" in recheck_title,
                        "please verify" in recheck_content,
                    ]
                    if any(recheck_indicators):
                        logger.warning(
                            "  [Cloudflare] Challenge reappeared after reload!"
                        )
                        clear_count = 0
                        # This counts as a new challenge cycle
                        challenge_solve_count += 1
                        if challenge_solve_count >= 2:
                            logger.warning(
                                "  [Cloudflare] INFINITE LOOP - will try manual download."
                            )
                            self._cloudflare_loop_detected = True
                            return False
                        continue
                except Exception:
                    pass

                logger.info(
                    "  [Cloudflare] Challenge passed! Content loaded successfully."
                )
                return True

            # Require multiple consecutive clear checks for stability
            # This prevents false positives during page transitions
            if clear_count >= 5:
                if cloudflare_detected:
                    # Extra verification: wait and re-check one more time
                    logger.info(
                        "  [Cloudflare] Appears clear, final verification (5s)..."
                    )
                    time.sleep(5)
                    try:
                        final_content = page.content().lower()
                        final_title = page.title().lower() if page.title() else ""
                        final_indicators = [
                            "checking your browser" in final_content,
                            "just a moment" in final_title,
                            "cf-browser-verification" in final_content,
                            "challenge-running" in final_content,
                            "turnstile" in final_content,
                            "challenge-form" in final_content,
                            "g-recaptcha" in final_content,
                            "recaptcha" in final_content
                            and "challenge" in final_content,
                            "h-captcha" in final_content,
                            "hcaptcha" in final_content,
                            "verify you are human" in final_content,
                            "verify you're human" in final_content,
                            "human verification" in final_content,
                            "security check" in final_title,
                            "please verify" in final_content,
                        ]
                        if any(final_indicators):
                            logger.warning(
                                "  [Cloudflare] Challenge reappeared on final check!"
                            )
                            clear_count = 0
                            challenge_solve_count += 1
                            if challenge_solve_count >= 2:
                                logger.warning(
                                    "  [Cloudflare] INFINITE LOOP - will try manual download."
                                )
                                self._cloudflare_loop_detected = True
                                return False
                            continue
                    except Exception:
                        pass

                if check_count > 1:
                    logger.info("  [Cloudflare] Passed or not present")
                self._cloudflare_loop_detected = False
                return True

            # Still checking for stability
            time.sleep(1)

        logger.warning(
            f"  [Cloudflare] Timed out after {max_wait}s waiting for challenge"
        )
        return False

    def _browser_download_pdf(self, url: str, pmid: str) -> Optional[Path]:
        """Download PDF using browser (handles JS redirects like PMC)."""
        before_files = self._get_existing_files()
        logger.info(f"  [Browser] Attempting download: {url[:80]}...")

        try:
            # Navigate and wait for download to start
            with self.page.expect_download(timeout=30000) as download_info:
                self.page.goto(url, wait_until="commit", timeout=self.timeout)

            download = download_info.value
            logger.info(
                f"  [Browser] Download triggered, filename: {download.suggested_filename}"
            )

            # Save with PMID prefix
            suggested = download.suggested_filename
            if not suggested.endswith(".pdf"):
                suggested = f"{pmid}_Main_Text.pdf"
            else:
                suggested = f"{pmid}_{suggested}"

            dest_path = self.target_dir / suggested
            download.save_as(str(dest_path))

            # Verify it's actually a PDF
            with open(dest_path, "rb") as f:
                header = f.read(4)
            if header == b"%PDF":
                logger.info(f"  [Browser] SUCCESS: {dest_path.name}")
                return dest_path
            else:
                logger.info(
                    f"  [Browser] Downloaded file is not a valid PDF (header: {header!r})"
                )
                dest_path.unlink()  # Delete invalid file
                return None

        except PlaywrightTimeout:
            logger.info(
                "  [Browser] Timeout waiting for download, checking for Cloudflare..."
            )

            # Check for Cloudflare challenge
            if self._wait_for_cloudflare(self.page):
                logger.info(
                    "  [Browser] Cloudflare passed or not present, retrying download..."
                )
                # Try to get the PDF now that Cloudflare is cleared
                try:
                    # Check if we're now on a PDF page
                    content_type = self.page.evaluate("() => document.contentType")
                    logger.info(f"  [Browser] Page content type: {content_type}")

                    if content_type == "application/pdf":
                        # The PDF is rendered inline in browser, try to save it
                        result = self._save_embedded_pdf(self.page, pmid)
                        if result:
                            return result

                        # Fallback: Try clicking any download button that might have appeared
                        result = self._click_download_button(self.page, pmid)
                        if result:
                            return result
                except Exception as e:
                    logger.debug(f"  [Browser] Post-Cloudflare check failed: {e}")

            # Check if a file appeared in downloads anyway
            new_file = self._wait_for_new_file(before_files, timeout=10)
            if new_file and new_file.suffix.lower() == ".pdf":
                logger.info(f"  [Browser] Found file after timeout: {new_file.name}")
                # Verify and move
                with open(new_file, "rb") as f:
                    header = f.read(4)
                if header == b"%PDF":
                    dest_path = self.target_dir / f"{pmid}_{new_file.name}"
                    new_file.rename(dest_path)
                    logger.info(f"  [Browser] SUCCESS: {dest_path.name}")
                    return dest_path
            logger.info("  [Browser] No valid PDF found after timeout")
            return None

        except Exception as e:
            logger.info(f"  [Browser] ERROR: {e}")
            return None

    def _download_pdf_direct(self, url: str, pmid: str) -> Optional[Path]:
        """Try to download a PDF directly using requests."""
        import requests as req

        logger.info(f"  [Direct] Attempting: {url[:80]}...")
        try:
            headers = {
                "User-Agent": (
                    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/120.0.0.0 Safari/537.36"
                ),
            }
            response = req.get(url, headers=headers, timeout=30)
            logger.info(
                f"  [Direct] Response: HTTP {response.status_code}, {len(response.content)} bytes"
            )

            if response.status_code == 200:
                content = response.content
                content_type = response.headers.get("Content-Type", "unknown")
                logger.info(f"  [Direct] Content-Type: {content_type}")

                # Check if it's actually a PDF (PDF files start with %PDF)
                if content[:4] == b"%PDF":
                    # Extract filename from URL or use default
                    filename = url.split("/")[-1]
                    if not filename.endswith(".pdf"):
                        filename = f"{pmid}_Main_Text.pdf"
                    else:
                        filename = f"{pmid}_{filename}"

                    dest_path = self.target_dir / filename
                    with open(dest_path, "wb") as f:
                        f.write(content)

                    logger.info(f"  [Direct] SUCCESS: {dest_path.name}")
                    return dest_path
                else:
                    header_preview = content[:50].decode("utf-8", errors="replace")
                    logger.info(f"  [Direct] Not a PDF. Header: {header_preview!r}")
            else:
                logger.info(f"  [Direct] FAILED: HTTP {response.status_code}")

        except Exception as e:
            logger.info(f"  [Direct] ERROR: {e}")

        return None

    def _download_supplements(self, page: Page, pmid: str) -> List[Path]:
        """
        Find and download supplement files from the current page.

        Uses SupplementScraper to find supplement links, then downloads them.
        """
        if not SCRAPER_AVAILABLE:
            logger.info(
                "  [Supplements] SupplementScraper not available (harvesting module not found)"
            )
            return []

        downloaded = []
        scraper = SupplementScraper()
        html = page.content()
        base_url = page.url
        logger.info(f"  [Supplements] Scanning page: {base_url[:60]}...")

        # Determine which scraper to use based on domain
        domain = base_url.lower()
        if "nature.com" in domain:
            logger.info("  [Supplements] Using Nature scraper")
            supplements = scraper.scrape_nature_supplements(html, base_url)
        elif any(
            d in domain for d in ["sciencedirect.com", "elsevier.com", "cell.com"]
        ):
            logger.info("  [Supplements] Using Elsevier scraper")
            supplements = scraper.scrape_elsevier_supplements(html, base_url)
        else:
            logger.info("  [Supplements] Using generic scraper")
            supplements = scraper.scrape_generic_supplements(html, base_url)

        if not supplements:
            logger.info("  [Supplements] No supplement links found on page")
            return []

        logger.info(f"  [Supplements] Found {len(supplements)} potential supplement(s)")

        for idx, supp in enumerate(supplements, 1):
            supp_url = supp.get("url", "")
            supp_name = supp.get("name", f"supplement_{idx}")

            if not supp_url:
                continue

            # Clean up filename
            if not any(
                supp_name.lower().endswith(ext)
                for ext in [".pdf", ".xlsx", ".xls", ".docx", ".doc", ".csv", ".zip"]
            ):
                # Try to get extension from URL
                url_path = supp_url.split("?")[0]
                if "." in url_path.split("/")[-1]:
                    supp_name = url_path.split("/")[-1]

            dest_name = f"{pmid}_Supp_{idx}_{supp_name}"
            dest_path = self.target_dir / dest_name

            # Try direct download first
            result = self._download_file_direct(supp_url, dest_path)
            if result:
                downloaded.append(result)
                logger.info(f"    Downloaded supplement: {dest_name}")
            else:
                logger.debug(f"    Failed to download: {supp_url}")

        return downloaded

    def _download_file_direct(self, url: str, dest_path: Path) -> Optional[Path]:
        """Download a file directly using requests."""
        import requests as req

        try:
            headers = {
                "User-Agent": (
                    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/120.0.0.0 Safari/537.36"
                ),
            }
            response = req.get(url, headers=headers, timeout=60, allow_redirects=True)

            if response.status_code == 200:
                content = response.content

                # Skip if it's an HTML error page
                if content[:5] == b"<!DOC" or content[:5] == b"<html":
                    logger.debug(f"    Got HTML instead of file for {url}")
                    return None

                dest_path.write_bytes(content)
                return dest_path

        except Exception as e:
            logger.debug(f"    Download failed: {e}")

        return None

    def _find_pdf_links_heuristic(self, page: Page) -> List[str]:
        """Find PDF links using heuristics (no Claude)."""
        links = []
        content = page.content()
        current_url = page.url
        logger.info(f"  [Heuristic] Scanning page: {current_url[:60]}...")

        # Try regex patterns
        for pattern in self.PDF_PATTERNS:
            matches = re.findall(pattern, content, re.IGNORECASE)
            if matches:
                logger.info(
                    f"  [Heuristic] Pattern '{pattern[:30]}...' found {len(matches)} match(es)"
                )
            links.extend(matches)

        # Try publisher-specific selectors
        for domain, selectors in self.PUBLISHER_SELECTORS.items():
            if domain in current_url:
                logger.info(f"  [Heuristic] Using selectors for {domain}")
                for selector in selectors:
                    try:
                        elements = page.query_selector_all(selector)
                        for el in elements:
                            href = el.get_attribute("href")
                            if href:
                                logger.info(
                                    f"  [Heuristic] Selector '{selector}' found: {href[:60]}..."
                                )
                                links.append(href)
                    except Exception as e:
                        logger.debug(f"  [Heuristic] Selector '{selector}' failed: {e}")

        # Make URLs absolute
        base_url = "/".join(current_url.split("/")[:3])
        absolute_links = []
        for link in links:
            if link.startswith("http"):
                absolute_links.append(link)
            elif link.startswith("//"):
                absolute_links.append("https:" + link)
            elif link.startswith("/"):
                absolute_links.append(base_url + link)
            else:
                absolute_links.append(current_url.rsplit("/", 1)[0] + "/" + link)

        unique_links = list(set(absolute_links))
        logger.info(f"  [Heuristic] Found {len(unique_links)} unique PDF link(s)")
        for i, link in enumerate(unique_links[:5]):
            logger.info(f"  [Heuristic]   {i + 1}. {link[:70]}...")

        return unique_links

    def _try_download_link(self, page: Page, url: str, pmid: str) -> Optional[Path]:
        """Try to download a file from a URL."""
        before_files = self._get_existing_files()

        try:
            # Set up download handler
            with page.expect_download(timeout=self.timeout) as download_info:
                page.goto(url, wait_until="commit", timeout=self.timeout)

            download = download_info.value

            # Wait for download to complete
            dest_path = self.target_dir / f"{pmid}_{download.suggested_filename}"
            download.save_as(str(dest_path))
            logger.info(f"  Downloaded: {dest_path.name}")
            return dest_path

        except PlaywrightTimeout:
            # Check if file appeared anyway
            new_file = self._wait_for_new_file(before_files, timeout=5)
            if new_file:
                # Move to target
                dest_path = self.target_dir / f"{pmid}_{new_file.name}"
                new_file.rename(dest_path)
                return dest_path
            return None

        except Exception as e:
            logger.debug(f"  Download attempt failed: {e}")
            return None

    def _click_download_button(self, page: Page, pmid: str) -> Optional[Path]:
        """Try to click a download button on the page."""
        before_files = self._get_existing_files()
        logger.info("  [Button] Scanning for download buttons...")

        # Common download button selectors
        button_selectors = [
            'a:has-text("PDF")',
            'a:has-text("Download PDF")',
            'a:has-text("Full Text PDF")',
            'button:has-text("PDF")',
            "a[download]",
            ".pdf-download",
            ".download-pdf",
        ]

        for selector in button_selectors:
            try:
                button = page.query_selector(selector)
                if button:
                    is_visible = button.is_visible()
                    button_text = (
                        button.text_content()[:30]
                        if button.text_content()
                        else "(no text)"
                    )
                    logger.info(
                        f"  [Button] Found '{selector}' - visible: {is_visible}, text: {button_text}"
                    )

                    if is_visible:
                        logger.info(f"  [Button] Clicking: {selector}")

                        with page.expect_download(
                            timeout=self.timeout
                        ) as download_info:
                            button.click()

                        download = download_info.value
                        logger.info(
                            f"  [Button] Download started: {download.suggested_filename}"
                        )
                        dest_path = (
                            self.target_dir / f"{pmid}_{download.suggested_filename}"
                        )
                        download.save_as(str(dest_path))
                        logger.info(f"  [Button] SUCCESS: {dest_path.name}")
                        return dest_path

            except PlaywrightTimeout:
                logger.info(
                    f"  [Button] Timeout after clicking '{selector}', checking for CAPTCHA..."
                )

                # Check if a CAPTCHA appeared after clicking
                if self._wait_for_cloudflare(page, max_wait=60):
                    logger.info(
                        "  [Button] CAPTCHA solved or not present, retrying download..."
                    )
                    # After CAPTCHA, check if we're on a PDF page or can retry
                    try:
                        content_type = page.evaluate("() => document.contentType")
                        if content_type == "application/pdf":
                            result = self._save_embedded_pdf(page, pmid)
                            if result:
                                return result
                    except Exception:
                        pass

                    # Try clicking the button again if we're back on the article page
                    try:
                        button = page.query_selector(selector)
                        if button and button.is_visible():
                            logger.info(
                                f"  [Button] Re-clicking after CAPTCHA: {selector}"
                            )
                            with page.expect_download(timeout=30000) as download_info:
                                button.click()
                            download = download_info.value
                            dest_path = (
                                self.target_dir
                                / f"{pmid}_{download.suggested_filename}"
                            )
                            download.save_as(str(dest_path))
                            logger.info(
                                f"  [Button] SUCCESS after CAPTCHA: {dest_path.name}"
                            )
                            return dest_path
                    except Exception as e:
                        logger.debug(f"  [Button] Retry after CAPTCHA failed: {e}")

                # Check for new files
                new_file = self._wait_for_new_file(before_files, timeout=5)
                if new_file:
                    dest_path = self.target_dir / f"{pmid}_{new_file.name}"
                    new_file.rename(dest_path)
                    logger.info(f"  [Button] SUCCESS (late): {dest_path.name}")
                    return dest_path
            except Exception as e:
                logger.debug(f"  [Button] Selector '{selector}' failed: {e}")
                continue

        logger.info("  [Button] No working download buttons found")
        return None

    def fetch_paper(self, pmid: str, url: str) -> DownloadResult:
        """
        Attempt to download a paper using browser automation.

        Args:
            pmid: PubMed ID
            url: URL to the paper (DOI or publisher page)

        Returns:
            DownloadResult with success status and downloaded files
        """
        logger.info(f"\n{'=' * 60}")
        logger.info(f"Fetching PMID {pmid}")
        logger.info(f"URL: {url}")
        logger.info(f"{'=' * 60}")
        downloaded_files = []

        # Check if we already have this paper from a previous run
        existing_md = self.target_dir / f"PMID_{pmid}_FULL_CONTEXT.md"
        if existing_md.exists() and existing_md.stat().st_size > 1000:
            logger.info(
                f">>> ALREADY EXISTS: {existing_md.name} ({existing_md.stat().st_size} bytes)"
            )
            logger.info(">>> Skipping download - marking as success")
            return DownloadResult(
                pmid=pmid,
                success=True,
                files_downloaded=[str(existing_md)],
                method="already_exists",
            )

        try:
            # Check if URL is a direct PDF link
            if url.lower().endswith(".pdf") or "/pdf/" in url.lower():
                logger.info(">>> STEP 1: Direct PDF URL detected")
                # Try direct download with requests first
                result = self._download_pdf_direct(url, pmid)
                if result:
                    return DownloadResult(
                        pmid=pmid,
                        success=True,
                        files_downloaded=[str(result)],
                        method="direct_pdf",
                    )

                # If direct download failed, try browser-based (handles JS redirects)
                logger.info(
                    ">>> STEP 1b: Trying browser-based download for JS-redirect PDF..."
                )
                result = self._browser_download_pdf(url, pmid)
                if result:
                    return DownloadResult(
                        pmid=pmid,
                        success=True,
                        files_downloaded=[str(result)],
                        method="browser_pdf",
                    )

            # Navigate to the paper - use 'load' instead of 'networkidle' for faster response
            logger.info(">>> STEP 2: Navigating to page...")
            try:
                self.page.goto(url, wait_until="load", timeout=self.timeout)
                logger.info(f"  Page loaded. Final URL: {self.page.url[:70]}...")
            except PlaywrightTimeout:
                # Page might still be usable even if not fully loaded
                logger.warning("  Page load timeout, continuing anyway...")

            # EARLY HTML CAPTURE: Scrape article content BEFORE cookies/Cloudflare
            # Some sites show content behind cookie dialogs that gets blocked later
            early_html_result = None
            try:
                early_content = self.page.content()
                if len(early_content) > 10000 and "abstract" in early_content.lower():
                    logger.info(
                        "  [EarlyCapture] Article content detected, saving as backup..."
                    )
                    early_html_result = self._scrape_article_html(self.page, pmid)
            except Exception as e:
                logger.debug(f"  [EarlyCapture] Failed: {e}")

            # Handle cookie consent dialogs (common on publisher sites)
            time.sleep(1)  # Brief wait for dialog to appear
            self._handle_cookie_consent(self.page)

            # Check for Cloudflare and wait if needed
            # Reset the loop flag before checking
            self._cloudflare_loop_detected = False
            has_early_content = early_html_result is not None
            if has_early_content:
                logger.info(
                    "  [EarlyCapture] Content already saved - will exit faster if CAPTCHA blocks us"
                )
            if not self._wait_for_cloudflare(
                self.page, has_early_content=has_early_content
            ):
                # Cloudflare blocked us - but check if we captured content early
                if early_html_result:
                    logger.info(
                        ">>> Cloudflare blocked, but early HTML capture succeeded!"
                    )
                    return DownloadResult(
                        pmid=pmid,
                        success=True,
                        files_downloaded=[str(early_html_result)],
                        method="early_html_capture",
                    )

                # Check if we're stuck in an infinite loop
                if (
                    self._cloudflare_loop_detected
                    and self.fallback_to_manual
                    and not self.headless
                ):
                    # Try manual download fallback
                    logger.info(
                        ">>> Cloudflare loop detected - switching to manual download mode"
                    )
                    result = self._manual_download_fallback(pmid)
                    if result:
                        return result
                    # If manual fallback failed/timed out, return failure
                    return DownloadResult(
                        pmid=pmid,
                        success=False,
                        files_downloaded=[],
                        error="Manual download fallback timed out",
                        method="manual_fallback_timeout",
                        failure_type="captcha_blocked",
                    )
                else:
                    # No fallback available - abort
                    logger.warning(
                        "  Cloudflare challenge not passed - aborting this paper"
                    )
                    return DownloadResult(
                        pmid=pmid,
                        success=False,
                        files_downloaded=[],
                        error="Cloudflare blocks automation. Use non-headless mode with --fallback-manual for manual download.",
                        method="cloudflare_blocked",
                        failure_type="captcha_blocked",
                    )

            time.sleep(3)  # Wait for dynamic content

            # Handle common redirects (PubMed -> publisher)
            current_url = self.page.url
            logger.info(f"  Current URL after redirect: {current_url[:70]}...")

            if "pubmed.ncbi.nlm.nih.gov" in current_url:
                logger.info("  On PubMed page, looking for full text link...")
                # Try to find and click through to full text
                full_text_link = self.page.query_selector("a.id-link")
                if full_text_link:
                    href = full_text_link.get_attribute("href")
                    logger.info(f"  Found full text link: {href}")
                    full_text_link.click()
                    try:
                        self.page.wait_for_load_state(
                            "networkidle", timeout=self.timeout
                        )
                    except PlaywrightTimeout:
                        logger.info(
                            "  Timeout waiting for publisher page, continuing..."
                        )
                    time.sleep(2)
                    logger.info(f"  Now on: {self.page.url[:70]}...")

            # Method 1: Try clicking download buttons for main PDF
            logger.info(">>> STEP 3: Trying download button click...")
            result = self._click_download_button(self.page, pmid)
            if result:
                downloaded_files.append(str(result))
                logger.info("  Downloaded main PDF via button click")

            # Method 2: Try heuristic PDF link detection (if no PDF yet)
            if not downloaded_files:
                logger.info(">>> STEP 4: Trying heuristic PDF link detection...")
                pdf_links = self._find_pdf_links_heuristic(self.page)
                for i, link in enumerate(pdf_links[:5]):  # Try first 5 links
                    logger.info(
                        f"  Trying link {i + 1}/{min(5, len(pdf_links))}: {link[:60]}..."
                    )
                    result = self._try_download_link(self.page, link, pmid)
                    if result:
                        downloaded_files.append(str(result))
                        logger.info("  Downloaded main PDF via heuristic")
                        break

            # Method 3: Download supplement files (regardless of main PDF success)
            logger.info(">>> STEP 5: Scanning for supplement files...")
            supp_files = self._download_supplements(self.page, pmid)
            for supp in supp_files:
                downloaded_files.append(str(supp))

            # If we got any files, consider it a success
            if downloaded_files:
                logger.info(f">>> SUCCESS: Downloaded {len(downloaded_files)} file(s)")
                return DownloadResult(
                    pmid=pmid,
                    success=True,
                    files_downloaded=downloaded_files,
                    method="auto_with_supplements",
                )

            # Method 4: Use Claude to analyze page (if enabled and still no files)
            if self.use_claude and self.analyzer:
                logger.info(">>> STEP 6: Using Claude to analyze page...")
                page_content = self.page.content()
                logger.info(f"  Page content length: {len(page_content)} chars")
                claude_links = self.analyzer.find_pdf_links(page_content, self.page.url)
                logger.info(f"  Claude found {len(claude_links)} link(s)")

                for link_info in claude_links:
                    confidence = link_info.get("confidence", 0)
                    logger.info(
                        f"  Link: {link_info.get('url', '')[:50]}... (confidence: {confidence})"
                    )
                    if confidence > 0.5:
                        result = self._try_download_link(
                            self.page, link_info["url"], pmid
                        )
                        if result:
                            downloaded_files.append(str(result))
                            return DownloadResult(
                                pmid=pmid,
                                success=True,
                                files_downloaded=downloaded_files,
                                method="claude_assisted",
                            )

            # If we get here, PDF download failed - try HTML scraping as fallback
            logger.info(">>> STEP 7: PDF not found, trying HTML scrape fallback...")

            # Navigate back to article page if we're on reader page
            current = self.page.url
            if "/reader/" in current or "/epub/" in current:
                # Go back to main article page
                article_url = current.replace("/reader/", "/").replace("/epub/", "/")
                logger.info(
                    f"  [HTMLScrape] Navigating to article page: {article_url[:60]}..."
                )
                try:
                    self.page.goto(article_url, wait_until="load", timeout=30000)
                    time.sleep(2)
                    self._handle_cookie_consent(self.page)
                except Exception as e:
                    logger.info(f"  [HTMLScrape] Navigation failed: {e}")

            html_result = self._scrape_article_html(self.page, pmid)
            if html_result:
                return DownloadResult(
                    pmid=pmid,
                    success=True,
                    files_downloaded=[str(html_result)],
                    method="html_scrape",
                )

            # If we get here, everything failed
            logger.info(">>> FAILED: Could not download PDF or scrape HTML")
            return DownloadResult(
                pmid=pmid,
                success=False,
                files_downloaded=[],
                error="Could not find download link or scrape article content",
                method="needs_manual",
                failure_type="paywall",
            )

        except Exception as e:
            logger.error(f">>> ERROR: {e}")
            import traceback

            error_msg = str(e)
            # Detect browser/page crash vs other errors
            if any(
                x in error_msg.lower()
                for x in ["closed", "crashed", "target page", "context"]
            ):
                failure_type = "crash"
            elif any(
                x in error_msg.lower() for x in ["captcha", "cloudflare", "challenge"]
            ):
                failure_type = "captcha_blocked"
            else:
                failure_type = "crash"

            logger.error(traceback.format_exc())
            return DownloadResult(
                pmid=pmid,
                success=False,
                files_downloaded=downloaded_files,
                error=error_msg,
                method="error",
                failure_type=failure_type,
            )


def _interactive_fetch(
    fetcher: BrowserFetcher,
    pmid: str,
    url: str,
    downloads_dir: Path,
    target_dir: Path,
) -> DownloadResult:
    """
    Interactive fetch: navigate to page, let user manually download, monitor files.

    This mode opens the browser to the paper's page, then waits for the user to
    manually download files (login if needed, click download buttons, etc.).
    Once the user confirms, it monitors the Downloads folder and moves files.
    """
    logger.info(f"  Opening: {url}")
    logger.info("  Browser will navigate to the paper. Download files manually.")

    # Take snapshot before navigation
    before_files = {f.name for f in downloads_dir.iterdir() if f.is_file()}

    # Navigate to page
    try:
        fetcher.page.goto(url, wait_until="load", timeout=fetcher.timeout)
    except PlaywrightTimeout:
        logger.warning("  Page load timeout, browser still usable")

    time.sleep(2)

    # Wait for user input
    print("\n" + "-" * 50)
    print(f"  PMID: {pmid}")
    print(f"  URL:  {fetcher.page.url}")
    print("-" * 50)
    print("  Download the PDF/supplements manually in the browser.")
    print("  Press ENTER when done, 's' to skip, 'q' to quit: ", end="", flush=True)

    try:
        response = input().strip().lower()
    except (EOFError, KeyboardInterrupt):
        response = "q"

    if response == "q":
        return DownloadResult(
            pmid=pmid,
            success=False,
            files_downloaded=[],
            error="User quit",
            method="manual",
        )

    if response == "s":
        return DownloadResult(
            pmid=pmid,
            success=False,
            files_downloaded=[],
            error="Skipped",
            method="manual",
        )

    # Small delay for downloads to complete
    time.sleep(1)

    # Check for new files
    after_files = {f.name for f in downloads_dir.iterdir() if f.is_file()}
    new_files = after_files - before_files

    # Filter out partial downloads
    complete_files = [
        f
        for f in new_files
        if not f.endswith((".crdownload", ".part", ".tmp", ".download"))
        and not f.startswith(".")
    ]

    if not complete_files:
        print(
            "  No new files detected. Mark as done anyway? (y/n): ", end="", flush=True
        )
        confirm = input().strip().lower()
        if confirm == "y":
            return DownloadResult(
                pmid=pmid, success=True, files_downloaded=[], method="manual_no_files"
            )
        return DownloadResult(
            pmid=pmid,
            success=False,
            files_downloaded=[],
            error="No files downloaded",
            method="manual",
        )

    # Move and rename files
    downloaded = []
    for idx, filename in enumerate(complete_files):
        src = downloads_dir / filename
        ext = src.suffix.lower()

        if idx == 0 and ext == ".pdf":
            new_name = f"{pmid}_Main_Text.pdf"
        elif ext == ".pdf":
            new_name = f"{pmid}_Supplement_{idx}.pdf"
        elif ext in (".xlsx", ".xls", ".csv"):
            new_name = f"{pmid}_Supp_Data_{idx}{ext}"
        else:
            new_name = f"{pmid}_file_{idx}{ext}"

        dest = target_dir / new_name
        try:
            src.rename(dest)
            logger.info(f"  Moved: {filename} -> {new_name}")
            downloaded.append(str(dest))
        except Exception as e:
            logger.error(f"  Failed to move {filename}: {e}")

    return DownloadResult(
        pmid=pmid,
        success=len(downloaded) > 0,
        files_downloaded=downloaded,
        method="manual",
    )


def run_browser_fetch(
    csv_file: Path,
    downloads_dir: Optional[Path] = None,
    target_dir: Optional[Path] = None,
    headless: bool = False,
    use_claude: bool = False,
    pmid_column: str = "PMID",
    status_column: str = "Status",
    specific_pmids: Optional[List[str]] = None,
    max_papers: Optional[int] = None,
    interactive: bool = False,
    wait_for_captcha: bool = False,
    retry_failures_only: bool = False,
    retry_captcha_only: bool = False,
    retry_paywall_only: bool = False,
    fallback_to_manual: bool = True,
    use_profile: bool = False,
    profile_path: Optional[Path] = None,
    slow_mo: int = 0,
):
    """
    Run browser-based fetching for papers in a CSV file.

    Args:
        csv_file: Path to CSV with paper list
        downloads_dir: Browser downloads directory
        target_dir: Where to save renamed files
        headless: Run browser without visible window
        use_claude: Use Claude to help find download links
        pmid_column: Name of PMID column
        status_column: Name of status column
        specific_pmids: Only process these PMIDs
        max_papers: Maximum number of papers to process
        interactive: Interactive mode - navigate and wait for manual download
        wait_for_captcha: Wait up to 5 minutes for manual CAPTCHA completion
        retry_failures_only: Only process papers with browser_failed status
        fallback_to_manual: Auto-switch to manual download when Cloudflare loops (default True)
    """
    if not PLAYWRIGHT_AVAILABLE:
        logger.error(
            "Playwright not installed. Run: pip install playwright && playwright install chromium"
        )
        sys.exit(1)

    # Set up directories
    downloads_dir = downloads_dir or Path.home() / "Downloads"

    # Auto-detect target directory from CSV location
    if target_dir is None:
        if csv_file.parent.name == "pmc_fulltext":
            target_dir = csv_file.parent
        else:
            target_dir = csv_file.parent / "browser_downloads"

    target_dir.mkdir(parents=True, exist_ok=True)

    # Load CSV
    logger.info(f"Loading CSV: {csv_file}")
    df = pd.read_csv(csv_file)

    if pmid_column not in df.columns:
        logger.error(f"Column '{pmid_column}' not found. Available: {list(df.columns)}")
        sys.exit(1)

    # Add status column if needed
    if status_column not in df.columns:
        df[status_column] = ""

    # Filter to pending papers
    df[status_column] = df[status_column].fillna("")

    if retry_captcha_only:
        # Only retry CAPTCHA-blocked papers
        pending_mask = df[status_column].str.lower() == "captcha_blocked"
        logger.info("Retry mode: processing only captcha_blocked papers")
    elif retry_paywall_only:
        # Only retry paywalled papers
        pending_mask = df[status_column].str.lower() == "paywall"
        logger.info("Retry mode: processing only paywall papers")
    elif retry_failures_only:
        # Retry any failed papers
        failed_statuses = [
            "browser_failed",
            "captcha_blocked",
            "paywall",
            "browser_crashed",
        ]
        pending_mask = df[status_column].str.lower().isin(failed_statuses)
        logger.info(f"Retry mode: processing papers with status in {failed_statuses}")
    else:
        # Process papers that haven't been attempted yet (exclude all terminal states)
        terminal_statuses = [
            "done",
            "skipped",
            "browser_done",
            "browser_failed",
            "captcha_blocked",
            "paywall",
            "browser_crashed",
        ]
        pending_mask = ~df[status_column].str.lower().isin(terminal_statuses)

    # Filter to specific PMIDs if provided
    if specific_pmids:
        pmid_mask = df[pmid_column].astype(str).isin(specific_pmids)
        pending_mask = pending_mask & pmid_mask

    papers_to_process = df[pending_mask].copy()

    if max_papers:
        papers_to_process = papers_to_process.head(max_papers)

    total = len(papers_to_process)
    logger.info(f"Processing {total} papers...")

    if total == 0:
        logger.info("No papers to process")
        return

    # Run browser automation
    # If wait_for_captcha is enabled, use 5 minute timeout for Cloudflare
    captcha_wait_time = 300 if wait_for_captcha else 30
    if wait_for_captcha:
        logger.info(
            "CAPTCHA wait mode: will wait up to 5 minutes for manual CAPTCHA completion"
        )
        if headless:
            logger.warning(
                "WARNING: wait_for_captcha works best with headless=False so you can see the browser"
            )

    results = []
    with BrowserFetcher(
        downloads_dir=downloads_dir,
        target_dir=target_dir,
        headless=headless,
        use_claude=use_claude,
        captcha_wait_time=captcha_wait_time,
        fallback_to_manual=fallback_to_manual,
        use_profile=use_profile,
        profile_path=profile_path,
        slow_mo=slow_mo,
    ) as fetcher:
        for idx, (row_idx, row) in enumerate(papers_to_process.iterrows(), 1):
            pmid = str(row[pmid_column]).strip()

            # Construct URL
            if "URL" in row and pd.notna(row["URL"]):
                url = str(row["URL"]).strip()
            elif "DOI" in row and pd.notna(row["DOI"]):
                doi = str(row["DOI"]).strip()
                if not doi.startswith("http"):
                    url = f"https://doi.org/{doi}"
                else:
                    url = doi
            else:
                url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            logger.info(f"\n[{idx}/{total}] Processing PMID {pmid}")

            if interactive:
                # Interactive mode: navigate, let user download, monitor files
                result = _interactive_fetch(
                    fetcher, pmid, url, downloads_dir, target_dir
                )
            else:
                # Fully automated mode with retry on crash
                max_retries = 3
                for attempt in range(max_retries):
                    try:
                        result = fetcher.fetch_paper(pmid, url)
                        break  # Success, exit retry loop
                    except Exception as e:
                        error_str = str(e).lower()
                        is_crash = any(
                            term in error_str
                            for term in [
                                "crash",
                                "target closed",
                                "connection closed",
                                "browser has been closed",
                                "page closed",
                                "protocol error",
                                "session closed",
                            ]
                        )

                        if is_crash and attempt < max_retries - 1:
                            logger.warning(
                                f"  Browser crash detected (attempt {attempt + 1}/{max_retries}): {e}"
                            )
                            logger.info("  Attempting to restart browser...")
                            try:
                                # Try to restart the browser
                                fetcher._restart_browser()
                                time.sleep(2)  # Give browser time to stabilize
                                logger.info("  Browser restarted, retrying...")
                                continue
                            except Exception as restart_error:
                                logger.error(
                                    f"  Failed to restart browser: {restart_error}"
                                )

                        # Final failure
                        result = DownloadResult(
                            pmid=pmid,
                            success=False,
                            files_downloaded=[],
                            error=f"Browser crash after {attempt + 1} attempts: {str(e)}",
                            method="crash_recovery_failed",
                            failure_type="crash",
                        )
                        break

            results.append(result)

            # Update status in DataFrame with granular failure tracking
            # Add Failure_Reason column if it doesn't exist
            if "Failure_Reason" not in df.columns:
                df["Failure_Reason"] = ""

            if result.success:
                df.at[row_idx, status_column] = "browser_done"
                df.at[row_idx, "Failure_Reason"] = ""
                logger.info(
                    f"  SUCCESS: Downloaded {len(result.files_downloaded)} file(s)"
                )
            else:
                # Use granular status based on failure type
                failure_type = result.failure_type or "unknown"
                if failure_type == "captcha_blocked":
                    df.at[row_idx, status_column] = "captcha_blocked"
                elif failure_type == "paywall":
                    df.at[row_idx, status_column] = "paywall"
                elif failure_type == "crash":
                    df.at[row_idx, status_column] = "browser_crashed"
                else:
                    df.at[row_idx, status_column] = "browser_failed"

                df.at[row_idx, "Failure_Reason"] = result.error or ""
                logger.warning(f"  FAILED [{failure_type}]: {result.error}")

            # Save progress
            df.to_csv(csv_file, index=False)

            # Small delay between papers
            time.sleep(1)

    # Print summary
    successful = sum(1 for r in results if r.success)
    failed = sum(1 for r in results if not r.success)

    print("\n" + "=" * 60)
    print("Browser Fetch Complete")
    print(f"  Successful: {successful}/{total}")
    print(f"  Failed: {failed}/{total}")
    print(f"  Files saved to: {target_dir}")
    print("=" * 60)

    # Save detailed results
    results_file = target_dir / "browser_fetch_results.json"
    with open(results_file, "w") as f:
        json.dump(
            [
                {
                    "pmid": r.pmid,
                    "success": r.success,
                    "files": r.files_downloaded,
                    "error": r.error,
                    "method": r.method,
                    "failure_type": r.failure_type,
                }
                for r in results
            ],
            f,
            indent=2,
        )
    logger.info(f"Results saved to: {results_file}")

    # Print failure breakdown
    if failed > 0:
        captcha_blocked = sum(1 for r in results if r.failure_type == "captcha_blocked")
        paywalled = sum(1 for r in results if r.failure_type == "paywall")
        crashed = sum(1 for r in results if r.failure_type == "crash")
        other = failed - captcha_blocked - paywalled - crashed
        print("\nFailure breakdown:")
        if captcha_blocked:
            print(f"  CAPTCHA blocked: {captcha_blocked}")
        if paywalled:
            print(f"  Paywall/no access: {paywalled}")
        if crashed:
            print(f"  Browser crashed: {crashed}")
        if other:
            print(f"  Other: {other}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Browser-based paper fetcher with optional Claude assistance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (recommended for paywalled papers)
  python browser_fetch.py paywalled_missing.csv --interactive

  # Fully automated mode - tries to find and click download links
  python browser_fetch.py paywalled_missing.csv

  # Use Claude to help find download links on complex sites
  python browser_fetch.py paywalled_missing.csv --use-claude

  # Run headless (no visible browser, automated only)
  python browser_fetch.py paywalled_missing.csv --headless

  # Process specific PMIDs
  python browser_fetch.py paywalled_missing.csv -i --pmids 12345678,87654321

  # Limit to first 10 papers
  python browser_fetch.py paywalled_missing.csv -i --max 10

  # Two-step workflow for sites with CAPTCHAs:
  python browser_fetch.py paywalled_missing.csv              # Step 1: automated run
  python browser_fetch.py paywalled_missing.csv --retry-captcha --wait-for-captcha  # Step 2: manually solve CAPTCHAs

  # Retry only paywalled papers (e.g., after getting VPN/library access)
  python browser_fetch.py paywalled_missing.csv --retry-paywall

Status tracking (in CSV Status column):
  - browser_done: Successfully downloaded
  - captcha_blocked: Blocked by CAPTCHA/Cloudflare
  - paywall: Requires subscription/login
  - browser_crashed: Browser error during fetch
        """,
    )

    parser.add_argument("csv_file", type=Path, help="CSV file with paper list")
    parser.add_argument("--downloads", type=Path, help="Browser downloads directory")
    parser.add_argument("--target-dir", type=Path, help="Target directory for files")
    parser.add_argument(
        "--headless", action="store_true", help="Run browser without window"
    )
    parser.add_argument(
        "--interactive",
        "-i",
        action="store_true",
        help="Interactive mode: navigate to page, let you download manually, auto-organize files",
    )
    parser.add_argument(
        "--use-claude", action="store_true", help="Use Claude to find download links"
    )
    parser.add_argument(
        "--wait-for-captcha",
        "-w",
        action="store_true",
        help="Wait up to 5 minutes for manual CAPTCHA completion (use with visible browser)",
    )
    parser.add_argument(
        "--retry-failures",
        "-r",
        action="store_true",
        help="Retry all previously failed papers (any failure type)",
    )
    parser.add_argument(
        "--retry-captcha",
        action="store_true",
        help="Retry only CAPTCHA-blocked papers (use with --wait-for-captcha for manual solving)",
    )
    parser.add_argument(
        "--retry-paywall",
        action="store_true",
        help="Retry only paywalled papers (useful if you now have institutional access)",
    )
    parser.add_argument("--pmid-column", default="PMID", help="PMID column name")
    parser.add_argument("--status-column", default="Status", help="Status column name")
    parser.add_argument(
        "--pmids", type=str, help="Comma-separated list of specific PMIDs"
    )
    parser.add_argument("--max", type=int, help="Maximum papers to process")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument(
        "--no-fallback-manual",
        action="store_true",
        help="Disable automatic fallback to manual download when Cloudflare loops (default: enabled)",
    )
    parser.add_argument(
        "--use-profile",
        action="store_true",
        help="Use your real Chrome profile (has cookies/history, better against CAPTCHAs)",
    )
    parser.add_argument(
        "--profile-path",
        type=Path,
        help="Custom Chrome profile path (default: ~/Library/Application Support/Google/Chrome)",
    )
    parser.add_argument(
        "--slow-mo",
        type=int,
        default=0,
        help="Slow down browser actions by N milliseconds (helps avoid detection)",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not args.csv_file.exists():
        logger.error(f"CSV file not found: {args.csv_file}")
        sys.exit(1)

    specific_pmids = None
    if args.pmids:
        specific_pmids = [p.strip() for p in args.pmids.split(",")]

    run_browser_fetch(
        csv_file=args.csv_file,
        downloads_dir=args.downloads,
        target_dir=args.target_dir,
        headless=args.headless,
        use_claude=args.use_claude,
        pmid_column=args.pmid_column,
        status_column=args.status_column,
        specific_pmids=specific_pmids,
        max_papers=args.max,
        interactive=args.interactive,
        wait_for_captcha=args.wait_for_captcha,
        retry_failures_only=args.retry_failures,
        retry_captcha_only=args.retry_captcha,
        retry_paywall_only=args.retry_paywall,
        fallback_to_manual=not args.no_fallback_manual,
        use_profile=args.use_profile,
        profile_path=args.profile_path,
        slow_mo=args.slow_mo,
    )


if __name__ == "__main__":
    main()
