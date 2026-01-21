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

import os
import sys
import time
import json
import logging
import re
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass
from datetime import datetime

import pandas as pd

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
        sync_playwright,
        Page,
        Browser,
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
    ):
        self.downloads_dir = Path(downloads_dir)
        self.target_dir = Path(target_dir)
        self.headless = headless
        self.timeout = timeout
        self.use_claude = use_claude
        self.captcha_wait_time = captcha_wait_time

        self.analyzer = PageAnalyzer() if use_claude else None
        self.browser: Optional[Browser] = None
        self.page: Optional[Page] = None

    def __enter__(self):
        if not PLAYWRIGHT_AVAILABLE:
            raise RuntimeError("Playwright not installed")

        self.playwright = sync_playwright().start()
        self.browser = self.playwright.chromium.launch(
            headless=self.headless,
            downloads_path=str(self.downloads_dir),
        )
        self.context = self.browser.new_context(
            accept_downloads=True,
            user_agent=(
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/120.0.0.0 Safari/537.36"
            ),
        )
        self.page = self.context.new_page()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.browser:
            self.browser.close()
        if hasattr(self, "playwright"):
            self.playwright.stop()

    def _get_existing_files(self) -> set:
        """Get set of existing files in downloads directory."""
        return {f.name for f in self.downloads_dir.iterdir() if f.is_file()}

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

            logger.info(f"  [EmbeddedPDF] Response is not a valid PDF")

        except Exception as e:
            logger.info(f"  [EmbeddedPDF] ERROR: {e}")

        return None

    def _wait_for_cloudflare(self, page: Page, max_wait: Optional[int] = None) -> bool:
        """
        Detect and wait for Cloudflare challenge to be completed.

        Returns True if page is ready (no Cloudflare or challenge passed).
        Returns False if still blocked after max_wait seconds.
        """
        if max_wait is None:
            max_wait = self.captcha_wait_time
        start = time.time()
        check_count = 0
        clear_count = 0  # Track consecutive clear checks for stability
        initial_url = page.url
        cloudflare_detected = False
        challenge_solve_count = 0  # Track how many times we've "solved" a challenge

        while time.time() - start < max_wait:
            check_count += 1
            # Check for common Cloudflare indicators
            try:
                page_content = page.content().lower()
                title = page.title().lower() if page.title() else ""
                current_url = page.url
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
            ]

            # Check for actual content indicators (means we're past Cloudflare)
            content_indicators = [
                "<article" in page_content,
                "<main" in page_content,
                "abstract" in page_content and len(page_content) > 5000,
                "pdf" in page_content and "download" in page_content,
                "full text" in page_content,
            ]

            # URL change after Cloudflare detection often means success
            url_changed = current_url != initial_url and cloudflare_detected

            if any(cloudflare_indicators):
                # If we were previously clear and now see challenge again, we're in a loop
                if clear_count > 0:
                    challenge_solve_count += 1
                    if challenge_solve_count >= 3:
                        elapsed = int(time.time() - start)
                        logger.warning(
                            f"  [Cloudflare] Infinite challenge loop detected "
                            f"({challenge_solve_count} cycles, {elapsed}s elapsed). "
                            f"Site may require manual intervention."
                        )
                        return False
                    logger.info(
                        f"  [Cloudflare] New challenge after clear (cycle {challenge_solve_count})"
                    )

                cloudflare_detected = True
                clear_count = 0  # Reset stability counter
                elapsed = int(time.time() - start)
                logger.info(
                    f"  [Cloudflare] Challenge detected, waiting... ({elapsed}s/{max_wait}s)"
                )
                time.sleep(2)
                continue

            # No Cloudflare indicators detected
            clear_count += 1

            # If we detected Cloudflare earlier and now see content, we passed
            if cloudflare_detected and any(content_indicators):
                # Wait a bit longer and re-verify to catch reload loops
                logger.info(f"  [Cloudflare] Content detected, verifying stability...")
                time.sleep(3)

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
                    ]
                    if any(recheck_indicators):
                        logger.info(f"  [Cloudflare] Challenge reappeared after reload")
                        clear_count = 0
                        continue
                except Exception:
                    pass

                logger.info(f"  [Cloudflare] Challenge passed! Content loaded.")
                return True

            # Require multiple consecutive clear checks for stability
            # This prevents false positives during page transitions
            if clear_count >= 5:  # Increased from 3 to 5 for more stability
                if cloudflare_detected:
                    # Extra verification: wait and re-check one more time
                    logger.info(f"  [Cloudflare] Appears clear, final verification...")
                    time.sleep(3)
                    try:
                        final_content = page.content().lower()
                        final_title = page.title().lower() if page.title() else ""
                        final_indicators = [
                            "checking your browser" in final_content,
                            "just a moment" in final_title,
                            "cf-browser-verification" in final_content,
                            "challenge-running" in final_content,
                            "turnstile" in final_content,
                        ]
                        if any(final_indicators):
                            logger.info(
                                f"  [Cloudflare] Challenge reappeared on final check"
                            )
                            clear_count = 0
                            continue
                    except Exception:
                        pass

                if check_count > 1:
                    logger.info(f"  [Cloudflare] Passed or not present")
                return True

            # Still checking for stability
            time.sleep(1)

        logger.info(
            f"  [Cloudflare] Still blocked after {max_wait}s, continuing anyway"
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
            logger.info(f"  [Heuristic]   {i+1}. {link[:70]}...")

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
                    f"  [Button] Timeout after clicking '{selector}', checking files..."
                )
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
        logger.info(f"\n{'='*60}")
        logger.info(f"Fetching PMID {pmid}")
        logger.info(f"URL: {url}")
        logger.info(f"{'='*60}")
        downloaded_files = []

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

            # Check for Cloudflare and wait if needed
            if not self._wait_for_cloudflare(self.page):
                logger.warning(
                    "  Cloudflare challenge not passed, may have limited success"
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
                logger.info(f"  Downloaded main PDF via button click")

            # Method 2: Try heuristic PDF link detection (if no PDF yet)
            if not downloaded_files:
                logger.info(">>> STEP 4: Trying heuristic PDF link detection...")
                pdf_links = self._find_pdf_links_heuristic(self.page)
                for i, link in enumerate(pdf_links[:5]):  # Try first 5 links
                    logger.info(
                        f"  Trying link {i+1}/{min(5, len(pdf_links))}: {link[:60]}..."
                    )
                    result = self._try_download_link(self.page, link, pmid)
                    if result:
                        downloaded_files.append(str(result))
                        logger.info(f"  Downloaded main PDF via heuristic")
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

            # If we get here, automatic download failed
            logger.info(">>> FAILED: Could not find download link automatically")
            return DownloadResult(
                pmid=pmid,
                success=False,
                files_downloaded=[],
                error="Could not find download link automatically",
                method="needs_manual",
            )

        except Exception as e:
            logger.error(f">>> ERROR: {e}")
            import traceback

            logger.error(traceback.format_exc())
            return DownloadResult(
                pmid=pmid,
                success=False,
                files_downloaded=downloaded_files,
                error=str(e),
                method="error",
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

    if retry_failures_only:
        # Only retry papers that previously failed
        pending_mask = df[status_column].str.lower() == "browser_failed"
        logger.info("Retry mode: processing only browser_failed papers")
    else:
        # Process papers that haven't been attempted yet
        pending_mask = ~df[status_column].str.lower().isin(
            ["done", "skipped", "browser_done", "browser_failed"]
        )

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
                # Fully automated mode
                result = fetcher.fetch_paper(pmid, url)

            results.append(result)

            # Update status in DataFrame
            if result.success:
                df.at[row_idx, status_column] = "browser_done"
                logger.info(
                    f"  SUCCESS: Downloaded {len(result.files_downloaded)} file(s)"
                )
            else:
                df.at[row_idx, status_column] = "browser_failed"
                logger.warning(f"  FAILED: {result.error}")

            # Save progress
            df.to_csv(csv_file, index=False)

            # Small delay between papers
            time.sleep(1)

    # Print summary
    successful = sum(1 for r in results if r.success)
    failed = sum(1 for r in results if not r.success)

    print("\n" + "=" * 60)
    print(f"Browser Fetch Complete")
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
                }
                for r in results
            ],
            f,
            indent=2,
        )
    logger.info(f"Results saved to: {results_file}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Browser-based paper fetcher with optional Claude assistance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (recommended for paywalled papers)
  # Opens browser, you download manually, files auto-organized
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

  # Two-step workflow: automated first, then retry failures with CAPTCHA wait
  python browser_fetch.py paywalled_missing.csv              # Step 1: automated
  python browser_fetch.py paywalled_missing.csv --retry-failures --wait-for-captcha  # Step 2: retry failures
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
        help="Only process papers with browser_failed status (retry previously failed papers)",
    )
    parser.add_argument("--pmid-column", default="PMID", help="PMID column name")
    parser.add_argument("--status-column", default="Status", help="Status column name")
    parser.add_argument(
        "--pmids", type=str, help="Comma-separated list of specific PMIDs"
    )
    parser.add_argument("--max", type=int, help="Maximum papers to process")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")

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
    )


if __name__ == "__main__":
    main()
