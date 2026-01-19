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

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

# Try to import Playwright
try:
    from playwright.sync_api import sync_playwright, Page, Browser, TimeoutError as PlaywrightTimeout
    PLAYWRIGHT_AVAILABLE = True
except ImportError:
    PLAYWRIGHT_AVAILABLE = False
    logger.warning("Playwright not installed. Run: pip install playwright && playwright install chromium")

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
                messages=[{"role": "user", "content": prompt}]
            )

            # Parse JSON from response
            text = response.content[0].text.strip()
            # Find JSON array in response
            match = re.search(r'\[[\s\S]*\]', text)
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
        'sciencedirect.com': [
            'a.pdf-download-btn-link',
            'a[data-aa-name="PDF link"]',
            '#pdfLink',
        ],
        'nature.com': [
            'a[data-article-pdf]',
            'a.c-pdf-download__link',
        ],
        'springer.com': [
            'a.c-pdf-download__link',
            'a[data-track-action="Download PDF"]',
        ],
        'wiley.com': [
            'a.pdf-download',
            'a[title="PDF"]',
        ],
        'cell.com': [
            'a.pdf-download',
            'a[data-action="download-pdf"]',
        ],
        'pnas.org': [
            'a.article-dl-pdf-link',
        ],
        'ahajournals.org': [
            'a.article__ctrl--pdf',
        ],
        'oup.com': [  # Oxford Academic
            'a.article-pdf-download',
            'a[data-track="pdf-download"]',
        ],
        'pubmed.ncbi.nlm.nih.gov': [
            'a.link-item.pmc',  # Link to PMC
            'a.id-link',  # DOI link
        ],
    }

    def __init__(
        self,
        downloads_dir: Path,
        target_dir: Path,
        headless: bool = False,
        use_claude: bool = False,
        timeout: int = 60000,
    ):
        self.downloads_dir = Path(downloads_dir)
        self.target_dir = Path(target_dir)
        self.headless = headless
        self.timeout = timeout
        self.use_claude = use_claude

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
        if hasattr(self, 'playwright'):
            self.playwright.stop()

    def _get_existing_files(self) -> set:
        """Get set of existing files in downloads directory."""
        return {f.name for f in self.downloads_dir.iterdir() if f.is_file()}

    def _wait_for_new_file(self, before_files: set, timeout: int = 60) -> Optional[Path]:
        """Wait for a new file to appear in downloads directory."""
        start = time.time()
        while time.time() - start < timeout:
            current_files = {f.name for f in self.downloads_dir.iterdir() if f.is_file()}
            new_files = current_files - before_files

            # Filter out partial downloads
            complete_files = [
                f for f in new_files
                if not f.endswith(('.crdownload', '.part', '.tmp'))
            ]

            if complete_files:
                # Return the newest file
                newest = max(complete_files, key=lambda f: (self.downloads_dir / f).stat().st_mtime)
                return self.downloads_dir / newest

            time.sleep(0.5)

        return None

    def _find_pdf_links_heuristic(self, page: Page) -> List[str]:
        """Find PDF links using heuristics (no Claude)."""
        links = []
        content = page.content()

        # Try regex patterns
        for pattern in self.PDF_PATTERNS:
            matches = re.findall(pattern, content, re.IGNORECASE)
            links.extend(matches)

        # Try publisher-specific selectors
        current_url = page.url
        for domain, selectors in self.PUBLISHER_SELECTORS.items():
            if domain in current_url:
                for selector in selectors:
                    try:
                        elements = page.query_selector_all(selector)
                        for el in elements:
                            href = el.get_attribute('href')
                            if href:
                                links.append(href)
                    except Exception:
                        pass

        # Make URLs absolute
        base_url = '/'.join(current_url.split('/')[:3])
        absolute_links = []
        for link in links:
            if link.startswith('http'):
                absolute_links.append(link)
            elif link.startswith('//'):
                absolute_links.append('https:' + link)
            elif link.startswith('/'):
                absolute_links.append(base_url + link)
            else:
                absolute_links.append(current_url.rsplit('/', 1)[0] + '/' + link)

        return list(set(absolute_links))

    def _try_download_link(self, page: Page, url: str, pmid: str) -> Optional[Path]:
        """Try to download a file from a URL."""
        before_files = self._get_existing_files()

        try:
            # Set up download handler
            with page.expect_download(timeout=self.timeout) as download_info:
                page.goto(url, wait_until='commit', timeout=self.timeout)

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

        # Common download button selectors
        button_selectors = [
            'a:has-text("PDF")',
            'a:has-text("Download PDF")',
            'a:has-text("Full Text PDF")',
            'button:has-text("PDF")',
            'a[download]',
            '.pdf-download',
            '.download-pdf',
        ]

        for selector in button_selectors:
            try:
                button = page.query_selector(selector)
                if button and button.is_visible():
                    logger.info(f"  Clicking download button: {selector}")

                    with page.expect_download(timeout=self.timeout) as download_info:
                        button.click()

                    download = download_info.value
                    dest_path = self.target_dir / f"{pmid}_{download.suggested_filename}"
                    download.save_as(str(dest_path))
                    return dest_path

            except PlaywrightTimeout:
                new_file = self._wait_for_new_file(before_files, timeout=5)
                if new_file:
                    dest_path = self.target_dir / f"{pmid}_{new_file.name}"
                    new_file.rename(dest_path)
                    return dest_path
            except Exception:
                continue

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
        logger.info(f"Fetching PMID {pmid}: {url}")
        downloaded_files = []

        try:
            # Navigate to the paper - use 'load' instead of 'networkidle' for faster response
            try:
                self.page.goto(url, wait_until='load', timeout=self.timeout)
            except PlaywrightTimeout:
                # Page might still be usable even if not fully loaded
                logger.warning("  Page load timeout, continuing anyway...")

            time.sleep(3)  # Wait for dynamic content

            # Handle common redirects (PubMed -> publisher)
            current_url = self.page.url
            if 'pubmed.ncbi.nlm.nih.gov' in current_url:
                # Try to find and click through to full text
                full_text_link = self.page.query_selector('a.id-link')
                if full_text_link:
                    full_text_link.click()
                    self.page.wait_for_load_state('networkidle', timeout=self.timeout)
                    time.sleep(2)

            # Method 1: Try clicking download buttons
            result = self._click_download_button(self.page, pmid)
            if result:
                downloaded_files.append(str(result))
                return DownloadResult(
                    pmid=pmid,
                    success=True,
                    files_downloaded=downloaded_files,
                    method="button_click"
                )

            # Method 2: Try heuristic PDF link detection
            pdf_links = self._find_pdf_links_heuristic(self.page)
            for link in pdf_links[:5]:  # Try first 5 links
                result = self._try_download_link(self.page, link, pmid)
                if result:
                    downloaded_files.append(str(result))
                    return DownloadResult(
                        pmid=pmid,
                        success=True,
                        files_downloaded=downloaded_files,
                        method="heuristic"
                    )

            # Method 3: Use Claude to analyze page (if enabled)
            if self.use_claude and self.analyzer:
                logger.info("  Using Claude to analyze page...")
                page_content = self.page.content()
                claude_links = self.analyzer.find_pdf_links(page_content, self.page.url)

                for link_info in claude_links:
                    if link_info.get('confidence', 0) > 0.5:
                        result = self._try_download_link(
                            self.page,
                            link_info['url'],
                            pmid
                        )
                        if result:
                            downloaded_files.append(str(result))
                            return DownloadResult(
                                pmid=pmid,
                                success=True,
                                files_downloaded=downloaded_files,
                                method="claude_assisted"
                            )

            # If we get here, automatic download failed
            return DownloadResult(
                pmid=pmid,
                success=False,
                files_downloaded=[],
                error="Could not find download link automatically",
                method="needs_manual"
            )

        except Exception as e:
            logger.error(f"  Error fetching {pmid}: {e}")
            return DownloadResult(
                pmid=pmid,
                success=False,
                files_downloaded=downloaded_files,
                error=str(e),
                method="error"
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
        fetcher.page.goto(url, wait_until='load', timeout=fetcher.timeout)
    except PlaywrightTimeout:
        logger.warning("  Page load timeout, browser still usable")

    time.sleep(2)

    # Wait for user input
    print("\n" + "-" * 50)
    print(f"  PMID: {pmid}")
    print(f"  URL:  {fetcher.page.url}")
    print("-" * 50)
    print("  Download the PDF/supplements manually in the browser.")
    print("  Press ENTER when done, 's' to skip, 'q' to quit: ", end='', flush=True)

    try:
        response = input().strip().lower()
    except (EOFError, KeyboardInterrupt):
        response = 'q'

    if response == 'q':
        return DownloadResult(pmid=pmid, success=False, files_downloaded=[], error="User quit", method="manual")

    if response == 's':
        return DownloadResult(pmid=pmid, success=False, files_downloaded=[], error="Skipped", method="manual")

    # Small delay for downloads to complete
    time.sleep(1)

    # Check for new files
    after_files = {f.name for f in downloads_dir.iterdir() if f.is_file()}
    new_files = after_files - before_files

    # Filter out partial downloads
    complete_files = [
        f for f in new_files
        if not f.endswith(('.crdownload', '.part', '.tmp', '.download'))
        and not f.startswith('.')
    ]

    if not complete_files:
        print("  No new files detected. Mark as done anyway? (y/n): ", end='', flush=True)
        confirm = input().strip().lower()
        if confirm == 'y':
            return DownloadResult(pmid=pmid, success=True, files_downloaded=[], method="manual_no_files")
        return DownloadResult(pmid=pmid, success=False, files_downloaded=[], error="No files downloaded", method="manual")

    # Move and rename files
    downloaded = []
    for idx, filename in enumerate(complete_files):
        src = downloads_dir / filename
        ext = src.suffix.lower()

        if idx == 0 and ext == '.pdf':
            new_name = f"{pmid}_Main_Text.pdf"
        elif ext == '.pdf':
            new_name = f"{pmid}_Supplement_{idx}.pdf"
        elif ext in ('.xlsx', '.xls', '.csv'):
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
        method="manual"
    )


def run_browser_fetch(
    csv_file: Path,
    downloads_dir: Optional[Path] = None,
    target_dir: Optional[Path] = None,
    headless: bool = False,
    use_claude: bool = False,
    pmid_column: str = 'PMID',
    status_column: str = 'Status',
    specific_pmids: Optional[List[str]] = None,
    max_papers: Optional[int] = None,
    interactive: bool = False,
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
    """
    if not PLAYWRIGHT_AVAILABLE:
        logger.error("Playwright not installed. Run: pip install playwright && playwright install chromium")
        sys.exit(1)

    # Set up directories
    downloads_dir = downloads_dir or Path.home() / "Downloads"

    # Auto-detect target directory from CSV location
    if target_dir is None:
        if csv_file.parent.name == 'pmc_fulltext':
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
        df[status_column] = ''

    # Filter to pending papers
    df[status_column] = df[status_column].fillna('')
    pending_mask = ~df[status_column].str.lower().isin(['done', 'skipped', 'browser_done', 'browser_failed'])

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
    results = []
    with BrowserFetcher(
        downloads_dir=downloads_dir,
        target_dir=target_dir,
        headless=headless,
        use_claude=use_claude,
    ) as fetcher:
        for idx, (row_idx, row) in enumerate(papers_to_process.iterrows(), 1):
            pmid = str(row[pmid_column]).strip()

            # Construct URL
            if 'URL' in row and pd.notna(row['URL']):
                url = str(row['URL']).strip()
            elif 'DOI' in row and pd.notna(row['DOI']):
                doi = str(row['DOI']).strip()
                if not doi.startswith('http'):
                    url = f"https://doi.org/{doi}"
                else:
                    url = doi
            else:
                url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            logger.info(f"\n[{idx}/{total}] Processing PMID {pmid}")

            if interactive:
                # Interactive mode: navigate, let user download, monitor files
                result = _interactive_fetch(fetcher, pmid, url, downloads_dir, target_dir)
            else:
                # Fully automated mode
                result = fetcher.fetch_paper(pmid, url)

            results.append(result)

            # Update status in DataFrame
            if result.success:
                df.at[row_idx, status_column] = 'browser_done'
                logger.info(f"  SUCCESS: Downloaded {len(result.files_downloaded)} file(s)")
            else:
                df.at[row_idx, status_column] = 'browser_failed'
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
    with open(results_file, 'w') as f:
        json.dump(
            [
                {
                    'pmid': r.pmid,
                    'success': r.success,
                    'files': r.files_downloaded,
                    'error': r.error,
                    'method': r.method,
                }
                for r in results
            ],
            f,
            indent=2
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
        """
    )

    parser.add_argument('csv_file', type=Path, help='CSV file with paper list')
    parser.add_argument('--downloads', type=Path, help='Browser downloads directory')
    parser.add_argument('--target-dir', type=Path, help='Target directory for files')
    parser.add_argument('--headless', action='store_true', help='Run browser without window')
    parser.add_argument('--interactive', '-i', action='store_true',
                        help='Interactive mode: navigate to page, let you download manually, auto-organize files')
    parser.add_argument('--use-claude', action='store_true', help='Use Claude to find download links')
    parser.add_argument('--pmid-column', default='PMID', help='PMID column name')
    parser.add_argument('--status-column', default='Status', help='Status column name')
    parser.add_argument('--pmids', type=str, help='Comma-separated list of specific PMIDs')
    parser.add_argument('--max', type=int, help='Maximum papers to process')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not args.csv_file.exists():
        logger.error(f"CSV file not found: {args.csv_file}")
        sys.exit(1)

    specific_pmids = None
    if args.pmids:
        specific_pmids = [p.strip() for p in args.pmids.split(',')]

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
    )


if __name__ == "__main__":
    main()
