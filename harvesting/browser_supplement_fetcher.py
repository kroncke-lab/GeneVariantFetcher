"""
Browser-based Supplement Fetcher

Uses Playwright for publishers with Cloudflare/bot protection.
Fallback when direct HTTP requests fail with 403.
"""

import asyncio
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import urljoin, urlparse

try:
    from playwright.async_api import Browser, Page, async_playwright

    PLAYWRIGHT_AVAILABLE = True
except ImportError:
    PLAYWRIGHT_AVAILABLE = False


async def fetch_karger_supplements_browser(
    doi: str, output_dir: Path, timeout_ms: int = 30000
) -> Tuple[List[Dict], Optional[str]]:
    """
    Fetch Karger supplements using browser automation.

    Args:
        doi: DOI of the article (e.g., "10.1159/000440608")
        output_dir: Directory to save downloaded supplements
        timeout_ms: Timeout in milliseconds

    Returns:
        Tuple of (list of downloaded files, error message or None)
    """
    if not PLAYWRIGHT_AVAILABLE:
        return [], "Playwright not installed"

    downloaded_files = []
    article_url = f"https://doi.org/{doi}"

    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        context = await browser.new_context(
            user_agent="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
        )
        page = await context.new_page()

        try:
            # Navigate to article page
            print(f"    🌐 Loading article page via browser: {article_url}")
            await page.goto(article_url, wait_until="networkidle", timeout=timeout_ms)

            # Wait for Cloudflare challenge if present
            await asyncio.sleep(3)

            final_url = page.url
            print(f"    → Resolved to: {final_url}")

            # Look for supplementary material section
            supp_links = []

            # Try common supplement link patterns
            selectors = [
                'a[href*="suppl"]',
                'a:has-text("Supplementary")',
                'a:has-text("supplement")',
                ".supplementary-material a",
                "#supplementary-information a",
            ]

            for selector in selectors:
                try:
                    links = await page.query_selector_all(selector)
                    for link in links:
                        href = await link.get_attribute("href")
                        text = await link.inner_text()
                        if href:
                            full_url = urljoin(final_url, href)
                            supp_links.append(
                                {
                                    "url": full_url,
                                    "name": text or Path(urlparse(full_url).path).name,
                                }
                            )
                except Exception:
                    continue

            # Also try constructing the supplement URL directly
            karger_supp_url = f"https://karger.com/doi/suppl/{doi}"
            print(f"    → Trying supplement URL: {karger_supp_url}")

            try:
                await page.goto(
                    karger_supp_url, wait_until="networkidle", timeout=timeout_ms
                )
                await asyncio.sleep(2)

                # If we got a page (not 404), look for download links
                if "404" not in await page.title():
                    file_links = await page.query_selector_all(
                        'a[href*=".pdf"], a[href*=".xls"], a[href*=".doc"], a[href*=".csv"]'
                    )
                    for link in file_links:
                        href = await link.get_attribute("href")
                        text = await link.inner_text()
                        if href:
                            full_url = urljoin(karger_supp_url, href)
                            if not any(s["url"] == full_url for s in supp_links):
                                supp_links.append(
                                    {
                                        "url": full_url,
                                        "name": text
                                        or Path(urlparse(full_url).path).name,
                                    }
                                )
            except Exception as e:
                print(f"    ⚠️ Supplement URL failed: {e}")

            # Download found supplements
            for supp in supp_links:
                try:
                    url = supp["url"]
                    filename = supp["name"] or "supplement"

                    # Ensure filename has extension
                    if not any(
                        filename.endswith(ext)
                        for ext in [
                            ".pdf",
                            ".xls",
                            ".xlsx",
                            ".doc",
                            ".docx",
                            ".csv",
                            ".zip",
                        ]
                    ):
                        # Try to get extension from URL
                        ext_match = re.search(
                            r"\.(pdf|xlsx?|docx?|csv|zip)(?:\?|$)", url, re.I
                        )
                        if ext_match:
                            filename = f"{filename}.{ext_match.group(1)}"

                    output_path = output_dir / filename

                    # Download using page context (handles cookies/auth)
                    async with page.expect_download() as download_info:
                        await page.goto(url)
                    download = await download_info.value
                    await download.save_as(output_path)

                    downloaded_files.append(
                        {"url": url, "name": filename, "path": str(output_path)}
                    )
                    print(f"    ✅ Downloaded: {filename}")

                except Exception as e:
                    print(
                        f"    ⚠️ Failed to download {supp.get('name', 'unknown')}: {e}"
                    )

            return downloaded_files, None

        except Exception as e:
            return downloaded_files, str(e)
        finally:
            await browser.close()


def fetch_karger_supplements_sync(
    doi: str, output_dir: Path, timeout_ms: int = 30000
) -> Tuple[List[Dict], Optional[str]]:
    """
    Synchronous wrapper for fetch_karger_supplements_browser.
    """
    return asyncio.run(fetch_karger_supplements_browser(doi, output_dir, timeout_ms))


if __name__ == "__main__":
    # Test with PMID 26496715
    import tempfile

    test_doi = "10.1159/000440608"
    output_dir = Path(tempfile.mkdtemp())

    print(f"Testing Karger supplement fetch for DOI: {test_doi}")
    print(f"Output directory: {output_dir}")

    files, error = fetch_karger_supplements_sync(test_doi, output_dir)

    if error:
        print(f"Error: {error}")
    else:
        print(f"Downloaded {len(files)} files:")
        for f in files:
            print(f"  - {f['name']}")
