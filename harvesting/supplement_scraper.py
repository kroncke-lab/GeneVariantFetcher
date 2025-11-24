"""
Supplement Scraper Module

Web scraping logic for extracting supplemental files from publisher websites.
Includes domain-specific scrapers for Nature, Elsevier, and a generic fallback.
"""

import re
import json
from typing import List, Dict
from pathlib import Path
from urllib.parse import urlparse, urljoin
from bs4 import BeautifulSoup


class SupplementScraper:
    """Scrapes supplemental files from various publisher websites."""

    def scrape_nature_supplements(self, html: str, base_url: str) -> List[Dict]:
        """
        Scrape supplemental files from a Nature journal page.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL for resolving relative links

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys
        """
        print("  Scraping with scrape_nature_supplements...")
        soup = BeautifulSoup(html, 'html.parser')
        found_files = []

        # Nature articles often have a "Supplementary Information" section
        supp_info_section = soup.find(id='supplementary-information')
        if not supp_info_section:
            # Fallback to searching for common heading patterns
            supp_info_section = soup.find('h2', string=re.compile(r'Supplementary (Information|Data)', re.IGNORECASE))
            if supp_info_section:
                supp_info_section = supp_info_section.parent

        if supp_info_section:
            print("    Found 'Supplementary Information' section.")
            for a in supp_info_section.find_all('a', href=True):
                href = a['href']
                # Filter for links that look like file downloads
                if '/articles/' in href and '/figures/' not in href and 'author-information' not in href:
                    url = urljoin(base_url, href)
                    filename = Path(urlparse(url).path).name
                    if filename and not any(f['name'] == filename for f in found_files):
                        found_files.append({'url': url, 'name': filename})
                        print(f"      Found potential supplement: {filename}")

        if not found_files:
            print("    'Supplementary Information' section not found or empty. Trying generic scan.")
            return self.scrape_generic_supplements(html, base_url)

        return found_files

    def scrape_elsevier_supplements(self, html: str, base_url: str) -> List[Dict]:
        """
        Scrape supplemental files from an Elsevier/GIM journal page.
        Uses a multi-step approach: JSON data, regex patterns, and HTML parsing.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL for resolving relative links

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys
        """
        print("  Scraping with scrape_elsevier_supplements...")
        soup = BeautifulSoup(html, 'html.parser')
        found_files = []

        # 1. Look for data in <script type="application/json"> blocks
        # This is often more reliable than parsing the HTML structure
        for script in soup.find_all('script', type='application/json'):
            try:
                data = json.loads(script.string)
                # This path can be very specific and fragile; may need adjustment
                supp_data = data.get('article', {}).get('supplementaryMaterials', {}).get('supplementaryMaterial', [])
                for item in supp_data:
                    url = item.get('downloadUrl')
                    filename = item.get('title')
                    if url and filename and not any(f['name'] == filename for f in found_files):
                        found_files.append({'url': url, 'name': filename})
                        print(f"    Found supplement in JSON data: {filename}")
                if found_files:
                    return found_files
            except (json.JSONDecodeError, AttributeError):
                continue

        # 2. Regex for "mmc" (multimedia component) links, a common Elsevier pattern
        mmc_links = re.findall(r'href="(/cms/attachment/[^"]+/mmc\d+\.(?:pdf|docx|xlsx|zip|csv|txt))"', html)
        for link in mmc_links:
            url = urljoin(base_url, link)
            filename = Path(urlparse(url).path).name
            if filename and not any(f['name'] == filename for f in found_files):
                found_files.append({'url': url, 'name': filename})
                print(f"    Found MMC link via regex: {filename}")
        if found_files:
            return found_files

        # 3. Look for specific "Download file" links in ScienceDirect as a fallback
        for a in soup.find_all('a', class_='S_C_9cf8451f', href=True):
            if 'Download file' in a.get_text():
                url = urljoin(base_url, a['href'])
                filename = a.get_text().replace('Download file', '').strip()
                if filename and not any(f['name'] == filename for f in found_files):
                    found_files.append({'url': url, 'name': filename})
                    print(f"    Found ScienceDirect download link: {filename}")
        if found_files:
            return found_files

        print("    No specific Elsevier/GIM supplements found. Trying generic scan.")
        return self.scrape_generic_supplements(html, base_url)

    def scrape_generic_supplements(self, html: str, base_url: str) -> List[Dict]:
        """
        A best-effort generic scraper for supplemental files.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL for resolving relative links

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys
        """
        print("  Scraping with scrape_generic_supplements...")
        soup = BeautifulSoup(html, 'html.parser')
        found_files = []

        keywords = ['supplement', 'supporting', 'appendix', 'additional file']
        file_extensions = ['.pdf', '.docx', '.xlsx', '.csv', '.zip', '.rar', '.gz', '.txt', '.doc']

        for a in soup.find_all('a', href=True):
            link_text = a.get_text().lower()
            href = a['href'].lower()

            # Check if link text or href contain any of the keywords or file extensions
            is_supplement_link = any(keyword in link_text for keyword in keywords)
            is_file_link = any(href.endswith(ext) for ext in file_extensions)

            if is_supplement_link or is_file_link:
                try:
                    url = urljoin(base_url, a['href'])
                    filename = Path(urlparse(url).path).name

                    # Basic filtering to avoid irrelevant links
                    if filename and not any(f['name'] == filename for f in found_files):
                        found_files.append({'url': url, 'name': filename})
                        print(f"    Found potential supplement: {filename}")
                except Exception:
                    continue  # Ignore malformed URLs

        return found_files
