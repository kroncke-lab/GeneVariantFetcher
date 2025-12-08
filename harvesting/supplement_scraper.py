"""
Supplement Scraper Module

Web scraping logic for extracting supplemental files and full-text content
from publisher websites.
Includes domain-specific scrapers for Nature, Elsevier, and a generic fallback.
Also includes full-text extraction for free articles without PMCIDs.
"""

import re
import json
from typing import List, Dict, Optional, Tuple
from pathlib import Path
from urllib.parse import urlparse, urljoin
from bs4 import BeautifulSoup


class SupplementScraper:
    """Scrapes supplemental files from various publisher websites."""

    def _normalize_pmc_url(self, url: str, base_url: str) -> str:
        """
        Normalize PMC supplement URLs to use the correct path format.

        NCBI PMC pages sometimes have links with /articles/instance/{id}/bin/
        but the actual files may be served from different URL patterns.
        This method returns a normalized URL, but the caller should be prepared
        to try multiple URL formats if the first one fails.

        Args:
            url: The URL to normalize
            base_url: The base PMC article URL

        Returns:
            Normalized URL
        """
        # Check if this is a PMC page and if the URL needs fixing
        if 'ncbi.nlm.nih.gov/pmc/articles/' in base_url or 'pmc.ncbi.nlm.nih.gov/articles/' in base_url:
            # Extract the PMCID from the base URL
            import re
            pmcid_match = re.search(r'/(?:pmc/)?articles/(PMC\d+)', base_url)
            if pmcid_match:
                pmcid = pmcid_match.group(1)
                numeric_id = pmcid.replace('PMC', '')

                # Fix URLs that use /articles/instance/{numeric_id}/bin/ format
                instance_match = re.search(r'/articles/instance/(\d+)/bin/(.+)$', url)
                if instance_match:
                    filename = instance_match.group(2)
                    # Use the new pmc.ncbi.nlm.nih.gov domain with instance path
                    fixed_url = f"https://pmc.ncbi.nlm.nih.gov/articles/instance/{numeric_id}/bin/{filename}"
                    print(f"    Normalized PMC URL: {url} -> {fixed_url}")
                    return fixed_url

        return url

    def get_pmc_supplement_url_variants(self, url: str, base_url: str) -> list:
        """
        Generate multiple URL variants to try for PMC supplement downloads.

        PMC has multiple URL formats and has migrated domains, so we need to
        try several patterns to find the working URL.

        Args:
            url: The original supplement URL
            base_url: The base PMC article URL

        Returns:
            List of URL variants to try (in order of preference)
        """
        variants = []

        # Extract PMCID and numeric ID from base URL
        pmcid_match = re.search(r'/(?:pmc/)?articles/(PMC\d+)', base_url)
        if not pmcid_match:
            return [url]

        pmcid = pmcid_match.group(1)
        numeric_id = pmcid.replace('PMC', '')

        # Extract filename from URL
        instance_match = re.search(r'/articles/instance/\d+/bin/(.+)$', url)
        bin_match = re.search(r'/bin/(.+)$', url)

        if instance_match:
            filename = instance_match.group(1)
        elif bin_match:
            filename = bin_match.group(1)
        else:
            # Can't extract filename, just return original
            return [url]

        # URL variants to try (in order of preference):
        # 1. New domain with instance path (most likely to work for newer structure)
        variants.append(f"https://pmc.ncbi.nlm.nih.gov/articles/instance/{numeric_id}/bin/{filename}")
        # 2. Old domain with instance path
        variants.append(f"https://www.ncbi.nlm.nih.gov/pmc/articles/instance/{numeric_id}/bin/{filename}")
        # 3. New domain with PMC path
        variants.append(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/bin/{filename}")
        # 4. Old domain with PMC path
        variants.append(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/bin/{filename}")
        # 5. Europe PMC (often bypasses PMC proof-of-work for downloads)
        variants.append(f"https://europepmc.org/articles/{pmcid}/bin/{filename}")
        if filename.lower().endswith(".pdf"):
            variants.append(f"https://europepmc.org/articles/{pmcid}/pdf/{filename}")
        # 6. Original URL as fallback
        if url not in variants:
            variants.append(url)

        return variants

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
                # Try multiple possible paths in the JSON structure
                supp_data = (
                    data.get('article', {}).get('supplementaryMaterials', {}).get('supplementaryMaterial', []) or
                    data.get('supplementaryMaterials', {}).get('supplementaryMaterial', []) or
                    data.get('supplementaryMaterial', []) or
                    []
                )
                for item in supp_data:
                    url = item.get('downloadUrl') or item.get('url') or item.get('href')
                    filename = item.get('title') or item.get('name') or item.get('filename')
                    if url and filename and not any(f['name'] == filename for f in found_files):
                        found_files.append({'url': url, 'name': filename})
                        print(f"    Found supplement in JSON data: {filename}")
                if found_files:
                    return found_files
            except (json.JSONDecodeError, AttributeError, TypeError):
                continue

        # 2. Look for sections with "supplement" or "supporting" in headings
        supp_keywords = ['supplement', 'supporting', 'additional', 'appendix']
        for heading in soup.find_all(['h1', 'h2', 'h3', 'h4', 'h5', 'h6']):
            heading_text = heading.get_text().lower()
            if any(keyword in heading_text for keyword in supp_keywords):
                # Look for links in the section containing this heading
                section = heading.find_next_sibling() or heading.parent
                if section:
                    for link in section.find_all('a', href=True):
                        href = link['href']
                        link_text = link.get_text().lower()
                        # Check if it's a file link
                        if any(ext in href.lower() for ext in ['.pdf', '.docx', '.xlsx', '.zip', '.csv', '.txt', '.doc']):
                            url = urljoin(base_url, href)
                            filename = link.get_text().strip() or Path(urlparse(url).path).name
                            if filename and not any(f['name'] == filename for f in found_files):
                                found_files.append({'url': url, 'name': filename})
                                print(f"    Found supplement in section: {filename}")

        # 3. Regex for "mmc" (multimedia component) links, a common Elsevier pattern
        # More flexible pattern to catch various MMC link formats
        mmc_patterns = [
            r'href="(/cms/attachment/[^"]+/mmc\d+\.(?:pdf|docx|xlsx|zip|csv|txt))"',
            r'href="([^"]*mmc\d+\.(?:pdf|docx|xlsx|zip|csv|txt))"',
            r'href="([^"]*attachment[^"]*mmc[^"]*\.(?:pdf|docx|xlsx|zip|csv|txt))"',
        ]
        for pattern in mmc_patterns:
            mmc_links = re.findall(pattern, html, re.IGNORECASE)
            for link in mmc_links:
                url = urljoin(base_url, link)
                filename = Path(urlparse(url).path).name
                if filename and not any(f['name'] == filename for f in found_files):
                    found_files.append({'url': url, 'name': filename})
                    print(f"    Found MMC link via regex: {filename}")
        if found_files:
            return found_files

        # 4. Look for links with "supplement" or "mmc" in the text or href
        for link in soup.find_all('a', href=True):
            link_text = link.get_text().lower()
            href = link['href'].lower()
            is_supplement = (
                any(keyword in link_text for keyword in ['supplement', 'supporting', 'additional', 'mmc']) or
                'mmc' in href or
                'supplement' in href or
                'attachment' in href
            )
            is_file = any(href.endswith(ext) for ext in ['.pdf', '.docx', '.xlsx', '.zip', '.csv', '.txt', '.doc'])
            if is_supplement and is_file:
                url = urljoin(base_url, link['href'])
                filename = link.get_text().strip() or Path(urlparse(url).path).name
                if filename and not any(f['name'] == filename for f in found_files):
                    found_files.append({'url': url, 'name': filename})
                    print(f"    Found supplement link: {filename}")

        # 5. Look for specific "Download file" links in ScienceDirect as a fallback
        for link in soup.find_all('a', class_='S_C_9cf8451f', href=True):
            if 'Download file' in link.get_text():
                url = urljoin(base_url, link['href'])
                filename = link.get_text().replace('Download file', '').strip()
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
            List of supplement file dictionaries with 'url', 'name', and optionally
            'base_url' (for PMC pages to enable URL variant generation) keys
        """
        print("  Scraping with scrape_generic_supplements...")
        soup = BeautifulSoup(html, 'html.parser')
        found_files = []

        keywords = ['supplement', 'supporting', 'appendix', 'additional file']
        file_extensions = ['.pdf', '.docx', '.xlsx', '.csv', '.zip', '.rar', '.gz', '.txt', '.doc']

        # Check if this is a PMC page
        is_pmc_page = 'ncbi.nlm.nih.gov/pmc/articles/' in base_url or 'pmc.ncbi.nlm.nih.gov/articles/' in base_url

        for link in soup.find_all('a', href=True):
            link_text = link.get_text().lower()
            href = link['href'].lower()

            # Check if link text or href contain any of the keywords or file extensions
            is_supplement_link = any(keyword in link_text for keyword in keywords)
            is_file_link = any(href.endswith(ext) for ext in file_extensions)

            if is_supplement_link or is_file_link:
                try:
                    original_url = urljoin(base_url, link['href'])
                    # Normalize PMC URLs to use correct path format
                    url = self._normalize_pmc_url(original_url, base_url)
                    filename = Path(urlparse(url).path).name

                    # Skip invalid filenames
                    if not filename:
                        continue

                    # Skip entries that look like PMCIDs (e.g., "PMC3049907")
                    if re.match(r'^PMC\d+$', filename, re.IGNORECASE):
                        continue

                    # Skip entries that don't have a file extension and aren't actual files
                    has_extension = any(filename.lower().endswith(ext) for ext in file_extensions)
                    if not has_extension and is_supplement_link:
                        # Only skip non-extension links if they also look like IDs
                        if re.match(r'^[A-Z]+\d+$', filename, re.IGNORECASE):
                            continue

                    # Basic filtering to avoid irrelevant links
                    if not any(f['name'] == filename for f in found_files):
                        file_info = {'url': url, 'name': filename}
                        # Store original URL and base URL for PMC pages to enable URL variant generation
                        if is_pmc_page:
                            file_info['original_url'] = original_url
                            file_info['base_url'] = base_url
                        found_files.append(file_info)
                        print(f"    Found potential supplement: {filename}")
                except Exception:
                    continue  # Ignore malformed URLs

        return found_files

    # ==========================================================================
    # FULL-TEXT EXTRACTION METHODS
    # These extract main article content from publisher pages for free articles
    # without PMCIDs.
    # ==========================================================================

    def extract_fulltext_nature(self, html: str, base_url: str) -> Tuple[Optional[str], str]:
        """
        Extract full-text content from a Nature journal page.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL of the article

        Returns:
            Tuple of (markdown_content, title)
        """
        print("  Extracting full text from Nature article...")
        soup = BeautifulSoup(html, 'html.parser')
        markdown = "# MAIN TEXT\n\n"

        # Extract title
        title = ""
        title_elem = soup.find('h1', class_='c-article-title') or soup.find('h1')
        if title_elem:
            title = title_elem.get_text().strip()
            markdown += f"## {title}\n\n"

        # Extract abstract
        abstract = soup.find('div', id='Abs1-content') or soup.find('section', {'data-title': 'Abstract'})
        if abstract:
            markdown += "### Abstract\n\n"
            markdown += abstract.get_text().strip() + "\n\n"

        # Extract main article content
        article_body = soup.find('div', class_='c-article-body') or soup.find('main')
        if article_body:
            # Find all sections
            sections = article_body.find_all(['section', 'div'], class_=re.compile(r'c-article-section'))
            if not sections:
                sections = article_body.find_all('section')

            for section in sections:
                # Get section title
                section_title = section.find(['h2', 'h3', 'h4'])
                if section_title:
                    markdown += f"### {section_title.get_text().strip()}\n\n"

                # Get paragraphs
                for p in section.find_all('p', recursive=False):
                    text = p.get_text().strip()
                    if text:
                        markdown += f"{text}\n\n"

        # Fallback: just extract all paragraph text
        if len(markdown) < 500:
            markdown = "# MAIN TEXT\n\n"
            if title:
                markdown += f"## {title}\n\n"

            for p in soup.find_all('p'):
                text = p.get_text().strip()
                if len(text) > 50:  # Skip very short paragraphs (likely nav elements)
                    markdown += f"{text}\n\n"

        return markdown if len(markdown) > 200 else None, title

    def extract_fulltext_elsevier(self, html: str, base_url: str) -> Tuple[Optional[str], str]:
        """
        Extract full-text content from an Elsevier/ScienceDirect page.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL of the article

        Returns:
            Tuple of (markdown_content, title)
        """
        print("  Extracting full text from Elsevier/ScienceDirect article...")
        soup = BeautifulSoup(html, 'html.parser')
        markdown = "# MAIN TEXT\n\n"

        # Extract title
        title = ""
        title_elem = soup.find('span', class_='title-text') or soup.find('h1', class_='svTitle')
        if title_elem:
            title = title_elem.get_text().strip()
            markdown += f"## {title}\n\n"

        # Try to extract from JSON data embedded in page (ScienceDirect often has this)
        for script in soup.find_all('script', type='application/json'):
            try:
                data = json.loads(script.string)
                # Look for article content in various possible JSON paths
                article_data = data.get('article', {})
                if article_data:
                    # Extract abstract
                    abstract_data = article_data.get('abstract', {})
                    if abstract_data:
                        markdown += "### Abstract\n\n"
                        if isinstance(abstract_data, dict):
                            abstract_text = abstract_data.get('content', '')
                        else:
                            abstract_text = str(abstract_data)
                        # Clean HTML tags from abstract
                        abstract_soup = BeautifulSoup(abstract_text, 'html.parser')
                        markdown += abstract_soup.get_text().strip() + "\n\n"

                    # Extract body sections
                    body_data = article_data.get('body', {})
                    if body_data and isinstance(body_data, dict):
                        content = body_data.get('content', [])
                        if isinstance(content, list):
                            for section in content:
                                if isinstance(section, dict):
                                    sec_title = section.get('label', '') or section.get('title', '')
                                    if sec_title:
                                        markdown += f"### {sec_title}\n\n"
                                    sec_content = section.get('content', '')
                                    if sec_content:
                                        content_soup = BeautifulSoup(sec_content, 'html.parser')
                                        markdown += content_soup.get_text().strip() + "\n\n"
                if len(markdown) > 500:
                    return markdown, title
            except (json.JSONDecodeError, AttributeError, TypeError):
                continue

        # Fallback: HTML-based extraction
        # Extract abstract
        abstract = soup.find('div', class_='abstract') or soup.find('div', id='abstracts')
        if abstract:
            markdown += "### Abstract\n\n"
            markdown += abstract.get_text().strip() + "\n\n"

        # Extract main body content
        body = soup.find('div', id='body') or soup.find('div', class_='Body')
        if body:
            for section in body.find_all(['section', 'div'], class_=re.compile(r'section')):
                section_title = section.find(['h2', 'h3', 'h4'])
                if section_title:
                    markdown += f"### {section_title.get_text().strip()}\n\n"

                for p in section.find_all('p'):
                    text = p.get_text().strip()
                    if text:
                        markdown += f"{text}\n\n"

        # Ultimate fallback: int all paragraphs
        if len(markdown) < 500:
            markdown = "# MAIN TEXT\n\n"
            if title:
                markdown += f"## {title}\n\n"

            content_divs = soup.find_all('div', class_=re.compile(r'(Body|content|article)', re.I))
            for div in content_divs:
                for p in div.find_all('p'):
                    text = p.get_text().strip()
                    if len(text) > 50:
                        markdown += f"{text}\n\n"

        return markdown if len(markdown) > 200 else None, title

    def extract_fulltext_wiley(self, html: str, base_url: str) -> Tuple[Optional[str], str]:
        """
        Extract full-text content from a Wiley Online Library page.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL of the article

        Returns:
            Tuple of (markdown_content, title)
        """
        print("  Extracting full text from Wiley article...")
        soup = BeautifulSoup(html, 'html.parser')
        markdown = "# MAIN TEXT\n\n"

        # Extract title
        title = ""
        title_elem = soup.find('h1', class_='citation__title') or soup.find('h1')
        if title_elem:
            title = title_elem.get_text().strip()
            markdown += f"## {title}\n\n"

        # Extract abstract
        abstract = soup.find('section', class_='article-section__abstract') or soup.find('div', class_='abstract')
        if abstract:
            markdown += "### Abstract\n\n"
            markdown += abstract.get_text().strip() + "\n\n"

        # Extract main body
        body = soup.find('section', class_='article-section__content') or soup.find('div', class_='article__body')
        if body:
            for section in body.find_all(['section', 'div']):
                section_title = section.find(['h2', 'h3', 'h4'])
                if section_title:
                    markdown += f"### {section_title.get_text().strip()}\n\n"

                for p in section.find_all('p', recursive=False):
                    text = p.get_text().strip()
                    if text:
                        markdown += f"{text}\n\n"

        # Fallback
        if len(markdown) < 500:
            for p in soup.find_all('p'):
                text = p.get_text().strip()
                if len(text) > 50:
                    markdown += f"{text}\n\n"

        return markdown if len(markdown) > 200 else None, title

    def extract_fulltext_generic(self, html: str, base_url: str) -> Tuple[Optional[str], str]:
        """
        Generic full-text extraction for any publisher page.

        Uses heuristics to find article content in HTML.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL of the article

        Returns:
            Tuple of (markdown_content, title)
        """
        print("  Extracting full text using generic scraper...")
        soup = BeautifulSoup(html, 'html.parser')
        markdown = "# MAIN TEXT\n\n"

        # Remove script, style, nav, footer, header elements
        for element in soup.find_all(['script', 'style', 'nav', 'footer', 'header', 'aside']):
            element.decompose()

        # Extract title
        title = ""
        title_elem = soup.find('h1') or soup.find('title')
        if title_elem:
            title = title_elem.get_text().strip()
            # Clean up title (remove site name etc)
            title = re.sub(r'\s*[-|].*$', '', title)
            markdown += f"## {title}\n\n"

        # Look for common article containers
        article_containers = [
            soup.find('article'),
            soup.find('div', class_=re.compile(r'(article|content|body|main)', re.I)),
            soup.find('main'),
            soup.find('div', id=re.compile(r'(article|content|body|main)', re.I)),
        ]

        content_found = False
        for container in article_containers:
            if container:
                # Extract sections
                for section in container.find_all(['section', 'div'], recursive=False):
                    section_title = section.find(['h2', 'h3', 'h4'])
                    if section_title:
                        markdown += f"### {section_title.get_text().strip()}\n\n"

                    for p in section.find_all('p'):
                        text = p.get_text().strip()
                        if len(text) > 30:
                            markdown += f"{text}\n\n"
                            content_found = True

                if content_found:
                    break

        # Ultimate fallback: int all paragraphs
        if not content_found or len(markdown) < 500:
            markdown = "# MAIN TEXT\n\n"
            if title:
                markdown += f"## {title}\n\n"

            for p in soup.find_all('p'):
                text = p.get_text().strip()
                # Filter out short paragraphs that are likely navigation or UI elements
                if len(text) > 100:
                    markdown += f"{text}\n\n"

        return markdown if len(markdown) > 200 else None, title

    def extract_fulltext(self, html: str, base_url: str) -> Tuple[Optional[str], str]:
        """
        Main entry point for full-text extraction.
        Routes to domain-specific extractors based on URL.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL of the article

        Returns:
            Tuple of (markdown_content, title)
        """
        domain = urlparse(base_url).netloc.lower()

        if 'nature.com' in domain:
            return self.extract_fulltext_nature(html, base_url)
        elif any(d in domain for d in ['sciencedirect.com', 'elsevier.com', 'gimjournal.org']):
            return self.extract_fulltext_elsevier(html, base_url)
        elif 'wiley.com' in domain or 'onlinelibrary.wiley.com' in domain:
            return self.extract_fulltext_wiley(html, base_url)
        else:
            return self.extract_fulltext_generic(html, base_url)
