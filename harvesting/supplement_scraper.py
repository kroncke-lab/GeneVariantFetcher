r"""
Supplement Scraper Module
========================

Web scraping logic for extracting supplemental files and full-text content
from publisher websites.

CRITICAL CONTEXT - 94.6% Variant Extraction Failure
----------------------------------------------------
The Gene Variant Fetcher pipeline fails to extract variants from ~94.6% of papers.
This module is a major bottleneck because most variant data is in supplementary
files (Excel tables, PDFs), NOT the main article text.

PUBLISHER COVERAGE SUMMARY
==========================

SPECIFIC HANDLERS IMPLEMENTED (4 publishers):
---------------------------------------------
1. Nature (scrape_nature_supplements)
   - Looks for: #supplementary-information section, h2 with "Supplementary Information/Data"
   - Coverage: Good for Nature, Nature Genetics, Nature Communications
   - Pattern: /articles/ URLs in the supplementary section

2. Elsevier (scrape_elsevier_supplements) - Also covers ScienceDirect, GIM Journal
   - Looks for:
     * JSON in <script type="application/json"> with supplementaryMaterials
     * Sections with supplement/supporting/additional/appendix headings
     * MMC (multimedia component) links via regex: mmc\d+\.(pdf|docx|xlsx|...)
     * Links with class "S_C_9cf8451f" (ScienceDirect specific)
   - Coverage: Good for most Elsevier journals

3. Karger (scrape_karger_supplements)
   - Looks for:
     * Supplement sections/divs with supplement|supp-material|additional classes
     * Constructs supplement URL from DOI: karger.com/doi/suppl/10.1159/XXXXXX
     * Makes HTTP request to supplement page to find files
   - Coverage: POOR - See "Why Karger Fails" below
   - Browser fallback: harvesting/browser_supplement_fetcher.py

4. SPRINGER/BMC (scrape_springer_supplements) - IMPLEMENTED 2026-02-01
   - Looks for:
     * Section classes: c-article-supplementary, SupplementaryMaterial, additional-files
     * Heading patterns: "Electronic Supplementary Material", "Additional files", "ESM"
     * /MediaObjects/ links (Springer's CDN for ESM files)
     * MOESM naming pattern in URLs and link text
     * data-track attributes for download links
   - Coverage: Good for link.springer.com, biomedcentral.com, springeropen.com
   - Filters out: citation links, anchor-only URLs, main article PDFs
   - Tested with: 10.1186/s12864-019-6413-7 (5 supplements found)

NO SPECIFIC HANDLERS (Remaining Gaps):
--------------------------------------
These publishers fall through to the generic scraper, which often fails:

5. OXFORD ACADEMIC (scrape_oxford_supplements) - IMPLEMENTED 2026-02-01
   - Looks for:
     * "Supplementary data" section/tab
     * /downloadSupplement?file=... endpoint URLs (Oxford's download pattern)
     * Links with supplement/supplementary in href or text
     * data-doi links in supplementary sections
   - Coverage: Good for academic.oup.com journals
   - URL pattern: /downloadSupplement?file=XXXXX_supplementary_data.xlsx

6. WILEY (onlinelibrary.wiley.com)
   - NO supplement handler (only full-text extractor exists)
   - Supplement location: "Supporting Information" section
   - Actual URLs: /action/downloadSupplement?doi=...&file=...
   - Similar issues to Oxford

7. AHA JOURNALS (ahajournals.org) - Circulation, etc.
   - NO dedicated handler
   - Supplements at: /doi/suppl/10.1161/XXXXX

8. PLOS, Frontiers, MDPI
   - NO dedicated handlers
   - Usually better handled by generic due to cleaner HTML structure

GENERIC SCRAPER (scrape_generic_supplements)
============================================
Keywords searched: "supplement", "supporting", "appendix", "additional file"
Extensions matched: .pdf, .docx, .xlsx, .csv, .zip, .rar, .gz, .txt, .doc

CRITICAL LIMITATION: The generic scraper only finds supplements when:
- Link text contains one of the keywords, OR
- Link href ends with a known file extension

This misses:
- JavaScript-rendered supplement sections
- Supplements on separate pages (requires navigation)
- Downloads via API endpoints (e.g., /downloadSupplement?...)
- Publisher-specific naming (ESM, MOESM, mmc, etc.)
- Links with generic text like "Download" or just icons

WHY KARGER FAILS
================
1. DOI extraction unreliable - DOI may not be on page or in expected format
2. Supplement page (karger.com/doi/suppl/...) often returns:
   - 403 Forbidden (anti-bot)
   - Redirect to login
   - Different HTML structure than expected
3. No browser fallback integration in main pipeline
4. HTTP requests from scraper blocked by Cloudflare

OXFORD ACADEMIC - HANDLER IMPLEMENTED (2026-02-01)
===================================================
scrape_oxford_supplements() now handles academic.oup.com:
- Detects "Supplementary data" sections
- Parses /downloadSupplement?file=... URLs
- Extracts file parameter for proper naming
Note: Some Oxford supplements may still require authentication.

RECOMMENDED FIXES (Priority Order)
==================================
1. ✅ DONE: scrape_springer_supplements handler implemented (2026-02-01)
   - Handles link.springer.com, biomedcentral.com, springeropen.com
   - Extracts /MediaObjects/ URLs with MOESM naming
   - Routing added to doi_resolver.py

2. Add scrape_oxford_supplements handler:
   - Look for "Supplementary data" section
   - Handle /downloadSupplement endpoint URLs
   - May need to follow "Supplementary data" tab links

3. Improve generic scraper:
   - Add more keywords: "ESM", "MOESM", "electronic supplementary", "supporting information"
   - Match /downloadSupplement, /MediaObjects/, /figures/ paths
   - Try following links that say "Supplementary" even without file extension

4. Browser fallback integration:
   - Karger browser fetcher exists but not integrated into main pipeline
   - Need similar for Oxford, Springer (Cloudflare protection common)

5. Consider Playwright for all publishers:
   - JavaScript rendering would solve most issues
   - Expensive in terms of time/resources
   - browser_supplement_fetcher.py has Karger-specific implementation

FULL-TEXT EXTRACTORS
====================
Also includes extractors for article main text (for free articles without PMCIDs):
- extract_fulltext_nature
- extract_fulltext_elsevier
- extract_fulltext_wiley
- extract_fulltext_generic

These route via extract_fulltext() based on domain.

FILE STRUCTURE
==============
- scrape_*_supplements: Returns List[Dict] with 'url' and 'name' keys
- extract_fulltext_*: Returns Tuple[Optional[str], str] - (markdown, title)
- _normalize_pmc_url: Handles PMC URL variants
- get_pmc_supplement_url_variants: Generates fallback URLs for PMC downloads

Related files:
- harvesting/doi_resolver.py - Routes DOIs to appropriate scrapers
- harvesting/browser_supplement_fetcher.py - Playwright-based Karger fetcher
- cli/browser_fetch.py - General browser-based PDF fetcher with selector patterns

Author: Gene Variant Fetcher Team
Last audit: 2026-02-01
"""

import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import urljoin, urlparse

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
        if (
            "ncbi.nlm.nih.gov/pmc/articles/" in base_url
            or "pmc.ncbi.nlm.nih.gov/articles/" in base_url
        ):
            # Extract the PMCID from the base URL
            import re

            pmcid_match = re.search(r"/(?:pmc/)?articles/(PMC\d+)", base_url)
            if pmcid_match:
                pmcid = pmcid_match.group(1)
                numeric_id = pmcid.replace("PMC", "")

                # Fix URLs that use /articles/instance/{numeric_id}/bin/ format
                instance_match = re.search(r"/articles/instance/(\d+)/bin/(.+)$", url)
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
        pmcid_match = re.search(r"/(?:pmc/)?articles/(PMC\d+)", base_url)
        if not pmcid_match:
            return [url]

        pmcid = pmcid_match.group(1)
        numeric_id = pmcid.replace("PMC", "")

        # Extract filename from URL
        instance_match = re.search(r"/articles/instance/\d+/bin/(.+)$", url)
        bin_match = re.search(r"/bin/(.+)$", url)

        if instance_match:
            filename = instance_match.group(1)
        elif bin_match:
            filename = bin_match.group(1)
        else:
            # Can't extract filename, just return original
            return [url]

        # URL variants to try (in order of preference):
        # 1. New domain with instance path (most likely to work for newer structure)
        variants.append(
            f"https://pmc.ncbi.nlm.nih.gov/articles/instance/{numeric_id}/bin/{filename}"
        )
        # 2. Old domain with instance path
        variants.append(
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/instance/{numeric_id}/bin/{filename}"
        )
        # 3. New domain with PMC path
        variants.append(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/bin/{filename}")
        # 4. Old domain with PMC path
        variants.append(
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/bin/{filename}"
        )
        # 5. Europe PMC (often bypasses PMC proof-of-work for downloads)
        variants.append(f"https://europepmc.org/articles/{pmcid}/bin/{filename}")
        if filename.lower().endswith(".pdf"):
            variants.append(f"https://europepmc.org/articles/{pmcid}/pdf/{filename}")
        # 6. Original URL as fallback
        if url not in variants:
            variants.append(url)

        return variants

    def scrape_springer_supplements(self, html: str, base_url: str) -> List[Dict]:
        """
        Scrape supplemental files from Springer/BMC journal pages.

        Springer/BMC (link.springer.com, biomedcentral.com) is one of the largest
        publishers. This handler addresses a major gap in the 94.6% failure rate.

        Key patterns:
        - Section heading: "Electronic Supplementary Material" or "Additional files"
        - Section classes: "c-article-supplementary", "Supplementary", "SupplementaryMaterial"
        - Link patterns:
          * /MediaObjects/XXXXX_YYYY_ZZZZ_MOESM1_ESM.pdf
          * /article/10.1007/.../figures/N (figure files)
        - Link text: "ESM 1", "MOESM1", "Additional file 1", etc.
        - data-track attributes for analytics (useful for finding links)

        Args:
            html: HTML content of the publisher page
            base_url: Base URL for resolving relative links

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys
        """
        print("  Scraping with scrape_springer_supplements...")
        soup = BeautifulSoup(html, "html.parser")
        found_files = []
        seen_urls = set()  # Track URLs to avoid duplicates

        def add_file(url: str, name: str) -> bool:
            """Add a file if not already seen and passes validation. Returns True if added."""
            # Normalize URL for deduplication (strip trailing slashes and anchors for comparison)
            parsed = urlparse(url)

            # Skip anchor-only links (e.g., article.html#MOESM1) - these aren't actual files
            if parsed.fragment and not any(
                ext in parsed.path.lower()
                for ext in [
                    ".pdf",
                    ".doc",
                    ".xls",
                    ".csv",
                    ".zip",
                    ".txt",
                    ".pptx",
                    ".xlsx",
                    ".docx",
                ]
            ):
                return False

            # Skip citation/reference links
            if "citation-needed" in url or "refman" in url or "format=refman" in url:
                return False

            # Skip main article PDF (we want supplements only, not the paper itself)
            # Main article PDFs typically have paths like /content/pdf/DOI.pdf
            if "/content/pdf/" in url and "MOESM" not in url and "ESM" not in url:
                return False

            # Create a normalized key for deduplication
            norm_url = url.split("#")[0].rstrip("/")  # Remove anchor and trailing slash
            if norm_url in seen_urls:
                return False
            seen_urls.add(norm_url)

            # Skip if name looks like a DOI fragment without extension
            if re.match(r"^s?\d{5}-\d{3}-\d{4}-\d+$", name) and not any(
                ext in name.lower()
                for ext in [
                    ".pdf",
                    ".doc",
                    ".xls",
                    ".csv",
                    ".zip",
                    ".txt",
                    ".pptx",
                    ".xlsx",
                    ".docx",
                ]
            ):
                return False

            found_files.append({"url": url, "name": name})
            print(f"      Found supplement: {name}")
            return True

        # 1. Look for ESM section by class (most reliable)
        # Springer uses various class names for supplementary sections
        esm_section = None
        section_classes = [
            r"c-article-supplementary",  # Nature/Springer style
            r"SupplementaryMaterial",
            r"Supplementary",
            r"additional-files",
            r"supplementary-content",
            r"esm-container",
        ]

        for class_pattern in section_classes:
            esm_section = soup.find(
                ["section", "div", "aside"],
                class_=re.compile(class_pattern, re.IGNORECASE),
            )
            if esm_section:
                print(
                    f"    Found supplementary section with class matching: {class_pattern}"
                )
                break

        # 2. Fallback: Find by heading text
        if not esm_section:
            # Try various heading patterns
            heading_patterns = [
                r"Electronic\s+Supplementary\s+Material",
                r"Additional\s+files?",
                r"Supplementary\s+(Information|Material|Data|Files?)",
                r"Supporting\s+Information",
                r"ESM",
            ]
            for pattern in heading_patterns:
                heading = soup.find(
                    ["h2", "h3", "h4", "h5"], string=re.compile(pattern, re.IGNORECASE)
                )
                if heading:
                    # Get the parent section or the next sibling containing links
                    esm_section = heading.find_parent(["section", "div", "article"])
                    if not esm_section:
                        esm_section = heading.find_next_sibling(
                            ["div", "section", "ul", "ol"]
                        )
                    if esm_section:
                        print(
                            f"    Found supplementary section via heading: {heading.get_text().strip()[:50]}"
                        )
                        break

        # 3. Extract links from ESM section
        if esm_section:
            for link in esm_section.find_all("a", href=True):
                href = link["href"]
                link_text = link.get_text().strip()

                # Check if this looks like a supplement download link
                is_media_object = "/MediaObjects/" in href or "/media/" in href.lower()
                is_esm_link = re.search(
                    r"(MOESM|ESM|supplementary|additional)",
                    href + link_text,
                    re.IGNORECASE,
                )
                has_extension = any(
                    ext in href.lower()
                    for ext in [
                        ".pdf",
                        ".doc",
                        ".xls",
                        ".csv",
                        ".zip",
                        ".txt",
                        ".pptx",
                        ".xlsx",
                        ".docx",
                    ]
                )

                if is_media_object or is_esm_link or has_extension:
                    url = urljoin(base_url, href)
                    # Extract filename from URL or use link text
                    filename = Path(urlparse(url).path).name
                    if not filename or filename == "":
                        filename = link_text if link_text else "supplement"
                    # Clean up filename
                    filename = re.sub(r"[^\w\-_\.\s]", "", filename).strip()
                    if not filename:
                        filename = "supplement"
                    add_file(url, filename)

        # 4. Scan entire page for /MediaObjects/ links (Springer's CDN for supplements)
        for link in soup.find_all(
            "a", href=re.compile(r"/MediaObjects/", re.IGNORECASE)
        ):
            href = link["href"]
            url = urljoin(base_url, href)
            filename = Path(urlparse(url).path).name
            if filename:
                add_file(url, filename)

        # 5. Look for links with data-track attributes (Springer analytics)
        for link in soup.find_all(
            "a",
            attrs={
                "data-track-action": re.compile(r"(download|supplement)", re.IGNORECASE)
            },
        ):
            href = link.get("href", "")
            if href:
                url = urljoin(base_url, href)
                link_text = link.get_text().strip()
                filename = Path(urlparse(url).path).name or link_text or "supplement"
                add_file(url, filename)

        # 6. Look for MOESM pattern in any link on the page
        for link in soup.find_all("a", href=re.compile(r"MOESM\d*", re.IGNORECASE)):
            href = link["href"]
            url = urljoin(base_url, href)
            filename = Path(urlparse(url).path).name
            link_text = link.get_text().strip()
            if filename:
                add_file(url, filename)
            elif link_text:
                add_file(url, link_text)

        # 7. Check for BMC-style "Additional file" links
        for link in soup.find_all(
            "a", string=re.compile(r"Additional\s+file\s*\d*", re.IGNORECASE)
        ):
            href = link.get("href", "")
            if href:
                url = urljoin(base_url, href)
                link_text = link.get_text().strip()
                filename = Path(urlparse(url).path).name or link_text
                add_file(url, filename)

        # 8. Look for links in figure/table supplement sections
        for container in soup.find_all(
            ["figure", "div"],
            class_=re.compile(r"(figure|table).*supplement", re.IGNORECASE),
        ):
            for link in container.find_all("a", href=True):
                href = link["href"]
                if any(
                    ext in href.lower()
                    for ext in [".pdf", ".doc", ".xls", ".csv", ".zip", ".txt"]
                ):
                    url = urljoin(base_url, href)
                    filename = Path(urlparse(url).path).name
                    add_file(url, filename)

        if found_files:
            print(f"    Total Springer supplements found: {len(found_files)}")
            return found_files

        print("    No specific Springer/BMC supplements found. Trying generic scan.")
        return self.scrape_generic_supplements(html, base_url)

    def scrape_oxford_supplements(self, html: str, base_url: str) -> List[Dict]:
        """
        Scrape supplemental files from Oxford Academic journal pages.

        Oxford Academic (academic.oup.com) uses a distinctive pattern:
        - Supplement section: "Supplementary data" or "Supplementary Material"
        - Download URLs: /downloadSupplement?file=XXXXX_supplementary_data.xlsx
        - The file parameter contains the actual filename

        Key patterns:
        - Section classes: "supplementary-data", "supplementary-material"
        - Link patterns: /downloadSupplement?file=..., /oup/backfile/...
        - data-doi attributes on supplement links

        Args:
            html: HTML content of the publisher page
            base_url: Base URL for resolving relative links

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys
        """
        print("  Scraping with scrape_oxford_supplements...")
        soup = BeautifulSoup(html, "html.parser")
        found_files = []
        seen_urls = set()  # Track URLs to avoid duplicates

        def add_file(url: str, name: str) -> bool:
            """Add a file if not already seen. Returns True if added."""
            # Normalize for deduplication
            norm_url = url.split("#")[0]
            if norm_url in seen_urls:
                return False
            seen_urls.add(norm_url)

            # Clean up the name
            name = re.sub(r"[^\w\-_\.\s]", "", name).strip()
            if not name:
                name = "supplement"

            found_files.append({"url": url, "name": name})
            print(f"      Found supplement: {name}")
            return True

        # 1. Look for /downloadSupplement links anywhere on page (most reliable pattern)
        for link in soup.find_all(
            "a", href=re.compile(r"downloadSupplement", re.IGNORECASE)
        ):
            href = link["href"]
            url = urljoin(base_url, href)

            # Extract filename from the 'file' parameter
            file_match = re.search(r"file=([^&]+)", href)
            if file_match:
                filename = file_match.group(1)
                # URL decode if needed
                from urllib.parse import unquote

                filename = unquote(filename)
            else:
                # Fallback to link text
                filename = link.get_text().strip() or "supplement"

            add_file(url, filename)

        # 2. Look for supplementary data section by class
        supp_section = None
        section_classes = [
            r"supplementary-data",
            r"supplementary-material",
            r"supp-data",
            r"article-supplementary",
        ]

        for class_pattern in section_classes:
            supp_section = soup.find(
                ["section", "div", "aside"],
                class_=re.compile(class_pattern, re.IGNORECASE),
            )
            if supp_section:
                print(
                    f"    Found supplementary section with class matching: {class_pattern}"
                )
                break

        # 3. Fallback: Find by heading text
        if not supp_section:
            heading_patterns = [
                r"Supplementary\s+[Dd]ata",
                r"Supplementary\s+[Mm]aterial",
                r"Supporting\s+[Ii]nformation",
                r"Additional\s+[Ff]iles?",
            ]
            for pattern in heading_patterns:
                heading = soup.find(
                    ["h2", "h3", "h4", "h5"], string=re.compile(pattern, re.IGNORECASE)
                )
                if heading:
                    supp_section = heading.find_parent(["section", "div", "article"])
                    if not supp_section:
                        supp_section = heading.find_next_sibling(
                            ["div", "section", "ul", "ol"]
                        )
                    if supp_section:
                        print(
                            f"    Found supplementary section via heading: {heading.get_text().strip()[:50]}"
                        )
                        break

        # 4. Extract links from supplementary section
        if supp_section:
            for link in supp_section.find_all("a", href=True):
                href = link["href"]
                link_text = link.get_text().strip()

                # Check if this looks like a download link
                is_download = (
                    "download" in href.lower() or "download" in link_text.lower()
                )
                has_extension = any(
                    ext in href.lower()
                    for ext in [
                        ".pdf",
                        ".doc",
                        ".xls",
                        ".csv",
                        ".zip",
                        ".txt",
                        ".pptx",
                        ".xlsx",
                        ".docx",
                    ]
                )
                is_supplement = "suppl" in href.lower() or "suppl" in link_text.lower()

                if is_download or has_extension or is_supplement:
                    url = urljoin(base_url, href)

                    # Try to get filename from file parameter or URL path
                    file_match = re.search(r"file=([^&]+)", href)
                    if file_match:
                        from urllib.parse import unquote

                        filename = unquote(file_match.group(1))
                    else:
                        filename = Path(urlparse(url).path).name
                        if not filename or filename == "":
                            filename = link_text if link_text else "supplement"

                    add_file(url, filename)

        # 5. Look for Oxford backfile links (another common pattern)
        for link in soup.find_all(
            "a", href=re.compile(r"/oup/backfile/", re.IGNORECASE)
        ):
            href = link["href"]
            url = urljoin(base_url, href)
            filename = Path(urlparse(url).path).name
            if filename:
                add_file(url, filename)

        # 6. Look for links with data-doi attribute in supplementary context
        for link in soup.find_all("a", attrs={"data-doi": True}):
            href = link.get("href", "")
            if href and ("suppl" in href.lower() or "download" in href.lower()):
                url = urljoin(base_url, href)
                filename = (
                    link.get_text().strip()
                    or Path(urlparse(url).path).name
                    or "supplement"
                )
                add_file(url, filename)

        # 7. Look for links in supplementary list items
        for li in soup.find_all("li", class_=re.compile(r"suppl", re.IGNORECASE)):
            for link in li.find_all("a", href=True):
                href = link["href"]
                url = urljoin(base_url, href)
                filename = link.get_text().strip() or Path(urlparse(url).path).name
                if filename:
                    add_file(url, filename)

        if found_files:
            print(f"    Total Oxford supplements found: {len(found_files)}")
            return found_files

        print("    No specific Oxford supplements found. Trying generic scan.")
        return self.scrape_generic_supplements(html, base_url)

    def scrape_wiley_supplements(self, html: str, base_url: str) -> List[Dict]:
        """
        Scrape supplemental files from Wiley Online Library journal pages.

        Wiley (onlinelibrary.wiley.com) is a major publisher for genetics/medical journals.
        This handler addresses a gap in the 94.6% failure rate.

        Key patterns:
        - Section heading: "Supporting Information" or "Data Availability"
        - Section classes: "article-section__supporting", "article-section__supplement"
        - Link patterns:
          * /action/downloadSupplement?doi=...&file=...
          * Direct file links ending in .pdf, .xlsx, .docx, .zip
        - Link text: "Filename:", "Supporting Information", "Appendix S1", etc.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL for resolving relative links

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys
        """
        print("  Scraping with scrape_wiley_supplements...")
        soup = BeautifulSoup(html, "html.parser")
        found_files = []
        seen_urls = set()  # Track URLs to avoid duplicates

        def add_file(url: str, name: str) -> bool:
            """Add a file if not already seen. Returns True if added."""
            # Normalize for deduplication
            norm_url = url.split("#")[0]
            if norm_url in seen_urls:
                return False
            seen_urls.add(norm_url)

            # Clean up the name
            name = re.sub(r"[^\w\-_\.\s]", "", name).strip()
            if not name:
                name = "supplement"

            found_files.append({"url": url, "name": name})
            print(f"      Found supplement: {name}")
            return True

        # 1. Look for /action/downloadSupplement links anywhere on page (most reliable)
        for link in soup.find_all(
            "a", href=re.compile(r"/action/downloadSupplement", re.IGNORECASE)
        ):
            href = link["href"]
            url = urljoin(base_url, href)

            # Extract filename from the 'file' parameter
            file_match = re.search(r"file=([^&]+)", href)
            if file_match:
                from urllib.parse import unquote

                filename = unquote(file_match.group(1))
            else:
                # Fallback to link text or inner content
                filename = link.get_text().strip()
                # Wiley sometimes has "Filename: xxx.pdf" format
                filename_match = re.search(r"Filename:\s*(.+)", filename, re.IGNORECASE)
                if filename_match:
                    filename = filename_match.group(1).strip()
                if not filename:
                    filename = "supplement"

            add_file(url, filename)

        # 2. Look for "Supporting Information" section by class
        supp_section = None
        section_classes = [
            r"article-section__supporting",
            r"article-section__supplement",
            r"article-supplementary",
            r"supporting-info",
            r"data-availability",
        ]

        for class_pattern in section_classes:
            supp_section = soup.find(
                ["section", "div", "aside"],
                class_=re.compile(class_pattern, re.IGNORECASE),
            )
            if supp_section:
                print(
                    f"    Found supplementary section with class matching: {class_pattern}"
                )
                break

        # 3. Fallback: Find by heading text
        if not supp_section:
            heading_patterns = [
                r"Supporting\s+Information",
                r"Supplementary\s+(Material|Data|Files?)",
                r"Data\s+Availability",
                r"Additional\s+Supporting\s+Information",
                r"Appendix",
            ]
            for pattern in heading_patterns:
                heading = soup.find(
                    ["h2", "h3", "h4", "h5"], string=re.compile(pattern, re.IGNORECASE)
                )
                if heading:
                    supp_section = heading.find_parent(["section", "div", "article"])
                    if not supp_section:
                        supp_section = heading.find_next_sibling(
                            ["div", "section", "ul", "ol"]
                        )
                    if supp_section:
                        print(
                            f"    Found supplementary section via heading: {heading.get_text().strip()[:50]}"
                        )
                        break

        # 4. Extract links from supplementary section
        if supp_section:
            for link in supp_section.find_all("a", href=True):
                href = link["href"]
                link_text = link.get_text().strip()

                # Check if this looks like a download link
                is_download = (
                    "download" in href.lower() or "download" in link_text.lower()
                )
                has_extension = any(
                    ext in href.lower()
                    for ext in [
                        ".pdf",
                        ".doc",
                        ".xls",
                        ".csv",
                        ".zip",
                        ".txt",
                        ".pptx",
                        ".xlsx",
                        ".docx",
                    ]
                )
                is_supplement = (
                    "suppl" in href.lower()
                    or "suppl" in link_text.lower()
                    or "appendix" in link_text.lower()
                )

                if is_download or has_extension or is_supplement:
                    url = urljoin(base_url, href)

                    # Try to get filename from file parameter, URL path, or link text
                    file_match = re.search(r"file=([^&]+)", href)
                    if file_match:
                        from urllib.parse import unquote

                        filename = unquote(file_match.group(1))
                    else:
                        filename = Path(urlparse(url).path).name
                        if not filename or filename == "":
                            # Look for "Filename: xxx" in link text
                            filename_match = re.search(
                                r"Filename:\s*(.+)", link_text, re.IGNORECASE
                            )
                            if filename_match:
                                filename = filename_match.group(1).strip()
                            elif link_text:
                                filename = link_text
                            else:
                                filename = "supplement"

                    add_file(url, filename)

        # 5. Look for links with "supporting" or "supplementary" in text
        for link in soup.find_all(
            "a",
            string=re.compile(r"(supporting|supplementary|appendix)", re.IGNORECASE),
        ):
            href = link.get("href", "")
            if href and any(
                ext in href.lower()
                for ext in [".pdf", ".doc", ".xls", ".csv", ".zip", ".txt"]
            ):
                url = urljoin(base_url, href)
                filename = link.get_text().strip()
                # Extract "Filename: xxx" pattern
                filename_match = re.search(r"Filename:\s*(.+)", filename, re.IGNORECASE)
                if filename_match:
                    filename = filename_match.group(1).strip()
                if not filename:
                    filename = Path(urlparse(url).path).name
                add_file(url, filename)

        # 6. Look for data-article-supporting-info attributes (Wiley specific)
        for link in soup.find_all("a", attrs={"data-article-supporting-info": True}):
            href = link.get("href", "")
            if href:
                url = urljoin(base_url, href)
                filename = (
                    link.get_text().strip()
                    or Path(urlparse(url).path).name
                    or "supplement"
                )
                add_file(url, filename)

        # 7. Look for list items in supporting info sections
        for li in soup.find_all(
            "li", class_=re.compile(r"(supporting|supplement)", re.IGNORECASE)
        ):
            for link in li.find_all("a", href=True):
                href = link["href"]
                if any(
                    ext in href.lower()
                    for ext in [".pdf", ".doc", ".xls", ".csv", ".zip", ".txt"]
                ):
                    url = urljoin(base_url, href)
                    filename = link.get_text().strip()
                    filename_match = re.search(
                        r"Filename:\s*(.+)", filename, re.IGNORECASE
                    )
                    if filename_match:
                        filename = filename_match.group(1).strip()
                    if not filename:
                        filename = Path(urlparse(url).path).name
                    if filename:
                        add_file(url, filename)

        if found_files:
            print(f"    Total Wiley supplements found: {len(found_files)}")
            return found_files

        print("    No specific Wiley supplements found. Trying generic scan.")
        return self.scrape_generic_supplements(html, base_url)

    def scrape_karger_supplements(self, html: str, base_url: str) -> List[Dict]:
        """
        Scrape supplemental files from a Karger journal page.

        Karger supplements are typically at a separate URL:
        - Article: karger.com/crd/article/133/2/73/...
        - Supplements: karger.com/doi/suppl/10.1159/XXXXXX

        KNOWN ISSUES (contributing to 94.6% failure):
        1. DOI extraction often fails if DOI not in expected format
        2. Supplement page (karger.com/doi/suppl/...) returns 403 (Cloudflare)
        3. HTTP requests blocked by anti-bot protection
        4. Returns empty list when HTTP fails, relies on browser fallback
           that isn't integrated into the main pipeline

        BROWSER FALLBACK: harvesting/browser_supplement_fetcher.py has
        fetch_karger_supplements_browser() but it's not called from the
        main orchestrator.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL for resolving relative links

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys.
            Often returns [] due to access issues.
        """
        print("  Scraping with scrape_karger_supplements...")
        soup = BeautifulSoup(html, "html.parser")
        found_files = []

        # 1. Look for supplement links on the page
        # Karger often has a "Supplementary Material" section
        supp_sections = soup.find_all(
            ["section", "div"],
            class_=re.compile(r"supplement|supp-material|additional", re.IGNORECASE),
        )

        for section in supp_sections:
            print("    Found supplementary section")
            for link in section.find_all("a", href=True):
                href = link["href"]
                if any(
                    ext in href.lower()
                    for ext in [".pdf", ".doc", ".xls", ".csv", ".zip"]
                ):
                    url = urljoin(base_url, href)
                    filename = Path(urlparse(url).path).name or link.get_text().strip()
                    if filename and not any(f["name"] == filename for f in found_files):
                        found_files.append({"url": url, "name": filename})
                        print(f"      Found supplement: {filename}")

        # 2. Look for explicit supplement links with common Karger patterns
        for link in soup.find_all("a", href=True):
            href = link["href"]
            link_text = link.get_text().lower()

            # Check for supplement-related links
            if "suppl" in href.lower() or "suppl" in link_text:
                url = urljoin(base_url, href)
                filename = Path(urlparse(url).path).name or "supplement"
                if filename and not any(f["name"] == filename for f in found_files):
                    found_files.append({"url": url, "name": filename})
                    print(f"      Found supplement link: {filename}")

        # 3. Try to construct the supplement URL from DOI
        # Karger pattern: karger.com/doi/suppl/10.1159/XXXXXX
        doi_match = re.search(r"10\.1159/\d+", base_url)
        if not doi_match:
            # Try to find DOI in the page content
            doi_match = re.search(r"10\.1159/\d+", html)

        if doi_match:
            doi = doi_match.group(0)
            supp_url = f"https://karger.com/doi/suppl/{doi}"
            print(f"    Constructed supplement URL from DOI: {supp_url}")

            # Try to fetch the supplement page
            try:
                import requests

                headers = {
                    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
                    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
                }
                supp_response = requests.get(
                    supp_url, headers=headers, timeout=30, allow_redirects=True
                )

                if supp_response.status_code == 200:
                    supp_soup = BeautifulSoup(supp_response.text, "html.parser")

                    # Find all downloadable files on the supplement page
                    for link in supp_soup.find_all("a", href=True):
                        href = link["href"]
                        # Look for file downloads
                        if any(
                            ext in href.lower()
                            for ext in [".pdf", ".doc", ".xls", ".csv", ".zip", ".txt"]
                        ):
                            file_url = urljoin(supp_url, href)
                            filename = Path(urlparse(file_url).path).name
                            if filename and not any(
                                f["name"] == filename for f in found_files
                            ):
                                found_files.append({"url": file_url, "name": filename})
                                print(
                                    f"      Found supplement from DOI page: {filename}"
                                )
                else:
                    print(
                        f"    Supplement page returned status {supp_response.status_code}"
                    )
            except Exception as e:
                print(f"    Could not fetch supplement page: {e}")

        if not found_files:
            print(
                "    No Karger supplements found via HTTP. Will try browser fallback if available."
            )
            # Return empty - the orchestrator should try browser-based fetching
            # We don't call generic here because Karger pages are different
            return []

        return found_files

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
        soup = BeautifulSoup(html, "html.parser")
        found_files = []

        # Nature articles often have a "Supplementary Information" section
        supp_info_section = soup.find(id="supplementary-information")
        if not supp_info_section:
            # Fallback to searching for common heading patterns
            supp_info_section = soup.find(
                "h2",
                string=re.compile(r"Supplementary (Information|Data)", re.IGNORECASE),
            )
            if supp_info_section:
                supp_info_section = supp_info_section.parent

        if supp_info_section:
            print("    Found 'Supplementary Information' section.")
            for a in supp_info_section.find_all("a", href=True):
                href = a["href"]
                # Filter for links that look like file downloads
                if (
                    "/articles/" in href
                    and "/figures/" not in href
                    and "author-information" not in href
                ):
                    url = urljoin(base_url, href)
                    filename = Path(urlparse(url).path).name
                    if filename and not any(f["name"] == filename for f in found_files):
                        found_files.append({"url": url, "name": filename})
                        print(f"      Found potential supplement: {filename}")

        if not found_files:
            print(
                "    'Supplementary Information' section not found or empty. Trying generic scan."
            )
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
        soup = BeautifulSoup(html, "html.parser")
        found_files = []

        # 1. Look for data in <script type="application/json"> blocks
        # This is often more reliable than parsing the HTML structure
        for script in soup.find_all("script", type="application/json"):
            try:
                data = json.loads(script.string)
                # Try multiple possible paths in the JSON structure
                supp_data = (
                    data.get("article", {})
                    .get("supplementaryMaterials", {})
                    .get("supplementaryMaterial", [])
                    or data.get("supplementaryMaterials", {}).get(
                        "supplementaryMaterial", []
                    )
                    or data.get("supplementaryMaterial", [])
                    or []
                )
                for item in supp_data:
                    url = item.get("downloadUrl") or item.get("url") or item.get("href")
                    filename = (
                        item.get("title") or item.get("name") or item.get("filename")
                    )
                    if (
                        url
                        and filename
                        and not any(f["name"] == filename for f in found_files)
                    ):
                        found_files.append({"url": url, "name": filename})
                        print(f"    Found supplement in JSON data: {filename}")
                if found_files:
                    return found_files
            except (json.JSONDecodeError, AttributeError, TypeError):
                continue

        # 2. Look for sections with "supplement" or "supporting" in headings
        supp_keywords = ["supplement", "supporting", "additional", "appendix"]
        for heading in soup.find_all(["h1", "h2", "h3", "h4", "h5", "h6"]):
            heading_text = heading.get_text().lower()
            if any(keyword in heading_text for keyword in supp_keywords):
                # Look for links in the section containing this heading
                section = heading.find_next_sibling() or heading.parent
                if section:
                    for link in section.find_all("a", href=True):
                        href = link["href"]
                        link_text = link.get_text().lower()
                        # Check if it's a file link
                        if any(
                            ext in href.lower()
                            for ext in [
                                ".pdf",
                                ".docx",
                                ".xlsx",
                                ".zip",
                                ".csv",
                                ".txt",
                                ".doc",
                            ]
                        ):
                            url = urljoin(base_url, href)
                            filename = (
                                link.get_text().strip() or Path(urlparse(url).path).name
                            )
                            if filename and not any(
                                f["name"] == filename for f in found_files
                            ):
                                found_files.append({"url": url, "name": filename})
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
                if filename and not any(f["name"] == filename for f in found_files):
                    found_files.append({"url": url, "name": filename})
                    print(f"    Found MMC link via regex: {filename}")
        if found_files:
            return found_files

        # 4. Look for links with "supplement" or "mmc" in the text or href
        for link in soup.find_all("a", href=True):
            link_text = link.get_text().lower()
            href = link["href"].lower()
            is_supplement = (
                any(
                    keyword in link_text
                    for keyword in ["supplement", "supporting", "additional", "mmc"]
                )
                or "mmc" in href
                or "supplement" in href
                or "attachment" in href
            )
            is_file = any(
                href.endswith(ext)
                for ext in [".pdf", ".docx", ".xlsx", ".zip", ".csv", ".txt", ".doc"]
            )
            if is_supplement and is_file:
                url = urljoin(base_url, link["href"])
                filename = link.get_text().strip() or Path(urlparse(url).path).name
                if filename and not any(f["name"] == filename for f in found_files):
                    found_files.append({"url": url, "name": filename})
                    print(f"    Found supplement link: {filename}")

        # 5. Look for specific "Download file" links in ScienceDirect as a fallback
        for link in soup.find_all("a", class_="S_C_9cf8451f", href=True):
            if "Download file" in link.get_text():
                url = urljoin(base_url, link["href"])
                filename = link.get_text().replace("Download file", "").strip()
                if filename and not any(f["name"] == filename for f in found_files):
                    found_files.append({"url": url, "name": filename})
                    print(f"    Found ScienceDirect download link: {filename}")

        if found_files:
            return found_files

        print("    No specific Elsevier/GIM supplements found. Trying generic scan.")
        return self.scrape_generic_supplements(html, base_url)

    def scrape_generic_supplements(self, html: str, base_url: str) -> List[Dict]:
        """
        A best-effort generic scraper for supplemental files.

        This is the fallback for publishers without specific handlers.

        CRITICAL LIMITATIONS (major contributor to 94.6% failure rate):

        Current keywords: "supplement", "supporting", "appendix", "additional file"

        MISSES these common patterns:
        - "Electronic Supplementary Material" / "ESM" (Springer)
        - "MOESM" file naming (Springer)
        - /MediaObjects/ URLs (Springer)
        - /downloadSupplement endpoint (Oxford, Wiley)
        - Links with icon-only text (no keyword match)
        - JavaScript-rendered supplement sections
        - Supplements on separate pages requiring navigation

        MATCHING LOGIC:
        - Link text contains keyword, OR
        - Link href ends with known extension (.pdf, .xlsx, etc.)

        This means a link like:
          <a href="/downloadSupplement?file=data.xlsx">Download</a>
        Will NOT match because:
        - "Download" doesn't contain keywords
        - href doesn't END with .xlsx (has query params)

        TODO: Expand keywords and match URL patterns not just extensions.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL for resolving relative links

        Returns:
            List of supplement file dictionaries with 'url', 'name', and optionally
            'base_url' (for PMC pages to enable URL variant generation) keys
        """
        print("  Scraping with scrape_generic_supplements...")
        soup = BeautifulSoup(html, "html.parser")
        found_files = []

        keywords = ["supplement", "supporting", "appendix", "additional file"]
        file_extensions = [
            ".pdf",
            ".docx",
            ".xlsx",
            ".csv",
            ".zip",
            ".rar",
            ".gz",
            ".txt",
            ".doc",
        ]

        # Check if this is a PMC page
        is_pmc_page = (
            "ncbi.nlm.nih.gov/pmc/articles/" in base_url
            or "pmc.ncbi.nlm.nih.gov/articles/" in base_url
        )

        for link in soup.find_all("a", href=True):
            link_text = link.get_text().lower()
            href = link["href"].lower()

            # Check if link text or href contain any of the keywords or file extensions
            is_supplement_link = any(keyword in link_text for keyword in keywords)
            is_file_link = any(href.endswith(ext) for ext in file_extensions)

            if is_supplement_link or is_file_link:
                try:
                    original_url = urljoin(base_url, link["href"])
                    # Normalize PMC URLs to use correct path format
                    url = self._normalize_pmc_url(original_url, base_url)
                    filename = Path(urlparse(url).path).name

                    # Skip invalid filenames
                    if not filename:
                        continue

                    # Skip entries that look like PMCIDs (e.g., "PMC3049907")
                    if re.match(r"^PMC\d+$", filename, re.IGNORECASE):
                        continue

                    # Skip entries that don't have a file extension and aren't actual files
                    has_extension = any(
                        filename.lower().endswith(ext) for ext in file_extensions
                    )
                    if not has_extension and is_supplement_link:
                        # Only skip non-extension links if they also look like IDs
                        if re.match(r"^[A-Z]+\d+$", filename, re.IGNORECASE):
                            continue

                    # Basic filtering to avoid irrelevant links
                    if not any(f["name"] == filename for f in found_files):
                        file_info = {"url": url, "name": filename}
                        # Store original URL and base URL for PMC pages to enable URL variant generation
                        if is_pmc_page:
                            file_info["original_url"] = original_url
                            file_info["base_url"] = base_url
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

    def extract_fulltext_nature(
        self, html: str, base_url: str
    ) -> Tuple[Optional[str], str]:
        """
        Extract full-text content from a Nature journal page.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL of the article

        Returns:
            Tuple of (markdown_content, title)
        """
        print("  Extracting full text from Nature article...")
        soup = BeautifulSoup(html, "html.parser")
        markdown = "# MAIN TEXT\n\n"

        # Extract title
        title = ""
        title_elem = soup.find("h1", class_="c-article-title") or soup.find("h1")
        if title_elem:
            title = title_elem.get_text().strip()
            markdown += f"## {title}\n\n"

        # Extract abstract
        abstract = soup.find("div", id="Abs1-content") or soup.find(
            "section", {"data-title": "Abstract"}
        )
        if abstract:
            markdown += "### Abstract\n\n"
            markdown += abstract.get_text().strip() + "\n\n"

        # Extract main article content
        article_body = soup.find("div", class_="c-article-body") or soup.find("main")
        if article_body:
            # Find all sections
            sections = article_body.find_all(
                ["section", "div"], class_=re.compile(r"c-article-section")
            )
            if not sections:
                sections = article_body.find_all("section")

            for section in sections:
                # Get section title
                section_title = section.find(["h2", "h3", "h4"])
                if section_title:
                    markdown += f"### {section_title.get_text().strip()}\n\n"

                # Get paragraphs
                for p in section.find_all("p", recursive=False):
                    text = p.get_text().strip()
                    if text:
                        markdown += f"{text}\n\n"

        # Fallback: just extract all paragraph text
        if len(markdown) < 500:
            markdown = "# MAIN TEXT\n\n"
            if title:
                markdown += f"## {title}\n\n"

            for p in soup.find_all("p"):
                text = p.get_text().strip()
                if len(text) > 50:  # Skip very short paragraphs (likely nav elements)
                    markdown += f"{text}\n\n"

        return markdown if len(markdown) > 200 else None, title

    def extract_fulltext_elsevier(
        self, html: str, base_url: str
    ) -> Tuple[Optional[str], str]:
        """
        Extract full-text content from an Elsevier/ScienceDirect page.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL of the article

        Returns:
            Tuple of (markdown_content, title)
        """
        print("  Extracting full text from Elsevier/ScienceDirect article...")
        soup = BeautifulSoup(html, "html.parser")
        markdown = "# MAIN TEXT\n\n"

        # Extract title
        title = ""
        title_elem = soup.find("span", class_="title-text") or soup.find(
            "h1", class_="svTitle"
        )
        if title_elem:
            title = title_elem.get_text().strip()
            markdown += f"## {title}\n\n"

        # Try to extract from JSON data embedded in page (ScienceDirect often has this)
        for script in soup.find_all("script", type="application/json"):
            try:
                data = json.loads(script.string)
                # Look for article content in various possible JSON paths
                article_data = data.get("article", {})
                if article_data:
                    # Extract abstract
                    abstract_data = article_data.get("abstract", {})
                    if abstract_data:
                        markdown += "### Abstract\n\n"
                        if isinstance(abstract_data, dict):
                            abstract_text = abstract_data.get("content", "")
                        else:
                            abstract_text = str(abstract_data)
                        # Clean HTML tags from abstract
                        abstract_soup = BeautifulSoup(abstract_text, "html.parser")
                        markdown += abstract_soup.get_text().strip() + "\n\n"

                    # Extract body sections
                    body_data = article_data.get("body", {})
                    if body_data and isinstance(body_data, dict):
                        content = body_data.get("content", [])
                        if isinstance(content, list):
                            for section in content:
                                if isinstance(section, dict):
                                    sec_title = section.get("label", "") or section.get(
                                        "title", ""
                                    )
                                    if sec_title:
                                        markdown += f"### {sec_title}\n\n"
                                    sec_content = section.get("content", "")
                                    if sec_content:
                                        content_soup = BeautifulSoup(
                                            sec_content, "html.parser"
                                        )
                                        markdown += (
                                            content_soup.get_text().strip() + "\n\n"
                                        )
                if len(markdown) > 500:
                    return markdown, title
            except (json.JSONDecodeError, AttributeError, TypeError):
                continue

        # Fallback: HTML-based extraction
        # Extract abstract
        abstract = soup.find("div", class_="abstract") or soup.find(
            "div", id="abstracts"
        )
        if abstract:
            markdown += "### Abstract\n\n"
            markdown += abstract.get_text().strip() + "\n\n"

        # Extract main body content
        body = soup.find("div", id="body") or soup.find("div", class_="Body")
        if body:
            for section in body.find_all(
                ["section", "div"], class_=re.compile(r"section")
            ):
                section_title = section.find(["h2", "h3", "h4"])
                if section_title:
                    markdown += f"### {section_title.get_text().strip()}\n\n"

                for p in section.find_all("p"):
                    text = p.get_text().strip()
                    if text:
                        markdown += f"{text}\n\n"

        # Ultimate fallback: more aggressive extraction
        if len(markdown) < 500:
            markdown = "# MAIN TEXT\n\n"
            if title:
                markdown += f"## {title}\n\n"

            # Try multiple content containers
            content_divs = soup.find_all(
                "div",
                class_=re.compile(
                    r"(Body|content|article|fulltext|full-text|text)", re.I
                ),
            )
            for div in content_divs:
                for p in div.find_all("p"):
                    text = p.get_text().strip()
                    if len(text) > 50:
                        markdown += f"{text}\n\n"

        # If still too short, extract all paragraphs from body
        if len(markdown) < 500:
            body = soup.find("body")
            if body:
                for p in body.find_all("p"):
                    text = p.get_text().strip()
                    if len(text) > 50:
                        # Skip common non-content patterns
                        skip_patterns = [
                            "cookie",
                            "privacy",
                            "sign in",
                            "subscribe",
                            "javascript",
                        ]
                        if not any(
                            pattern in text.lower() for pattern in skip_patterns
                        ):
                            markdown += f"{text}\n\n"

        return markdown if len(markdown) > 200 else None, title

    def extract_fulltext_wiley(
        self, html: str, base_url: str
    ) -> Tuple[Optional[str], str]:
        """
        Extract full-text content from a Wiley Online Library page.

        Args:
            html: HTML content of the publisher page
            base_url: Base URL of the article

        Returns:
            Tuple of (markdown_content, title)
        """
        print("  Extracting full text from Wiley article...")
        soup = BeautifulSoup(html, "html.parser")
        markdown = "# MAIN TEXT\n\n"

        # Extract title
        title = ""
        title_elem = soup.find("h1", class_="citation__title") or soup.find("h1")
        if title_elem:
            title = title_elem.get_text().strip()
            markdown += f"## {title}\n\n"

        # Extract abstract
        abstract = soup.find(
            "section", class_="article-section__abstract"
        ) or soup.find("div", class_="abstract")
        if abstract:
            markdown += "### Abstract\n\n"
            markdown += abstract.get_text().strip() + "\n\n"

        # Extract main body
        body = soup.find("section", class_="article-section__content") or soup.find(
            "div", class_="article__body"
        )
        if body:
            for section in body.find_all(["section", "div"]):
                section_title = section.find(["h2", "h3", "h4"])
                if section_title:
                    markdown += f"### {section_title.get_text().strip()}\n\n"

                for p in section.find_all("p", recursive=False):
                    text = p.get_text().strip()
                    if text:
                        markdown += f"{text}\n\n"

        # Fallback
        if len(markdown) < 500:
            for p in soup.find_all("p"):
                text = p.get_text().strip()
                if len(text) > 50:
                    markdown += f"{text}\n\n"

        return markdown if len(markdown) > 200 else None, title

    def extract_fulltext_generic(
        self, html: str, base_url: str
    ) -> Tuple[Optional[str], str]:
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
        soup = BeautifulSoup(html, "html.parser")
        markdown = "# MAIN TEXT\n\n"

        # Remove script, style, nav, footer, header elements
        for element in soup.find_all(
            ["script", "style", "nav", "footer", "header", "aside", "noscript"]
        ):
            element.decompose()

        # Extract title
        title = ""
        title_elem = soup.find("h1") or soup.find("title")
        if title_elem:
            title = title_elem.get_text().strip()
            # Clean up title (remove site name etc)
            title = re.sub(r"\s*[-|].*$", "", title)
            markdown += f"## {title}\n\n"

        # Look for common article containers - expanded list for better coverage
        article_containers = [
            soup.find("article"),
            soup.find(
                "div",
                class_=re.compile(
                    r"(article|content|body|main|fulltext|full-text)", re.I
                ),
            ),
            soup.find("main"),
            soup.find(
                "div",
                id=re.compile(r"(article|content|body|main|fulltext|full-text)", re.I),
            ),
            # Publisher-specific containers
            soup.find(
                "div",
                class_=re.compile(
                    r"(hlFld-Fulltext|article-content|ArticleBody)", re.I
                ),
            ),
            soup.find("section", class_=re.compile(r"(article|content|body)", re.I)),
            # AHA journals specific
            soup.find("div", class_="article__body"),
            soup.find("div", class_="article-full-text"),
            # Oxford Academic specific
            soup.find("div", class_="article-body"),
            # Springer specific
            soup.find("div", class_="c-article-body"),
            # PNAS specific
            soup.find("div", class_="article fulltext-view"),
        ]

        content_found = False
        extracted_text = []

        for container in article_containers:
            if container:
                # Try to extract from sections first
                sections = container.find_all(["section", "div"], recursive=True)
                for section in sections:
                    section_title = section.find(["h2", "h3", "h4"], recursive=False)
                    if section_title:
                        title_text = section_title.get_text().strip()
                        if (
                            title_text and len(title_text) < 100
                        ):  # Skip overly long "titles"
                            extracted_text.append(f"### {title_text}\n\n")

                    for p in section.find_all("p", recursive=False):
                        text = p.get_text().strip()
                        if len(text) > 30:
                            extracted_text.append(f"{text}\n\n")
                            content_found = True

                # If no sections, try getting paragraphs directly
                if not content_found:
                    for p in container.find_all("p"):
                        text = p.get_text().strip()
                        if len(text) > 30:
                            extracted_text.append(f"{text}\n\n")
                            content_found = True

                if content_found:
                    break

        # If containers didn't work, try extracting from the whole page more aggressively
        if not content_found or len("".join(extracted_text)) < 500:
            extracted_text = []

            # Look for abstract first
            abstract = soup.find(
                ["div", "section"], class_=re.compile(r"abstract", re.I)
            )
            if abstract:
                abstract_text = abstract.get_text().strip()
                if len(abstract_text) > 100:
                    extracted_text.append(f"### Abstract\n\n{abstract_text}\n\n")

            # Get all paragraphs from body
            body = soup.find("body")
            if body:
                for p in body.find_all("p"):
                    text = p.get_text().strip()
                    # Filter out short paragraphs that are likely navigation or UI elements
                    # But use a lower threshold to capture more content
                    if len(text) > 50:
                        # Skip cookie notices, login prompts, etc.
                        skip_patterns = [
                            "cookie",
                            "privacy policy",
                            "sign in",
                            "log in",
                            "subscribe",
                            "register",
                            "your browser",
                            "javascript",
                        ]
                        if not any(
                            pattern in text.lower() for pattern in skip_patterns
                        ):
                            extracted_text.append(f"{text}\n\n")
                            content_found = True

        if extracted_text:
            markdown += "".join(extracted_text)

        # Final fallback: just get all text if nothing worked
        if len(markdown) < 500:
            body = soup.find("body")
            if body:
                all_text = body.get_text(separator="\n", strip=True)
                # Clean up excessive whitespace
                all_text = re.sub(r"\n{3,}", "\n\n", all_text)
                if len(all_text) > 500:
                    markdown = f"# MAIN TEXT\n\n## {title}\n\n{all_text}\n"

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

        if "nature.com" in domain:
            return self.extract_fulltext_nature(html, base_url)
        elif any(
            d in domain for d in ["sciencedirect.com", "elsevier.com", "gimjournal.org"]
        ):
            return self.extract_fulltext_elsevier(html, base_url)
        elif "wiley.com" in domain or "onlinelibrary.wiley.com" in domain:
            return self.extract_fulltext_wiley(html, base_url)
        else:
            return self.extract_fulltext_generic(html, base_url)


# =============================================================================
# TODO: REMAINING IMPLEMENTATIONS TO IMPROVE FAILURE RATE
# =============================================================================
#
# PRIORITY 1: ✅ DONE - scrape_springer_supplements() implemented (2026-02-01)
# See the actual implementation above in the SupplementScraper class.
# Handles: link.springer.com, biomedcentral.com, springeropen.com
# Tested with: 10.1186/s12864-019-6413-7 (found 5 MOESM PDFs)
#
# PRIORITY 2: Add scrape_oxford_supplements()
# -------------------------------------------
# Oxford Academic (academic.oup.com) is critical for genetics journals.
#
# def scrape_oxford_supplements(self, html: str, base_url: str) -> List[Dict]:
#     """
#     Scrape supplements from Oxford Academic (academic.oup.com)
#
#     Key patterns:
#     - Section: "Supplementary data" tab/section
#     - Link class: "supplementary-data" or similar
#     - URL patterns:
#       * /downloadSupplement?file=XXXXX_supplementary_data.xlsx
#       * /article/XXXXX/XXXXX/supplementary-data
#
#     CHALLENGE: Often requires JavaScript to render supplement tab.
#     May need to construct supplement URL from article URL pattern.
#     """
#     soup = BeautifulSoup(html, "html.parser")
#     found_files = []
#
#     # 1. Look for direct download links
#     for link in soup.find_all('a', href=re.compile(r'downloadSupplement')):
#         href = link['href']
#         url = urljoin(base_url, href)
#         filename = link.get_text().strip() or 'supplement'
#         found_files.append({'url': url, 'name': filename})
#
#     # 2. Look for supplementary-data section
#     supp_section = soup.find(id='supplementary-data') or \
#                    soup.find('section', class_=re.compile(r'supplementary'))
#
#     # 3. Try constructing supplement page URL from article URL
#     # Oxford pattern: /article/doi/xxxxx -> /article/doi/xxxxx/supplementary-data
#
#     return found_files if found_files else self.scrape_generic_supplements(html, base_url)
#
#
# PRIORITY 3: Improve scrape_generic_supplements()
# ------------------------------------------------
# Current keywords miss too much. Add:
#   - Keywords: "ESM", "MOESM", "electronic supplementary", "supporting information"
#   - URL patterns: /downloadSupplement, /MediaObjects/, /figures/, /suppl/
#   - Don't require file extension if URL contains clear supplement indicators
#
#
# PRIORITY 4: Add routing in doi_resolver.py
# ------------------------------------------
# Once handlers are implemented, add to resolve_and_scrape_supplements():
#
#     elif "springer.com" in domain or "biomedcentral.com" in domain:
#         return scraper.scrape_springer_supplements(response.text, final_url)
#     elif "academic.oup.com" in domain or "oup.com" in domain:
#         return scraper.scrape_oxford_supplements(response.text, final_url)
#
#
# PRIORITY 5: Browser fallback integration
# ----------------------------------------
# browser_supplement_fetcher.py has Karger implementation.
# Need to:
# 1. Integrate browser fallback into main pipeline (when HTTP fails)
# 2. Add browser handlers for Oxford, Springer (similar to Karger)
# 3. Consider making browser the default for known-problematic publishers
#
# =============================================================================
