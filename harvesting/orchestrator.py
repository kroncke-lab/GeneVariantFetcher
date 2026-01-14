"""
Orchestrator Module

Main PMCHarvester class that coordinates all harvesting operations:
- Converts PMIDs to PMCIDs
- Downloads full-text and supplemental files
- Creates unified markdown files for LLM processing
"""

import os
import re
import time
import csv
import requests
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import urlparse

from .pmc_api import PMCAPIClient
from .doi_resolver import DOIResolver
from .supplement_scraper import SupplementScraper
from .format_converters import FormatConverter
from .elsevier_api import ElsevierAPIClient
from .wiley_api import WileyAPIClient

# Import scout components (with fallback for import errors)
try:
    from pipeline.data_scout import GeneticDataScout
    from config.settings import get_settings
    SCOUT_AVAILABLE = True
except ImportError:
    SCOUT_AVAILABLE = False
    get_settings = None


class PMCHarvester:
    """Harvests full-text and supplemental materials from PubMed Central."""

    SUSPICIOUS_FREE_URL_DOMAINS = {"antibodies.cancer.gov"}

    def __init__(self, output_dir: str = "pmc_harvest", gene_symbol: str = None):
        """
        Initialize PMC Harvester.

        Args:
            output_dir: Directory to save harvested files
            gene_symbol: Target gene symbol for data scout analysis
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.paywalled_log = self.output_dir / "paywalled_missing.csv"
        self.success_log = self.output_dir / "successful_downloads.csv"
        self.gene_symbol = gene_symbol

        # Initialize session with browser-like headers
        self.session = requests.Session()
        HEADERS = {
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
            'Accept-Language': 'en-US,en;q=0.9',
            'Accept-Encoding': 'gzip, deflate, br',
            'Referer': 'https://pubmed.ncbi.nlm.nih.gov/',
            'DNT': '1',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1'
        }
        self.session.headers.update(HEADERS)

        # Initialize component modules
        self.pmc_api = PMCAPIClient(session=self.session)
        self.doi_resolver = DOIResolver(session=self.session, paywalled_log=self.paywalled_log)
        self.scraper = SupplementScraper()
        self.converter = FormatConverter()

        # Initialize Elsevier API client (optional - uses API key from settings if available)
        elsevier_api_key = None
        elsevier_insttoken = None
        wiley_api_key = None
        if get_settings is not None:
            try:
                settings = get_settings()
                elsevier_api_key = settings.elsevier_api_key
                elsevier_insttoken = settings.elsevier_insttoken
                wiley_api_key = settings.wiley_api_key
            except Exception:
                pass  # Settings validation may fail if other keys are missing
        self.elsevier_api = ElsevierAPIClient(api_key=elsevier_api_key, insttoken=elsevier_insttoken, session=self.session)
        self.wiley_api = WileyAPIClient(api_key=wiley_api_key, session=self.session)

        # Initialize log files
        with open(self.paywalled_log, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['PMID', 'Reason', 'URL'])

        with open(self.success_log, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['PMID', 'PMCID', 'Supplements_Downloaded'])

    def _run_data_scout(self, pmid: str, unified_content: str) -> bool:
        """
        Run the Genetic Data Scout to identify high-value data zones.

        Creates two additional files:
        - {PMID}_DATA_ZONES.json: Zone metadata for debugging/analysis
        - {PMID}_DATA_ZONES.md: Condensed markdown with only high-value zones

        Args:
            pmid: PubMed ID
            unified_content: Full markdown content to analyze

        Returns:
            True if scout ran successfully, False otherwise
        """
        if not SCOUT_AVAILABLE:
            return False

        if not self.gene_symbol:
            print(f"  - Skipping data scout: no gene symbol provided")
            return False

        try:
            settings = get_settings()
            if not settings.scout_enabled:
                return False

            scout = GeneticDataScout(
                gene_symbol=self.gene_symbol,
                min_relevance_score=settings.scout_min_relevance,
                max_zones=settings.scout_max_zones,
            )

            report = scout.scan(unified_content, pmid=pmid)

            # Write zone metadata JSON
            zones_json_path = self.output_dir / f"{pmid}_DATA_ZONES.json"
            zones_json_path.write_text(scout.to_json(report))

            # Write condensed markdown with high-value zones only
            zones_md_path = self.output_dir / f"{pmid}_DATA_ZONES.md"
            zones_md_path.write_text(scout.format_markdown(report, unified_content))

            print(f"  ✓ Data Scout: {report.zones_kept}/{report.total_zones_found} zones kept ({report.compression_ratio:.0%} of original)")
            return True

        except Exception as e:
            print(f"  - Data scout failed: {e}")
            return False

    @staticmethod
    def _markdown_missing_body(markdown: str) -> bool:
        """
        Detect whether the converted markdown is missing main body text.

        Returns True when only the abstract is present or the content
        is suspiciously short, signalling that we should try an HTML scrape.
        """
        if not markdown:
            return True
        headings = [line for line in markdown.splitlines() if line.startswith("### ")]
        non_abstract = [h for h in headings if not h.startswith("### Abstract")]
        return (len(non_abstract) == 0) or (len(markdown) < 2000)

    def get_supplemental_files(self, pmcid: str, pmid: str, doi: str) -> List[Dict[str, Any]]:
        """
        Orchestrates fetching supplemental files, first via API, then by scraping.

        Args:
            pmcid: PubMed Central ID
            pmid: PubMed ID (for logging)
            doi: Digital Object Identifier (for fallback scraping)

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys,
            plus 'base_url' and 'original_url' for PMC supplements
        """
        pmc_base_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/" if pmcid else None

        # 1. First, try the EuropePMC API
        api_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
        try:
            response = self.session.get(api_url, timeout=30)
            response.raise_for_status()

            # Check if response body is not empty
            if not response.text.strip():
                raise ValueError("Empty response from API")

            content_type = response.headers.get('Content-Type', '').lower()

            # Try to parse as JSON first
            if 'application/json' in content_type:
                data = response.json()
                # The API can return a success response with an empty list of files
                if data.get('result', {}).get('supplementaryFiles'):
                    print("  ✓ Found supplemental files via EuropePMC API (JSON)")
                    supp_files = data['result']['supplementaryFiles']
                    # Add PMC base_url to enable URL variant generation
                    if pmc_base_url:
                        for f in supp_files:
                            f['base_url'] = pmc_base_url
                            f['original_url'] = f.get('url', '')
                    return supp_files
            # Try to parse as XML if JSON fails
            elif 'application/xml' in content_type or 'text/xml' in content_type:
                try:
                    from xml.etree import ElementTree as ET
                    root = ET.fromstring(response.text)
                    # Look for supplementary file elements in the XML
                    supp_files = []
                    # Common XML paths for supplementary files
                    for elem in root.iter():
                        if 'supplement' in elem.tag.lower() or 'supplementary' in elem.tag.lower():
                            url = elem.get('url') or elem.get('href') or elem.text
                            name = elem.get('name') or elem.get('title') or elem.get('filename')
                            if url and name:
                                file_info = {'url': url, 'name': name}
                                if pmc_base_url:
                                    file_info['base_url'] = pmc_base_url
                                    file_info['original_url'] = url
                                supp_files.append(file_info)
                    if supp_files:
                        print(f"  ✓ Found {len(supp_files)} supplemental files via EuropePMC API (XML)")
                        return supp_files
                    # If XML parsing succeeded but no files found, continue to fallback
                except ET.ParseError as e:
                    raise ValueError(f"Could not parse XML response: {e}")
                # If we get here, XML was parsed but no files found - continue to fallback
            else:
                raise ValueError(f"Expected JSON or XML but got Content-Type: {content_type}")
        except requests.exceptions.RequestException as e:
            print(f"  - EuropePMC supplemental files API request failed for {pmcid}: {e}")
        except ValueError as e:
            print(f"  - EuropePMC supplemental files API returned invalid response for {pmcid}: {e}")
        except Exception as e:
            print(f"  - EuropePMC supplemental files API failed for {pmcid}: {e}")

        # 2. If API fails, try scraping PMC page directly
        if pmcid:
            try:
                pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                print(f"  - Trying to scrape PMC page directly: {pmc_url}")
                pmc_response = self.session.get(pmc_url, timeout=30, allow_redirects=True)
                pmc_response.raise_for_status()

                # Use the final URL after redirects as the base URL
                final_pmc_url = pmc_response.url

                # Use generic scraper for PMC pages
                pmc_files = self.scraper.scrape_generic_supplements(pmc_response.text, final_pmc_url)
                if pmc_files:
                    print(f"  ✓ Found {len(pmc_files)} supplemental files via PMC page scraping")
                    return pmc_files
            except Exception as e:
                print(f"  - PMC page scraping failed: {e}")

        # 3. If PMC scraping fails, fall back to DOI-based scraping
        if doi:
            print(f"  - Falling back to DOI scraping for {doi}")
            return self.doi_resolver.resolve_and_scrape_supplements(doi, pmid, self.scraper)
        else:
            print("  - API failed and no DOI available. Cannot scrape.")
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/" if pmcid else "N/A"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Supplemental files API failed, no DOI', pmc_url])
            return []

    def download_supplement(self, url: str, output_path: Path, pmid: str, filename: str,
                            base_url: str = None, original_url: str = None) -> bool:
        """
        Download a supplemental file, trying multiple URL variants for PMC supplements.

        Args:
            url: URL of the supplement file
            output_path: Path to save the file
            pmid: PubMed ID (for logging)
            filename: Name of the file (for logging)
            base_url: Base URL of the article page (for PMC URL variant generation)
            original_url: Original URL before normalization (for PMC URL variant generation)

        Returns:
            True if download succeeded, False otherwise
        """
        # Generate URL variants for PMC supplements
        urls_to_try = [url]
        if base_url and ('ncbi.nlm.nih.gov/pmc/' in base_url or 'pmc.ncbi.nlm.nih.gov/' in base_url):
            # Use the scraper's method to generate URL variants
            source_url = original_url if original_url else url
            variants = self.scraper.get_pmc_supplement_url_variants(source_url, base_url)
            # Put primary URL first, then add unique variants
            urls_to_try = [url] + [v for v in variants if v != url]

        last_error = None
        for try_url in urls_to_try:
            try:
                # Use the session to download, which includes our headers
                response = self.session.get(try_url, timeout=60, stream=True, allow_redirects=True)
                response.raise_for_status()

                # Collect content to validate before writing
                content_chunks = []
                for chunk in response.iter_content(chunk_size=8192):
                    content_chunks.append(chunk)

                content = b''.join(content_chunks)

                # Validate the downloaded content
                ext = output_path.suffix.lower()

                # Check for HTML error pages disguised as other file types
                # HTML pages typically start with <!DOCTYPE, <html, or have these early in the content
                content_start = content[:4096].lower() if len(content) >= 4096 else content.lower()
                is_html_page = (
                    content_start.startswith(b'<!doctype') or
                    content_start.startswith(b'<html') or
                    b'<!doctype html' in content_start or
                    b'<html' in content_start[:500] or
                    # Some servers return XHTML
                    b'<?xml' in content_start[:100] and b'<html' in content_start[:500]
                )

                # Detect specific access denial patterns in HTML
                access_denied_patterns = [
                    # HTTP error codes
                    b'403 forbidden', b'401 unauthorized', b'access denied',
                    # Paywall messages
                    b'subscription required', b'purchase this article',
                    b'institutional access', b'sign in to access',
                    b'log in to access', b'login required',
                    # Content not available
                    b'pdf not available', b'full text not available',
                    b'content not available', b'article not found',
                    # ScienceDirect specific
                    b'sciencedirect', b'elsevier', b'get access',
                    # Wiley specific
                    b'wiley online library',
                    # Generic paywall indicators
                    b'buy this article', b'rent this article',
                    b'get full access', b'view full text',
                ]
                is_access_denied = is_html_page and any(
                    pattern in content_start for pattern in access_denied_patterns
                )

                # Check for specific publisher error page signatures
                is_publisher_error = is_html_page and (
                    b'error' in content_start[:200] or
                    b'not found' in content_start[:500] or
                    b'unavailable' in content_start[:500]
                )

                # PDF validation: PDFs should start with %PDF
                if ext == '.pdf':
                    if not content.startswith(b'%PDF'):
                        if is_access_denied:
                            last_error = f"{filename}: ACCESS DENIED (paywall/login required)"
                            print(f"    ⚠ {filename}: Access denied - likely paywall or login required")
                            continue  # Try next URL variant
                        elif is_publisher_error:
                            last_error = f"{filename}: Publisher error page received instead of PDF"
                            print(f"    ⚠ {filename}: Received publisher error page instead of PDF")
                            continue  # Try next URL variant
                        elif is_html_page:
                            last_error = f"{filename}: HTML page received instead of PDF (file may be paywalled)"
                            print(f"    ⚠ {filename}: Received HTML page instead of PDF - file may be paywalled or URL expired")
                            continue  # Try next URL variant
                        else:
                            # Content doesn't start with %PDF but isn't HTML either
                            # Could be corrupted or wrong content type
                            content_preview = content[:50].hex() if len(content) > 0 else "(empty)"
                            last_error = f"{filename}: Invalid PDF (starts with: {content_preview[:20]}...)"
                            print(f"    ⚠ {filename}: Invalid PDF file (wrong magic bytes)")
                            continue  # Try next URL variant

                # For other file types, check if we got an HTML error page
                elif ext in ['.docx', '.xlsx', '.xls', '.doc', '.zip']:
                    if is_access_denied:
                        last_error = f"{filename}: ACCESS DENIED - file appears paywalled"
                        print(f"    ⚠ {filename}: Access denied - file may be paywalled")
                        continue  # Try next URL variant
                    elif is_html_page:
                        last_error = f"{filename}: HTML page received instead of {ext} file"
                        print(f"    ⚠ {filename}: Received HTML page instead of {ext} file - may be paywalled")
                        continue  # Try next URL variant

                # Write the validated content
                with open(output_path, 'wb') as f:
                    f.write(content)

                if try_url != url:
                    print(f"    ✓ Downloaded from alternate URL: {try_url}")
                return True

            except Exception as e:
                last_error = str(e)
                continue  # Try next URL variant

        # All URL variants failed
        print(f"    Error downloading {filename}: {last_error}")
        # Log failed supplemental download
        with open(self.paywalled_log, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([pmid, f'Supplemental file download failed: {filename}', url])
        return False

    def _process_supplements(self, pmid: str, pmcid: str, doi: str) -> Tuple[str, int]:
        """
        Download supplemental files and convert them to markdown.

        Args:
            pmid: PubMed ID
            pmcid: PubMed Central ID
            doi: Digital Object Identifier for fallback scraping

        Returns:
            Tuple containing the combined supplement markdown and number of downloaded files
        """
        supplements_dir = self.output_dir / f"{pmid}_supplements"
        supplements_dir.mkdir(exist_ok=True)

        supp_files = self.get_supplemental_files(pmcid, pmid, doi)
        print(f"  Found {len(supp_files)} supplemental files")

        supplement_markdown = ""
        downloaded_count = 0

        for idx, supp in enumerate(supp_files, 1):
            url = supp.get('url', '')
            filename = supp.get('name', f'supplement_{idx}')
            # Get PMC-specific URL info for trying multiple variants
            base_url = supp.get('base_url')
            original_url = supp.get('original_url')

            if not url:
                continue

            file_path = supplements_dir / filename
            print(f"    Downloading: {filename}")

            if self.download_supplement(url, file_path, pmid, filename, base_url, original_url):
                downloaded_count += 1

                ext = file_path.suffix.lower()
                supplement_markdown += f"\n\n# SUPPLEMENTAL FILE {idx}: {filename}\n\n"

                if ext in ['.xlsx', '.xls']:
                    supplement_markdown += self.converter.excel_to_markdown(file_path)
                elif ext == '.docx':
                    supplement_markdown += self.converter.docx_to_markdown(file_path)
                elif ext == '.doc':
                    supplement_markdown += self.converter.doc_to_markdown(file_path)
                elif ext == '.pdf':
                    supplement_markdown += self.converter.pdf_to_markdown(file_path)
                elif ext in ['.txt', '.csv']:
                    try:
                        text = file_path.read_text(encoding='utf-8', errors='ignore')
                        supplement_markdown += text + "\n\n"
                    except Exception as e:
                        supplement_markdown += f"[Error reading text file: {e}]\n\n"
                else:
                    supplement_markdown += f"[File available at: {file_path}]\n\n"

            time.sleep(0.5)

        return supplement_markdown, downloaded_count

    def process_pmid(self, pmid: str) -> Tuple[bool, str]:
        """
        Process a single PMID: convert to PMCID, download content, create unified markdown.

        For papers without PMCIDs but marked as "Free Full Text" on PubMed, this method
        will attempt to fetch full text directly from the publisher's website.

        Args:
            pmid: PubMed ID to process

        Returns:
            Tuple of (success: bool, result: str) where result is output file path or error message
        """
        print(f"\nProcessing PMID: {pmid}")

        # Get DOI first (needed for supplement fallback and free text retrieval)
        doi = self.pmc_api.get_doi_from_pmid(pmid)
        if doi:
            print(f"  ✓ DOI: {doi}")
        else:
            print("  - No DOI found for this PMID.")

        # Convert PMID to PMCID
        pmcid = self.pmc_api.pmid_to_pmcid(pmid)

        if not pmcid:
            # No PMCID - check if this is a free full text article via publisher
            print(f"  - No PMCID found, checking for free full text via publisher...")
            return self._process_free_text_pmid(pmid, doi)

        print(f"  ✓ PMCID: {pmcid}")

        # Get full-text XML from PMC
        xml_content = self.pmc_api.get_fulltext_xml(pmcid)

        if not xml_content:
            print(f"  ❌ Full-text not available from PMC")
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Full-text not available', pmc_url])
            return False, "No full-text"

        print(f"  ✓ Full-text XML retrieved from PMC")

        # Convert main text to markdown
        main_markdown = self.converter.xml_to_markdown(xml_content)

        # Some PMC entries block XML full text (only abstract). Fall back to scraping
        # the PMC HTML page when the markdown looks incomplete.
        if self._markdown_missing_body(main_markdown):
            pmc_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
            try:
                html_response = self.session.get(pmc_url, timeout=30, allow_redirects=True)
                html_response.raise_for_status()
                html_markdown = self.converter.pmc_html_to_markdown(html_response.text)
                if html_markdown and len(html_markdown) > len(main_markdown):
                    main_markdown = html_markdown
                    print("  ✓ Recovered full text from PMC HTML (XML body unavailable)")
            except Exception as e:
                print(f"  - PMC HTML fallback failed: {e}")

        supplement_markdown, downloaded_count = self._process_supplements(pmid, pmcid, doi)

        # Create unified markdown file
        unified_content = main_markdown + supplement_markdown

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(unified_content)

        print(f"  ✅ Created: {output_file.name} ({downloaded_count} supplements)")

        # Run data scout to create condensed DATA_ZONES.md
        self._run_data_scout(pmid, unified_content)

        # Log success
        with open(self.success_log, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([pmid, pmcid, downloaded_count])

        return True, str(output_file)

    def _try_elsevier_api(self, doi: str, pmid: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Try to fetch full text via Elsevier API if applicable.

        Args:
            doi: Digital Object Identifier
            pmid: PubMed ID (for logging)

        Returns:
            Tuple of (markdown_content, error_message)
            - markdown_content: Full text as markdown if successful, None otherwise
            - error_message: Error description if failed, None if successful
        """
        if not self.elsevier_api.is_available:
            return None, "Elsevier API key not configured"

        if not doi or not self.elsevier_api.is_elsevier_doi(doi):
            return None, "Not an Elsevier DOI"

        print(f"  Trying Elsevier API for DOI: {doi}")
        markdown, error = self.elsevier_api.fetch_fulltext(doi=doi)

        if markdown:
            print(f"  ✓ Full text retrieved via Elsevier API ({len(markdown)} characters)")
            return markdown, None
        else:
            print(f"  - Elsevier API: {error}")
            return None, error

    def _try_wiley_api(self, doi: str, pmid: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Try to fetch full text via Wiley TDM API if applicable.

        Args:
            doi: Digital Object Identifier
            pmid: PubMed ID (for logging)

        Returns:
            Tuple of (markdown_content, error_message)
            - markdown_content: Full text as markdown if successful, None otherwise
            - error_message: Error description if failed, None if successful
        """
        if not self.wiley_api.is_available:
            return None, "Wiley API key not configured"

        if not doi or not self.wiley_api.is_wiley_doi(doi):
            return None, "Not a Wiley DOI"

        print(f"  Trying Wiley API for DOI: {doi}")
        markdown, error = self.wiley_api.fetch_fulltext(doi=doi)

        if markdown:
            print(f"  ✓ Full text retrieved via Wiley API ({len(markdown)} characters)")
            return markdown, None
        else:
            print(f"  - Wiley API: {error}")
            return None, error

    def _process_free_text_pmid(self, pmid: str, doi: str) -> Tuple[bool, str]:
        """
        Process a PMID that has no PMCID but may have free full text via publisher.

        This method checks if the article is marked as "Free Full Text" on PubMed
        and attempts to fetch the content from the publisher's website.
        For Elsevier articles, it first tries the Elsevier API if an API key is configured.

        Args:
            pmid: PubMed ID to process
            doi: DOI for the article (may be None)

        Returns:
            Tuple of (success: bool, result: str) where result is output file path or error message
        """
        # Check if article is marked as free full text
        is_free, free_url = self.pmc_api.is_free_full_text(pmid)

        if not is_free:
            # Not a free article - log and skip
            print(f"  ❌ No PMCID and not marked as free full text (likely paywalled)")
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'No PMCID found, not free full text', pubmed_url])
            return False, "No PMCID"

        print(f"  ✓ Article marked as free full text on PubMed")

        # Initialize variables for tracking content source
        main_markdown = None
        final_url = None
        supp_files = []
        used_elsevier_api = False
        used_wiley_api = False

        # Track whether we need to try web scraping as fallback
        api_insufficient_content = False

        # Try Elsevier API first if this is an Elsevier article
        if doi and self.elsevier_api.is_elsevier_doi(doi):
            elsevier_markdown, elsevier_error = self._try_elsevier_api(doi, pmid)
            if elsevier_markdown:
                main_markdown = elsevier_markdown
                final_url = f"https://doi.org/{doi}"
                used_elsevier_api = True
                # Still need to get supplements via scraping
                supp_files = self.doi_resolver.resolve_and_scrape_supplements(doi, pmid, self.scraper)
            elif elsevier_error and "insufficient content" in elsevier_error.lower():
                api_insufficient_content = True
                print(f"  → Elsevier API returned abstract only, falling back to web scraping...")

        # Try Wiley API if this is a Wiley article and Elsevier didn't work
        if not main_markdown and doi and self.wiley_api.is_wiley_doi(doi):
            wiley_markdown, wiley_error = self._try_wiley_api(doi, pmid)
            if wiley_markdown:
                main_markdown = wiley_markdown
                final_url = f"https://doi.org/{doi}"
                used_wiley_api = True
                # Still need to get supplements via scraping
                supp_files = self.doi_resolver.resolve_and_scrape_supplements(doi, pmid, self.scraper)
            elif wiley_error and "insufficient content" in wiley_error.lower():
                api_insufficient_content = True
                print(f"  → Wiley API returned abstract only, falling back to web scraping...")

        # Fall back to DOI resolver if publisher APIs didn't work or returned insufficient content
        if not main_markdown and doi:
            # Try to fetch full text and supplements from publisher using DOI
            main_markdown, final_url, supp_files = self.doi_resolver.resolve_and_fetch_fulltext(
                doi, pmid, self.scraper
            )

            # Validate that we got actual article content, not just an abstract
            MIN_ARTICLE_LENGTH = 1000  # Minimum characters for a valid full-text article
            if main_markdown and len(main_markdown) < MIN_ARTICLE_LENGTH:
                print(f"  ⚠ DOI content too short ({len(main_markdown)} chars) - likely abstract only")
                print(f"  ❌ Skipping - does not contain valid full article text")
                with open(self.paywalled_log, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([pmid, 'DOI content too short (abstract only)', f"https://doi.org/{doi}"])
                return False, "DOI content too short (abstract only)"

        # If no DOI or DOI-based methods failed, try free URL
        if not main_markdown and free_url:
            # No DOI, but we have a direct URL to the free full text
            parsed_url = urlparse(free_url)
            if parsed_url.netloc in self.SUSPICIOUS_FREE_URL_DOMAINS:
                print(f"  - Skipping suspicious free URL on {parsed_url.netloc} (likely non-article content)")
                with open(self.paywalled_log, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([pmid, 'Suspicious free URL skipped', free_url])
                return False, "Suspicious free URL"

            print(f"  - No DOI, attempting to fetch from free URL: {free_url}")

            # Check if this is an Elsevier URL and try API first
            if self.elsevier_api.is_available and self.elsevier_api.is_elsevier_url(free_url):
                pii = self.elsevier_api.extract_pii_from_url(free_url)
                if pii:
                    print(f"  Trying Elsevier API for PII: {pii}")
                    elsevier_markdown, elsevier_error = self.elsevier_api.fetch_fulltext(pii=pii)
                    if elsevier_markdown:
                        main_markdown = elsevier_markdown
                        final_url = free_url
                        used_elsevier_api = True
                        print(f"  ✓ Full text retrieved via Elsevier API ({len(main_markdown)} characters)")
                        # Get supplements via web scraping
                        try:
                            response = self.session.get(free_url, allow_redirects=True, timeout=10)
                            response.raise_for_status()
                            supp_files = self.scraper.scrape_elsevier_supplements(response.text, response.url)
                        except Exception:
                            supp_files = []

            # Check if this is a Wiley URL and try API
            if not main_markdown and self.wiley_api.is_available and self.wiley_api.is_wiley_url(free_url):
                extracted_doi = self.wiley_api.extract_doi_from_url(free_url)
                if extracted_doi:
                    print(f"  Trying Wiley API for DOI: {extracted_doi}")
                    wiley_markdown, wiley_error = self.wiley_api.fetch_fulltext(doi=extracted_doi)
                    if wiley_markdown:
                        main_markdown = wiley_markdown
                        final_url = free_url
                        used_wiley_api = True
                        print(f"  ✓ Full text retrieved via Wiley API ({len(main_markdown)} characters)")
                        # Get supplements via web scraping
                        try:
                            response = self.session.get(free_url, allow_redirects=True, timeout=10)
                            response.raise_for_status()
                            supp_files = self.scraper.scrape_generic_supplements(response.text, response.url)
                        except Exception:
                            supp_files = []

            # Fall back to web scraping if API didn't work
            if not main_markdown:
                try:
                    response = self.session.get(free_url, allow_redirects=True, timeout=10)
                    response.raise_for_status()
                    final_url = response.url
                    print(f"  ✓ Retrieved free full text page")

                    # Handle Elsevier linkinghub redirects
                    domain = urlparse(final_url).netloc
                    if "linkinghub.elsevier.com" in domain:
                        try:
                            pii_match = re.search(r'/pii/([^/?]+)', final_url)
                            if pii_match:
                                pii = pii_match.group(1)
                                sciencedirect_url = f"https://www.sciencedirect.com/science/article/pii/{pii}"
                                print(f"  → Attempting to access ScienceDirect page: {sciencedirect_url}")
                                redirect_response = self.session.get(sciencedirect_url, allow_redirects=True, timeout=30)
                                redirect_response.raise_for_status()
                                final_url = redirect_response.url
                                response = redirect_response
                                domain = urlparse(final_url).netloc
                        except Exception as e:
                            print(f"  - Could not follow redirect from linkinghub: {e}")

                    # Extract full text and supplements from the page
                    html_content = response.text
                    main_markdown, title = self.scraper.extract_fulltext(html_content, final_url)

                    # Validate that we actually got article content (not a database page or error)
                    # Article text should be reasonably long and contain typical article indicators
                    MIN_ARTICLE_LENGTH = 1000  # Minimum characters for a valid article
                    if main_markdown and len(main_markdown) >= MIN_ARTICLE_LENGTH:
                        print(f"  ✓ Extracted full text ({len(main_markdown)} characters)")
                    elif main_markdown:
                        # Got some content but it's too short - likely not an article
                        print(f"  ⚠ Extracted content too short ({len(main_markdown)} chars) - likely not full article text")
                        print(f"  ❌ Skipping this URL as it doesn't contain valid article content")
                        with open(self.paywalled_log, 'a', newline='') as f:
                            writer = csv.writer(f)
                            writer.writerow([pmid, 'Free URL content too short (not an article)', free_url])
                        return False, "Free URL content invalid"
                    else:
                        print(f"  ❌ Could not extract full text from page")

                    # Get supplements from the page
                    domain = urlparse(final_url).netloc

                    # Route to domain-specific scraper
                    if "nature.com" in domain:
                        supp_files = self.scraper.scrape_nature_supplements(html_content, final_url)
                    elif any(d in domain for d in ["gimjournal.org", "sciencedirect.com", "elsevier.com"]):
                        supp_files = self.scraper.scrape_elsevier_supplements(html_content, final_url)
                    else:
                        supp_files = self.scraper.scrape_generic_supplements(html_content, final_url)

                except requests.exceptions.RequestException as e:
                    print(f"  ❌ Failed to fetch free full text from {free_url}: {e}")
                    with open(self.paywalled_log, 'a', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([pmid, f'Free full text fetch failed: {e}', free_url])
                    return False, "Free text fetch failed"
        else:
            # No DOI and no free URL
            print(f"  ❌ No DOI or free URL available to fetch full text")
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Free full text but no DOI or URL available', pubmed_url])
            return False, "No DOI or URL for free text"

        if not main_markdown:
            print(f"  ❌ Could not retrieve full text from publisher")
            fallback_url = final_url or (f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Free full text extraction failed', fallback_url])
            return False, "Free text extraction failed"

        print(f"  ✓ Full text retrieved from publisher ({len(main_markdown)} characters)")

        # Create supplements directory
        supplements_dir = self.output_dir / f"{pmid}_supplements"
        supplements_dir.mkdir(exist_ok=True)

        print(f"  Found {len(supp_files)} supplemental files")

        supplement_markdown = ""
        downloaded_count = 0

        # Download and convert each supplement
        for idx, supp in enumerate(supp_files, 1):
            url = supp.get('url', '')
            filename = supp.get('name', f'supplement_{idx}')

            if not url:
                continue

            file_path = supplements_dir / filename
            print(f"    Downloading: {filename}")

            if self.download_supplement(url, file_path, pmid, filename):
                downloaded_count += 1

                ext = file_path.suffix.lower()

                supplement_markdown += f"\n\n# SUPPLEMENTAL FILE {idx}: {filename}\n\n"

                # Convert supplement to markdown based on file type
                if ext in ['.xlsx', '.xls']:
                    supplement_markdown += self.converter.excel_to_markdown(file_path)
                elif ext == '.docx':
                    supplement_markdown += self.converter.docx_to_markdown(file_path)
                elif ext == '.doc':
                    supplement_markdown += self.converter.doc_to_markdown(file_path)
                elif ext == '.pdf':
                    supplement_markdown += self.converter.pdf_to_markdown(file_path)
                elif ext in ['.txt', '.csv']:
                    try:
                        text = file_path.read_text(encoding='utf-8', errors='ignore')
                        supplement_markdown += text + "\n\n"
                    except Exception as e:
                        supplement_markdown += f"[Error reading text file: {e}]\n\n"
                else:
                    supplement_markdown += f"[File available at: {file_path}]\n\n"

            time.sleep(0.5)

        # Create unified markdown file
        unified_content = main_markdown + supplement_markdown

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(unified_content)

        if used_elsevier_api:
            source_tag = "[via Elsevier API]"
        elif used_wiley_api:
            source_tag = "[via Wiley API]"
        else:
            source_tag = "[from publisher]"
        print(f"  ✅ Created: {output_file.name} ({downloaded_count} supplements) {source_tag}")

        # Run data scout to create condensed DATA_ZONES.md
        self._run_data_scout(pmid, unified_content)

        # Log success with special marker for publisher-sourced content
        if used_elsevier_api:
            source_marker = 'ELSEVIER_API'
        elif used_wiley_api:
            source_marker = 'WILEY_API'
        else:
            source_marker = 'PUBLISHER_FREE'
        with open(self.success_log, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([pmid, source_marker, downloaded_count])

        return True, str(output_file)

    def harvest(self, pmids: List[str], delay: float = 2.0):
        """
        Harvest full-text and supplements for a list of PMIDs.

        Args:
            pmids: List of PubMed IDs to process
            delay: Delay in seconds between processing each PMID
        """
        print(f"Starting harvest for {len(pmids)} PMIDs")
        print(f"Output directory: {self.output_dir.absolute()}\n")

        successful = 0
        failed = 0

        for idx, pmid in enumerate(pmids, 1):
            print(f"[{idx}/{len(pmids)}]", end=" ")

            success, result = self.process_pmid(pmid)

            if success:
                successful += 1
            else:
                failed += 1

            if idx < len(pmids):
                time.sleep(delay)

        print(f"\n{'='*60}")
        print(f"Harvest complete!")
        print(f"  ✅ Successful: {successful}")
        print(f"  ❌ Failed: {failed}")
        print(f"  Output directory: {self.output_dir.absolute()}")
        print(f"  Success log: {self.success_log}")
        print(f"  Paywalled log: {self.paywalled_log}")
        print(f"{'='*60}")

    # Backward-compatible methods for tests and legacy code
    def pmid_to_pmcid(self, pmid: str):
        """Backward-compatible wrapper for PMC API."""
        return self.pmc_api.pmid_to_pmcid(pmid)

    def get_doi_from_pmid(self, pmid: str):
        """Backward-compatible wrapper for PMC API."""
        return self.pmc_api.get_doi_from_pmid(pmid)

    def get_fulltext_xml(self, pmcid: str):
        """Backward-compatible wrapper for PMC API."""
        return self.pmc_api.get_fulltext_xml(pmcid)

    def _get_supplemental_files_from_doi(self, doi: str, pmid: str):
        """Backward-compatible wrapper for DOI resolver."""
        return self.doi_resolver.resolve_and_scrape_supplements(doi, pmid, self.scraper)
