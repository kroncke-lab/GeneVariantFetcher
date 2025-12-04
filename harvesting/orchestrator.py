"""
Orchestrator Module

Main PMCHarvester class that coordinates all harvesting operations:
- Converts PMIDs to PMCIDs
- Downloads full-text and supplemental files
- Creates unified markdown files for LLM processing
"""

import os
import time
import csv
from pathlib import Path
from typing import List, Tuple
import requests

from .pmc_api import PMCAPIClient
from .doi_resolver import DOIResolver
from .supplement_scraper import SupplementScraper
from .format_converters import FormatConverter


class PMCHarvester:
    """Harvests full-text and supplemental materials from PubMed Central."""

    def __init__(self, output_dir: str = "pmc_harvest"):
        """
        Initialize PMC Harvester.

        Args:
            output_dir: Directory to save harvested files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.paywalled_log = self.output_dir / "paywalled_missing.csv"
        self.success_log = self.output_dir / "successful_downloads.csv"

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

        # Initialize log files
        with open(self.paywalled_log, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['PMID', 'Reason', 'URL'])

        with open(self.success_log, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['PMID', 'PMCID', 'Supplements_Downloaded'])

    def get_supplemental_files(self, pmcid: str, pmid: str, doi: str) -> List[dict]:
        """
        Orchestrates fetching supplemental files, first via API, then by scraping.

        Args:
            pmcid: PubMed Central ID
            pmid: PubMed ID (for logging)
            doi: Digital Object Identifier (for fallback scraping)

        Returns:
            List of supplement file dictionaries with 'url' and 'name' keys
        """
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
                    return data['result']['supplementaryFiles']
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
                                supp_files.append({'url': url, 'name': name})
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
                pmc_response = self.session.get(pmc_url, timeout=30)
                pmc_response.raise_for_status()
                
                # Use generic scraper for PMC pages
                pmc_files = self.scraper.scrape_generic_supplements(pmc_response.text, pmc_url)
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

    def download_supplement(self, url: str, output_path: Path, pmid: str, filename: str) -> bool:
        """
        Download a supplemental file.

        Args:
            url: URL of the supplement file
            output_path: Path to save the file
            pmid: PubMed ID (for logging)
            filename: Name of the file (for logging)

        Returns:
            True if download succeeded, False otherwise
        """
        try:
            # Use the session to download, which includes our headers
            response = self.session.get(url, timeout=60, stream=True)
            response.raise_for_status()

            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            return True
        except Exception as e:
            print(f"    Error downloading {url}: {e}")
            # Log failed supplemental download
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, f'Supplemental file download failed: {filename}', url])
            return False

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

        # Create supplements directory
        supplements_dir = self.output_dir / f"{pmid}_supplements"
        supplements_dir.mkdir(exist_ok=True)

        # Get supplemental files
        supp_files = self.get_supplemental_files(pmcid, pmid, doi)
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
                elif ext in ['.docx']:
                    supplement_markdown += self.converter.docx_to_markdown(file_path)
                elif ext == '.pdf':
                    supplement_markdown += self.converter.pdf_to_markdown(file_path)
                else:
                    supplement_markdown += f"[File available at: {file_path}]\n\n"

            time.sleep(0.5)

        # Create unified markdown file
        unified_content = main_markdown + supplement_markdown

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(unified_content)

        print(f"  ✅ Created: {output_file.name} ({downloaded_count} supplements)")

        # Log success
        with open(self.success_log, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([pmid, pmcid, downloaded_count])

        return True, str(output_file)

    def _process_free_text_pmid(self, pmid: str, doi: str) -> Tuple[bool, str]:
        """
        Process a PMID that has no PMCID but may have free full text via publisher.

        This method checks if the article is marked as "Free Full Text" on PubMed
        and attempts to fetch the content from the publisher's website.

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

        # Need DOI to fetch from publisher
        if not doi:
            print(f"  ❌ No DOI available to fetch free full text from publisher")
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Free full text but no DOI available', pubmed_url])
            return False, "No DOI for free text"

        # Try to fetch full text and supplements from publisher
        main_markdown, final_url, supp_files = self.doi_resolver.resolve_and_fetch_fulltext(
            doi, pmid, self.scraper
        )

        if not main_markdown:
            print(f"  ❌ Could not retrieve full text from publisher")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Free full text extraction failed', final_url or f"https://doi.org/{doi}"])
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
                elif ext in ['.docx']:
                    supplement_markdown += self.converter.docx_to_markdown(file_path)
                elif ext == '.pdf':
                    supplement_markdown += self.converter.pdf_to_markdown(file_path)
                else:
                    supplement_markdown += f"[File available at: {file_path}]\n\n"

            time.sleep(0.5)

        # Create unified markdown file
        unified_content = main_markdown + supplement_markdown

        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(unified_content)

        print(f"  ✅ Created: {output_file.name} ({downloaded_count} supplements) [from publisher]")

        # Log success with special marker for publisher-sourced content
        with open(self.success_log, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([pmid, 'PUBLISHER_FREE', downloaded_count])

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
