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
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
            'Referer': 'https://pubmed.ncbi.nlm.nih.gov/'
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
            data = response.json()
            # The API can return a success response with an empty list of files
            if data.get('result', {}).get('supplementaryFiles'):
                print("  ✓ Found supplemental files via EuropePMC API")
                return data['result']['supplementaryFiles']
        except Exception as e:
            print(f"  - EuropePMC supplemental files API failed for {pmcid}: {e}")

        # 2. If API fails or returns no files, fall back to DOI-based scraping
        if doi:
            print(f"  - API failed or returned no files. Falling back to DOI scraping for {doi}")
            return self.doi_resolver.resolve_and_scrape_supplements(doi, pmid, self.scraper)
        else:
            print("  - API failed and no DOI available. Cannot scrape.")
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
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

        Args:
            pmid: PubMed ID to process

        Returns:
            Tuple of (success: bool, result: str) where result is output file path or error message
        """
        print(f"\nProcessing PMID: {pmid}")

        # Get DOI first (needed for supplement fallback)
        doi = self.pmc_api.get_doi_from_pmid(pmid)
        if doi:
            print(f"  ✓ DOI: {doi}")
        else:
            print("  - No DOI found for this PMID.")

        # Convert PMID to PMCID
        pmcid = self.pmc_api.pmid_to_pmcid(pmid)

        if not pmcid:
            print(f"  ❌ No PMCID found (likely paywalled)")
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'No PMCID found', pubmed_url])
            return False, "No PMCID"

        print(f"  ✓ PMCID: {pmcid}")

        # Get full-text XML
        xml_content = self.pmc_api.get_fulltext_xml(pmcid)

        if not xml_content:
            print(f"  ❌ Full-text not available")
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Full-text not available', pmc_url])
            return False, "No full-text"

        print(f"  ✓ Full-text XML retrieved")

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
