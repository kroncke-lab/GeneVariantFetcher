#!/usr/bin/env python3
"""
PubMed Central Full-Text & Supplemental Materials Harvester

This script converts PMIDs to PMCIDs, downloads full-text XML and all supplemental files,
and converts everything to unified markdown files for LLM processing.

Usage:
    python harvest_pmc_fulltext.py
"""

import os
import time
import json
import csv
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import xml.etree.ElementTree as ET
from urllib.parse import urlparse, urljoin

import requests
from bs4 import BeautifulSoup
from Bio import Entrez
import pandas as pd

try:
    from markitdown import MarkItDown
    MARKITDOWN_AVAILABLE = True
except ImportError:
    MARKITDOWN_AVAILABLE = False
    print("WARNING: markitdown not available. Will use basic conversion for docx/xlsx.")


Entrez.email = "brett.kroncke@gmail.com"
Entrez.tool = "PMCHarvester"


class PMCHarvester:
    """Harvests full-text and supplemental materials from PubMed Central."""
    
    def __init__(self, output_dir: str = "pmc_harvest"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.paywalled_log = self.output_dir / "paywalled_missing.csv"
        self.success_log = self.output_dir / "successful_downloads.csv"
        self.markitdown = MarkItDown() if MARKITDOWN_AVAILABLE else None
        
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.9',
            'Accept-Encoding': 'gzip, deflate, br',
            'Referer': 'https://pubmed.ncbi.nlm.nih.gov/',
            'DNT': '1',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1'
        })

        with open(self.paywalled_log, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['PMID', 'Reason', 'URL'])
        
        with open(self.success_log, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['PMID', 'PMCID', 'Supplements_Downloaded'])
    
    def pmid_to_pmcid(self, pmid: str) -> Optional[str]:
        """Convert PMID to PMCID using NCBI E-utilities."""
        try:
            handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid, linkname="pubmed_pmc")
            record = Entrez.read(handle)
            handle.close()
            
            if record[0]["LinkSetDb"]:
                pmc_ids = record[0]["LinkSetDb"][0]["Link"]
                if pmc_ids:
                    pmcid = pmc_ids[0]["Id"]
                    return f"PMC{pmcid}"
            
            return None
        except Exception as e:
            print(f"  Error converting PMID {pmid} to PMCID: {e}")
            return None

    def get_doi_from_pmid(self, pmid: str) -> Optional[str]:
        """Fetch the DOI for a given PMID."""
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            article = records['PubmedArticle'][0]['MedlineCitation']['Article']
            for item in article['ELocationID']:
                if item.attributes['EIdType'] == 'doi':
                    return str(item)
            # Fallback for older records
            if 'ArticleIdList' in article:
                 for identifier in article['ArticleIdList']:
                    if identifier.attributes.get('IdType') == 'doi':
                        return str(identifier)
            return None
        except Exception as e:
            print(f"  Error fetching DOI for PMID {pmid}: {e}")
            return None

    def get_fulltext_xml(self, pmcid: str) -> Optional[str]:
        """Download full-text XML from EuropePMC."""
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
        
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            return response.text
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 404:
                print(f"  Full-text not available for {pmcid} (404)")
            else:
                print(f"  HTTP error fetching {pmcid}: {e}")
            return None
        except Exception as e:
            print(f"  Error fetching full-text for {pmcid}: {e}")
            return None
    
    def get_supplemental_files(self, pmcid: str, pmid: str, doi: Optional[str]) -> List[Dict]:
        """
        Orchestrates fetching supplemental files, first via API, then by scraping.
        """
        # 1. First, try the EuropePMC API
        api_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
        try:
            response = self.session.get(api_url, timeout=30)
            response.raise_for_status()
            data = response.json()
            if 'result' in data and 'supplementaryFiles' in data['result']:
                print("  ✓ Found supplemental files via EuropePMC API")
                return data['result']['supplementaryFiles']
        except Exception as e:
            print(f"  - EuropePMC supplemental files API failed for {pmcid}: {e}")

        # 2. If API fails or returns no files, fall back to DOI-based scraping
        if doi:
            print(f"  - API failed. Falling back to DOI scraping for {doi}")
            return self._get_supplemental_files_from_doi(doi, pmid)
        else:
            print("  - API failed and no DOI available. Cannot scrape.")
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Supplemental files API failed, no DOI', pmc_url])
            return []

    def _get_supplemental_files_from_doi(self, doi: str, pmid: str) -> List[Dict]:
        """Resolve DOI to final URL and route to the appropriate scraper."""
        # Try multiple times with delays to handle rate limiting
        max_retries = 3
        for attempt in range(max_retries):
            try:
                # Add a delay before each attempt (except the first)
                if attempt > 0:
                    delay = 2 ** attempt  # Exponential backoff: 2s, 4s
                    print(f"  ⏳ Retry {attempt + 1}/{max_retries} after {delay}s delay...")
                    time.sleep(delay)

                response = self.session.get(f"https://doi.org/{doi}", allow_redirects=True, timeout=30)

                # Check if we got a successful response
                if response.status_code == 200:
                    final_url = response.url
                    domain = urlparse(final_url).netloc
                    print(f"  ✓ DOI resolved to: {final_url}")
                    break
                elif response.status_code == 403:
                    print(f"  ⚠️  403 Forbidden (attempt {attempt + 1}/{max_retries})")
                    if attempt == max_retries - 1:
                        raise requests.exceptions.HTTPError(f"403 Forbidden after {max_retries} attempts")
                    continue
                else:
                    response.raise_for_status()

            except requests.exceptions.RequestException as e:
                if attempt == max_retries - 1:
                    print(f"  ❌ DOI resolution failed for {doi} after {max_retries} attempts: {e}")
                    # Try constructing URL directly as a fallback
                    constructed_url = self._construct_url_from_doi(doi)
                    if constructed_url:
                        print(f"  ⚡ Trying direct URL construction: {constructed_url}")
                        try:
                            response = self.session.get(constructed_url, timeout=30)
                            print(f"  📊 Response status: {response.status_code}")
                            if response.status_code == 200:
                                final_url = constructed_url
                                domain = urlparse(final_url).netloc
                                print(f"  ✓ Successfully accessed via constructed URL")
                                break
                            else:
                                print(f"  ⚠️  Constructed URL returned status {response.status_code}")
                        except requests.exceptions.RequestException as url_error:
                            print(f"  ❌ Constructed URL also failed: {url_error}")

                    # If all attempts failed, log and return empty
                    with open(self.paywalled_log, 'a', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([pmid, f'DOI resolution failed: {doi}', f"https://doi.org/{doi}"])
                    return []
                continue
            except Exception as e:
                print(f"  ❌ Unexpected error during DOI resolution: {e}")
                with open(self.paywalled_log, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([pmid, f'DOI resolution error: {doi}', f"https://doi.org/{doi}"])
                return []

        # Route to specific scraper based on resolved domain
        if "nature.com" in domain:
            return self._scrape_nature_supplements(response.text, final_url)
        elif "gimjournal.org" in domain or "sciencedirect.com" in domain:
            return self._scrape_elsevier_supplements(response.text, final_url)
        else:
            print(f"  - No specific scraper for domain: {domain}. Using generic scraper.")
            return self._scrape_generic_supplements(response.text, final_url)

    def _construct_url_from_doi(self, doi: str) -> Optional[str]:
        """Construct publisher URL directly from DOI when resolution fails."""
        # Nature journals: 10.1038/...
        if doi.startswith('10.1038/'):
            return f"https://www.nature.com/articles/{doi}"
        # Elsevier/GIM: 10.1016/...
        elif doi.startswith('10.1016/'):
            # Try both ScienceDirect and GIM journal
            # Note: This may not always work as Elsevier uses PIIs internally
            return f"https://www.gimjournal.org/article/{doi}/fulltext"
        return None

    def _scrape_nature_supplements(self, html: str, base_url: str) -> List[Dict]:
        """Scrape supplemental files from a Nature journal page."""
        print("  Scraping with _scrape_nature_supplements...")
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
                if '/articles/' in href and '/figures/' not in href:
                    url = urljoin(base_url, href)
                    filename = Path(urlparse(url).path).name
                    if filename:
                        found_files.append({'url': url, 'name': filename})
                        print(f"      Found potential supplement: {filename}")
        
        if not found_files:
            print("    'Supplementary Information' section not found or empty. Trying generic scan.")
            return self._scrape_generic_supplements(html, base_url)

        return found_files

    def _scrape_elsevier_supplements(self, html: str, base_url: str) -> List[Dict]:
        """Scrape supplemental files from an Elsevier/GIM journal page."""
        print("  Scraping with _scrape_elsevier_supplements...")
        soup = BeautifulSoup(html, 'html.parser')
        found_files = []

        # 1. Look for specific "Download file" links in ScienceDirect
        for a in soup.find_all('a', class_='S_C_9cf8451f', href=True):
            if 'Download file' in a.get_text():
                url = urljoin(base_url, a['href'])
                filename = a.get_text().replace('Download file', '').strip()
                if filename:
                    found_files.append({'url': url, 'name': filename})
                    print(f"    Found ScienceDirect download link: {filename}")
        if found_files:
            return found_files

        # 2. Look for data in <script type="application/json">
        for script in soup.find_all('script', type='application/json'):
            try:
                data = json.loads(script.string)
                # This path can be very specific and fragile
                supp_data = data.get('article', {}).get('supplementaryMaterials', {}).get('supplementaryMaterial', [])
                for item in supp_data:
                    url = item.get('downloadUrl')
                    filename = item.get('title')
                    if url and filename:
                        found_files.append({'url': url, 'name': filename})
                        print(f"    Found supplement in JSON data: {filename}")
                if found_files:
                    return found_files
            except (json.JSONDecodeError, AttributeError):
                continue
        
        # 3. Regex for "mmc" (multimedia component) links
        mmc_links = re.findall(r'href="(/cms/attachment/[^"]+/mmc\d+\.(?:pdf|docx|xlsx|zip))"', html)
        for link in mmc_links:
            url = urljoin(base_url, link)
            filename = Path(urlparse(url).path).name
            if filename and not any(f['name'] == filename for f in found_files):
                found_files.append({'url': url, 'name': filename})
                print(f"    Found MMC link via regex: {filename}")
        if found_files:
            return found_files

        print("    No specific Elsevier/GIM supplements found. Trying generic scan.")
        return self._scrape_generic_supplements(html, base_url)

    def _scrape_generic_supplements(self, html: str, base_url: str) -> List[Dict]:
        """A best-effort generic scraper for supplemental files."""
        print("  Scraping with _scrape_generic_supplements...")
        soup = BeautifulSoup(html, 'html.parser')
        found_files = []
        
        keywords = ['supplement', 'supporting', 'appendix']
        file_extensions = ['.pdf', '.docx', '.xlsx', '.csv', '.zip', '.rar', '.gz']
        
        for a in soup.find_all('a', href=True):
            link_text = a.get_text().lower()
            href = a['href'].lower()
            
            # Check if link text or href contain any of the keywords or file extensions
            if any(keyword in link_text for keyword in keywords) or \
               any(href.endswith(ext) for ext in file_extensions):
                
                url = urljoin(base_url, a['href'])
                filename = Path(urlparse(url).path).name
                
                if filename:
                    found_files.append({'url': url, 'name': filename})
                    print(f"    Found potential supplement: {filename}")

        return found_files

    def download_supplement(self, url: str, output_path: Path, pmid: str, filename: str) -> bool:
        """Download a supplemental file."""
        try:
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
                writer.writerow([pmid, f'Supplemental file failed: {filename}', url])
            return False
    
    def xml_to_markdown(self, xml_content: str) -> str:
        """Convert PubMed Central XML to markdown."""
        try:
            root = ET.fromstring(xml_content)
            
            markdown = "# MAIN TEXT\n\n"
            
            title_elem = root.find(".//article-title")
            if title_elem is not None:
                title = ''.join(title_elem.itertext()).strip()
                markdown += f"## {title}\n\n"
            
            abstract_elem = root.find(".//abstract")
            if abstract_elem is not None:
                markdown += "### Abstract\n\n"
                abstract_text = ''.join(abstract_elem.itertext()).strip()
                markdown += f"{abstract_text}\n\n"
            
            body_elem = root.find(".//body")
            if body_elem is not None:
                for sec in body_elem.findall(".//sec"):
                    title_elem = sec.find("title")
                    if title_elem is not None:
                        sec_title = ''.join(title_elem.itertext()).strip()
                        markdown += f"### {sec_title}\n\n"
                    
                    for p in sec.findall(".//p"):
                        para_text = ''.join(p.itertext()).strip()
                        if para_text:
                            markdown += f"{para_text}\n\n"
            
            return markdown
        except Exception as e:
            print(f"  Error parsing XML: {e}")
            return "# MAIN TEXT\n\n[Error parsing XML content]\n\n"
    
    def excel_to_markdown(self, file_path: Path) -> str:
        """Convert Excel file to markdown tables."""
        try:
            xl_file = pd.ExcelFile(file_path)
            markdown = ""
            
            for sheet_name in xl_file.sheet_names:
                df = pd.read_excel(file_path, sheet_name=sheet_name)
                
                if df.empty:
                    continue
                
                markdown += f"#### Sheet: {sheet_name}\n\n"
                
                if len(df) > 100:
                    markdown += f"*Note: Showing first 100 rows of {len(df)} total rows*\n\n"
                    df_display = df.head(100)
                else:
                    df_display = df
                
                markdown += df_display.to_markdown(index=False)
                markdown += "\n\n"
            
            return markdown
        except Exception as e:
            print(f"    Error converting Excel {file_path}: {e}")
            return f"[Error converting Excel file: {e}]\n\n"
    
    def docx_to_markdown(self, file_path: Path) -> str:
        """Convert Word document to markdown."""
        if self.markitdown:
            try:
                result = self.markitdown.convert(str(file_path))
                return result.text_content
            except Exception as e:
                print(f"    Error converting DOCX with markitdown {file_path}: {e}")
                return f"[Error converting DOCX file: {e}]\n\n"
        else:
            try:
                from docx import Document
                doc = Document(file_path)
                text = "\n\n".join([para.text for para in doc.paragraphs if para.text.strip()])
                return text + "\n\n"
            except Exception as e:
                print(f"    Error converting DOCX {file_path}: {e}")
                return f"[Error converting DOCX file: {e}]\n\n"
    
    def pdf_to_markdown(self, file_path: Path) -> str:
        """Convert PDF to markdown."""
        if self.markitdown:
            try:
                result = self.markitdown.convert(str(file_path))
                return result.text_content
            except Exception as e:
                print(f"    Error converting PDF with markitdown {file_path}: {e}")
                return f"[PDF file available at: {file_path.name}]\n\n"
        else:
            return f"[PDF file available at: {file_path.name}]\n\n"
    
    def process_pmid(self, pmid: str) -> Tuple[bool, str]:
        """Process a single PMID: convert to PMCID, download content, create unified markdown."""
        print(f"\nProcessing PMID: {pmid}")
        
        doi = self.get_doi_from_pmid(pmid)
        if doi:
            print(f"  ✓ DOI: {doi}")
        else:
            print("  - No DOI found for this PMID.")

        pmcid = self.pmid_to_pmcid(pmid)
        
        if not pmcid:
            print(f"  ❌ No PMCID found (likely paywalled)")
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'No PMCID found', pubmed_url])
            return False, "No PMCID"
        
        print(f"  ✓ PMCID: {pmcid}")
        
        xml_content = self.get_fulltext_xml(pmcid)
        
        if not xml_content:
            print(f"  ❌ Full-text not available")
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Full-text not available', pmc_url])
            return False, "No full-text"
        
        print(f"  ✓ Full-text XML retrieved")
        
        main_markdown = self.xml_to_markdown(xml_content)
        
        supplements_dir = self.output_dir / f"{pmid}_supplements"
        supplements_dir.mkdir(exist_ok=True)
        
        supp_files = self.get_supplemental_files(pmcid, pmid, doi)
        print(f"  Found {len(supp_files)} supplemental files")
        
        supplement_markdown = ""
        downloaded_count = 0
        
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
                
                if ext in ['.xlsx', '.xls']:
                    supplement_markdown += self.excel_to_markdown(file_path)
                elif ext in ['.docx']:
                    supplement_markdown += self.docx_to_markdown(file_path)
                elif ext == '.pdf':
                    supplement_markdown += self.pdf_to_markdown(file_path)
                else:
                    supplement_markdown += f"[File available at: {file_path}]\n\n"
            
            time.sleep(0.5)
        
        unified_content = main_markdown + supplement_markdown
        
        output_file = self.output_dir / f"{pmid}_FULL_CONTEXT.md"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(unified_content)
        
        print(f"  ✅ Created: {output_file.name} ({downloaded_count} supplements)")
        
        with open(self.success_log, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([pmid, pmcid, downloaded_count])
        
        return True, str(output_file)
    
    def harvest(self, pmids: List[str], delay: float = 2.0):
        """Harvest full-text and supplements for a list of PMIDs."""
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


def main():
    """Main entry point."""
    
    pmids = [
        '34931732', # Nature
        '35443093', # GIM
        '33442691', # Some other
    ]
    
    print("PubMed Central Full-Text & Supplemental Materials Harvester")
    print("="*60)
    print(f"\nNOTE: Make sure to set your email in Entrez.email")
    print(f"Current email: {Entrez.email}\n")
    
    harvester = PMCHarvester(output_dir="pmc_harvest")
    
    harvester.harvest(pmids, delay=2.0)


if __name__ == "__main__":
    main()
