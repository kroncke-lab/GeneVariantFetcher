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

import requests
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
        
        with open(self.paywalled_log, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['PMID', 'Reason'])
        
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
    
    def get_fulltext_xml(self, pmcid: str) -> Optional[str]:
        """Download full-text XML from EuropePMC."""
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
        
        try:
            response = requests.get(url, timeout=30)
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
    
    def get_supplemental_files(self, pmcid: str) -> List[Dict]:
        """Get list of supplemental files from EuropePMC."""
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if 'result' in data and 'supplementaryFiles' in data['result']:
                return data['result']['supplementaryFiles']
            return []
        except Exception as e:
            print(f"  Error fetching supplemental files for {pmcid}: {e}")
            return []
    
    def download_supplement(self, url: str, output_path: Path) -> bool:
        """Download a supplemental file."""
        try:
            response = requests.get(url, timeout=60, stream=True)
            response.raise_for_status()
            
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            return True
        except Exception as e:
            print(f"    Error downloading {url}: {e}")
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
        
        pmcid = self.pmid_to_pmcid(pmid)
        
        if not pmcid:
            print(f"  ❌ No PMCID found (likely paywalled)")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'No PMCID found'])
            return False, "No PMCID"
        
        print(f"  ✓ PMCID: {pmcid}")
        
        xml_content = self.get_fulltext_xml(pmcid)
        
        if not xml_content:
            print(f"  ❌ Full-text not available")
            with open(self.paywalled_log, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([pmid, 'Full-text not available'])
            return False, "No full-text"
        
        print(f"  ✓ Full-text XML retrieved")
        
        main_markdown = self.xml_to_markdown(xml_content)
        
        supplements_dir = self.output_dir / f"{pmid}_supplements"
        supplements_dir.mkdir(exist_ok=True)
        
        supp_files = self.get_supplemental_files(pmcid)
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
            
            if self.download_supplement(url, file_path):
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
        '34931732',
        '35443093', 
        '33442691',
    ]
    
    print("PubMed Central Full-Text & Supplemental Materials Harvester")
    print("="*60)
    print(f"\nNOTE: Make sure to set your email in Entrez.email")
    print(f"Current email: {Entrez.email}\n")
    
    harvester = PMCHarvester(output_dir="pmc_harvest")
    
    harvester.harvest(pmids, delay=2.0)


if __name__ == "__main__":
    main()
