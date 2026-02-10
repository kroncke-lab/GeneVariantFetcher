"""
Europe PMC Full-Text API Handler

This module provides specialized access to Europe PMC's full-text articles and supplementary files.
Based on the official API documentation at: https://europepmc.org/RestfulWebService
"""

import json
import logging
import requests
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from xml.etree import ElementTree as ET

from utils.retry_utils import api_retry

logger = logging.getLogger(__name__)

class EuropePMCClient:
    """
    Client for Europe PMC RESTful Web Service APIs.
    
    Provides access to:
    - Full-text OA articles (over 6.5 million available)
    - XML and HTML format options
    - Supplementary files
    - Citation data
    - Cross-references to other databases
    
    API base URLs:
    - Production: https://www.ebi.ac.uk/europepmc/webservices/rest/
    - Test: https://www.ebi.ac.uk/europepmc/webservices/test/rest/
    """
    
    BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest"
    
    def __init__(self, timeout: int = 30):
        """
        Initialize Europe PMC client.
        
        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout
        self.session = requests.Session()
        
        # Set appropriate headers for Europe PMC
        self.session.headers.update({
            'User-Agent': 'GeneVariantFetcher/1.0 (Bot for clinical variant extraction)',
            'Accept': 'application/json, application/xml, text/html',
            'Accept-Language': 'en-US,en;q=0.9',
        })
    
    def search_papers(self, query: str, max_results: int = 100) -> List[Dict]:
        """Search for papers using the search API."""
        url = f"{self.BASE_URL}/search"
        params = {
            'query': query,
            'format': 'json',
            'pageSize': max_results,
            'resultType': 'core'
        }
        
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            
            data = response.json()
            results = data.get('resultList', {}).get('result', [])
            
            logger.info(f"Search query '{query}' returned {len(results)} results")
            return results
            
        except Exception as e:
            logger.error(f"Europe PMC search failed for query '{query}': {e}")
            return []
    
    def get_paper_metadata(self, pmid: str) -> Optional[Dict]:
        """
        Get comprehensive metadata for a specific paper by PMID.
        
        Args:
            pmid: PubMed ID
            
        Returns:
            Dictionary with complete paper metadata or None if not found
        """
        url = f"{self.BASE_URL}/search"
        params = {
            'query': f'EXT_ID:{pmid}',
            'format': 'json',
            'resultType': 'core'
        }
        
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            
            data = response.json()
            results = data.get('resultList', {}).get('result', [])
            
            if not results:
                logger.warning(f"No metadata found for PMID {pmid}")
                return None
                
            result = results[0]
            
            # Parse authors
            authors = []
            for author in result.get('authorList', {}).get('author', []):
                if isinstance(author, dict):
                    authors.append({
                        'forename': author.get('firstName', ''),
                        'lastname': author.get('lastName', ''),
                        'affiliation': author.get('affiliation', '')
                    })
            
            # Escape pmcid
            pmcid = result.get('pmcid', None)
            if pmcid and not pmcid.startswith('PMC'):
                pmcid = f'PMC{pmcid}'
            
            metadata = {
                'pmid': str(result.get('pmid', pmid)),
                'pmcid': pmcid,
                'doi': result.get('doi', ''),
                'title': result.get('title', ''),
                'abstract': result.get('abstractText', ''),
                'authors': authors,
                'journal_info': {
                    'title': result.get('journalTitle', ''),
                    'year': result.get('pubYear', ''),
                    'month': result.get('pubMonth', ''),
                    'volume': result.get('journalVolume', ''),
                    'issue': result.get('issue', ''),
                    'pages': result.get('pageInfo', '')
                },
                'fulltext_frequency': result.get('hasFullText', 'N'),
                'is_open_access': result.get('isOpenAccess', 'N') == 'Y',
                'citation_count': int(result.get('citedByCount', 0)),
                'links': {
                    'europepmc_url': f"https://europepmc.org/article/MED/{pmid}",
                    'fulltext_html_url': f"https://europepmc.org/articles/{pmcid}" if pmcid else None,
                    'fulltext_pdf_url': f"https://europepmc.org/backend/ptpmcrender.fcgi?accid={pmcid}&blobtype=pdf" if pmcid else None
                },
                'textmining': {
                    'has_textmined_terms': result.get('hasTextMinedTerms', 'N') == 'Y',
                    'has_references': result.get('hasReferences', 'N') == 'Y',
                    'has_annotations': result.get('hasAnnotations', 'N') == 'Y'
                }
            }
            
            return metadata
            
        except Exception as e:
            logger.error(f"Failed to fetch metadata for PMID {pmid}: {e}")
            return None
    
    @api_retry
    def get_fulltext_xml(self, pmcid: str) -> Optional[str]:
        """
        Download full-text XML for a PMC article from Europe PMC.
        
        Args:
            pmcid: PubMed Central ID (with or without PMC prefix)
            
        Returns:
            XML content as string or None if not available
        """
        if not pmcid:
            logger.warning("No PMCID provided")
            return None
            
        # Ensure PMCID has proper format
        pmcid = pmcid.strip()
        if not pmcid.startswith('PMC'):
            pmcid = f'PMC{pmcid}'
        
        url = f"{self.BASE_URL}/{pmcid}/fullTextXML"
        
        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            
            xml_content = response.text
            
            # Validate XML content
            if not xml_content or "<error>" in xml_content.lower():
                logger.warning(f"Invalid XML response for {pmcid}")
                return None
                
            logger.info(f"Successfully retrieved XML for {pmcid}")
            return xml_content
            
        except Exception as e:
            logger.error(f"Failed to retrieve XML for {pmcid}: {e}")
            return None
    
    @api_retry
    def get_fulltext_html(self, pmcid: str) -> Optional[str]:
        """
        Download full-text HTML for a PMC article from Europe PMC.
        
        Args:
            pmcid: PubMed Central ID
            
        Returns:
            HTML content as string or None if not available
        """
        if not pmcid:
            logger.warning("No PMCID provided")
            return None
            
        pmcid = pmcid.strip()
        if not pmcid.startswith('PMC'):
            pmcid = f'PMC{pmcid}'
        
        url = f"https://europepmc.org/articles/{pmcid}"
        
        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            
            html_content = response.text
            
            if not html_content or "404" in html_content.lower():
                logger.warning(f"Invalid HTML response for {pmcid}")
                return None
                
            logger.info(f"Successfully retrieved HTML for {pmcid}")
            return html_content
            
        except Exception as e:
            logger.error(f"Failed to retrieve HTML for {pmcid}: {e}")
            return None
    
    @api_retry
    def get_supplementary_files(self, pmcid: str) -> List[Dict]:
        """
        Retrieve supplementary files for a PMC article.
        
        Args:
            pmcid: PubMed Central ID
            
        Returns:
            List of supplementary file information dictionaries
        """
        if not pmcid:
            return []
            
        pmcid = pmcid.strip()
        if not pmcid.startswith('PMC'):
            pmcid = f'PMC{pmcid}'
        
        url = f"{self.BASE_URL}/{pmcid}/supplementaryFiles"
        
        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            
            data = response.json()
            supplementary_files = data.get('result', {}).get('supplementaryFiles', [])
            
            logger.info(f"Found {len(supplementary_files)} supplementary files for {pmcid}")
            return supplementary_files
            
        except Exception as e:
            logger.error(f"Failed to retrieve supplementary files for {pmcid}: {e}")
            return []
    
    def download_full_paper(self, pmid: str, output_dir: Path) -> Dict:
        """
        Comprehensive method to download a complete paper from Europe PMC.
        
        Args:
            pmid: PubMed ID
            output_dir: Directory to save files
            
        Returns:
            Dictionary with download results and metadata
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        results = {
            'pmid': pmid,
            'success': False,
            'metadata': None,
            'files_downloaded': [],
            'errors': []
        }
        
        # Get metadata first
        metadata = self.get_paper_metadata(pmid)
        if not metadata:
            error_msg = f"Could not retrieve metadata for PMID {pmid}"
            results['errors'].append(error_msg)
            logger.error(error_msg)
            return results
            
        results['metadata'] = metadata
        
        # Save metadata as JSON
        metadata_file = output_dir / f"{pmid}_europepmc_metadata.json"
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        results['files_downloaded'].append(str(metadata_file))
        
        # Check if we have PMCID for full text access
        pmcid = metadata.get('pmcid')
        if not pmcid:
            error_msg = f"No PMCID found for PMID {pmid} - full text not available"
            results['errors'].append(error_msg)
            logger.warning(error_msg)
            return results
        
        # Download full-text XML
        xml_content = self.get_fulltext_xml(pmcid)
        if xml_content:
            xml_file = output_dir / f"{pmid}_europepmc_full.xml"
            with open(xml_file, 'w', encoding='utf-8') as f:
                f.write(xml_content)
            results['files_downloaded'].append(str(xml_file))
            
            # Also try PDF for easier reading
            pdf_url = metadata.get('links', {}).get('fulltext_pdf_url')
            if pdf_url:
                try:
                    pdf_response = self.session.get(pdf_url, timeout=self.timeout)
                    pdf_response.raise_for_status()
                    
                    pdf_file = output_dir / f"{pmid}_europepmc_full.pdf"
                    with open(pdf_file, 'wb') as f:
                        f.write(pdf_response.content)
                    results['files_downloaded'].append(str(pdf_file))
                    
                except Exception as e:
                    logger.warning(f"Could not download PDF for {pmid}: {e}")
        
        # Download supplementary files
        supplements = self.get_supplementary_files(pmcid)
        if supplements:
            supplements_dir = output_dir / f"{pmid}_supplements" 
            supplements_dir.mkdir(exist_ok=True)
            
            for idx, supp in enumerate(supplements, 1):
                try:
                    supp_url = supp.get('downloadUrl', '') or supp.get('url', '')
                    supp_filename = supp.get('filename', '') or supp.get('name', f"supplement_{idx}")
                    
                    if not supp_url or not supp_filename:
                        continue
                    
                    # Download the file
                    response = self.session.get(supp_url, timeout=self.timeout)
                    response.raise_for_status()
                    
                    supp_file = supplements_dir / supp_filename
                    with open(supp_file, 'wb') as f:
                        f.write(response.content)
                    
                    results['files_downloaded'].append(str(supp_file))
                    logger.info(f"Downloaded supplementary file: {supp_filename}")
                    
                except Exception as e:
                    logger.warning(f"Failed to download supplement {idx}: {e}")
                    results['errors'].append(f"Failed to download supplement {idx}: {e}")
        
        results['success'] = len(results['files_downloaded']) > 0
        return results
    
    def get_citations(self, pmid: str) -> List[Dict]:
        """
        Get papers that cite the specified paper.
        
        Args:
            pmid: PubMed ID
            
        Returns:
            List of citing papers with metadata
        """
        url = f"{self.BASE_URL}/MED/{pmid}/citations"
        params = {
            'format': 'json',
            'pageSize': 100
        }
        
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            
            data = response.json()
            citations = data.get('resultList', {}).get('result', [])
            
            logger.info(f"Found {len(citations)} citations for PMID {pmid}")
            return citations
            
        except Exception as e:
            logger.error(f"Failed to retrieve citations for PMID {pmid}: {e}")
            return []
    
    def get_references(self, pmid: str) -> List[Dict]:
        """
        Get references cited by the specified paper.
        
        Args:
            pmid: PubMed ID
            
        Returns:
            List of referenced papers with metadata
        """
        url = f"{self.BASE_URL}/MED/{pmid}/references"
        params = {
            'format': 'json',
            'pageSize': 500
        }
        
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            
            data = response.json()
            references = data.get('resultList', {}).get('result', [])
            
            logger.info(f"Found {len(references)} references for PMID {pmid}")
            return references
            
        except Exception as e:
            logger.error(f"Failed to retrieve references for PMID {pmid}: {e}")
            return []