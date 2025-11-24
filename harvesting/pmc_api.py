"""
PMC API Module

Handles interactions with NCBI PubMed Central APIs:
- PMID to PMCID conversion
- Full-text XML retrieval from EuropePMC
- DOI fetching from PubMed records
"""

import os
from typing import Optional
import requests
from Bio import Entrez


# Configure Entrez
Entrez.email = os.getenv("ENTREZ_EMAIL", "your.email@example.com")
Entrez.tool = "PMCHarvester"


class PMCAPIClient:
    """Client for PubMed Central API operations."""

    def __init__(self, session: Optional[requests.Session] = None):
        """
        Initialize PMC API client.

        Args:
            session: Optional requests session with configured headers
        """
        self.session = session or requests.Session()

    def pmid_to_pmcid(self, pmid: str) -> Optional[str]:
        """
        Convert PMID to PMCID using NCBI E-utilities.

        Args:
            pmid: PubMed ID

        Returns:
            PMCID string (e.g., "PMC1234567") or None if not found
        """
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
        """
        Fetch the DOI for a given PMID.

        Args:
            pmid: PubMed ID

        Returns:
            DOI string or None if not found
        """
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
        """
        Download full-text XML from EuropePMC.

        Args:
            pmcid: PubMed Central ID (e.g., "PMC1234567")

        Returns:
            XML content as string or None if not available
        """
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
