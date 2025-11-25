"""
PMC API Module

Handles interactions with NCBI PubMed Central APIs:
- PMID to PMCID conversion
- Full-text XML retrieval from EuropePMC
- DOI fetching from PubMed records
"""

import os
import time
from typing import Optional
import requests
from Bio import Entrez


# Configure Entrez
Entrez.email = os.getenv("ENTREZ_EMAIL", "your.email@example.com")
Entrez.tool = "PMCHarvester"
Entrez.api_key = os.getenv("NCBI_API_KEY")  # Optional: get API key from https://www.ncbi.nlm.nih.gov/account/


# Rate limiting: NCBI allows 3 req/sec without API key, 10 req/sec with API key
# We'll be conservative and use 0.34 seconds between requests (3 req/sec) without key
# and 0.11 seconds (9 req/sec) with key to stay safely under the limit
RATE_LIMIT_DELAY = 0.11 if Entrez.api_key else 0.34
_last_request_time = 0


class PMCAPIClient:
    """Client for PubMed Central API operations."""

    def __init__(self, session: Optional[requests.Session] = None):
        """
        Initialize PMC API client.

        Args:
            session: Optional requests session with configured headers
        """
        self.session = session or requests.Session()

    def _rate_limit(self):
        """Enforce NCBI rate limiting by adding delays between requests."""
        global _last_request_time
        current_time = time.time()
        time_since_last = current_time - _last_request_time

        if time_since_last < RATE_LIMIT_DELAY:
            sleep_time = RATE_LIMIT_DELAY - time_since_last
            time.sleep(sleep_time)

        _last_request_time = time.time()

    def pmid_to_pmcid(self, pmid: str) -> Optional[str]:
        """
        Convert PMID to PMCID using NCBI E-utilities.

        Args:
            pmid: PubMed ID

        Returns:
            PMCID string (e.g., "PMC1234567") or None if not found
        """
        try:
            self._rate_limit()
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
            self._rate_limit()
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
