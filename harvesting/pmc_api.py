"""
PMC API Module

Handles interactions with NCBI PubMed Central APIs:
- PMID to PMCID conversion
- Full-text XML retrieval from NCBI
- DOI fetching from PubMed records
"""

import os
import time
from typing import List, Optional, Tuple
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
        Download full-text XML from NCBI using E-utilities, with fallback to Europe PMC.

        Args:
            pmcid: PubMed Central ID (e.g., "PMC1234567")

        Returns:
            XML content as string or None if not available
        """
        if not pmcid or not pmcid.upper().startswith("PMC"):
            print(f"  Invalid PMCID format: {pmcid}")
            return None

        # E-utilities require the numeric part of the PMCID
        numeric_pmcid = pmcid[3:]

        # Try NCBI E-utilities first
        try:
            self._rate_limit()
            handle = Entrez.efetch(db="pmc", id=numeric_pmcid, rettype="full", retmode="xml")
            xml_content = handle.read()
            handle.close()

            # The response is bytes, so decode it to a string
            xml_string = xml_content.decode('utf-8')

            # Check for an empty or error response from NCBI
            if not xml_string or "<error>" in xml_string.lower():
                 print(f"  NCBI E-utilities returned empty or error for {pmcid}, trying Europe PMC...")
            else:
                return xml_string
        except Exception as e:
            # This can include HTTP errors from Entrez or parsing errors
            print(f"  NCBI E-utilities failed for {pmcid}: {e}")
            print(f"  Trying Europe PMC fallback...")

        # Fallback to Europe PMC full-text API
        try:
            europe_pmc_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
            response = self.session.get(europe_pmc_url, timeout=30)
            response.raise_for_status()

            xml_string = response.text

            # Check for valid XML response
            if xml_string and "<article" in xml_string.lower() and not "<error>" in xml_string.lower():
                print(f"  ✓ Retrieved full-text from Europe PMC")
                return xml_string
            else:
                print(f"  Europe PMC returned invalid or empty response for {pmcid}")
                return None

        except Exception as e:
            print(f"  Europe PMC fallback also failed for {pmcid}: {e}")
            return None

    def get_pubmed_linkout_urls(self, pmid: str) -> List[dict]:
        """
        Retrieve LinkOut URLs associated with a PubMed record.

        Args:
            pmid: PubMed ID

        Returns:
            List of dictionaries containing provider, url, and attributes metadata
        """
        try:
            self._rate_limit()
            handle = Entrez.elink(dbfrom="pubmed", id=pmid, cmd="llinks", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            linkouts = []
            for link_set in records:
                for link_db in link_set.get("LinkSetDb", []):
                    provider = link_db.get("LinkName", "")
                    for link in link_db.get("Link", []):
                        url = link.get("Url") or link.get("Id")
                        attributes = link.get("Attr") if isinstance(link.get("Attr"), list) else link.get("Attr")
                        attributes = attributes or []
                        if url:
                            linkouts.append({
                                "provider": provider,
                                "url": url,
                                "attributes": attributes,
                            })
            return linkouts
        except Exception as e:
            print(f"  Error fetching LinkOut URLs for PMID {pmid}: {e}")
            return []

    def is_free_full_text(self, pmid: str) -> Tuple[bool, Optional[str]]:
        """
        Determine whether a PubMed article has free full text available.

        Args:
            pmid: PubMed ID

        Returns:
            Tuple of (is_free: bool, free_url: Optional[str])
        """
        # If a PMCID exists, the article is available in PMC and therefore free
        pmcid = self.pmid_to_pmcid(pmid)
        if pmcid:
            return True, f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"

        linkouts = self.get_pubmed_linkout_urls(pmid)

        free_indicators = [
            "free full text",
            "free article",
            "free",
            "publisher free",
            "open access"
        ]

        for link in linkouts:
            provider = (link.get("provider") or "").lower()
            attributes = [attr.lower() for attr in link.get("attributes", []) if isinstance(attr, str)]
            url = link.get("url")

            if any(indicator in provider for indicator in free_indicators):
                return True, url

            if any(indicator in attr for attr in attributes for indicator in free_indicators):
                return True, url

        return False, None
