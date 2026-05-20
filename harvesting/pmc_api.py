"""
PMC API Module

Handles interactions with NCBI PubMed Central APIs:
- PMID to PMCID conversion
- Full-text XML retrieval from NCBI
- DOI fetching from PubMed records
- Free full-text status detection from PubMed
"""

import logging
import os
import threading
import time
from typing import List, Optional, Tuple
from urllib.parse import urlparse

import requests
from Bio import Entrez

from utils.retry_utils import api_retry

logger = logging.getLogger(__name__)

# Configure Entrez - uses NCBI_EMAIL from environment
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.tool = "PMCHarvester"
Entrez.api_key = os.getenv("NCBI_API_KEY")


# Rate limiting: NCBI allows 3 req/sec without API key, 10 req/sec with API key
# We'll be conservative and use 0.34 seconds between requests (3 req/sec) without key
# and 0.11 seconds (9 req/sec) with key to stay safely under the limit
RATE_LIMIT_DELAY = 0.11 if Entrez.api_key else 0.34

# Thread-safe rate limiting
_rate_limit_lock = threading.Lock()
_last_request_time = 0


NON_ARTICLE_LINKOUT_DOMAINS = {
    "antibodies.cancer.gov",
    "antibodies-online.com",
    "assays.cancer.gov",
    "biocyc.org",
    "clinicaltrials.gov",
    "disgenet.org",
    "genecards.org",
    "glygen.org",
    "hgmd.cf.ac.uk",
    "lens.org",
    "lovd.nl",
    "malacards.org",
    "medlineplus.gov",
    "ncbi.nlm.nih.gov/clinvar",
    "ncbi.nlm.nih.gov/gene",
    "ncbi.nlm.nih.gov/nuccore",
    "ncbi.nlm.nih.gov/protein",
    "omim.org",
    "orphanet.org",
    "proteinatlas.org",
    "scite.ai",
    "uniprot.org",
}

PUBLISHER_LINKOUT_DOMAINS = {
    "academic.oup.com",
    "ahajournals.org",
    "biomedcentral.com",
    "bmj.com",
    "cell.com",
    "frontiersin.org",
    "gimjournal.org",
    "hindawi.com",
    "jamanetwork.com",
    "journals.lww.com",
    "karger.com",
    "mdpi.com",
    "nature.com",
    "nejm.org",
    "onlinelibrary.wiley.com",
    "plos.org",
    "pnas.org",
    "science.org",
    "sciencedirect.com",
    "springer.com",
    "tandfonline.com",
    "thelancet.com",
    "wiley.com",
}


def _host_and_path(url: str) -> tuple[str, str]:
    parsed = urlparse(url)
    host = parsed.netloc.lower()
    if host.startswith("www."):
        host = host[4:]
    return host, parsed.path.lower()


def _domain_rule_matches(url: str, rule: str) -> bool:
    host, path = _host_and_path(url)
    rule = rule.lower()
    if "/" in rule:
        rule_host, rule_path = rule.split("/", 1)
        rule_path = "/" + rule_path
        return (
            host == rule_host or host.endswith(f".{rule_host}")
        ) and path.startswith(rule_path)
    return host == rule or host.endswith(f".{rule}")


def is_non_article_linkout_url(url: str) -> bool:
    """Return True for PubMed LinkOuts that are databases/tools, not papers."""
    return any(
        _domain_rule_matches(url, domain) for domain in NON_ARTICLE_LINKOUT_DOMAINS
    )


def is_publisher_linkout_url(url: str) -> bool:
    """Return True for LinkOut URLs that are likely publisher article pages."""
    return any(
        _domain_rule_matches(url, domain) for domain in PUBLISHER_LINKOUT_DOMAINS
    )


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
        """Enforce NCBI rate limiting by adding delays between requests (thread-safe)."""
        global _last_request_time
        with _rate_limit_lock:
            current_time = time.time()
            time_since_last = current_time - _last_request_time

            if time_since_last < RATE_LIMIT_DELAY:
                sleep_time = RATE_LIMIT_DELAY - time_since_last
                time.sleep(sleep_time)

            _last_request_time = time.time()

    @api_retry
    def pmid_to_pmcid(self, pmid: str) -> Optional[str]:
        """
        Convert PMID to PMCID using NCBI E-utilities.

        Args:
            pmid: PubMed ID

        Returns:
            PMCID string (e.g., "PMC1234567") or None if not found
        """
        last_error: Exception | None = None
        for attempt in range(1, 4):
            handle = None
            try:
                self._rate_limit()
                handle = Entrez.elink(
                    dbfrom="pubmed", db="pmc", id=pmid, linkname="pubmed_pmc"
                )
                record = Entrez.read(handle)

                if record[0]["LinkSetDb"]:
                    pmc_ids = record[0]["LinkSetDb"][0]["Link"]
                    if pmc_ids:
                        pmcid = pmc_ids[0]["Id"]
                        return f"PMC{pmcid}"

                return None
            except Exception as e:
                last_error = e
                if attempt < 3:
                    time.sleep(2**attempt)
                    continue
            finally:
                if handle is not None:
                    try:
                        handle.close()
                    except Exception:
                        pass

        print(f"  Error converting PMID {pmid} to PMCID after retries: {last_error}")
        return None

    @api_retry
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

            pubmed_article = records["PubmedArticle"][0]
            article = pubmed_article["MedlineCitation"]["Article"]

            # Check ELocationID in Article (newer records)
            for item in article.get("ELocationID", []):
                if item.attributes.get("EIdType") == "doi":
                    return str(item)

            # Fallback: Check ArticleIdList in PubmedData (older records)
            # Note: ArticleIdList is in PubmedData, not Article
            pubmed_data = pubmed_article.get("PubmedData", {})
            for identifier in pubmed_data.get("ArticleIdList", []):
                if identifier.attributes.get("IdType") == "doi":
                    return str(identifier)

            return None
        except Exception as e:
            print(f"  Error fetching DOI for PMID {pmid}: {e}")
            return None

    @api_retry
    def _fetch_xml_from_ncbi(self, numeric_pmcid: str) -> str:
        """Helper to fetch XML from NCBI with retries."""
        self._rate_limit()
        handle = Entrez.efetch(
            db="pmc", id=numeric_pmcid, rettype="full", retmode="xml"
        )
        xml_content = handle.read()
        handle.close()
        return xml_content.decode("utf-8")

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

        numeric_pmcid = pmcid[3:]

        # Try NCBI E-utilities first
        try:
            xml_string = self._fetch_xml_from_ncbi(numeric_pmcid)

            if not xml_string or "<error>" in xml_string.lower():
                print(
                    f"  NCBI E-utilities returned empty or error for {pmcid}, trying Europe PMC..."
                )
            else:
                return xml_string
        except Exception as e:
            # This can include HTTP errors from Entrez or parsing errors
            print(f"  NCBI E-utilities failed for {pmcid}: {e}")
            print("  Trying Europe PMC fallback...")

        # Fallback to Europe PMC full-text API
        try:
            europe_pmc_url = (
                f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
            )
            response = self.session.get(europe_pmc_url, timeout=30)
            response.raise_for_status()

            xml_string = response.text

            # Check for valid XML response
            if (
                xml_string
                and "<article" in xml_string.lower()
                and "<error>" not in xml_string.lower()
            ):
                print("  ✓ Retrieved full-text from Europe PMC")
                return xml_string
            else:
                print(f"  Europe PMC returned invalid or empty response for {pmcid}")
                return None

        except Exception as e:
            print(f"  Europe PMC fallback also failed for {pmcid}: {e}")
            return None

    @api_retry
    def is_free_full_text(self, pmid: str) -> Tuple[bool, Optional[str]]:
        """
        Check if a PubMed article is marked as free full text via the publisher.

        This detects articles that have "Free article" or "Free full text" indicators
        on PubMed but do NOT have a PMCID (i.e., not in PMC). These articles provide
        free access through the publisher's website.

        Args:
            pmid: PubMed ID

        Returns:
            Tuple of (is_free: bool, publisher_url: Optional[str])
            - is_free: True if article is marked as free full text
            - publisher_url: Direct URL to free full text if available
        """
        try:
            self._rate_limit()

            # Use elink to get LinkOut information for this PMID
            # cmd="llinks" returns external links including free full text sources
            free_indicators = [
                "free full text",
                "free article",
                "free",
                "publisher free",
                "open access",
            ]

            # Fetch LinkOut information via elink
            handle = Entrez.elink(dbfrom="pubmed", id=pmid, cmd="llinks")
            records = Entrez.read(handle)
            handle.close()

            # First try structured LinkOut helper (easier to monkeypatch in tests)
            free_text_url = None
            prioritized_url = None  # For publisher domains
            is_free = False

            if records and len(records) > 0:
                # Check IdUrlList for LinkOut URLs
                id_url_list = records[0].get("IdUrlList", {})
                if id_url_list:
                    id_url_set = id_url_list.get("IdUrlSet", [])
                    if id_url_set and len(id_url_set) > 0:
                        obj_urls = id_url_set[0].get("ObjUrl", [])
                        for obj_url in obj_urls:
                            # Check for free full text attributes
                            attributes = obj_url.get("Attribute", [])
                            url = obj_url.get("Url", "")

                            # PubMed uses these attributes to indicate free access
                            free_indicators = [
                                "free full text",
                                "free article",
                                "free",
                                "full text",
                                "publisher free",
                                "open access",
                            ]

                            attr_text = " ".join(str(a).lower() for a in attributes)
                            if any(
                                indicator in attr_text for indicator in free_indicators
                            ):
                                is_free = True
                                if url:
                                    url_str = str(url)

                                    # Skip irrelevant domains
                                    if is_non_article_linkout_url(url_str):
                                        print(
                                            f"  - Skipping irrelevant LinkOut URL: {url_str}"
                                        )
                                        continue

                                    # Prioritize known publisher domains
                                    if is_publisher_linkout_url(url_str):
                                        if not prioritized_url:
                                            prioritized_url = url_str
                                            print(f"  ✓ Found publisher URL: {url_str}")
                                    elif not free_text_url:
                                        # Only use as fallback if no prioritized URL found yet
                                        free_text_url = url_str
                                        print(f"  ~ Found generic free URL: {url_str}")

            # Prefer prioritized URL (from known publishers) over generic free URL
            final_url = prioritized_url or free_text_url

            # If elink didn't find free text, check the PubMed record itself
            # Some records have the "Free" indicator in different locations
            if not is_free:
                try:
                    self._rate_limit()
                    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
                    records = Entrez.read(handle)
                    handle.close()

                    if records.get("PubmedArticle"):
                        article = records["PubmedArticle"][0]
                        pubmed_data = article.get("PubmedData", {})

                        # Check ArticleIdList for free full text URL
                        article_ids = pubmed_data.get("ArticleIdList", [])
                        for aid in article_ids:
                            id_type = aid.attributes.get("IdType", "")
                            # Check for PMC free article status
                            if id_type == "pmc" and "free" in str(aid).lower():
                                is_free = True

                        # Check history for free access indicators
                        history = pubmed_data.get("History", [])
                        for status in history:
                            pub_status = status.attributes.get("PubStatus", "")
                            if "free" in pub_status.lower():
                                is_free = True

                except Exception as e:
                    print(f"    - Secondary free text check failed: {e}")

            return is_free, final_url

        except Exception as e:
            print(f"  Error checking free full text status for PMID {pmid}: {e}")
            return False, None

    @api_retry
    def get_pubmed_linkout_urls(self, pmid: str) -> List[str]:
        """
        Get all LinkOut URLs for a PubMed article.

        This retrieves external links to publisher websites, which can be used
        to access free full text articles.

        Args:
            pmid: PubMed ID

        Returns:
            List of dictionaries with 'url', 'provider', and 'category' keys
        """
        try:
            self._rate_limit()
            handle = Entrez.elink(dbfrom="pubmed", id=pmid, cmd="llinks")
            record = Entrez.read(handle)
            handle.close()

            links = []
            if record and len(record) > 0:
                id_url_list = record[0].get("IdUrlList", {})
                if id_url_list:
                    id_url_set = id_url_list.get("IdUrlSet", [])
                    if id_url_set and len(id_url_set) > 0:
                        obj_urls = id_url_set[0].get("ObjUrl", [])
                        for obj_url in obj_urls:
                            url = obj_url.get("Url", "")
                            provider = obj_url.get("Provider", {}).get(
                                "Name", "Unknown"
                            )
                            category = (
                                obj_url.get("Category", ["Unknown"])[0]
                                if obj_url.get("Category")
                                else "Unknown"
                            )
                            attributes = obj_url.get("Attribute", [])

                            if url:
                                links.append(
                                    {
                                        "url": str(url),
                                        "provider": str(provider),
                                        "category": str(category),
                                        "attributes": [str(a) for a in attributes],
                                    }
                                )

            return links

        except Exception as e:
            print(f"  Error getting LinkOut URLs for PMID {pmid}: {e}")
            return []
