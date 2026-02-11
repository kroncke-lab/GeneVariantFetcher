"""PMC supplement fetcher (Tier 1 - free, no API key required).

Retrieves supplementary files from PubMed Central using:
1. Europe PMC supplementary files API (primary)
2. NCBI eUtils OA service for PMC XML with embedded supplement links
3. FTP fallback: ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/

PMID 24667783 is a known good test case (multiple supplements).
"""

import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urljoin
from xml.etree import ElementTree as ET

import requests
from Bio import Entrez

from utils.retry_utils import api_retry

from .base import SupplementFetcher, SupplementFile

logger = logging.getLogger(__name__)

# Configure Entrez from environment
Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.tool = "GeneVariantFetcher"
Entrez.api_key = os.getenv("NCBI_API_KEY")


class PMCSupplementFetcher(SupplementFetcher):
    """Fetch supplements from PMC (free tier).

    Strategy:
        1. Resolve PMID -> PMCID via NCBI eLink
        2. Query Europe PMC supplementary files endpoint
        3. If no results, parse PMC full-text XML for supplement links
        4. FTP fallback for OA subset articles
    """

    EUROPEPMC_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
    PMC_OA_BASE = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"
    PMC_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/pmc"

    def __init__(self, timeout: int = 30):
        super().__init__(timeout=timeout)
        # Europe PMC prefers JSON
        self.session.headers.update(
            {
                "Accept": "application/json, application/xml, text/html",
            }
        )

    def fetch(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Fetch supplement metadata from PMC sources.

        Tries Europe PMC API first, then falls back to XML parsing.
        """
        pmcid = self._resolve_pmcid(pmid)
        if not pmcid:
            logger.info(f"PMID {pmid} has no PMCID - not in PMC")
            return []

        # Tier 1a: Europe PMC supplementary files endpoint
        supplements = self._fetch_europepmc_supplements(pmcid)
        if supplements:
            logger.info(
                f"Europe PMC returned {len(supplements)} supplements for {pmcid}"
            )
            return supplements

        # Tier 1b: Parse PMC full-text XML for supplement links
        supplements = self._fetch_from_pmc_xml(pmcid)
        if supplements:
            logger.info(f"PMC XML yielded {len(supplements)} supplements for {pmcid}")
            return supplements

        # Tier 1c: OA service FTP links
        supplements = self._fetch_from_oa_service(pmcid)
        if supplements:
            logger.info(
                f"OA service returned {len(supplements)} supplements for {pmcid}"
            )
            return supplements

        logger.info(f"No PMC supplements found for PMID {pmid} ({pmcid})")
        return []

    def list_supplements(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Alias for fetch()."""
        return self.fetch(pmid, doi)

    # ------------------------------------------------------------------
    # PMID -> PMCID resolution
    # ------------------------------------------------------------------

    @api_retry
    def _resolve_pmcid(self, pmid: str) -> Optional[str]:
        """Convert PMID to PMCID using NCBI ID converter API."""
        url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
        params = {
            "ids": pmid,
            "format": "json",
            "tool": "GeneVariantFetcher",
            "email": os.getenv("NCBI_EMAIL", ""),
        }
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()

            records = data.get("records", [])
            if records:
                pmcid = records[0].get("pmcid", "")
                if pmcid:
                    logger.debug(f"PMID {pmid} -> {pmcid}")
                    return pmcid

            return None
        except Exception as e:
            logger.warning(f"PMCID resolution failed for PMID {pmid}: {e}")
            return None

    # ------------------------------------------------------------------
    # Europe PMC supplementary files
    # ------------------------------------------------------------------

    @api_retry
    def _fetch_europepmc_supplements(self, pmcid: str) -> List[SupplementFile]:
        """Query Europe PMC for supplementary files.

        The endpoint returns different content types:
        - application/zip for OA articles (the actual supplement archive)
        - application/xml with an error for non-OA articles
        - application/json with file listings (rare)
        """
        if not pmcid.startswith("PMC"):
            pmcid = f"PMC{pmcid}"

        url = f"{self.EUROPEPMC_BASE}/{pmcid}/supplementaryFiles"

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            if e.response is not None and e.response.status_code == 404:
                logger.debug(f"No Europe PMC supplements for {pmcid}")
                return []
            raise
        except Exception as e:
            logger.warning(
                f"Europe PMC supplementary files request failed for {pmcid}: {e}"
            )
            return []

        content_type = response.headers.get("Content-Type", "")

        # ZIP response means OA article with downloadable supplements
        if "application/zip" in content_type:
            return [
                SupplementFile(
                    url=url,
                    name=f"{pmcid}_supplements.zip",
                    source="pmc_europepmc",
                    pmcid=pmcid,
                    mime_type="application/zip",
                    size_bytes=len(response.content),
                    description="Europe PMC supplementary files archive",
                )
            ]

        # XML error response (non-OA article)
        if "xml" in content_type:
            if (
                "not open access" in response.text.lower()
                or "<error" in response.text.lower()
            ):
                logger.debug(f"{pmcid} is not open access in Europe PMC")
                return []

        # JSON response with file listings
        try:
            data = response.json()
        except ValueError:
            logger.debug(f"Non-JSON response from Europe PMC supplements for {pmcid}")
            return []

        files_data = data.get("supplementaryFiles", [])
        if not files_data:
            files_data = data.get("result", {}).get("supplementaryFiles", [])

        results = []
        seen_urls = set()

        for entry in files_data:
            file_url = entry.get("url", "") or entry.get("downloadUrl", "")
            filename = (
                entry.get("fileName", "")
                or entry.get("name", "")
                or entry.get("label", "")
            )

            if not file_url:
                continue

            normalized = file_url.split("#")[0].rstrip("/")
            if normalized in seen_urls:
                continue
            seen_urls.add(normalized)

            if not filename:
                filename = Path(file_url).name or "supplement"

            results.append(
                SupplementFile(
                    url=file_url,
                    name=self._clean_filename(filename),
                    source="pmc_europepmc",
                    pmcid=pmcid,
                    mime_type=entry.get("mimeType", ""),
                    description=entry.get("description", "") or entry.get("title", ""),
                )
            )

        return results

    # ------------------------------------------------------------------
    # PMC full-text XML parsing
    # ------------------------------------------------------------------

    @api_retry
    def _fetch_from_pmc_xml(self, pmcid: str) -> List[SupplementFile]:
        """Parse PMC full-text XML to extract supplement links."""
        if not pmcid.startswith("PMC"):
            pmcid = f"PMC{pmcid}"

        url = f"{self.EUROPEPMC_BASE}/{pmcid}/fullTextXML"

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            xml_text = response.text
        except Exception as e:
            logger.debug(f"Could not fetch XML for {pmcid}: {e}")
            return []

        if not xml_text or "<error>" in xml_text.lower():
            return []

        return self._parse_supplement_links_from_xml(xml_text, pmcid)

    def _parse_supplement_links_from_xml(
        self, xml_text: str, pmcid: str
    ) -> List[SupplementFile]:
        """Extract supplement URLs from PMC JATS XML.

        Handles two common patterns in JATS XML:
        1. <supplementary-material xlink:href="filename.pdf">
        2. <supplementary-material><media xlink:href="filename.pdf"/></supplementary-material>
        """
        results = []
        seen_urls = set()
        xlink = "{http://www.w3.org/1999/xlink}"

        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError as e:
            logger.warning(f"XML parse error for {pmcid}: {e}")
            return []

        def _add_supplement(href: str, label: str, caption: str, mime: str) -> None:
            if not href:
                return
            # Make absolute URL
            if not href.startswith("http"):
                href = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/bin/{href}"

            normalized = href.split("#")[0].rstrip("/")
            if normalized in seen_urls:
                return
            seen_urls.add(normalized)

            filename = Path(href).name or label or "supplement"
            results.append(
                SupplementFile(
                    url=href,
                    name=self._clean_filename(filename),
                    source="pmc_xml",
                    pmcid=pmcid,
                    mime_type=mime,
                    description=f"{label}: {caption}".strip(": ")
                    if (label or caption)
                    else "",
                )
            )

        for supp in root.iter("supplementary-material"):
            # Get label/caption from the supplementary-material element
            label_el = supp.find(".//label")
            caption_el = supp.find(".//caption")
            label = ""
            if label_el is not None and label_el.text:
                label = label_el.text.strip()
            caption = ""
            if caption_el is not None:
                caption = " ".join(caption_el.itertext()).strip()

            # Pattern 1: href directly on supplementary-material
            href = supp.get(f"{xlink}href", "")
            if href:
                _add_supplement(href, label, caption, "")

            # Pattern 2: <media> children with xlink:href
            for media in supp.iter("media"):
                media_href = media.get(f"{xlink}href", "")
                if media_href:
                    mimetype = media.get("mimetype", "")
                    subtype = media.get("mime-subtype", "")
                    mime = f"{mimetype}/{subtype}" if mimetype and subtype else ""
                    _add_supplement(media_href, label, caption, mime)

        return results

    # ------------------------------------------------------------------
    # NCBI OA service (FTP links)
    # ------------------------------------------------------------------

    @api_retry
    def _fetch_from_oa_service(self, pmcid: str) -> List[SupplementFile]:
        """Use NCBI OA service to find FTP download links for supplements."""
        if not pmcid.startswith("PMC"):
            pmcid = f"PMC{pmcid}"

        params = {
            "id": pmcid,
            "format": "json",
        }

        try:
            response = self.session.get(
                self.PMC_OA_BASE, params=params, timeout=self.timeout
            )
            response.raise_for_status()
            data = response.json()
        except Exception as e:
            logger.debug(f"OA service failed for {pmcid}: {e}")
            return []

        records = data.get("records", [])
        if not records:
            return []

        results = []
        seen_urls = set()

        for record in records:
            # OA service returns links to tar.gz packages containing supplements
            for link_key in ("ftp", "https"):
                link = record.get(link_key, "")
                if not link:
                    continue

                normalized = link.split("#")[0].rstrip("/")
                if normalized in seen_urls:
                    continue
                seen_urls.add(normalized)

                filename = Path(link).name or f"{pmcid}_package"
                results.append(
                    SupplementFile(
                        url=link,
                        name=self._clean_filename(filename),
                        source="pmc_oa_ftp",
                        pmcid=pmcid,
                        description="OA package (may contain supplements)",
                    )
                )

        return results
