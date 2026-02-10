"""Elsevier supplement fetcher (Tier 2 - requires API key).

Uses the Elsevier Article Retrieval API to discover and download
supplementary materials. Requires ELSEVIER_API_KEY in .env.

Elsevier supplements typically follow the "mmc" naming pattern:
    mmc1.pdf, mmc2.xlsx, mmc3.docx, etc.

API docs: https://dev.elsevier.com/documentation/ArticleRetrievalAPI.wadl
"""

import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urljoin, urlparse

import requests

from utils.retry_utils import api_retry
from .base import SupplementFetcher, SupplementFile

logger = logging.getLogger(__name__)


class ElsevierSupplementFetcher(SupplementFetcher):
    """Fetch supplements from Elsevier via their API.

    Requires ELSEVIER_API_KEY environment variable.

    Strategy:
        1. Use DOI to query Article Retrieval API for full article metadata
        2. Extract supplement references from the response
        3. Build download URLs using ScienceDirect supplement patterns
    """

    API_BASE = "https://api.elsevier.com/content"
    SCIENCEDIRECT_BASE = "https://www.sciencedirect.com"

    # Common Elsevier supplement filename pattern
    MMC_PATTERN = re.compile(r"mmc\d+\.(pdf|docx?|xlsx?|csv|pptx?|zip|txt)", re.IGNORECASE)

    def __init__(self, timeout: int = 30, api_key: str = ""):
        super().__init__(timeout=timeout)
        self.api_key = api_key or os.getenv("ELSEVIER_API_KEY", "")
        if self.api_key:
            self.session.headers.update({
                "X-ELS-APIKey": self.api_key,
                "Accept": "application/json",
            })

    @property
    def available(self) -> bool:
        """Check if API key is configured."""
        return bool(self.api_key)

    def fetch(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Fetch Elsevier supplement metadata.

        Requires a DOI to query the Elsevier API effectively.
        """
        if not self.available:
            logger.debug("Elsevier API key not configured, skipping")
            return []

        if not doi:
            logger.debug(f"No DOI for PMID {pmid}, cannot query Elsevier API")
            return []

        # Try article metadata endpoint for supplement references
        supplements = self._fetch_from_article_api(doi, pmid)
        if supplements:
            return supplements

        # Try to discover mmc-pattern supplements from the PII
        pii = self._resolve_pii(doi)
        if pii:
            supplements = self._discover_mmc_supplements(pii, doi, pmid)
            if supplements:
                return supplements

        return []

    def list_supplements(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Alias for fetch()."""
        return self.fetch(pmid, doi)

    # ------------------------------------------------------------------
    # Article Retrieval API
    # ------------------------------------------------------------------

    @api_retry
    def _fetch_from_article_api(self, doi: str, pmid: str) -> List[SupplementFile]:
        """Query Elsevier Article Retrieval API for supplement info."""
        url = f"{self.API_BASE}/article/doi/{doi}"
        params = {"view": "FULL"}

        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            if response.status_code == 404:
                logger.debug(f"Article not found in Elsevier for DOI {doi}")
                return []
            response.raise_for_status()
        except requests.exceptions.HTTPError:
            return []
        except Exception as e:
            logger.warning(f"Elsevier API request failed for DOI {doi}: {e}")
            return []

        try:
            data = response.json()
        except ValueError:
            return []

        return self._extract_supplements_from_response(data, doi, pmid)

    def _extract_supplements_from_response(
        self, data: dict, doi: str, pmid: str
    ) -> List[SupplementFile]:
        """Parse Elsevier API response for supplement references."""
        results = []
        seen_urls = set()

        # Navigate the nested Elsevier response
        article = data.get("full-text-retrieval-response", {})

        # Check for explicit supplement references in coredata
        coredata = article.get("coredata", {})
        pii = coredata.get("pii", "")

        # Look in the object list for supplementary materials
        objects = article.get("objects", {}).get("object", [])
        if isinstance(objects, dict):
            objects = [objects]

        for obj in objects:
            ref = obj.get("ref", {})
            href = ""
            if isinstance(ref, dict):
                href = ref.get("$", "")
            elif isinstance(ref, str):
                href = ref

            if not href:
                continue

            # Only include actual supplement files
            if not any(kw in href.lower() for kw in ("mmc", "suppl", "appendix", "si_")):
                continue

            if not href.startswith("http"):
                href = urljoin(self.SCIENCEDIRECT_BASE, href)

            normalized = href.split("#")[0].rstrip("/")
            if normalized in seen_urls:
                continue
            seen_urls.add(normalized)

            filename = Path(urlparse(href).path).name or "supplement"
            results.append(SupplementFile(
                url=href,
                name=self._clean_filename(filename),
                source="elsevier_api",
                pmid=pmid,
                mime_type=obj.get("mime-type", ""),
                description=obj.get("caption", ""),
            ))

        # Also check for attachment references
        attachments = article.get("attachment", [])
        if isinstance(attachments, dict):
            attachments = [attachments]

        for att in attachments:
            href = att.get("$", "") or att.get("url", "")
            if not href:
                continue

            if not href.startswith("http"):
                href = urljoin(self.SCIENCEDIRECT_BASE, href)

            normalized = href.split("#")[0].rstrip("/")
            if normalized in seen_urls:
                continue
            seen_urls.add(normalized)

            filename = att.get("filename", "") or Path(urlparse(href).path).name or "supplement"
            results.append(SupplementFile(
                url=href,
                name=self._clean_filename(filename),
                source="elsevier_api",
                pmid=pmid,
                mime_type=att.get("mime-type", ""),
            ))

        return results

    # ------------------------------------------------------------------
    # PII resolution and MMC discovery
    # ------------------------------------------------------------------

    @api_retry
    def _resolve_pii(self, doi: str) -> Optional[str]:
        """Resolve DOI to Elsevier PII (article identifier)."""
        url = f"{self.API_BASE}/article/doi/{doi}"
        params = {"view": "META"}

        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            if response.status_code == 404:
                return None
            response.raise_for_status()
            data = response.json()

            coredata = data.get("full-text-retrieval-response", {}).get("coredata", {})
            if not coredata:
                coredata = data.get("abstracts-retrieval-response", {}).get("coredata", {})

            return coredata.get("pii", None)
        except Exception as e:
            logger.debug(f"PII resolution failed for DOI {doi}: {e}")
            return None

    def _discover_mmc_supplements(
        self, pii: str, doi: str, pmid: str
    ) -> List[SupplementFile]:
        """Try common mmc1..mmcN supplement URLs for an Elsevier article.

        Elsevier supplements typically follow:
            https://ars.els-cdn.com/content/image/1-s2.0-{PII}-mmc1.pdf
        """
        results = []
        clean_pii = re.sub(r"[^A-Za-z0-9]", "", pii)

        for i in range(1, 11):  # Try mmc1 through mmc10
            for ext in ("pdf", "xlsx", "docx", "csv", "zip"):
                file_url = (
                    f"https://ars.els-cdn.com/content/image/"
                    f"1-s2.0-{clean_pii}-mmc{i}.{ext}"
                )

                try:
                    resp = self.session.head(file_url, timeout=10, allow_redirects=True)
                    if resp.status_code == 200:
                        results.append(SupplementFile(
                            url=file_url,
                            name=f"mmc{i}.{ext}",
                            source="elsevier_mmc",
                            pmid=pmid,
                        ))
                        break  # Found this mmc number, try next
                except Exception:
                    continue

            # If we didn't find mmc{i} in any extension, stop probing
            if not any(s.name.startswith(f"mmc{i}.") for s in results):
                break

        if results:
            logger.info(f"Discovered {len(results)} mmc supplements for PII {pii}")

        return results
