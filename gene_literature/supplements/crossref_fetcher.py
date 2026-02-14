"""CrossRef supplement fetcher (Tier 2.5 - free, no API key required).

Uses the CrossRef REST API to discover supplementary materials via
DOI metadata. CrossRef indexes supplement links from most major publishers.

Polite pool: requests with a mailto: User-Agent get routed to faster servers
and higher rate limits (~50 req/s vs ~1 req/s).

API docs: https://api.crossref.org/swagger-ui/index.html
"""

import logging
import os
import time
from pathlib import Path
from typing import List, Optional
from urllib.parse import urlparse

import requests

from utils.retry_utils import api_retry

from .base import SupplementFetcher, SupplementFile

logger = logging.getLogger(__name__)

# URL path keywords that indicate supplementary material
_SUPPLEMENT_KEYWORDS = {
    "supplement", "supplementary", "suppl", "supp",
    "mmc", "appendix", "additional", "supporting",
    "si_", "s1_", "s2_", "s3_", "table_s", "figure_s",
}

# Intended-application values to skip (not supplement files)
_SKIP_APPLICATIONS = {"similarity-checking", "text-mining"}


class CrossRefSupplementFetcher(SupplementFetcher):
    """Fetch supplement metadata from CrossRef via DOI lookup.

    Strategy:
        1. Query CrossRef works API with the paper's DOI
        2. Parse ``message.link`` array for supplement URLs
           (filtering out similarity-checking and text-mining entries)
        3. Parse ``message.relation.has-supplement`` for URI entries
        4. Deduplicate by normalized URL
    """

    API_BASE = "https://api.crossref.org/works"

    def __init__(self, timeout: int = 30, email: str = ""):
        super().__init__(timeout=timeout)
        self._email = email or os.getenv("NCBI_EMAIL", "")
        # CrossRef polite pool: include mailto in User-Agent
        ua = "GeneVariantFetcher/1.0 (https://github.com/KronckeLab/GeneVariantFetcher"
        if self._email:
            ua += f"; mailto:{self._email}"
        ua += ")"
        self.session.headers.update({"User-Agent": ua})
        self._last_request_time: float = 0.0

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def fetch(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Fetch CrossRef supplement metadata for a paper.

        Requires a DOI — returns [] immediately without one.
        """
        if not doi:
            logger.debug("No DOI for PMID %s, skipping CrossRef lookup", pmid)
            return []

        return self._fetch_from_crossref(doi, pmid)

    def list_supplements(self, pmid: str, doi: str = "") -> List[SupplementFile]:
        """Alias for fetch()."""
        return self.fetch(pmid, doi)

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    @api_retry
    def _fetch_from_crossref(self, doi: str, pmid: str) -> List[SupplementFile]:
        """Query CrossRef works endpoint and extract supplement links."""
        self._rate_limit()
        url = f"{self.API_BASE}/{doi}"

        try:
            response = self.session.get(url, timeout=self.timeout)
            if response.status_code == 404:
                logger.debug("DOI %s not found in CrossRef", doi)
                return []
            response.raise_for_status()
        except requests.exceptions.HTTPError:
            return []
        except Exception as e:
            logger.warning("CrossRef request failed for DOI %s: %s", doi, e)
            return []

        try:
            data = response.json()
        except ValueError:
            return []

        message = data.get("message", {})
        return self._extract_supplements(message, doi, pmid)

    def _extract_supplements(
        self, message: dict, doi: str, pmid: str
    ) -> List[SupplementFile]:
        """Extract supplement URLs from a CrossRef works response message."""
        results: List[SupplementFile] = []
        seen_urls: set = set()

        # Strategy 1: Parse message["link"] array
        for link in message.get("link", []):
            intended = link.get("intended-application", "").lower()
            if intended in _SKIP_APPLICATIONS:
                continue

            href = link.get("URL", "")
            if not href:
                continue

            if self._is_supplement_url(href):
                norm = href.split("#")[0].rstrip("/")
                if norm in seen_urls:
                    continue
                seen_urls.add(norm)

                content_type = link.get("content-type", "")
                filename = Path(urlparse(href).path).name or "supplement"
                results.append(
                    SupplementFile(
                        url=href,
                        name=self._clean_filename(filename),
                        source="crossref",
                        pmid=pmid,
                        mime_type=content_type,
                    )
                )

        # Strategy 2: Parse message["relation"]["has-supplement"]
        relations = message.get("relation", {})
        for supp_rel in relations.get("has-supplement", []):
            if supp_rel.get("id-type", "") != "uri":
                continue
            href = supp_rel.get("id", "")
            if not href:
                continue

            norm = href.split("#")[0].rstrip("/")
            if norm in seen_urls:
                continue
            seen_urls.add(norm)

            filename = Path(urlparse(href).path).name or "supplement"
            results.append(
                SupplementFile(
                    url=href,
                    name=self._clean_filename(filename),
                    source="crossref",
                    pmid=pmid,
                )
            )

        if results:
            logger.info(
                "CrossRef: found %d supplement(s) for DOI %s", len(results), doi
            )

        return results

    def _is_supplement_url(self, url: str) -> bool:
        """Check if a URL looks like it points to supplementary material."""
        url_lower = url.lower()

        # Check for supplement keywords in URL path
        for kw in _SUPPLEMENT_KEYWORDS:
            if kw in url_lower:
                return True

        # Check for valid file extension (data files are likely supplements)
        if self._has_valid_extension(url):
            path = urlparse(url).path.lower()
            # Exclude the main article PDF (often the only non-supplement link)
            if not path.endswith(".pdf") or any(
                kw in path for kw in ("suppl", "mmc", "appendix", "table_s", "si_")
            ):
                return True

        return False

    def _rate_limit(self) -> None:
        """Enforce polite rate limiting (min 0.1s between requests)."""
        now = time.time()
        elapsed = now - self._last_request_time
        if elapsed < 0.1:
            time.sleep(0.1 - elapsed)
        self._last_request_time = time.time()
