"""Client for the Unpaywall API to find open access versions of papers.

Unpaywall (https://unpaywall.org) indexes legal OA copies of scholarly articles.
The API is free and requires only an email address for identification.
"""

from __future__ import annotations

import json
import logging
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

UNPAYWALL_BASE = "https://api.unpaywall.org/v2"
EMAIL = "brett.kroncke@gmail.com"


class UnpaywallError(Exception):
    """Raised when communication with the Unpaywall API fails."""


class UnpaywallClient:
    """Client for querying the Unpaywall REST API.

    Parameters
    ----------
    email : str
        Email address used for API identification (required by Unpaywall TOS).
    timeout : float
        HTTP request timeout in seconds.
    """

    def __init__(
        self,
        email: str = EMAIL,
        *,
        timeout: float = 15.0,
    ) -> None:
        self.email = email
        self.timeout = timeout
        self._cache: Dict[str, dict] = {}
        self._last_request_time: float = 0.0

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def find_oa_url(self, doi: str) -> Optional[str]:
        """Return the best open-access URL for *doi*, or ``None``.

        Prefers the ``best_oa_location`` field returned by Unpaywall, falling
        back to the first entry in ``oa_locations`` if needed.
        """
        data = self._fetch(doi)
        if data is None:
            return None

        best = data.get("best_oa_location")
        if best:
            url = best.get("url_for_pdf") or best.get("url_for_landing_page") or best.get("url")
            if url:
                return url

        # Fallback: iterate oa_locations
        for loc in data.get("oa_locations", []):
            url = loc.get("url_for_pdf") or loc.get("url_for_landing_page") or loc.get("url")
            if url:
                return url

        return None

    def get_oa_locations(self, doi: str) -> List[dict]:
        """Return all OA locations for *doi* as a list of dicts.

        Each dict contains keys like ``url``, ``url_for_pdf``,
        ``url_for_landing_page``, ``host_type``, ``version``, etc.
        Returns an empty list when the DOI is not found or has no OA copies.
        """
        data = self._fetch(doi)
        if data is None:
            return []
        return data.get("oa_locations", [])

    def is_open_access(self, doi: str) -> bool:
        """Return ``True`` if *doi* has at least one open-access copy."""
        data = self._fetch(doi)
        if data is None:
            return False
        return bool(data.get("is_oa"))

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _fetch(self, doi: str) -> Optional[dict]:
        """Fetch and cache the Unpaywall record for *doi*."""
        doi = doi.strip()
        if doi in self._cache:
            return self._cache[doi]

        self._rate_limit()

        url = f"{UNPAYWALL_BASE}/{urllib.parse.quote(doi, safe='')}?email={urllib.parse.quote(self.email)}"
        req = urllib.request.Request(url)
        req.add_header("User-Agent", "GeneVariantFetcher/0.1")

        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                data = json.loads(resp.read().decode())
            self._cache[doi] = data
            logger.debug("Unpaywall hit for DOI %s – is_oa=%s", doi, data.get("is_oa"))
            return data
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                logger.debug("Unpaywall: DOI not found – %s", doi)
                self._cache[doi] = None  # type: ignore[assignment]
                return None
            logger.warning("Unpaywall HTTP %s for DOI %s: %s", exc.code, doi, exc.reason)
            raise UnpaywallError(f"HTTP {exc.code}: {exc.reason}") from exc
        except urllib.error.URLError as exc:
            logger.warning("Unpaywall network error for DOI %s: %s", doi, exc.reason)
            raise UnpaywallError(str(exc.reason)) from exc

    def _rate_limit(self) -> None:
        """Enforce a minimum 1-second gap between API requests."""
        now = time.monotonic()
        elapsed = now - self._last_request_time
        if elapsed < 1.0:
            time.sleep(1.0 - elapsed)
        self._last_request_time = time.monotonic()
