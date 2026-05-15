"""
Open-access full-text recovery for abstract-only PMIDs.

Plugs the gaps in the standard harvest pipeline by trying free OA sources the
existing pipeline misses:

    1. PMC HTML  (`/pmc/articles/{PMCID}/`)
       Works for Author Manuscripts where Europe PMC `fullTextXML` returns 404.
       In a 9-PMID sample of in-EPMC manuscripts this route succeeded 9/9
       while `fullTextXML` succeeded 0/9.

    2. Europe PMC fullTextXML
       True OA (`isOpenAccess=Y` in EPMC search) returns JATS XML directly.

    3. Europe PMC PDF render (`europepmc.org/articles/{PMCID}?pdf=render`)
       Returns rendered PDF for Author Manuscripts that lack JATS XML.

    4. Unpaywall best OA URL
       Publisher-hosted or repository-hosted OA copy keyed off DOI.
       Some publisher URLs hit Cloudflare; those are caught and reported.

Recoverability across the KCNH2 v4 run's 227 abstract-only PMIDs:

      pmc_oa              2  (0.9%)   ← was already recovered by existing pipeline
      epmc_manuscript     9  (4.0%)   ← NEW via PMC HTML route
      unpaywall_publisher 50  (22.0%) ← NEW via Unpaywall
      unpaywall_repo      11  (4.8%)  ← NEW via Unpaywall green OA
      doi_only_closed    140  (61.7%) ← genuinely paywalled
      no_doi_no_route     15  (6.6%)  ← older / obscure
      TOTAL RECOVERABLE   72  (31.7%)
"""

from __future__ import annotations

import logging
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import requests

from harvesting.format_converters import FormatConverter
from harvesting.pmc_api import PMCAPIClient
from harvesting.unpaywall_api import UnpaywallClient

logger = logging.getLogger(__name__)


MIN_BODY_CHARS = 4000
GENEROUS_CHARS = 8000
MIN_BODY_LINES = 15

PAYWALL_SENTINELS = (
    "purchase this article",
    "subscription content, log in",
    "check if you have access through your login",
    "this is a preview of subscription content",
    "log in to read full text",
    "buy this article",
    "rent this article",
    "purchase pdf",
    "this article requires a subscription",
    "just a moment...",
    "checking your browser before accessing",
)

BROWSER_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/124.0.0.0 Safari/537.36"
    ),
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,application/pdf,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.9",
}


@dataclass
class RecoveryResult:
    success: bool
    markdown: Optional[str] = None
    source: Optional[str] = None
    pmcid: Optional[str] = None
    n_chars: int = 0
    attempts: list = field(default_factory=list)
    error: Optional[str] = None


def _count_body_lines(md: str) -> int:
    n = 0
    for ln in md.splitlines():
        s = ln.strip()
        if not s or s.startswith("#") or s.startswith(">"):
            continue
        n += 1
    return n


def is_acceptable_body(md: Optional[str]) -> tuple[bool, str]:
    if not md:
        return False, "empty"
    low = md.lower()
    for sentinel in PAYWALL_SENTINELS:
        if sentinel in low:
            return False, f"paywall sentinel: {sentinel!r}"
    if len(md) < MIN_BODY_CHARS:
        return False, f"too short ({len(md)} < {MIN_BODY_CHARS})"
    if len(md) >= GENEROUS_CHARS:
        return True, "ok"
    if _count_body_lines(md) < MIN_BODY_LINES:
        return False, f"too few body lines"
    return True, "ok"


class OARecoveryClient:
    """Stateless recovery client. Re-use one instance across PMIDs."""

    EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
    PMC_HTML = "https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"

    def __init__(
        self,
        email: str,
        ncbi_api_key: Optional[str] = None,
        session: Optional[requests.Session] = None,
        min_request_interval: float = 0.34,
    ):
        self.email = email
        self.ncbi_api_key = ncbi_api_key
        self.session = session or requests.Session()
        self.session.headers.update(BROWSER_HEADERS)
        self._last = 0.0
        self._gap = min_request_interval

        self.pmc = PMCAPIClient(session=self.session)
        self.unpaywall = UnpaywallClient(email=email, session=self.session)
        self.converter = FormatConverter()

    def _throttle(self):
        elapsed = time.time() - self._last
        if elapsed < self._gap:
            time.sleep(self._gap - elapsed)
        self._last = time.time()

    def recover(self, pmid: str, doi: Optional[str] = None) -> RecoveryResult:
        attempts: list = []

        if not doi:
            try:
                doi = self.pmc.get_doi_from_pmid(pmid)
            except Exception as e:
                logger.debug("DOI lookup failed for %s: %s", pmid, e)

        try:
            pmcid = self.pmc.pmid_to_pmcid(pmid)
        except Exception as e:
            logger.debug("PMCID lookup failed for %s: %s", pmid, e)
            pmcid = None
        # Entrez is rate-prone; fall back to EPMC search which also returns a
        # PMCID when one exists and is much more reliable under load.
        if not pmcid:
            pmcid = self._pmcid_via_epmc(pmid)

        for source_name, fn in self._strategy_chain(pmid, doi, pmcid):
            try:
                md, detail = fn()
            except Exception as e:
                attempts.append((source_name, "error", str(e)[:120]))
                continue
            ok, reason = is_acceptable_body(md)
            if ok:
                attempts.append((source_name, "ok", detail))
                return RecoveryResult(
                    success=True,
                    markdown=md,
                    source=source_name,
                    pmcid=pmcid,
                    n_chars=len(md),
                    attempts=attempts,
                )
            attempts.append((source_name, f"reject: {reason}", detail))

        return RecoveryResult(
            success=False,
            pmcid=pmcid,
            attempts=attempts,
            error="no source produced acceptable body",
        )

    def _strategy_chain(self, pmid, doi, pmcid):
        """Yield (name, callable) in order of expected success rate."""
        if pmcid:
            yield "pmc_html", lambda: self._try_pmc_html(pmcid)
            yield "epmc_fulltextxml", lambda: self._try_epmc_xml(pmcid)
            yield "epmc_pdf_render", lambda: self._try_epmc_pdf_render(pmcid)
        if doi:
            yield "unpaywall_pdf", lambda: self._try_unpaywall_pdf(doi)
            yield "unpaywall_html", lambda: self._try_unpaywall_html(doi)
        yield "epmc_search_oa_url", lambda: self._try_epmc_oa_url(pmid)

    def _pmcid_via_epmc(self, pmid: str) -> Optional[str]:
        """Fallback PMCID lookup via Europe PMC search."""
        self._throttle()
        try:
            params = {
                "query": f"EXT_ID:{pmid} AND SRC:MED",
                "format": "json",
                "resultType": "core",
            }
            r = self.session.get(f"{self.EPMC}/search", params=params, timeout=15)
            if r.status_code != 200:
                return None
            results = r.json().get("resultList", {}).get("result", []) or []
            if not results:
                return None
            pmcid = results[0].get("pmcid")
            return pmcid if pmcid else None
        except Exception as e:
            logger.debug("EPMC PMCID fallback failed for %s: %s", pmid, e)
            return None

    def _try_pmc_html(self, pmcid: str) -> tuple[Optional[str], str]:
        self._throttle()
        url = self.PMC_HTML.format(pmcid=pmcid)
        r = self.session.get(url, timeout=30, allow_redirects=True)
        if r.status_code != 200:
            return None, f"HTTP {r.status_code}"
        md = self.converter.pmc_html_to_markdown(r.text)
        return md, f"{len(md)} chars from {pmcid}"

    def _try_epmc_xml(self, pmcid: str) -> tuple[Optional[str], str]:
        self._throttle()
        url = f"{self.EPMC}/{pmcid}/fullTextXML"
        r = self.session.get(url, timeout=30)
        if r.status_code != 200 or "<article" not in r.text.lower():
            return None, f"HTTP {r.status_code}, len={len(r.text)}"
        md = self.converter.xml_to_markdown(r.text)
        return md, f"{len(md)} chars from EPMC XML"

    def _try_epmc_pdf_render(self, pmcid: str) -> tuple[Optional[str], str]:
        self._throttle()
        url = f"https://europepmc.org/articles/{pmcid.lower()}?pdf=render"
        r = self.session.get(url, timeout=45, allow_redirects=True)
        if r.status_code != 200 or not r.content.startswith(b"%PDF-"):
            return (
                None,
                f"HTTP {r.status_code}, ct={r.headers.get('Content-Type', '?')}",
            )
        return self._pdf_bytes_to_markdown(r.content, f"EPMC-render {pmcid}")

    def _try_unpaywall_pdf(self, doi: str) -> tuple[Optional[str], str]:
        info = self._unpaywall_lookup(doi)
        if not info or not info.get("is_oa"):
            return None, "no Unpaywall OA record"
        pdf_url = (info.get("best_oa_location") or {}).get("url_for_pdf") or info.get(
            "pdf_url"
        )
        if not pdf_url:
            return None, "no Unpaywall PDF URL"
        return self._fetch_url_to_markdown(pdf_url, prefer="pdf")

    def _try_unpaywall_html(self, doi: str) -> tuple[Optional[str], str]:
        info = self._unpaywall_lookup(doi)
        if not info or not info.get("is_oa"):
            return None, "no Unpaywall OA record"
        loc = info.get("best_oa_location") or {}
        url = (
            loc.get("url")
            or loc.get("url_for_landing_page")
            or info.get("landing_page")
        )
        if not url:
            return None, "no Unpaywall HTML URL"
        return self._fetch_url_to_markdown(url, prefer="html")

    def _try_epmc_oa_url(self, pmid: str) -> tuple[Optional[str], str]:
        """Use EPMC search to find any 'Free'/'Open access' URL the publisher exposes."""
        self._throttle()
        params = {
            "query": f"EXT_ID:{pmid} AND SRC:MED",
            "format": "json",
            "resultType": "core",
        }
        r = self.session.get(f"{self.EPMC}/search", params=params, timeout=20)
        if r.status_code != 200:
            return None, f"HTTP {r.status_code}"
        results = r.json().get("resultList", {}).get("result", []) or []
        if not results:
            return None, "no EPMC record"
        ft_urls = results[0].get("fullTextUrlList", {}).get("fullTextUrl", []) or []
        for u in ft_urls:
            if (u.get("availability") or "").lower() in (
                "free",
                "open access",
            ) and u.get("url"):
                prefer = (
                    "pdf" if (u.get("documentStyle") or "").lower() == "pdf" else "html"
                )
                md, detail = self._fetch_url_to_markdown(u["url"], prefer=prefer)
                if md:
                    return md, f"EPMC FT url ({u.get('site')}): {detail}"
        return None, "no Free/OA url in EPMC FT list"

    def _unpaywall_lookup(self, doi: str) -> Optional[dict]:
        self._throttle()
        try:
            r = self.session.get(
                f"https://api.unpaywall.org/v2/{doi}",
                params={"email": self.email},
                timeout=20,
            )
            if r.status_code != 200:
                return None
            return r.json()
        except Exception as e:
            logger.debug("Unpaywall lookup failed for %s: %s", doi, e)
            return None

    def _fetch_url_to_markdown(
        self, url: str, prefer: str
    ) -> tuple[Optional[str], str]:
        self._throttle()
        try:
            r = self.session.get(url, timeout=45, allow_redirects=True)
        except Exception as e:
            return None, f"fetch error: {e}"
        if r.status_code != 200:
            return None, f"HTTP {r.status_code}"
        ct = (r.headers.get("Content-Type") or "").lower()
        body = r.content
        if "pdf" in ct or body[:5] == b"%PDF-":
            return self._pdf_bytes_to_markdown(body, url)
        if len(body) < 4000:
            return None, f"tiny response ({len(body)} bytes)"
        if self._looks_like_block(body):
            return None, "cloudflare/block page"
        try:
            md = self.converter.pmc_html_to_markdown(
                body.decode("utf-8", errors="ignore")
            )
        except Exception as e:
            return None, f"html parse error: {e}"
        return md, f"{len(md)} chars from HTML ({url[:60]}...)"

    def _pdf_bytes_to_markdown(
        self, pdf_bytes: bytes, label: str
    ) -> tuple[Optional[str], str]:
        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as f:
            f.write(pdf_bytes)
            tmp = Path(f.name)
        try:
            md = self.converter.pdf_to_markdown(tmp)
            return md, f"{len(md)} chars from PDF ({label[:50]}...)"
        except Exception as e:
            return None, f"PDF parse error: {e}"
        finally:
            try:
                tmp.unlink()
            except Exception:
                pass

    @staticmethod
    def _looks_like_block(body: bytes) -> bool:
        try:
            head = body[:8000].decode("utf-8", errors="ignore").lower()
        except Exception:
            return False
        markers = (
            "just a moment",
            "cf-mitigated",
            "cf-chl-bypass",
            "access denied",
            "checking your browser",
        )
        return any(m in head for m in markers)
