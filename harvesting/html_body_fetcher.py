"""
HTML body fetcher — multi-source fallback for abstract-only / paywall-scrape PMIDs.

When the standard harvest pipeline (PMC XML → DOI scrape → browser fallback) gives
up and writes an `# ABSTRACT-ONLY FALLBACK` FULL_CONTEXT.md, this module tries a
fresh set of OA-friendly sources to recover real body text:

    1. Europe PMC fullTextXML       (independent index, often has BMC/Wellcome/etc.)
    2. NCBI ELink → PMC OA          (find PMCID even when PubMed metadata hides it)
    3. Unpaywall OA HTML landing    (preprint / accepted-MS HTML hosted by author/archive)
    4. Unpaywall OA PDF             (last resort — PDF → markdown via FormatConverter)

Each candidate is gated by `is_acceptable_body()`:
    - ≥ MIN_BODY_LINES non-heading lines
    - ≥ MIN_BODY_CHARS characters
    - no paywall sentinel ("Purchase this article", "subscription content, log in",
      "Check if you have access", etc.)

Designed as a sidecar: results overwrite FULL_CONTEXT.md in place, leaving the
rest of the pipeline (extraction, aggregation, DB build) untouched.
"""

from __future__ import annotations

import logging
import re
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import requests

from harvesting.format_converters import FormatConverter
from harvesting.pmc_api import PMCAPIClient
from harvesting.supplement_scraper import SupplementScraper
from harvesting.unpaywall_api import UnpaywallClient

logger = logging.getLogger(__name__)


# Quality gates — calibrated so a typical journal article (Methods + Results +
# Discussion) clears them but a landing-page abstract does not. Body lines and
# char count are checked as alternatives because many JATS XML papers render
# whole paragraphs as a single long line (low line count, plenty of content).
MIN_BODY_LINES = 15
MIN_BODY_CHARS = 4000
GENEROUS_CHARS = 8000  # If above this, line-count gate is bypassed

PAYWALL_SENTINELS = (
    "purchase this article",
    "subscription content, log in",
    "check if you have access through your login",
    "to read this article in full",
    "this is a preview of subscription content",
    "log in to read full text",
    "buy this article",
    "rent this article",
    "get access",
    "purchase pdf",
    "this article requires a subscription",
)

USER_AGENT = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/120.0.0.0 Safari/537.36 "
    "GeneVariantFetcher/2.1 (mailto:{email})"
)


@dataclass
class FetchResult:
    """Result of an HTML-body fetch attempt."""

    success: bool
    markdown: Optional[str] = None
    source: Optional[str] = None  # which sub-fetcher won
    n_body_lines: int = 0
    n_chars: int = 0
    attempts: list = None  # list of (source, status, detail)
    error: Optional[str] = None

    def __post_init__(self):
        if self.attempts is None:
            self.attempts = []


def count_body_lines(markdown: str) -> int:
    """Lines that are neither blank nor a markdown heading."""
    n = 0
    for ln in markdown.splitlines():
        s = ln.strip()
        if not s or s.startswith("#") or s.startswith(">"):
            continue
        n += 1
    return n


def has_paywall_sentinel(markdown: str) -> bool:
    low = markdown.lower()
    return any(s in low for s in PAYWALL_SENTINELS)


def is_acceptable_body(markdown: Optional[str]) -> tuple[bool, str]:
    """Quality gate. Returns (acceptable, reason)."""
    if not markdown:
        return False, "empty markdown"
    if has_paywall_sentinel(markdown):
        return False, "paywall sentinel detected"
    if len(markdown) < MIN_BODY_CHARS:
        return False, f"too short ({len(markdown)} < {MIN_BODY_CHARS} chars)"
    n_lines = count_body_lines(markdown)
    # Long-paragraph JATS papers can have few "lines" but lots of chars; accept
    # them on char count alone when body is generous.
    if len(markdown) >= GENEROUS_CHARS:
        return True, "ok"
    if n_lines < MIN_BODY_LINES:
        return False, f"too few body lines ({n_lines} < {MIN_BODY_LINES})"
    return True, "ok"


class HtmlBodyFetcher:
    """Multi-source body fetcher. Stateless across calls; reuse via a single instance."""

    EPMC_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
    NCBI_ELINK = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

    def __init__(
        self,
        email: str,
        ncbi_api_key: Optional[str] = None,
        session: Optional[requests.Session] = None,
        min_request_interval: float = 0.34,  # NCBI is 3 req/s with key; play safe
    ):
        self.email = email
        self.ncbi_api_key = ncbi_api_key
        self.session = session or requests.Session()
        self.session.headers.update({"User-Agent": USER_AGENT.format(email=email)})

        self._last_request_time = 0.0
        self._min_interval = min_request_interval

        # Reuse existing clients where possible
        self.pmc = PMCAPIClient(session=self.session)
        self.unpaywall = UnpaywallClient(email=email, session=self.session)
        self.scraper = SupplementScraper()
        self.converter = FormatConverter()

    # -------------------------------------------------------------------------
    # rate limiting
    # -------------------------------------------------------------------------
    def _rate_limit(self):
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_request_time = time.time()

    # -------------------------------------------------------------------------
    # public entry point
    # -------------------------------------------------------------------------
    def fetch_body(self, pmid: str, doi: Optional[str] = None) -> FetchResult:
        """
        Try every available source in order. Returns the first acceptable result.

        Args:
            pmid: PubMed ID
            doi:  Optional DOI hint; if missing, will be looked up via NCBI

        Returns:
            FetchResult with `success=True` and `markdown` set if any source
            cleared the quality gate.
        """
        attempts: list = []

        # Look up DOI if not provided (also seeds the NCBI ELink step)
        if not doi:
            try:
                doi = self.pmc.get_doi_from_pmid(pmid)
            except Exception as exc:
                logger.debug(f"PMID {pmid}: DOI lookup failed: {exc}")

        # Source 1: Europe PMC fullTextXML
        md, detail = self._try_europepmc_xml(pmid)
        attempts.append(("europepmc_xml", "ok" if md else "miss", detail))
        ok, reason = is_acceptable_body(md)
        if ok:
            return FetchResult(
                success=True,
                markdown=md,
                source="europepmc_xml",
                n_body_lines=count_body_lines(md),
                n_chars=len(md),
                attempts=attempts,
            )
        elif md:
            attempts[-1] = ("europepmc_xml", f"rejected: {reason}", detail)

        # Source 2: NCBI ELink → PMC OA XML (find PMCID even when PubMed hides it)
        md, detail = self._try_ncbi_elink_pmc(pmid)
        attempts.append(("ncbi_elink_pmc", "ok" if md else "miss", detail))
        ok, reason = is_acceptable_body(md)
        if ok:
            return FetchResult(
                success=True,
                markdown=md,
                source="ncbi_elink_pmc",
                n_body_lines=count_body_lines(md),
                n_chars=len(md),
                attempts=attempts,
            )
        elif md:
            attempts[-1] = ("ncbi_elink_pmc", f"rejected: {reason}", detail)

        # Source 3 + 4: Unpaywall (HTML landing first, then PDF)
        if doi:
            md, detail = self._try_unpaywall_html(doi)
            attempts.append(("unpaywall_html", "ok" if md else "miss", detail))
            ok, reason = is_acceptable_body(md)
            if ok:
                return FetchResult(
                    success=True,
                    markdown=md,
                    source="unpaywall_html",
                    n_body_lines=count_body_lines(md),
                    n_chars=len(md),
                    attempts=attempts,
                )
            elif md:
                attempts[-1] = ("unpaywall_html", f"rejected: {reason}", detail)

            md, detail = self._try_unpaywall_pdf(doi)
            attempts.append(("unpaywall_pdf", "ok" if md else "miss", detail))
            ok, reason = is_acceptable_body(md)
            if ok:
                return FetchResult(
                    success=True,
                    markdown=md,
                    source="unpaywall_pdf",
                    n_body_lines=count_body_lines(md),
                    n_chars=len(md),
                    attempts=attempts,
                )
            elif md:
                attempts[-1] = ("unpaywall_pdf", f"rejected: {reason}", detail)
        else:
            attempts.append(("unpaywall_html", "skip", "no DOI"))
            attempts.append(("unpaywall_pdf", "skip", "no DOI"))

        return FetchResult(
            success=False, attempts=attempts, error="no source produced acceptable body"
        )

    # -------------------------------------------------------------------------
    # Source 1: Europe PMC fullTextXML
    # -------------------------------------------------------------------------
    def _try_europepmc_xml(self, pmid: str) -> tuple[Optional[str], str]:
        """Look up PMID at Europe PMC, fetch JATS XML if available, convert to markdown."""
        self._rate_limit()
        try:
            search_url = f"{self.EPMC_BASE}/search"
            params = {
                "query": f"ext_id:{pmid} AND src:MED",
                "resultType": "core",
                "format": "json",
                "pageSize": 1,
            }
            r = self.session.get(search_url, params=params, timeout=30)
            r.raise_for_status()
            data = r.json()
        except Exception as exc:
            return None, f"search error: {exc}"

        results = (data.get("resultList") or {}).get("result") or []
        if not results:
            return None, "not found at europepmc"
        rec = results[0]

        src = rec.get("source")
        pmcid = rec.get("pmcid")
        has_xml = (
            rec.get("hasTextMinedTerms") == "Y"
            or rec.get("inEPMC") == "Y"
            or rec.get("isOpenAccess") == "Y"
        )

        # Decide which id Europe PMC will accept for fullTextXML
        # For PMC-hosted articles: src=PMC, id=PMCID (numeric)
        # For preprints: src=PPR, id=full id
        # For MEDLINE-only OA: src=MED, id=PMID (rarely works)
        candidates = []
        if pmcid:
            num = pmcid.replace("PMC", "").strip()
            if num:
                candidates.append(("PMC", num))
        if src == "PPR" and rec.get("id"):
            candidates.append((src, rec["id"]))
        # Always try MED:PMID as last fallback
        candidates.append(("MED", pmid))

        last_detail = "no candidate worked"
        for cand_src, cand_id in candidates:
            self._rate_limit()
            try:
                xml_url = f"{self.EPMC_BASE}/{cand_src}/{cand_id}/fullTextXML"
                r = self.session.get(xml_url, timeout=60)
                if (
                    r.status_code != 200
                    or not r.text
                    or "<error" in r.text.lower()[:200]
                ):
                    last_detail = f"{cand_src}/{cand_id} → HTTP {r.status_code}, len={len(r.text)}"
                    continue
                xml_content = r.text
                if len(xml_content) < 1000:
                    last_detail = (
                        f"{cand_src}/{cand_id} → tiny XML ({len(xml_content)} bytes)"
                    )
                    continue
                markdown = self.converter.xml_to_markdown(xml_content)
                if markdown and len(markdown) > 500:
                    return (
                        markdown,
                        f"epmc {cand_src}/{cand_id} → {len(markdown)} chars",
                    )
                last_detail = f"{cand_src}/{cand_id} → xml_to_markdown empty"
            except Exception as exc:
                last_detail = f"{cand_src}/{cand_id} → {exc}"
                continue

        return None, last_detail

    # -------------------------------------------------------------------------
    # Source 2: NCBI ELink → PMC
    # -------------------------------------------------------------------------
    def _try_ncbi_elink_pmc(self, pmid: str) -> tuple[Optional[str], str]:
        """
        Use NCBI ELink to find a PMC entry for this PMID. This catches articles
        where PubMed's per-record `ArticleIdList` does NOT list a PMCID but
        a PMC entry exists (common for delayed-deposit OA articles).
        """
        self._rate_limit()
        try:
            params = {
                "dbfrom": "pubmed",
                "db": "pmc",
                "id": pmid,
                "tool": "GeneVariantFetcher",
                "email": self.email,
            }
            if self.ncbi_api_key:
                params["api_key"] = self.ncbi_api_key
            r = self.session.get(self.NCBI_ELINK, params=params, timeout=30)
            r.raise_for_status()
        except Exception as exc:
            return None, f"elink error: {exc}"

        # ELink returns multiple LinkSetDb blocks. We want ONLY the one whose
        # <LinkName> is `pubmed_pmc` (the direct PMID→PMC mapping). Other blocks
        # (e.g. `pubmed_pmc_refs`, `pubmed_pmc_local`) point at related/citing
        # papers — picking the first <Id> in the document yields wrong content.
        text_resp = r.text
        pmc_block_match = re.search(
            r"<LinkSetDb>\s*<DbTo>pmc</DbTo>\s*<LinkName>pubmed_pmc</LinkName>(.*?)</LinkSetDb>",
            text_resp,
            re.DOTALL,
        )
        if not pmc_block_match:
            return None, "elink: no pubmed_pmc linkset for this PMID"
        id_match = re.search(r"<Link>\s*<Id>(\d+)</Id>", pmc_block_match.group(1))
        if not id_match:
            return None, "elink: pubmed_pmc linkset empty"

        numeric_pmcid = id_match.group(1)
        pmcid = f"PMC{numeric_pmcid}"

        # Try the existing PMC client first (handles BioC + JATS retries)
        try:
            xml_content = self.pmc.get_fulltext_xml(pmcid)
            if xml_content:
                markdown = self.converter.xml_to_markdown(xml_content)
                if markdown and len(markdown) > 500:
                    return markdown, f"elink → {pmcid} → {len(markdown)} chars"
        except Exception as exc:
            return None, f"elink {pmcid}: pmc xml fetch error: {exc}"

        return None, f"elink {pmcid}: xml empty/short"

    # -------------------------------------------------------------------------
    # Source 3: Unpaywall HTML landing
    # -------------------------------------------------------------------------
    def _try_unpaywall_html(self, doi: str) -> tuple[Optional[str], str]:
        """
        Pull all OA locations Unpaywall knows about, fetch the *HTML* landing
        pages (skip pure PDF locations), and run each through the generic
        scraper. Return first acceptable markdown.
        """
        self._rate_limit()
        result, err = self.unpaywall.find_open_access(doi)
        if not result or not result.get("is_oa"):
            return None, f"unpaywall not OA: {err or 'no OA location'}"

        # Walk every OA location, not just the "best" one — author-archive
        # HTML often has fuller text than the publisher landing.
        locations = []
        best = result.get("best_oa_location")
        if best:
            locations.append(best)
        # All other locations
        # find_open_access only returns "best_oa_location" — we need the raw
        # response. Refetch the raw record once to enumerate.
        raw = self._unpaywall_raw(doi)
        if raw:
            for loc in raw.get("oa_locations", []) or []:
                if loc not in locations:
                    locations.append(loc)

        if not locations:
            return None, "unpaywall: no OA locations"

        last_detail = "no HTML location worked"
        for loc in locations:
            url = loc.get("url_for_landing_page") or loc.get("url")
            if not url:
                continue
            # Skip explicit PDF locations (handled by _try_unpaywall_pdf)
            if loc.get("url_for_pdf") and url == loc.get("url_for_pdf"):
                continue
            if url.lower().endswith(".pdf"):
                continue

            self._rate_limit()
            try:
                hr = self.session.get(url, timeout=60, allow_redirects=True)
                if hr.status_code != 200:
                    last_detail = f"{url[:60]}… → HTTP {hr.status_code}"
                    continue
                ctype = hr.headers.get("content-type", "").lower()
                if "html" not in ctype:
                    last_detail = f"{url[:60]}… → content-type {ctype}"
                    continue
                markdown, _title = self.scraper.extract_fulltext(hr.text, hr.url)
                if markdown and len(markdown) > 500:
                    return (
                        markdown,
                        f"unpaywall html {url[:60]}… → {len(markdown)} chars",
                    )
                last_detail = f"{url[:60]}… → empty extraction"
            except Exception as exc:
                last_detail = f"{url[:60]}… → {exc}"
                continue

        return None, last_detail

    def _unpaywall_raw(self, doi: str) -> Optional[dict]:
        """Direct Unpaywall fetch to get the full oa_locations list."""
        self._rate_limit()
        try:
            from urllib.parse import quote

            encoded = quote(doi.strip(), safe="/")
            url = f"https://api.unpaywall.org/v2/{encoded}"
            r = self.session.get(url, params={"email": self.email}, timeout=30)
            if r.status_code == 200 and r.text.strip():
                return r.json()
        except Exception:
            pass
        return None

    # -------------------------------------------------------------------------
    # Source 4: Unpaywall PDF → markdown
    # -------------------------------------------------------------------------
    def _try_unpaywall_pdf(self, doi: str) -> tuple[Optional[str], str]:
        """Fetch the best OA PDF and convert via FormatConverter.pdf_to_markdown."""
        self._rate_limit()
        result, err = self.unpaywall.find_open_access(doi)
        if not result or not result.get("is_oa"):
            return None, f"unpaywall not OA: {err or 'no OA'}"

        pdf_url = result.get("pdf_url")
        if not pdf_url:
            return None, "unpaywall: no PDF URL"

        self._rate_limit()
        try:
            hr = self.session.get(pdf_url, timeout=120, allow_redirects=True)
            if hr.status_code != 200:
                return None, f"pdf HTTP {hr.status_code}"
            ctype = hr.headers.get("content-type", "").lower()
            if "pdf" not in ctype and not pdf_url.lower().endswith(".pdf"):
                return None, f"not a PDF (content-type {ctype})"

            # Write to a temp file so the existing converter can stream it
            with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
                tmp.write(hr.content)
                tmp_path = Path(tmp.name)
            try:
                markdown = self.converter.pdf_to_markdown(tmp_path)
            finally:
                try:
                    tmp_path.unlink(missing_ok=True)
                except Exception:
                    pass

            if markdown and len(markdown) > 500:
                return markdown, f"unpaywall pdf → {len(markdown)} chars"
            return None, "pdf_to_markdown empty/short"
        except Exception as exc:
            return None, f"pdf fetch/convert error: {exc}"


# -----------------------------------------------------------------------------
# Convenience: classify an existing FULL_CONTEXT.md as needing re-harvest
# -----------------------------------------------------------------------------


def needs_reharvest(full_context_path: Path) -> tuple[bool, str]:
    """
    Decide whether a FULL_CONTEXT.md needs re-harvest.

    Returns (needs_reharvest, reason).
    """
    if not full_context_path.exists():
        return True, "missing"
    try:
        text = full_context_path.read_text(errors="ignore")
    except Exception as exc:
        return True, f"unreadable: {exc}"

    if "# ABSTRACT-ONLY FALLBACK" in text[:500]:
        return True, "abstract_only_fallback"
    if "# REHARVESTED FULL TEXT" in text[:200]:
        # Previously recovered — leave it alone to avoid undoing prior work.
        return False, "already_reharvested"
    size = full_context_path.stat().st_size
    if size < 4000:
        return True, f"tiny ({size} bytes)"
    if has_paywall_sentinel(text):
        return True, "paywall_sentinel"
    # Use the same gate as is_acceptable_body so we only re-harvest files that
    # would NOT pass the gate today.
    n_lines = count_body_lines(text)
    if len(text) < GENEROUS_CHARS and n_lines < MIN_BODY_LINES:
        return True, f"too few body lines ({n_lines})"
    return False, "ok"
