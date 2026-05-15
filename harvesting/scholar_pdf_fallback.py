"""Google Scholar title-based PDF recovery.

This is a last-resort OA path for author/lab-hosted PDFs that are not indexed
by PMC, Europe PMC, or Unpaywall. It intentionally does not try to bypass
Scholar CAPTCHA or bot checks; those are reported as clean misses.
"""

from __future__ import annotations

import logging
import os
import re
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional
from urllib.parse import urlparse

import requests
from bs4 import BeautifulSoup

from harvesting.format_converters import FormatConverter

logger = logging.getLogger(__name__)

SCHOLAR_URL = "https://scholar.google.com/scholar"
PUBMED_SUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
DEFAULT_SCHOLAR_MIN_REQUEST_INTERVAL = 10.0
DEFAULT_SCHOLAR_BLOCK_BACKOFF_SECONDS = 15 * 60
MAX_PDF_CANDIDATES = 4

BLOCK_MARKERS = (
    "unusual traffic",
    "not a robot",
    "recaptcha",
    "captcha",
    "sorry/index",
    "please show you're not a robot",
    "automated queries",
)

BLOCKED_PDF_HOST_MARKERS = (
    "sci-hub",
    "researchgate.net",
    "academia.edu",
)

CONVERTER_FAILURE_PREFIXES = (
    "[PDF file available",
    "[Invalid PDF file",
    "[Error reading PDF",
)

MIN_PDF_BYTES = 50_000

DEFAULT_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/131.0.0.0 Safari/537.36"
    ),
    "Accept": "text/html,application/xhtml+xml,application/pdf;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.9",
}

_THROTTLE_LOCK = threading.Lock()
_LAST_REQUEST = 0.0
_SCHOLAR_BLOCKED_UNTIL = 0.0


@dataclass
class ScholarPDFResult:
    success: bool
    markdown: Optional[str] = None
    source_url: Optional[str] = None
    title: Optional[str] = None
    attempts: list[tuple[str, str, str]] = field(default_factory=list)
    error: Optional[str] = None

    @property
    def detail(self) -> str:
        if self.success and self.source_url:
            return self.source_url
        return self.error or ""


class ScholarPDFRecovery:
    """Recover a likely article PDF from Google Scholar's title search."""

    def __init__(
        self,
        *,
        session: requests.Session,
        converter: FormatConverter,
        email: str,
        quality_gate: Callable[[str], tuple[bool, str]],
        api_key: Optional[str] = None,
        min_request_interval: float = DEFAULT_SCHOLAR_MIN_REQUEST_INTERVAL,
    ):
        self.session = session
        self.converter = converter
        self.email = email
        self.quality_gate = quality_gate
        self.api_key = api_key
        self.min_request_interval = min_request_interval

    def recover(
        self,
        *,
        pmid: str,
        title: Optional[str] = None,
    ) -> ScholarPDFResult:
        attempts: list[tuple[str, str, str]] = []
        clean_title = _clean_title(title) if title else self._pubmed_title(pmid)
        if not clean_title:
            return ScholarPDFResult(
                success=False,
                attempts=attempts,
                error="no title available for Scholar search",
            )

        pdf_urls, detail = self._scholar_pdf_urls(clean_title)
        attempts.append(("scholar_search", "ok" if pdf_urls else "miss", detail))
        if not pdf_urls:
            return ScholarPDFResult(
                success=False,
                title=clean_title,
                attempts=attempts,
                error=detail,
            )

        for url in pdf_urls[:MAX_PDF_CANDIDATES]:
            markdown, detail = self._pdf_url_to_markdown(url)
            ok, reason = self.quality_gate(markdown or "")
            if markdown and ok:
                attempts.append(("scholar_pdf", "ok", f"{url}: {reason}"))
                return ScholarPDFResult(
                    success=True,
                    markdown=markdown,
                    source_url=url,
                    title=clean_title,
                    attempts=attempts,
                )
            status = f"reject: {reason}" if markdown else "miss"
            attempts.append(("scholar_pdf", status, f"{url}: {detail}"))

        return ScholarPDFResult(
            success=False,
            title=clean_title,
            attempts=attempts,
            error="no Scholar PDF produced acceptable body",
        )

    def _throttle(self) -> None:
        global _LAST_REQUEST
        with _THROTTLE_LOCK:
            elapsed = time.time() - _LAST_REQUEST
            if elapsed < self.min_request_interval:
                time.sleep(self.min_request_interval - elapsed)
            _LAST_REQUEST = time.time()

    def _mark_blocked(self) -> None:
        global _SCHOLAR_BLOCKED_UNTIL
        if DEFAULT_SCHOLAR_BLOCK_BACKOFF_SECONDS <= 0:
            return
        with _THROTTLE_LOCK:
            _SCHOLAR_BLOCKED_UNTIL = max(
                _SCHOLAR_BLOCKED_UNTIL,
                time.time() + DEFAULT_SCHOLAR_BLOCK_BACKOFF_SECONDS,
            )

    def _pubmed_title(self, pmid: str) -> Optional[str]:
        try:
            params = {
                "db": "pubmed",
                "id": pmid,
                "retmode": "json",
                "tool": "GeneVariantFetcher",
                "email": self.email,
            }
            if self.api_key:
                params["api_key"] = self.api_key
            response = self.session.get(
                PUBMED_SUMMARY_URL,
                params=params,
                timeout=20,
            )
            if response.status_code != 200:
                return None
            item = response.json().get("result", {}).get(str(pmid), {})
            return _clean_title(item.get("title"))
        except Exception as exc:
            logger.debug("PubMed title lookup failed for PMID %s: %s", pmid, exc)
            return None

    def _scholar_pdf_urls(self, title: str) -> tuple[list[str], str]:
        with _THROTTLE_LOCK:
            blocked_remaining = _SCHOLAR_BLOCKED_UNTIL - time.time()
        if blocked_remaining > 0:
            return (
                [],
                f"Scholar temporarily disabled after block ({blocked_remaining:.0f}s remaining)",
            )

        self._throttle()
        try:
            response = self.session.get(
                SCHOLAR_URL,
                params={"q": f'"{title}"', "hl": "en"},
                headers={**DEFAULT_HEADERS, "Referer": "https://scholar.google.com/"},
                timeout=30,
                allow_redirects=True,
            )
        except Exception as exc:
            return [], f"Scholar search error: {exc}"

        text = response.text or ""
        if response.status_code in (403, 429):
            self._mark_blocked()
            return [], f"Scholar blocked HTTP {response.status_code}"
        if _looks_blocked(text, response.url):
            self._mark_blocked()
            return [], "Scholar captcha/block page"
        if response.status_code != 200:
            return [], f"Scholar HTTP {response.status_code}"

        urls = _extract_pdf_urls(text)
        return urls, f"{len(urls)} PDF candidate(s)"

    def _pdf_url_to_markdown(self, url: str) -> tuple[Optional[str], str]:
        try:
            response = self.session.get(
                url,
                headers={**DEFAULT_HEADERS, "Referer": SCHOLAR_URL},
                timeout=90,
                allow_redirects=True,
            )
        except Exception as exc:
            return None, f"PDF fetch error: {exc}"

        if response.status_code != 200:
            return None, f"PDF HTTP {response.status_code}"
        body = response.content or b""
        content_length = response.headers.get("content-length")
        if (
            content_length
            and content_length.isdigit()
            and int(content_length) < MIN_PDF_BYTES
        ):
            return None, f"PDF too small ({content_length} bytes)"
        if len(body) < MIN_PDF_BYTES:
            return None, f"PDF too small ({len(body)} bytes)"
        content_type = (response.headers.get("content-type") or "").lower()
        if not (body.startswith(b"%PDF-") or "pdf" in content_type):
            return None, f"not PDF content ({content_type or 'unknown content-type'})"

        with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
            tmp.write(body)
            tmp_path = Path(tmp.name)
        try:
            markdown = self.converter.pdf_to_markdown(tmp_path)
            if not markdown:
                return None, "PDF conversion returned empty markdown"
            if markdown.startswith(CONVERTER_FAILURE_PREFIXES):
                return None, "PDF conversion failed"
            markdown = re.sub(r"(?<!\n)\n(?!\n)", " ", markdown)
            return markdown, f"{len(markdown)} chars from PDF"
        except Exception as exc:
            return None, f"PDF conversion error: {exc}"
        finally:
            try:
                tmp_path.unlink(missing_ok=True)
            except Exception:
                pass


def _clean_title(title: Optional[str]) -> Optional[str]:
    if not title:
        return None
    cleaned = re.sub(r"\s+", " ", str(title)).strip()
    cleaned = cleaned.strip(" .")
    return cleaned or None


def try_scholar_pdf(
    *,
    title: Optional[str],
    pmid: str,
    session: requests.Session,
    converter: FormatConverter,
    quality_gate: Callable[[str], tuple[bool, str]],
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) -> ScholarPDFResult:
    """Convenience wrapper used by the existing recovery flows."""
    client = ScholarPDFRecovery(
        session=session,
        converter=converter,
        email=email
        or os.environ.get("NCBI_EMAIL")
        or os.environ.get("ENTREZ_EMAIL")
        or "gvf@example.com",
        quality_gate=quality_gate,
        api_key=api_key
        or os.environ.get("NCBI_API_KEY")
        or os.environ.get("ENTREZ_API_KEY"),
    )
    return client.recover(pmid=pmid, title=title)


def _looks_blocked(html: str, final_url: str) -> bool:
    low = f"{final_url}\n{html[:8000]}".lower()
    return any(marker in low for marker in BLOCK_MARKERS)


def _extract_pdf_urls(html: str) -> list[str]:
    soup = BeautifulSoup(html, "html.parser")
    urls: list[str] = []
    seen: set[str] = set()

    for anchor in soup.select("a[href]"):
        href = anchor.get("href") or ""
        text = anchor.get_text(" ", strip=True).lower()
        if not _is_pdf_candidate(href, text):
            continue
        if href in seen:
            continue
        seen.add(href)
        urls.append(href)

    return urls


def _is_pdf_candidate(href: str, text: str) -> bool:
    if not href.startswith(("http://", "https://")):
        return False
    host = urlparse(href).netloc.lower()
    if "google." in host or host.endswith("gstatic.com"):
        return False
    if any(marker in host for marker in BLOCKED_PDF_HOST_MARKERS):
        return False
    href_low = href.lower()
    return (
        "[pdf]" in text
        or " pdf" in f" {text}"
        or href_low.endswith(".pdf")
        or ".pdf?" in href_low
        or "/pdf/" in href_low
    )
