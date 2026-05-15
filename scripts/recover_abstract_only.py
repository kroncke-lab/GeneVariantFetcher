"""Recover full text for abstract-only papers via Unpaywall + preprint servers.

Scans a harvest directory for FULL_CONTEXT.md files that are abstract-only,
resolves DOIs (parsed from the abstract text, or via NCBI ESummary as
fallback), queries Unpaywall for OA copies, also probes bioRxiv/medRxiv,
downloads PDFs/HTML, converts to markdown, and rewrites the FULL_CONTEXT.md
with the recovered body. The original abstract-only file is preserved at
<PMID>_FULL_CONTEXT.ABSTRACT_ONLY.md.

Run:
    .venv/bin/python scripts/recover_abstract_only.py \\
        --harvest-dir results/KCNH2/20260506_102238/pmc_fulltext \\
        --email brett.kroncke@gmail.com \\
        --out-report recovery_report.json \\
        [--limit N] [--dry-run]
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import urlparse, quote

import requests

# Project imports
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from harvesting.unpaywall_api import UnpaywallClient  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s %(message)s",
)
log = logging.getLogger("recover")

ABSTRACT_MARKER = "ABSTRACT-ONLY FALLBACK"
ABSTRACT_FILE_SIZE_THRESHOLD = 6000

# PubMed-style abstracts label DOIs as "doi: " followed by the registrant/suffix,
# often wrapped onto the next line. The suffix can contain dots, parens, angle
# brackets, semicolons (e.g. legacy SICI DOIs like
# 10.1002/(SICI)1098-1004(1999)13:4<301::AID-HUMU7>3.0.CO;2-V), so we grab
# everything up to the next whitespace and clean trailing punctuation.
DOI_KEYWORD_RE = re.compile(
    r"doi:\s*\n?\s*(10\.\d{4,9}/\S+)",
    re.IGNORECASE,
)
# Fallback for DOIs that appear without a "doi:" label
DOI_BARE_RE = re.compile(
    r"\b(10\.\d{4,9}/\S+)",
    re.IGNORECASE,
)
# Strip trailing punctuation noise like "." ";" "," ")" "]"
DOI_TRAILING_JUNK = re.compile(r"[.;,\]\)\>]+$")

UA = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)
# Hosts known to require interactive browser fetch (Cloudflare / strict bot
# detection on OA content). When we hit these and Unpaywall says is_oa=True,
# we surface them as a distinct category so they can be sent to the browser
# fetch workflow rather than being lumped in with "paywalled".
BROWSER_BLOCKED_HOSTS = {
    "karger.com",
    "www.karger.com",
    "ahajournals.org",
    "www.ahajournals.org",
    # Linkinghub is Elsevier's JS-redirect shell. The HTTP response is a 2KB
    # stub that requires a browser to follow into ScienceDirect.
    "linkinghub.elsevier.com",
    # OUP and JBC return 403 to bot UAs even when Unpaywall says is_oa
    "academic.oup.com",
    "www.jbc.org",
    # Wiley OnlineLibrary 403s bot requests; the Wiley TDM API key in .env
    # would be the proper route but is not used by this script.
    "onlinelibrary.wiley.com",
}

# Patterns that suggest a link is the article PDF on an institutional repo.
# We follow these from a landing page if direct downloads have failed.
IR_PDF_LINK_PATTERNS = (
    re.compile(r"\.pdf(?:[?#]|$)", re.IGNORECASE),
    re.compile(r"/(?:download|fulltext|content|pdf|attachment)/", re.IGNORECASE),
    re.compile(r"viewcontent\.cgi", re.IGNORECASE),  # bepress / Digital Commons
)
# Anchor text that frequently labels the article PDF download
IR_PDF_LINK_TEXT_RE = re.compile(
    r"(?:download|full[\s_-]*text|view\s*pdf|article\s*pdf|preprint|accepted\s*manuscript|aam)",
    re.IGNORECASE,
)

# Hosts that index OA copies but don't actually serve the article body.
# Visiting them returns a catalog/listing page, not paper text.
CATALOG_HOSTS = {
    "doaj.org",
    "www.doaj.org",
    "openalex.org",
    "api.openalex.org",
    "scholar.google.com",
    "core.ac.uk",  # often only metadata for closed papers
    "semanticscholar.org",
    "www.semanticscholar.org",
    "researchgate.net",  # paywalled / login-walled
    "www.researchgate.net",
}

# Section-heading hints that a converted document is real article content,
# not a nav/landing page. Case-insensitive substring match on the rendered
# markdown.
ARTICLE_SECTION_HINTS = (
    "introduction",
    "methods",
    "materials and methods",
    "results",
    "discussion",
    "conclusion",
    "references",
    "patients and methods",
    "subjects and methods",
    "abstract",
)
MIN_ARTICLE_CHARS = 6000
MIN_ARTICLE_SECTION_HITS = 3

# Strong paywall / stub-page indicators. If any of these phrases appears in
# the converted text, the page is rejected even if it nominally has section
# headings — those headings come from the publisher's site chrome, not the
# article body. Tuned against false positives observed on Nature, Maastricht
# CRIS, NCBI bookshelf landing pages.
PAYWALL_STUB_PHRASES = (
    "log in or create a free account",
    "log in or create an account",
    "access through your institution",
    "buy this article",
    "rent or buy this article",
    "subscribe to journal",
    "purchase this article",
    "access denied",
    "full text availability",  # PMC "no full text" landing
    "this is an excerpt",
    "request full-text pdf",
    "request full text",
    "sign in to view",
)
# A real article should mention research-paper-style phrases. We require at
# least N matches to consider the body genuine (defensive against pages
# where the abstract alone gives 1-2 hits).
ARTICLE_DEEP_PHRASES = (
    " p < ",
    " p <",
    " p = ",
    " p =",
    " n = ",
    " n =",
    " ± ",
    " +/- ",
    "patients with",
    "we identified",
    "we report",
    "we found",
    "in this study",
    "table 1",
    "table 2",
    "figure 1",
    "figure 2",
    "et al.",
)
MIN_DEEP_PHRASE_HITS = 4

# ------------------------------------------------------------------ helpers


def find_abstract_only_pmids(harvest_dir: Path) -> List[str]:
    """Return sorted PMIDs that look abstract-only in the harvest dir."""
    pmids: list[str] = []
    for f in sorted(harvest_dir.glob("*_FULL_CONTEXT.md")):
        pmid = f.name.split("_")[0]
        try:
            head = f.read_text(encoding="utf-8", errors="ignore")[:8192]
        except OSError:
            continue
        size = f.stat().st_size
        if ABSTRACT_MARKER in head or size < ABSTRACT_FILE_SIZE_THRESHOLD:
            pmids.append(pmid)
    return pmids


def extract_doi_from_context(full_context: str) -> Optional[str]:
    """Pull the first DOI-looking token from the FULL_CONTEXT text.

    PubMed's abstract block usually contains `doi: 10.xxxx/yyy.` near the top,
    sometimes wrapped onto the next line. Prefer the labeled form when present.

    Strip publisher-image asset paths (e.g. ".../asset/.../graphic/...jpeg")
    that occasionally appear concatenated with the canonical DOI in
    PubMed-formatted abstracts (notably AHA/Circulation).
    """
    m = DOI_KEYWORD_RE.search(full_context)
    if not m:
        m = DOI_BARE_RE.search(full_context)
    if not m:
        return None
    doi = m.group(1).strip()
    # Truncate publisher-page extensions that got concatenated to the DOI
    # (AHA pages sometimes embed image-asset paths inside the DOI string).
    for marker in ("/asset/", "/graphic/", "/figure/", "/full/", "/pdf/", "/-/"):
        idx = doi.lower().find(marker)
        if idx > 0:
            doi = doi[:idx]
    # Strip trailing punctuation iteratively (e.g. ".)" or ";.")
    while True:
        cleaned = DOI_TRAILING_JUNK.sub("", doi)
        if cleaned == doi:
            break
        doi = cleaned
    return doi if doi else None


def fetch_doi_via_entrez(
    pmid: str, email: str, session: requests.Session
) -> Optional[str]:
    """Resolve a DOI for a PMID via NCBI ESummary as a fallback."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "pubmed",
        "id": pmid,
        "retmode": "json",
        "email": email,
        "tool": "GeneVariantFetcher",
    }
    try:
        r = session.get(url, params=params, timeout=20)
        if r.status_code != 200:
            return None
        data = r.json()
        result = data.get("result", {}).get(pmid, {})
        for entry in result.get("articleids", []):
            if entry.get("idtype") == "doi":
                return entry.get("value")
    except Exception as e:
        log.debug("Entrez DOI lookup failed for %s: %s", pmid, e)
    return None


def publisher_from_doi(doi: str) -> str:
    """Crude publisher tag from DOI prefix (registrant)."""
    prefix = doi.split("/", 1)[0]
    return prefix


# ------------------------------------------------------------------ Unpaywall


def query_unpaywall(
    client: UnpaywallClient, doi: str
) -> Tuple[Optional[dict], Optional[str]]:
    """Return the full Unpaywall payload (not the trimmed result)."""
    if not doi:
        return None, "no-doi"
    encoded = quote(doi.strip(), safe="/")
    url = f"https://api.unpaywall.org/v2/{encoded}"
    try:
        r = client.session.get(url, params={"email": client.email}, timeout=30)
        if r.status_code == 404:
            return None, "doi-not-found"
        if r.status_code == 422:
            return None, "invalid-doi"
        if r.status_code != 200:
            return None, f"http-{r.status_code}"
        if not r.text.strip():
            time.sleep(1.0)
            r = client.session.get(url, params={"email": client.email}, timeout=30)
            if r.status_code != 200 or not r.text.strip():
                return None, "empty-response"
        try:
            return r.json(), None
        except ValueError:
            return None, "invalid-json"
    except requests.RequestException as e:
        return None, f"request-error:{e}"


# ------------------------------------------------------------------ preprints


def query_biorxiv(doi: str, session: requests.Session) -> Optional[dict]:
    """Check bioRxiv/medRxiv for a preprint of this DOI."""
    if not doi:
        return None
    for server in ("biorxiv", "medrxiv"):
        try:
            r = session.get(
                f"https://api.biorxiv.org/details/{server}/{doi}",
                timeout=20,
                headers={"User-Agent": UA},
            )
            if r.status_code == 200:
                data = r.json()
                if isinstance(data, dict) and data.get("collection"):
                    return {"server": server, "records": data["collection"]}
        except Exception as e:
            log.debug("%s lookup failed for %s: %s", server, doi, e)
    return None


def biorxiv_pdf_url(record: dict) -> Optional[str]:
    """Construct a PDF URL for a bioRxiv/medRxiv record."""
    doi = record.get("doi")
    server = record.get("server")
    if not doi or not server:
        return None
    return (
        f"https://www.{server}.org/content/{doi}v{record.get('version', '1')}.full.pdf"
    )


# ------------------------------------------------------------------ download


def download_binary(
    url: str, session: requests.Session, timeout: int = 60
) -> Tuple[Optional[bytes], Optional[int]]:
    """Download a URL; return (bytes, status_code). bytes is None on failure."""
    try:
        r = session.get(
            url,
            timeout=timeout,
            headers={
                "User-Agent": UA,
                "Accept": "text/html,application/pdf,application/xhtml+xml,*/*",
                "Accept-Language": "en-US,en;q=0.9",
            },
            allow_redirects=True,
        )
        if r.status_code != 200:
            return None, r.status_code
        return r.content, 200
    except Exception as e:
        log.debug("download failed for %s: %s", url, e)
        return None, None


def _is_browser_blocked_host(url: str) -> bool:
    try:
        host = (urlparse(url).hostname or "").lower()
    except ValueError:
        return False
    return host in BROWSER_BLOCKED_HOSTS


def looks_like_pdf(data: bytes) -> bool:
    return bool(data) and data[:5] == b"%PDF-"


def convert_pdf_bytes(data: bytes) -> Optional[str]:
    """Convert PDF bytes to markdown text. Tries markitdown then pdftotext."""
    if not looks_like_pdf(data):
        return None

    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tf:
        tf.write(data)
        tf_path = tf.name
    try:
        # markitdown
        try:
            from markitdown import MarkItDown  # noqa: WPS433

            md = MarkItDown()
            result = md.convert(tf_path)
            if (
                result
                and result.text_content
                and len(result.text_content.strip()) > 200
            ):
                return result.text_content
        except Exception as e:
            log.debug("markitdown failed: %s", e)

        # pdftotext fallback
        import subprocess

        try:
            out = subprocess.run(
                ["pdftotext", "-layout", "-q", tf_path, "-"],
                capture_output=True,
                timeout=60,
            )
            if out.returncode == 0 and out.stdout and len(out.stdout.strip()) > 200:
                return out.stdout.decode("utf-8", errors="ignore")
        except FileNotFoundError:
            log.debug("pdftotext not installed")
        except Exception as e:
            log.debug("pdftotext failed: %s", e)
    finally:
        try:
            os.unlink(tf_path)
        except OSError:
            pass
    return None


def convert_html_bytes(data: bytes, url: str) -> Optional[str]:
    """Convert HTML bytes to markdown. Tries markitdown directly."""
    if not data:
        return None
    text_head = data[:512].lower()
    if (
        b"<html" not in text_head
        and b"<!doctype html" not in text_head
        and b"<body" not in text_head
    ):
        # Allow XML/JATS too — markitdown may still cope
        if not (b"<article" in text_head or b"<?xml" in text_head):
            return None

    try:
        from markitdown import MarkItDown  # noqa: WPS433

        # Hint markitdown about html
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as tf:
            tf.write(data)
            tf_path = tf.name
        try:
            md = MarkItDown()
            result = md.convert(tf_path)
            if (
                result
                and result.text_content
                and len(result.text_content.strip()) > 300
            ):
                return result.text_content
        finally:
            try:
                os.unlink(tf_path)
            except OSError:
                pass
    except Exception as e:
        log.debug("html→md via markitdown failed for %s: %s", url, e)

    # Fallback: strip tags crudely. Better than nothing for HTML landing pages
    try:
        from html.parser import HTMLParser  # noqa: WPS433

        class _Stripper(HTMLParser):
            def __init__(self) -> None:
                super().__init__()
                self.parts: list[str] = []
                self._skip = 0

            def handle_starttag(self, tag, attrs):
                if tag in ("script", "style", "noscript"):
                    self._skip += 1

            def handle_endtag(self, tag):
                if tag in ("script", "style", "noscript") and self._skip:
                    self._skip -= 1

            def handle_data(self, data):
                if not self._skip:
                    s = data.strip()
                    if s:
                        self.parts.append(s)

        stripper = _Stripper()
        stripper.feed(data.decode("utf-8", errors="ignore"))
        text = "\n".join(stripper.parts)
        if len(text) > 1000:
            return text
    except Exception as e:
        log.debug("html strip fallback failed for %s: %s", url, e)
    return None


# ------------------------------------------------------------------ pipeline


@dataclass
class RecoveryRecord:
    pmid: str
    doi: Optional[str] = None
    doi_source: Optional[str] = None  # context|entrez|none
    is_oa: Optional[bool] = None
    oa_status: Optional[str] = None
    publisher: Optional[str] = None
    journal: Optional[str] = None
    unpaywall_error: Optional[str] = None
    biorxiv_hit: bool = False
    oa_locations_tried: List[dict] = field(default_factory=list)
    recovered: bool = False
    recovered_via: Optional[str] = None  # url that worked
    recovered_kind: Optional[str] = None  # pdf|html
    recovered_chars: int = 0
    final_reason: Optional[str] = None
    browser_blocked_hosts: List[str] = field(default_factory=list)


def candidate_urls_from_unpaywall(payload: dict) -> List[Tuple[str, str]]:
    """Return list of (url, kind_hint) candidates from Unpaywall payload.

    kind_hint is 'pdf' or 'landing'. PDF candidates come first.
    """
    cands: list[tuple[str, str]] = []
    seen: set[str] = set()

    def add(url: Optional[str], kind: str) -> None:
        if not url:
            return
        if url in seen:
            return
        seen.add(url)
        cands.append((url, kind))

    best = payload.get("best_oa_location") or {}
    add(best.get("url_for_pdf"), "pdf")
    add(best.get("url_for_landing_page"), "landing")
    add(best.get("url"), "landing")
    for loc in payload.get("oa_locations") or []:
        add(loc.get("url_for_pdf"), "pdf")
        add(loc.get("url_for_landing_page"), "landing")
        add(loc.get("url"), "landing")
    # Sort: PDFs first, but preserve original order within each kind
    pdfs = [c for c in cands if c[1] == "pdf"]
    rest = [c for c in cands if c[1] != "pdf"]
    return pdfs + rest


def try_recover(
    pmid: str,
    doi: Optional[str],
    session: requests.Session,
    unpaywall: UnpaywallClient,
) -> Tuple[RecoveryRecord, Optional[str]]:
    """Try to recover full text for one paper. Returns (record, body_md_or_None)."""
    rec = RecoveryRecord(pmid=pmid, doi=doi)

    if doi:
        rec.publisher = publisher_from_doi(doi)
        payload, err = query_unpaywall(unpaywall, doi)
        if err:
            rec.unpaywall_error = err
        if payload:
            rec.is_oa = payload.get("is_oa")
            rec.oa_status = payload.get("oa_status")
            rec.journal = payload.get("journal_name")
            if payload.get("publisher"):
                rec.publisher = payload.get("publisher")

            candidates = candidate_urls_from_unpaywall(payload)
            for url, kind in candidates:
                body, status, blocked = _attempt_url(url, kind, session)
                tried = {
                    "url": url,
                    "kind": kind,
                    "ok": bool(body),
                    "status": status,
                    "browser_blocked": blocked,
                }
                rec.oa_locations_tried.append(tried)
                if blocked:
                    host = urlparse(url).hostname or url
                    if host not in rec.browser_blocked_hosts:
                        rec.browser_blocked_hosts.append(host)
                if body:
                    rec.recovered = True
                    rec.recovered_via = url
                    rec.recovered_kind = kind
                    rec.recovered_chars = len(body)
                    return rec, body

    # Preprint servers
    if doi:
        bx = query_biorxiv(doi, session)
        if bx:
            rec.biorxiv_hit = True
            for record in bx["records"][:3]:
                record["server"] = bx["server"]
                url = biorxiv_pdf_url(record)
                if not url:
                    continue
                body, status, _ = _attempt_url(url, "pdf", session)
                rec.oa_locations_tried.append(
                    {
                        "url": url,
                        "kind": "pdf-preprint",
                        "ok": bool(body),
                        "status": status,
                    }
                )
                if body:
                    rec.recovered = True
                    rec.recovered_via = url
                    rec.recovered_kind = "pdf-preprint"
                    rec.recovered_chars = len(body)
                    return rec, body

    if not rec.recovered:
        if not doi:
            rec.final_reason = "no-doi"
        elif rec.is_oa is False:
            rec.final_reason = "not-oa"
        elif rec.unpaywall_error:
            rec.final_reason = f"unpaywall:{rec.unpaywall_error}"
        elif rec.browser_blocked_hosts and rec.is_oa:
            rec.final_reason = "oa-but-browser-blocked"
        elif rec.oa_locations_tried:
            rec.final_reason = "oa-but-all-downloads-failed"
        else:
            rec.final_reason = "no-oa-locations"
    return rec, None


def _is_catalog_host(url: str) -> bool:
    try:
        host = urlparse(url).hostname or ""
    except ValueError:
        return False
    return host.lower() in CATALOG_HOSTS


def _looks_like_article(text: str) -> bool:
    """Heuristic: is this converted markdown a real paper body?

    The bar is "this looks like article body text", not "this is well-formed
    markdown". We need to reject Nature / Maastricht CRIS / NCBI landing
    pages that are essentially the abstract plus navigation.
    """
    if not text or len(text) < MIN_ARTICLE_CHARS:
        return False
    lower = text.lower()

    # Hard reject: paywall / stub phrases. A page that explicitly tells the
    # user to log in is not full text, regardless of length.
    for phrase in PAYWALL_STUB_PHRASES:
        if phrase in lower:
            return False

    # Require article section hints (Methods/Results/Discussion etc.)
    section_hits = sum(1 for h in ARTICLE_SECTION_HINTS if h in lower)
    if section_hits < MIN_ARTICLE_SECTION_HITS:
        return False

    # Require deep article phrases (research-paper language). Filters pages
    # where the section headings are present in nav links but the body text
    # is just the abstract.
    deep_hits = sum(1 for p in ARTICLE_DEEP_PHRASES if p in lower)
    if deep_hits < MIN_DEEP_PHRASE_HITS:
        return False

    # Reject pages dominated by site-chrome noise.
    noise = sum(lower.count(t) for t in ("cookie", "subscribe", "log in", "sign in"))
    if noise > 50 and len(text) < 12000:
        return False
    return True


def _find_ir_pdf_links(html: str, base_url: str) -> List[str]:
    """Find candidate PDF links inside an institutional-repo HTML landing page.

    Returns absolute URLs ordered by likelihood. Looks at <a href> with
    matching patterns and at the anchor text. De-duplicates while preserving
    order.
    """
    from html.parser import HTMLParser
    from urllib.parse import urljoin

    class _Collector(HTMLParser):
        def __init__(self) -> None:
            super().__init__()
            self.candidates: list[tuple[int, str]] = []
            self._open_a: Optional[str] = None
            self._a_text: list[str] = []

        def handle_starttag(self, tag, attrs):
            if tag == "a":
                ad = dict(attrs)
                href = ad.get("href")
                if href:
                    self._open_a = href
                    self._a_text = []

        def handle_data(self, data):
            if self._open_a is not None:
                self._a_text.append(data)

        def handle_endtag(self, tag):
            if tag == "a" and self._open_a is not None:
                href = self._open_a
                text = " ".join(self._a_text).strip()
                score = 0
                for pat in IR_PDF_LINK_PATTERNS:
                    if pat.search(href):
                        score += 2
                        break
                if IR_PDF_LINK_TEXT_RE.search(text):
                    score += 1
                if score:
                    self.candidates.append((score, href))
                self._open_a = None
                self._a_text = []

    parser = _Collector()
    try:
        parser.feed(html)
    except Exception:
        return []
    # Sort by score (desc) but stable
    parser.candidates.sort(key=lambda x: -x[0])
    seen: set[str] = set()
    out: list[str] = []
    for _, href in parser.candidates:
        absu = urljoin(base_url, href)
        # Skip obvious junk: privacy policies, terms, citations export
        low = absu.lower()
        if any(
            x in low
            for x in (
                "privacy",
                "/terms",
                "cookies",
                "citation",
                "/export",
                ".ris",
                ".bib",
                "twitter",
                "facebook",
            )
        ):
            continue
        if absu in seen:
            continue
        seen.add(absu)
        out.append(absu)
        if len(out) >= 5:
            break
    return out


def _attempt_url(
    url: str,
    kind_hint: str,
    session: requests.Session,
    follow_ir_links: bool = True,
) -> Tuple[Optional[str], Optional[int], bool]:
    """Download and convert one URL.

    Returns (markdown_body_or_None, http_status, was_blocked_bot_host).
    For HTML landing pages on institutional repositories, scrape for
    candidate PDF links and try those too (one level deep).
    """
    if _is_catalog_host(url):
        return None, None, False
    data, status = download_binary(url, session)
    blocked = status in (403, 401, 429) and _is_browser_blocked_host(url)
    if not data:
        return None, status, blocked

    if looks_like_pdf(data):
        body = convert_pdf_bytes(data)
        if body and _looks_like_article(body):
            return body, status, False
        return None, status, blocked

    # HTML path: try direct conversion first
    body = convert_html_bytes(data, url)
    if body and _looks_like_article(body):
        return body, status, False

    # Institutional-repo HTML landing pages: look for an actual PDF link
    if follow_ir_links:
        try:
            html = data.decode("utf-8", errors="ignore")
        except Exception:
            html = ""
        for pdf_url in _find_ir_pdf_links(html, url):
            sub, sub_status, sub_blocked = _attempt_url(
                pdf_url, "pdf", session, follow_ir_links=False
            )
            if sub:
                return sub, sub_status, False
    return None, status, blocked


# ------------------------------------------------------------------ output


def write_recovered_full_context(
    harvest_dir: Path,
    pmid: str,
    body_md: str,
    rec: RecoveryRecord,
) -> Path:
    """Write a new FULL_CONTEXT.md from the recovered body and preserve old."""
    target = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
    backup = harvest_dir / f"{pmid}_FULL_CONTEXT.ABSTRACT_ONLY.md"

    if target.exists() and not backup.exists():
        backup.write_bytes(target.read_bytes())

    header = (
        f"# RECOVERED VIA UNPAYWALL/PREPRINT\n\n"
        f"> Source URL: {rec.recovered_via}\n"
        f"> Kind: {rec.recovered_kind}\n"
        f"> Publisher: {rec.publisher}\n"
        f"> Journal: {rec.journal}\n"
        f"> OA status: {rec.oa_status}\n"
        f"> DOI: {rec.doi}\n"
        f"> Chars: {rec.recovered_chars}\n"
        f"> Original abstract-only stub preserved at {backup.name}\n\n---\n\n"
    )
    target.write_text(header + body_md, encoding="utf-8")
    return target


# ------------------------------------------------------------------ main


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--harvest-dir",
        required=True,
        type=Path,
        help="Path to pmc_fulltext harvest directory",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL", "brett.kroncke@gmail.com"),
    )
    parser.add_argument(
        "--out-report",
        type=Path,
        default=Path("recovery_report.json"),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Process only first N PMIDs (0 = all)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Do everything except writing recovered FULL_CONTEXT files",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=6,
        help="Concurrent recovery workers",
    )
    parser.add_argument(
        "--pmids",
        type=str,
        default="",
        help="Comma-separated subset of PMIDs to process (overrides discovery)",
    )
    args = parser.parse_args()

    harvest_dir: Path = args.harvest_dir.resolve()
    if not harvest_dir.is_dir():
        log.error("harvest dir not found: %s", harvest_dir)
        return 2

    if args.pmids:
        pmids = [p.strip() for p in args.pmids.split(",") if p.strip()]
    else:
        pmids = find_abstract_only_pmids(harvest_dir)
    if args.limit:
        pmids = pmids[: args.limit]
    log.info("processing %d PMIDs from %s", len(pmids), harvest_dir)

    session = requests.Session()
    session.headers.update({"User-Agent": UA})
    unpaywall = UnpaywallClient(email=args.email, session=session)

    # Step 1: resolve DOIs (sequential — cheap, hits same tools)
    doi_map: dict[str, tuple[Optional[str], str]] = {}
    for pmid in pmids:
        path = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
        doi: Optional[str] = None
        source = "none"
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            text = ""
        if text:
            doi = extract_doi_from_context(text)
            if doi:
                source = "context"
        if not doi:
            doi = fetch_doi_via_entrez(pmid, args.email, session)
            if doi:
                source = "entrez"
        doi_map[pmid] = (doi, source)

    with_doi = sum(1 for d, _ in doi_map.values() if d)
    log.info("resolved DOIs: %d/%d", with_doi, len(pmids))

    # Step 2: recover (parallel)
    records: list[RecoveryRecord] = []
    futures = {}
    body_by_pmid: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for pmid in pmids:
            doi, src = doi_map[pmid]
            futures[ex.submit(try_recover, pmid, doi, session, unpaywall)] = (pmid, src)
        for i, fut in enumerate(as_completed(futures), 1):
            pmid, src = futures[fut]
            try:
                rec, body = fut.result()
            except Exception as e:
                log.exception("recovery raised for %s", pmid)
                rec = RecoveryRecord(
                    pmid=pmid, doi=doi_map[pmid][0], final_reason=f"exception:{e}"
                )
                body = None
            rec.doi_source = src
            records.append(rec)
            if body:
                body_by_pmid[pmid] = body
            if i % 10 == 0 or i == len(pmids):
                log.info(
                    "progress %d/%d (recovered so far: %d)",
                    i,
                    len(pmids),
                    sum(1 for r in records if r.recovered),
                )

    # Step 3: write recovered bodies
    written = 0
    if not args.dry_run:
        for rec in records:
            if rec.recovered:
                body = body_by_pmid.get(rec.pmid)
                if body:
                    write_recovered_full_context(harvest_dir, rec.pmid, body, rec)
                    written += 1

    # Step 4: aggregate stats
    recovered = sum(1 for r in records if r.recovered)
    not_oa = sum(1 for r in records if r.is_oa is False)
    is_oa_failed = sum(1 for r in records if r.is_oa and not r.recovered)
    browser_blocked = sum(
        1 for r in records if r.browser_blocked_hosts and not r.recovered
    )
    no_doi = sum(1 for r in records if not r.doi)
    by_final_reason = Counter(r.final_reason for r in records if not r.recovered)
    unpaywall_errs = Counter(r.unpaywall_error for r in records if r.unpaywall_error)
    publisher_blocked: Dict[str, int] = defaultdict(int)
    for r in records:
        if not r.recovered:
            key = r.publisher or "unknown"
            publisher_blocked[key] += 1

    summary = {
        "harvest_dir": str(harvest_dir),
        "total_pmids": len(pmids),
        "doi_resolved": with_doi,
        "doi_missing": no_doi,
        "recovered": recovered,
        "oa_but_failed": is_oa_failed,
        "oa_but_browser_blocked": browser_blocked,
        "not_oa": not_oa,
        "files_written": written,
        "dry_run": args.dry_run,
        "unpaywall_errors": dict(unpaywall_errs),
        "publisher_blocked": dict(
            sorted(publisher_blocked.items(), key=lambda kv: -kv[1])
        ),
        "final_reason_counts": dict(by_final_reason),
    }

    args.out_report.write_text(
        json.dumps(
            {"summary": summary, "records": [asdict(r) for r in records]},
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    print("\n=== RECOVERY SUMMARY ===")
    for k, v in summary.items():
        if isinstance(v, dict):
            print(f"{k}:")
            for kk, vv in v.items():
                print(f"  {kk}: {vv}")
        else:
            print(f"{k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
