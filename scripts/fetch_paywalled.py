#!/usr/bin/env python3
"""Fetch paywalled biomedical articles using the local Chrome session.

Drives the Tier 3.5 browser stack with an ``AuthenticatedBrowserPool`` that
inherits cookies from the user's Chrome profile, so publishers that gate
content behind institutional SSO (e.g. VUMC's Microsoft sign-in) serve full
text instead of a paywall stub.

Usage::

    # Custom PMID list.
    python -m scripts.fetch_paywalled --pmid 15840476 --pmid 10973849 \\
        --output ./out --no-headless

    # Pipeline recovery CSV.
    python -m scripts.fetch_paywalled --input results/KCNH2/run/pmc_fulltext/paywalled_missing.csv \\
        --output results/KCNH2/run/pmc_fulltext

    # Pre-resolved DOI (skips the NCBI roundtrip).
    python -m scripts.fetch_paywalled --pmid-doi 15840476=10.1161/01.CIR.0000164255.06478.96 \\
        --output ./out

The script writes per-PMID artifacts under ``{output_dir}/{PMID}/`` (page
HTML, result.json, FULL_CONTEXT.md) *and* a canonical flat mirror at
``{output_dir}/{PMID}_FULL_CONTEXT.md`` so the main extraction discovery
path (``cli.extract.find_input_files`` with ``--full-text``) finds it
without per-PMID-dir glob plumbing. The two files always share content
(enriched unified markdown when available, body-only otherwise). A summary
line per PMID is printed at the end (publisher / chars / paywall-marker
count / outcome).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple
from urllib.parse import urlparse

# Path bootstrap so `python scripts/fetch_paywalled.py` works without -m.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import requests  # noqa: E402

from harvesting.browser_html import (  # noqa: E402
    AuthenticatedBrowserPool,
    BrowserHTMLFetcher,
    cookie_domain_summary,
    load_chrome_cookies,
    validate_article_content,
)
from harvesting.browser_html.dom_extract import extract_body_markdown  # noqa: E402
from harvesting.browser_html.strategies import (  # noqa: E402
    find_strategy,
    registered_names,
)
from harvesting.elsevier_api import ElsevierAPIClient  # noqa: E402
from harvesting.format_converters import FormatConverter  # noqa: E402
from harvesting.paywall_context_enrichment import (  # noqa: E402
    EnrichmentResult,
    enrich_paywall_full_context,
)
from harvesting.scholar_pdf_fallback import try_scholar_pdf  # noqa: E402
from harvesting.springer_api import SpringerAPIClient  # noqa: E402
from harvesting.supplement_scraper import SupplementScraper  # noqa: E402
from harvesting.wiley_api import WileyAPIClient  # noqa: E402


# Body selectors used by the PMC HTML fallback. PMC articles render the
# body inside <article>, <main>, or a #mc/.tsec container depending on the
# theme.
_PMC_BODY_SELECTORS = (
    "article",
    "main",
    "div.tsec",
    "#mc",
    "div#maincontent",
    "div.article-page",
    "body",
)

FETCH_SUCCESS_OUTCOMES = {
    "success",
    "success_via_pmc",
    "success_via_scholar_pdf",
    "success_via_elsevier_api",
    "success_via_publisher_api",
    "success_via_springer_api",
    "success_via_wiley_api",
    "success_supplement_only",
}

_BROWSER_SUPPLEMENT_SIZE_LIMIT_BYTES = 25 * 1024 * 1024
_HTMLISH_SUPPLEMENT_EXTS = {".html", ".htm", ".xhtml", ".xml"}


def _canonical_full_context_path(output_dir: Path, pmid: str) -> Path:
    """Path to the flat mirror that ``cli.extract`` discovers via glob."""
    return output_dir / f"{pmid}_FULL_CONTEXT.md"


def write_canonical_mirror(output_dir: Path, pmid: str, content: str) -> Path:
    """Write ``{output_dir}/{pmid}_FULL_CONTEXT.md`` and return its path.

    The flat layout matches ``cli.extract.find_input_files`` so the rescued
    paper is picked up by ``gvf extract --full-text`` without any per-PMID-
    directory awareness.
    """
    canonical = _canonical_full_context_path(output_dir, pmid)
    canonical.parent.mkdir(parents=True, exist_ok=True)
    canonical.write_text(content, encoding="utf-8")
    return canonical


def europepmc_lookup_pmcid(pmid: str, session: requests.Session) -> Optional[str]:
    """Return PMCID for a PMID if Europe PMC has one, else None.

    Europe PMC's search API returns the PMCID directly when an article has
    been deposited in PMC. NIH-funded papers usually do; subscription papers
    usually don't.
    """
    url = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        f"?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json&resultType=lite"
    )
    try:
        r = session.get(url, timeout=20)
        if r.status_code != 200:
            return None
        results = r.json().get("resultList", {}).get("result", [])
        if not results:
            return None
        pmcid = results[0].get("pmcid")
        if pmcid and not pmcid.startswith("PMC"):
            pmcid = f"PMC{pmcid}"
        return pmcid
    except Exception:
        return None


def try_pmc_fallback(
    pmid: str,
    output_dir: Path,
    session: requests.Session,
    converter: Optional[FormatConverter] = None,
    scraper: Optional[SupplementScraper] = None,
) -> Optional[Dict]:
    """Attempt to recover a stub via the PMC HTML deposit.

    Many subscription papers (especially NIH-funded ones) have a PMC version
    that's free to read after the embargo. We look up the PMCID via Europe
    PMC, fetch the PMC HTML, run our DOM extractor on it, and pull supplement
    links from the same HTML so captions AND downloadable supplements both
    land in the rescued FULL_CONTEXT.md.

    Returns a result row dict on success, or None if no PMC version exists
    or the extraction failed quality gates. On success the row's ``path``
    points at the rewritten FULL_CONTEXT.md.
    """
    pmcid = europepmc_lookup_pmcid(pmid, session)
    if not pmcid:
        return None

    pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
    try:
        r = session.get(
            pmc_url,
            timeout=30,
            headers={
                "User-Agent": (
                    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"
                ),
                "Accept": "text/html,application/xhtml+xml",
            },
        )
    except Exception as e:
        LOG.warning("PMC fallback fetch failed for PMID %s: %s", pmid, e)
        return None
    if not r.ok or len(r.content) < 5000:
        return None

    html = r.text
    markdown = extract_body_markdown(html, _PMC_BODY_SELECTORS)
    if not markdown:
        return None

    ok, reason = validate_article_content(markdown)
    if not ok:
        # PMC content didn't pass the gate either — could be an abstract-
        # only PMC record. Don't overwrite the stub from Tier 3.5.
        LOG.info("PMC fallback for %s extracted but failed gate: %s", pmid, reason)
        return None

    # Write a PMC FULL_CONTEXT.md and source marker.
    pmid_dir = output_dir / pmid
    pmid_dir.mkdir(parents=True, exist_ok=True)
    (pmid_dir / "page.html").write_text(html, encoding="utf-8")
    (pmid_dir / "source.txt").write_text(
        f"Recovered via PMC fallback: {pmc_url}\n", encoding="utf-8"
    )

    # Scrape supplement links directly from the PMC HTML so the enricher can
    # download convertible PMC supplements (xlsx/pdf/docx variant tables) into
    # the rescued FULL_CONTEXT.md. The orchestrator's PMC path normally does
    # this; in fallback mode we replicate it here rather than ship a body-
    # only context that silently drops supplement variants.
    scraper = scraper or SupplementScraper()
    try:
        supp_files = scraper.scrape_generic_supplements(html, pmc_url) or []
    except Exception as exc:
        LOG.info("PMC supplement scrape failed for %s: %s", pmid, exc)
        supp_files = []

    enrichment = enrich_paywall_full_context(
        body_markdown=markdown,
        html=html,
        supp_files=supp_files,
        pmid=pmid,
        output_dir=output_dir,
        converter=converter or FormatConverter(),
        session=session,
        download_supplements=True,
        source_url=pmc_url,
    )
    unified = enrichment.unified_markdown or markdown
    per_pmid_full_ctx = pmid_dir / "FULL_CONTEXT.md"
    per_pmid_full_ctx.write_text(unified, encoding="utf-8")
    canonical_path = write_canonical_mirror(output_dir, pmid, unified)

    return {
        "pmid": pmid,
        "strategy": "pmc_fallback",
        "outcome": "success",
        "reason": reason,
        "chars": len(unified),
        "supp_files": len(supp_files),
        "final_url": pmc_url,
        "path": str(canonical_path),
        "canonical_path": str(canonical_path),
        "per_pmid_path": str(per_pmid_full_ctx),
        "pmcid": pmcid,
        "figure_captions": enrichment.figure_caption_count,
        "table_captions": enrichment.table_caption_count,
        "supplements_downloaded": enrichment.supplement_count,
    }


def try_scholar_pdf_fallback(
    pmid: str,
    output_dir: Path,
    session: requests.Session,
    converter: Optional[FormatConverter] = None,
) -> Optional[Dict]:
    """Try Google Scholar title search for an author/lab-hosted PDF."""
    result = try_scholar_pdf(
        title=None,
        pmid=pmid,
        session=session,
        converter=converter or FormatConverter(),
        quality_gate=validate_article_content,
        email=_ncbi_email(),
    )
    if not result.success or not result.markdown:
        return None

    pmid_dir = output_dir / pmid
    pmid_dir.mkdir(parents=True, exist_ok=True)
    per_pmid_full_ctx = pmid_dir / "FULL_CONTEXT.md"
    per_pmid_full_ctx.write_text(result.markdown, encoding="utf-8")
    canonical_path = write_canonical_mirror(output_dir, pmid, result.markdown)
    (pmid_dir / "source.txt").write_text(
        f"Recovered via Google Scholar PDF: {result.source_url}\n",
        encoding="utf-8",
    )

    return {
        "reason": f"Google Scholar PDF ({result.title or 'title lookup'}): {result.detail}",
        "chars": len(result.markdown),
        "final_url": result.source_url or "",
        "path": str(canonical_path),
        "canonical_path": str(canonical_path),
        "per_pmid_path": str(per_pmid_full_ctx),
    }


def mark_scholar_pdf_success(row: Dict, scholar_row: Dict) -> None:
    """Merge a successful Scholar PDF fallback into a summary row."""
    base_strategy = row.get("strategy") or ""
    row["outcome"] = "success_via_scholar_pdf"
    row["strategy"] = (
        f"{base_strategy}+google_scholar_pdf"
        if base_strategy and base_strategy != "(none)"
        else "google_scholar_pdf"
    )
    row["reason"] = scholar_row["reason"]
    row["chars"] = scholar_row["chars"]
    row["final_url"] = scholar_row["final_url"]
    row["path"] = scholar_row["path"]
    row["canonical_path"] = scholar_row["canonical_path"]
    row["per_pmid_path"] = scholar_row["per_pmid_path"]
    row["supp_files"] = 0


def _write_publisher_api_result(
    *,
    pmid: str,
    doi: str,
    output_dir: Path,
    final_url: str,
    markdown: str,
    reason: str,
    publisher: str,
    source_label: str,
    outcome: str,
    supplement_markdown: str,
    supplements_downloaded: int,
) -> Dict:
    recovered_markdown = _append_recovered_supplement_markdown(
        markdown,
        supplement_markdown,
    )

    pmid_dir = output_dir / pmid
    pmid_dir.mkdir(parents=True, exist_ok=True)
    per_pmid_full_ctx = pmid_dir / "FULL_CONTEXT.md"
    per_pmid_full_ctx.write_text(recovered_markdown, encoding="utf-8")
    canonical_path = write_canonical_mirror(output_dir, pmid, recovered_markdown)
    (pmid_dir / "source.txt").write_text(
        f"Recovered via {source_label}: {doi}\n",
        encoding="utf-8",
    )
    meta = {
        "pmid": pmid,
        "publisher": publisher,
        "doi": doi,
        "final_url": final_url,
        "markdown_chars": len(recovered_markdown),
        "api_markdown_chars": len(markdown),
        "supplements_downloaded": supplements_downloaded,
        "canonical_full_context_path": str(canonical_path),
        "per_pmid_full_context_path": str(per_pmid_full_ctx),
        "notes": [],
        "error": None,
    }
    (pmid_dir / "result.json").write_text(
        json.dumps(meta, indent=2, default=str),
        encoding="utf-8",
    )

    return {
        "pmid": pmid,
        "strategy": publisher,
        "outcome": outcome,
        "api_label": source_label,
        "reason": reason,
        "chars": len(recovered_markdown),
        "supp_files": 0,
        "supplements_downloaded": supplements_downloaded,
        "final_url": final_url,
        "path": str(canonical_path),
        "canonical_path": str(canonical_path),
        "per_pmid_path": str(per_pmid_full_ctx),
    }


def try_elsevier_api_fallback(
    pmid: str,
    doi: Optional[str],
    output_dir: Path,
    *,
    final_url: str = "",
    session: Optional[requests.Session] = None,
    supplement_markdown: str = "",
    supplements_downloaded: int = 0,
) -> Optional[Dict]:
    """Recover Elsevier full text via the official API + insttoken.

    Browser recovery for Elsevier imprint sites can land on paywall stubs even
    when the Article Retrieval API returns the full XML, including table
    bodies. This path writes the same flat/per-PMID FULL_CONTEXT artifacts as
    the browser and PMC fallbacks.
    """
    if not doi or not ElsevierAPIClient.is_elsevier_doi(doi):
        return None

    client = ElsevierAPIClient(
        api_key=os.environ.get("ELSEVIER_API_KEY"),
        insttoken=os.environ.get("ELSEVIER_INSTTOKEN"),
        session=session if isinstance(session, requests.Session) else None,
    )
    if not client.is_available:
        return None

    markdown, error = client.fetch_fulltext(doi=doi, url=final_url or None)
    if not markdown:
        LOG.info("Elsevier API fallback failed for PMID %s: %s", pmid, error)
        return None

    ok, reason = validate_article_content(markdown)
    if not ok:
        LOG.info("Elsevier API fallback for PMID %s failed gate: %s", pmid, reason)
        return None

    return _write_publisher_api_result(
        pmid=pmid,
        doi=doi,
        output_dir=output_dir,
        final_url=final_url,
        markdown=markdown,
        reason=reason,
        publisher="elsevier_api",
        source_label="Elsevier Article Retrieval API",
        outcome="success_via_elsevier_api",
        supplement_markdown=supplement_markdown,
        supplements_downloaded=supplements_downloaded,
    )


def try_wiley_api_fallback(
    pmid: str,
    doi: Optional[str],
    output_dir: Path,
    *,
    final_url: str = "",
    session: Optional[requests.Session] = None,
) -> Optional[Dict]:
    """Recover Wiley full text via the Wiley TDM API."""
    if not doi or not WileyAPIClient.is_wiley_doi(doi):
        return None

    client = WileyAPIClient(
        api_key=os.environ.get("WILEY_API_KEY"),
        session=session if isinstance(session, requests.Session) else None,
    )
    if not client.is_available:
        return None

    markdown, error = client.fetch_fulltext(
        doi=doi,
        url=final_url or None,
        try_web_scraping=True,
    )
    if not markdown:
        LOG.info("Wiley API fallback failed for PMID %s: %s", pmid, error)
        return None

    ok, reason = validate_article_content(markdown)
    if not ok:
        LOG.info("Wiley API fallback for PMID %s failed gate: %s", pmid, reason)
        return None

    return _write_publisher_api_result(
        pmid=pmid,
        doi=doi,
        output_dir=output_dir,
        final_url=final_url,
        markdown=markdown,
        reason=reason,
        publisher="wiley_api",
        source_label="Wiley TDM API",
        outcome="success_via_wiley_api",
        supplement_markdown="",
        supplements_downloaded=0,
    )


def try_springer_api_fallback(
    pmid: str,
    doi: Optional[str],
    output_dir: Path,
    *,
    final_url: str = "",
    session: Optional[requests.Session] = None,
) -> Optional[Dict]:
    """Recover Springer/Nature/BMC full text via Springer OpenAccess API."""
    if not doi or not SpringerAPIClient.is_springer_doi(doi):
        return None

    client = SpringerAPIClient(
        api_key=os.environ.get("SPRINGER_API_KEY"),
        session=session if isinstance(session, requests.Session) else None,
    )
    if not client.is_available:
        return None

    markdown, _metadata, error = client.fetch_article(doi)
    if not markdown:
        LOG.info("Springer API fallback failed for PMID %s: %s", pmid, error)
        return None

    ok, reason = validate_article_content(markdown)
    if not ok:
        LOG.info("Springer API fallback for PMID %s failed gate: %s", pmid, reason)
        return None

    return _write_publisher_api_result(
        pmid=pmid,
        doi=doi,
        output_dir=output_dir,
        final_url=final_url,
        markdown=markdown,
        reason=reason,
        publisher="springer_api",
        source_label="Springer OpenAccess API",
        outcome="success_via_springer_api",
        supplement_markdown="",
        supplements_downloaded=0,
    )


def try_publisher_api_fallback(
    pmid: str,
    doi: Optional[str],
    output_dir: Path,
    *,
    final_url: str = "",
    session: Optional[requests.Session] = None,
    supplement_markdown: str = "",
    supplements_downloaded: int = 0,
) -> Optional[Dict]:
    """Try publisher API fallbacks in DOI-prefix order."""
    api_row = try_elsevier_api_fallback(
        pmid,
        doi,
        output_dir,
        final_url=final_url,
        session=session,
        supplement_markdown=supplement_markdown,
        supplements_downloaded=supplements_downloaded,
    )
    if api_row is not None:
        return api_row
    api_row = try_wiley_api_fallback(
        pmid,
        doi,
        output_dir,
        final_url=final_url,
        session=session,
    )
    if api_row is not None:
        return api_row
    return try_springer_api_fallback(
        pmid,
        doi,
        output_dir,
        final_url=final_url,
        session=session,
    )


LOG = logging.getLogger("fetch_paywalled")


def _ncbi_email() -> str:
    email = os.environ.get("ENTREZ_EMAIL") or os.environ.get("NCBI_EMAIL")
    if not email:
        raise RuntimeError("Set ENTREZ_EMAIL or NCBI_EMAIL before querying NCBI.")
    return email


def pubmed_resolve_doi(
    pmid: str, session: requests.Session, max_attempts: int = 4
) -> Optional[str]:
    """Resolve a PMID to its DOI via NCBI esummary.

    NCBI E-utilities rate-limits unauthenticated callers to 3 req/s; we retry
    with exponential backoff on 429. Returns None when no DOI is registered
    or the API repeatedly fails.
    """
    import time

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=pubmed&id={pmid}&retmode=json"
    )
    delay = 0.5
    for attempt in range(1, max_attempts + 1):
        try:
            r = session.get(url, timeout=20)
            if r.status_code == 429:
                time.sleep(delay)
                delay *= 2
                continue
            r.raise_for_status()
            item = r.json().get("result", {}).get(str(pmid), {})
            for artid in item.get("articleids", []) or []:
                if (artid.get("idtype") or "").lower() == "doi":
                    return artid.get("value")
            return None
        except Exception as e:
            if attempt == max_attempts:
                LOG.warning("DOI resolution failed for PMID %s: %s", pmid, e)
                return None
            time.sleep(delay)
            delay *= 2
    return None


def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update(
        {
            "User-Agent": (f"GVF-PaywalledFetch/1.0 (mailto:{_ncbi_email()})"),
        }
    )
    # Auto-route Cloudflare-blocked publisher URLs (Wiley/Karger/Sage) through the
    # institutional EZproxy when configured (GVF_EZPROXY_PREFIX/HOST). No-op
    # otherwise. Article HTML and /action/downloadSupplement files both route.
    from harvesting.browser_html import ezproxy

    ezproxy.install_on_session(s)
    if ezproxy.is_configured():
        LOG.info("EZproxy routing active for CF-blocked publishers")
    return s


def prime_authenticated_browser(
    pool: AuthenticatedBrowserPool,
    auth_urls: List[str],
    *,
    timeout_s: int,
) -> None:
    """Open authorization URLs and pause for a human-controlled browser session.

    This is for authorized publisher/institutional access, not unattended
    challenge solving. The caller can use a dedicated persistent profile so SSO
    cookies and publisher access state survive into later recovery batches.
    """
    urls = [u.strip() for u in auth_urls or [] if u and u.strip()]
    if not urls:
        return
    if not sys.stdin.isatty():
        raise RuntimeError("--auth-url requires an interactive terminal")
    if pool.headless:
        print("WARNING: --auth-url is intended for use with --no-headless.")

    timeout_ms = max(30, int(timeout_s or 180)) * 1000
    for url in urls:
        print(f"Opening auth URL: {url}")
        with pool.page() as page:
            try:
                page.set_default_timeout(timeout_ms)
            except Exception:
                pass
            try:
                page.goto(url, wait_until="domcontentloaded", timeout=timeout_ms)
            except Exception as exc:
                LOG.info("Auth URL load did not complete for %s: %s", url, exc)
            print(
                "Complete institutional login/access checks in the browser, "
                "then press Enter to continue."
            )
            input()


def hydrate_session_with_browser_cookies(
    session: requests.Session, cookies: List[dict]
) -> int:
    """Copy Playwright-format browser cookies into a requests session.

    The Playwright browser pool uses these cookies for publisher page access,
    but supplement downloads run through ``requests`` in the enrichment step.
    Mirroring the same cookie jar lets authenticated supplement files use the
    same institutional session as the rendered article page.
    """
    added = 0
    for cookie in cookies or []:
        name = cookie.get("name")
        value = cookie.get("value")
        if not name or value is None:
            continue

        domain = cookie.get("domain") or ""
        if not domain and cookie.get("url"):
            domain = urlparse(str(cookie["url"])).hostname or ""
        if not domain:
            continue

        session.cookies.set(
            str(name),
            str(value),
            domain=str(domain),
            path=str(cookie.get("path") or "/"),
        )
        added += 1
    return added


def _payload_looks_like_html_error(
    payload: bytes, *, content_type: str, filename: str
) -> bool:
    suffix = Path(filename).suffix.lower()
    if suffix in _HTMLISH_SUPPLEMENT_EXTS:
        return False
    content_type_l = content_type.lower()
    head = payload[:2048].lstrip().lower()
    if "text/html" in content_type_l:
        return True
    return (
        head.startswith(b"<!doctype html")
        or head.startswith(b"<html")
        or b"cf-mitigated" in head
        or b"cloudflare" in head
        or b"just a moment" in head
    )


def _browser_response_download(
    page: Any,
    *,
    url: str,
    file_path: Path,
    filename: str,
    source_url: str,
    timeout_ms: int,
) -> bool:
    headers = {"Referer": source_url} if source_url else None
    request = getattr(getattr(page, "context", None), "request", None)
    if request is None:
        return False

    kwargs: Dict[str, Any] = {"timeout": timeout_ms}
    if headers:
        kwargs["headers"] = headers
    try:
        response = request.get(url, max_redirects=10, **kwargs)
    except TypeError:
        response = request.get(url, **kwargs)
    except Exception as exc:
        LOG.info("Browser supplement request failed for %s: %s", url, exc)
        return False

    if not bool(getattr(response, "ok", False)):
        LOG.info(
            "Browser supplement request non-200 for %s: %s",
            url,
            getattr(response, "status", "(unknown)"),
        )
        return False

    raw_headers = getattr(response, "headers", {}) or {}
    headers_l = {str(k).lower(): str(v) for k, v in raw_headers.items()}
    content_length = headers_l.get("content-length")
    if content_length:
        try:
            if int(content_length) > _BROWSER_SUPPLEMENT_SIZE_LIMIT_BYTES:
                LOG.info("Browser supplement too large by header for %s", url)
                return False
        except ValueError:
            pass

    try:
        payload = response.body()
    except Exception as exc:
        LOG.info("Browser supplement response body failed for %s: %s", url, exc)
        return False

    if not payload or len(payload) > _BROWSER_SUPPLEMENT_SIZE_LIMIT_BYTES:
        return False
    if _payload_looks_like_html_error(
        payload,
        content_type=headers_l.get("content-type", ""),
        filename=filename,
    ):
        LOG.info("Browser supplement response looked like HTML challenge for %s", url)
        return False

    file_path.parent.mkdir(parents=True, exist_ok=True)
    file_path.write_bytes(payload)
    return file_path.exists() and file_path.stat().st_size > 0


def _save_playwright_download(download: Any, file_path: Path) -> bool:
    file_path.parent.mkdir(parents=True, exist_ok=True)
    download.save_as(str(file_path))
    if not file_path.exists() or file_path.stat().st_size == 0:
        return False
    if file_path.stat().st_size > _BROWSER_SUPPLEMENT_SIZE_LIMIT_BYTES:
        file_path.unlink(missing_ok=True)
        return False
    return True


def _browser_navigation_download(
    page: Any,
    *,
    url: str,
    file_path: Path,
    timeout_ms: int,
) -> bool:
    try:
        with page.expect_download(timeout=timeout_ms) as download_info:
            page.goto(url, wait_until="commit", timeout=timeout_ms)
        return _save_playwright_download(download_info.value, file_path)
    except Exception as exc:
        LOG.info("Browser supplement navigation download failed for %s: %s", url, exc)
        return False


def _browser_click_download(
    page: Any,
    *,
    url: str,
    source_url: str,
    file_path: Path,
    timeout_ms: int,
) -> bool:
    if not source_url:
        return False
    try:
        page.goto(source_url, wait_until="domcontentloaded", timeout=timeout_ms)
    except Exception as exc:
        LOG.info("Browser source navigation failed for %s: %s", source_url, exc)
        return False

    find_script = """
    (target) => {
      return Array.from(document.querySelectorAll('a[href]')).some((a) => {
        try {
          return new URL(a.getAttribute('href'), document.baseURI).href === target;
        } catch (e) {
          return false;
        }
      });
    }
    """
    try:
        if not page.evaluate(find_script, url):
            return False
    except Exception as exc:
        LOG.info("Browser supplement link lookup failed for %s: %s", url, exc)
        return False

    click_script = """
    (target) => {
      const links = Array.from(document.querySelectorAll('a[href]'));
      const match = links.find((a) => {
        try {
          return new URL(a.getAttribute('href'), document.baseURI).href === target;
        } catch (e) {
          return false;
        }
      });
      if (!match) {
        return false;
      }
      match.target = '_self';
      match.click();
      return true;
    }
    """
    try:
        with page.expect_download(timeout=timeout_ms) as download_info:
            if not page.evaluate(click_script, url):
                return False
        return _save_playwright_download(download_info.value, file_path)
    except Exception as exc:
        LOG.info("Browser supplement click download failed for %s: %s", url, exc)
        return False


def make_browser_supplement_download_fallback(
    fetcher: BrowserHTMLFetcher,
) -> Optional[Callable[[str, Path, str, str, Dict[str, Any]], bool]]:
    """Build a supplement downloader backed by the authenticated browser pool."""
    pool = getattr(fetcher, "_pool", None)
    if pool is None or not hasattr(pool, "page"):
        return None

    timeout_s = max(5, min(int(getattr(fetcher, "_timeout_s", 90) or 90), 30))
    timeout_ms = timeout_s * 1000

    def _fallback(
        url: str,
        file_path: Path,
        pmid: str,
        filename: str,
        supp: Dict[str, Any],
    ) -> bool:
        source_url = str(supp.get("source_url") or "")
        try:
            with pool.page() as page:
                try:
                    page.set_default_timeout(timeout_ms)
                except Exception:
                    pass
                if source_url:
                    try:
                        page.goto(
                            source_url,
                            wait_until="domcontentloaded",
                            timeout=timeout_ms,
                        )
                    except Exception as exc:
                        LOG.info(
                            "PMID %s browser supplement source load failed: %s",
                            pmid,
                            exc,
                        )
                if _browser_response_download(
                    page,
                    url=url,
                    file_path=file_path,
                    filename=filename,
                    source_url=source_url,
                    timeout_ms=timeout_ms,
                ):
                    return True
                if _browser_navigation_download(
                    page,
                    url=url,
                    file_path=file_path,
                    timeout_ms=timeout_ms,
                ):
                    return True
                return _browser_click_download(
                    page,
                    url=url,
                    source_url=source_url,
                    file_path=file_path,
                    timeout_ms=timeout_ms,
                )
        except Exception as exc:
            LOG.info(
                "PMID %s browser supplement fallback failed for %s: %s", pmid, url, exc
            )
            try:
                file_path.unlink(missing_ok=True)
            except OSError:
                pass
            return False

    return _fallback


def _append_recovered_supplement_markdown(
    markdown: str,
    supplement_markdown: str,
) -> str:
    supplement_markdown = (supplement_markdown or "").strip()
    if not supplement_markdown:
        return markdown
    return (
        markdown.rstrip()
        + "\n\n# RECOVERED PUBLISHER SUPPLEMENTS\n\n"
        + supplement_markdown
        + "\n"
    )


def parse_pmid_doi_overrides(raw: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for entry in raw or []:
        if "=" not in entry:
            continue
        k, v = entry.split("=", 1)
        k, v = k.strip(), v.strip()
        if k and v:
            out[k] = v
    return out


def read_pmid_csv(path: Path) -> Tuple[List[str], Dict[str, str]]:
    """Read PMID and optional DOI columns from a recovery CSV."""
    pmids: List[str] = []
    doi_overrides: Dict[str, str] = {}
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return pmids, doi_overrides
        lower_to_field = {name.strip().lower(): name for name in reader.fieldnames}
        pmid_col = None
        for candidate in ("pmid", "pubmed_id", "pubmed", "id"):
            if candidate in lower_to_field:
                pmid_col = lower_to_field[candidate]
                break
        if pmid_col is None:
            raise ValueError(
                f"No PMID column found in {path}; columns={reader.fieldnames}"
            )
        doi_col = None
        for candidate in ("doi", "article_doi"):
            if candidate in lower_to_field:
                doi_col = lower_to_field[candidate]
                break
        for row in reader:
            pmid = (row.get(pmid_col) or "").strip()
            if not pmid:
                continue
            pmids.append(pmid)
            doi = (row.get(doi_col) or "").strip() if doi_col else ""
            if doi:
                doi_overrides[pmid] = doi
    return pmids, doi_overrides


def write_outputs(
    pmid: str,
    result,
    output_dir: Path,
    converter: Optional[FormatConverter] = None,
    session: Optional[requests.Session] = None,
    enrich: bool = True,
    supplement_download_fallback: Optional[
        Callable[[str, Path, str, str, Dict[str, Any]], bool]
    ] = None,
) -> Tuple[Path, Path, Optional[str], Optional[EnrichmentResult]]:
    """Write FULL_CONTEXT.md (per-PMID dir + flat canonical mirror) + raw HTML.

    Returns ``(canonical_path, per_pmid_path, body_for_gating, enrichment)``.
    The canonical path is the flat mirror ``{output_dir}/{pmid}_FULL_CONTEXT.md``
    that ``cli.extract.find_input_files`` discovers; the per-PMID path is the
    artifact-bundle copy under ``{output_dir}/{pmid}/FULL_CONTEXT.md``.

    The body returned for gating is the *un-enriched* article text — that's
    what the quality validator should judge. Captions and supplements are
    appended afterward, so the gate's verdict isn't polluted by figure
    legends or supplement padding.
    """
    pmid_dir = output_dir / pmid
    pmid_dir.mkdir(parents=True, exist_ok=True)

    body = result.main_markdown or ""
    body_for_enrichment = body
    if not body and result.supp_files:
        body_for_enrichment = (
            "# MAIN TEXT\n\n"
            "[No usable article body recovered; supplement-only recovery.]\n\n"
        )

    if result.main_html:
        (pmid_dir / "page.html").write_text(result.main_html, encoding="utf-8")

    enrichment: Optional[EnrichmentResult] = None
    if enrich and body_for_enrichment:
        try:
            enrichment = enrich_paywall_full_context(
                body_markdown=body_for_enrichment,
                html=result.main_html or "",
                supp_files=list(result.supp_files or []),
                pmid=pmid,
                output_dir=output_dir,
                converter=converter or FormatConverter(),
                session=session,
                source_url=result.final_url,
                supplement_download_fallback=supplement_download_fallback,
            )
        except Exception as exc:
            LOG.warning("Enrichment failed for PMID %s: %s", pmid, exc)
            enrichment = None

    if enrichment is not None and enrichment.unified_markdown:
        full_text = enrichment.unified_markdown
    else:
        full_text = body

    per_pmid_full_ctx = pmid_dir / "FULL_CONTEXT.md"
    per_pmid_full_ctx.write_text(full_text, encoding="utf-8")
    canonical_path = write_canonical_mirror(output_dir, pmid, full_text)

    meta = {
        "pmid": pmid,
        "publisher": result.publisher,
        "final_url": result.final_url,
        "supp_files": result.supp_files,
        "figure_paths": [str(p) for p in (result.figure_paths or [])],
        "notes": result.notes,
        "error": result.error,
        "markdown_chars": len(body),
        "canonical_full_context_path": str(canonical_path),
        "per_pmid_full_context_path": str(per_pmid_full_ctx),
    }
    if enrichment is not None:
        meta["unified_chars"] = len(enrichment.unified_markdown or "")
        meta["figure_captions"] = enrichment.figure_caption_count
        meta["table_captions"] = enrichment.table_caption_count
        meta["supplements_downloaded"] = enrichment.supplement_count
    (pmid_dir / "result.json").write_text(
        json.dumps(meta, indent=2, default=str), encoding="utf-8"
    )

    return canonical_path, per_pmid_full_ctx, (body or None), enrichment


def fetch_one(
    fetcher: BrowserHTMLFetcher,
    pmid: str,
    doi: Optional[str],
    output_dir: Path,
    pmc_session: Optional[requests.Session] = None,
) -> Dict:
    """Run Tier 3.5 for a single PMID; return a result row.

    On stub/empty outcomes from Tier 3.5, attempt the PMC fallback (Europe
    PMC PMCID lookup + PMC HTML fetch) so subscription papers with NIH PMC
    deposits are recovered automatically.
    """
    row: Dict = {
        "pmid": pmid,
        "doi": doi or "",
        "strategy": "",
        "outcome": "skipped",
        "reason": "",
        "chars": 0,
        "supp_files": 0,
        "supplements_downloaded": 0,
        "final_url": "",
    }
    recovery_session = pmc_session or getattr(fetcher, "session", None)
    publisher_supp_file_count = 0

    def run_scholar_fallback() -> Optional[Dict]:
        if recovery_session is None:
            return None
        scholar_row = try_scholar_pdf_fallback(
            pmid,
            output_dir,
            recovery_session,
            converter=fetcher.converter,
        )
        if scholar_row is not None:
            mark_scholar_pdf_success(row, scholar_row)
        return scholar_row

    def promote_pmc_fallback(pmc_row: Dict) -> None:
        print(
            f"  PMC fallback succeeded: {pmc_row['pmcid']} -> {pmc_row['chars']} chars"
        )
        row["outcome"] = "success_via_pmc"
        row["strategy"] = f"{row['strategy']}+pmc_fallback"
        row["reason"] = f"PMC fallback ({pmc_row['pmcid']}): {pmc_row['reason']}"
        row["chars"] = pmc_row["chars"]
        row["final_url"] = pmc_row["final_url"]
        row["path"] = pmc_row["path"]
        row["canonical_path"] = pmc_row.get("canonical_path", pmc_row["path"])
        row["per_pmid_path"] = pmc_row.get("per_pmid_path", pmc_row["path"])
        row["pmcid"] = pmc_row["pmcid"]
        row["supp_files"] = pmc_row.get("supp_files", 0)
        row["figure_captions"] = pmc_row.get("figure_captions", 0)
        row["table_captions"] = pmc_row.get("table_captions", 0)
        row["supplements_downloaded"] = pmc_row.get("supplements_downloaded", 0)

    def promote_publisher_api_fallback(api_row: Dict) -> None:
        label = api_row.get("api_label") or "Publisher API"
        strategy_name = api_row.get("strategy") or "publisher_api"
        outcome = api_row.get("outcome") or "success_via_publisher_api"
        print(f"  {label} fallback succeeded: {api_row['chars']} chars")
        base_strategy = row.get("strategy") or ""
        row["outcome"] = outcome
        row["strategy"] = (
            f"{base_strategy}+{strategy_name}"
            if base_strategy and base_strategy != "(none)"
            else strategy_name
        )
        row["reason"] = f"{label} fallback: {api_row['reason']}"
        row["chars"] = api_row["chars"]
        row["final_url"] = api_row.get("final_url") or row.get("final_url") or ""
        row["path"] = api_row["path"]
        row["canonical_path"] = api_row.get("canonical_path", api_row["path"])
        row["per_pmid_path"] = api_row.get("per_pmid_path", api_row["path"])
        row["supp_files"] = max(
            int(row.get("supp_files") or 0), int(api_row.get("supp_files") or 0)
        )
        row["supplements_downloaded"] = max(
            int(row.get("supplements_downloaded") or 0),
            int(api_row.get("supplements_downloaded") or 0),
        )

    if not doi:
        row["reason"] = "no DOI"
        run_scholar_fallback()
        return row

    # Tell the user which strategy will run, before we burn time on the fetch.
    strategy = find_strategy(doi=doi, allowlist=None)
    row["strategy"] = strategy.NAME if strategy else "(none)"
    if strategy is None:
        row["reason"] = "no matching strategy"
        api_row = try_publisher_api_fallback(
            pmid,
            doi,
            output_dir,
            session=recovery_session,
        )
        if api_row is not None:
            promote_publisher_api_fallback(api_row)
            return row
        run_scholar_fallback()
        return row

    result = fetcher.fetch(pmid=pmid, doi=doi, pub_date=None)
    if result is None:
        row["outcome"] = "skipped"
        row["reason"] = "fetcher returned None"
        api_row = try_publisher_api_fallback(
            pmid,
            doi,
            output_dir,
            session=recovery_session,
        )
        if api_row is not None:
            promote_publisher_api_fallback(api_row)
            return row
        run_scholar_fallback()
        return row

    row["final_url"] = result.final_url or ""
    row["supp_files"] = len(result.supp_files or [])
    publisher_supp_file_count = int(row["supp_files"] or 0)

    supplement_download_fallback = make_browser_supplement_download_fallback(fetcher)
    canonical_path, per_pmid_path, body, enrichment = write_outputs(
        pmid,
        result,
        output_dir,
        session=pmc_session,
        supplement_download_fallback=supplement_download_fallback,
    )
    # The gate validates only the body, but the chars we report should
    # reflect what actually landed in FULL_CONTEXT.md so downstream
    # observers see whether captions/supplements were appended.
    if enrichment is not None and enrichment.unified_markdown:
        row["chars"] = len(enrichment.unified_markdown)
        row["figure_captions"] = enrichment.figure_caption_count
        row["table_captions"] = enrichment.table_caption_count
        row["supplements_downloaded"] = enrichment.supplement_count
    else:
        row["chars"] = len(body or "")
    # ``path`` is the canonical flat mirror — that's the file
    # ``cli.extract.find_input_files`` discovers via glob. ``per_pmid_path``
    # is preserved for diagnostics so summary.json carries both locations.
    row["path"] = str(canonical_path)
    row["canonical_path"] = str(canonical_path)
    row["per_pmid_path"] = str(per_pmid_path)

    # Apply quality gate (the fetcher itself also runs the validator, but we
    # surface the reason here for the per-PMID summary table).
    ok, reason = validate_article_content(body or "")
    if result.error:
        row["outcome"] = "error"
        row["reason"] = result.error
    elif not body:
        row["outcome"] = "empty"
        row["reason"] = ";".join(result.notes) or "no markdown"
    elif not ok:
        row["outcome"] = "paywall_or_stub"
        row["reason"] = reason
    else:
        row["outcome"] = "success"
        row["reason"] = reason

    if (
        row["outcome"] in ("empty", "paywall_or_stub")
        and enrichment is not None
        and enrichment.supplement_count > 0
        and enrichment.unified_markdown
    ):
        row["outcome"] = "success_supplement_only"
        row["reason"] = (
            f"supplement-only recovery: downloaded "
            f"{enrichment.supplement_count} supplement(s)"
        )

    if row["outcome"] in (
        "paywall_or_stub",
        "empty",
        "error",
        "success_supplement_only",
    ):
        api_row = try_publisher_api_fallback(
            pmid,
            doi,
            output_dir,
            final_url=row.get("final_url") or "",
            session=recovery_session,
            supplement_markdown=(
                enrichment.supplement_markdown
                if enrichment is not None and enrichment.supplement_count > 0
                else ""
            ),
            supplements_downloaded=(
                enrichment.supplement_count if enrichment is not None else 0
            ),
        )
        if api_row is not None:
            promote_publisher_api_fallback(api_row)

    # PMC fallback. If Tier 3.5 left us with a stub/empty body, or recovered
    # a publisher page whose supplement downloads all failed, try the Europe
    # PMC route. Many NIH-funded subscription papers have a free PMC deposit,
    # and PMC supplement URLs can be easier to download than publisher copies.
    current_supplements = int(row.get("supplements_downloaded") or 0)
    publisher_supplements_failed = (
        row["outcome"]
        in {
            "success",
            "success_supplement_only",
            "success_via_elsevier_api",
            "success_via_springer_api",
            "success_via_wiley_api",
        }
        and (publisher_supp_file_count > 0 or int(row.get("supp_files") or 0) > 0)
        and current_supplements == 0
    )
    if (
        row["outcome"] in ("paywall_or_stub", "empty", "error")
        or publisher_supplements_failed
    ) and recovery_session is not None:
        pmc_row = try_pmc_fallback(
            pmid,
            output_dir,
            recovery_session,
            converter=fetcher.converter,
            scraper=fetcher.scraper,
        )
        if pmc_row is not None:
            pmc_supplements = int(pmc_row.get("supplements_downloaded") or 0)
            if (
                row["outcome"] in ("paywall_or_stub", "empty", "error")
                or pmc_supplements > current_supplements
                or int(pmc_row.get("chars") or 0) > int(row.get("chars") or 0) + 2048
            ):
                promote_pmc_fallback(pmc_row)

    if (
        row["outcome"] in ("paywall_or_stub", "empty", "error")
        and recovery_session is not None
    ):
        scholar_row = run_scholar_fallback()
        if scholar_row is not None:
            print(
                f"  Google Scholar PDF fallback succeeded: {scholar_row['chars']} chars"
            )

    return row


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--pmid",
        action="append",
        default=None,
        help="PMID to fetch (repeat for multiple). Defaults to the 4 top-value test PMIDs.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="CSV with a PMID column and optional DOI column, e.g. paywalled_missing.csv.",
    )
    parser.add_argument(
        "--pmid-doi",
        action="append",
        default=[],
        help="PMID=DOI override (repeat). Skips the NCBI esummary lookup.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("./paywall_test"),
        help="Output directory (default: ./paywall_test)",
    )
    parser.add_argument(
        "--output-dir",
        dest="output",
        type=Path,
        help="Alias for --output.",
    )
    parser.add_argument(
        "--profile",
        default=None,
        help="Chrome profile name (e.g. 'Default'). Default: merge all profiles.",
    )
    parser.add_argument(
        "--no-cookies",
        action="store_true",
        help="Skip Chrome cookie loading. Useful when macOS Keychain is locked; publisher SSO may not work.",
    )
    parser.add_argument(
        "--cookie-timeout-s",
        type=float,
        default=8.0,
        help="Per-domain Chrome cookie loading timeout in seconds (default: 8).",
    )
    parser.add_argument(
        "--headless",
        dest="headless",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "--no-headless",
        dest="headless",
        action="store_false",
        help="Run with a visible browser window (useful for debugging SSO challenges).",
    )
    parser.add_argument(
        "--allow",
        action="append",
        default=None,
        help="Strategy NAME to allow (repeat). Default: all registered strategies.",
    )
    parser.add_argument(
        "--no-chrome-channel",
        dest="use_chrome_channel",
        action="store_false",
        default=True,
        help="Use bundled Chromium instead of the locally-installed Chrome binary.",
    )
    parser.add_argument(
        "--browser-profile-dir",
        type=Path,
        default=None,
        help=(
            "Dedicated persistent browser user-data directory for authorized "
            "publisher/SSO sessions. Use with --no-headless and do not point "
            "at a daily Chrome profile."
        ),
    )
    parser.add_argument(
        "--stealth",
        action="store_true",
        help=(
            "Use the patchright stealth backend (real Chrome, forced headful) to "
            "clear Cloudflare managed challenges (e.g. Wiley Online Library) "
            "without a proxy. Requires `pip install patchright && patchright "
            "install chromium`. Best combined with --browser-profile-dir."
        ),
    )
    parser.add_argument(
        "--auth-url",
        action="append",
        default=[],
        help=(
            "URL to open before fetching, then pause for Enter after manual "
            "institutional login/access checks. Repeat for multiple publishers."
        ),
    )
    parser.add_argument(
        "--auth-timeout-s",
        type=int,
        default=180,
        help="Navigation timeout for each --auth-url page (default: 180).",
    )
    parser.add_argument(
        "--timeout-s",
        type=int,
        default=120,
        help="Per-paper timeout in seconds (default: 120).",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    input_pmids: List[str] = []
    input_overrides: Dict[str, str] = {}
    if args.input:
        input_pmids, input_overrides = read_pmid_csv(args.input.expanduser())

    pmids: List[str] = []
    pmids.extend(input_pmids)
    if args.pmid:
        pmids.extend(args.pmid)
    if not pmids:
        parser.error("Provide at least one PMID via --input or --pmid.")
    pmids = list(dict.fromkeys(str(p).strip() for p in pmids if str(p).strip()))

    overrides = parse_pmid_doi_overrides(args.pmid_doi)
    overrides = {**input_overrides, **overrides}
    output_dir: Path = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    browser_profile_dir = (
        args.browser_profile_dir.expanduser().resolve()
        if args.browser_profile_dir
        else None
    )

    print(f"PMIDs:       {', '.join(pmids)}")
    print(f"Output dir:  {output_dir}")
    if browser_profile_dir:
        print(f"Profile dir: {browser_profile_dir}")
    print(f"Strategies:  {registered_names()}")
    print(f"Headless:    {args.headless}")
    print()

    # ---- 1. Load cookies from local Chrome ----
    print("Loading Chrome cookies…")
    if args.no_cookies:
        cookies = []
        print("  skipped by --no-cookies")
    else:
        cookies = load_chrome_cookies(
            profile_name=args.profile,
            timeout_seconds=args.cookie_timeout_s,
        )
    summary = cookie_domain_summary(cookies)
    print(f"  total cookies: {len(cookies)}")
    for d, n in list(summary.items())[:20]:
        print(f"    {d:40s} {n}")
    print()

    # ---- 2. Build the authenticated pool ----
    pool = AuthenticatedBrowserPool(
        cookies=cookies,
        headless=args.headless,
        use_chrome_channel=args.use_chrome_channel,
        persistent_profile_path=str(browser_profile_dir)
        if browser_profile_dir
        else None,
        use_stealth=args.stealth,
    )

    prime_authenticated_browser(
        pool,
        args.auth_url,
        timeout_s=args.auth_timeout_s,
    )

    # ---- 3. Build minimal settings for the fetcher ----
    class _S:
        enable_browser_html_fallback = True
        browser_html_publisher_allowlist = args.allow or registered_names()
        browser_html_headless = args.headless
        browser_html_max_per_run = max(50, len(pmids) * 2)
        browser_html_per_paper_timeout_s = args.timeout_s
        browser_html_min_embargo_months = None

    scraper = SupplementScraper()
    converter = FormatConverter()
    session = make_session()
    session_cookie_count = hydrate_session_with_browser_cookies(session, cookies)
    print(f"Hydrated requests session with {session_cookie_count} browser cookies.")
    print()

    fetcher = BrowserHTMLFetcher(
        scraper=scraper,
        converter=converter,
        session=session,
        output_dir=output_dir,
        settings=_S(),
        validate_content_quality=validate_article_content,
        pool=pool,
        bypass_embargo=True,
        enabled_override=True,
    )

    # ---- 4. Resolve DOIs ----
    import time as _time

    print("Resolving DOIs…")
    doi_map: Dict[str, Optional[str]] = {}
    for pmid in pmids:
        if pmid in overrides:
            doi_map[pmid] = overrides[pmid]
            print(f"  {pmid}  (override) -> {doi_map[pmid]}")
            continue
        doi = pubmed_resolve_doi(pmid, session)
        doi_map[pmid] = doi
        print(f"  {pmid}  -> {doi or '(no DOI)'}")
        _time.sleep(0.4)  # NCBI E-utilities: 3 req/s cap without API key.
    print()

    # ---- 5. Fetch ----
    import time as _t

    def _domain_for(doi: Optional[str]) -> str:
        # Map DOI prefixes to the publisher domain so we can pace per-host.
        if not doi:
            return ""
        if doi.startswith("10.1161"):
            return "ahajournals.org"
        if doi.startswith("10.1016"):
            return "sciencedirect.com"
        if doi.startswith("10.4065"):
            return "mayoclinicproceedings.org"
        if doi.startswith("10.1093"):
            return "academic.oup.com"
        if doi.startswith("10.1002") or doi.startswith("10.1111"):
            return "onlinelibrary.wiley.com"
        if doi.startswith("10.1159"):
            return "karger.com"
        if doi.startswith("10.1212"):
            return "neurology.org"
        return "other"

    last_domain_at: Dict[str, float] = {}
    rows: List[Dict] = []

    def _run_with_pacing(pmid: str, gap_seconds: float) -> Dict:
        dom = _domain_for(doi_map[pmid])
        now = _t.time()
        wait_for = gap_seconds - (now - last_domain_at.get(dom, 0.0))
        if wait_for > 0 and last_domain_at:
            print(f"  pacing: sleeping {wait_for:.1f}s before {dom}")
            _t.sleep(wait_for)
        r = fetch_one(fetcher, pmid, doi_map[pmid], output_dir, pmc_session=session)
        last_domain_at[dom] = _t.time()
        return r

    try:
        for i, pmid in enumerate(pmids):
            print(f"--- PMID {pmid} ---")
            row = _run_with_pacing(pmid, gap_seconds=15.0)
            rows.append(row)
            print(
                f"  strategy={row['strategy']:14s} outcome={row['outcome']:18s} "
                f"chars={row['chars']:>7d} "
                f"supp_links={int(row.get('supp_files') or 0):>2d} "
                f"supp_downloaded={int(row.get('supplements_downloaded') or 0):>2d} "
                f"url={row['final_url'][:80]}"
            )
            if row.get("reason"):
                print(f"  reason: {row['reason']}")
            print()

        # Retry sweep — anything that failed with a Cloudflare interstitial
        # often resolves with a longer cool-down. CF escalates difficulty on
        # rapid same-domain hits, so wait substantially longer before retrying.
        retry_idxs = [
            i
            for i, r in enumerate(rows)
            if r["outcome"] != "success"
            and "cloudflare" in (r.get("reason") or "").lower()
        ]
        if retry_idxs:
            print(f"\nRetrying {len(retry_idxs)} CF-blocked PMID(s) after cool-down…\n")
            for i in retry_idxs:
                pmid = rows[i]["pmid"]
                print(f"--- retry PMID {pmid} ---")
                new_row = _run_with_pacing(pmid, gap_seconds=45.0)
                if new_row["outcome"] == "success":
                    rows[i] = new_row
                    print(
                        f"  retry succeeded: chars={new_row['chars']} "
                        f"reason={new_row['reason']}"
                    )
                else:
                    rows[i]["notes"] = (
                        rows[i].get("notes") or ""
                    ) + " retry_also_failed"
                    print(f"  retry failed: {new_row.get('reason')}")
                print()
    finally:
        try:
            pool.close()
        except Exception:
            pass

    # ---- 6. Final summary ----
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    by_publisher: Dict[str, Dict[str, int]] = {}
    for r in rows:
        p = r["strategy"] or "(none)"
        by_publisher.setdefault(p, {"success": 0, "fail": 0})
        if r["outcome"] in FETCH_SUCCESS_OUTCOMES:
            by_publisher[p]["success"] += 1
        else:
            by_publisher[p]["fail"] += 1
    for p, c in sorted(by_publisher.items()):
        print(f"  {p:14s}  success={c['success']}  fail={c['fail']}")

    success = sum(1 for r in rows if r["outcome"] in FETCH_SUCCESS_OUTCOMES)
    print(f"\nOverall: {success}/{len(rows)} PMIDs succeeded.")

    # Persist a CSV-ish JSON for downstream inspection.
    (output_dir / "summary.json").write_text(
        json.dumps(rows, indent=2, default=str), encoding="utf-8"
    )
    print(f"Wrote {output_dir / 'summary.json'}")

    return 0 if success == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
