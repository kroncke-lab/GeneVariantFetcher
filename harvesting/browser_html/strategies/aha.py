"""AHA Journals (Circulation, JAHA, Stroke, etc.) — DOI prefix 10.1161.

AHA HTML is freely readable 12 months after publication. The body has lived
in ``#bodymatter`` since AHA's 2024 redesign — the historical
``.hlFld-Fulltext`` wrapper now contains only the eLetters block (~8 chars).
Supplements are linked under ``/doi/suppl/10.1161/...`` and not exposed via
any AHA API.

AHA fronts content with Cloudflare, which fires a JS interstitial when it
sees rapid same-session requests. This strategy polls for the body
container across multiple cycles and reloads when the challenge sticks.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List
from urllib.parse import urljoin, urlparse

from bs4 import BeautifulSoup

from ..base import FetchContext, FetchResult, PublisherStrategy
from ..dom_extract import (
    extract_body_markdown,
    looks_like_cloudflare_challenge,
    looks_like_paywall_stub,
    pick_better_markdown,
)
from . import register

logger = logging.getLogger(__name__)


@register
class AHAStrategy(PublisherStrategy):
    NAME = "aha"
    DOI_PREFIXES = ("10.1161",)
    DOMAINS = ("ahajournals.org",)
    EMBARGO_MONTHS = 12

    # Current AHA selector is ``#bodymatter``; the legacy wrappers are kept
    # as fallbacks for older article variants. ``article`` / ``main`` are
    # last-resort generic containers that the DOM walker can still mine.
    BODY_SELECTORS = (
        "#bodymatter",
        ".article__body",
        ".hlFld-Fulltext",
        "article",
        "main",
    )

    # Subset used as wait-for hints. We only need one of these to appear
    # before we trust that the CF challenge cleared and the page has DOM.
    WAIT_SELECTORS = ("#bodymatter", "article", ".article__body")

    def fetch(self, page: Any, ctx: FetchContext) -> FetchResult:
        result = FetchResult(publisher=self.NAME)
        if not ctx.doi:
            result.error = "no doi"
            return result

        target = f"https://www.ahajournals.org/doi/{self.encode_doi_for_path(ctx.doi)}"
        try:
            page.goto(target, wait_until="load", timeout=ctx.timeout_s * 1000)
        except Exception as e:
            result.error = f"navigation failed: {e}"
            return result

        # Cloudflare interstitial loop. CF's JS challenge clears itself by
        # running JS, setting cookies, and reloading. We poll for the body
        # container and reissue the navigation on cycles 1 and 2 — that
        # covers the common case where CF expects a fresh GET after the JS
        # challenge sets its clearance cookie.
        for cycle in range(4):
            try:
                html_probe = page.content()
            except Exception:
                html_probe = ""
            if not looks_like_cloudflare_challenge(html_probe):
                break
            page.wait_for_timeout(5000)
            try:
                page.wait_for_selector(
                    "#bodymatter, article, .article__body", timeout=8000
                )
            except Exception:
                if cycle in (1, 2):
                    try:
                        page.goto(
                            target,
                            wait_until="domcontentloaded",
                            timeout=ctx.timeout_s * 1000,
                        )
                    except Exception:
                        pass

        self.accept_cookies(page)

        # Wait for any of the body selectors. Don't fail if missing —
        # the DOM extractor still has wider fallbacks below.
        for sel in self.WAIT_SELECTORS:
            if self.wait_for_body(page, sel, timeout_ms=8000):
                break

        # AHA lazy-loads figures and section anchors. A short scroll
        # triggers deferred renders before we snapshot HTML.
        try:
            for _ in range(2):
                page.evaluate("window.scrollBy(0, document.body.scrollHeight)")
                page.wait_for_timeout(500)
            page.evaluate("window.scrollTo(0, 0)")
        except Exception:
            pass

        try:
            html = page.content()
            final_url = page.url
        except Exception as e:
            result.error = f"could not read page: {e}"
            return result

        result.main_html = html
        result.final_url = final_url

        # Surface CF and paywall signals before extraction so the caller's
        # retry sweep can route on the reason string (it looks for the
        # literal substring "cloudflare").
        if looks_like_cloudflare_challenge(html):
            result.notes.append("cloudflare challenge unresolved")
            return result
        if looks_like_paywall_stub(html):
            result.notes.append("paywall stub detected")
            # Continue — the page may still carry an abstract we want to
            # preserve for downstream debugging.

        # Dual extraction: legacy publisher-aware path + DOM walker.
        # ``pick_better_markdown`` prefers DOM when the legacy path's
        # class-name targets have rotted (the AHA 2024 redesign case).
        primary_md = self.extract_via_scraper(
            html, final_url, ctx, selectors=list(self.BODY_SELECTORS)
        )
        dom_md = extract_body_markdown(html, self.BODY_SELECTORS)
        result.main_markdown = pick_better_markdown(primary_md, dom_md)

        result.supp_files = self._scrape_aha_supplements(html, final_url)

        try:
            result.figure_paths = self.download_figures(page, ctx)
        except Exception as e:
            result.notes.append(f"figure download error: {e}")

        return result

    # ------------------------------------------------------------------
    # Self-contained AHA supplement scraping
    # ------------------------------------------------------------------

    def _scrape_aha_supplements(self, html: str, base_url: str) -> List[Dict]:
        """Find supplement links on AHA article pages.

        AHA links supplements via URLs containing ``/doi/suppl/`` and direct
        links to ``-supplemental-`` named PDFs/Excel files.
        """
        soup = BeautifulSoup(html, "html.parser")
        results: List[Dict] = []
        seen: set = set()
        valid_exts = (".pdf", ".docx", ".doc", ".xlsx", ".xls", ".csv", ".zip")

        # Look in supplement-specific sections first.
        candidates = []
        for selector in (
            "section[id*='supplementary' i] a[href]",
            "div[class*='supplementary' i] a[href]",
            "div[class*='supplement' i] a[href]",
            "a[href*='/doi/suppl/']",
            "a[href*='supplemental-material']",
        ):
            try:
                candidates.extend(soup.select(selector))
            except Exception:
                continue

        for a in candidates:
            href = (a.get("href") or "").strip()
            if not href or href.startswith("#"):
                continue
            full = urljoin(base_url, href)
            text = a.get_text(" ", strip=True)
            href_l = full.lower()

            looks_like_file = href_l.endswith(valid_exts)
            looks_like_supp = (
                "/suppl/" in href_l
                or "supplement" in (text + " " + href_l).lower()
                or "supporting" in text.lower()
            )
            if not (looks_like_file or looks_like_supp):
                continue
            if full in seen:
                continue
            seen.add(full)

            name = self._filename_from(full, text)
            results.append(
                {
                    "url": full,
                    "name": name,
                    "base_url": base_url,
                    "original_url": href,
                }
            )

        return results

    @staticmethod
    def _filename_from(url: str, text: str) -> str:
        path = urlparse(url).path
        last = path.rsplit("/", 1)[-1] if path else ""
        if last and "." in last:
            return last
        if text:
            cleaned = text.strip().replace("/", "_")[:80]
            return cleaned or "aha_supplement"
        return "aha_supplement"
