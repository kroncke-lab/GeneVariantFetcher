"""AHA Journals (Circulation, JAHA, Stroke, etc.) — DOI prefix 10.1161.

AHA HTML is freely readable 12 months after publication. Body text lives in
``.hlFld-Fulltext``; supplements are linked under ``/doi/suppl/10.1161/...``
and not exposed via any AHA API.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List
from urllib.parse import urljoin, urlparse

from bs4 import BeautifulSoup

from ..base import FetchContext, FetchResult, PublisherStrategy
from . import register

logger = logging.getLogger(__name__)


@register
class AHAStrategy(PublisherStrategy):
    NAME = "aha"
    DOI_PREFIXES = ("10.1161",)
    DOMAINS = ("ahajournals.org",)
    EMBARGO_MONTHS = 12

    BODY_SELECTORS = (".hlFld-Fulltext", ".article__body", "article")

    def fetch(self, page: Any, ctx: FetchContext) -> FetchResult:
        result = FetchResult(publisher=self.NAME)
        if not ctx.doi:
            result.error = "no doi"
            return result

        target = f"https://www.ahajournals.org/doi/{ctx.doi}"
        try:
            page.goto(target, wait_until="load", timeout=ctx.timeout_s * 1000)
        except Exception as e:
            result.error = f"navigation failed: {e}"
            return result

        self.accept_cookies(page)

        # Wait for any of the known body selectors. Don't fail if missing —
        # some article variants render outside the standard wrapper and the
        # generic extractor can still recover text.
        for sel in self.BODY_SELECTORS:
            if self.wait_for_body(page, sel, timeout_ms=8000):
                break

        try:
            html = page.content()
            final_url = page.url
        except Exception as e:
            result.error = f"could not read page: {e}"
            return result

        result.main_html = html
        result.final_url = final_url

        markdown = self.extract_via_scraper(html, final_url, ctx)
        if markdown:
            result.main_markdown = markdown

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
