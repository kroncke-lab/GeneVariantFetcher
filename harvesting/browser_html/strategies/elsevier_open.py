"""Elsevier "Open Archive" via ScienceDirect — DOI prefix 10.1016.

A subset of Elsevier journals make older content freely readable. Body
content is heavily JS-rendered; we wait for ``article#main-content`` or the
older ``.Body`` selector before scraping.
"""

from __future__ import annotations

from typing import Any

from ..base import FetchContext, FetchResult, PublisherStrategy
from . import register


@register
class ElsevierOpenStrategy(PublisherStrategy):
    NAME = "elsevier_open"
    DOI_PREFIXES = ("10.1016",)
    DOMAINS = ("sciencedirect.com", "linkinghub.elsevier.com")
    # Open Archive policies vary by journal — leave embargo unset so the
    # strategy attempts and bails out fast on access denial.
    EMBARGO_MONTHS = None

    BODY_SELECTORS = (
        "article#main-content",
        "div.Body",
        "section[aria-labelledby='section-cited-by']",
        "article",
    )

    def fetch(self, page: Any, ctx: FetchContext) -> FetchResult:
        result = FetchResult(publisher=self.NAME)
        if not ctx.doi:
            result.error = "no doi"
            return result

        target = f"https://doi.org/{self.encode_doi_for_path(ctx.doi)}"
        try:
            page.goto(target, wait_until="load", timeout=ctx.timeout_s * 1000)
        except Exception as e:
            result.error = f"navigation failed: {e}"
            return result

        self.accept_cookies(page)

        for sel in self.BODY_SELECTORS:
            if self.wait_for_body(page, sel, timeout_ms=10000):
                break
        # ScienceDirect lazy-loads sections on scroll.
        try:
            for _ in range(3):
                page.evaluate("window.scrollBy(0, document.body.scrollHeight)")
                page.wait_for_timeout(800)
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

        markdown = self.extract_via_scraper(
            html, final_url, ctx, selectors=list(self.BODY_SELECTORS)
        )
        if markdown:
            result.main_markdown = markdown

        try:
            result.supp_files = ctx.scraper.scrape_elsevier_supplements(html, final_url)
        except Exception as e:
            result.notes.append(f"elsevier supp scrape error: {e}")

        try:
            result.figure_paths = self.download_figures(page, ctx)
        except Exception as e:
            result.notes.append(f"figure download error: {e}")

        return result
