"""Wiley Online Library — DOI prefixes 10.1002, 10.1111.

Some Wiley titles (notably Human Mutation under 10.1002) are open access
with no embargo; others gate by subscription. We try and let the content
validator decide.
"""

from __future__ import annotations

from typing import Any

from ..base import FetchContext, FetchResult, PublisherStrategy
from . import register


@register
class WileyStrategy(PublisherStrategy):
    NAME = "wiley"
    DOI_PREFIXES = ("10.1002", "10.1111")
    DOMAINS = ("onlinelibrary.wiley.com",)
    EMBARGO_MONTHS = 0  # try immediately; OA Wiley journals have no embargo

    BODY_SELECTORS = (
        ".article-section__full",
        "section.article__body",
        "article",
    )

    def fetch(self, page: Any, ctx: FetchContext) -> FetchResult:
        result = FetchResult(publisher=self.NAME)
        if not ctx.doi:
            result.error = "no doi"
            return result

        target = f"https://onlinelibrary.wiley.com/doi/full/{ctx.doi}"
        try:
            page.goto(target, wait_until="load", timeout=ctx.timeout_s * 1000)
        except Exception as e:
            result.error = f"navigation failed: {e}"
            return result

        self.accept_cookies(page)

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

        try:
            result.supp_files = ctx.scraper.scrape_wiley_supplements(html, final_url)
        except Exception as e:
            result.notes.append(f"wiley supp scrape error: {e}")

        try:
            result.figure_paths = self.download_figures(page, ctx)
        except Exception as e:
            result.notes.append(f"figure download error: {e}")

        return result
