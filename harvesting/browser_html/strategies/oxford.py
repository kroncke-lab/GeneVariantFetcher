"""Oxford Academic (OUP) — Europace, EHJ, Human Molecular Genetics, etc.

DOI prefix 10.1093. Many OUP journals make full HTML free 12 months after
publication. Body lives in ``.widget-ArticleFulltext``; supplements use
``/downloadSupplement?file=...``.
"""

from __future__ import annotations

from typing import Any

from ..base import FetchContext, FetchResult, PublisherStrategy
from . import register


@register
class OxfordStrategy(PublisherStrategy):
    NAME = "oxford"
    DOI_PREFIXES = ("10.1093",)
    DOMAINS = ("academic.oup.com",)
    EMBARGO_MONTHS = 12

    BODY_SELECTORS = (
        ".widget-ArticleFulltext",
        ".article-body",
        ".content-box-content",
    )

    def fetch(self, page: Any, ctx: FetchContext) -> FetchResult:
        result = FetchResult(publisher=self.NAME)
        if not ctx.doi:
            result.error = "no doi"
            return result

        target = (
            "https://academic.oup.com/lookup/doi/"
            f"{self.encode_doi_for_path(ctx.doi)}"
        )
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

        markdown = self.extract_via_scraper(
            html, final_url, ctx, selectors=list(self.BODY_SELECTORS)
        )
        if markdown:
            result.main_markdown = markdown

        # OUP's supplement scraping is already handled in the existing scraper.
        try:
            result.supp_files = ctx.scraper.scrape_oxford_supplements(html, final_url)
        except Exception as e:
            result.notes.append(f"oxford supp scrape error: {e}")

        try:
            result.figure_paths = self.download_figures(page, ctx)
        except Exception as e:
            result.notes.append(f"figure download error: {e}")

        return result
