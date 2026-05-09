"""Generic fallback strategy for publishers without a dedicated handler.

Navigates to ``https://doi.org/{doi}``, follows redirects, runs the existing
``SupplementScraper.extract_fulltext`` and ``scrape_generic_supplements`` on
the rendered HTML.
"""

from __future__ import annotations

from typing import Any

from ..base import FetchContext, FetchResult, PublisherStrategy
from . import register


@register
class GenericStrategy(PublisherStrategy):
    NAME = "generic"
    DOI_PREFIXES = ()
    DOMAINS = ()
    EMBARGO_MONTHS = None  # don't gate — try and let the scraper validate

    def fetch(self, page: Any, ctx: FetchContext) -> FetchResult:
        result = FetchResult(publisher=self.NAME)

        target = f"https://doi.org/{ctx.doi}" if ctx.doi else None
        if not target:
            result.error = "no doi"
            return result

        try:
            page.goto(target, wait_until="load", timeout=ctx.timeout_s * 1000)
        except Exception as e:
            result.error = f"navigation failed: {e}"
            return result

        self.accept_cookies(page)

        try:
            page.wait_for_timeout(1500)  # let lazy content settle
            html = page.content()
            final_url = page.url
        except Exception as e:
            result.error = f"could not read page: {e}"
            return result

        result.main_html = html
        result.final_url = final_url

        # Reuse existing publisher-aware extractor — it routes by domain.
        markdown = self.extract_via_scraper(html, final_url, ctx)
        if markdown:
            result.main_markdown = markdown
        else:
            result.notes.append(
                f"extract_via_scraper returned empty markdown for {final_url}"
            )

        # Try the existing publisher supplement scrapers based on domain.
        try:
            domain = final_url.lower()
            if "nature.com" in domain:
                result.supp_files = ctx.scraper.scrape_nature_supplements(
                    html, final_url
                )
            elif any(
                d in domain
                for d in ("sciencedirect.com", "elsevier.com", "gimjournal.org")
            ):
                result.supp_files = ctx.scraper.scrape_elsevier_supplements(
                    html, final_url
                )
            elif "wiley.com" in domain or "onlinelibrary.wiley.com" in domain:
                result.supp_files = ctx.scraper.scrape_wiley_supplements(
                    html, final_url
                )
            elif "academic.oup.com" in domain:
                result.supp_files = ctx.scraper.scrape_oxford_supplements(
                    html, final_url
                )
            elif any(
                d in domain
                for d in ("link.springer.com", "biomedcentral.com", "springeropen.com")
            ):
                result.supp_files = ctx.scraper.scrape_springer_supplements(
                    html, final_url
                )
            else:
                result.supp_files = ctx.scraper.scrape_generic_supplements(
                    html, final_url
                )
        except Exception as e:
            result.notes.append(f"supplement scrape error: {e}")

        # Best-effort figure capture.
        try:
            result.figure_paths = self.download_figures(page, ctx)
        except Exception as e:
            result.notes.append(f"figure download error: {e}")

        return result
