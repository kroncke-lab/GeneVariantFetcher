"""Wiley Online Library — DOI prefixes 10.1002, 10.1111.

Some Wiley titles (notably Human Mutation under 10.1002) are open access
with no embargo; others gate by subscription. Wiley fronts content with
Cloudflare, so we run the same CF wait/reload loop the AHA strategy uses.
"""

from __future__ import annotations

from typing import Any

from ..base import FetchContext, FetchResult, PublisherStrategy
from ..dom_extract import (
    extract_body_markdown,
    looks_like_cloudflare_challenge,
    looks_like_paywall_stub,
    pick_better_markdown,
)
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
        "main",
    )

    def fetch(self, page: Any, ctx: FetchContext) -> FetchResult:
        result = FetchResult(publisher=self.NAME)
        if not ctx.doi:
            result.error = "no doi"
            return result

        target = (
            "https://onlinelibrary.wiley.com/doi/full/"
            f"{self.encode_doi_for_path(ctx.doi)}"
        )
        # Route through the institutional EZproxy when configured so the request
        # egresses from the publisher-allowlisted subscriber IP (clears the
        # Cloudflare managed challenge). No-op when GVF_EZPROXY_* is unset.
        from harvesting.browser_html import ezproxy

        target = ezproxy.wrap(target)
        try:
            page.goto(target, wait_until="load", timeout=ctx.timeout_s * 1000)
        except Exception as e:
            result.error = f"navigation failed: {e}"
            return result

        # Cloudflare interstitial loop (Wiley uses CF aggressively).
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
                    ".article-section__full, section.article__body, article",
                    timeout=8000,
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

        # CF + paywall signals (retry sweep keys on "cloudflare").
        if looks_like_cloudflare_challenge(html):
            result.notes.append("cloudflare challenge unresolved")
            return result
        if looks_like_paywall_stub(html):
            result.notes.append("paywall stub detected")

        # Dual extraction: legacy publisher-aware path + DOM walker.
        primary_md = self.extract_via_scraper(
            html, final_url, ctx, selectors=list(self.BODY_SELECTORS)
        )
        dom_md = extract_body_markdown(html, self.BODY_SELECTORS)
        result.main_markdown = pick_better_markdown(primary_md, dom_md)

        try:
            result.supp_files = ctx.scraper.scrape_wiley_supplements(html, final_url)
        except Exception as e:
            result.notes.append(f"wiley supp scrape error: {e}")

        try:
            result.figure_paths = self.download_figures(page, ctx)
        except Exception as e:
            result.notes.append(f"figure download error: {e}")

        return result
