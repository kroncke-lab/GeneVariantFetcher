"""Generic fallback strategy for publishers without a dedicated handler.

Navigates to ``https://doi.org/{doi}``, follows redirects, then runs both
the legacy ``SupplementScraper.extract_fulltext`` chain and the DOM-walker
``extract_body_markdown``. ``pick_better_markdown`` chooses between them so
publishers whose HTML class names have rotted (Sage, Karger, smaller titles)
still recover article body via the selector-agnostic DOM path.
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

# Wide body-selector set covering publishers that fall through to this
# strategy — Karger, Springer, BMC, Sage, smaller Wiley imprints, etc.
# ``extract_body_markdown`` rejects containers under 200 chars so a stub
# matching here doesn't produce false positives.
GENERIC_BODY_SELECTORS = (
    "div.article-body",
    "div.fulltext-view",
    "section.article-section",
    "article#article",
    "article.article",
    "main#main-content",
    "#bodymatter",
    "article",
    "main",
)


@register
class GenericStrategy(PublisherStrategy):
    NAME = "generic"
    DOI_PREFIXES = ()
    DOMAINS = ()
    EMBARGO_MONTHS = None  # don't gate — try and let the scraper validate

    def fetch(self, page: Any, ctx: FetchContext) -> FetchResult:
        result = FetchResult(publisher=self.NAME)

        target = (
            f"https://doi.org/{self.encode_doi_for_path(ctx.doi)}" if ctx.doi else None
        )
        if not target:
            result.error = "no doi"
            return result

        try:
            page.goto(target, wait_until="load", timeout=ctx.timeout_s * 1000)
        except Exception as e:
            result.error = f"navigation failed: {e}"
            return result

        self.accept_cookies(page)

        # Cloudflare interstitial handling. Many smaller publishers (Sage,
        # Karger, BMJ, etc.) front content with CF. Wait for an article-
        # level container to appear before snapshotting; reload on cycles
        # 1 and 2 the way the AHA strategy does.
        for cycle in range(3):
            try:
                html_probe = page.content()
            except Exception:
                html_probe = ""
            if not looks_like_cloudflare_challenge(html_probe):
                break
            page.wait_for_timeout(5000)
            try:
                page.wait_for_selector(
                    "article, main, #bodymatter, div.article-body", timeout=8000
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

        try:
            page.wait_for_timeout(1500)  # let lazy content settle
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

        # Dual extraction: legacy domain-routed scraper + DOM walker.
        # The DOM walker wins on Sage and Karger where the legacy
        # publisher selectors don't match.
        primary_md = self.extract_via_scraper(html, final_url, ctx)
        dom_md = extract_body_markdown(html, GENERIC_BODY_SELECTORS)
        chosen = pick_better_markdown(primary_md, dom_md)
        if chosen:
            result.main_markdown = chosen
        else:
            result.notes.append(f"no body container matched for {final_url}")

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
            elif "mayoclinicproceedings.org" in domain:
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
            elif "karger.com" in domain:
                result.supp_files = ctx.scraper.scrape_karger_supplements(
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
