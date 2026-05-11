"""Elsevier "Open Archive" via ScienceDirect — DOI prefix 10.1016.

A subset of Elsevier journals make older content freely readable. Body
content is heavily JS-rendered; we wait for ``article#main-content`` or the
older ``.Body`` selector before scraping.

Imprint sites (Heart Rhythm, JACC, Cell, Lancet, …) share Elsevier
infrastructure but land on ``/article/{pii}/abstract`` by default; the
fulltext URL is the same path with ``/fulltext``. We rewrite once after the
initial redirect so we don't waste cycles parsing an abstract.
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
class ElsevierOpenStrategy(PublisherStrategy):
    NAME = "elsevier_open"
    DOI_PREFIXES = ("10.1016",)
    DOMAINS = (
        "sciencedirect.com",
        "linkinghub.elsevier.com",
        "heartrhythmjournal.com",
        "jacc.org",
        "cell.com",
        "thelancet.com",
    )
    # Open Archive policies vary by journal — leave embargo unset so the
    # strategy attempts and bails out fast on access denial.
    EMBARGO_MONTHS = None

    BODY_SELECTORS = (
        "article#main-content",
        "div.Body",
        "article#article",
        "section.article__sections",
        "section[aria-labelledby='section-cited-by']",
        "article",
        "main",
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

        # Imprint sites (Heart Rhythm, JACC, Cell, Lancet) land on
        # ``/article/.../abstract`` after the DOI redirect — that page has
        # only ~500 chars of body. Rewriting to ``/fulltext`` retrieves the
        # full article when the user has institutional access.
        try:
            current = page.url
            if "/abstract" in current and "/fulltext" not in current:
                fulltext_url = current.replace("/abstract", "/fulltext")
                try:
                    page.goto(
                        fulltext_url,
                        wait_until="load",
                        timeout=ctx.timeout_s * 1000,
                    )
                except Exception:
                    # Some titles only honor /fulltext when authorized;
                    # fall back to the abstract page (still useful for diag).
                    pass
        except Exception:
            pass

        for sel in self.BODY_SELECTORS:
            if self.wait_for_body(page, sel, timeout_ms=10000):
                break
        # ScienceDirect lazy-loads sections on scroll.
        try:
            for _ in range(3):
                page.evaluate("window.scrollBy(0, document.body.scrollHeight)")
                page.wait_for_timeout(800)
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

        # CF / paywall signals (retry sweep keys on "cloudflare" substring).
        if looks_like_cloudflare_challenge(html):
            result.notes.append("cloudflare challenge unresolved")
            return result
        if "/abstract" in (final_url or "") and looks_like_paywall_stub(html):
            result.notes.append("redirected to /abstract — institutional access denied")
        elif looks_like_paywall_stub(html):
            result.notes.append("paywall stub detected")

        # Dual extraction; DOM walker wins when legacy class targets miss.
        primary_md = self.extract_via_scraper(
            html, final_url, ctx, selectors=list(self.BODY_SELECTORS)
        )
        dom_md = extract_body_markdown(html, self.BODY_SELECTORS)
        result.main_markdown = pick_better_markdown(primary_md, dom_md)

        try:
            result.supp_files = ctx.scraper.scrape_elsevier_supplements(html, final_url)
        except Exception as e:
            result.notes.append(f"elsevier supp scrape error: {e}")

        try:
            result.figure_paths = self.download_figures(page, ctx)
        except Exception as e:
            result.notes.append(f"figure download error: {e}")

        return result
