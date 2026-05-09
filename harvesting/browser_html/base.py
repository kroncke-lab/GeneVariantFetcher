"""Abstract base for Tier 3.5 publisher strategies and shared dataclasses."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class FetchResult:
    """Result of a Tier 3.5 browser HTML fetch.

    Shape matches what the orchestrator expects to plug into the existing
    write_free_text_output / process_supplement_files pipeline.
    """

    main_markdown: Optional[str] = None
    main_html: Optional[str] = None
    supp_files: List[dict] = field(default_factory=list)
    figure_paths: List[Path] = field(default_factory=list)
    final_url: Optional[str] = None
    publisher: str = ""
    embargo_passed: bool = True
    notes: List[str] = field(default_factory=list)
    error: Optional[str] = None

    def is_usable(self) -> bool:
        """Did we actually get content worth writing?"""
        return bool(self.main_markdown and len(self.main_markdown) > 500)


@dataclass
class FetchContext:
    """Per-fetch context handed to a strategy.

    Carries the shared session, scraper, converter, and output dir so
    strategies don't have to know how to instantiate them.
    """

    pmid: str
    doi: str
    output_dir: Path
    scraper: Any  # harvesting.supplement_scraper.SupplementScraper
    converter: Any  # harvesting.format_converters.FormatConverter
    session: Any  # requests.Session (for figure downloads)
    timeout_s: int = 90


class PublisherStrategy(ABC):
    """One publisher's fetch policy.

    Subclasses declare:

    - ``DOI_PREFIXES``: tuple of DOI prefixes this strategy claims
      (e.g. ``("10.1161",)`` for AHA journals)
    - ``DOMAINS``: tuple of substrings; if a final URL contains any, this
      strategy can also claim it post-redirect
    - ``EMBARGO_MONTHS``: months that must elapse between pub date and now
      before the publisher's HTML is freely readable. ``0`` means no embargo;
      ``None`` means "unknown — always try"
    - ``NAME``: short slug used in logs and the publisher allowlist setting
    """

    NAME: str = "base"
    DOI_PREFIXES: tuple = ()
    DOMAINS: tuple = ()
    EMBARGO_MONTHS: Optional[int] = None

    def matches(self, doi: str, url: str = "") -> bool:
        """Does this strategy claim the given DOI/URL?"""
        if doi:
            doi_l = doi.lower().strip()
            for prefix in self.DOI_PREFIXES:
                if doi_l.startswith(prefix.lower()):
                    return True
        if url:
            url_l = url.lower()
            for dom in self.DOMAINS:
                if dom.lower() in url_l:
                    return True
        return False

    @abstractmethod
    def fetch(self, page: Any, ctx: FetchContext) -> FetchResult:
        """Drive the browser to retrieve full text + supplements + figures.

        Args:
            page: A Playwright ``Page`` ready for navigation.
            ctx: Shared fetch context.

        Returns:
            FetchResult populated with whatever could be retrieved.
        """

    # ------------------------------------------------------------------
    # Helpers concrete strategies can use
    # ------------------------------------------------------------------

    def accept_cookies(self, page: Any, selectors: Optional[List[str]] = None) -> bool:
        """Click a cookie-consent button if present. Returns True on click."""
        default_selectors = [
            "#onetrust-accept-btn-handler",
            "button#onetrust-accept-btn-handler",
            "button[aria-label*='Accept' i]",
            "button:has-text('Accept all cookies')",
            "button:has-text('Accept All Cookies')",
            "button:has-text('Accept all')",
            "button:has-text('I Accept')",
            "button:has-text('Got it')",
        ]
        for sel in (selectors or []) + default_selectors:
            try:
                el = page.query_selector(sel)
                if el and el.is_visible():
                    el.click(timeout=2000)
                    page.wait_for_timeout(500)
                    return True
            except Exception:
                continue
        return False

    def wait_for_body(self, page: Any, selector: str, timeout_ms: int = 15000) -> bool:
        """Wait for an article body selector. Returns True if it appeared."""
        try:
            page.wait_for_selector(selector, timeout=timeout_ms)
            return True
        except Exception:
            return False

    def extract_via_scraper(
        self, html: str, final_url: str, ctx: FetchContext
    ) -> Optional[str]:
        """Run the existing publisher-aware scraper on rendered HTML.

        Returns markdown or None if extraction failed. Logs the failure cause
        so silent 0-byte FULL_CONTEXT.md outcomes (observed on cohort-paper
        URLs whose page loads but body extraction fails) are visible.
        """
        try:
            markdown, _title = ctx.scraper.extract_fulltext(html, final_url)
            if not markdown:
                logger.info(
                    "extract_via_scraper: PMID %s scraper returned empty markdown for %s",
                    ctx.pmid,
                    final_url,
                )
            return markdown or None
        except Exception as e:
            logger.warning(
                "extract_via_scraper: PMID %s exception for %s: %s",
                ctx.pmid,
                final_url,
                e,
            )
            return None

    def download_figures(
        self,
        page: Any,
        ctx: FetchContext,
        img_selector: str = "figure img",
        max_figures: int = 30,
    ) -> List[Path]:
        """Download figure images via the shared requests session.

        We use the session (not Playwright) so cookies/headers established
        during page load can be reused, and so PNG/JPG validation matches
        what the rest of the pipeline expects.
        """
        figures_dir = ctx.output_dir / f"{ctx.pmid}_figures"
        figures_dir.mkdir(exist_ok=True)
        downloaded: List[Path] = []

        try:
            elements = page.query_selector_all(img_selector)
        except Exception:
            return downloaded

        seen_urls: set = set()
        for idx, el in enumerate(elements[:max_figures], 1):
            try:
                src = (
                    el.get_attribute("src")
                    or el.get_attribute("data-src")
                    or el.get_attribute("data-lazy-src")
                    or ""
                )
                if not src or src in seen_urls:
                    continue
                # Skip data URIs and tiny icons
                if src.startswith("data:") or "icon" in src.lower():
                    continue
                seen_urls.add(src)

                # Resolve relative URLs
                if src.startswith("//"):
                    src = "https:" + src
                elif src.startswith("/"):
                    from urllib.parse import urlparse

                    parsed = urlparse(page.url)
                    src = f"{parsed.scheme}://{parsed.netloc}{src}"

                ext = ".jpg"
                for candidate in (".png", ".gif", ".jpeg", ".jpg", ".svg"):
                    if candidate in src.lower():
                        ext = candidate
                        break
                # Skip SVG icons
                if ext == ".svg":
                    continue

                target = figures_dir / f"fig_browser_{idx}{ext}"
                resp = ctx.session.get(src, timeout=30)
                if resp.status_code != 200 or len(resp.content) < 2048:
                    continue
                magic = resp.content[:4]
                if magic[:3] in (b"\x89PN", b"GIF") or magic[:2] in (b"\xff\xd8",):
                    target.write_bytes(resp.content)
                    downloaded.append(target)
            except Exception:
                continue

        return downloaded
