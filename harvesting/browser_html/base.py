"""Abstract base for Tier 3.5 publisher strategies and shared dataclasses."""

from __future__ import annotations

import logging
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Optional
from urllib.parse import quote

from bs4 import BeautifulSoup

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

    @staticmethod
    def encode_doi_for_path(doi: str) -> str:
        """Percent-encode a DOI for publisher path URLs without escaping slashes."""
        return quote((doi or "").strip(), safe="/:")

    @staticmethod
    def _clean_text(value: str) -> str:
        text = re.sub(r"\s+", " ", value or "").strip()
        return text.replace("|", "\\|")

    def _table_to_markdown(self, table: Any) -> str:
        """Convert an HTML table to markdown while preserving row/cell order."""
        rows: List[List[str]] = []
        for tr in table.find_all("tr"):
            cells = [
                self._clean_text(cell.get_text(" ", strip=True))
                for cell in tr.find_all(["th", "td"])
            ]
            if any(cells):
                rows.append(cells)
        if not rows:
            return ""

        max_cols = max(len(row) for row in rows)
        padded = [row + [""] * (max_cols - len(row)) for row in rows]
        header = padded[0]
        lines = [
            "| " + " | ".join(header) + " |",
            "| " + " | ".join(["---"] * max_cols) + " |",
        ]
        for row in padded[1:]:
            lines.append("| " + " | ".join(row) + " |")
        return "\n".join(lines)

    def extract_readable_html(
        self,
        html: str,
        selectors: Optional[List[str]] = None,
    ) -> Optional[str]:
        """Fallback HTML-to-markdown pass that keeps article tables.

        The legacy publisher scrapers mainly preserve paragraphs. Cohort
        papers often put the useful variants in HTML tables, so this pass
        walks the rendered article container and emits headings, paragraphs,
        lists, and tables in document order.
        """
        if not html:
            return None

        soup = BeautifulSoup(html, "html.parser")
        for tag in soup(["script", "style", "noscript", "nav", "footer", "header"]):
            tag.decompose()

        title = ""
        title_elem = soup.find("h1") or soup.find("title")
        if title_elem:
            title = self._clean_text(title_elem.get_text(" ", strip=True))

        candidates = []
        for selector in selectors or []:
            try:
                candidates.extend(soup.select(selector))
            except Exception:
                continue
        for selector in (
            "article",
            "main",
            "[role='main']",
            ".article__body",
            ".article-body",
            ".article-section__full",
            ".article-section__content",
            ".hlFld-Fulltext",
            ".widget-ArticleFulltext",
            "#body",
            ".Body",
        ):
            try:
                candidates.extend(soup.select(selector))
            except Exception:
                continue
        if not candidates and soup.body:
            candidates = [soup.body]

        if not candidates:
            return None

        def score(node: Any) -> int:
            text_len = len(node.get_text(" ", strip=True))
            table_bonus = 1200 * len(node.find_all("table"))
            return text_len + table_bonus

        root = max(candidates, key=score)
        parts: List[str] = ["# MAIN TEXT"]
        if title:
            parts.append(f"## {title}")

        skip_phrases = (
            "cookie",
            "privacy policy",
            "sign in",
            "log in",
            "subscribe",
            "register",
            "your browser",
            "javascript",
            "advertisement",
        )

        for node in root.find_all(["h2", "h3", "h4", "h5", "h6", "p", "li", "table"]):
            if node.find_parent("table") and node.name != "table":
                continue
            if node.name == "table":
                caption = node.find("caption")
                cap_text = (
                    self._clean_text(caption.get_text(" ", strip=True))
                    if caption
                    else ""
                )
                table_md = self._table_to_markdown(node)
                if table_md:
                    if cap_text:
                        parts.append(f"### Table: {cap_text}")
                    parts.append(table_md)
                continue

            text = self._clean_text(node.get_text(" ", strip=True))
            if not text:
                continue
            lower = text.lower()
            if node.name == "p" and len(text) < 30:
                continue
            if any(phrase in lower for phrase in skip_phrases) and len(text) < 300:
                continue
            if node.name and node.name.startswith("h"):
                if len(text) < 140:
                    parts.append(f"### {text}")
            elif node.name == "li":
                if len(text) > 10:
                    parts.append(f"- {text}")
            else:
                parts.append(text)

        markdown = "\n\n".join(parts).strip() + "\n"
        return markdown if len(markdown) > 300 else None

    def extract_via_scraper(
        self,
        html: str,
        final_url: str,
        ctx: FetchContext,
        selectors: Optional[List[str]] = None,
    ) -> Optional[str]:
        """Run the existing publisher-aware scraper on rendered HTML.

        Returns markdown or None if extraction failed. Logs the failure cause
        so silent 0-byte FULL_CONTEXT.md outcomes (observed on cohort-paper
        URLs whose page loads but body extraction fails) are visible.
        """
        fallback_markdown = self.extract_readable_html(html, selectors)
        try:
            markdown, _title = ctx.scraper.extract_fulltext(html, final_url)
            if not markdown:
                logger.info(
                    "extract_via_scraper: PMID %s scraper returned empty markdown for %s",
                    ctx.pmid,
                    final_url,
                )
            if markdown and fallback_markdown:
                if len(fallback_markdown) > max(
                    len(markdown) * 1.2, len(markdown) + 1000
                ):
                    logger.info(
                        "extract_via_scraper: PMID %s using table-preserving HTML fallback "
                        "(%d chars vs %d scraper chars)",
                        ctx.pmid,
                        len(fallback_markdown),
                        len(markdown),
                    )
                    return fallback_markdown
            return markdown or fallback_markdown or None
        except Exception as e:
            logger.warning(
                "extract_via_scraper: PMID %s exception for %s: %s",
                ctx.pmid,
                final_url,
                e,
            )
            return fallback_markdown

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
