"""Lightweight DOM → markdown helper for Tier 3.5 strategies.

Publisher HTML structures shift often; rather than chase every CSS class
change in the downstream ``SupplementScraper`` extractors, this module gives
each strategy a self-contained way to pick a body container and emit
reasonably clean markdown.

Kept small on purpose — `markitdown` and `pdfminer.six` already cover the
heavy lifting elsewhere; here we just need to lift article text out of a
loaded DOM with structure intact (headings, paragraphs, lists, tables).
"""

from __future__ import annotations

import re
from typing import Iterable, Optional, Sequence

from bs4 import BeautifulSoup, NavigableString, Tag

# Tags whose contents are not article body. Stripped before extraction so
# nav menus, citation popups, and figure-share widgets don't pollute output.
_STRIP_TAGS = (
    "script",
    "style",
    "noscript",
    "nav",
    "header",
    "footer",
    "form",
    "iframe",
    "svg",
    "button",
    "aside",
)

_STRIP_CLASS_CONTAINING = (
    "advert",
    "advertising",
    "altmetric",
    "social-share",
    "share-tools",
    "article-tools",
    "core-fv",
    "cookie",
    "menu",
    "skip-link",
    "back-to-top",
    "comments",
    "related-content",
)


def _is_strippable(el: Tag) -> bool:
    if el.name in _STRIP_TAGS:
        return True
    cls = " ".join(el.get("class") or []).lower()
    if not cls:
        return False
    return any(token in cls for token in _STRIP_CLASS_CONTAINING)


def _strip_chrome(container: Tag) -> None:
    for el in list(container.find_all(True)):
        try:
            if _is_strippable(el):
                el.decompose()
        except Exception:
            continue


def _text_of(el) -> str:
    if el is None:
        return ""
    return el.get_text(" ", strip=True)


def _table_to_markdown(table: Tag) -> str:
    """Quick markdown rendering of a <table>. Falls back to flat text on weird DOMs."""
    rows: list[list[str]] = []
    for tr in table.find_all("tr"):
        cells = [_text_of(c) for c in tr.find_all(["th", "td"])]
        if cells:
            rows.append(cells)
    if not rows:
        return _text_of(table)
    width = max(len(r) for r in rows)
    rows = [r + [""] * (width - len(r)) for r in rows]
    lines = ["| " + " | ".join(rows[0]) + " |"]
    lines.append("| " + " | ".join("---" for _ in range(width)) + " |")
    for r in rows[1:]:
        lines.append("| " + " | ".join(c.replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _render(el: Tag, depth: int = 0) -> str:
    """Recursive markdown render for a body container."""
    parts: list[str] = []
    for child in el.children:
        if isinstance(child, NavigableString):
            txt = str(child)
            if txt.strip():
                parts.append(txt)
            continue
        if not isinstance(child, Tag):
            continue
        name = child.name
        if name in _STRIP_TAGS:
            continue
        if name in ("h1", "h2", "h3", "h4", "h5", "h6"):
            level = int(name[1])
            parts.append("\n\n" + "#" * level + " " + _text_of(child) + "\n")
        elif name == "p":
            txt = _text_of(child)
            if txt:
                parts.append("\n\n" + txt + "\n")
        elif name in ("ul", "ol"):
            for li in child.find_all("li", recursive=False):
                parts.append("\n- " + _text_of(li))
            parts.append("\n")
        elif name == "table":
            parts.append("\n\n" + _table_to_markdown(child) + "\n")
        elif name == "figure":
            cap = child.find("figcaption")
            if cap:
                parts.append("\n\n*Figure: " + _text_of(cap) + "*\n")
        elif name == "section" or name == "div" or name == "article":
            parts.append(_render(child, depth + 1))
        elif name == "br":
            parts.append("\n")
        else:
            # Inline elements — concatenate their text content.
            txt = _text_of(child)
            if txt:
                parts.append(txt)
    return "".join(parts)


def extract_body_markdown(
    html: str,
    selectors: Sequence[str],
    title_selectors: Sequence[str] = ("h1.article-title", "h1.citation__title", "h1"),
) -> Optional[str]:
    """Pick the first matching body selector and render it to markdown.

    Returns None if no selector matches a non-empty container.
    """
    if not html:
        return None
    soup = BeautifulSoup(html, "html.parser")

    title = ""
    for ts in title_selectors:
        t = soup.select_one(ts)
        if t and _text_of(t):
            title = _text_of(t)
            break

    container: Optional[Tag] = None
    chosen_selector = ""
    for sel in selectors:
        for cand in soup.select(sel):
            if not isinstance(cand, Tag):
                continue
            text_len = len(_text_of(cand))
            if text_len < 200:
                continue
            container = cand
            chosen_selector = sel
            break
        if container is not None:
            break
    if container is None:
        return None

    _strip_chrome(container)
    body_md = _render(container).strip()

    # Collapse runs of >2 blank lines to keep the markdown readable.
    body_md = re.sub(r"\n{3,}", "\n\n", body_md)

    parts = ["# MAIN TEXT", ""]
    if title:
        parts.extend(["## " + title, ""])
    parts.append(body_md)
    return "\n".join(parts).strip()


def pick_better_markdown(
    primary: Optional[str],
    dom: Optional[str],
    threshold: float = 0.80,
) -> Optional[str]:
    """Choose between two markdown extractions of the same page.

    Strategies typically have two extraction paths: the legacy
    publisher-aware ``SupplementScraper.extract_fulltext`` chain (high
    quality when its class-name targets match the current HTML, broken when
    they don't), and the DOM-walker ``extract_body_markdown`` (selector-
    agnostic, sometimes pulls duplicate citation lists or nav residue).

    Picking the longer string indiscriminately is wrong — primary's chrome-
    stripping is often better. But when primary returns substantially less
    text than DOM, the legacy selectors have rotted and DOM is the only
    path that recovered the article body. Default threshold: prefer DOM
    when primary is under 80% of DOM's length.
    """
    if not primary:
        return dom
    if not dom:
        return primary
    if len(primary) < threshold * len(dom):
        return dom
    return primary


def looks_like_cloudflare_challenge(html: str) -> bool:
    """Quick check whether the page is the CF JS-challenge interstitial."""
    if not html:
        return False
    h = html.lower()
    return (
        "checking your browser" in h
        or ("cf-browser-verification" in h)
        or ("lds-ring" in h and "loading-verifying" in h)
        or ("ray-id" in h and "challenge-platform" in h)
    )


def looks_like_paywall_stub(html: str) -> bool:
    """Quick check whether the page is a "we couldn't authenticate you" stub."""
    if not html:
        return False
    h = html.lower()
    indicators = (
        "this content requires a subscription",
        "purchase this article",
        "buy this article",
        "you do not have access",
        "your access options",
    )
    hits = sum(1 for ind in indicators if ind in h)
    return hits >= 2
