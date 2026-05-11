"""Content quality gate for Tier 3.5 results.

Authenticated fetches can return things that *look* like a full article but
are actually paywall stubs, "access denied" pages, or institutional-login
landing pages. The pipeline cannot tell the difference downstream, so we
reject obvious failures here.

Heuristics:

- Length floor: usable articles have at least a few KB of body text. The
  default 5_000-char floor is below the shortest real KCNH2 paper in our
  corpus but above any paywall stub we've seen.
- Negative markers: text phrases that only appear on paywall / SSO pages.
  We require at least 2 distinct markers AND a short body to reject — one
  marker alone shouldn't be enough since some real reviews quote them.
- Positive markers: phrases that strongly indicate full text (Methods,
  References section, etc.). If present, we don't reject for short length.
"""

from __future__ import annotations

import re
from typing import Tuple

# Phrases that, in clusters, indicate the page is a paywall/access stub
# rather than full text.
PAYWALL_MARKERS = (
    "purchase this article",
    "buy this article",
    "rent this article",
    "subscribe to access",
    "get access",
    "institutional sign in",
    "institutional login",
    "log in to your personal account",
    "your access options",
    "this content requires a subscription",
    "please sign in",
    "you do not have access",
    "access through your institution",
    "no access",
    "checking your browser before accessing",  # Cloudflare interstitial
    "verify you are human",
    "ddos protection by cloudflare",
)

# Phrases that strongly suggest full text was actually rendered.
FULLTEXT_MARKERS = (
    "references",
    "materials and methods",
    "patients and methods",
    "study population",
    "results",
    "discussion",
    "supplementary",
    "supporting information",
    "acknowledgments",
    "conflicts of interest",
)


def _count_distinct_markers(text_lower: str, markers: Tuple[str, ...]) -> int:
    return sum(1 for m in markers if m in text_lower)


def validate_article_content(
    text: str,
    min_chars: int = 5_000,
    min_words: int = 800,
) -> Tuple[bool, str]:
    """Return (ok, reason). ``reason`` describes why we accepted or rejected.

    Designed to be pluggable into ``BrowserHTMLFetcher(validate_content_quality=...)``.
    """
    if not text:
        return False, "empty body"

    body = text.strip()
    char_count = len(body)
    word_count = len(re.findall(r"\b\w+\b", body))

    text_lower = body.lower()
    paywall_hits = _count_distinct_markers(text_lower, PAYWALL_MARKERS)
    fulltext_hits = _count_distinct_markers(text_lower, FULLTEXT_MARKERS)

    # Hard reject: dense paywall markers and not enough fulltext signal.
    # We require 2+ paywall hits AND <3 fulltext hits, so a real article
    # that quotes "purchase" once in its discussion isn't killed.
    if paywall_hits >= 2 and fulltext_hits < 3:
        return False, (
            f"paywall markers={paywall_hits}, fulltext markers={fulltext_hits} "
            f"({char_count} chars)"
        )

    # Cloudflare interstitial is unambiguous on its own.
    if "checking your browser" in text_lower and char_count < 20_000:
        return False, "cloudflare interstitial"
    if "verify you are human" in text_lower and char_count < 20_000:
        return False, "captcha challenge"

    # Length floor. Fulltext markers alone are NOT enough to waive the
    # floor — section headings ("Methods", "Results", "Discussion") show up
    # in nav menus and citation lists even on a stub. Real full text always
    # has at least a few KB of body. We only allow a soft floor (half of
    # min_chars) when the page has BOTH fulltext markers AND substantial
    # body length.
    soft_floor = max(2_000, min_chars // 2)
    if char_count < soft_floor:
        return False, f"body too short: {char_count} chars (floor {soft_floor})"
    if char_count < min_chars and fulltext_hits < 4:
        return (
            False,
            f"body too short: {char_count} chars, {fulltext_hits} fulltext hits",
        )

    if word_count < min_words and fulltext_hits < 4:
        return (
            False,
            f"body too thin: {word_count} words, {fulltext_hits} fulltext hits",
        )

    return True, (
        f"ok: {char_count} chars, {word_count} words, "
        f"fulltext={fulltext_hits}, paywall={paywall_hits}"
    )
