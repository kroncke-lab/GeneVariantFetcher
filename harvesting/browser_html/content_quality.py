"""Content quality gate for Tier 3.5 results.

Authenticated fetches can return things that *look* like a full article but
are actually paywall stubs, "access denied" pages, or institutional-login
landing pages. The pipeline cannot tell the difference downstream, so we
reject obvious failures here.

Hardened after a 2026-05-11 audit of 9 "successful" paywall recoveries
showed that 5 of 9 were abstract-only stubs that the prior gate let
through. The fix layers four hard rejects on top of the existing length
heuristics:

1. **Explicit paywall sentences** — exact phrases that only appear on
   paywall stubs (e.g. "Available to Purchase", "You do not currently have
   access to this content", "Log in, subscribe or purchase for full
   access"). One hit ⇒ reject.
2. **Body content under section headers** — fulltext markers like
   "Methods"/"Results"/"Discussion" appear as nav menu items and structured-
   abstract subsection headers on stubs too. We require the *combined*
   character count of paragraph text under top-level Methods/Results/
   Discussion/Conclusion headers to clear a threshold; abstract-only stubs
   max out around 1.5–2 KB in those sections, real bodies clear 3 KB
   trivially.
3. **Repeated-line padding** — Wiley/Elsevier stubs that pad to ~20 KB
   by repeating the same author affiliation line 5–10 times.
4. **Variant-token count** — surfaced as a soft signal in the reason
   string; useful for callers (paywall recovery) that expect variants.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Dict, Iterable, List, Tuple

# Phrases that, in clusters, indicate the page is a paywall/access stub.
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

# Hard-reject phrases. Each is a *full-sentence* phrase that only appears
# on publisher paywall stubs — not in real biomedical article body or in
# the publisher's nav/footer chrome. Generic short phrases like "purchase
# this article" are intentionally NOT here because they appear verbatim in
# the access-options footer/sidebar of real AHA/Elsevier articles. The
# phrases below were verified against 5 confirmed stubs and 4 real papers
# from the 2026-05-11 audit; each hits all stubs and no real papers.
HARD_REJECT_PHRASES = (
    "available to purchase",  # Karger stub title
    "you do not currently have access to this content",  # Karger stub
    "log in, subscribe or purchase for full access",  # Elsevier stub
    "get full access to this article",  # AAN stub
    "view all available purchase options and get full access",  # AAN variant
    "checking your browser before accessing",  # Cloudflare hard
    "this content requires a subscription",  # generic stub
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

# Section names treated as "real body" indicators. We look for these as
# top-level (## or ###) markdown headers and require accumulated paragraph
# text under them to clear a threshold.
BODY_SECTION_NAMES = (
    "methods",
    "materials and methods",
    "patients and methods",
    "study population",
    "subjects",
    "results",
    "findings",
    "discussion",
    "conclusion",
    "conclusions",
)

# Pattern for variant-like tokens. Covers single-letter notation (A123V,
# R273X, K897T) and HGVS p./c. notation. 3-letter notation (Ala123Val) is
# handled by a separate regex below.
_SHORT_VARIANT_RE = re.compile(r"\b[A-Z][0-9]{2,4}[A-Z*X]\b")
_HGVS_RE = re.compile(
    r"\bp\.[A-Za-z][a-z]{2}\d+[A-Za-z][a-z]{2}\b|\bc\.\d+[ACGT]>[ACGT]\b"
)
_THREE_LETTER_RE = re.compile(
    r"\b(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|"
    r"Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)\d+"
    r"(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|"
    r"Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter|Stop|fs|\*)\b"
)


def _count_distinct_markers(text_lower: str, markers: Tuple[str, ...]) -> int:
    return sum(1 for m in markers if m in text_lower)


def _find_hard_reject(text_lower: str) -> str:
    """Return the first hard-reject phrase found, or empty string."""
    for phrase in HARD_REJECT_PHRASES:
        if phrase in text_lower:
            return phrase
    return ""


def _body_content_chars(text: str, min_para_line_chars: int = 80) -> int:
    """Count chars in substantive body content: paragraphs OR table cells.

    Real biomedical papers have body content distributed across:
    - Long narrative paragraphs (Methods/Results/Discussion text)
    - Markdown tables (cohort papers put their variant data here — Moss
      2002 has a 14-variant pore table; Splawski 2000 has 60+ row variant
      tables; Horigome 2010 has a 27-case cohort table).

    Stubs have neither — only short affiliation/nav lines and abstract.

    A real cohort paper might have only a few KB of narrative but 10+ KB
    of table content. Counting both keeps cohort papers from being killed
    while still flagging abstract-only stubs (which have no tables).
    """
    total = 0
    for line in text.split("\n"):
        s = line.strip()
        if not s or s.startswith("#") or s.startswith(">"):
            continue
        if s.startswith("- ") or s.startswith("* ") or s.startswith("+ "):
            continue
        # Table rows: count cell content if at least 2 non-empty, non-separator
        # cells. ``| --- | --- |`` separator rows and single-cell rows don't
        # count as body — those exist in stubs too.
        if s.startswith("|"):
            cells = [c.strip() for c in s.strip("|").split("|")]
            substantive = [
                c for c in cells if c and c != "---" and not re.fullmatch(r"-{3,}", c)
            ]
            if len(substantive) >= 2:
                total += sum(len(c) for c in substantive)
            continue
        if len(s) < min_para_line_chars:
            continue
        total += len(s)
    return total


# Kept for back-compat with any external caller that imported the old name.


def _detect_repeated_line_padding(
    text: str, min_line_chars: int = 20, min_repeats: int = 5
) -> Tuple[bool, int, str]:
    """Detect the affiliation-padding pattern (same line repeated 5+ times).

    Returns ``(triggered, count, line)``. A real biomedical article rarely
    repeats any single >20-char line more than 3–4 times; affiliation-padded
    stubs hit 10+ easily.

    Markdown furniture (table separator rows like ``| --- | --- |``, list
    markers, blank-cell rows) is excluded because cohort papers with wide
    tables have many separator-row repeats by construction.
    """
    counts: Counter = Counter()
    for line in text.split("\n"):
        s = line.strip()
        if len(s) < min_line_chars:
            continue
        # Skip markdown table separator rows (only |, -, : and whitespace).
        if re.fullmatch(r"[\s\|:\-]+", s):
            continue
        # Skip rows that are mostly empty cells (table-furniture).
        if s.startswith("|"):
            cells = [c.strip() for c in s.strip("|").split("|")]
            non_empty = [c for c in cells if c and c not in ("---", "-")]
            if len(non_empty) < 2:
                continue
        counts[s] += 1
    if not counts:
        return False, 0, ""
    line, n = counts.most_common(1)[0]
    if n >= min_repeats:
        return True, n, line
    return False, n, line


def _variant_token_count(text: str) -> int:
    """Count distinct variant-like tokens across the standard regexes."""
    seen: set = set()
    seen.update(_SHORT_VARIANT_RE.findall(text))
    seen.update(_HGVS_RE.findall(text))
    seen.update(_THREE_LETTER_RE.findall(text))
    return len(seen)


def validate_article_content(
    text: str,
    min_chars: int = 5_000,
    min_words: int = 800,
    min_body_content_chars: int = 3_000,
) -> Tuple[bool, str]:
    """Return ``(ok, reason)``. ``reason`` describes why we accepted or rejected.

    Designed to be pluggable into
    ``BrowserHTMLFetcher(validate_content_quality=...)``.

    Args:
        text: Extracted markdown body.
        min_chars: Hard char-count floor for total body.
        min_words: Soft word-count floor (waived when many fulltext markers).
        min_body_content_chars: Floor for chars in body paragraphs
            (lines ≥80 chars that aren't headers/tables/lists). Abstract-
            only stubs cap at ~1500 chars here; real papers clear 10 000+.
    """
    if not text:
        return False, "empty body"

    body = text.strip()
    char_count = len(body)
    word_count = len(re.findall(r"\b\w+\b", body))
    text_lower = body.lower()

    # FIX #1 — explicit paywall sentences. Hard reject on the first match.
    hard = _find_hard_reject(text_lower)
    if hard:
        return False, f"explicit paywall phrase: '{hard}' ({char_count} chars)"

    # Cloudflare interstitial / captcha are unambiguous on their own.
    if "checking your browser" in text_lower and char_count < 20_000:
        return False, "cloudflare interstitial"
    if "verify you are human" in text_lower and char_count < 20_000:
        return False, "captcha challenge"

    # FIX #2 — require substantive body paragraph content, not just section
    # header names in nav or affiliation boilerplate. ``_body_content_chars``
    # sums chars in lines ≥80 chars that aren't headers/tables/lists. Real
    # bodies clear this trivially; stubs cap below 2 KB.
    body_para_chars = _body_content_chars(body)
    paywall_hits = _count_distinct_markers(text_lower, PAYWALL_MARKERS)
    fulltext_hits = _count_distinct_markers(text_lower, FULLTEXT_MARKERS)
    variant_count = _variant_token_count(body)

    if body_para_chars < min_body_content_chars:
        return False, (
            f"insufficient body paragraph content: {body_para_chars} chars "
            f"(floor {min_body_content_chars}); "
            f"paywall={paywall_hits} fulltext={fulltext_hits} "
            f"variants={variant_count} total={char_count}"
        )

    # FIX #3 — affiliation/boilerplate padding. Only triggers when the
    # body-paragraph floor is also marginal. BMC papers repeat the
    # "Article CAS PubMed Google Scholar" reference-list furniture 40+
    # times but have 50+ KB of real body — those shouldn't be killed.
    if body_para_chars < 10_000:
        padded, n_repeats, padded_line = _detect_repeated_line_padding(body)
        if padded:
            snippet = (padded_line or "")[:60]
            return False, (
                f"repeated-line padding: '{snippet}…' x{n_repeats} "
                f"with thin body ({body_para_chars} chars)"
            )

    # Legacy dense-paywall-marker check.
    if paywall_hits >= 2 and fulltext_hits < 3:
        return False, (
            f"paywall markers={paywall_hits}, fulltext markers={fulltext_hits} "
            f"({char_count} chars)"
        )

    # Length floor (kept as defense in depth — real papers always clear
    # the soft floor easily).
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
        f"body_para={body_para_chars}, variants={variant_count}, "
        f"fulltext={fulltext_hits}, paywall={paywall_hits}"
    )
