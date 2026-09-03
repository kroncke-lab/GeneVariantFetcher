"""The one canonical normalization layer for harvested source text.

Every extraction route used to clean text its own way. Measured divergence on
cells that appear verbatim in converted publisher sources:

===========================  ==================  ======================
raw cell                     ``table_router``    ``variant_normalizer``
===========================  ==================  ======================
``c.1471C\\u2009>T``           ``None``            ``c.1471C>T``
``c.722\\u00a0T>C``            ``None``            ``c.722T>C``
``c.2550\\u20112551insTG``     ``None``            ``c.2550_2551insTG``
``p.Arg491Trp\\u200b``         ``None``            ``ARG491TRP\\u200b``
``c.354T\\uff1eG``             ``None``            ``c.354T\\uff1eG``
===========================  ==================  ======================

So one allele was a parsed identity on one route, a dropped row on another and
a junk identity on a third. This module is the single place that decision is
made; the per-route cleaners delegate here.

Deliberately **not** ``unicodedata.normalize("NFKC")``
------------------------------------------------------
NFKC is a semantic rewrite. It folds primes, superscripts, subscripts, Greek
letters and the micro sign, any of which can be meaning-bearing in a residue
label or an OCR'd figure dump. This module instead enumerates a closed set of
*presentational* substitutions, so the mapping is auditable, testable and
reversible in review. Characters outside that set are passed through unchanged.

Two rules exist because they are the ones that silently corrupt identities:

**Zero-width characters become a space, never nothing.** A publisher PDF that
extracts as ``c.1471C>T<ZWSP>p.Arg491Trp`` holds *two* alleles. Deleting the
zero-width joiner fuses them into one chimeric token that parses as neither and
can be persisted as a fabricated identity. Substituting a space splits them.
The cost is that a zero-width character *inside* one allele turns that allele
into two fragments and it is missed. A miss is recoverable; a fabricated
identity in a curated database is not.

**Dashes become ASCII hyphen-minus, never underscore.** In HGVS a hyphen is an
intronic acceptor offset (``c.123-2A>G``), a 5' UTR coordinate (``c.-15C>T``)
or a range separator (``c.2550-2551insTG``), and only context tells them apart.
Collapsing the Unicode dash family to one ASCII character is presentation.
Deciding that a given hyphen means a range is notation semantics, and stays
where it already lives, behind the strict adjacency guard in
``utils.variant_normalizer._repair_hyphen_range``.
"""

from __future__ import annotations

import re
from typing import Optional

__all__ = [
    "normalize_source_text",
    "normalize_notation_token",
    "normalize_header_text",
    "split_glued_notation",
    "NULL_CELL_VALUES",
]


#: Cell payloads that mean "no value". Compared after normalization/casefold.
NULL_CELL_VALUES = frozenset(
    {"", "-", "--", ".", "na", "n/a", "nan", "none", "nd", "n.d.", "null"}
)

# --- Closed substitution tables (presentation only) -------------------------

#: Space-like characters publishers emit instead of U+0020. Includes the
#: NBSP/thin/hair spaces seen in real converted tables plus the rest of the
#: Unicode space family, so the set cannot drift as new sources arrive.
_SPACE_LIKE = (
    " "  # no-break space
    " "  # ogham space mark
    "           "
    " "  # narrow no-break space
    " "  # medium mathematical space
    "　"  # ideographic space
)

#: Zero-width characters. Mapped to a space, never removed (see module docs).
_ZERO_WIDTH = (
    "​"  # zero width space
    "‌"  # zero width non-joiner
    "‍"  # zero width joiner
    "⁠"  # word joiner
    "﻿"  # zero width no-break space / BOM
)

#: Dash family mapped to ASCII hyphen-minus. Range-versus-offset semantics are
#: decided downstream, never here.
_DASH_LIKE = (
    "‐"  # hyphen
    "‑"  # non-breaking hyphen
    "‒"  # figure dash
    "–"  # en dash
    "—"  # em dash
    "―"  # horizontal bar
    "−"  # minus sign
    "˗"  # modifier letter minus sign
    "﹣"  # small hyphen-minus
    "－"  # fullwidth hyphen-minus
)

#: Arrow glyphs used for substitutions in mutation tables ("1682C→T").
_ARROW_LIKE = "→⇒⟶➔➜➡"

#: Characters deliberately NOT folded, with the reason. Referenced by tests so
#: the exclusion is a pinned decision rather than an oversight.
PRESERVED_CHARACTERS = {
    "′": "prime: can be a residue/chain marker, not punctuation",
    "″": "double prime: same",
    "µ": "micro sign: unit, meaning-bearing in assay text",
    "μ": "greek mu: unit",
    "×": "multiplication sign: magnification/fold-change, not 'x'",
    "°": "degree sign: unit",
    "⁰": "superscript zero: exponent, not a digit in-line",
    "²": "superscript two: exponent",
}

_TRANSLATION: dict[int, str] = {}
for _ch in _SPACE_LIKE:
    _TRANSLATION[ord(_ch)] = " "
for _ch in _ZERO_WIDTH:
    _TRANSLATION[ord(_ch)] = " "
for _ch in _DASH_LIKE:
    _TRANSLATION[ord(_ch)] = "-"
for _ch in _ARROW_LIKE:
    _TRANSLATION[ord(_ch)] = ">"
# Fullwidth ASCII block U+FF01..U+FF5E maps 1:1 onto U+0021..U+007E. This is
# what turns "c.354T＞G" into "c.354T>G". Fullwidth digits and letters fold
# with it, which is intended: a fullwidth digit in a coordinate is the same
# coordinate.
for _code in range(0xFF01, 0xFF5F):
    _TRANSLATION[_code] = chr(_code - 0xFEE0)
# The fullwidth block contains a hyphen at U+FF0D; the dash rule above already
# claimed it and must win, so re-assert after the block loop.
_TRANSLATION[0xFF0D] = "-"
del _ch, _code

#: Footnote/marker glyphs stripped only from the END of an isolated notation
#: token ("c.1149insT*"), never from inside it and never from prose.
_TRAILING_MARKERS = ("*", "†", "‡", "§", "¶", "#", "^")

#: A notation prefix. Used to refuse gluing two alleles into one token.
_NOTATION_PREFIX_RE = re.compile(r"(?:\b|^)[cpgmnr]\.", re.IGNORECASE)

#: A standalone allele-shaped token that carries NO prefix: a three-letter or
#: one-letter residue, a position, and a consequence ("Arg491Trp", "Q403*",
#: "R147Ter", "Gly262Alafs*98"). Prefix-based splitting alone missed these, so
#: a zero-width joiner between two bare residue labels still fused them.
_BARE_NOTATION_RE = re.compile(
    r"^(?:[A-Z][a-z]{2}|[A-Z])\d{1,5}"
    r"(?:[A-Z][a-z]{2}|[A-Z]|\*|X|Ter|=|fs\w*|del\w*|dup\w*|ins\w*)",
    re.IGNORECASE,
)


def _looks_like_notation(part: str) -> bool:
    return bool(_NOTATION_PREFIX_RE.match(part) or _BARE_NOTATION_RE.match(part))


def normalize_source_text(text: Optional[str]) -> str:
    """Canonicalize harvested text while preserving token boundaries.

    This is the document- and cell-level entry point: the form fed to the
    regex scanner, the table parsers, and LLM prompts. Line structure and word
    boundaries survive, so a downstream tokenizer sees the same words it always
    did, just spelled in ASCII where the difference was presentational.
    """
    if not text:
        return ""
    normalized = str(text).translate(_TRANSLATION)
    # Normalize line endings so line-oriented parsers agree across converters.
    normalized = normalized.replace("\r\n", "\n").replace("\r", "\n")
    return normalized


def normalize_header_text(value: Optional[str]) -> str:
    """Canonical key for deterministic table-header matching.

    Applies the shared character policy first, so a header that differs only by
    a no-break space collapses onto the same key as its ASCII twin.
    """
    return re.sub(r"[^a-z0-9]+", "", normalize_source_text(value).lower())


def split_glued_notation(value: Optional[str]) -> list[str]:
    """Split a cell that holds several whitespace-separated notation tokens.

    ``"c.1471C>T p.Arg491Trp"`` is two identities. Returning them separately is
    what lets :func:`normalize_notation_token` refuse to fuse them while the
    caller can still parse each half.
    """
    normalized = normalize_source_text(value).strip()
    if not normalized:
        return []
    parts = [part for part in re.split(r"\s+", normalized) if part]
    if len(parts) < 2:
        return parts
    # Only treat whitespace as a separator when more than one part is itself
    # allele-shaped; otherwise "c.1471C >T" is one token wrapped by the
    # converter. Bare residue labels count: two unprefixed protein calls glued
    # by a zero-width joiner ("Arg491Trp<ZWSP>Gln403Ter") are still two
    # alleles, and matching only on a "c."/"p." prefix let those fuse.
    if sum(1 for part in parts if _looks_like_notation(part)) < 2:
        return [normalized]
    return parts


def normalize_notation_token(
    value: Optional[str], *, strip_trailing_markers: bool = True
) -> Optional[str]:
    """Canonicalize one isolated variant-notation token.

    Beyond the shared character policy this removes internal ASCII whitespace,
    because a table cell holding one allele may have been wrapped mid-token by
    the converter (``"c.1471C >T"``, ``"c.722\\u00a0T>C"``).

    Returns ``None`` for an empty or null-marker cell, and — the guard that
    matters — returns ``None`` rather than a fused token when the cell holds
    two whitespace-separated notation prefixes. Gluing those would mint an
    identity that appears in no source.
    """
    if value is None:
        return None
    normalized = normalize_source_text(value).strip()
    if not normalized:
        return None
    if normalized.casefold() in NULL_CELL_VALUES:
        return None

    parts = split_glued_notation(normalized)
    if len(parts) > 1:
        # Two alleles in one cell: the caller must split first. Never fuse.
        return None

    collapsed = re.sub(r"\s+", "", normalized)
    if strip_trailing_markers:
        for marker in _TRAILING_MARKERS:
            if collapsed.endswith(marker) and len(collapsed) > len(marker):
                collapsed = collapsed[: -len(marker)]
                break
    collapsed = collapsed.rstrip(",;")
    if not collapsed or collapsed.casefold() in NULL_CELL_VALUES:
        return None
    return collapsed
