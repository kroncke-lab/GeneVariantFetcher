"""The one canonical DOI extraction and cleaning path.

Four modules carried their own copy of this logic (`harvesting/orchestrator`,
`harvesting/supplement_processing_service`, and two under
`scripts/recall_audit/`), each with the same greedy pattern::

    10\\.\\d{4,9}/[^\\s"'<>]+

That pattern has no word boundary, which is fine for a DOI sitting in a URL or
a metadata field and wrong for one harvested from *rendered page text*. A real
publisher landing page renders::

    https://doi.org/10.3390/ijms25052734Submission received: 6 December 2023

with no delimiter between the DOI and the next word, so the extracted DOI
became ``10.3390/ijms25052734Submission``. Downstream that is not a resolution
failure with a clear message; it is a fetch attempt against an identifier that
cannot exist, and it failed a production run's source-recovery stage.

The repair is a boundary rule, not a publisher rule: a trailing run of
capitalised words that follows a digit is prose that got glued on, because a
DOI suffix does not end that way. ``Submission``, ``Received``, ``Published``
and ``Abstract`` are all handled by the same rule without an allowlist — this
codebase has already recorded that ignore-lists of this shape do not converge.

Trimming is conservative in the direction that matters. Over-trimming yields a
shorter DOI that still resolves or fails loudly; under-trimming yields a
silently unresolvable identifier.
"""

from __future__ import annotations

import re
from typing import Any, Optional
from urllib.parse import unquote

__all__ = ["DOI_RE", "clean_doi", "extract_doi", "doi_from_text"]

#: A DOI token. Deliberately unchanged from the historical pattern so this is a
#: drop-in for every call site; the boundary repair happens in `clean_doi`.
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\"'<>]+", re.IGNORECASE)

_DOI_URL_PREFIX_RE = re.compile(
    r"^(?:https?://)?(?:dx\.)?doi\.org/", flags=re.IGNORECASE
)

#: One or more capitalised words glued directly onto a digit at the end of the
#: token. Requires the uppercase letter to be followed by at least two
#: lowercase letters, so a legitimate uppercase suffix segment (``.../CIRCEP``,
#: ``.../S1``) is never touched.
_GLUED_TRAILING_WORDS_RE = re.compile(r"(?<=\d)(?:[A-Z][a-z]{2,})+$")

#: Trailing punctuation that ordinary prose leaves attached.
_TRAILING_PUNCTUATION = ".,;:)]}>"


def clean_doi(value: Any) -> str:
    """Normalize one DOI string: strip a resolver prefix, punctuation and glue."""
    doi = unquote(str(value or "").strip())
    doi = _DOI_URL_PREFIX_RE.sub("", doi)
    doi = doi.strip().strip(_TRAILING_PUNCTUATION)
    # Prose glued onto the end by a page render with no delimiter.
    repaired = _GLUED_TRAILING_WORDS_RE.sub("", doi)
    # Never trim away the whole suffix: a DOI must keep something after "/".
    if repaired and repaired.rstrip("/").count("/") >= 1:
        prefix, _, suffix = repaired.partition("/")
        if suffix:
            doi = repaired
    return doi.strip().strip(_TRAILING_PUNCTUATION)


def extract_doi(value: Any) -> str:
    """Return the first DOI in ``value``, cleaned. Empty string when absent."""
    if value is None:
        return ""
    match = DOI_RE.search(unquote(str(value)))
    return clean_doi(match.group(0)) if match else ""


def doi_from_text(text: Optional[str]) -> str:
    """Pull the article's own DOI out of free text.

    Prefers a ``doi:``-labelled match — the NLM citation line names the article
    itself — before falling back to the first DOI anywhere, so a DOI belonging
    to a cited reference is not mistaken for the paper's.
    """
    if not text:
        return ""
    labelled = re.search(
        r"doi:\s*(10\.\d{4,9}/[^\s\"'<>]+)", str(text), flags=re.IGNORECASE
    )
    if labelled:
        return clean_doi(labelled.group(1))
    return extract_doi(text)
