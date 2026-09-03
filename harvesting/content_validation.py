"""Content-quality validation for harvested full text."""

from __future__ import annotations

import re
from urllib.parse import urlparse

JUNK_CONTENT_PATTERNS = [
    "we've detected you are running an older version",
    "browsehappy.com",
    "phosphositeplus",
    "uniprot database entry",
    "sign in to access",
    "subscribe to read",
    "please enable javascript",
    "your browser does not support",
    "access denied",
    "403 forbidden",
    "404 not found",
    "page not found",
    "cookies must be enabled",
    "this site requires javascript",
]

JUNK_CONTENT_DOMAINS = {
    "assays.cancer.gov",
    "antibodies.cancer.gov",
    "biocyc.org",
    "glygen.org",
    "malacards.org",
    "lens.org",
    "clinicaltrials.gov",
    "medlineplus.gov",
}

PAPER_INDICATORS = [
    "abstract",
    "introduction",
    "methods",
    "results",
    "discussion",
    "conclusion",
    "references",
    "materials and methods",
    "patients and methods",
    "study population",
]

VARIANT_PATTERNS = [
    r"[A-Z]\d{2,4}[A-Z]",
    r"p\.[A-Z][a-z]{2}\d+",
    r"c\.\d+[ACGT]>[ACGT]",
    r"mutation",
    r"variant",
]

_ABSTRACT_HEADING_RE = re.compile(r"(?im)^\s{0,3}(?:#{1,6}\s*)?abstract\s*$")
_REFERENCES_HEADING_RE = re.compile(
    r"(?im)^\s{0,3}(?:#{1,6}\s*)?(?:references|bibliography)\s*$"
)
_BODY_HEADING_RE = re.compile(
    r"(?im)^\s{0,3}(?:#{1,6}\s*)?(?:"
    r"introduction|background|methods?|materials(?:\s+and\s+methods)?|"
    r"patients?(?:\s+and\s+methods)?|study\s+population|results?|findings|"
    r"case(?:\s+presentation|\s+report)?|discussion|conclusions?|"
    r"main\s+text|letter|correspondence"
    r")\s*$"
)


def is_abstract_reference_shell(content: str) -> bool:
    """Return True when a purported article is only abstract plus references.

    Large reference lists made byte-size and generic paper-indicator checks
    classify publisher metadata shells as full text.  Require a substantive
    body heading between the Abstract and References headings.  Letters and
    case reports are included in the accepted body vocabulary so this does not
    impose an IMRAD-only article shape.
    """
    if not content:
        return False
    abstract = _ABSTRACT_HEADING_RE.search(content)
    if abstract is None:
        return False
    references = _REFERENCES_HEADING_RE.search(content, abstract.end())
    if references is None:
        return False
    between = content[abstract.end() : references.start()]
    return _BODY_HEADING_RE.search(between) is None


# Folded-supplement evidence. ``harvesting/supplement_fold.py`` appends converted
# supplement files to FULL_CONTEXT under these markers; a paywalled publisher
# page (abstract + navigation + reference list) that carries them is not a
# shell for extraction purposes: the supplement tables are the source.
FOLDED_SUPPLEMENT_MARKERS = (
    "<!-- GVF_FOLDED_SUPPLEMENTS_BEGIN -->",
    "# SUPPLEMENTAL FILE",
)
SUPPLEMENT_CONTENT_MIN_CHARS = 1000
SUPPLEMENT_CONTENT_MIN_TABLE_ROWS = 3
SUPPLEMENT_CONTENT_MIN_VARIANT_TOKENS = 3
_SUPPLEMENT_VARIANT_TOKEN_RE = re.compile(
    r"\b(?:p\.)?[A-Z][a-z]{2}\d{1,4}(?:[A-Z][a-z]{2}|Ter|X|\*|fs)"
    r"|\b[ACDEFGHIKLMNPQRSTVWY]\d{1,4}(?:[ACDEFGHIKLMNPQRSTVWY]|X|\*)\b"
    r"|\bc\.\d+[+-]?\d*(?:_\d+)?(?:[ACGT]>[ACGT]|del|dup|ins)"
)
_PIPE_ROW_RE = re.compile(r"(?m)^\s{0,3}\|.*\|\s*$")


def has_substantive_supplement_content(content: str) -> bool:
    """Return True when folded supplement text large enough to extract from is present.

    RYR2 22222782 is a Springer access shell (abstract, ``Access this article``,
    ``Buy Now``, reference list) followed by two folded ``.doc`` supplements
    holding 235 table rows and every gold variant. ``is_abstract_reference_shell``
    is right that the article body is missing, but refusing the document sent
    extraction a PubMed abstract instead of those tables. 22 of the 70 corpus
    full texts the reuse gate refused on 2026-09-03 were this shape.
    """
    if not content:
        return False
    starts = [i for i in (content.find(m) for m in FOLDED_SUPPLEMENT_MARKERS) if i >= 0]
    if not starts:
        return False
    tail = content[min(starts) :]
    if len(tail) < SUPPLEMENT_CONTENT_MIN_CHARS:
        return False
    # Length alone is fold chrome plus leftover prose; require the fold to
    # carry table rows or variant tokens before it can stand in for a body.
    if len(_PIPE_ROW_RE.findall(tail)) >= SUPPLEMENT_CONTENT_MIN_TABLE_ROWS:
        return True
    return (
        len(_SUPPLEMENT_VARIANT_TOKEN_RE.findall(tail))
        >= SUPPLEMENT_CONTENT_MIN_VARIANT_TOKENS
    )


def _is_binary_content(content: str) -> tuple[bool, str]:
    """Check if content appears to be binary data rather than text.

    PDF conversion failures can produce files with raw binary data.
    This catches those cases before they pollute FULL_CONTEXT.md files.
    """
    if not content:
        return False, ""

    # Check for PDF magic bytes at start
    if content.startswith("%PDF"):
        return True, "Raw PDF binary (starts with %PDF)"

    # Check first 1000 chars for high ratio of non-printable characters
    sample = content[:1000]
    non_printable = sum(1 for c in sample if ord(c) < 32 and c not in "\n\r\t")
    ratio = non_printable / len(sample) if sample else 0

    if ratio > 0.1:  # More than 10% non-printable = likely binary
        return True, f"Binary content detected ({ratio:.1%} non-printable chars)"

    # Check for common binary patterns (null bytes, control chars)
    if "\x00" in sample:
        return True, "Contains null bytes (binary data)"

    # Check for excessive high-byte chars (common in corrupted PDF extraction)
    high_bytes = sum(1 for c in sample if ord(c) > 127)
    high_ratio = high_bytes / len(sample) if sample else 0

    if high_ratio > 0.3:  # More than 30% high-byte chars
        return (
            True,
            f"Likely binary/corrupted content ({high_ratio:.1%} high-byte chars)",
        )

    return False, ""


def validate_content_quality(
    content: str, source_url: str | None = None
) -> tuple[bool, str]:
    """Validate harvested content is likely a real paper."""
    if not content:
        return False, "Empty content"

    # Check for binary content first (PDF conversion failures)
    is_binary, binary_reason = _is_binary_content(content)
    if is_binary:
        return False, binary_reason

    content_lower = content.lower()

    if is_abstract_reference_shell(content) and not has_substantive_supplement_content(
        content
    ):
        return False, "Abstract/reference shell: no substantive article body"

    for pattern in JUNK_CONTENT_PATTERNS:
        if pattern in content_lower:
            return False, f"Junk content detected: '{pattern}'"

    if source_url:
        try:
            domain = urlparse(source_url).netloc.lower()
            for junk_domain in JUNK_CONTENT_DOMAINS:
                if junk_domain in domain:
                    return False, f"Content from non-article domain: {junk_domain}"
        except Exception:
            pass

    if len(content) < 1500:
        return False, f"Content too short ({len(content)} chars)"

    indicator_count = sum(1 for ind in PAPER_INDICATORS if ind in content_lower)
    if indicator_count < 2:
        variant_matches = sum(
            1 for pat in VARIANT_PATTERNS if re.search(pat, content, re.IGNORECASE)
        )
        if variant_matches < 2:
            return (
                False,
                f"Missing paper structure (only {indicator_count} indicators) and no variant content",
            )

    return True, "Valid content"
