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


def validate_content_quality(content: str, source_url: str | None = None) -> tuple[bool, str]:
    """Validate harvested content is likely a real paper."""
    if not content:
        return False, "Empty content"

    content_lower = content.lower()

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
