"""Helpers for PMC's lightweight proof-of-work supplement gate."""

from __future__ import annotations

import hashlib
import re
from typing import Any
from urllib.parse import urlparse


PMC_POW_COOKIE = "cloudpmc-viewer-pow"
_CHALLENGE_RE = re.compile(r'POW_CHALLENGE\s*=\s*"([^"]+)"')
_DIFFICULTY_RE = re.compile(r'POW_DIFFICULTY\s*=\s*"(\d+)"')


def is_pmc_pow_challenge(content: bytes) -> bool:
    """Return True when a PMC supplement response is the download POW page."""
    sample = content[:4096].lower()
    return (
        b"preparing to download" in sample
        and b"cloudpmc-viewer-pow" in sample
        and b"pow_challenge" in sample
    )


def solve_pmc_pow_cookie(
    html: str, *, max_nonce: int = 2_000_000, max_difficulty: int = 5
) -> str | None:
    """Solve PMC's SHA256 nonce challenge and return the cookie value."""
    challenge_match = _CHALLENGE_RE.search(html)
    if not challenge_match:
        return None
    difficulty_match = _DIFFICULTY_RE.search(html)
    difficulty = int(difficulty_match.group(1)) if difficulty_match else 4
    if difficulty > max_difficulty:
        return None
    prefix = "0" * difficulty
    challenge = challenge_match.group(1)

    for nonce in range(max_nonce + 1):
        digest = hashlib.sha256(f"{challenge}{nonce}".encode("utf-8")).hexdigest()
        if digest.startswith(prefix):
            return f"{challenge},{nonce}"
    return None


def attach_pmc_pow_cookie(session: Any, *, html: str, url: str) -> bool:
    """Solve a PMC POW page and attach its cookie to the requests session."""
    cookie_value = solve_pmc_pow_cookie(html)
    if not cookie_value:
        return False

    hostname = urlparse(url).hostname or "pmc.ncbi.nlm.nih.gov"
    session.cookies.set(PMC_POW_COOKIE, cookie_value, domain=hostname, path="/")
    return True
