"""Shared authenticated HTTP session helpers for paywall recovery.

This module lives in a shipped package so both the repository recovery CLI and
installed ``gvf`` preflight use identical EZproxy and browser-cookie behavior.
"""

from __future__ import annotations

import os
from urllib.parse import urlparse

import requests


def _ncbi_email() -> str:
    email = os.environ.get("ENTREZ_EMAIL") or os.environ.get("NCBI_EMAIL")
    if not email:
        raise RuntimeError("Set ENTREZ_EMAIL or NCBI_EMAIL before querying NCBI.")
    return email


def make_session() -> requests.Session:
    """Build the requests session used for publisher and EZproxy traffic."""
    session = requests.Session()
    session.headers.update(
        {"User-Agent": f"GVF-PaywalledFetch/1.0 (mailto:{_ncbi_email()})"}
    )

    from harvesting.browser_html import ezproxy

    ezproxy.install_on_session(session)
    return session


def hydrate_session_with_browser_cookies(
    session: requests.Session, cookies: list[dict]
) -> int:
    """Copy Playwright-format browser cookies into a requests session."""
    added = 0
    for cookie in cookies or []:
        name = cookie.get("name")
        value = cookie.get("value")
        if not name or value is None:
            continue

        domain = cookie.get("domain") or ""
        if not domain and cookie.get("url"):
            domain = urlparse(str(cookie["url"])).hostname or ""
        if not domain:
            continue

        session.cookies.set(
            str(name),
            str(value),
            domain=str(domain),
            path=str(cookie.get("path") or "/"),
        )
        added += 1
    return added
