"""Load cookies from a local Chrome profile and convert to Playwright format.

Used by AuthenticatedBrowserPool so headless Playwright inherits the user's
publisher SSO sessions (e.g. VUMC EZproxy / institutional access).

Reads Chrome's cookie store via `browser_cookie3` — no need to attach to a
running Chrome process or copy the user-data-dir, which avoids profile-lock
conflicts when Chrome is open.
"""

from __future__ import annotations

import logging
from typing import Iterable, List, Optional

logger = logging.getLogger(__name__)


# Publisher domains we care about for paywalled-paper recovery.
# Subdomains are matched by `browser_cookie3` (a parent-domain query returns
# cookies set on any subdomain), so listing the bare registrable domain is
# enough.
DEFAULT_PUBLISHER_DOMAINS: tuple = (
    # AHA Journals (Circulation, JAHA, Stroke, Hypertension, Circ Res, ...)
    "ahajournals.org",
    # Elsevier (Heart Rhythm, JACC, Cell, Lancet, ...) — content is on
    # sciencedirect.com plus a CDN at sciencedirectassets.com.
    "sciencedirect.com",
    "sciencedirectassets.com",
    "elsevier.com",
    "heartrhythmjournal.com",
    "jacc.org",
    "cell.com",
    "thelancet.com",
    "mayoclinicproceedings.org",
    # Wiley (Human Mutation, etc.)
    "onlinelibrary.wiley.com",
    "wiley.com",
    # Oxford Academic (Europace, EHJ, HMG)
    "academic.oup.com",
    "oup.com",
    # Springer / BMC / Nature
    "springer.com",
    "link.springer.com",
    "biomedcentral.com",
    "nature.com",
    # Karger
    "karger.com",
    # Sage / Mary Ann Liebert (DOI 10.1089 hosts on journals.sagepub.com
    # despite the Liebert imprint — Sage acquired Liebert content delivery).
    "sagepub.com",
    "journals.sagepub.com",
    "liebertpub.com",
    # BMJ
    "bmj.com",
    # Wolters Kluwer / Neurology journals (e.g. DOI 10.1212 supplements)
    "neurology.org",
    # Lippincott Williams & Wilkins / Wolters Kluwer
    "lww.com",
    # PMC (cookies aren't needed but harmless)
    "ncbi.nlm.nih.gov",
    "pmc.ncbi.nlm.nih.gov",
    # Institutional SSO / proxy — needed for the auth handshake redirect chain.
    "vumc.org",
    "vanderbilt.edu",
    "ezproxy.library.vanderbilt.edu",
    "microsoftonline.com",
    "login.microsoft.com",
    "okta.com",
    "shibboleth.net",
)


def _coerce_same_site(rest: dict) -> Optional[str]:
    """Best-effort SameSite extraction from cookiejar `rest` attrs."""
    raw = (rest or {}).get("SameSite") or (rest or {}).get("samesite")
    if not raw:
        return None
    raw_l = str(raw).strip().lower()
    if raw_l in ("strict",):
        return "Strict"
    if raw_l in ("lax",):
        return "Lax"
    if raw_l in ("none", "no_restriction", "unspecified"):
        return "None"
    return None


def chrome_cookie_to_playwright(cookie) -> Optional[dict]:
    """Convert a `http.cookiejar.Cookie` to the dict shape Playwright expects.

    Returns None if the cookie is unusable (no name/value or expired session
    cookie with no expiry that Playwright would reject).
    """
    if not cookie.name:
        return None

    # Domain: cookiejar uses a leading dot when the cookie applies to
    # subdomains; Playwright accepts either, but a leading dot is the
    # safer signal of "all subdomains".
    domain = cookie.domain or ""
    if not domain:
        return None

    out: dict = {
        "name": cookie.name,
        "value": cookie.value or "",
        "domain": domain,
        "path": cookie.path or "/",
        "secure": bool(cookie.secure),
        # cookiejar stores HttpOnly in `_rest` / `rest` under "HttpOnly".
        "httpOnly": bool(
            (getattr(cookie, "_rest", None) or {}).get("HttpOnly")
            or (getattr(cookie, "_rest", None) or {}).get("httponly")
        ),
    }

    # Playwright requires expires to be a positive int (Unix seconds) or -1
    # for session cookies. cookiejar uses None for session cookies.
    if cookie.expires:
        try:
            out["expires"] = int(cookie.expires)
        except (TypeError, ValueError):
            out["expires"] = -1
    else:
        out["expires"] = -1

    same_site = _coerce_same_site(getattr(cookie, "_rest", None))
    if same_site:
        out["sameSite"] = same_site

    return out


def load_chrome_cookies(
    domains: Optional[Iterable[str]] = None,
    profile_name: Optional[str] = None,
) -> List[dict]:
    """Load Chrome cookies for the given domains, in Playwright format.

    Args:
        domains: registrable domains to fetch cookies for. Defaults to the
            publisher + SSO list we maintain in this module.
        profile_name: optional Chrome profile folder name (e.g. ``"Default"``,
            ``"Profile 1"``). When omitted, browser_cookie3 walks all profiles
            and merges results — that's usually what we want, since a user may
            be signed into the publisher in any profile.

    Returns:
        List of cookie dicts ready for ``BrowserContext.add_cookies()``.
        Failures to read a given domain are logged and skipped (we don't want
        one bad domain to disable the whole fetcher).
    """
    try:
        import browser_cookie3 as bc
    except ImportError as e:
        raise RuntimeError(
            "browser_cookie3 is required for authenticated fetching. "
            "Install with `pip install browser-cookie3`."
        ) from e

    domains_iter = (
        list(domains) if domains is not None else list(DEFAULT_PUBLISHER_DOMAINS)
    )
    seen_keys: set = set()
    out: List[dict] = []

    for d in domains_iter:
        try:
            if profile_name:
                jar = bc.chrome(domain_name=d, profile=profile_name)
            else:
                jar = bc.chrome(domain_name=d)
        except Exception as e:
            logger.warning("cookie load failed for %s: %s", d, e)
            continue

        count_before = len(out)
        for c in jar:
            pw = chrome_cookie_to_playwright(c)
            if pw is None:
                continue
            # Dedup across overlapping domain queries (e.g. wiley.com vs
            # onlinelibrary.wiley.com both pulling the same cookie).
            key = (pw["name"], pw["domain"], pw["path"])
            if key in seen_keys:
                continue
            seen_keys.add(key)
            out.append(pw)
        logger.debug("cookie load %s: +%d cookies", d, len(out) - count_before)

    logger.info(
        "Loaded %d Chrome cookies across %d domains", len(out), len(domains_iter)
    )
    return out


def cookie_domain_summary(cookies: List[dict]) -> dict:
    """Count cookies per registrable-ish domain. For diagnostics only."""
    counts: dict = {}
    for c in cookies:
        d = (c.get("domain") or "").lstrip(".")
        counts[d] = counts.get(d, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))
