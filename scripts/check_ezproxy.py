#!/usr/bin/env python3
"""Verify institutional EZproxy routing clears the Cloudflare wall (Wiley etc.).

Mission-critical sanity check: confirms that, with ``GVF_EZPROXY_PREFIX``/``HOST``
set and a valid EZproxy session cookie loaded, a Wiley Online Library request
returns licensed content instead of the Cloudflare 403 "Just a moment…" page.

  GVF_EZPROXY_PREFIX="https://ezproxy.library.vanderbilt.edu/login?url=" \\
      python scripts/check_ezproxy.py --doi 10.1111/jce.14865

Reports, for the direct and the EZproxy-routed request: HTTP status, whether a
Cloudflare challenge was detected, and the response size. A PASS is a non-403
response with no CF challenge and a real body.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from harvesting.browser_html import ezproxy  # noqa: E402

_CF_MARKERS = ("just a moment", "cf-chl", "challenge-platform", "/cdn-cgi/challenge")


def _looks_cf(text: str, status: int, headers) -> bool:
    low = (text or "")[:4000].lower()
    server = str(headers.get("server", "")).lower() if headers else ""
    return (
        status == 403
        or "cloudflare" in server
        and any(m in low for m in _CF_MARKERS)
        or any(m in low for m in _CF_MARKERS)
    )


def _probe(session, url: str, label: str) -> bool:
    try:
        r = session.get(url, timeout=45, allow_redirects=True)
    except Exception as exc:  # noqa: BLE001
        print(f"  {label:18s}: ERROR {exc}")
        return False
    cf = _looks_cf(r.text, r.status_code, r.headers)
    ok = (r.status_code == 200) and not cf and len(r.text) > 4000
    print(
        f"  {label:18s}: status={r.status_code} cf_challenge={cf} "
        f"bytes={len(r.text)} -> {'PASS' if ok else 'blocked'}"
    )
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--doi", default="10.1111/jce.14865", help="Wiley DOI to probe.")
    args = ap.parse_args()

    print(
        f"EZproxy configured: {ezproxy.is_configured()}  prefix={ezproxy.proxy_prefix() or '(none)'}"
    )
    if not ezproxy.is_configured():
        print(
            "\nSet GVF_EZPROXY_PREFIX (or GVF_EZPROXY_HOST) and ensure an EZproxy\n"
            "session cookie is loaded (log into the library proxy once in Chrome;\n"
            "cookie_loader reads ezproxy.library.vanderbilt.edu). Then re-run."
        )
        return 2

    # Build a session with browser cookies (incl. the EZproxy session cookie) and
    # the EZproxy URL rewriter installed — exactly what fetch_paywalled.py uses.
    from scripts.fetch_paywalled import (
        make_session,
        hydrate_session_with_browser_cookies,
    )
    from harvesting.browser_html.cookie_loader import load_chrome_cookies

    raw = make_session()  # already has EZproxy installed
    try:
        cookies = load_chrome_cookies()
        n = hydrate_session_with_browser_cookies(raw, cookies)
        ez = sum(1 for c in cookies if "ezproxy" in (c.get("domain") or "").lower())
        print(f"loaded {n} browser cookies into the session ({ez} EZproxy)")
        if ez == 0:
            print(
                "  WARNING: no EZproxy session cookie found — log into the proxy in Chrome."
            )
    except Exception as exc:  # noqa: BLE001
        print(f"  (could not load browser cookies: {exc})")

    url = f"https://onlinelibrary.wiley.com/doi/full/{args.doi}"
    print(f"\nProbing Wiley DOI {args.doi}:")
    # direct (rewriter no-ops only if unconfigured; here it WILL route — so to show
    # the contrast we hit the raw URL with a plain session too)
    import requests

    plain = requests.Session()
    plain.headers.update(raw.headers)
    _probe(plain, url, "direct (no proxy)")
    routed_ok = _probe(raw, url, "via EZproxy")

    print(
        "\nRESULT: EZproxy routing "
        + ("WORKS — Cloudflare cleared." if routed_ok else "did NOT clear CF.")
    )
    if not routed_ok:
        print(
            "Most likely the EZproxy session cookie is missing/expired. Log into\n"
            "the library proxy in Chrome once, then re-run. If it still fails, the\n"
            "proxy may require host-rewrite form instead of login?url= — set\n"
            "GVF_EZPROXY_PREFIX to the working form from a manual browser session."
        )
    return 0 if routed_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
