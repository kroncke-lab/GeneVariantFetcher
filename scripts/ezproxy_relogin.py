#!/usr/bin/env python3
"""Refresh the EZproxy session cookie without re-doing SSO by hand.

The recurring preflight exit-5 is the *short-lived* EZproxy session cookie
dying server-side while the *long-lived* VUMC Microsoft SSO session stays
valid. That asymmetry makes the refresh silently automatable: a dedicated
persistent browser profile holds the SSO tokens, so driving it through
``login.proxy.library.vanderbilt.edu/login?url=...`` completes the whole
redirect chain without a human and mints a fresh EZproxy session cookie.
(The ``login?url=`` starting-point form is deliberate here — it triggers the
SSO chain. Unattended *fetching* keeps using the host-rewritten form; see
``harvesting/browser_html/ezproxy.py::wrap``.)

The fresh cookies are merged into ``GVF_COOKIE_FILE`` (Netscape cookies.txt —
the one file every GVF cookie consumer reads: the requests session, the
Playwright pool, and the gvf-run preflight probe), then verified with the same
live probe gvf-run uses.

Operator flow:

  one-time bootstrap (headed browser; complete the VUMC Microsoft MFA once):
      .venv/bin/python scripts/ezproxy_relogin.py --bootstrap

  every later refresh (headless, no interaction):
      .venv/bin/python scripts/ezproxy_relogin.py

Exit codes:
  0  session refreshed and the live full-text probe PASSED
  1  unexpected error
  2  not configured (no GVF_EZPROXY_PREFIX/HOST, or no cookie-file path)
  3  SSO interaction required — run again with --bootstrap
  4  cookies refreshed but the live probe is still blocked

The persistent profile holds live SSO tokens: it is created ``chmod 700``,
must stay untracked, and defaults to ``~/.gvf/ezproxy_profile`` (override with
``GVF_EZPROXY_PROFILE_DIR`` or ``--profile-dir``).
"""

from __future__ import annotations

import argparse
import os
import stat
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from harvesting.browser_html import ezproxy  # noqa: E402
from harvesting.browser_html.authenticated_pool import (  # noqa: E402
    AuthenticatedBrowserPool,
)
from harvesting.browser_html.cookie_loader import _cookies_from_file  # noqa: E402

DEFAULT_PROFILE_DIR = "~/.gvf/ezproxy_profile"
PROFILE_DIR_ENV = "GVF_EZPROXY_PROFILE_DIR"

# How long a headless refresh waits for the silent SSO redirect chain to mint
# a session cookie, and how long a headed bootstrap waits for the human to
# finish MFA.
DEFAULT_REFRESH_WAIT_S = 60
DEFAULT_BOOTSTRAP_WAIT_S = 300

# A final URL on one of these hosts means the flow is parked at an identity
# provider waiting for credentials/MFA — headless cannot proceed. Narrower than
# the preflight's _LOGIN_URL_MARKERS on purpose: the EZproxy *login host* itself
# contains "login", and landing there with a fresh session cookie is success.
_IDP_URL_MARKERS = (
    "microsoftonline",
    "login.microsoft",
    "adfs",
    "okta",
    "shibboleth",
    "wayf",
    "openathens",
    "duosecurity",
)


def effective_wait_s(headed: bool, wait_s: Optional[int]) -> int:
    """A human doing MFA gets minutes; a silent redirect chain gets seconds."""
    if wait_s is not None:
        return wait_s
    return DEFAULT_BOOTSTRAP_WAIT_S if headed else DEFAULT_REFRESH_WAIT_S


def resolve_profile_dir(explicit: Optional[str] = None) -> Path:
    """The persistent-profile directory (not created here)."""
    raw = explicit or os.environ.get(PROFILE_DIR_ENV) or DEFAULT_PROFILE_DIR
    return Path(raw).expanduser()


def profile_is_bootstrapped(profile_dir: Path) -> bool:
    """True when the profile exists and has ever been used (holds SSO state)."""
    try:
        return profile_dir.is_dir() and any(profile_dir.iterdir())
    except OSError:
        return False


def resolve_cookie_file(explicit: Optional[str] = None) -> Optional[Path]:
    raw = explicit or os.environ.get("GVF_COOKIE_FILE")
    if not raw:
        return None
    return Path(raw).expanduser()


def _is_proxy_scoped_cookie(cookie: dict) -> bool:
    """Broad scope used for FILE selection: anything on the proxy base domain
    (the apex session cookie plus host-rewritten publisher subdomain cookies)
    or any ``ezproxy*``-named cookie."""
    name = (cookie.get("name") or "").lower()
    domain = (cookie.get("domain") or "").lower().lstrip(".")
    if name.startswith("ezproxy"):
        return True
    base = ezproxy.proxy_base()
    return bool(base) and (domain == base or domain.endswith("." + base))


def _is_session_cookie(cookie: dict) -> bool:
    """STRICT success detection: only the ``ezproxy*``-NAMED cookie counts.

    Deliberately narrower than :func:`_is_proxy_scoped_cookie` — after one
    bootstrap the persistent profile durably holds publisher/analytics cookies
    on host-rewritten ``*.proxy...`` subdomains, and the SSO chain can set
    pre-auth infrastructure cookies (Shibboleth state etc.) on proxy hosts.
    Counting any of those as "the session" made a dead session look refreshed
    (found in adversarial review). If an institution renames its session
    cookie away from ``ezproxy*`` the refresh times out and falls back to the
    bootstrap guidance — a safe failure, arbitrated by the live probe anyway.
    """
    return (cookie.get("name") or "").lower().startswith("ezproxy")


def _looks_parked_at_idp(url: str) -> bool:
    low = (url or "").lower()
    return any(m in low for m in _IDP_URL_MARKERS)


def _url_host_path(url: str) -> tuple[str, str]:
    try:
        from urllib.parse import urlparse

        parts = urlparse(url)
        return (parts.hostname or "").lower(), (parts.path or "").lower()
    except ValueError:
        return "", ""


def _on_proxy_host(url: str) -> bool:
    base = ezproxy.proxy_base()
    if not base:
        return False
    host, _ = _url_host_path(url)
    return host == base or host.endswith("." + base)


def _is_proxy_login_page(url: str) -> bool:
    """The proxy's own /login endpoint — where an EXPIRED session lands as a
    plain HTTP-200 SAML/SSO handoff page. Being here is not success."""
    _, path = _url_host_path(url)
    return _on_proxy_host(url) and path.startswith("/login")


def _at_proxy_content(url: str) -> bool:
    """The chain finished on proxy-served CONTENT (menu or a host-rewritten
    resource), not the login endpoint and not an identity provider — i.e.
    EZproxy actually accepted the session."""
    return (
        _on_proxy_host(url)
        and not _is_proxy_login_page(url)
        and not _looks_parked_at_idp(url)
    )


def select_relevant_cookies(
    cookies: List[dict], *, all_domains: bool = False
) -> List[dict]:
    """Cookies worth persisting to the shared cookie file.

    By default only proxy-domain cookies (the EZproxy session plus anything on
    the host-rewritten publisher subdomains) and ``ezproxy*``-named cookies are
    kept. Identity-provider tokens (Microsoft SSO etc.) deliberately stay in
    the browser profile — they are useless to a requests session and live
    credentials do not belong in a plain-text file.
    """
    if all_domains:
        return list(cookies)
    return [c for c in cookies if _is_proxy_scoped_cookie(c)]


def _cookie_key(cookie: dict) -> tuple:
    return (
        (cookie.get("domain") or "").lower().lstrip("."),
        cookie.get("name") or "",
        cookie.get("path") or "/",
    )


def merge_cookies(existing: List[dict], fresh: List[dict]) -> List[dict]:
    """Merge *fresh* cookies over *existing* ones by (domain, name, path).

    Existing entries the refresh did not touch are preserved — once
    ``GVF_COOKIE_FILE`` exists it is the ONLY cookie source (the loader never
    falls back to Chrome for an empty match), so dropping unrelated entries
    would silently revoke publisher sessions.
    """
    merged = {_cookie_key(c): c for c in existing}
    for c in fresh:
        merged[_cookie_key(c)] = c
    return list(merged.values())


def _sanitize(value: str) -> str:
    """Strip characters that would corrupt the tab-delimited line format."""
    return (value or "").replace("\n", "").replace("\r", "")


def playwright_cookies_to_netscape_lines(cookies: List[dict]) -> List[str]:
    """Render Playwright cookie dicts as Netscape ``cookies.txt`` lines.

    Exact inverse of ``cookie_loader._cookies_from_file``: 7 tab-separated
    fields, ``#HttpOnly_`` line prefix for HttpOnly cookies (the EZproxy
    session cookie usually is one — stdlib parsers drop those lines, ours
    must not), and expires ``0`` for session cookies (parsed back as ``-1``).
    """
    lines: List[str] = []
    for c in cookies:
        name = _sanitize(str(c.get("name") or ""))
        if not name:
            continue
        domain = _sanitize(str(c.get("domain") or ""))
        path = _sanitize(str(c.get("path") or "/")) or "/"
        include_subdomains = "TRUE" if domain.startswith(".") else "FALSE"
        secure = "TRUE" if c.get("secure") else "FALSE"
        try:
            expires = int(float(c.get("expires") or 0))
        except (TypeError, ValueError):
            expires = 0
        if expires < 0:
            expires = 0
        value = _sanitize(str(c.get("value") or ""))
        prefix = "#HttpOnly_" if c.get("httpOnly") else ""
        lines.append(
            f"{prefix}{domain}\t{include_subdomains}\t{path}\t{secure}"
            f"\t{expires}\t{name}\t{value}"
        )
    return lines


def write_cookie_file(path: Path, cookies: List[dict]) -> None:
    """Atomically write *cookies* as a Netscape file readable only by the user."""
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join(
        ["# Netscape HTTP Cookie File", *playwright_cookies_to_netscape_lines(cookies)]
    )
    fd, tmp_name = tempfile.mkstemp(
        dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(body + "\n")
        os.chmod(tmp_name, stat.S_IRUSR | stat.S_IWUSR)  # 0600
        os.replace(tmp_name, path)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise


def read_existing_cookies(path: Path) -> List[dict]:
    """Parse the current cookie file with the SAME parser the pipeline uses."""
    try:
        if not path.expanduser().is_file():
            return []
        return _cookies_from_file(str(path.expanduser()))
    except Exception:  # noqa: BLE001 - an unparseable file is replaced, not fatal
        return []


@dataclass
class RefreshResult:
    """Outcome of one browser pass through the EZproxy login chain."""

    ok: bool
    needs_bootstrap: bool = False
    final_url: str = ""
    cookies: List[dict] = field(default_factory=list)
    detail: str = ""


def refresh_session(
    profile_dir: Path,
    *,
    headed: bool = False,
    wait_s: Optional[int] = None,
    target_url: Optional[str] = None,
) -> RefreshResult:
    """Drive the persistent profile through ``login?url=...`` and harvest cookies.

    Success here means the browser flow completed and an EZproxy-domain cookie
    exists afterwards — whether that session actually unlocks licensed full
    text is decided by the caller via ``probe_institutional_access()``.
    """
    prefix = ezproxy.proxy_prefix()
    if not prefix:
        return RefreshResult(
            ok=False, detail="EZproxy is not configured (GVF_EZPROXY_PREFIX/HOST)."
        )
    wait_s = effective_wait_s(headed, wait_s)
    target = target_url or (
        "https://onlinelibrary.wiley.com/doi/full/"
        + (os.environ.get("GVF_PREFLIGHT_DOI") or "10.1111/jce.14865")
    )
    login_url = prefix + target

    profile_dir = profile_dir.expanduser()
    profile_dir.mkdir(parents=True, exist_ok=True)
    os.chmod(profile_dir, stat.S_IRWXU)  # 0700 — the profile holds live SSO tokens

    pool = AuthenticatedBrowserPool(
        headless=not headed,
        persistent_profile_path=str(profile_dir),
        use_chrome_channel=True,
    )
    try:
        pool.start()
        # The persistent profile RETAINS the previous (possibly dead) session
        # cookie, so a session cookie existing is NOT success on its own, and
        # neither is being on a proxy-host URL (the /login handoff page is a
        # proxy-host HTTP 200). Success = an ezproxy*-NAMED cookie exists AND
        # the redirect chain finished on proxy CONTENT (menu / rewritten
        # resource). Whether that session truly unlocks full text is decided
        # afterwards by probe_institutional_access().
        stale_values = {
            (c.get("name"), c.get("value"))
            for c in pool.context_cookies()
            if _is_session_cookie(c)
        }
        with pool.page() as page:
            page.goto(login_url, wait_until="domcontentloaded", timeout=60_000)
            deadline = time.monotonic() + wait_s
            final_url = ""
            while True:
                final_url = str(getattr(page, "url", "") or "")
                cookies = pool.context_cookies()
                session_cookies = [c for c in cookies if _is_session_cookie(c)]
                if session_cookies and _at_proxy_content(final_url):
                    minted = any(
                        (c.get("name"), c.get("value")) not in stale_values
                        for c in session_cookies
                    )
                    return RefreshResult(
                        ok=True,
                        final_url=final_url,
                        cookies=cookies,
                        detail=(
                            "fresh EZproxy session cookie minted."
                            if minted
                            else "existing EZproxy session accepted by the proxy."
                        ),
                    )
                if time.monotonic() >= deadline:
                    break
                time.sleep(1.0)
            cookies = pool.context_cookies()
            # The proxy's own /login page is where an expired session parks
            # when the silent SSO chain cannot finish — interactive territory,
            # same as an identity-provider URL.
            parked = _looks_parked_at_idp(final_url) or _is_proxy_login_page(final_url)
            return RefreshResult(
                ok=False,
                needs_bootstrap=parked or not headed,
                final_url=final_url,
                cookies=cookies,
                detail=(
                    "no fresh EZproxy session appeared within "
                    f"{wait_s}s (final URL: {final_url[:120] or 'unknown'})"
                    + (" — parked at a login/identity-provider page" if parked else "")
                ),
            )
    finally:
        pool.close()


def _print_next_steps() -> None:
    print(
        "\nSession refreshed. Typical next steps:\n"
        "  - resume paywall recovery on the latest run:\n"
        "      .venv/bin/python scripts/fetch_paywalled.py --run-dir results/<GENE>/<ts>\n"
        "      .venv/bin/python scripts/refresh_run_db.py --run-dir results/<GENE>/<ts>\n"
        "  - or just re-run gvf-run; the preflight now self-heals from this profile."
    )


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--bootstrap",
        action="store_true",
        help="open a VISIBLE browser for the one-time SSO/MFA login",
    )
    ap.add_argument(
        "--profile-dir",
        help=f"persistent browser profile (default: ${PROFILE_DIR_ENV} or {DEFAULT_PROFILE_DIR})",
    )
    ap.add_argument(
        "--cookie-file",
        help="Netscape cookies.txt to update (default: $GVF_COOKIE_FILE)",
    )
    ap.add_argument(
        "--wait-s",
        type=int,
        help=(
            "seconds to wait for the session cookie "
            f"(default: {DEFAULT_REFRESH_WAIT_S} headless, "
            f"{DEFAULT_BOOTSTRAP_WAIT_S} with --bootstrap)"
        ),
    )
    ap.add_argument(
        "--all-domains",
        action="store_true",
        help="persist every profile cookie, not just proxy-domain ones "
        "(identity-provider tokens included — normally leave this off)",
    )
    ap.add_argument(
        "--skip-verify",
        action="store_true",
        help="skip the live full-text probe after writing the cookie file",
    )
    args = ap.parse_args(argv)

    if not ezproxy.is_configured():
        print(
            "ERROR: EZproxy is not configured. Set GVF_EZPROXY_HOST "
            "(e.g. login.proxy.library.vanderbilt.edu) or GVF_EZPROXY_PREFIX in .env.",
            file=sys.stderr,
        )
        return 2
    cookie_file = resolve_cookie_file(args.cookie_file)
    if cookie_file is None:
        print(
            "ERROR: no cookie-file path. Set GVF_COOKIE_FILE in .env (the pipeline "
            "reads cookies from that file) or pass --cookie-file.",
            file=sys.stderr,
        )
        return 2

    profile_dir = resolve_profile_dir(args.profile_dir)
    if not args.bootstrap and not profile_is_bootstrapped(profile_dir):
        print(
            f"ERROR: profile {profile_dir} has never been bootstrapped.\n"
            "Run once with --bootstrap and complete the VUMC Microsoft login "
            "(@vumc.org identity) in the browser window.",
            file=sys.stderr,
        )
        return 3

    print(
        f"{'Bootstrapping' if args.bootstrap else 'Refreshing'} EZproxy session "
        f"(profile: {profile_dir})…"
    )
    if args.bootstrap:
        print(
            "A browser window will open. Complete the library/Microsoft login "
            "with your @vumc.org identity; the window closes by itself when the "
            "proxy session lands."
        )
    try:
        result = refresh_session(profile_dir, headed=args.bootstrap, wait_s=args.wait_s)
    except RuntimeError as exc:  # Playwright missing, browser launch failure
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if not result.ok:
        print(f"Refresh did not complete: {result.detail}", file=sys.stderr)
        if result.needs_bootstrap and not args.bootstrap:
            print(
                "The SSO session needs interactive login — run:\n"
                "  .venv/bin/python scripts/ezproxy_relogin.py --bootstrap",
                file=sys.stderr,
            )
            return 3
        return 1

    fresh = select_relevant_cookies(result.cookies, all_domains=args.all_domains)
    merged = merge_cookies(read_existing_cookies(cookie_file), fresh)
    write_cookie_file(cookie_file, merged)
    print(
        f"Wrote {len(fresh)} refreshed cookie(s) (file now holds {len(merged)}) "
        f"to {cookie_file}."
    )

    if args.skip_verify:
        return 0

    print("Verifying with the live full-text probe…")
    from cli.institutional_preflight import probe_institutional_access

    report = probe_institutional_access()
    for line in report.lines:
        print(f"  {line}")
    if report.viable:
        print("✅ EZproxy session refreshed and live full-text access confirmed.")
        _print_next_steps()
        return 0
    print(
        f"⚠️  cookies were refreshed but the probe is still blocked: {report.reason}",
        file=sys.stderr,
    )
    return 4


if __name__ == "__main__":
    raise SystemExit(main())
