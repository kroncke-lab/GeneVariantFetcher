"""Load cookies from a local Chrome profile and convert to Playwright format.

Used by AuthenticatedBrowserPool so headless Playwright inherits the user's
publisher SSO sessions (e.g. VUMC EZproxy / institutional access).

Reads Chrome's cookie store via `browser_cookie3` — no need to attach to a
running Chrome process or copy the user-data-dir, which avoids profile-lock
conflicts when Chrome is open.
"""

from __future__ import annotations

import logging
import multiprocessing
import os
import queue
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
    # Institutional SSO / proxy — needed for the auth handshake redirect chain
    # AND for the EZproxy session cookie that authorizes proxied publisher
    # access (Wiley/Karger/Sage routing in harvesting/browser_html/ezproxy.py).
    "vumc.org",
    "vanderbilt.edu",
    # The live Vanderbilt library proxy. EZproxy host-rewriting drops the session
    # cookie on the apex `.proxy.library.vanderbilt.edu`, so querying the apex
    # catches it; the actual *login* host is the `login.` subdomain (the bare apex
    # has no matching cert). Both listed so host-only cookies are caught too.
    "proxy.library.vanderbilt.edu",
    "login.proxy.library.vanderbilt.edu",
    "ezproxy.library.vanderbilt.edu",  # legacy alias (no longer resolves), kept for safety
    "microsoftonline.com",
    "login.microsoft.com",
    "okta.com",
    "shibboleth.net",
)


def _env_cookie_domains() -> list[str]:
    """Return user-configured cookie domains for institution-specific SSO.

    Honors ``GVF_COOKIE_DOMAINS`` / ``GVF_SSO_COOKIE_DOMAINS`` (space/comma list)
    and the plain ``COOKIE_DOMAIN`` (e.g. ``proxy.library.vanderbilt.edu``).
    """
    out: list[str] = []
    raw = os.environ.get("GVF_COOKIE_DOMAINS") or os.environ.get(
        "GVF_SSO_COOKIE_DOMAINS"
    )
    if raw:
        out.extend(d.strip() for d in _split_cookie_domains(raw) if d.strip())
    single = (os.environ.get("COOKIE_DOMAIN") or "").strip()
    if single:
        out.append(single)
    return out


def _split_cookie_domains(raw: str) -> list[str]:
    return [part for part in raw.replace(",", " ").split(" ") if part]


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


def _cookies_from_file(path: str) -> List[dict]:
    """Load cookies from a Netscape/Mozilla ``cookies.txt`` export, bypassing
    browser_cookie3 and the macOS Keychain entirely.

    Recent Chrome builds harden the "Chrome Safe Storage" keychain item so only
    Chrome itself can read the cookie-encryption key; browser_cookie3 then fails
    with "Unable to get key for cookie decryption" even in an interactive
    Terminal (and always in a headless / agent context). Exporting the cookies
    from *inside* Chrome — e.g. a "Get cookies.txt" browser extension — sidesteps
    that: Chrome decrypts its own cookies and writes a plain tab-delimited file
    we can read directly.

    Handles the ``#HttpOnly_`` line prefix that stdlib cookiejar parsers treat as
    a comment and drop — the EZproxy session cookie is typically HttpOnly, so
    losing it would defeat the purpose.
    """
    out: List[dict] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            http_only = False
            if line.startswith("#HttpOnly_"):
                http_only = True
                line = line[len("#HttpOnly_") :]
            elif line.startswith("#"):
                continue  # a genuine comment line
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            domain, _flag, cpath, secure, expires, name = parts[:6]
            value = "\t".join(parts[6:])  # tolerate a tab inside the value
            if not name:
                continue
            try:
                exp = int(expires)
            except (TypeError, ValueError):
                exp = 0
            out.append(
                {
                    "name": name,
                    "value": value,
                    "domain": domain,
                    "path": cpath or "/",
                    "secure": str(secure).upper() == "TRUE",
                    "httpOnly": http_only,
                    "expires": exp if exp > 0 else -1,
                }
            )
    return out


def _load_domain_cookie_dicts(
    domain: str,
    profile_name: Optional[str],
) -> List[dict]:
    """Load and convert cookies for one domain.

    Kept at module scope so it can run in a short-lived child process. Chrome
    cookie decryption can hang indefinitely on macOS when Keychain is locked or
    denying non-interactive access; isolating one domain lets the caller kill
    that attempt without losing the entire recovery run.
    """
    import browser_cookie3 as bc

    if profile_name:
        jar = bc.chrome(domain_name=domain, profile=profile_name)
    else:
        jar = bc.chrome(domain_name=domain)
    return [pw for c in jar if (pw := chrome_cookie_to_playwright(c)) is not None]


def _load_domain_cookie_dicts_worker(
    domain: str,
    profile_name: Optional[str],
    out_queue,
) -> None:
    try:
        out_queue.put({"cookies": _load_domain_cookie_dicts(domain, profile_name)})
    except Exception as exc:  # pragma: no cover - exercised through parent process
        out_queue.put({"error": str(exc)})


def _cookie_domain_matches(cookie_domain: str, domain_list) -> bool:
    """True if a cookie's domain equals, or is a subdomain of, any requested domain."""
    cd = (cookie_domain or "").lower().lstrip(".")
    for d in domain_list:
        dl = (d or "").lower().lstrip(".")
        if dl and (cd == dl or cd.endswith("." + dl)):
            return True
    return False


def load_chrome_cookies(
    domains: Optional[Iterable[str]] = None,
    profile_name: Optional[str] = None,
    timeout_seconds: Optional[float] = None,
) -> List[dict]:
    """Load Chrome cookies for the given domains, in Playwright format.

    Args:
        domains: registrable domains to fetch cookies for. Defaults to the
            publisher + SSO list we maintain in this module.
        profile_name: optional Chrome profile folder name (e.g. ``"Default"``,
            ``"Profile 1"``). When omitted, browser_cookie3 walks all profiles
            and merges results — that's usually what we want, since a user may
            be signed into the publisher in any profile.
        timeout_seconds: optional per-domain timeout. When set, each domain is
            loaded in an isolated child process so a stuck macOS Keychain
            lookup can be terminated instead of hanging the full long run.

    Returns:
        List of cookie dicts ready for ``BrowserContext.add_cookies()``.
        Failures to read a given domain are logged and skipped (we don't want
        one bad domain to disable the whole fetcher).
    """
    # Resolve the effective domain set once; both the cookie-file path and the
    # browser_cookie3 path filter to it.
    if domains is not None:
        domains_iter = list(domains)
    else:
        domains_iter = list(
            dict.fromkeys((*DEFAULT_PUBLISHER_DOMAINS, *_env_cookie_domains()))
        )

    # A GVF_COOKIE_FILE (Netscape cookies.txt exported from inside Chrome) takes
    # precedence and skips browser_cookie3/Keychain entirely — the only reliable
    # path on Chrome builds that lock the Safe Storage keychain item to Chrome,
    # and the only path that works headless / through an agent. Filter to the
    # requested domains so an untrimmed export can't inject unrelated cookies
    # (banking, mail, ...) into the recovery session / browser pool.
    cookie_file = os.environ.get("GVF_COOKIE_FILE")
    if cookie_file:
        p = os.path.expanduser(cookie_file)
        if os.path.exists(p):
            try:
                cookies = [
                    c
                    for c in _cookies_from_file(p)
                    if _cookie_domain_matches(c.get("domain", ""), domains_iter)
                ]
                logger.info(
                    "Loaded %d cookies from GVF_COOKIE_FILE=%s (Keychain bypassed)",
                    len(cookies),
                    p,
                )
                return cookies
            except Exception as e:  # noqa: BLE001 - fall back to Chrome on parse error
                logger.warning(
                    "GVF_COOKIE_FILE=%s could not be parsed (%s); falling back to Chrome",
                    p,
                    e,
                )
        else:
            logger.warning(
                "GVF_COOKIE_FILE=%s does not exist; falling back to Chrome cookies", p
            )

    # An explicit empty domain set has nothing to load from Chrome. Return before
    # importing the optional browser_cookie3 dependency so this no-op path stays
    # usable in base/CI installs that intentionally omit the browser extra.
    if not domains_iter:
        return []

    try:
        import browser_cookie3 as bc
    except ImportError as e:
        raise RuntimeError(
            "browser_cookie3 is required for authenticated fetching. "
            "Install with `pip install browser-cookie3`."
        ) from e

    seen_keys: set = set()
    out: List[dict] = []

    for d in domains_iter:
        try:
            if timeout_seconds and timeout_seconds > 0:
                ctx = multiprocessing.get_context("spawn")
                q = ctx.Queue()
                proc = ctx.Process(
                    target=_load_domain_cookie_dicts_worker,
                    args=(d, profile_name, q),
                )
                proc.start()
                proc.join(timeout_seconds)
                if proc.is_alive():
                    proc.terminate()
                    proc.join(timeout=2)
                    logger.warning(
                        "cookie load timed out for %s after %.1fs",
                        d,
                        timeout_seconds,
                    )
                    continue
                try:
                    payload = q.get_nowait()
                except queue.Empty:
                    logger.warning("cookie load failed for %s: no child result", d)
                    continue
                if payload.get("error"):
                    logger.warning("cookie load failed for %s: %s", d, payload["error"])
                    continue
                cookie_dicts = payload.get("cookies") or []
            else:
                cookie_dicts = _load_domain_cookie_dicts(d, profile_name)
        except Exception as e:
            logger.warning("cookie load failed for %s: %s", d, e)
            continue

        count_before = len(out)
        for pw in cookie_dicts:
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
