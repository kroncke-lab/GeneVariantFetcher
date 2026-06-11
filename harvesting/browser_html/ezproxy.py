"""Route Cloudflare-fronted publisher requests through an institutional EZproxy.

Wiley Online Library (and Karger / Sage) sit behind a Cloudflare *managed
challenge* that blocks plain `requests` and headless browsers (HTTP 403 "Just a
moment…"). The robust, authorized way past it for a **licensed institution** is
to route the request through the library EZproxy: its egress IP is allowlisted by
the publisher as a paying subscriber, so the proxied request returns the licensed
full text / supplement without a CF challenge.

This module is **inert unless configured** (no hardcoded institution is ever
activated). Set one of:

  GVF_EZPROXY_PREFIX   full login-rewrite prefix, e.g.
                       "https://login.proxy.library.vanderbilt.edu/login?url="
  GVF_EZPROXY_HOST     just the host, e.g. "login.proxy.library.vanderbilt.edu"
                       (the "/login?url=" prefix is built for you)

Use the institution's actual *login* host. Vanderbilt's EZproxy uses host-based
rewriting behind a wildcard cert (`*.proxy.library.vanderbilt.edu`) that does NOT
match the bare apex, so `proxy.library.vanderbilt.edu` throws a TLS error — the
working host is `login.proxy.library.vanderbilt.edu` (verified live 2026-06-08).
The session cookie still lands on the apex `.proxy.library.vanderbilt.edu`, which
is what `cookie_loader.py` reads.

Optional:
  GVF_EZPROXY_ALL=1    proxy *every* publisher URL, not just the CF-blocked set.

You also need a valid EZproxy **session cookie** in the requests session /
browser profile (one-time SSO; cookies for `proxy.library.vanderbilt.edu` are
already loaded by ``cookie_loader.py``). With the cookie present, routing is
fully automated until the SSO session expires.
"""

from __future__ import annotations

import os
from urllib.parse import urlparse, urlsplit, urlunsplit

# Registrable domains that are Cloudflare-fronted / subscription-walled and worth
# routing through the institutional proxy. Subdomains match by suffix.
CF_BLOCKED_DOMAINS: tuple[str, ...] = (
    "onlinelibrary.wiley.com",
    "wiley.com",
    "karger.com",
    "sagepub.com",
    "journals.sagepub.com",
    "liebertpub.com",
    # AHA Journals (Circulation, Circ Cardiovasc Genet, JAHA, Stroke, ...) —
    # DOI 10.1161, Cloudflare-fronted and subscription-walled for <12-month and
    # many legacy articles. Proxy so the subscriber IP returns licensed full text
    # instead of a "get full access to this article" stub.
    "ahajournals.org",
)

# Resolution / API / CDN / infrastructure hosts that must NEVER be routed through
# EZproxy — they are not subscription resources, and proxying them breaks them.
# In particular GVF_EZPROXY_ALL=1 would otherwise send NCBI E-utilities DOI
# lookups through the proxy, which returns a non-JSON page and silently drops
# every PMID as "no DOI". Checked before both the ALL=1 and CF-list rules.
NEVER_PROXY_DOMAINS: tuple[str, ...] = (
    "ncbi.nlm.nih.gov",  # E-utilities (esummary/efetch), PMC, FTP
    "nih.gov",
    "doi.org",  # DOI resolver — follow the redirect directly, then proxy the target
    "crossref.org",
    "els-cdn.com",  # Elsevier open supplement CDN (ars.els-cdn.com) — served openly
    "api.elsevier.com",  # Elsevier full-text API (uses insttoken, not the proxy)
    "id.elsevier.com",  # Elsevier SSO
    "unpaywall.org",
    "googleapis.com",
)


def proxy_prefix() -> str:
    """The configured EZproxy login-rewrite prefix, or '' when unconfigured.

    Accepts ``GVF_EZPROXY_PREFIX`` or the plainer ``PROXY_LOGIN_PREFIX``; or a
    bare host via ``GVF_EZPROXY_HOST`` / ``PROXY_HOST`` (the ``/login?url=``
    prefix is built for you). e.g. PROXY_HOST=proxy.library.vanderbilt.edu.
    """
    raw = (
        os.environ.get("GVF_EZPROXY_PREFIX")
        or os.environ.get("PROXY_LOGIN_PREFIX")
        or ""
    ).strip()
    if raw:
        return raw
    host = (
        (os.environ.get("GVF_EZPROXY_HOST") or os.environ.get("PROXY_HOST") or "")
        .strip()
        .rstrip("/")
    )
    if host:
        if "://" not in host:
            host = f"https://{host}"
        return f"{host}/login?url="
    return ""


def is_configured() -> bool:
    return bool(proxy_prefix())


def proxy_base() -> str:
    """The proxy's base host used for host-rewriting, or '' when unconfigured.

    EZproxy host-based rewriting turns a resource into a subdomain of this base
    (e.g. ``onlinelibrary.wiley.com`` -> ``onlinelibrary-wiley-com.<base>``). For
    Vanderbilt the base is ``proxy.library.vanderbilt.edu`` while the *login* host
    is ``login.proxy.library.vanderbilt.edu`` — so a leading ``login.`` label is
    stripped off the configured host. Override directly with ``GVF_EZPROXY_BASE``
    / ``PROXY_BASE`` for libraries whose base differs from ``host`` minus
    ``login.``.
    """
    explicit = (
        os.environ.get("GVF_EZPROXY_BASE") or os.environ.get("PROXY_BASE") or ""
    ).strip()
    if explicit:
        return _host(explicit if "://" in explicit else f"https://{explicit}")
    host = _host(proxy_prefix() + "x")
    if host.startswith("login."):
        host = host[len("login.") :]
    return host


def _host(url: str) -> str:
    try:
        return (urlparse(url).hostname or "").lower()
    except ValueError:
        return ""


def should_proxy(url: str) -> bool:
    """True if *url* should be routed through the configured EZproxy."""
    if not is_configured() or not url:
        return False
    host = _host(url)
    if not host:
        return False
    # never double-wrap a URL already on the proxy (login host OR a rewritten
    # resource subdomain both end in the base, e.g. *.proxy.library.vanderbilt.edu)
    base = proxy_base()
    if base and (host == base or host.endswith("." + base)):
        return False
    # resolution/API/CDN/infra hosts are never proxied (even with GVF_EZPROXY_ALL)
    if any(host == d or host.endswith("." + d) for d in NEVER_PROXY_DOMAINS):
        return False
    if (os.environ.get("GVF_EZPROXY_ALL") or "").strip() in ("1", "true", "yes"):
        return True
    return any(host == d or host.endswith("." + d) for d in CF_BLOCKED_DOMAINS)


def _hostify(host: str, base: str) -> str:
    """EZproxy host-based rewrite of a single hostname.

    Standard EZproxy encoding: literal dashes are doubled, then dots become
    dashes, then the proxy base is appended. ``onlinelibrary.wiley.com`` ->
    ``onlinelibrary-wiley-com.<base>``; ``my-host.com`` -> ``my--host-com.<base>``.
    """
    label = host.replace("-", "--").replace(".", "-")
    return f"{label}.{base}"


def wrap(url: str) -> str:
    """Return the EZproxy host-rewritten URL when proxying applies, else *url*.

    Uses EZproxy *host-based rewriting* (the resource becomes a subdomain of the
    proxy base), NOT the ``login?url=`` starting-point form. On an
    already-authenticated session the ``login?url=`` form redirects to the proxy
    *menu* instead of the resource (observed live on Vanderbilt's EZproxy
    2026-06-08), so it silently fetches the menu page and is unusable for
    unattended harvesting. The host-rewritten form returns the resource directly,
    cold, with just the session cookie.
    """
    if not should_proxy(url):
        return url
    base = proxy_base()
    if not base:
        return url
    parts = urlsplit(url)
    if not parts.hostname:
        return url
    new_netloc = _hostify(parts.hostname, base)
    if parts.port:
        new_netloc = f"{new_netloc}:{parts.port}"
    return urlunsplit(
        (parts.scheme or "https", new_netloc, parts.path, parts.query, parts.fragment)
    )


def install_on_session(session):
    """Make a ``requests.Session`` auto-route CF-blocked URLs through EZproxy.

    Idempotent and a no-op when unconfigured. Wraps ``Session.request`` so every
    ``.get``/``.post``/streamed download of a CF-blocked publisher URL (article
    HTML *and* `/action/downloadSupplement` files) is rewritten to the proxy form
    — no per-call-site changes needed.
    """
    if session is None or getattr(session, "_gvf_ezproxy_installed", False):
        return session
    if not is_configured():
        return session
    original_request = session.request

    def _request(method, url, *args, **kwargs):
        return original_request(method, wrap(url), *args, **kwargs)

    session.request = _request
    session._gvf_ezproxy_installed = True
    return session
