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
                       "https://proxy.library.vanderbilt.edu/login?url="
  GVF_EZPROXY_HOST     just the host, e.g. "proxy.library.vanderbilt.edu"
                       (the "/login?url=" prefix is built for you)

Optional:
  GVF_EZPROXY_ALL=1    proxy *every* publisher URL, not just the CF-blocked set.

You also need a valid EZproxy **session cookie** in the requests session /
browser profile (one-time SSO; cookies for `proxy.library.vanderbilt.edu` are
already loaded by ``cookie_loader.py``). With the cookie present, routing is
fully automated until the SSO session expires.
"""

from __future__ import annotations

import os
from urllib.parse import quote, urlparse

# Registrable domains that are Cloudflare-fronted / subscription-walled and worth
# routing through the institutional proxy. Subdomains match by suffix.
CF_BLOCKED_DOMAINS: tuple[str, ...] = (
    "onlinelibrary.wiley.com",
    "wiley.com",
    "karger.com",
    "sagepub.com",
    "journals.sagepub.com",
    "liebertpub.com",
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
    # never double-wrap a URL already pointing at the proxy host
    prefix_host = _host(proxy_prefix() + "x")
    if prefix_host and prefix_host in host:
        return False
    if (os.environ.get("GVF_EZPROXY_ALL") or "").strip() in ("1", "true", "yes"):
        return True
    return any(host == d or host.endswith("." + d) for d in CF_BLOCKED_DOMAINS)


def wrap(url: str) -> str:
    """Return the EZproxy-rewritten URL when proxying applies, else *url*."""
    if not should_proxy(url):
        return url
    return proxy_prefix() + quote(url, safe="")


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
