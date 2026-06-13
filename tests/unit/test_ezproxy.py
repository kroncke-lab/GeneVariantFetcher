"""Unit tests for institutional EZproxy routing of CF-blocked publishers.

``wrap()`` uses EZproxy *host-based rewriting* (resource -> subdomain of the proxy
base), NOT the ``login?url=`` starting-point form: on an authenticated session the
latter bounces to the proxy menu (verified live on Vanderbilt's EZproxy), so it
can't be used for unattended fetch.
"""

import requests

from harvesting.browser_html import ezproxy

# Configured value is the *login* host; the rewrite base is that minus "login.".
VU = "https://login.proxy.library.vanderbilt.edu/login?url="
BASE = "proxy.library.vanderbilt.edu"


def test_inert_when_unconfigured(monkeypatch):
    monkeypatch.delenv("GVF_EZPROXY_PREFIX", raising=False)
    monkeypatch.delenv("GVF_EZPROXY_HOST", raising=False)
    monkeypatch.delenv("PROXY_LOGIN_PREFIX", raising=False)
    monkeypatch.delenv("PROXY_HOST", raising=False)
    assert not ezproxy.is_configured()
    url = "https://onlinelibrary.wiley.com/doi/10.1111/jce.14865"
    assert ezproxy.wrap(url) == url  # unchanged
    assert ezproxy.should_proxy(url) is False


def test_proxy_base_strips_login_label(monkeypatch):
    monkeypatch.setenv("GVF_EZPROXY_PREFIX", VU)
    monkeypatch.delenv("GVF_EZPROXY_BASE", raising=False)
    assert ezproxy.proxy_base() == BASE


def test_wraps_cf_blocked_to_host_rewrite(monkeypatch):
    monkeypatch.setenv("GVF_EZPROXY_PREFIX", VU)
    monkeypatch.delenv("GVF_EZPROXY_ALL", raising=False)
    monkeypatch.delenv("GVF_EZPROXY_BASE", raising=False)
    wiley = "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fx&file=y"
    assert ezproxy.should_proxy(wiley)
    wrapped = ezproxy.wrap(wiley)
    # host-rewritten onto the proxy base, NOT a login?url= prefix
    assert wrapped.startswith(f"https://onlinelibrary-wiley-com.{BASE}/")
    assert "/login?url=" not in wrapped
    # path + query are preserved verbatim (supplement downloads depend on this)
    assert wrapped.endswith("/action/downloadSupplement?doi=10.1111%2Fx&file=y")
    # an open / non-blocked host is left alone unless GVF_EZPROXY_ALL
    cdn = "https://ars.els-cdn.com/content/image/1-s2.0-x-mmc1.docx"
    assert ezproxy.should_proxy(cdn) is False
    # never double-wrap a URL already on the proxy host
    assert ezproxy.should_proxy(wrapped) is False
    assert ezproxy.wrap(wrapped) == wrapped


def test_host_only_config_rewrites(monkeypatch):
    monkeypatch.delenv("GVF_EZPROXY_PREFIX", raising=False)
    monkeypatch.setenv("GVF_EZPROXY_HOST", "login.proxy.library.vanderbilt.edu")
    assert ezproxy.wrap("https://karger.com/x") == f"https://karger-com.{BASE}/x"


def test_dashes_in_host_are_doubled(monkeypatch):
    # standard EZproxy encoding: literal dashes double, then dots -> dashes
    monkeypatch.setenv("GVF_EZPROXY_PREFIX", VU)
    monkeypatch.setenv("GVF_EZPROXY_ALL", "1")
    assert ezproxy.wrap("https://my-host.example.com/p") == (
        f"https://my--host-example-com.{BASE}/p"
    )


def test_ezproxy_all_proxies_everything(monkeypatch):
    monkeypatch.setenv("GVF_EZPROXY_PREFIX", VU)
    monkeypatch.setenv("GVF_EZPROXY_ALL", "1")
    assert ezproxy.should_proxy("https://www.sciencedirect.com/science/article/pii/X")
    assert ezproxy.wrap("https://www.sciencedirect.com/science/article/pii/X") == (
        f"https://www-sciencedirect-com.{BASE}/science/article/pii/X"
    )


def test_install_on_session_rewrites_requests(monkeypatch):
    monkeypatch.setenv("GVF_EZPROXY_PREFIX", VU)
    monkeypatch.delenv("GVF_EZPROXY_ALL", raising=False)
    seen = {}

    class _S(requests.Session):
        def request(self, method, url, *a, **kw):
            seen["url"] = url
            return "ok"

    s = _S()
    ezproxy.install_on_session(s)
    s.get("https://onlinelibrary.wiley.com/doi/10.1111/x")
    assert seen["url"] == f"https://onlinelibrary-wiley-com.{BASE}/doi/10.1111/x"
    # idempotent: a second install doesn't double-wrap
    ezproxy.install_on_session(s)
    s.get("https://onlinelibrary.wiley.com/doi/10.1111/y")
    assert seen["url"].count(BASE) == 1


def test_honors_plain_proxy_env_var_aliases(monkeypatch):
    # the user's .env may use PROXY_LOGIN_PREFIX / PROXY_HOST (not GVF_-prefixed)
    monkeypatch.delenv("GVF_EZPROXY_PREFIX", raising=False)
    monkeypatch.delenv("GVF_EZPROXY_HOST", raising=False)
    monkeypatch.delenv("GVF_EZPROXY_BASE", raising=False)
    monkeypatch.setenv(
        "PROXY_LOGIN_PREFIX", "https://login.proxy.library.vanderbilt.edu/login?url="
    )
    assert ezproxy.is_configured()
    assert ezproxy.wrap("https://onlinelibrary.wiley.com/doi/x") == (
        f"https://onlinelibrary-wiley-com.{BASE}/doi/x"
    )
    monkeypatch.delenv("PROXY_LOGIN_PREFIX", raising=False)
    monkeypatch.setenv("PROXY_HOST", "login.proxy.library.vanderbilt.edu")
    assert ezproxy.proxy_base() == BASE


def test_never_proxy_infra_hosts_even_with_all(monkeypatch):
    # GVF_EZPROXY_ALL must NOT route NCBI/DOI/API/CDN hosts (would break DOI lookup)
    monkeypatch.setenv("GVF_EZPROXY_PREFIX", VU)
    monkeypatch.setenv("GVF_EZPROXY_ALL", "1")
    for u in (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=1",
        "https://doi.org/10.1111/jce.14865",
        "https://api.elsevier.com/content/article/doi/10.1016/x",
        "https://ars.els-cdn.com/content/image/1-s2.0-x-mmc1.docx",
    ):
        assert ezproxy.should_proxy(u) is False
        assert ezproxy.wrap(u) == u
    # but a real publisher host still proxies under ALL=1
    assert ezproxy.should_proxy("https://www.nature.com/articles/x")


def test_explicit_base_override(monkeypatch):
    # GVF_EZPROXY_BASE wins for libraries whose base != host minus "login."
    monkeypatch.setenv("GVF_EZPROXY_HOST", "signin.example.edu")
    monkeypatch.setenv("GVF_EZPROXY_BASE", "ezp.example.edu")
    monkeypatch.setenv("GVF_EZPROXY_ALL", "1")
    assert ezproxy.proxy_base() == "ezp.example.edu"
    assert (
        ezproxy.wrap("https://karger.com/x") == "https://karger-com.ezp.example.edu/x"
    )
