"""Unit tests for institutional EZproxy routing of CF-blocked publishers."""

import requests

from harvesting.browser_html import ezproxy

VU = "https://ezproxy.library.vanderbilt.edu/login?url="


def test_inert_when_unconfigured(monkeypatch):
    monkeypatch.delenv("GVF_EZPROXY_PREFIX", raising=False)
    monkeypatch.delenv("GVF_EZPROXY_HOST", raising=False)
    assert not ezproxy.is_configured()
    url = "https://onlinelibrary.wiley.com/doi/10.1111/jce.14865"
    assert ezproxy.wrap(url) == url  # unchanged
    assert ezproxy.should_proxy(url) is False


def test_wraps_cf_blocked_only(monkeypatch):
    monkeypatch.setenv("GVF_EZPROXY_PREFIX", VU)
    monkeypatch.delenv("GVF_EZPROXY_ALL", raising=False)
    wiley = "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fx&file=y"
    assert ezproxy.should_proxy(wiley)
    wrapped = ezproxy.wrap(wiley)
    assert wrapped.startswith(VU) and "onlinelibrary.wiley.com" in wrapped
    # an open / non-blocked host is left alone unless GVF_EZPROXY_ALL
    cdn = "https://ars.els-cdn.com/content/image/1-s2.0-x-mmc1.docx"
    assert ezproxy.should_proxy(cdn) is False
    # never double-wrap a URL already on the proxy host
    assert ezproxy.should_proxy(wrapped) is False
    assert ezproxy.wrap(wrapped) == wrapped


def test_host_only_builds_prefix(monkeypatch):
    monkeypatch.delenv("GVF_EZPROXY_PREFIX", raising=False)
    monkeypatch.setenv("GVF_EZPROXY_HOST", "ezproxy.library.vanderbilt.edu")
    assert ezproxy.wrap("https://karger.com/x").startswith(
        "https://ezproxy.library.vanderbilt.edu/login?url="
    )


def test_ezproxy_all_proxies_everything(monkeypatch):
    monkeypatch.setenv("GVF_EZPROXY_PREFIX", VU)
    monkeypatch.setenv("GVF_EZPROXY_ALL", "1")
    assert ezproxy.should_proxy("https://www.sciencedirect.com/science/article/pii/X")


def test_install_on_session_rewrites_requests(monkeypatch):
    monkeypatch.setenv("GVF_EZPROXY_PREFIX", VU)
    seen = {}

    class _S(requests.Session):
        def request(self, method, url, *a, **kw):
            seen["url"] = url
            return "ok"

    s = _S()
    ezproxy.install_on_session(s)
    s.get("https://onlinelibrary.wiley.com/doi/10.1111/x")
    assert seen["url"].startswith(VU)
    # idempotent: a second install doesn't double-wrap
    ezproxy.install_on_session(s)
    s.get("https://onlinelibrary.wiley.com/doi/10.1111/y")
    assert seen["url"].count("ezproxy.library.vanderbilt.edu") == 1
