import sys
from types import SimpleNamespace

from harvesting.browser_html import cookie_loader


def test_load_chrome_cookies_can_use_non_timeout_helper(monkeypatch):
    monkeypatch.setitem(sys.modules, "browser_cookie3", SimpleNamespace())

    def fake_load(domain, profile_name):
        assert domain == "academic.oup.com"
        assert profile_name == "Default"
        return [
            {
                "name": "session",
                "value": "abc",
                "domain": ".academic.oup.com",
                "path": "/",
                "secure": True,
                "httpOnly": True,
                "expires": -1,
            }
        ]

    monkeypatch.setattr(cookie_loader, "_load_domain_cookie_dicts", fake_load)

    cookies = cookie_loader.load_chrome_cookies(
        domains=["academic.oup.com"],
        profile_name="Default",
        timeout_seconds=None,
    )

    assert cookies == [
        {
            "name": "session",
            "value": "abc",
            "domain": ".academic.oup.com",
            "path": "/",
            "secure": True,
            "httpOnly": True,
            "expires": -1,
        }
    ]


def test_load_chrome_cookies_adds_env_sso_domains(monkeypatch):
    monkeypatch.setitem(sys.modules, "browser_cookie3", SimpleNamespace())
    monkeypatch.setenv("GVF_SSO_COOKIE_DOMAINS", "proxy.example.edu, login.example.edu")
    monkeypatch.setattr(cookie_loader, "DEFAULT_PUBLISHER_DOMAINS", ("publisher.org",))

    seen = []

    def fake_load(domain, profile_name):
        seen.append((domain, profile_name))
        return []

    monkeypatch.setattr(cookie_loader, "_load_domain_cookie_dicts", fake_load)

    cookies = cookie_loader.load_chrome_cookies(profile_name="Default")

    assert cookies == []
    assert seen == [
        ("publisher.org", "Default"),
        ("proxy.example.edu", "Default"),
        ("login.example.edu", "Default"),
    ]
