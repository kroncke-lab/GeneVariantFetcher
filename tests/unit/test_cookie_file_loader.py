"""Offline tests for GVF_COOKIE_FILE loading — the Keychain-bypass cookie source.

Newer Chrome locks its Safe Storage keychain item to Chrome, so browser_cookie3
fails ("Unable to get key") in a plain Terminal and always fails headless /
through an agent. A cookies.txt exported from inside Chrome is then the reliable
source. These tests are pure filesystem — no browser_cookie3, no Keychain.
"""

from harvesting.browser_html.cookie_loader import (
    _cookies_from_file,
    load_chrome_cookies,
)

# A minimal Netscape cookies.txt: a normal cookie, an #HttpOnly_ cookie (the
# prefix stdlib cookiejar parsers drop), a real comment, and a value with a tab.
_SAMPLE = "\n".join(
    [
        "# Netscape HTTP Cookie File",
        "# a genuine comment line",
        "\t".join(
            [
                "proxy.library.vanderbilt.edu",
                "TRUE",
                "/",
                "TRUE",
                "1999999999",
                "ezproxy",
                "sess-abc",
            ]
        ),
        "#HttpOnly_.proxy.library.vanderbilt.edu\t"
        + "\t".join(["TRUE", "/", "TRUE", "0", "ezproxyl", "ho-xyz"]),
        "\t".join(
            [
                "onlinelibrary.wiley.com",
                "FALSE",
                "/",
                "TRUE",
                "1999999999",
                "SESSION",
                "val\twith-tab",
            ]
        ),
        "",
    ]
)


def test_cookies_from_file_parses_httponly_and_tabs(tmp_path):
    f = tmp_path / "cookies.txt"
    f.write_text(_SAMPLE, encoding="utf-8")
    by = {c["name"]: c for c in _cookies_from_file(str(f))}

    assert set(by) == {"ezproxy", "ezproxyl", "SESSION"}  # comment lines skipped
    assert by["ezproxy"]["domain"] == "proxy.library.vanderbilt.edu"
    assert by["ezproxy"]["value"] == "sess-abc"
    assert by["ezproxy"]["secure"] is True
    assert by["ezproxyl"]["httpOnly"] is True  # #HttpOnly_ prefix honored, not dropped
    assert by["ezproxyl"]["expires"] == -1  # 0 -> session cookie (-1 for Playwright)
    assert by["SESSION"]["value"] == "val\twith-tab"  # tab inside the value tolerated


def test_load_chrome_cookies_uses_file_env(tmp_path, monkeypatch):
    f = tmp_path / "cookies.txt"
    f.write_text(_SAMPLE, encoding="utf-8")
    monkeypatch.setenv("GVF_COOKIE_FILE", str(f))
    # Must return the file cookies without ever touching browser_cookie3/Keychain.
    names = {c["name"] for c in load_chrome_cookies()}
    assert {"ezproxy", "ezproxyl", "SESSION"} <= names


def test_file_cookies_are_recognized_as_ezproxy(tmp_path, monkeypatch):
    from cli.institutional_preflight import _is_ez_cookie

    f = tmp_path / "cookies.txt"
    f.write_text(_SAMPLE, encoding="utf-8")
    monkeypatch.setenv("GVF_COOKIE_FILE", str(f))
    cookies = load_chrome_cookies()
    # The preflight probe counts these -> a file-sourced session unblocks a run.
    assert any(_is_ez_cookie(c) for c in cookies)


def test_missing_cookie_file_does_not_crash(tmp_path, monkeypatch):
    monkeypatch.setenv("GVF_COOKIE_FILE", str(tmp_path / "nope.txt"))
    # Non-existent file: warn + fall back (returns a list, never raises here).
    assert isinstance(load_chrome_cookies(domains=[]), list)
