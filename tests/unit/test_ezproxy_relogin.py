"""EZproxy relogin: the Netscape writer, cookie merge, refresh flow, and self-heal.

The refresh exists because the short-lived EZproxy session cookie dies while
the long-lived VUMC SSO session stays valid (TASKS.md design, 2026-08-10).
These tests pin the pieces that must not drift:

- the cookie WRITER must round-trip through the pipeline's own Netscape
  reader (``cookie_loader._cookies_from_file``), including the ``#HttpOnly_``
  prefix that stdlib parsers drop — the EZproxy session cookie is HttpOnly;
- merging must never drop unrelated entries: once ``GVF_COOKIE_FILE`` exists
  it is the ONLY cookie source (the loader does not fall back to Chrome on an
  empty match), so a destructive write would revoke publisher sessions;
- identity-provider tokens stay in the browser profile by default;
- gvf-run's self-heal must respect the offline local-data guard and only fire
  on failure classes a new cookie can fix.
"""

from __future__ import annotations

import stat
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

import pytest

import cli.gvf_run as gvf_run
import cli.institutional_preflight as preflight
import scripts.ezproxy_relogin as relogin
from harvesting.browser_html.cookie_loader import _cookies_from_file

EZPROXY_ENV_VARS = (
    "GVF_EZPROXY_PREFIX",
    "GVF_EZPROXY_HOST",
    "PROXY_LOGIN_PREFIX",
    "PROXY_HOST",
    "GVF_EZPROXY_BASE",
    "PROXY_BASE",
)


@pytest.fixture(autouse=True)
def _clean_ezproxy_env(monkeypatch):
    """EZproxy env vars are not in the suite's hermetic strip list."""
    for var in EZPROXY_ENV_VARS + ("GVF_EZPROXY_AUTOHEAL", "GVF_PREFLIGHT_DOI"):
        monkeypatch.delenv(var, raising=False)


def _configure_proxy(monkeypatch):
    monkeypatch.setenv("GVF_EZPROXY_HOST", "login.proxy.library.vanderbilt.edu")


SESSION_COOKIE = {
    "name": "ezproxy",
    "value": "s3ssion",
    "domain": ".proxy.library.vanderbilt.edu",
    "path": "/",
    "secure": True,
    "httpOnly": True,
    "expires": -1,
}


# ---------------------------------------------------------------------------
# Netscape writer round-trip
# ---------------------------------------------------------------------------


def test_netscape_lines_round_trip_through_the_pipeline_reader(tmp_path):
    cookies = [
        SESSION_COOKIE,
        {
            "name": "_c",
            "value": "abc\tdef",  # tab inside value must survive
            "domain": "onlinelibrary-wiley-com.proxy.library.vanderbilt.edu",
            "path": "/",
            "secure": False,
            "httpOnly": False,
            "expires": 1812463329.5,
        },
    ]
    out = tmp_path / "cookies.txt"
    relogin.write_cookie_file(out, cookies)

    text = out.read_text()
    assert text.startswith("# Netscape HTTP Cookie File")
    assert "#HttpOnly_.proxy.library.vanderbilt.edu\t" in text

    parsed = {c["name"]: c for c in _cookies_from_file(str(out))}
    assert parsed["ezproxy"]["httpOnly"] is True
    assert parsed["ezproxy"]["expires"] == -1  # session cookie: 0 on disk, -1 parsed
    assert parsed["ezproxy"]["secure"] is True
    assert parsed["_c"]["value"] == "abc\tdef"
    assert parsed["_c"]["expires"] == 1812463329


def test_netscape_writer_neutralizes_newlines_and_skips_nameless():
    lines = relogin.playwright_cookies_to_netscape_lines(
        [
            {"name": "a", "value": "evil\ninjected", "domain": "x.com", "path": "/"},
            {"name": "", "value": "nameless", "domain": "x.com", "path": "/"},
        ]
    )
    assert len(lines) == 1
    assert "\n" not in lines[0]
    assert "evilinjected" in lines[0]


def test_write_cookie_file_is_owner_only(tmp_path):
    out = tmp_path / "cookies.txt"
    relogin.write_cookie_file(out, [SESSION_COOKIE])
    mode = stat.S_IMODE(out.stat().st_mode)
    assert mode == 0o600


# ---------------------------------------------------------------------------
# Merge + relevance filtering
# ---------------------------------------------------------------------------


def test_merge_prefers_fresh_and_keeps_unrelated():
    existing = [
        {**SESSION_COOKIE, "value": "stale"},
        {
            "name": "wiley_sess",
            "domain": "onlinelibrary.wiley.com",
            "path": "/",
            "value": "keepme",
        },
    ]
    fresh = [{**SESSION_COOKIE, "value": "minty"}]
    merged = {
        (c["name"], c["domain"]): c for c in relogin.merge_cookies(existing, fresh)
    }
    assert merged[("ezproxy", ".proxy.library.vanderbilt.edu")]["value"] == "minty"
    assert merged[("wiley_sess", "onlinelibrary.wiley.com")]["value"] == "keepme"


def test_select_relevant_keeps_proxy_cookies_and_drops_idp_tokens(monkeypatch):
    _configure_proxy(monkeypatch)
    cookies = [
        SESSION_COOKIE,
        {
            "name": "_c",
            "domain": "onlinelibrary-wiley-com.proxy.library.vanderbilt.edu",
            "path": "/",
        },
        {"name": "ESTSAUTH", "domain": ".login.microsoftonline.com", "path": "/"},
        {"name": "ezproxyl", "domain": "somewhere.else.edu", "path": "/"},
    ]
    kept = {c["name"] for c in relogin.select_relevant_cookies(cookies)}
    assert kept == {"ezproxy", "_c", "ezproxyl"}  # MS SSO token stays in the profile
    all_kept = {
        c["name"] for c in relogin.select_relevant_cookies(cookies, all_domains=True)
    }
    assert "ESTSAUTH" in all_kept


# ---------------------------------------------------------------------------
# refresh_session against a fake browser pool
# ---------------------------------------------------------------------------


class _FakePage:
    def __init__(self, final_url, pool):
        self.url = final_url
        self.goto_calls = []
        self._pool = pool

    def goto(self, url, **kwargs):
        self.goto_calls.append(url)
        self._pool.navigated = True


class _FakePool:
    """Stands in for AuthenticatedBrowserPool; records construction kwargs.

    ``cookies_before`` models the persistent profile's retained (possibly
    stale) cookies; ``cookies_after`` what the context holds once navigation
    has happened.
    """

    last = None

    def __init__(
        self, *, cookies_before=None, cookies_after=None, final_url="", **kwargs
    ):
        self.kwargs = kwargs
        self._cookies_before = cookies_before or []
        self._cookies_after = cookies_after or []
        self.navigated = False
        self._page = _FakePage(final_url, self)
        self.closed = False
        _FakePool.last = self

    def start(self):
        pass

    @contextmanager
    def page(self):
        yield self._page

    def context_cookies(self):
        return list(self._cookies_after if self.navigated else self._cookies_before)

    def close(self):
        self.closed = True


def test_refresh_session_success_dumps_cookies(tmp_path, monkeypatch):
    _configure_proxy(monkeypatch)
    monkeypatch.setattr(
        relogin,
        "AuthenticatedBrowserPool",
        lambda **kw: _FakePool(
            cookies_after=[SESSION_COOKIE],
            final_url="https://login.proxy.library.vanderbilt.edu/menu",
            **kw,
        ),
    )
    result = relogin.refresh_session(tmp_path / "profile", headed=False, wait_s=5)
    assert result.ok is True
    assert result.needs_bootstrap is False
    assert any(c["name"] == "ezproxy" for c in result.cookies)
    pool = _FakePool.last
    assert pool.closed is True
    assert pool.kwargs["headless"] is True
    # the login?url= starting-point form, not the host-rewritten fetch form
    assert pool._page.goto_calls[0].startswith(
        "https://login.proxy.library.vanderbilt.edu/login?url="
    )
    # the profile dir was created private: it holds live SSO tokens
    assert stat.S_IMODE((tmp_path / "profile").stat().st_mode) == 0o700


def test_refresh_session_parked_at_idp_needs_bootstrap(tmp_path, monkeypatch):
    _configure_proxy(monkeypatch)
    monkeypatch.setattr(relogin.time, "sleep", lambda s: None)
    monkeypatch.setattr(
        relogin,
        "AuthenticatedBrowserPool",
        lambda **kw: _FakePool(
            cookies_after=[],
            final_url="https://login.microsoftonline.com/common/oauth2/authorize",
            **kw,
        ),
    )
    result = relogin.refresh_session(tmp_path / "profile", headed=False, wait_s=1)
    assert result.ok is False
    assert result.needs_bootstrap is True
    assert "login/identity-provider" in result.detail


def test_refresh_session_stale_cookie_alone_is_not_success(tmp_path, monkeypatch):
    """The profile retains the previous dead session cookie; parked at the IdP
    with only that stale cookie must NOT read as success."""
    _configure_proxy(monkeypatch)
    monkeypatch.setattr(relogin.time, "sleep", lambda s: None)
    monkeypatch.setattr(
        relogin,
        "AuthenticatedBrowserPool",
        lambda **kw: _FakePool(
            cookies_before=[SESSION_COOKIE],
            cookies_after=[SESSION_COOKIE],  # same value: nothing minted
            final_url="https://login.microsoftonline.com/common/oauth2/authorize",
            **kw,
        ),
    )
    result = relogin.refresh_session(tmp_path / "profile", headed=False, wait_s=1)
    assert result.ok is False
    assert result.needs_bootstrap is True


def test_refresh_session_minted_value_wins_even_off_proxy(tmp_path, monkeypatch):
    _configure_proxy(monkeypatch)
    fresh = {**SESSION_COOKIE, "value": "minty-new"}
    monkeypatch.setattr(
        relogin,
        "AuthenticatedBrowserPool",
        lambda **kw: _FakePool(
            cookies_before=[SESSION_COOKIE],
            cookies_after=[fresh],
            final_url="https://onlinelibrary-wiley-com.proxy.library.vanderbilt.edu/doi/full/x",
            **kw,
        ),
    )
    result = relogin.refresh_session(tmp_path / "profile", headed=False, wait_s=5)
    assert result.ok is True
    assert "minted" in result.detail


def test_refresh_session_proxy_login_page_is_not_success(tmp_path, monkeypatch):
    """The expired-session landing page IS a proxy-host URL (an HTTP-200 SAML
    handoff on /login). With the profile's retained stale cookie present, this
    exact shape used to read as success and close the browser before the SSO
    chain ran (adversarial-review finding)."""
    _configure_proxy(monkeypatch)
    monkeypatch.setattr(relogin.time, "sleep", lambda s: None)
    monkeypatch.setattr(
        relogin,
        "AuthenticatedBrowserPool",
        lambda **kw: _FakePool(
            cookies_before=[SESSION_COOKIE],
            cookies_after=[SESSION_COOKIE],
            final_url="https://login.proxy.library.vanderbilt.edu/login?auth=shibboleth",
            **kw,
        ),
    )
    result = relogin.refresh_session(tmp_path / "profile", headed=False, wait_s=1)
    assert result.ok is False
    assert result.needs_bootstrap is True


def test_refresh_session_preauth_infra_cookie_is_not_a_session(tmp_path, monkeypatch):
    """A pre-auth cookie minted on a proxy host (Shibboleth state, LB affinity)
    must not count as the session while the flow is parked at the IdP."""
    _configure_proxy(monkeypatch)
    monkeypatch.setattr(relogin.time, "sleep", lambda s: None)
    infra = {
        "name": "_shibstate_1",
        "value": "xyz",
        "domain": "login.proxy.library.vanderbilt.edu",
        "path": "/",
    }
    monkeypatch.setattr(
        relogin,
        "AuthenticatedBrowserPool",
        lambda **kw: _FakePool(
            cookies_before=[],
            cookies_after=[infra],
            final_url="https://login.microsoftonline.com/common/oauth2/authorize",
            **kw,
        ),
    )
    result = relogin.refresh_session(tmp_path / "profile", headed=False, wait_s=1)
    assert result.ok is False
    assert result.needs_bootstrap is True


def test_refresh_session_unconfigured_fails_cleanly(tmp_path):
    result = relogin.refresh_session(tmp_path / "profile", wait_s=1)
    assert result.ok is False
    assert "not configured" in result.detail


def test_effective_wait_defaults():
    assert relogin.effective_wait_s(True, None) == relogin.DEFAULT_BOOTSTRAP_WAIT_S
    assert relogin.effective_wait_s(False, None) == relogin.DEFAULT_REFRESH_WAIT_S
    assert relogin.effective_wait_s(True, 45) == 45


# ---------------------------------------------------------------------------
# main(): exit codes and the cookie-file write
# ---------------------------------------------------------------------------


def _viable_report(viable=True, reason=""):
    return preflight.AccessReport(
        viable=viable, should_block=not viable, reason=reason, lines=[]
    )


def test_main_without_proxy_config_exits_2(tmp_path, capsys):
    rc = relogin.main(["--cookie-file", str(tmp_path / "c.txt")])
    assert rc == 2
    assert "not configured" in capsys.readouterr().err


def test_main_without_cookie_file_exits_2(monkeypatch, capsys):
    _configure_proxy(monkeypatch)
    monkeypatch.delenv("GVF_COOKIE_FILE", raising=False)
    rc = relogin.main([])
    assert rc == 2
    assert "GVF_COOKIE_FILE" in capsys.readouterr().err


def test_main_unbootstrapped_profile_exits_3(tmp_path, monkeypatch, capsys):
    _configure_proxy(monkeypatch)
    rc = relogin.main(
        [
            "--cookie-file",
            str(tmp_path / "c.txt"),
            "--profile-dir",
            str(tmp_path / "never_used"),
        ]
    )
    assert rc == 3
    assert "--bootstrap" in capsys.readouterr().err


def _bootstrapped_profile(tmp_path) -> Path:
    profile = tmp_path / "profile"
    profile.mkdir()
    (profile / "Default").mkdir()
    return profile


def test_main_success_merges_into_cookie_file_and_exits_0(tmp_path, monkeypatch):
    _configure_proxy(monkeypatch)
    profile = _bootstrapped_profile(tmp_path)
    cookie_file = tmp_path / "c.txt"
    relogin.write_cookie_file(
        cookie_file,
        [
            {
                "name": "wiley_sess",
                "value": "keepme",
                "domain": "onlinelibrary.wiley.com",
                "path": "/",
            }
        ],
    )
    monkeypatch.setattr(
        relogin,
        "refresh_session",
        lambda *a, **k: relogin.RefreshResult(ok=True, cookies=[SESSION_COOKIE]),
    )
    monkeypatch.setattr(
        preflight, "probe_institutional_access", lambda **k: _viable_report(True)
    )
    rc = relogin.main(
        ["--cookie-file", str(cookie_file), "--profile-dir", str(profile)]
    )
    assert rc == 0
    names = {c["name"] for c in _cookies_from_file(str(cookie_file))}
    assert names == {"wiley_sess", "ezproxy"}  # merged, not replaced


def test_main_probe_still_blocked_exits_4(tmp_path, monkeypatch):
    _configure_proxy(monkeypatch)
    profile = _bootstrapped_profile(tmp_path)
    monkeypatch.setattr(
        relogin,
        "refresh_session",
        lambda *a, **k: relogin.RefreshResult(ok=True, cookies=[SESSION_COOKIE]),
    )
    monkeypatch.setattr(
        preflight,
        "probe_institutional_access",
        lambda **k: _viable_report(False, reason="still walled"),
    )
    rc = relogin.main(
        ["--cookie-file", str(tmp_path / "c.txt"), "--profile-dir", str(profile)]
    )
    assert rc == 4


def test_main_bootstrap_passes_the_profile_gate_and_runs_headed(tmp_path, monkeypatch):
    """--bootstrap on a NEVER-bootstrapped profile must not bounce off the
    'run --bootstrap first' gate, and must launch the browser headed — the
    one flow that cannot be exercised live in CI (adversarial-review finding)."""
    _configure_proxy(monkeypatch)
    seen = {}

    def _fake_refresh(profile_dir, *, headed, wait_s=None, **kw):
        seen["headed"] = headed
        return relogin.RefreshResult(ok=True, cookies=[SESSION_COOKIE])

    monkeypatch.setattr(relogin, "refresh_session", _fake_refresh)
    monkeypatch.setattr(
        preflight, "probe_institutional_access", lambda **k: _viable_report(True)
    )
    rc = relogin.main(
        [
            "--bootstrap",
            "--cookie-file",
            str(tmp_path / "c.txt"),
            "--profile-dir",
            str(tmp_path / "brand_new_profile"),
        ]
    )
    assert rc == 0
    assert seen["headed"] is True


def test_main_bootstrap_timeout_exits_1_not_3(tmp_path, monkeypatch, capsys):
    """A bootstrap that still ends needing interaction must not tell the
    operator to run --bootstrap (they just did) — it is a plain failure."""
    _configure_proxy(monkeypatch)
    monkeypatch.setattr(
        relogin,
        "refresh_session",
        lambda *a, **k: relogin.RefreshResult(
            ok=False, needs_bootstrap=True, detail="timed out at MFA"
        ),
    )
    rc = relogin.main(
        [
            "--bootstrap",
            "--cookie-file",
            str(tmp_path / "c.txt"),
            "--profile-dir",
            str(tmp_path / "brand_new_profile"),
        ]
    )
    assert rc == 1


def test_main_headless_needing_interaction_exits_3(tmp_path, monkeypatch, capsys):
    _configure_proxy(monkeypatch)
    profile = _bootstrapped_profile(tmp_path)
    monkeypatch.setattr(
        relogin,
        "refresh_session",
        lambda *a, **k: relogin.RefreshResult(
            ok=False, needs_bootstrap=True, detail="parked"
        ),
    )
    rc = relogin.main(
        ["--cookie-file", str(tmp_path / "c.txt"), "--profile-dir", str(profile)]
    )
    assert rc == 3
    assert "--bootstrap" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# The probe's structured login_redirect field
# ---------------------------------------------------------------------------


class _FakeResponse:
    def __init__(self, url):
        self.url = url
        self.status_code = 200
        self.text = "x" * 5000
        self.headers = {}


def test_probe_sets_login_redirect_field(tmp_path, monkeypatch):
    _configure_proxy(monkeypatch)
    cookie_file = tmp_path / "c.txt"
    relogin.write_cookie_file(cookie_file, [SESSION_COOKIE])
    monkeypatch.setenv("GVF_COOKIE_FILE", str(cookie_file))

    import scripts.fetch_paywalled as fetch_paywalled

    class _FakeSession:
        def get(self, url, **kwargs):
            return _FakeResponse(
                "https://login.proxy.library.vanderbilt.edu/login?auth=shibboleth"
            )

    monkeypatch.setattr(fetch_paywalled, "make_session", lambda: _FakeSession())
    monkeypatch.setattr(
        fetch_paywalled, "hydrate_session_with_browser_cookies", lambda s, c: None
    )
    report = preflight.probe_institutional_access(timeout_s=1)
    assert report.viable is False
    assert report.login_redirect is True
    assert report.should_block is True


# ---------------------------------------------------------------------------
# gvf-run self-heal
# ---------------------------------------------------------------------------


def _blocked_report(**overrides):
    kwargs = dict(
        viable=False,
        should_block=True,
        ez_cookies=2,
        login_redirect=True,
        reason="expired",
    )
    kwargs.update(overrides)
    return preflight.AccessReport(**kwargs)


def test_self_heal_is_inert_under_the_local_data_guard(monkeypatch):
    monkeypatch.setenv("GVF_DISABLE_LOCAL_DATA", "1")

    def _boom(*a, **k):  # pragma: no cover - must never run
        raise AssertionError("self-heal must not touch subprocess under the guard")

    monkeypatch.setattr(gvf_run.subprocess, "run", _boom)
    assert gvf_run._attempt_ezproxy_self_heal(_blocked_report()) is None


def test_self_heal_skips_unhealable_failures(monkeypatch):
    monkeypatch.delenv("GVF_DISABLE_LOCAL_DATA", raising=False)
    report = _blocked_report(login_redirect=False, ez_cookies=3, reason="cloudflare")
    assert gvf_run._attempt_ezproxy_self_heal(report) is None


def test_self_heal_respects_opt_out(monkeypatch):
    monkeypatch.delenv("GVF_DISABLE_LOCAL_DATA", raising=False)
    monkeypatch.setenv("GVF_EZPROXY_AUTOHEAL", "0")
    assert gvf_run._attempt_ezproxy_self_heal(_blocked_report()) is None


def test_self_heal_requires_a_bootstrapped_profile(tmp_path, monkeypatch):
    monkeypatch.delenv("GVF_DISABLE_LOCAL_DATA", raising=False)
    monkeypatch.setenv("GVF_EZPROXY_PROFILE_DIR", str(tmp_path / "empty"))

    def _boom(*a, **k):  # pragma: no cover
        raise AssertionError("must not launch a refresh without a profile")

    monkeypatch.setattr(gvf_run.subprocess, "run", _boom)
    assert gvf_run._attempt_ezproxy_self_heal(_blocked_report()) is None


def test_self_heal_refreshes_and_reprobes(tmp_path, monkeypatch):
    monkeypatch.delenv("GVF_DISABLE_LOCAL_DATA", raising=False)
    profile = tmp_path / "profile"
    profile.mkdir()
    (profile / "Default").mkdir()
    monkeypatch.setenv("GVF_EZPROXY_PROFILE_DIR", str(profile))

    calls = {}

    def _fake_run(cmd, **kwargs):
        calls["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")

    monkeypatch.setattr(gvf_run.subprocess, "run", _fake_run)
    healed_report = _viable_report(True)
    monkeypatch.setattr(
        preflight, "probe_institutional_access", lambda **k: healed_report
    )

    healed = gvf_run._attempt_ezproxy_self_heal(_blocked_report())
    assert healed is healed_report
    assert calls["cmd"][0] == sys.executable
    assert any(str(part).endswith("ezproxy_relogin.py") for part in calls["cmd"])
    assert "--skip-verify" in calls["cmd"]


def test_self_heal_keeps_original_on_nonzero_exit(tmp_path, monkeypatch):
    monkeypatch.delenv("GVF_DISABLE_LOCAL_DATA", raising=False)
    profile = tmp_path / "profile"
    profile.mkdir()
    (profile / "Default").mkdir()
    monkeypatch.setenv("GVF_EZPROXY_PROFILE_DIR", str(profile))
    monkeypatch.setattr(
        gvf_run.subprocess,
        "run",
        lambda cmd, **k: subprocess.CompletedProcess(
            cmd, 3, stdout="", stderr="needs mfa"
        ),
    )
    assert gvf_run._attempt_ezproxy_self_heal(_blocked_report()) is None
