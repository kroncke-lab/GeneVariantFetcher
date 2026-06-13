import sys
from types import ModuleType, SimpleNamespace

from harvesting.browser_html.authenticated_pool import AuthenticatedBrowserPool


def test_authenticated_pool_uses_persistent_profile(monkeypatch, tmp_path):
    calls = {}

    class _Context:
        def __init__(self):
            self.cookies = []
            self.init_scripts = []
            self.closed = False

        def add_init_script(self, script):
            self.init_scripts.append(script)

        def add_cookies(self, cookies):
            self.cookies.extend(cookies)

        def close(self):
            self.closed = True

    context = _Context()

    class _Chromium:
        def launch_persistent_context(self, user_data_dir, **kwargs):
            calls["persistent"] = {
                "user_data_dir": user_data_dir,
                "kwargs": kwargs,
            }
            return context

    class _Playwright:
        chromium = _Chromium()

        def stop(self):
            calls["stopped"] = True

    def _sync_playwright():
        return SimpleNamespace(start=lambda: _Playwright())

    playwright_module = ModuleType("playwright")
    playwright_module.__path__ = []
    sync_api_module = ModuleType("playwright.sync_api")
    sync_api_module.sync_playwright = _sync_playwright
    monkeypatch.setitem(sys.modules, "playwright", playwright_module)
    monkeypatch.setitem(sys.modules, "playwright.sync_api", sync_api_module)

    profile_dir = tmp_path / "publisher-profile"
    cookie = {
        "name": "PublisherSession",
        "value": "abc123",
        "domain": ".onlinelibrary.wiley.com",
        "path": "/",
    }
    pool = AuthenticatedBrowserPool(
        cookies=[cookie],
        headless=False,
        persistent_profile_path=str(profile_dir),
    )

    pool.start()
    pool.close()

    assert calls["persistent"]["user_data_dir"] == str(profile_dir)
    assert calls["persistent"]["kwargs"]["headless"] is False
    assert calls["persistent"]["kwargs"]["channel"] == "chrome"
    assert context.cookies == [cookie]
    assert context.init_scripts
    assert context.closed is True
    assert calls["stopped"] is True


def test_use_stealth_forces_headful_chrome_and_skips_init_script(monkeypatch, tmp_path):
    """Stealth backend: patchright provider, real Chrome, forced headful, no
    custom init-script (which would re-introduce a detectable fingerprint)."""
    calls = {}

    class _Context:
        def __init__(self):
            self.init_scripts = []
            self.cookies = []
            self.closed = False

        def add_init_script(self, script):
            self.init_scripts.append(script)

        def add_cookies(self, cookies):
            self.cookies.extend(cookies)

        def close(self):
            self.closed = True

    context = _Context()

    class _Browser:
        def new_context(self, **kwargs):
            calls["context_kwargs"] = kwargs
            return context

        def close(self):
            calls["browser_closed"] = True

    class _Chromium:
        def launch(self, **kwargs):
            calls["launch"] = kwargs
            return _Browser()

    class _Playwright:
        chromium = _Chromium()

        def stop(self):
            calls["stopped"] = True

    pw_mod = ModuleType("patchright")
    pw_mod.__path__ = []
    sync_mod = ModuleType("patchright.sync_api")
    sync_mod.sync_playwright = lambda: SimpleNamespace(start=lambda: _Playwright())
    monkeypatch.setitem(sys.modules, "patchright", pw_mod)
    monkeypatch.setitem(sys.modules, "patchright.sync_api", sync_mod)

    # headless=True must be overridden to headful under stealth
    pool = AuthenticatedBrowserPool(headless=True, use_stealth=True)
    pool.start()
    pool.close()

    assert pool.headless is False  # forced headful
    assert calls["launch"]["channel"] == "chrome"
    assert calls["launch"]["headless"] is False
    # minimal context: no custom user_agent / extra headers under stealth
    assert "user_agent" not in calls["context_kwargs"]
    assert "extra_http_headers" not in calls["context_kwargs"]
    # the hand-rolled stealth init-script is NOT injected under patchright
    assert context.init_scripts == []
