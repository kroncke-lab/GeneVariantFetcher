"""Tests for BrowserHTMLFetcher orchestration.

We mock the browser pool and strategies so these tests don't require
Playwright. Focus: enable/disable behavior, embargo gating, allowlist
filtering, log writing, and max-per-run cap.
"""

from __future__ import annotations

import csv
import datetime
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from harvesting.browser_html.base import FetchResult, PublisherStrategy
from harvesting.browser_html.fetcher import BrowserHTMLFetcher
from harvesting.browser_html.strategies import _REGISTRY


@pytest.fixture
def settings_enabled():
    return SimpleNamespace(
        enable_browser_html_fallback=True,
        browser_html_publisher_allowlist=["aha", "wiley", "generic"],
        browser_html_min_embargo_months=None,
        browser_html_headless=True,
        browser_html_max_per_run=50,
        browser_html_per_paper_timeout_s=10,
    )


@pytest.fixture
def settings_disabled():
    return SimpleNamespace(
        enable_browser_html_fallback=False,
        browser_html_publisher_allowlist=["generic"],
        browser_html_min_embargo_months=None,
        browser_html_headless=True,
        browser_html_max_per_run=50,
        browser_html_per_paper_timeout_s=10,
    )


@pytest.fixture
def fetcher_args(tmp_path):
    return {
        "scraper": MagicMock(),
        "converter": MagicMock(),
        "session": MagicMock(),
        "output_dir": tmp_path,
    }


def _read_log(log_path: Path) -> list[dict]:
    with open(log_path, newline="") as f:
        return list(csv.DictReader(f))


# ----------------------------------------------------------------------
# Enable/disable behavior
# ----------------------------------------------------------------------


def test_disabled_fetcher_returns_none(settings_disabled, fetcher_args):
    fetcher = BrowserHTMLFetcher(settings=settings_disabled, **fetcher_args)
    assert fetcher.is_enabled() is False
    result = fetcher.fetch(pmid="1", doi="10.1161/test")
    assert result is None


def test_enabled_fetcher_with_no_doi_skips(settings_enabled, fetcher_args):
    fetcher = BrowserHTMLFetcher(settings=settings_enabled, **fetcher_args)
    # No DOI → no strategy matches except generic (always claims via fallback).
    # Generic doesn't gate by DOI but its fetch needs one. We expect it to
    # attempt and the strategy itself returns an error.
    fake_strategy = _FakeStrategy(name="generic", embargo_months=None)
    with _replace_strategies({"generic": _FakeStrategy}), _stub_pool_returning(
        FetchResult(publisher="generic", error="no doi")
    ):
        result = fetcher.fetch(pmid="1", doi="")

    # Even with empty DOI we may attempt — outcome is "failed".
    assert result is None or not result.is_usable()


# ----------------------------------------------------------------------
# Embargo gating
# ----------------------------------------------------------------------


def test_embargo_blocks_recent_aha_paper(settings_enabled, fetcher_args, tmp_path):
    fetcher = BrowserHTMLFetcher(settings=settings_enabled, **fetcher_args)
    six_months_ago = datetime.date.today() - datetime.timedelta(days=180)

    with _replace_strategies({"aha": _FakeAHAStrategy}):
        result = fetcher.fetch(
            pmid="1", doi="10.1161/CIRCRESAHA.118.999", pub_date=six_months_ago
        )

    assert result is None
    rows = _read_log(fetcher.log_path)
    assert len(rows) == 1
    assert rows[0]["outcome"] == "skipped"
    assert "embargo" in rows[0]["reason"]


def test_embargo_passes_old_aha_paper(settings_enabled, fetcher_args):
    fetcher = BrowserHTMLFetcher(settings=settings_enabled, **fetcher_args)
    two_years_ago = datetime.date.today() - datetime.timedelta(days=730)

    expected_result = FetchResult(
        main_markdown="x" * 1000,
        publisher="aha",
        final_url="https://www.ahajournals.org/doi/...",
    )
    with _replace_strategies({"aha": _FakeAHAStrategy}), _stub_pool_returning(
        expected_result
    ):
        result = fetcher.fetch(
            pmid="1", doi="10.1161/CIRCRESAHA.118.999", pub_date=two_years_ago
        )

    assert result is not None
    assert result.is_usable()
    rows = _read_log(fetcher.log_path)
    assert rows[-1]["outcome"] == "success"


# ----------------------------------------------------------------------
# Allowlist filtering
# ----------------------------------------------------------------------


def test_allowlist_skips_non_matching_strategy(settings_enabled, fetcher_args):
    settings_enabled.browser_html_publisher_allowlist = ["wiley"]
    fetcher = BrowserHTMLFetcher(settings=settings_enabled, **fetcher_args)
    two_years_ago = datetime.date.today() - datetime.timedelta(days=730)

    with _replace_strategies({"aha": _FakeAHAStrategy, "wiley": _FakeWileyStrategy}):
        result = fetcher.fetch(
            pmid="1",
            doi="10.1161/CIRCRESAHA.118.999",
            pub_date=two_years_ago,
        )

    # AHA filtered, no other strategy claims, no generic in allowlist → skip.
    assert result is None
    rows = _read_log(fetcher.log_path)
    assert rows[-1]["outcome"] == "skipped"
    assert "no strategy matched" in rows[-1]["reason"]


# ----------------------------------------------------------------------
# Max-per-run cap
# ----------------------------------------------------------------------


def test_max_per_run_caps_attempts(settings_enabled, fetcher_args):
    settings_enabled.browser_html_max_per_run = 2
    fetcher = BrowserHTMLFetcher(settings=settings_enabled, **fetcher_args)
    two_years_ago = datetime.date.today() - datetime.timedelta(days=730)

    expected = FetchResult(main_markdown="x" * 1000, publisher="aha", final_url="u")
    with _replace_strategies({"aha": _FakeAHAStrategy}), _stub_pool_returning(expected):
        for _ in range(3):
            fetcher.fetch(
                pmid="1",
                doi="10.1161/CIRCRESAHA.118.999",
                pub_date=two_years_ago,
            )

    assert fetcher.attempts_made == 2
    rows = _read_log(fetcher.log_path)
    cap_rows = [r for r in rows if "max_per_run" in r["reason"]]
    assert len(cap_rows) == 1


# ----------------------------------------------------------------------
# Embargo override
# ----------------------------------------------------------------------


def test_global_embargo_override_can_widen(settings_enabled, fetcher_args):
    settings_enabled.browser_html_min_embargo_months = 24  # stricter than AHA's 12
    fetcher = BrowserHTMLFetcher(settings=settings_enabled, **fetcher_args)
    pub_date = datetime.date.today() - datetime.timedelta(days=14 * 30)  # 14 mo

    with _replace_strategies({"aha": _FakeAHAStrategy}):
        result = fetcher.fetch(
            pmid="1",
            doi="10.1161/CIRCRESAHA.118.999",
            pub_date=pub_date,
        )

    # 14 mo > AHA's own 12-mo embargo (would pass) but < global 24-mo (blocks).
    assert result is None
    rows = _read_log(fetcher.log_path)
    assert "embargo" in rows[-1]["reason"]


# ======================================================================
# Helpers
# ======================================================================


class _FakeAHAStrategy(PublisherStrategy):
    NAME = "aha"
    DOI_PREFIXES = ("10.1161",)
    DOMAINS = ("ahajournals.org",)
    EMBARGO_MONTHS = 12

    def fetch(self, page, ctx):
        return FetchResult(publisher=self.NAME)


class _FakeWileyStrategy(PublisherStrategy):
    NAME = "wiley"
    DOI_PREFIXES = ("10.1002", "10.1111")
    DOMAINS = ("onlinelibrary.wiley.com",)
    EMBARGO_MONTHS = 0

    def fetch(self, page, ctx):
        return FetchResult(publisher=self.NAME)


class _FakeStrategy(PublisherStrategy):
    NAME = "generic"
    DOI_PREFIXES = ()
    DOMAINS = ()
    EMBARGO_MONTHS = None

    def __init__(self, name="generic", embargo_months=None):
        super().__init__()
        self.NAME = name
        self.EMBARGO_MONTHS = embargo_months

    def fetch(self, page, ctx):
        return FetchResult(publisher=self.NAME)


class _replace_strategies:
    """Context manager that swaps the registry to a fixed mapping.

    Restores both the registry and the autodiscovery flag so subsequent
    tests see the real registry state.
    """

    def __init__(self, mapping):
        self.mapping = mapping
        self.saved = None
        self.saved_loaded = None

    def __enter__(self):
        from harvesting.browser_html import strategies as strategies_mod

        # Force autodiscovery first so we capture the real registry state
        # rather than whatever empty/partial dict happens to exist.
        strategies_mod._autodiscover()
        self.saved = dict(_REGISTRY)
        self.saved_loaded = strategies_mod._LOADED
        _REGISTRY.clear()
        _REGISTRY.update(self.mapping)
        strategies_mod._LOADED = True  # don't re-discover during test

    def __exit__(self, exc_type, exc_val, exc_tb):
        from harvesting.browser_html import strategies as strategies_mod

        _REGISTRY.clear()
        _REGISTRY.update(self.saved)
        strategies_mod._LOADED = self.saved_loaded


class _stub_pool_returning:
    """Context manager that patches BrowserPool.page() to yield a no-op page,
    then patches the strategy's fetch to return a fixed FetchResult.

    We do this by patching BrowserPool to a stub that produces a context
    manager yielding an object the strategy will never inspect, and patching
    each registered strategy class's fetch method to return ``result``.
    """

    def __init__(self, result):
        self.result = result
        self._patchers = []
        self._fetch_patches = []

    def __enter__(self):
        # Patch the BrowserPool used by the fetcher.
        from harvesting.browser_html import fetcher as fetcher_module

        class _StubPool:
            def __init__(self, *a, **kw):
                pass

            def page(self):
                from contextlib import contextmanager

                @contextmanager
                def _cm():
                    page = MagicMock()
                    page.set_default_timeout = MagicMock()
                    yield page

                return _cm()

            def close(self):
                pass

        p = patch.object(fetcher_module, "BrowserPool", _StubPool)
        p.start()
        self._patchers.append(p)

        # Patch every registered strategy's fetch to return our result.
        for cls in _REGISTRY.values():
            patcher = patch.object(cls, "fetch", return_value=self.result)
            patcher.start()
            self._fetch_patches.append(patcher)

    def __exit__(self, exc_type, exc_val, exc_tb):
        for p in self._fetch_patches:
            p.stop()
        for p in self._patchers:
            p.stop()
