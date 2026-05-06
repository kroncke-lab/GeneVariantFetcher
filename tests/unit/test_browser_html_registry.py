"""Tests for the Tier 3.5 strategy registry.

Verify that strategies auto-discover, that DOI prefixes route correctly,
and that the registry supports allowlist filtering.
"""

from __future__ import annotations

from harvesting.browser_html.strategies import (
    all_strategies,
    find_strategy,
    register,
    registered_names,
)
from harvesting.browser_html.base import FetchResult, PublisherStrategy


def test_known_strategies_register_themselves():
    names = registered_names()
    for required in ("aha", "oxford", "wiley", "elsevier_open", "generic"):
        assert required in names, f"missing strategy: {required}"


def test_aha_doi_matches_aha_strategy():
    s = find_strategy("10.1161/CIRCRESAHA.118.123456")
    assert s is not None
    assert s.NAME == "aha"
    assert s.EMBARGO_MONTHS == 12


def test_oxford_doi_matches_oxford_strategy():
    s = find_strategy("10.1093/europace/euaa067")
    assert s is not None
    assert s.NAME == "oxford"
    assert s.EMBARGO_MONTHS == 12


def test_wiley_doi_matches_wiley_strategy():
    s = find_strategy("10.1002/humu.12345")
    assert s is not None
    assert s.NAME == "wiley"
    # Wiley/Human Mutation are OA — no embargo gate.
    assert s.EMBARGO_MONTHS == 0


def test_elsevier_doi_matches_elsevier_open_strategy():
    s = find_strategy("10.1016/j.cell.2020.01.001")
    assert s is not None
    assert s.NAME == "elsevier_open"
    # No embargo declared — strategy attempts and bails out fast.
    assert s.EMBARGO_MONTHS is None


def test_unknown_doi_falls_back_to_generic():
    s = find_strategy("10.9999/unknown.publisher.123")
    assert s is not None
    assert s.NAME == "generic"


def test_allowlist_excludes_a_strategy():
    s = find_strategy(
        "10.1161/CIRCRESAHA.118.123456",
        allowlist=["oxford", "wiley", "generic"],
    )
    # AHA filtered out → falls through to generic.
    assert s is not None
    assert s.NAME == "generic"


def test_allowlist_excluding_generic_returns_none_for_unknown():
    s = find_strategy(
        "10.9999/unknown",
        allowlist=["aha", "oxford"],
    )
    assert s is None


def test_registry_picks_up_new_strategies_dropped_into_package(monkeypatch):
    """Adding a new strategy class should require zero registry edits.

    We register a synthetic strategy via the @register decorator and
    verify it routes correctly. This proves the "drop a file in
    strategies/" extensibility claim.
    """

    @register
    class _SyntheticPLOS(PublisherStrategy):
        NAME = "test_plos_synthetic"
        DOI_PREFIXES = ("10.1371",)
        DOMAINS = ("journals.plos.org",)
        EMBARGO_MONTHS = 0

        def fetch(self, page, ctx):
            return FetchResult(publisher=self.NAME)

    try:
        s = find_strategy("10.1371/journal.pgen.1009000")
        assert s is not None
        assert s.NAME == "test_plos_synthetic"
        assert "test_plos_synthetic" in registered_names()
    finally:
        # Clean up so the synthetic strategy doesn't leak into other tests.
        from harvesting.browser_html.strategies import _REGISTRY

        _REGISTRY.pop("test_plos_synthetic", None)


def test_all_strategies_returns_instances():
    instances = all_strategies()
    assert all(isinstance(s, PublisherStrategy) for s in instances)
    assert len(instances) >= 5  # at least the five we ship
