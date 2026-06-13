"""T7: per-paper supplement cap is env-overridable and defaults higher than 12."""

import importlib


def test_default_max_supplements_env_override(monkeypatch):
    mod = importlib.import_module("harvesting.paywall_context_enrichment")
    monkeypatch.delenv("GVF_MAX_SUPPLEMENTS", raising=False)
    assert mod._default_max_supplements() == 40  # raised from the old 12

    monkeypatch.setenv("GVF_MAX_SUPPLEMENTS", "100")
    assert mod._default_max_supplements() == 100

    # junk / non-positive falls back to the default
    monkeypatch.setenv("GVF_MAX_SUPPLEMENTS", "0")
    assert mod._default_max_supplements() == 40
    monkeypatch.setenv("GVF_MAX_SUPPLEMENTS", "notanint")
    assert mod._default_max_supplements() == 40
