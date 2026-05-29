"""Tests for the count-guard / count-classifier wiring in pipeline.steps.

These cover the FLAG-ONLY (and off) wiring that pipeline.steps applies to
freshly extracted variants right before they are persisted to per-PMID JSON.

The wiring is exercised through the private helper ``_apply_count_hygiene``,
which is exactly what runs immediately before each ``json.dump`` site in
``extract_variants_step``. Tests are hermetic: synthetic ``extracted_data``
dicts, no network, no DB. The policy is driven by monkeypatching
``config.settings.get_settings`` (the symbol the helper imports lazily) to a
tiny stub so we never depend on the developer's ``.env``.
"""

import copy

import pytest

import pipeline.steps as steps


class _StubSettings:
    """Minimal stand-in exposing only the two count-policy attributes."""

    def __init__(self, guard="off", classifier="off"):
        self.count_guard_policy = guard
        self.count_classifier_policy = classifier


@pytest.fixture
def patch_policies(monkeypatch):
    """Return a setter that points steps' get_settings() at a stub."""

    def _set(guard="off", classifier="off"):
        monkeypatch.setattr(
            "config.settings.get_settings",
            lambda: _StubSettings(guard=guard, classifier=classifier),
        )

    return _set


def _make_extracted_data():
    """Synthetic per-PMID extraction with one study-wide-N count outlier.

    Five small per-variant carrier counts plus one value (453) that is both
    >50 and >10x the per-paper median, matching the canonical KCNQ1 29622001
    G589D 7-vs-453 outlier pattern from the MAE audit.
    """
    variants = [
        {"variant": "p.A1", "patients": {"count": 3}},
        {"variant": "p.A2", "patients": {"count": 5}},
        {"variant": "p.A3", "patients": {"count": 4}},
        {"variant": "p.A4", "patients": {"count": 6}},
        {"variant": "p.G589D", "patients": {"count": 453}},  # study-wide N
    ]
    return {"variants": variants, "extraction_metadata": {"model_used": "stub"}}


def test_flag_policy_annotates_outlier_and_preserves_raw_count(patch_policies):
    """COUNT_GUARD_POLICY=flag annotates the outlier but keeps the raw count."""
    patch_policies(guard="flag", classifier="off")
    data = _make_extracted_data()

    steps._apply_count_hygiene(data)

    outlier = data["variants"][4]
    # Raw count is PRESERVED in flag mode.
    assert outlier["patients"]["count"] == 453
    # The outlier carries a count_outlier_flags annotation recording the raw value.
    assert "count_outlier_flags" in outlier
    assert "carriers" in outlier["count_outlier_flags"]
    flag = outlier["count_outlier_flags"]["carriers"]
    assert flag["raw"] == 453
    assert flag["policy"] == "flag"

    # Non-outlier rows are untouched: no annotation key added.
    for v in data["variants"][:4]:
        assert "count_outlier_flags" not in v

    # No counts were cleared by the classifier (it was off).
    assert "count_classifier_flags" not in outlier

    # Run-level telemetry recorded under extraction_metadata.
    guard_meta = data["extraction_metadata"]["count_outlier_guard"]
    assert guard_meta["policy"] == "flag"
    assert guard_meta["flagged"] == 1
    assert guard_meta["cleared"] == 0


def test_off_policy_is_a_strict_noop(patch_policies):
    """Both knobs off (the default) leaves extracted_data byte-identical."""
    patch_policies(guard="off", classifier="off")
    data = _make_extracted_data()
    before = copy.deepcopy(data)

    steps._apply_count_hygiene(data)

    # Object is unchanged: no annotations, no cleared counts, no telemetry.
    assert data == before
    for v in data["variants"]:
        assert "count_outlier_flags" not in v
        assert "count_classifier_flags" not in v
    assert "count_outlier_guard" not in data["extraction_metadata"]
    assert "count_classifier" not in data["extraction_metadata"]


def test_resolve_count_policies_defaults_off_without_settings(monkeypatch):
    """If Settings can't load, the helper falls back to ('off', 'off')."""

    def _boom():
        raise RuntimeError("no env")

    monkeypatch.setattr("config.settings.get_settings", _boom)
    assert steps._resolve_count_policies() == ("off", "off")


def _baseline_env(monkeypatch):
    """Minimum env so Settings loads independent of the developer's .env."""
    monkeypatch.setenv("NCBI_EMAIL", "user@example.org")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-anthropic")


def test_settings_count_policies_default_off(monkeypatch):
    """The two new knobs default to 'off' on a real Settings instance."""
    from config.settings import get_settings

    _baseline_env(monkeypatch)
    monkeypatch.delenv("COUNT_GUARD_POLICY", raising=False)
    monkeypatch.delenv("COUNT_CLASSIFIER_POLICY", raising=False)
    get_settings.cache_clear()
    try:
        settings = get_settings()
        assert settings.count_guard_policy == "off"
        assert settings.count_classifier_policy == "off"
    finally:
        get_settings.cache_clear()


def test_settings_count_policy_accepts_flag_and_clear(monkeypatch):
    """flag/clear parse and normalize (case-insensitive) on real Settings."""
    from config.settings import get_settings

    _baseline_env(monkeypatch)
    monkeypatch.setenv("COUNT_GUARD_POLICY", "FLAG")
    monkeypatch.setenv("COUNT_CLASSIFIER_POLICY", "clear")
    get_settings.cache_clear()
    try:
        settings = get_settings()
        assert settings.count_guard_policy == "flag"
        assert settings.count_classifier_policy == "clear"
    finally:
        get_settings.cache_clear()


def test_settings_count_policy_rejects_unknown(monkeypatch):
    """An unknown policy value is rejected by the validator."""
    from config.settings import get_settings

    _baseline_env(monkeypatch)
    monkeypatch.setenv("COUNT_GUARD_POLICY", "nonsense")
    get_settings.cache_clear()
    try:
        with pytest.raises(ValueError):
            get_settings()
    finally:
        get_settings.cache_clear()
