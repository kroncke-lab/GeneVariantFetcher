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
    """Minimal stand-in exposing only count-hygiene settings."""

    def __init__(self, guard="off", classifier="off", classifier_fields="all"):
        self.count_guard_policy = guard
        self.count_classifier_policy = classifier
        self.count_classifier_fields = classifier_fields


@pytest.fixture
def patch_policies(monkeypatch):
    """Return a setter that points steps' get_settings() at a stub."""

    def _set(guard="off", classifier="off", classifier_fields="all"):
        monkeypatch.setattr(
            "config.settings.get_settings",
            lambda: _StubSettings(
                guard=guard,
                classifier=classifier,
                classifier_fields=classifier_fields,
            ),
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


def test_nonhuman_source_guard_clears_clinical_counts():
    data = {
        "variants": [
            {
                "gene_symbol": "MYBPC3",
                "protein_notation": "p.A31P",
                "patients": {"count": 1089},
                "penetrance_data": {
                    "total_carriers_observed": 1089,
                    "affected_count": 344,
                    "unaffected_count": 745,
                },
            }
        ],
        "extraction_metadata": {},
    }
    source_text = (
        "## A cardiac myosin binding protein C mutation in Maine Coon cats\n\n"
        "Animals: 3310 DNA samples from cats were evaluated for feline "
        "hypertrophic cardiomyopathy."
    )

    steps._apply_nonhuman_clinical_count_guard(data, source_text)

    variant = data["variants"][0]
    assert variant["patients"]["count"] is None
    assert variant["penetrance_data"]["total_carriers_observed"] is None
    assert variant["penetrance_data"]["affected_count"] is None
    assert variant["penetrance_data"]["unaffected_count"] is None
    assert variant["nonhuman_source_flags"]["carriers"]["raw"] == 1089
    assert data["extraction_metadata"]["nonhuman_count_guard"]["cleared"] == 3


def test_nonhuman_source_guard_skips_mixed_human_cat_reviews():
    data = {
        "variants": [
            {
                "gene_symbol": "MYBPC3",
                "protein_notation": "p.Arg502Trp",
                "patients": {"count": 25},
                "penetrance_data": {
                    "total_carriers_observed": 25,
                    "affected_count": 15,
                    "unaffected_count": 10,
                },
            }
        ],
        "extraction_metadata": {},
    }
    source_text = (
        "## The Genetic Basis of Hypertrophic Cardiomyopathy in Cats and Humans\n\n"
        "### HUMAN HCM\n\n"
        "Of the 25 individuals that carried this mutation, 15 had HCM.\n\n"
        "### FELINE HCM\n\n"
        "Maine Coon cats and Ragdoll cats carry other MYBPC3 mutations."
    )

    steps._apply_nonhuman_clinical_count_guard(data, source_text)

    variant = data["variants"][0]
    assert variant["patients"]["count"] == 25
    assert variant["penetrance_data"]["affected_count"] == 15
    assert "nonhuman_source_flags" not in variant
    assert "nonhuman_count_guard" not in data["extraction_metadata"]


def test_resolve_count_policies_defaults_off_without_settings(monkeypatch):
    """If Settings can't load, the helper falls back to ('off', 'off')."""

    def _boom():
        raise RuntimeError("no env")

    monkeypatch.setattr("config.settings.get_settings", _boom)
    assert steps._resolve_count_policies() == ("off", "off", None)


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
        assert settings.count_classifier_fields == "all"
    finally:
        get_settings.cache_clear()


def test_settings_count_policy_accepts_flag_clear_and_field_scope(monkeypatch):
    """flag/clear parse and normalize (case-insensitive) on real Settings."""
    from config.settings import get_settings

    _baseline_env(monkeypatch)
    monkeypatch.setenv("COUNT_GUARD_POLICY", "FLAG")
    monkeypatch.setenv("COUNT_CLASSIFIER_POLICY", "clear")
    monkeypatch.setenv("COUNT_CLASSIFIER_FIELDS", "Affected, Unaffected")
    get_settings.cache_clear()
    try:
        settings = get_settings()
        assert settings.count_guard_policy == "flag"
        assert settings.count_classifier_policy == "clear"
        assert settings.count_classifier_fields == "affected,unaffected"
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


def test_settings_count_classifier_fields_rejects_unknown(monkeypatch):
    """An unknown classifier field is rejected by the validator."""
    from config.settings import get_settings

    _baseline_env(monkeypatch)
    monkeypatch.setenv("COUNT_CLASSIFIER_FIELDS", "carriers,status")
    get_settings.cache_clear()
    try:
        with pytest.raises(ValueError):
            get_settings()
    finally:
        get_settings.cache_clear()


def test_classifier_field_scope_clears_only_requested_fields(patch_policies):
    """COUNT_CLASSIFIER_FIELDS can limit clearing to affected-status counts."""
    patch_policies(
        guard="off",
        classifier="clear",
        classifier_fields="affected,unaffected",
    )
    data = {
        "variants": [
            {
                "patients": {"count": 453},
                "penetrance_data": {
                    "total_carriers_observed": 453,
                    "affected_count": 200,
                    "unaffected_count": 253,
                },
                "count_provenance": {
                    "carriers_count_type": "cohort_total",
                    "affected_count_type": "cohort_total",
                    "unaffected_count_type": "control",
                },
            }
        ],
        "extraction_metadata": {},
    }

    steps._apply_count_hygiene(data)

    variant = data["variants"][0]
    assert variant["patients"]["count"] == 453
    assert variant["penetrance_data"]["total_carriers_observed"] == 453
    assert variant["penetrance_data"]["affected_count"] is None
    assert variant["penetrance_data"]["unaffected_count"] is None
    assert set(variant["count_classifier_flags"]) == {"affected", "unaffected"}
    assert data["extraction_metadata"]["count_classifier"]["fields"] == [
        "affected",
        "unaffected",
    ]
