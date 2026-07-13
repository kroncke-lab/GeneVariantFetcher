"""Tests for the refresh_recall --land MAE non-regression guard.

The land gate must promote a re-extraction into the canonical DB only when
unique-variant recall holds AND rows-mode count accuracy does not regress.
Recall-only gating let the 2026-06-06 KCNQ1 carriers MAE 0.882->1.274 spike land
silently; these tests pin the MAE/coverage half of the gate and the promotion
decision.
"""

from scripts.refresh_recall import (
    _corpus_bridge_cmd,
    _env_float,
    _mae,
    _mae_n,
    _mae_regressions,
    _promotion_decision,
    _supplement_fetch_cmd,
)


def test_refresh_commands_forward_custom_corpus(tmp_path):
    corpus = tmp_path / "custom-corpus"
    harvest = tmp_path / "harvest"

    fetch = _supplement_fetch_cmd("SCN5A", corpus)
    bridge = _corpus_bridge_cmd("SCN5A", corpus, harvest)

    assert fetch[fetch.index("--corpus") + 1] == str(corpus)
    assert bridge[bridge.index("--corpus") + 1] == str(corpus)
    assert bridge[bridge.index("--out") + 1] == str(harvest)


def _summary(carriers=None, affected=None, unaffected=None, n=1):
    """A minimal run_recall_suite-shaped summary with a rows-mode MAE block.

    Each field value is either a float MAE (matched-row count defaults to ``n``)
    or a ``(mae, n_matched)`` tuple when a test needs to control coverage.
    """
    mae = {}
    for field, val in (
        ("carriers", carriers),
        ("affected", affected),
        ("unaffected", unaffected),
    ):
        if val is None:
            continue
        m, nn = val if isinstance(val, tuple) else (val, n)
        mae[field] = {"sum_abs_error": 0, "n_matched": nn, "mae": m}
    return {"mae": mae}


# --- _mae / _mae_n extractors ----------------------------------------------


def test_mae_reads_carriers_by_default():
    assert _mae(_summary(carriers=0.614)) == 0.614


def test_mae_reads_named_field():
    s = _summary(affected=1.2, unaffected=0.3)
    assert _mae(s, "affected") == 1.2
    assert _mae(s, "unaffected") == 0.3


def test_mae_returns_none_when_unavailable():
    assert _mae(None) is None
    assert _mae({}) is None
    assert _mae(_summary(affected=1.0), "carriers") is None  # block absent
    assert _mae({"mae": {"carriers": {"mae": None}}}) is None  # no matched rows


def test_mae_n_reads_matched_count():
    assert _mae_n(_summary(carriers=(0.6, 120))) == 120
    assert _mae_n(None) == 0
    assert _mae_n(_summary(), "carriers") == 0  # block absent
    assert _mae_n({"mae": {"carriers": {"n_matched": "bad"}}}) == 0


# --- MAE-value regressions --------------------------------------------------


def test_mae_regressions_flags_the_documented_spike():
    # The 2026-06-06 KCNQ1 carriers MAE 0.882 -> 1.274 regression must be caught.
    regs = _mae_regressions(_summary(carriers=0.882), _summary(carriers=1.274))
    assert regs == ["carriers 0.882->1.274"]


def test_mae_regressions_tolerates_small_noise():
    # A re-extraction nudging carriers MAE 0.614 -> 0.620 is within tolerance.
    assert _mae_regressions(_summary(carriers=0.614), _summary(carriers=0.620)) == []


def test_mae_regressions_allows_improvement():
    before = _summary(carriers=0.882, affected=0.74, unaffected=1.2)
    after = _summary(carriers=0.614, affected=0.50, unaffected=0.9)
    assert _mae_regressions(before, after) == []


def test_mae_regressions_checks_all_three_count_fields():
    before = _summary(carriers=0.5, affected=0.5, unaffected=0.5)
    after = _summary(carriers=0.5, affected=0.5, unaffected=2.0)  # only unaffected
    assert _mae_regressions(before, after) == ["unaffected 0.500->2.000"]


def test_mae_regressions_carriers_tolerance_is_more_generous():
    # Carriers gets a modestly more generous absolute tolerance (0.10 default) so
    # recovering real supplement-table rows is not blocked by a small wiggle.
    # A +0.08 carriers increase passes (0.60 -> 0.68, within 0.10)...
    assert _mae_regressions(_summary(carriers=0.60), _summary(carriers=0.68)) == []
    # ...but the same +0.08 increase on unaffected (stricter 0.05 tol) is flagged.
    assert _mae_regressions(_summary(unaffected=0.60), _summary(unaffected=0.68)) == [
        "unaffected 0.600->0.680"
    ]
    # A genuine carriers blow-up (+0.30 > 0.10) is still caught.
    assert _mae_regressions(_summary(carriers=0.60), _summary(carriers=0.90)) == [
        "carriers 0.600->0.900"
    ]


# --- count-coverage regressions (Codex finding: can't game MAE by dropping
#     the very counts it is scored on) ------------------------------------


def test_mae_regressions_flags_count_coverage_collapse():
    # Baseline scored 120 carrier rows; the re-extraction drops them entirely.
    before = _summary(carriers=(0.6, 120))
    after = _summary()  # carriers block gone -> n_matched 0
    assert _mae_regressions(before, after) == ["carriers count-coverage 120->0"]


def test_mae_regressions_tolerates_small_coverage_drop():
    # A modest matched-set shift (100 -> 95, within 10%) is not a regression.
    before = _summary(carriers=(0.6, 100))
    after = _summary(carriers=(0.6, 95))
    assert _mae_regressions(before, after) == []


def test_mae_regressions_ignores_missing_baseline_block():
    # Nothing to regress from when the baseline never scored the field; gaining
    # counts (or both-empty) must not false-block.
    assert _mae_regressions(None, None) == []
    assert _mae_regressions(_summary(), _summary(carriers=99.0)) == []


# --- _promotion_decision (the actual land branch, made unit-testable) -------


def test_promotion_decision_promotes_when_clean():
    promote, msg = _promotion_decision(lm=100, bm=90, mae_regressions=[])
    assert promote is True
    assert msg == ""


def test_promotion_decision_blocks_recall_regression():
    promote, msg = _promotion_decision(lm=80, bm=90, mae_regressions=[])
    assert promote is False
    assert "unique-variant regression" in msg


def test_promotion_decision_blocks_mae_regression_even_when_recall_holds():
    # Recall held (100 >= 90) but a count field regressed -> must NOT promote,
    # so the caller never copies the landed DB over the canonical one.
    promote, msg = _promotion_decision(
        lm=100, bm=90, mae_regressions=["carriers 0.882->1.274"]
    )
    assert promote is False
    assert "carriers 0.882->1.274" in msg


def test_promotion_decision_blocks_variant_row_regression():
    # Unique recall held (100 >= 90) and MAE clean, but matched variant ROWS
    # dropped (4500 < 4600) -> must NOT promote (Codex review guardrail (d)).
    promote, msg = _promotion_decision(
        lm=100, bm=90, mae_regressions=[], lr=4500, br=4600
    )
    assert promote is False
    assert "variant-row regression" in msg


def test_promotion_decision_skips_row_check_when_not_supplied():
    # Back-compat: callers that don't pass row recall (br defaults to 0) skip
    # the row-recall guard entirely.
    promote, msg = _promotion_decision(lm=100, bm=90, mae_regressions=[])
    assert promote is True
    assert msg == ""


# --- _env_float robustness (must not crash the module at import) ------------


def test_env_float_uses_default_when_absent(monkeypatch):
    monkeypatch.delenv("GVF_LAND_MAE_ABS_TOL", raising=False)
    assert _env_float("GVF_LAND_MAE_ABS_TOL", 0.05) == 0.05


def test_env_float_parses_valid_override(monkeypatch):
    monkeypatch.setenv("GVF_LAND_MAE_ABS_TOL", "0.2")
    assert _env_float("GVF_LAND_MAE_ABS_TOL", 0.05) == 0.2


def test_env_float_falls_back_on_garbage(monkeypatch):
    monkeypatch.setenv("GVF_LAND_MAE_ABS_TOL", "foo")
    assert _env_float("GVF_LAND_MAE_ABS_TOL", 0.05) == 0.05
