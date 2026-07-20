"""Pure triage predicate for the per-paper final check."""

from pipeline.paper_final_check_triage import (
    PaperRiskView,
    decide_final_check_tier,
    policy_hash,
    should_escalate_cheap,
)


def _view(**kw) -> PaperRiskView:
    base = dict(
        pmid="1",
        n_count_facts=3,
        any_quarantine=False,
        any_trust_reason=False,
        frac_facts_with_quote=1.0,
        n_distinct_count_roles=1,
        extraction_confidence="medium",
        has_usable_source=True,
        empty_no_source=False,
    )
    base.update(kw)
    return PaperRiskView(**base)


# ---- skip lane ----
def test_empty_no_source_skips():
    d = decide_final_check_tier(
        _view(n_count_facts=0, has_usable_source=False, empty_no_source=True)
    )
    assert d.tier == "skip" and "empty_no_source" in d.reasons


def test_clean_db_only_skips_no_completeness_lane():
    d = decide_final_check_tier(_view(has_usable_source=False), audit_sample_rate=0.0)
    assert d.tier == "skip" and "deterministic_clean_db_only" in d.reasons


# ---- full lane (completeness-critical + hard risk) ----
def test_zero_counts_with_source_is_full():
    d = decide_final_check_tier(_view(n_count_facts=0))
    assert d.tier == "full" and "zero_counts_with_source" in d.reasons


def test_trust_gate_hit_is_full():
    assert decide_final_check_tier(_view(any_quarantine=True)).tier == "full"
    assert decide_final_check_tier(_view(any_trust_reason=True)).tier == "full"


def test_low_or_unknown_confidence_is_full():
    assert decide_final_check_tier(_view(extraction_confidence="low")).tier == "full"
    assert decide_final_check_tier(_view(extraction_confidence=None)).tier == "full"


def test_medium_confidence_is_not_a_risk_signal():
    # "medium" is the hardcoded default on several extraction paths — noise.
    assert (
        decide_final_check_tier(
            _view(extraction_confidence="medium"), audit_sample_rate=0.0
        ).tier
        == "cheap"
    )


def test_thin_provenance_is_full():
    d = decide_final_check_tier(_view(frac_facts_with_quote=0.1))
    assert d.tier == "full" and "thin_provenance" in d.reasons


def test_complex_count_surface_is_full():
    assert decide_final_check_tier(_view(n_count_facts=30)).tier == "full"


def test_mixed_count_roles_is_full():
    assert decide_final_check_tier(_view(n_distinct_count_roles=3)).tier == "full"


# ---- cheap lane ----
def test_clean_with_source_is_cheap_not_skip():
    # The core rule: clean != skip when source exists (completeness preserved).
    d = decide_final_check_tier(_view(), audit_sample_rate=0.0)
    assert d.tier == "cheap" and "clean_with_source" in d.reasons


def test_audit_sample_forces_full_for_a_stable_slice():
    # With rate=1.0 every otherwise-cheap paper is audited at full.
    d = decide_final_check_tier(_view(), audit_sample_rate=1.0)
    assert d.tier == "full" and "audit_sample" in d.reasons
    # ...and rate=0.0 lets it stay cheap (deterministic, no RNG).
    assert decide_final_check_tier(_view(), audit_sample_rate=0.0).tier == "cheap"


# ---- escalation ----
def test_escalate_on_flag_miss_and_low_confidence():
    assert should_escalate_cheap({"verdict": "flag", "confidence": 0.99})[0]
    assert should_escalate_cheap(
        {"verdict": "ok", "completeness": {"n_missing": 2}, "confidence": 0.99}
    )[0]
    assert should_escalate_cheap({"verdict": "ok", "confidence": 0.5})[0]
    assert should_escalate_cheap("not a dict")[0]


def test_no_escalate_on_clean_confident_complete():
    ok, reasons = should_escalate_cheap(
        {
            "verdict": "ok",
            "confidence": 0.97,
            "completeness": {"n_missing": 0, "status": "complete"},
        }
    )
    assert not ok and reasons == ()


def test_policy_hash_stable_and_tagged():
    assert policy_hash().startswith("tri1-")
    assert policy_hash() == policy_hash()
