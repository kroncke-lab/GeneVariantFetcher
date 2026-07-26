"""The per-paper final check (Steps 3.8/3.9) is parked -- default OFF.

Parked on 2026-07-26: one @xhigh reasoning call per paper cost more time and
money than its measured effect justified for a step that only RECORDS findings.
These tests exist so the expensive step cannot be re-enabled by accident, and so
the 3.8/3.9 pairing stays intact.
"""

from config.settings import Settings


def test_final_check_and_its_composer_default_off():
    settings = Settings()
    assert settings.paper_final_check_enabled is False
    assert settings.paper_final_check_gate_enabled is False


def test_reviving_the_composer_alone_is_documented_as_a_trap():
    """3.9 without 3.8 fails acceptance on stale findings, unsatisfiably.

    ``cli/gvf_run.py`` raises a stage failure demanding a "reviewer replay" when
    the composer refuses stale findings. With no live reviewer that demand can
    never be met, so the two flags must be revived together. Keep the warning
    attached to the settings a reader will actually consult.
    """
    fields = Settings.model_fields
    gate_doc = fields["paper_final_check_gate_enabled"].description or ""
    check_doc = fields["paper_final_check_enabled"].description or ""
    assert "PARKED" in check_doc and "PARKED" in gate_doc
    assert "PAPER_FINAL_CHECK_GATE_ENABLED" in check_doc
    assert "reviewer replay" in gate_doc


def test_parked_step_keeps_its_implementation_and_reason_codes():
    """Parked means dormant, not deleted -- revival must not need a rewrite."""
    from pipeline.paper_final_check import VALID_FLAG_REASONS
    from pipeline.paper_final_check_gate import ENFORCEABLE_REASON_CODES

    assert "attributed_to_other_study" in VALID_FLAG_REASONS
    assert "attributed_to_other_study" in ENFORCEABLE_REASON_CODES
