"""Tests for debate escalation queue construction."""

from scripts.recall_audit.build_claim_debate_escalation_queue import build_queue


def _verification_record(model: str, value: int, *, wrong: bool = True):
    return {
        "model": model,
        "gene": "RYR2",
        "pmid": "33686871",
        "variant": "p.Gly3118Arg",
        "failure_class": "count_semantics",
        "source_context": "/tmp/33686871_FULL_CONTEXT.md",
        "card": {
            "evidence": (
                "This study involves a family with homozygous G3118R carriers. "
                "Figure 1 pedigree shows affected and unaffected relatives."
            )
        },
        "field_evaluations": [
            {
                "field": "affected",
                "trusted_value": value,
                "after_wrong": wrong,
                "flagged_not_autopopulated": False,
            }
        ],
    }


def test_high_risk_direct_consensus_is_queued_without_debate_disagreement():
    rows = build_queue(
        records=[],
        escalation_model="anthropic/claude-opus-4-7",
        verification_records=[
            _verification_record("azure_ai/grok-4-20-reasoning", 4),
            _verification_record("azure_ai/gpt-5.4", 4),
            _verification_record("azure_ai/DeepSeek-V4-Pro", 4),
        ],
    )

    assert len(rows) == 1
    assert rows[0]["severity"] == "high"
    assert rows[0]["agreement_labels"] == "high_risk_consensus"
    assert rows[0]["disagreement_fields"] == "affected"
    assert "pedigree_family" in rows[0]["reasons"]


def test_low_risk_direct_consensus_is_not_queued():
    rows = build_queue(
        records=[],
        escalation_model="anthropic/claude-opus-4-7",
        verification_records=[
            {
                **_verification_record("azure_ai/grok-4-20-reasoning", 1),
                "card": {"evidence": "One patient carried p.R420W."},
            },
            {
                **_verification_record("azure_ai/gpt-5.4", 1),
                "card": {"evidence": "One patient carried p.R420W."},
            },
        ],
    )

    assert rows == []
