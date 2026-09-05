import json

from scripts.recall_audit.fn_root_cause import response_variant_notations


def trace(payload, stage="paper_variant_extraction"):
    return {
        "record_type": "llm_call",
        "context": {"stage": stage},
        "response": {"output_text": json.dumps(payload)},
    }


def test_narrative_mentions_are_not_parser_losses():
    assert (
        response_variant_notations(
            trace({"variants": [], "notes": "R594Q is not current-study evidence"})
        )
        == []
    )
    assert response_variant_notations(trace({"variants": [{"notes": "R594Q"}]})) == []


def test_actual_source_only_rows_are_still_emitted_identities():
    assert response_variant_notations(
        trace({"variants": [{"source_notation": "N1325S"}]})
    ) == ["N1325S"]


def test_router_verifier_and_nonjson_outputs_are_not_extraction_rows():
    assert (
        response_variant_notations(
            trace({"variants": [{"protein_notation": "R594Q"}]}, "claim_verification")
        )
        == []
    )
    assert (
        response_variant_notations(
            {
                "record_type": "llm_call",
                "context": {"stage": "paper_variant_extraction"},
                "response": {"output_text": "R594Q was not extracted"},
            }
        )
        == []
    )
