"""Tests for conditional Tier 3 ensemble adjudication."""

from pipeline.claim_verifier import (
    VariantClaimCard,
    VariantClaimVerifier,
    build_claim_card,
    build_claim_verification_prompt,
    normalize_verification,
)
from pipeline.extraction import ExpertExtractor
from utils.models import Paper


def _variant(name: str, total: int = 43, affected: int = 28, unaffected: int = 15):
    return {
        "gene_symbol": "RYR2",
        "protein_notation": name,
        "patients": {"count": total, "phenotype": "arrhythmia"},
        "penetrance_data": {
            "total_carriers_observed": total,
            "affected_count": affected,
            "unaffected_count": unaffected,
        },
        "source_location": "Results",
        "additional_notes": "LLM extraction from cohort summary",
    }


def test_risk_flags_repeated_large_non_deterministic_counts():
    extractor = ExpertExtractor(models=["primary-model"], tier_threshold=0)
    data = {
        "variants": [_variant(f"p.Arg{i}Trp") for i in range(1, 5)],
        "extraction_metadata": {},
    }

    risk = extractor._assess_extraction_risk(
        extracted_data=data,
        source_text="RYR2 mutation carriers were described. " * 20,
        estimated_variants=4,
        scanner_variant_count=4,
    )

    assert risk["score"] >= 2
    assert any(
        reason.startswith("repeated_large_count_tuple") for reason in risk["reasons"]
    )
    assert risk["requires_adjudication"] is True


def test_evidence_packet_is_compact_and_contains_count_lines():
    extractor = ExpertExtractor(models=["primary-model"], tier_threshold=0)
    extractor.evidence_packet_max_chars = 2200
    data = {
        "variants": [_variant("p.Arg420Trp")],
        "extraction_metadata": {},
    }
    paper = Paper(pmid="1", title="RYR2 families", gene_symbol="RYR2")
    source_text = "\n".join(
        [
            "Background line not useful.",
            "Among 43 RYR2 mutation carriers, 28 were affected and 15 were unaffected.",
            "The R420W variant was detected in family 125.",
        ]
        + ["irrelevant filler"] * 500
    )
    risk = {
        "score": 2,
        "reasons": ["repeated_large_count_tuple"],
        "source_blockers": [],
    }

    prompt = extractor._build_adjudication_prompt(
        paper=paper,
        primary_model="azure_ai/grok-4-20-reasoning",
        extracted_data=data,
        source_text=source_text,
        risk=risk,
    )

    assert len(prompt) <= extractor.evidence_packet_max_chars
    assert "43 RYR2 mutation carriers" in prompt
    assert "p.Arg420Trp" in prompt


def test_high_risk_extraction_uses_adjudicator_model(monkeypatch):
    captured = {}

    def fake_verify(self, card):
        captured["model"] = self.model
        captured.setdefault("cards", []).append(card)
        return {
            "verdict": "ambiguous",
            "field_verdicts": {
                "variant": "directly_supported",
                "total_carriers": "ambiguous",
                "affected": "ambiguous",
                "unaffected": "ambiguous",
            },
            "corrected_values": {
                "total_carriers": None,
                "affected": None,
                "unaffected": None,
            },
            "reason": "Aggregate count only; variant-level count not supported.",
            "evidence_quote": "Among 43 RYR2 mutation carriers...",
        }

    monkeypatch.setattr(VariantClaimVerifier, "verify", fake_verify)

    extractor = ExpertExtractor(
        models=["azure_ai/grok-4-20-reasoning"], tier_threshold=0
    )
    extractor.adjudicator_models = ["anthropic/claude-sonnet-4-6"]
    extractor.adjudication_risk_threshold = 1
    extractor.evidence_packet_max_chars = 4000
    data = {
        "variants": [_variant(f"p.Arg{i}Trp") for i in range(1, 5)],
        "extraction_metadata": {},
    }

    adjudicated = extractor._maybe_adjudicate_extraction(
        paper=Paper(
            pmid="12106942",
            title="RYR2 families",
            gene_symbol="RYR2",
            disease="catecholaminergic polymorphic ventricular tachycardia",
        ),
        primary_model="azure_ai/grok-4-20-reasoning",
        extracted_data=data,
        source_text=(
            "p.Arg1Trp and p.Arg2Trp were observed. Among 43 RYR2 mutation "
            "carriers, 28 were affected and 15 were unaffected."
        ),
        estimated_variants=4,
        scanner_variant_count=4,
    )

    assert captured["model"] == "anthropic/claude-sonnet-4-6"
    assert captured["cards"][0].variant == "p.Arg1Trp"
    assert (
        captured["cards"][0].disease
        == "catecholaminergic polymorphic ventricular tachycardia"
    )
    assert adjudicated["extraction_metadata"]["claim_verification_applied"] is True
    assert (
        adjudicated["extraction_metadata"]["claim_verification_model"]
        == captured["model"]
    )
    assert (
        adjudicated["variants"][0]["penetrance_data"]["total_carriers_observed"] is None
    )


def test_low_risk_extraction_skips_adjudicator(monkeypatch):
    def fail_if_called(self, prompt):
        raise AssertionError("adjudicator should not run for low-risk extraction")

    monkeypatch.setattr(ExpertExtractor, "call_llm_json", fail_if_called)
    extractor = ExpertExtractor(
        models=["azure_ai/grok-4-20-reasoning"], tier_threshold=0
    )
    extractor.adjudicator_models = ["anthropic/claude-sonnet-4-6"]
    extractor.adjudication_risk_threshold = 2
    data = {
        "variants": [
            {
                "gene_symbol": "RYR2",
                "protein_notation": "p.Arg420Trp",
                "patients": {"count": 1},
                "penetrance_data": {
                    "total_carriers_observed": 1,
                    "affected_count": 1,
                    "unaffected_count": 0,
                },
            }
        ],
        "extraction_metadata": {},
    }

    out = extractor._maybe_adjudicate_extraction(
        paper=Paper(pmid="1", title="Single case", gene_symbol="RYR2"),
        primary_model="azure_ai/grok-4-20-reasoning",
        extracted_data=data,
        source_text="A patient with R420W was affected.",
        estimated_variants=1,
        scanner_variant_count=1,
    )

    assert (
        out["extraction_metadata"]["adjudication_skipped_reason"]
        == "risk_below_threshold"
    )


def test_claim_prompt_distinguishes_disease_affected_from_symptomatic_subset():
    card = VariantClaimCard(
        gene="KCNH2",
        disease="Long QT syndrome",
        pmid="19160088",
        title="Founder mutations",
        variant="p.Arg176Trp",
        extracted={"total_carriers": 112, "affected": 18, "unaffected": 94},
        evidence=(
            "The KCNH2 R176W mutation was identified in 112 LQTS patients, "
            "of which 18 were symptomatic."
        ),
    )

    prompt = build_claim_verification_prompt(card)

    assert "affected=N even if only a smaller subset is symptomatic" in prompt
    assert "Long QT syndrome" in prompt


def test_claim_verification_guard_corrects_symptom_subset_as_affected_count():
    card = VariantClaimCard(
        gene="KCNH2",
        disease="Long QT syndrome",
        pmid="19160088",
        title="Founder mutations",
        variant="p.Arg176Trp",
        extracted={"total_carriers": 112, "affected": 18, "unaffected": 94},
        evidence="KCNH2 R176W mutation in 112 LQTS patients, of which 18 were symptomatic.",
    )
    raw = {
        "verdict": "ambiguous",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "directly_supported",
            "unaffected": "ambiguous",
        },
        "corrected_values": {
            "total_carriers": 112,
            "affected": 18,
            "unaffected": None,
        },
        "reason": "Model treated symptomatic subset as affected.",
        "evidence_quote": "112 LQTS patients, of which 18 were symptomatic",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"]["affected"] == 112
    assert normalized["corrected_values"]["unaffected"] == 0
    assert normalized["field_verdicts"]["unaffected"] == "inferred_supported"


def test_claim_verification_promotes_supported_corrected_values():
    card = VariantClaimCard(
        gene="KCNH2",
        disease="Long QT syndrome",
        pmid="19160088",
        title="Founder mutations",
        variant="p.Arg176Trp",
        extracted={"total_carriers": 112, "affected": 18, "unaffected": 94},
        evidence="KCNH2 R176W mutation in 112 LQTS patients.",
    )
    raw = {
        "verdict": "inferred_supported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "unsupported",
            "unaffected": "unsupported",
        },
        "corrected_values": {
            "total_carriers": 112,
            "affected": 112,
            "unaffected": 0,
        },
        "reason": "Original affected/unaffected were unsupported, but corrected values are supported.",
        "evidence_quote": "112 LQTS patients",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["field_verdicts"]["affected"] == "inferred_supported"
    assert normalized["field_verdicts"]["unaffected"] == "inferred_supported"


def test_claim_verification_infers_missing_count_from_supported_identity():
    raw = {
        "verdict": "inferred_supported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "unsupported",
            "unaffected": "directly_supported",
        },
        "corrected_values": {
            "total_carriers": 2,
            "affected": None,
            "unaffected": 2,
        },
        "reason": "Unselected cohort carriers are unaffected.",
        "evidence_quote": "KCNH2 L552S (n=2)",
    }

    normalized = normalize_verification(raw)

    assert normalized["corrected_values"]["affected"] == 0
    assert normalized["field_verdicts"]["affected"] == "inferred_supported"


def test_claim_card_table_evidence_includes_header_context():
    text = "\n".join(
        [
            "### Table 4",
            "",
            "In silico functional analysis of missense variants",
            "",
            "| Gene | Protein | Polyphen-2 | SIFT | Mutation assessor |",
            "|---|---|---|---|---|",
            "| KCNH2 | p.A561V | 2 | 1 | 2 |",
        ]
    )
    card = build_claim_card(
        source_text=text,
        gene="KCNH2",
        disease="Long QT syndrome",
        pmid="24606995",
        title="Danish LQTS mutations",
        variant={
            "protein_notation": "p.A561V",
            "patients": {"count": 1},
            "penetrance_data": {
                "total_carriers_observed": 1,
                "affected_count": 1,
                "unaffected_count": 0,
            },
        },
    )

    assert card is not None
    assert "Polyphen-2" in card.evidence
    assert "In silico functional analysis" in card.evidence
