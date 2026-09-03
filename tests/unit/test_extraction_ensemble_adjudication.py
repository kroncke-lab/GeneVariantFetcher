"""Tests for conditional Tier 3 ensemble adjudication."""

from pipeline.claim_verifier import (
    VariantClaimCard,
    CLAIM_VERIFICATION_SYSTEM_PROMPT,
    FIELD_NAMES,
    VariantClaimVerifier,
    apply_verification_to_variant,
    build_claim_card,
    build_claim_verification_prompt,
    normalize_verification,
    source_verified_claims,
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


def _derived_row_card(*, table_unaffected: int = 62) -> VariantClaimCard:
    return VariantClaimCard(
        gene="RYR2",
        disease="catecholaminergic polymorphic ventricular tachycardia",
        pmid="25814417",
        title="Founder mutation",
        variant="p.Gly357Ser",
        extracted={
            "total_carriers": 185,
            "affected": 97,
            "unaffected": 62,
            "uncertain": 26,
        },
        evidence=(
            "The p.G357S variant was identified in 179 living carriers and "
            "6 genotyped sudden-cardiac-death cases. Supplementary Table 3 "
            "enumerates the living carriers."
        ),
        source_location="Supplementary Table 3 and Results",
        derivation={
            "count_provenance": {
                "carriers_count_type": "per_variant_carrier",
                "affected_column_label": "Previous symptoms OR VA in basal test",
                "affected_count_type": "derived_from_patient_rows",
                "affected_source": "patient_row_phenotype_v2",
                "unaffected_column_label": "Previous symptoms AND VA in basal test",
                "unaffected_count_type": "derived_from_patient_rows",
                "unaffected_source": "patient_row_phenotype_v2",
            },
            "phenotype_derivation": {
                "method": "derived_from_patient_rows",
                "source_table": "Supplementary Table 3",
                "operational_rule": (
                    "affected = symptoms positive OR basal VA yes; unaffected = "
                    "symptoms no AND basal VA no; otherwise uncertain"
                ),
                "complete_table": True,
                "table_total": 179,
                "table_affected": 91,
                "table_unaffected": table_unaffected,
                "table_uncertain": 26,
                "additional_carriers": 6,
                "additional_affected": 6,
                "additional_unaffected": 0,
                "additional_uncertain": 0,
                "predicate_tallies": {
                    "previous_symptoms_positive": 45,
                    "basal_va_positive": 69,
                    "positive_overlap": 23,
                    "dual_negative": 62,
                    "incomplete": 26,
                },
            },
        },
    )


def test_luna_xhigh_claim_verifier_has_reasoning_headroom():
    verifier = VariantClaimVerifier(
        model="azure_ai/gpt-5.6-luna",
        max_tokens=2500,
        reasoning_effort="max",
    )

    assert verifier.max_tokens == 64_000
    assert verifier.reasoning_effort == "max"


def test_claim_verifier_keeps_compact_budget_without_xhigh():
    verifier = VariantClaimVerifier(
        model="azure_ai/gpt-5.6-luna",
        max_tokens=2500,
        reasoning_effort="high",
    )

    assert verifier.max_tokens == 2_500


def test_non_gpt56_xhigh_claim_verifier_keeps_compact_budget():
    verifier = VariantClaimVerifier(
        model="azure_ai/grok-4.5",
        max_tokens=2500,
        reasoning_effort="xhigh",
    )

    assert verifier.max_tokens == 2_500


def test_claim_verification_preserves_typed_family_count_observation():
    variant = {
        "gene_symbol": "BRCA2",
        "cdna_notation": "c.1813insA",
        "patients": {"count": 1},
        "penetrance_data": {"total_carriers_observed": 1},
        "count_provenance": {
            "carriers_column_label": "Number of families",
            "carriers_count_type": "family_count",
        },
    }
    verification = {
        "verdict": "unsupported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "unsupported",
            "affected": "unsupported",
            "unaffected": "unsupported",
        },
        "corrected_values": {
            "total_carriers": None,
            "affected": None,
            "unaffected": None,
        },
        "reason": "One family is not necessarily one carrier.",
        "evidence_quote": "c.1813insA | 1",
    }

    updated, changes = apply_verification_to_variant(variant, verification)

    assert updated["patients"]["count"] == 1
    assert updated["penetrance_data"]["total_carriers_observed"] == 1
    assert updated["claim_verification"] == verification
    assert changes == {}


def test_claim_verification_clears_typed_case_identifier_from_affected_count():
    variant = {
        "gene_symbol": "BRCA2",
        "cdna_notation": "c.3031G>A",
        "patients": {"count": 1},
        "penetrance_data": {
            "total_carriers_observed": 1,
            "affected_count": 6497,
            "unaffected_count": 1,
        },
        "count_provenance": {
            "carriers_count_type": "per_variant_carrier",
            "affected_count_type": "case",
            "unaffected_count_type": "per_variant_carrier",
        },
    }
    verification = {
        "verdict": "inferred_supported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "ambiguous",
            "unaffected": "directly_supported",
        },
        "corrected_values": {
            "total_carriers": 1,
            "affected": None,
            "unaffected": 1,
        },
        "reason": "6497 is the case identifier, not an affected count.",
        "evidence_quote": "case 6497 was a healthy sister",
    }

    updated, changes = apply_verification_to_variant(variant, verification)

    assert updated["penetrance_data"]["affected_count"] is None
    assert changes["affected"] == {
        "old": 6497,
        "new": None,
        "verdict": "ambiguous",
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


def test_claim_verification_ranks_count_risk_before_row_order(monkeypatch):
    captured = {"cards": []}

    def fake_verify(self, card):
        captured["model"] = self.model
        captured["cards"].append(card)
        return {
            "verdict": "directly_supported",
            "field_verdicts": {
                "variant": "directly_supported",
                "total_carriers": "directly_supported",
                "affected": "directly_supported",
                "unaffected": "directly_supported",
            },
            "corrected_values": card.extracted,
            "reason": "Local table row supports the claim.",
            "evidence_quote": card.variant,
        }

    monkeypatch.setattr(VariantClaimVerifier, "verify", fake_verify)

    extractor = ExpertExtractor(models=["primary-model"], tier_threshold=0)
    extractor.max_verifier_cards = 2
    extractor.evidence_packet_max_chars = 6000
    data = {
        "variants": [
            _variant("p.LowRiskAla", total=1, affected=1, unaffected=0),
            {
                **_variant("p.RegexHotArg", total=120, affected=75, unaffected=45),
                "source_location": "Table (regex extraction)",
                "source_layer": "regex_table",
            },
            _variant("p.LowRiskGly", total=2, affected=2, unaffected=0),
            {
                **_variant("p.RegexHotGly", total=80, affected=50, unaffected=30),
                "source_location": "Supplement Table 2 (regex extraction)",
                "source_layer": "regex_table",
            },
        ],
        "extraction_metadata": {},
    }
    source_text = "\n".join(
        [
            "p.LowRiskAla was observed in one affected carrier.",
            "Table row: p.RegexHotArg had 120 carriers, 75 affected, and 45 unaffected.",
            "p.LowRiskGly was observed in two affected carriers.",
            "Table row: p.RegexHotGly had 80 carriers, 50 affected, and 30 unaffected.",
        ]
    )

    verified = extractor._verify_claim_cards_for_extraction(
        paper=Paper(pmid="2", title="Dense table", gene_symbol="RYR2"),
        primary_model="primary-model",
        extracted_data=data,
        source_text=source_text,
        verifier_model="frontier-verifier",
    )

    assert captured["model"] == "frontier-verifier"
    assert [card.variant for card in captured["cards"]] == [
        "p.RegexHotArg",
        "p.RegexHotGly",
    ]
    assert [v["protein_notation"] for v in verified["variants"]] == [
        "p.LowRiskAla",
        "p.RegexHotArg",
        "p.LowRiskGly",
        "p.RegexHotGly",
    ]
    metadata = verified["extraction_metadata"]
    assert metadata["claim_verification_candidate_policy"] == "risk_ranked"
    assert metadata["claim_verification_cards"] == 2
    assert metadata["claim_verification_results"][0]["variant"] == "p.RegexHotArg"


def test_risk_flags_count_bearing_regex_table_rows():
    extractor = ExpertExtractor(models=["primary-model"], tier_threshold=0)
    data = {
        "variants": [
            {
                **_variant("p.RegexHotArg", total=120, affected=75, unaffected=45),
                "source_location": "Table (regex extraction)",
                "source_layer": "regex_table",
            },
            {
                **_variant("p.RegexHotGly", total=80, affected=50, unaffected=30),
                "source_location": "Supplement Table 2 (regex extraction)",
                "source_layer": "regex_table",
            },
            _variant("p.LowRiskAla", total=1, affected=1, unaffected=0),
        ],
        "extraction_metadata": {},
    }

    risk = extractor._assess_extraction_risk(
        extracted_data=data,
        source_text="Table rows mention RYR2 variants and carrier counts.",
        estimated_variants=3,
        scanner_variant_count=3,
    )

    assert risk["score"] >= 2
    assert "count_bearing_high_risk_source_layer:2" in risk["reasons"]
    assert risk["requires_adjudication"] is True


def test_low_risk_extraction_skips_adjudicator(monkeypatch):
    def fail_if_called(self, prompt, **_kw):
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

    # Rules moved to the byte-stable system prompt (cacheable prefix); the
    # per-card user prompt carries only the card and the task recap.
    assert "affected/unaffected split is null" in CLAIM_VERIFICATION_SYSTEM_PROMPT
    assert "not automatically\n  unaffected" in CLAIM_VERIFICATION_SYSTEM_PROMPT
    assert "Long QT syndrome" in prompt
    for field in ("{gene", "{pmid", "{card"):
        assert field not in CLAIM_VERIFICATION_SYSTEM_PROMPT


def test_claim_card_binds_paper_title_phenotype_before_gene_aliases():
    source = """# MAIN TEXT

## NOVEL SCN5A MUTATION ASSOCIATED WITH ARRHYTHMIC STORM DURING ISCHEMIA

One G400A carrier developed an arrhythmic storm during acute ischemia.
Challenge tests for Brugada and long QT syndromes were negative.
"""
    card = build_claim_card(
        source_text=source,
        gene="SCN5A",
        disease="long QT syndrome, Brugada syndrome, cardiac rhythm disorder",
        pmid="example",
        title="Paper 123",
        variant={
            "protein_notation": "p.Gly400Ala",
            "source_notation": "G400A",
            "penetrance_data": {
                "total_carriers_observed": 1,
                "affected_count": 1,
                "unaffected_count": 0,
            },
        },
    )

    assert card is not None
    assert card.paper_title_scope == (
        "NOVEL SCN5A MUTATION ASSOCIATED WITH ARRHYTHMIC STORM DURING ISCHEMIA"
    )
    assert card.paper_target_phenotypes == ["arrhythmic storm"]
    prompt = build_claim_verification_prompt(card)
    assert '"paper_target_phenotypes": [' in prompt
    assert "negative test for syndrome A" in CLAIM_VERIFICATION_SYSTEM_PROMPT


def test_paper_target_guard_rejects_off_target_negative_reversal():
    source = """# MAIN TEXT

## SCN5A MUTATION ASSOCIATED WITH ARRHYTHMIC STORM DURING ISCHEMIA

The G400A mutation carrier developed an arrhythmic storm during acute ischemia.
Challenge tests to unmask Brugada and long QT syndromes were negative.
"""
    card = build_claim_card(
        source_text=source,
        gene="SCN5A",
        disease="long QT syndrome, Brugada syndrome, cardiac rhythm disorder",
        pmid="example",
        title="Paper 123",
        variant={
            "protein_notation": "p.Gly400Ala",
            "source_notation": "G400A",
            "penetrance_data": {
                "total_carriers_observed": 1,
                "affected_count": 1,
                "unaffected_count": 0,
            },
        },
    )
    raw = {
        "verdict": "directly_supported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "directly_supported",
            "unaffected": "directly_supported",
        },
        "corrected_values": {
            "total_carriers": 1,
            "affected": 0,
            "unaffected": 1,
        },
        "reason": "Brugada and LQTS challenge tests were negative.",
        "evidence_quote": "Challenge tests for Brugada and long QT were negative.",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"]["affected"] == 1
    assert normalized["corrected_values"]["unaffected"] is None
    assert normalized["paper_target_scope_overrides"] == ["affected"]
    assert normalized["paper_target_scope"]["phenotypes"] == ["arrhythmic storm"]


def test_paper_target_guard_allows_negative_for_same_target():
    card = VariantClaimCard(
        gene="SCN5A",
        disease="Brugada syndrome",
        pmid="example",
        title="Brugada syndrome challenge testing",
        variant="p.Gly400Ala",
        extracted={
            "total_carriers": 1,
            "affected": 1,
            "unaffected": 0,
        },
        evidence=(
            "The p.Gly400Ala carrier was evaluated for Brugada syndrome. "
            "The Brugada syndrome challenge test was negative."
        ),
        paper_target_phenotypes=["Brugada syndrome"],
        paper_title_scope="Brugada syndrome challenge testing",
    )
    raw = {
        "verdict": "directly_supported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "directly_supported",
            "unaffected": "directly_supported",
        },
        "corrected_values": {
            "total_carriers": 1,
            "affected": 0,
            "unaffected": 1,
        },
        "reason": "The Brugada syndrome challenge was negative.",
        "evidence_quote": "The Brugada syndrome challenge test was negative.",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"]["affected"] == 0
    assert normalized["corrected_values"]["unaffected"] == 1
    assert "paper_target_scope_overrides" not in normalized


def test_verifier_cannot_introduce_unaffected_zero_complement():
    card = VariantClaimCard(
        gene="KCNH2",
        disease="Long QT syndrome",
        pmid="example",
        title="Variant-associated long QT syndrome",
        variant="p.Arg176Trp",
        extracted={"total_carriers": 1, "affected": 0, "unaffected": 1},
        evidence="One carrier was diagnosed with latent LQTS.",
    )
    raw = {
        "verdict": "directly_supported",
        "field_verdicts": {
            "total_carriers": "directly_supported",
            "affected": "directly_supported",
            "unaffected": "directly_supported",
        },
        "corrected_values": {
            "total_carriers": 1,
            "affected": 1,
            "unaffected": 0,
        },
        "reason": "The carrier was diagnosed with latent LQTS.",
        "evidence_quote": "One carrier was diagnosed with latent LQTS.",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"] == {
        "total_carriers": 1,
        "affected": 1,
        "unaffected": None,
    }
    assert normalized["field_verdicts"]["unaffected"] == "ambiguous"


def test_verifier_cannot_support_zero_for_open_carrier_set():
    card = VariantClaimCard(
        gene="KCNH2",
        disease="Long QT syndrome",
        pmid="example",
        title="Familial long QT syndrome",
        variant="p.Gly572Ser",
        extracted={"total_carriers": None, "affected": None, "unaffected": 0},
        evidence="QT prolongation was observed in all mutation carriers.",
    )
    raw = {
        "verdict": "directly_supported",
        "field_verdicts": {
            "total_carriers": "source_missing",
            "affected": "source_missing",
            "unaffected": "directly_supported",
        },
        "corrected_values": {
            "total_carriers": None,
            "affected": None,
            "unaffected": 0,
        },
        "reason": "All carriers had QT prolongation.",
        "evidence_quote": "QT prolongation was observed in all mutation carriers.",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"]["unaffected"] is None
    assert normalized["field_verdicts"]["unaffected"] == "ambiguous"


def test_claim_verification_guard_clears_ambiguous_symptom_partition():
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

    assert normalized["corrected_values"]["total_carriers"] == 112
    assert normalized["corrected_values"]["affected"] is None
    assert normalized["corrected_values"]["unaffected"] is None
    assert normalized["field_verdicts"]["affected"] == "ambiguous"
    assert normalized["field_verdicts"]["unaffected"] == "ambiguous"


def test_claim_verification_does_not_promote_unsupported_partition_values():
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

    assert normalized["corrected_values"]["affected"] is None
    assert normalized["corrected_values"]["unaffected"] is None
    assert normalized["field_verdicts"]["affected"] == "unsupported"
    assert normalized["field_verdicts"]["unaffected"] == "unsupported"


def test_claim_verification_accepts_audited_patient_row_derivation():
    card = _derived_row_card()
    raw = {
        "verdict": "inferred_supported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "inferred_supported",
            "unaffected": "inferred_supported",
        },
        "corrected_values": {
            "total_carriers": 185,
            "affected": 97,
            "unaffected": 62,
        },
        "reason": "The complete row audit reconciles with 26 uncertain carriers.",
        "evidence_quote": "Supplementary Table 3",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"] == {
        "total_carriers": 185,
        "affected": 97,
        "unaffected": 62,
    }
    assert normalized["field_verdicts"]["affected"] == "inferred_supported"
    assert normalized["field_verdicts"]["unaffected"] == "inferred_supported"


def test_claim_verification_rejects_mismatched_patient_row_derivation():
    card = _derived_row_card(table_unaffected=63)
    raw = {
        "verdict": "inferred_supported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "inferred_supported",
            "unaffected": "inferred_supported",
        },
        "corrected_values": {
            "total_carriers": 185,
            "affected": 97,
            "unaffected": 62,
        },
        "reason": "Claimed row derivation",
        "evidence_quote": "Supplementary Table 3",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"]["affected"] is None
    assert normalized["corrected_values"]["unaffected"] is None
    assert normalized["field_verdicts"]["affected"] == "ambiguous"
    assert normalized["field_verdicts"]["unaffected"] == "ambiguous"


def test_claim_verification_rejects_model_declared_patient_row_type_without_stamps():
    card = _derived_row_card()
    provenance = card.derivation["count_provenance"]
    provenance.pop("affected_source")
    provenance.pop("unaffected_source")
    raw = {
        "verdict": "inferred_supported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "inferred_supported",
            "unaffected": "inferred_supported",
        },
        "corrected_values": {
            "total_carriers": 185,
            "affected": 97,
            "unaffected": 62,
        },
        "reason": "Model-declared row derivation",
        "evidence_quote": "Supplementary Table 3",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"]["affected"] is None
    assert normalized["corrected_values"]["unaffected"] is None
    assert normalized["field_verdicts"]["affected"] == "ambiguous"
    assert normalized["field_verdicts"]["unaffected"] == "ambiguous"


def test_claim_verification_rejects_correction_that_breaks_valid_row_audit():
    card = _derived_row_card()
    raw = {
        "verdict": "inferred_supported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "inferred_supported",
            "unaffected": "inferred_supported",
        },
        "corrected_values": {
            "total_carriers": 185,
            "affected": 98,
            "unaffected": 61,
        },
        "reason": "Changed without changing the row audit",
        "evidence_quote": "Supplementary Table 3",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"]["affected"] is None
    assert normalized["corrected_values"]["unaffected"] is None


def test_claim_verification_does_not_complete_partition_arithmetically():
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

    assert normalized["corrected_values"]["affected"] is None
    assert normalized["field_verdicts"]["affected"] == "unsupported"


def test_claim_verification_does_not_invent_partition_for_nested_cohort():
    card = VariantClaimCard(
        gene="KCNQ1",
        disease="Long QT syndrome type 1",
        pmid="33141630",
        title="Amish founder variant",
        variant="T224M",
        extracted={"total_carriers": 124, "affected": None, "unaffected": None},
        evidence=(
            "The variant was found in 124 carriers. Of those, 88 consented to "
            "phenotype follow-up and 34/88 had clinical evidence of LQTS."
        ),
    )
    raw = {
        "verdict": "unsupported",
        "field_verdicts": {
            "variant": "directly_supported",
            "total_carriers": "directly_supported",
            "affected": "directly_supported",
            "unaffected": "unsupported",
        },
        "corrected_values": {
            "total_carriers": 124,
            "affected": 34,
            "unaffected": None,
        },
        "reason": "Disease status is not established for every carrier.",
        "evidence_quote": "124 carriers; 34/88 had clinical evidence",
    }

    normalized = normalize_verification(raw, card=card)

    assert normalized["corrected_values"] == {
        "total_carriers": 124,
        "affected": 34,
        "unaffected": None,
    }
    assert normalized["field_verdicts"]["unaffected"] == "unsupported"


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


def test_claim_card_carries_structured_patient_row_derivation():
    card = build_claim_card(
        source_text="The p.G357S variant had one complete patient table.",
        gene="RYR2",
        disease="CPVT",
        pmid="25814417",
        title="Founder mutation",
        variant={
            "protein_notation": "p.G357S",
            "patients": {"count": 185},
            "penetrance_data": {
                "total_carriers_observed": 185,
                "affected_count": 97,
                "unaffected_count": 62,
                "uncertain_count": 26,
            },
            "count_provenance": {
                "affected_count_type": "derived_from_patient_rows",
                "unaffected_count_type": "derived_from_patient_rows",
            },
            "phenotype_derivation": {"method": "derived_from_patient_rows"},
        },
    )

    assert card is not None
    assert card.extracted["uncertain"] == 26
    assert card.derivation["count_provenance"]["affected_count_type"] == (
        "derived_from_patient_rows"
    )
    assert card.derivation["phenotype_derivation"]["method"] == (
        "derived_from_patient_rows"
    )


def test_claim_card_finds_one_letter_alias_for_three_letter_variant():
    text = "\n".join(
        [
            "A study-wide group included 63 participants.",
            ("Long converted paragraph filler. " * 40)
            + "One variant, SCN5A-R1193Q, was found in 19 participants; "
            "7 of 19 had arrhythmia phenotypes.",
        ]
    )

    for notation in ("p.Arg1193Gln", "Arg1193Gln", "ARG1193GLN"):
        card = build_claim_card(
            source_text=text,
            gene="SCN5A",
            disease="Long QT syndrome",
            pmid="26746457",
            title="Incidental variants",
            variant=_variant(notation, total=19, affected=7, unaffected=12),
        )

        assert card is not None
        assert "SCN5A-R1193Q" in card.evidence


def test_claim_card_bridges_stop_codon_star_and_x_aliases():
    card = build_claim_card(
        source_text="The SCN5A-R1193X variant was present in one participant.",
        gene="SCN5A",
        disease="Long QT syndrome",
        pmid="1",
        title="Stop alias",
        variant=_variant("p.Arg1193Ter", total=1, affected=1, unaffected=0),
    )

    assert card is not None
    assert "SCN5A-R1193X" in card.evidence


def test_source_verified_alias_partition_overrides_false_verifier_rejection():
    source = (
        "The index case carried p.R360_Q361dupQKQR.\n"
        "Overall, 24 out of 29 KCNQ1 dup12 mutation carriers were affected "
        "by an LQT syndrome, suggesting incomplete penetrance."
    )
    variant = {
        "gene_symbol": "KCNQ1",
        "protein_notation": "p.R360_Q361dupQKQR",
        "source_notation": "KCNQ1 dup12; p.R360_Q361dupQKQR",
        "patients": {"count": 29},
        "penetrance_data": {
            "total_carriers_observed": 29,
            "affected_count": 24,
            "unaffected_count": 5,
        },
        "fact_provenance": [
            {
                "fact_type": "affected_count",
                "fact_value": "24",
                "evidence_quote": (
                    "24 out of 29 KCNQ1 dup12 mutation carriers were affected "
                    "by an LQT syndrome"
                ),
            }
        ],
    }
    card = build_claim_card(
        source_text=source,
        gene="KCNQ1",
        disease="Long QT syndrome",
        pmid="1",
        title="Alias partition",
        variant=variant,
        max_evidence_chars=1000,
    )
    raw = {
        "verdict": "unsupported",
        "field_verdicts": {field: "unsupported" for field in FIELD_NAMES},
        "corrected_values": {
            "total_carriers": None,
            "affected": None,
            "unaffected": None,
        },
        "reason": "Canonical-token window missed the alias count line.",
        "evidence_quote": "",
    }

    assert card is not None
    assert "24 out of 29 KCNQ1 dup12" in card.evidence
    normalized = normalize_verification(raw, card=card)
    assert normalized["corrected_values"] == {
        "total_carriers": 29,
        "affected": 24,
        "unaffected": 5,
    }
    assert normalized["source_verified_overrides"] == [
        "total_carriers",
        "affected",
        "unaffected",
    ]


def test_source_verified_same_population_coreference_keeps_affected_count():
    source = (
        "The current study population consisted of 30 carriers of the D1790G "
        "SCN5A mutation.\n"
        "The study population consisted of 30 LQT3 patients."
    )
    variant = {
        "protein_notation": "p.Asp1790Gly",
        "source_notation": "D1790G",
        "patients": {"count": 30},
        "penetrance_data": {
            "total_carriers_observed": 30,
            "affected_count": 30,
            "unaffected_count": None,
        },
        "fact_provenance": [
            {
                "fact_type": "total_carriers_observed",
                "fact_value": "30",
                "evidence_quote": (
                    "The current study population consisted of 30 carriers of "
                    "the D1790G SCN5A mutation"
                ),
            },
            {
                "fact_type": "affected_count",
                "fact_value": "30",
                "evidence_quote": "The study population consisted of 30 LQT3 patients",
            },
        ],
    }

    claims = source_verified_claims(variant, source, "Long QT syndrome")

    assert claims["affected"]["value"] == 30
    assert claims["affected"]["method"] == (
        "same_population_identity_disease_coreference"
    )


def test_disease_cohort_without_equal_identity_population_is_not_verified():
    source = "We enrolled 30 LQT3 patients. The D1790G variant was identified in 4."
    variant = {
        "source_notation": "D1790G",
        "patients": {"count": 4},
        "penetrance_data": {
            "total_carriers_observed": 4,
            "affected_count": 30,
            "unaffected_count": None,
        },
        "fact_provenance": [
            {
                "fact_type": "affected_count",
                "fact_value": "30",
                "evidence_quote": "We enrolled 30 LQT3 patients",
            }
        ],
    }

    assert "affected" not in source_verified_claims(variant, source, "Long QT syndrome")


def test_symptom_subset_is_not_a_closed_disease_partition():
    source = "24 out of 29 KCNQ1 dup12 mutation carriers had syncope."
    variant = {
        "source_notation": "KCNQ1 dup12",
        "patients": {"count": 29},
        "penetrance_data": {
            "total_carriers_observed": 29,
            "affected_count": 24,
            "unaffected_count": 5,
        },
    }

    assert source_verified_claims(variant, source, "Long QT syndrome") == {}
