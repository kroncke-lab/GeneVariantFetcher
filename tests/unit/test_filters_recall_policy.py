from pipeline.filters import ClinicalDataTriageFilter, InternFilter
from utils.models import FilterDecision, Paper


def test_intern_filter_missing_metadata_fails_open():
    result = InternFilter(model="stub").filter(Paper(pmid="1"))

    assert result.decision == FilterDecision.PASS
    assert result.metadata["fail_open"] is True


def test_clinical_triage_missing_metadata_fails_open():
    result = ClinicalDataTriageFilter(model="stub").triage(
        "", "", gene="BRCA2", pmid="1"
    )

    assert result["decision"] == "KEEP"
    assert result["fail_open"] is True


def test_intern_filter_prompt_uses_target_gene_without_kcnh2_examples(monkeypatch):
    captured = {}

    def fake_call_llm_json(self, prompt):
        captured["prompt"] = prompt
        return {"decision": "PASS", "reason": "cohort", "confidence": 0.9}

    monkeypatch.setattr(InternFilter, "call_llm_json", fake_call_llm_json)

    result = InternFilter(model="stub").filter(
        Paper(
            pmid="2",
            title="Genotyped hereditary cancer cohort",
            abstract="Patients underwent panel sequencing and variant analysis.",
            gene_symbol="BRCA2",
        )
    )

    assert result.decision == FilterDecision.PASS
    assert "BRCA2" in captured["prompt"]
    assert "KCNH2" not in captured["prompt"]
    assert "HERG" not in captured["prompt"]
    assert "LQTS" not in captured["prompt"]


def test_intern_filter_fails_open_when_llm_rejects_target_gene_screening(monkeypatch):
    def fake_call_llm_json(self, prompt):
        return {
            "decision": "FAIL",
            "reason": "KCNQ1-focused cohort",
            "confidence": 0.82,
        }

    monkeypatch.setattr(InternFilter, "call_llm_json", fake_call_llm_json)

    result = InternFilter(model="stub").filter(
        Paper(
            pmid="26669661",
            title="Genotype-phenotype relationship in LQT1, LQT2 and LQT3 patients",
            abstract=(
                "We screened patients with long QT syndrome for mutations in "
                "KCNQ1, KCNH2 and SCN5A and studied carriers in LQT1, LQT2 and LQT3."
            ),
            gene_symbol="SCN5A",
        )
    )

    assert result.decision == FilterDecision.PASS
    assert result.metadata["fail_open_target_gene_signal"] is True
    assert "Original classifier reason" in result.reason


def test_intern_filter_does_not_fail_open_unrelated_reviews(monkeypatch):
    def fake_call_llm_json(self, prompt):
        return {"decision": "FAIL", "reason": "review only", "confidence": 0.95}

    monkeypatch.setattr(InternFilter, "call_llm_json", fake_call_llm_json)

    result = InternFilter(model="stub").filter(
        Paper(
            pmid="3",
            title="Systematic review of KCNH2 channel structure",
            abstract="This review summarizes KCNH2 electrophysiology guidelines.",
            gene_symbol="KCNH2",
        )
    )

    assert result.decision == FilterDecision.FAIL
    assert result.metadata["fail_open_target_gene_signal"] is False
