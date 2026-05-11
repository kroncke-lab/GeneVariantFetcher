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
