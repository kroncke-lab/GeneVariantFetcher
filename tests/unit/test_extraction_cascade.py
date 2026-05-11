from pipeline.extraction import ExpertExtractor
from utils.models import ExtractionResult, Paper


def _result(pmid: str, model: str, variant_count: int) -> ExtractionResult:
    return ExtractionResult(
        pmid=pmid,
        success=True,
        model_used=model,
        extracted_data={
            "paper_metadata": {"pmid": pmid},
            "variants": [
                {
                    "gene_symbol": "KCNH2",
                    "protein_notation": f"Ala{idx + 1}Val",
                    "penetrance_data": {"total_carriers_observed": 1},
                }
                for idx in range(variant_count)
            ],
            "extraction_metadata": {"total_variants_found": variant_count},
        },
    )


def _patch_extractor(monkeypatch, results):
    monkeypatch.setattr(
        ExpertExtractor, "_prepare_full_text", lambda self, paper: "usable full text"
    )
    monkeypatch.setattr(
        ExpertExtractor, "_assess_input_quality", lambda self, text, gene: (True, "ok")
    )
    monkeypatch.setattr(ExpertExtractor, "_estimate_table_rows", lambda self, text: 0)

    calls = []

    def fake_attempt(self, paper, model, prepared_full_text, estimated_variants):
        calls.append(model)
        result = results[model]
        return result(paper.pmid) if callable(result) else result

    monkeypatch.setattr(ExpertExtractor, "_attempt_extraction", fake_attempt)
    return calls


def test_cascade_keeps_primary_when_fallback_finds_fewer_variants(monkeypatch):
    calls = _patch_extractor(
        monkeypatch,
        {
            "primary": lambda pmid: _result(pmid, "primary", 1),
            "fallback": lambda pmid: _result(pmid, "fallback", 0),
        },
    )
    extractor = ExpertExtractor(models=["primary", "fallback"], tier_threshold=2)

    result = extractor.extract(Paper(pmid="1", full_text="text", gene_symbol="KCNH2"))

    assert calls == ["primary", "fallback"]
    assert result.model_used == "primary"
    assert len(result.extracted_data["variants"]) == 1


def test_cascade_uses_fallback_when_it_finds_more_variants(monkeypatch):
    calls = _patch_extractor(
        monkeypatch,
        {
            "primary": lambda pmid: _result(pmid, "primary", 0),
            "fallback": lambda pmid: _result(pmid, "fallback", 1),
        },
    )
    extractor = ExpertExtractor(models=["primary", "fallback"], tier_threshold=2)

    result = extractor.extract(Paper(pmid="1", full_text="text", gene_symbol="KCNH2"))

    assert calls == ["primary", "fallback"]
    assert result.model_used == "fallback"
    assert len(result.extracted_data["variants"]) == 1


def test_cascade_uses_fallback_after_primary_failure(monkeypatch):
    failure = ExtractionResult(
        pmid="1",
        success=False,
        model_used="primary",
        error="rate limit",
    )
    calls = _patch_extractor(
        monkeypatch,
        {
            "primary": failure,
            "fallback": lambda pmid: _result(pmid, "fallback", 1),
        },
    )
    extractor = ExpertExtractor(models=["primary", "fallback"], tier_threshold=2)

    result = extractor.extract(Paper(pmid="1", full_text="text", gene_symbol="KCNH2"))

    assert calls == ["primary", "fallback"]
    assert result.model_used == "fallback"
    assert len(result.extracted_data["variants"]) == 1
