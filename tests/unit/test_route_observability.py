"""Route-observability regressions for the 2026-07-28 correction pass.

Each test in this file fails against the implementation it corrects; the
docstrings name the specific defect.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from utils.llm_utils import litellm_completion
from utils.llm_trace import (
    TRACE_MANIFEST_NAME,
    build_trace_manifest,
    configure_llm_tracing,
    reset_llm_tracing,
)


@pytest.fixture(autouse=True)
def _reset_trace_configuration():
    reset_llm_tracing()
    yield
    reset_llm_tracing()


def _records(trace_root: Path) -> list[dict]:
    return [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(trace_root.rglob("*.json"))
        if path.name != TRACE_MANIFEST_NAME
    ]


def _decision(trace_root: Path, event_type: str) -> dict:
    for record in _records(trace_root):
        if (record.get("event") or {}).get("type") == event_type:
            return record["event"]["data"]
    raise AssertionError(
        f"no {event_type} recorded; got "
        f"{[(r.get('event') or {}).get('type') or r.get('record_type') for r in _records(trace_root)]}"
    )


def _decisions(trace_root: Path, event_type: str) -> list[dict]:
    return [
        record["event"]["data"]
        for record in _records(trace_root)
        if (record.get("event") or {}).get("type") == event_type
    ]


def _chat_response(text: str):
    return SimpleNamespace(
        choices=[
            SimpleNamespace(message=SimpleNamespace(content=text), finish_reason="stop")
        ]
    )


def _fake_provider(monkeypatch, *texts_or_excs):
    """Patch litellm's `completion` — the layer BELOW `capture_llm_call`.

    Patching `litellm_completion` itself would replace the traced wrapper, so no
    `llm_call` record would be written and the test would prove nothing about
    linkage.
    """
    import utils.llm_utils as llm_utils

    served = iter(texts_or_excs)
    seen = []

    def completion(**kwargs):
        seen.append(kwargs)
        item = next(served)
        if isinstance(item, BaseException):
            raise item
        return _chat_response(item)

    monkeypatch.setattr(llm_utils, "completion", completion)
    return seen


# ---------------------------------------------------------------------------
# 1. Multi-model extraction: the retained result's call must be the accepted one.
# ---------------------------------------------------------------------------


class TestExtractionSelectionLinkage:
    def _extractor(self, monkeypatch, responses: list[str]):
        from pipeline.extraction import ExpertExtractor

        extractor = ExpertExtractor.__new__(ExpertExtractor)
        extractor.models = ["model-a", "model-b"]
        extractor.model = "model-a"
        extractor.temperature = 0.0
        extractor.max_tokens = 2048
        extractor.reasoning_effort = None
        extractor.tier_threshold = 100  # force trying the next model
        extractor._candidate_trace_ids = {}
        extractor.fulltext_dir = None

        served = iter(responses)

        def fake_attempt(paper, model, **kwargs):
            from utils.llm_trace import (
                capture_llm_call,
                note_llm_attempt,
                note_llm_outcome,
            )
            from utils.models import ExtractionResult

            payload = json.loads(next(served))
            _, summary = capture_llm_call(
                provider="fixture",
                requested_model=model,
                resolved_model=model,
                request={
                    "model": model,
                    "messages": [{"role": "user", "content": "x"}],
                },
                call=lambda: _chat_response(json.dumps(payload)),
            )
            trace_id = note_llm_attempt(summary, role="primary")
            note_llm_outcome(trace_id, "parsed")
            return ExtractionResult(
                pmid=paper.pmid,
                success=True,
                extracted_data=payload,
                model_used=model,
            )

        monkeypatch.setattr(extractor, "_attempt_extraction", fake_attempt)
        monkeypatch.setattr(
            extractor, "_prepare_full_text", lambda paper: paper.full_text or ""
        )
        monkeypatch.setattr(
            extractor, "_assess_input_quality", lambda text, gene: (True, "ok")
        )
        monkeypatch.setattr(extractor, "_estimate_table_rows", lambda text: 0)
        return extractor

    def test_first_model_wins_on_yield_and_is_the_accepted_link(
        self, tmp_path: Path, monkeypatch
    ):
        """Model A keeps the result; model B parsed LAST and used to claim accepted."""
        from utils.models import Paper

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="sel")
        # model-a: 3 variants, model-b: 1 variant. Threshold forces both to run.
        extractor = self._extractor(
            monkeypatch,
            [
                json.dumps(
                    {"variants": [{"protein_notation": f"p.A{i}B"} for i in range(3)]}
                ),
                json.dumps({"variants": [{"protein_notation": "p.C9D"}]}),
            ],
        )
        paper = Paper(pmid="12345678", title="T", abstract="A", full_text="body")
        paper.gene_symbol = "KCNH2"

        result = extractor.extract(paper)

        assert result.model_used == "model-a"
        assert len(result.extracted_data["variants"]) == 3
        calls = [r for r in _records(trace_root) if r["record_type"] == "llm_call"]
        assert len(calls) == 2
        by_model = {r["request"]["requested_model"]: r["trace_id"] for r in calls}
        data = _decision(trace_root, "paper_extraction_selection")
        assert data["selected_model"] == "model-a"
        assert data["accepted_response_trace_id"] == by_model["model-a"]
        assert data["accepted_response_trace_ids"] == [by_model["model-a"]]
        # The later successful call must be DISCARDED, not accepted.
        assert data["discarded_trace_ids"] == [by_model["model-b"]]
        outcomes = {
            link["trace_id"]: link["outcome"] for link in data["attempt_trace_links"]
        }
        assert outcomes[by_model["model-a"]] == "accepted"
        assert outcomes[by_model["model-b"]] == "discarded"
        # Attempt numbers are monotonic across the two candidates.
        assert [link["attempt"] for link in data["attempt_trace_links"]] == [1, 2]

    def test_second_model_wins_when_it_yields_more(self, tmp_path: Path, monkeypatch):
        from utils.models import Paper

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="sel2")
        extractor = self._extractor(
            monkeypatch,
            [
                json.dumps({"variants": [{"protein_notation": "p.A1B"}]}),
                json.dumps(
                    {"variants": [{"protein_notation": f"p.C{i}D"} for i in range(5)]}
                ),
            ],
        )
        paper = Paper(pmid="12345678", title="T", abstract="A", full_text="body")
        paper.gene_symbol = "KCNH2"

        result = extractor.extract(paper)

        assert result.model_used == "model-b"
        by_model = {
            r["request"]["requested_model"]: r["trace_id"]
            for r in _records(trace_root)
            if r["record_type"] == "llm_call"
        }
        data = _decision(trace_root, "paper_extraction_selection")
        assert data["accepted_response_trace_id"] == by_model["model-b"]
        assert data["discarded_trace_ids"] == [by_model["model-a"]]


# ---------------------------------------------------------------------------
# 2. Table routing.
# ---------------------------------------------------------------------------


ROUTER_TABLE = "| Variant | Something |\n| --- | --- |\n| p.Leu552Ser | 3 |\n"
#: A response the router will actually accept: `variant_tables`, with one
#: notation column and one count column (see parse_router_response).
_ROUTED_JSON = json.dumps(
    {
        "variant_tables": [
            {
                "table_id": "table_1",
                "column_mapping": {"protein": 0, "patient_count": 1},
                "confidence": 0.9,
                "notes": "carrier table",
            }
        ]
    }
)


class TestTableRoutingLinkage:
    def test_router_call_is_linked_as_accepted(self, tmp_path: Path, monkeypatch):
        """The ledger was opened but never fed: no accepted link existed at all."""
        from pipeline.table_router import route_tables

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="router")
        routed_json = _ROUTED_JSON
        _fake_provider(monkeypatch, routed_json)
        route_tables(
            ROUTER_TABLE,
            "KCNH2",
            model="fixture",
            llm_caller=litellm_completion,
            pmid="12345678",
        )

        calls = [r for r in _records(trace_root) if r["record_type"] == "llm_call"]
        assert len(calls) == 1
        data = _decision(trace_root, "table_routing_decision")
        assert data["decision_source"] == "model_routed"
        assert data["accepted_response_trace_id"] == calls[0]["trace_id"]
        assert data["provider_attempts"] == 1

    def test_empty_then_successful_retry_is_distinguishable(
        self, tmp_path: Path, monkeypatch
    ):
        from pipeline.table_router import route_tables

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="router-retry")
        _fake_provider(monkeypatch, "", _ROUTED_JSON)
        route_tables(
            ROUTER_TABLE,
            "KCNH2",
            model="fixture",
            llm_caller=litellm_completion,
            pmid="12345678",
        )

        data = _decision(trace_root, "table_routing_decision")
        assert data["provider_attempts"] == 2
        roles = [link["role"] for link in data["attempt_trace_links"]]
        assert roles == ["table_router", "empty_content_retry"]
        outcomes = [link["outcome"] for link in data["attempt_trace_links"]]
        assert outcomes == ["discarded", "accepted"]
        assert len(data["discarded_trace_ids"]) == 1

    def test_unusable_response_is_a_parse_failure_not_an_accepted_call(
        self, tmp_path: Path, monkeypatch
    ):
        from pipeline.table_router import route_tables

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="router-bad")
        _fake_provider(monkeypatch, "this is not json")
        route_tables(
            ROUTER_TABLE,
            "KCNH2",
            model="fixture",
            llm_caller=litellm_completion,
            pmid="12345678",
        )

        data = _decision(trace_root, "table_routing_decision")
        assert data["decision_source"] == "router_response_unusable"
        assert data["accepted_response_trace_id"] is None
        assert len(data["failed_trace_ids"]) == 1
        assert data["parse_failures"] == 1

    def test_provider_failure_is_recorded_as_a_failed_call(
        self, tmp_path: Path, monkeypatch
    ):
        from pipeline.table_router import route_tables

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="router-fail")
        _fake_provider(monkeypatch, ValueError("router down"))

        result = route_tables(
            ROUTER_TABLE,
            "KCNH2",
            model="fixture",
            llm_caller=litellm_completion,
            pmid="12345678",
        )

        assert result.error
        data = _decision(trace_root, "table_routing_decision")
        assert data["decision_source"] == "router_call_failed"
        assert data["accepted_response_trace_id"] is None
        assert data["failed_trace_ids"]

    def test_deterministic_route_claims_no_model_call(self, tmp_path: Path):
        """A fully deterministic route must not imply a model was consulted."""
        from pipeline.table_router import route_tables

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="router-det")
        deterministic_table = (
            "| Variant | Carriers |\n| --- | --- |\n| p.Leu552Ser | 3 |\n"
        )

        def must_not_be_called(**kwargs):
            raise AssertionError("router model called for a deterministic table")

        result = route_tables(
            deterministic_table,
            "KCNH2",
            model="fixture",
            llm_caller=must_not_be_called,
            pmid="12345678",
        )

        assert result.routed_tables
        data = _decision(trace_root, "table_routing_decision")
        assert data["decision_source"] == "deterministic"
        assert data["provider_attempts"] == 0
        assert data["accepted_response_trace_id"] is None
        assert data["model"] is None
        assert not [r for r in _records(trace_root) if r["record_type"] == "llm_call"]


# ---------------------------------------------------------------------------
# 3. Vision routes.
# ---------------------------------------------------------------------------


def _png(tmp_path: Path, name: str = "fig1.png") -> Path:
    path = tmp_path / name
    # Name-dependent bytes: the variant reader now skips byte-identical
    # duplicates per paper, and these tests are about trace linkage, not dedup.
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + name.encode() + b"0" * 32)
    return path


class TestFigureTextRoute:
    def test_successful_ocr_emits_a_linked_decision(self, tmp_path: Path, monkeypatch):
        """The route had scope but no decision event and no accepted link."""
        import harvesting.figure_text_extractor as mod

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="ocr")
        _fake_provider(monkeypatch, "## Table\n| a | b |")
        image = _png(tmp_path)

        out = mod.extract_images_to_markdown(
            [image], "anthropic/claude-fixture", gene="KCNH2", pmid="12345678"
        )

        assert "| a | b |" in out
        calls = [r for r in _records(trace_root) if r["record_type"] == "llm_call"]
        data = _decision(trace_root, "figure_text_extraction_decision")
        assert data["decision_source"] == "text_extracted"
        assert data["accepted_response_trace_id"] == calls[0]["trace_id"]
        assert data["image"] == "fig1.png"
        # Scope carries gene + PMID + image.
        decision_record = next(
            r
            for r in _records(trace_root)
            if (r.get("event") or {}).get("type") == "figure_text_extraction_decision"
        )
        assert decision_record["context"]["gene"] == "KCNH2"
        assert decision_record["context"]["pmid"] == "12345678"
        assert decision_record["context"]["stage"] == "figure_text_extraction"

    def test_empty_output_is_not_a_failure_and_not_accepted(
        self, tmp_path: Path, monkeypatch
    ):
        import harvesting.figure_text_extractor as mod

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="ocr-empty")
        _fake_provider(monkeypatch, "  ")

        out = mod.extract_images_to_markdown(
            [_png(tmp_path)], "anthropic/claude-fixture", gene="KCNH2", pmid="1"
        )

        assert out == ""  # biomedical outcome unchanged
        data = _decision(trace_root, "figure_text_extraction_decision")
        assert data["decision_source"] == "empty_output"
        assert data["text_extracted"] is False
        assert data["accepted_response_trace_id"] is None
        assert data["failed_trace_ids"] == []  # it did not FAIL, it found nothing

    def test_provider_failure_is_recorded_and_skipped(
        self, tmp_path: Path, monkeypatch
    ):
        import harvesting.figure_text_extractor as mod

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="ocr-fail")

        _fake_provider(monkeypatch, RuntimeError("vision down"))

        out = mod.extract_images_to_markdown(
            [_png(tmp_path)], "anthropic/claude-fixture", gene="KCNH2", pmid="1"
        )

        assert out == ""  # individual failures are skipped, as before
        data = _decision(trace_root, "figure_text_extraction_decision")
        assert data["decision_source"] == "call_failed"
        assert "vision down" in data["error"]
        assert data["failed_trace_ids"]


class TestFigureVariantRoute:
    def test_variants_read_emits_a_linked_decision(self, tmp_path: Path, monkeypatch):
        import harvesting.figure_variant_reader as mod

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="figvar")
        _fake_provider(
            monkeypatch, json.dumps({"variants": [{"variant": "p.Leu552Ser"}]})
        )

        report = mod.read_images(
            [_png(tmp_path)], "KCNH2", pmid="12345678", model="anthropic/x"
        )

        assert len(report.per_figure[0].variants) == 1
        calls = [r for r in _records(trace_root) if r["record_type"] == "llm_call"]
        data = _decision(trace_root, "figure_variant_read_decision")
        assert data["decision_source"] == "variants_read"
        assert data["variant_count"] == 1
        assert data["variants"] == ["p.Leu552Ser"]
        assert data["accepted_response_trace_id"] == calls[0]["trace_id"]

    def test_no_variants_and_unparseable_are_distinguishable(
        self, tmp_path: Path, monkeypatch
    ):
        """Both yield zero variants; a curator must be able to tell them apart."""
        import harvesting.figure_variant_reader as mod

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="figvar2")
        _fake_provider(
            monkeypatch,
            json.dumps({"variants": []}),
            "I could not read the figure.",
        )

        report = mod.read_images(
            [_png(tmp_path, "a.png"), _png(tmp_path, "b.png")],
            "KCNH2",
            pmid="12345678",
            model="anthropic/x",
        )

        assert [len(r.variants) for r in report.per_figure] == [0, 0]
        sources = [
            d["decision_source"]
            for d in _decisions(trace_root, "figure_variant_read_decision")
        ]
        assert sources == ["no_variants_reported", "unparseable_response"]
        clean, broken = _decisions(trace_root, "figure_variant_read_decision")
        assert clean["accepted_response_trace_id"] is not None
        assert broken["accepted_response_trace_id"] is None
        assert broken["parse_failures"] == 1

    def test_provider_failure_is_recorded(self, tmp_path: Path, monkeypatch):
        import harvesting.figure_variant_reader as mod

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="figvar3")

        _fake_provider(monkeypatch, RuntimeError("reader down"))

        report = mod.read_images(
            [_png(tmp_path)], "KCNH2", pmid="1", model="anthropic/x"
        )

        assert report.per_figure[0].error
        data = _decision(trace_root, "figure_variant_read_decision")
        assert data["decision_source"] == "call_failed"
        assert data["failed_trace_ids"]


class TestPedigreeRoute:
    def _extractor(self, monkeypatch, text_or_exc):
        from pipeline.pedigree_extractor import PedigreeExtractor

        extractor = PedigreeExtractor.__new__(PedigreeExtractor)
        extractor.model = "anthropic/claude-fixture"
        extractor.detection_threshold = 0.7
        extractor.gene = "KCNH2"
        extractor.pmid = "12345678"
        extractor._last_vision_status = None
        extractor._last_vision_error = None

        _fake_provider(monkeypatch, text_or_exc)
        monkeypatch.setattr(
            extractor, "_image_to_base64_url", lambda path: "data:image/png;base64,QQ=="
        )
        return extractor

    def test_detection_and_extraction_are_separate_linked_stages(
        self, tmp_path: Path, monkeypatch
    ):
        """One `pedigree_vision` stage could not tell detection from extraction."""
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="ped")
        image = _png(tmp_path)

        detector = self._extractor(
            monkeypatch,
            json.dumps({"is_pedigree": True, "confidence": 0.93, "reason": "symbols"}),
        )
        is_ped, confidence, _ = detector.is_pedigree(image)
        assert is_ped and confidence == 0.93

        extractor = self._extractor(
            monkeypatch,
            json.dumps(
                {
                    "individuals": [{"id": "II-1"}, {"id": "II-2"}],
                    "total_generations": 2,
                    "inheritance_pattern": "autosomal dominant",
                }
            ),
        )
        result = extractor.extract_pedigree(image)
        assert len(result["individuals"]) == 2

        stages = {
            r["context"]["stage"]
            for r in _records(trace_root)
            if r["record_type"] == "llm_call"
        }
        assert stages == {"pedigree_detection", "pedigree_extraction"}

        detect = _decision(trace_root, "pedigree_detection_decision")
        assert detect["decision_source"] == "model_detection"
        assert detect["is_pedigree"] is True
        assert detect["passes_threshold"] is True
        assert detect["accepted_response_trace_id"]

        extract = _decision(trace_root, "pedigree_extraction_decision")
        assert extract["decision_source"] == "model_extraction"
        assert extract["individual_count"] == 2
        assert extract["inheritance_pattern"] == "autosomal dominant"
        assert extract["accepted_response_trace_id"]
        assert (
            extract["accepted_response_trace_id"]
            != detect["accepted_response_trace_id"]
        )

    def test_unparseable_detection_is_not_an_accepted_response(
        self, tmp_path: Path, monkeypatch
    ):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="ped-bad")
        detector = self._extractor(monkeypatch, "not json at all")

        is_ped, confidence, reason = detector.is_pedigree(_png(tmp_path))

        assert (is_ped, confidence) == (False, 0.0)
        data = _decision(trace_root, "pedigree_detection_decision")
        assert data["decision_source"] == "unparseable_response"
        assert data["accepted_response_trace_id"] is None
        assert data["parse_failures"] == 1

    def test_failed_vision_call_is_recorded(self, tmp_path: Path, monkeypatch):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="ped-fail")
        detector = self._extractor(monkeypatch, RuntimeError("vision down"))

        assert detector.is_pedigree(_png(tmp_path))[0] is False

        data = _decision(trace_root, "pedigree_detection_decision")
        assert data["decision_source"] == "call_failed"
        assert "vision down" in data["error"]
        assert data["failed_trace_ids"]

    def test_extraction_missing_individuals_is_reported(
        self, tmp_path: Path, monkeypatch
    ):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="ped-noind")
        extractor = self._extractor(monkeypatch, json.dumps({"total_generations": 2}))

        assert extractor.extract_pedigree(_png(tmp_path)) is None

        data = _decision(trace_root, "pedigree_extraction_decision")
        assert data["decision_source"] == "missing_individuals_field"
        assert data["individual_count"] == 0


# ---------------------------------------------------------------------------
# 4. Claim debate pilot.
# ---------------------------------------------------------------------------


class TestClaimDebateRoute:
    def _card(self):
        from pipeline.claim_verifier import VariantClaimCard

        return VariantClaimCard(
            gene="KCNH2",
            disease="Long QT syndrome",
            pmid="12345678",
            title="A traced paper",
            variant="p.Leu552Ser",
            extracted={"total_carriers": 44, "affected": 12, "unaffected": 25},
            evidence="The p.Leu552Ser variant was identified in 44 carriers.",
        )

    def _verifier(self, monkeypatch, text_or_exc):
        from scripts.recall_audit.run_claim_debate_pilot import ClaimDebateVerifier

        verifier = ClaimDebateVerifier(model="fixture", max_tokens=512)
        # One failure per retry attempt: @llm_retry may re-enter the call.
        items = (
            [text_or_exc] * 12
            if isinstance(text_or_exc, BaseException)
            else [text_or_exc]
        )
        _fake_provider(monkeypatch, *items)
        return verifier

    def test_debate_is_scoped_and_emits_a_linked_decision(
        self, tmp_path: Path, monkeypatch
    ):
        """This route was fully unscoped: no gene, PMID, claim, stage or decision."""
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="debate")
        verifier = self._verifier(
            monkeypatch,
            json.dumps(
                {
                    "agreement": "disagree",
                    "field_agreement": {"total_carriers": "disagree"},
                    "reason": "The quote counts probands, not carriers.",
                }
            ),
        )

        out = verifier.debate(
            card=self._card(),
            baseline_model="baseline-model",
            baseline_verification={"verdict": "supported"},
        )

        assert out["agreement"] == "disagree"
        call = next(r for r in _records(trace_root) if r["record_type"] == "llm_call")
        assert call["context"]["gene"] == "KCNH2"
        assert call["context"]["pmid"] == "12345678"
        assert call["context"]["variant"] == "p.Leu552Ser"
        assert call["context"]["stage"] == "variant_claim_debate"
        data = _decision(trace_root, "claim_debate_decision")
        assert data["decision_source"] == "model_debate"
        assert data["agreement"] == "disagree"
        assert data["baseline_model"] == "baseline-model"
        assert data["accepted_response_trace_id"] == call["trace_id"]

    def test_debate_failure_is_linked_to_the_failed_call(
        self, tmp_path: Path, monkeypatch
    ):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="debate-fail")
        verifier = self._verifier(monkeypatch, RuntimeError("debater down"))

        with pytest.raises(Exception):
            verifier.debate(
                card=self._card(),
                baseline_model="baseline-model",
                baseline_verification={},
            )

        data = _decision(trace_root, "claim_debate_decision")
        assert data["decision_source"] == "call_failed"
        assert data["accepted_response_trace_id"] is None
        assert data["failed_trace_ids"]


# ---------------------------------------------------------------------------
# 5. Count recovery: provider success is not acceptance.
# ---------------------------------------------------------------------------


class TestCountRecoveryAcceptTiming:
    def _db(self, tmp_path: Path) -> Path:
        import sqlite3

        db = tmp_path / "T.db"
        con = sqlite3.connect(db)
        con.executescript(
            """
            CREATE TABLE variants(variant_id INTEGER PRIMARY KEY, gene_symbol TEXT,
                cdna_notation TEXT, protein_notation TEXT);
            CREATE TABLE variant_papers(variant_id INTEGER, pmid TEXT,
                source_location TEXT, key_quotes TEXT, source_layer TEXT);
            CREATE TABLE penetrance_data(penetrance_id INTEGER PRIMARY KEY,
                variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
                affected_count INTEGER, unaffected_count INTEGER);
            """
        )
        con.execute("INSERT INTO variants VALUES (1,'KCNH2',NULL,'p.Leu552Ser')")
        con.execute(
            "INSERT INTO variant_papers VALUES (1,'111','Table 2','','llm_table')"
        )
        con.execute("INSERT INTO penetrance_data (variant_id,pmid) VALUES (1,'111')")
        con.commit()
        con.close()
        return db

    SOURCE = "The p.Leu552Ser variant was identified in 44 carriers."

    def test_unparseable_response_cannot_claim_an_accepted_response(
        self, tmp_path: Path, monkeypatch
    ):
        """_call_text marked the call accepted BEFORE parse_response ran."""
        from pipeline.count_recovery import recover_counts

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="cr-bad")
        _fake_provider(monkeypatch, "the model rambled, no JSON here")
        stats = recover_counts(
            self._db(tmp_path),
            "KCNH2",
            source_for_pmid=lambda pmid: self.SOURCE,
            llm_caller=litellm_completion,
            model="fixture",
            dry_run=True,
        )

        assert stats["batch_failures"] == 1
        data = _decision(trace_root, "count_recovery_decision")
        assert data["decision_source"] == "batch_failed"
        assert data["accepted_response_trace_id"] is None, (
            "a failure decision must not claim an accepted response"
        )
        assert data["parse_failures"] == 1
        assert data["failed_trace_ids"]

    def test_parsed_response_with_zero_grounded_counts_is_still_accepted(
        self, tmp_path: Path, monkeypatch
    ):
        """Parsing succeeded; the decision reports zero grounded counts, not failure."""
        from pipeline.count_recovery import recover_counts

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="cr-zero")
        payload = json.dumps(
            {
                "variants": [
                    {
                        "variant": "p.Leu552Ser",
                        "carriers": 812,  # a denominator; the role gate refuses it
                        "quote": "We screened 812 probands for KCNH2 variants.",
                    }
                ]
            }
        )
        _fake_provider(monkeypatch, payload)
        stats = recover_counts(
            self._db(tmp_path),
            "KCNH2",
            source_for_pmid=lambda pmid: (
                "We screened 812 probands for KCNH2 variants. " + self.SOURCE
            ),
            llm_caller=litellm_completion,
            model="fixture",
            dry_run=True,
        )

        assert stats["batch_failures"] == 0
        assert stats["counts_accepted"] == 0
        assert stats["counts_rejected"] == 1
        data = _decision(trace_root, "count_recovery_decision")
        assert data["decision_source"] == "quote_grounded_validation"
        assert data["counts_accepted"] == 0
        # The provider response WAS usable input to the decision.
        assert data["accepted_response_trace_id"] is not None
        assert data["parse_failures"] == 0

    def test_provider_failure_has_no_accepted_response(
        self, tmp_path: Path, monkeypatch
    ):
        from pipeline.count_recovery import recover_counts

        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="cr-fail")
        _fake_provider(monkeypatch, RuntimeError("provider down"))

        recover_counts(
            self._db(tmp_path),
            "KCNH2",
            source_for_pmid=lambda pmid: self.SOURCE,
            llm_caller=litellm_completion,
            model="fixture",
            dry_run=True,
        )

        data = _decision(trace_root, "count_recovery_decision")
        assert data["accepted_response_trace_id"] is None
        assert data["parse_failures"] == 0
        assert data["failed_trace_ids"]


# ---------------------------------------------------------------------------
# 7. Route coverage must match within the stage.
# ---------------------------------------------------------------------------


class TestStageScopedCoverage:
    def test_an_event_from_another_stage_does_not_satisfy_a_stage(self, tmp_path: Path):
        """Coverage compared against a repo-wide event set, so gaps disappeared."""
        from utils.llm_trace import capture_llm_call, llm_trace_scope
        from utils.llm_trace_html import collect_trace_report_data

        run_dir = tmp_path / "run"
        trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="cov")

        # Paper A: count recovery WITH its decision.
        with llm_trace_scope(gene="KCNH2", pmid="111", stage="count_recovery"):
            capture_llm_call(
                provider="fixture",
                requested_model="m",
                resolved_model="m",
                request={"input": "a"},
                call=lambda: SimpleNamespace(output_text="{}"),
            )
        with llm_trace_scope(gene="KCNH2", pmid="111"):
            from utils.llm_trace import record_trace_event as rte

            rte("count_recovery_decision", {"accepted_response_trace_id": "x"})

        # Paper B: count recovery with NO decision. This is a real gap.
        with llm_trace_scope(gene="KCNH2", pmid="222", stage="count_recovery"):
            capture_llm_call(
                provider="fixture",
                requested_model="m",
                resolved_model="m",
                request={"input": "b"},
                call=lambda: SimpleNamespace(output_text="{}"),
            )

        build_trace_manifest(trace_root, run_id="cov")
        data = collect_trace_report_data(trace_root, run_dir=run_dir, run_id="cov")

        # The count_recovery stage has 2 calls and 1 decision, so it is satisfied
        # at run level; what must NOT happen is a stage with calls and zero
        # decisions being satisfied by another stage's event of the same name.
        stages = {s["stage"]: s for s in data["coverage"]["stages"]}
        assert stages["count_recovery"]["llm_calls"] == 2
        assert stages["count_recovery"]["decisions"] == 1
        assert "count_recovery_decision" in stages["count_recovery"]["event_types"]

    def test_stage_with_calls_and_no_event_of_its_own_is_a_gap(self, tmp_path: Path):
        """A different stage's decision must not satisfy this stage."""
        from utils.llm_trace import capture_llm_call, llm_trace_scope
        from utils.llm_trace import record_trace_event as rte
        from utils.llm_trace_html import collect_trace_report_data

        run_dir = tmp_path / "run"
        trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="cov2")

        with llm_trace_scope(gene="KCNH2", pmid="111", stage="count_recovery"):
            capture_llm_call(
                provider="fixture",
                requested_model="m",
                resolved_model="m",
                request={"input": "a"},
                call=lambda: SimpleNamespace(output_text="{}"),
            )
        with llm_trace_scope(gene="KCNH2", pmid="111"):
            rte("tier2_relevance_decision", {"accepted_response_trace_id": "x"})

        build_trace_manifest(trace_root, run_id="cov2")
        data = collect_trace_report_data(trace_root, run_dir=run_dir, run_id="cov2")

        gaps = {g["stage"]: g for g in data["coverage"]["missing_decision_links"]}
        assert "count_recovery" in gaps
        assert gaps["count_recovery"]["expected_event"] == "count_recovery_decision"

    def test_a_satisfied_stage_always_shows_its_own_decisions(self, tmp_path: Path):
        """The invariant global matching broke: 'satisfied' with `decisions: 0`.

        A decision emitted OUTSIDE any stage scope takes the event type as its
        stage. Matching a stage's expected event against a repo-wide set of every
        event type then reported the stage satisfied while its own row still
        showed zero decisions — a table that contradicted its own gap column.
        Decisions are now filed under the stage they DESCRIBE, so the row and the
        verdict agree.
        """
        from utils.llm_trace import capture_llm_call, llm_trace_scope
        from utils.llm_trace import record_trace_event as rte
        from utils.llm_trace_html import collect_trace_report_data

        run_dir = tmp_path / "run"
        trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="cov2b")

        with llm_trace_scope(gene="KCNH2", pmid="111", stage="count_recovery"):
            capture_llm_call(
                provider="fixture",
                requested_model="m",
                resolved_model="m",
                request={"input": "a"},
                call=lambda: SimpleNamespace(output_text="{}"),
            )
        # No stage in scope: record_trace_event defaults `stage` to the event type.
        with llm_trace_scope(gene="KCNH2", pmid="111"):
            rte("count_recovery_decision", {"accepted_response_trace_id": "x"})

        build_trace_manifest(trace_root, run_id="cov2b")
        data = collect_trace_report_data(trace_root, run_dir=run_dir, run_id="cov2b")

        stages = {s["stage"]: s for s in data["coverage"]["stages"]}
        gap_stages = {g["stage"] for g in data["coverage"]["missing_decision_links"]}
        # The decision is credited to the stage it describes...
        assert stages["count_recovery"]["decisions"] == 1
        assert "count_recovery_decision" in stages["count_recovery"]["event_types"]
        # ...and no stage is ever "satisfied" while showing zero decisions.
        for stage, entry in stages.items():
            if (
                entry["expected_event"]
                and entry["llm_calls"]
                and stage not in gap_stages
            ):
                assert entry["decisions"] > 0, (
                    f"{stage} reported satisfied with decisions=0"
                )
        # There is no phantom stage named after the event type.
        assert "count_recovery_decision" not in stages

    def test_continuation_and_adjudication_are_not_reported_as_gaps(
        self, tmp_path: Path
    ):
        """Their decision is folded into paper_extraction_selection by design."""
        from utils.llm_trace import capture_llm_call, llm_trace_scope
        from utils.llm_trace import record_trace_event as rte
        from utils.llm_trace_html import collect_trace_report_data

        run_dir = tmp_path / "run"
        trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="cov3")
        for stage in (
            "paper_variant_extraction",
            "paper_extraction_continuation",
            "paper_extraction_adjudication",
        ):
            with llm_trace_scope(gene="KCNH2", pmid="111", stage=stage):
                capture_llm_call(
                    provider="fixture",
                    requested_model="m",
                    resolved_model="m",
                    request={"input": stage},
                    call=lambda: SimpleNamespace(output_text="{}"),
                )
        with llm_trace_scope(gene="KCNH2", pmid="111"):
            rte("paper_extraction_selection", {"accepted_response_trace_id": "x"})

        build_trace_manifest(trace_root, run_id="cov3")
        data = collect_trace_report_data(trace_root, run_dir=run_dir, run_id="cov3")

        gap_stages = {g["stage"] for g in data["coverage"]["missing_decision_links"]}
        assert "paper_extraction_continuation" not in gap_stages
        assert "paper_extraction_adjudication" not in gap_stages
        assert "paper_variant_extraction" not in gap_stages

    def test_every_registered_vision_and_debate_stage_is_expected(self):
        from utils.llm_trace_html import EXPECTED_DECISION_EVENTS

        for stage in (
            "figure_text_extraction",
            "figure_variant_read",
            "pedigree_detection",
            "pedigree_extraction",
            "variant_claim_debate",
        ):
            assert stage in EXPECTED_DECISION_EVENTS, stage


# ---------------------------------------------------------------------------
# 8. HTML linkage safety.
# ---------------------------------------------------------------------------


class TestReportLinkSafety:
    def test_no_jump_button_for_a_trace_id_not_on_this_page(self, tmp_path: Path):
        """A dead button told the curator a record was reachable when it was not."""
        from utils.llm_trace import llm_trace_scope
        from utils.llm_trace import record_trace_event as rte
        from utils.llm_trace_html import collect_trace_report_data, render_trace_report

        run_dir = tmp_path / "run"
        trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="links")
        with llm_trace_scope(gene="KCNH2", pmid="111", stage="count_recovery"):
            rte(
                "count_recovery_decision",
                {
                    "accepted_response_trace_id": "not-embedded-anywhere",
                    "discarded_trace_ids": ["also-missing"],
                    "attempt_trace_links": [
                        {
                            "attempt": 1,
                            "role": "primary",
                            "trace_id": None,
                            "outcome": "error",
                        },
                    ],
                },
            )
        build_trace_manifest(trace_root, run_id="links")
        data = collect_trace_report_data(trace_root, run_dir=run_dir, run_id="links")

        rendered = render_trace_report(data)

        embedded = json.loads(
            rendered.split("const DATA = ", 1)[1].split(";\n    const state", 1)[0]
        )["embedded_trace_ids"]
        # Nothing referenced by the decision is embedded, so no jump target is
        # offered for any of them.
        assert "not-embedded-anywhere" not in embedded
        assert "also-missing" not in embedded
        # The viewer degrades to labelled text instead of a dead button.
        assert "not on this page" in rendered
        assert "(untraced)" in rendered

    def test_accepted_trace_ids_array_is_rendered_when_several_contribute(
        self, tmp_path: Path
    ):
        from utils.llm_trace import capture_llm_call, llm_trace_scope
        from utils.llm_trace import record_trace_event as rte
        from utils.llm_trace_html import collect_trace_report_data, render_trace_report

        run_dir = tmp_path / "run"
        trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="multi")
        ids = []
        with llm_trace_scope(
            gene="KCNH2", pmid="111", stage="paper_variant_extraction"
        ):
            for label in ("primary", "continuation"):
                _, summary = capture_llm_call(
                    provider="fixture",
                    requested_model="m",
                    resolved_model="m",
                    request={"input": label},
                    call=lambda: SimpleNamespace(output_text="{}"),
                )
                ids.append(summary["trace_id"])
            rte(
                "paper_extraction_selection",
                {
                    "accepted_response_trace_id": ids[1],
                    "accepted_response_trace_ids": ids,
                },
            )
        build_trace_manifest(trace_root, run_id="multi")
        data = collect_trace_report_data(trace_root, run_dir=run_dir, run_id="multi")

        rendered = render_trace_report(data)

        assert "all accepted trace ids" in rendered
        # The builder tells the viewer which targets exist on this page, so the
        # jump buttons it renders at runtime cannot be dead.
        embedded = json.loads(
            rendered.split("const DATA = ", 1)[1].split(";\n    const state", 1)[0]
        )["embedded_trace_ids"]
        assert set(ids).issubset(set(embedded))

    def test_same_tree_layout_gets_a_relative_record_link(self, tmp_path: Path):
        from utils.llm_trace import capture_llm_call, llm_trace_scope
        from utils.llm_trace_html import build_trace_html_report

        run_dir = tmp_path / "run"
        trace_root = configure_llm_tracing(run_dir / "llm_traces", run_id="href")
        with llm_trace_scope(gene="KCNH2", pmid="111", stage="count_recovery"):
            capture_llm_call(
                provider="fixture",
                requested_model="m",
                resolved_model="m",
                request={"input": "x" * 5000},
                call=lambda: SimpleNamespace(output_text="ok"),
            )
        build_trace_manifest(trace_root, run_id="href")

        data = build_trace_html_report(
            trace_root,
            output_path=run_dir / "llm_trace_report.html",
            run_dir=run_dir,
            run_id="href",
            max_field_chars=100,
        )

        assert data["trace_href"] == "llm_traces"
        html = (run_dir / "llm_trace_report.html").read_text(encoding="utf-8")
        assert '"trace_href":"llm_traces"' in html
        assert "This card shows a bounded copy" in html

    def test_cross_volume_layout_shows_the_path_as_text(self, tmp_path: Path):
        from utils.llm_trace import capture_llm_call, llm_trace_scope
        from utils.llm_trace_html import build_trace_html_report

        # Trace root far from the report: no sane relative path.
        trace_root = configure_llm_tracing(
            tmp_path / "vol" / "a" / "b" / "c" / "traces", run_id="cross"
        )
        report_dir = tmp_path / "reports" / "deep"
        report_dir.mkdir(parents=True)
        with llm_trace_scope(gene="KCNH2", pmid="111", stage="count_recovery"):
            capture_llm_call(
                provider="fixture",
                requested_model="m",
                resolved_model="m",
                request={"input": "y"},
                call=lambda: SimpleNamespace(output_text="ok"),
            )
        build_trace_manifest(trace_root, run_id="cross")

        data = build_trace_html_report(
            trace_root,
            output_path=report_dir / "report.html",
            run_id="cross",
        )

        assert data["trace_href"] is None
        html = (report_dir / "report.html").read_text(encoding="utf-8")
        assert '"trace_href":null' in html
