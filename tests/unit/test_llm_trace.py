from __future__ import annotations

import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import pytest

from gene_literature.llm_relevance import RelevanceChecker
from utils.llm_trace import (
    TRACE_INDEX_NAME,
    TRACE_MANIFEST_NAME,
    TraceRunMismatchError,
    attempt_link_summary,
    build_trace_manifest,
    capture_llm_call,
    configure_llm_tracing,
    exported_trace_identity,
    finalize_attempt_selection,
    integrity_level,
    json_safe,
    last_llm_trace,
    ledger_mark,
    ledger_trace_ids_since,
    llm_attempt_ledger,
    llm_trace_scope,
    note_llm_attempt,
    note_llm_outcome,
    reasoning_capture,
    record_trace_event,
    reset_llm_tracing,
    resolve_trace_location,
    trace_lock_targets,
    validate_trace_manifest,
)


@pytest.fixture(autouse=True)
def _reset_trace_configuration():
    reset_llm_tracing()
    yield
    reset_llm_tracing()


def test_call_trace_keeps_prompt_response_and_redacts_secrets(tmp_path: Path):
    trace_root = configure_llm_tracing(tmp_path / "traces", run_id="run-1")
    response = SimpleNamespace(
        id="resp-1",
        output_text='{"decision":"KEEP","reason":"patient table"}',
        usage=SimpleNamespace(input_tokens=12, output_tokens=4, total_tokens=16),
        reasoning=SimpleNamespace(summary="The table is original patient data."),
    )
    request = {
        "model": "model-a",
        "messages": [
            {
                "role": "user",
                "content": "Target gene: KCNH2\nPMID: 12345678\nExact prompt",
            }
        ],
        "api_key": "must-not-leak",
        "max_tokens": 200,
        "image_url": "data:image/png;base64,QUJD",
    }

    with llm_trace_scope(
        gene="KCNH2",
        pmid="12345678",
        stage="paper_curation",
        component="fixture",
    ):
        returned, summary = capture_llm_call(
            provider="fixture_provider",
            requested_model="model-a",
            resolved_model="provider/model-a",
            request=request,
            call=lambda: response,
        )

    assert returned is response
    assert summary is not None
    trace = json.loads((trace_root / summary["path"]).read_text())
    payload = trace["request"]["payload"]
    assert payload["messages"][0]["content"].endswith("Exact prompt")
    assert payload["api_key"] == "<redacted>"
    assert payload["image_url"]["inline_data_omitted"] is True
    assert trace["response"]["output_text"].startswith('{"decision"')
    assert trace["response"]["usage"]["total_tokens"] == 16
    assert trace["reasoning_capture"]["provider_exposed_reasoning_available"] is True
    assert "must-not-leak" not in json.dumps(trace)


def test_failed_call_is_traced_before_error_is_raised(tmp_path: Path):
    trace_root = configure_llm_tracing(tmp_path / "traces")

    with pytest.raises(RuntimeError, match="provider unavailable"):
        capture_llm_call(
            provider="fixture",
            requested_model="model-a",
            resolved_model=None,
            request={"input": "PMID: 12345678"},
            call=lambda: (_ for _ in ()).throw(RuntimeError("provider unavailable")),
        )

    traces = list(trace_root.rglob("*.json"))
    assert len(traces) == 1
    payload = json.loads(traces[0].read_text())
    assert payload["response"]["success"] is False
    assert payload["response"]["error"]["type"] == "RuntimeError"


def test_manifest_detects_trace_tampering_and_extra_files(tmp_path: Path):
    trace_root = configure_llm_tracing(tmp_path / "traces", run_id="run-2")
    with llm_trace_scope(gene="SCN5A", pmid="27566755", stage="route"):
        record_trace_event(
            "representation_route_decision",
            {"selected_tool": "text", "rationale": "only complete representation"},
        )
    manifest_path = trace_root / TRACE_MANIFEST_NAME
    manifest = build_trace_manifest(
        trace_root, output_path=manifest_path, run_id="run-2"
    )

    assert validate_trace_manifest(trace_root, manifest) == []
    record_path = trace_root / manifest["records"][0]["path"]
    record_path.write_text(record_path.read_text() + " ")
    assert any(
        "changed after manifest" in error
        for error in validate_trace_manifest(trace_root, manifest)
    )


def test_direct_anthropic_relevance_call_and_decision_are_traced(tmp_path: Path):
    trace_root = configure_llm_tracing(tmp_path / "traces", run_id="run-3")
    seen_request = {}
    response_text = (
        '{"is_relevant": false, "confidence": 0.93, '
        '"reasoning": "No variant-level evidence."}'
    )
    message = SimpleNamespace(
        content=[SimpleNamespace(text=response_text)],
        model_dump=lambda mode=None: {
            "id": "msg-1",
            "content": [{"type": "text", "text": response_text}],
            "usage": {"input_tokens": 41, "output_tokens": 18},
        },
    )

    def create_message(**request):
        seen_request.update(request)
        return message

    checker = RelevanceChecker(api_key="fixture-key", model="claude-fixture")
    checker._client = SimpleNamespace(messages=SimpleNamespace(create=create_message))

    decision = checker.check_relevance(
        gene_name="SCN5A",
        title="Expression-only study",
        abstract="No patient variants were studied.",
        pmid="12345678",
    )

    assert decision.is_relevant is False
    assert seen_request["model"] == "claude-fixture"
    manifest = build_trace_manifest(trace_root, run_id="run-3")
    assert manifest["llm_call_count"] == 1
    assert manifest["decision_event_count"] == 1
    records = [
        json.loads((trace_root / record["path"]).read_text())
        for record in manifest["records"]
    ]
    call = next(record for record in records if record["record_type"] == "llm_call")
    final = next(
        record for record in records if record["record_type"] == "decision_event"
    )
    assert (
        "Expression-only study" in call["request"]["payload"]["messages"][0]["content"]
    )
    assert call["response"]["envelope"]["content"][0]["text"] == response_text
    assert final["event"]["data"]["is_relevant"] is False
    assert final["event"]["data"]["accepted_response_trace_id"] == call["trace_id"]


# ---------------------------------------------------------------------------
# Concurrency: one shared caller must never cross-attribute another paper.
# The prior implementation kept {gene, pmid} on the BaseLLMCaller instance; an
# 8-thread probe bound 7 of 8 calls to another thread's PMID, so a curator
# asking "why was PMID X dropped" read another paper's evidence.
# ---------------------------------------------------------------------------


def test_shared_caller_under_concurrency_keeps_per_thread_attribution(tmp_path: Path):
    from pipeline.filters import InternFilter
    from utils.models import Paper

    trace_root = configure_llm_tracing(tmp_path / "traces", run_id="race")
    pmids = [f"3000000{index}" for index in range(8)]

    # One shared filter instance, exactly as pipeline/steps.py builds it.
    shared = InternFilter.__new__(InternFilter)
    shared.model = "fixture-model"
    shared.temperature = 0.0
    shared.max_tokens = 256
    shared.reasoning_effort = None
    shared.confidence_threshold = 0.5
    shared.disease = None
    shared.gene_aliases = ()

    def classify(pmid: str) -> None:
        paper = Paper(pmid=pmid, title=f"T{pmid}", abstract=f"A{pmid} variant carrier")
        paper.gene_symbol = "KCNH2"
        with (
            llm_trace_scope(
                gene=paper.gene_symbol,
                pmid=paper.pmid,
                stage="tier2_relevance_filter",
                component=type(shared).__name__,
            ),
            llm_attempt_ledger(),
        ):
            # Deterministic stand-in for the provider call: capture_llm_call is
            # the thing under test, not litellm.
            with shared._trace_scope("json", attempt=1):
                capture_llm_call(
                    provider="fixture",
                    requested_model=shared.model,
                    resolved_model=shared.model,
                    request={"model": shared.model, "input": f"classify {pmid}"},
                    call=lambda: SimpleNamespace(output_text='{"decision":"PASS"}'),
                )
            shared.record_llm_decision(
                "tier2_relevance_decision", {"final_decision": "PASS", "pmid": pmid}
            )

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(classify, pmids))

    paths = [
        path for path in trace_root.rglob("*.json") if path.name != TRACE_MANIFEST_NAME
    ]
    assert len(paths) == 16  # 8 calls + 8 decisions
    seen_pmids: set[str] = set()
    for path in paths:
        record = json.loads(path.read_text(encoding="utf-8"))
        pmid = record["context"]["pmid"]
        assert pmid in pmids
        seen_pmids.add(pmid)
        # Every record lands under its OWN gene/PMID directory.
        assert path.parent == trace_root / "KCNH2" / pmid
        blob = json.dumps(record)
        assert pmid in blob
        for other in pmids:
            if other != pmid:
                assert other not in blob, (
                    f"record for {pmid} leaked {other}: cross-thread attribution"
                )
    assert seen_pmids == set(pmids)


def test_trace_scope_accepts_caller_component_without_crashing():
    """`_trace_scope` used to raise TypeError on a duplicate `component` kwarg."""
    from utils.llm_utils import BaseLLMCaller

    caller = BaseLLMCaller(model="fixture", max_tokens=16)
    with llm_trace_scope(component="Outer", operation="outer"):
        with caller._trace_scope("json", attempt=2):
            from utils.llm_trace import current_trace_context

            context = current_trace_context()
    assert context["component"] == "BaseLLMCaller"
    assert context["operation"] == "json"
    assert context["attempt"] == 2


# ---------------------------------------------------------------------------
# Redaction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "key",
    [
        "api_key",
        "api-key",
        "API-Key",
        "x-api-key",
        "azure-ai-api-key",
        "openai_api_key",
        "authorization",
        "Authorization",
        "proxy-authorization",
        "bearer",
        "cookie",
        "set-cookie",
        "token",
        "refresh_token",
        "insttoken",
        "X-ELS-Insttoken",
        "client_secret",
        "password",
        "subscription-key",
    ],
)
def test_secret_key_spellings_are_redacted(key: str):
    assert json_safe({key: "sk-SECRET"})[key] == "<redacted>"


def test_header_maps_are_redacted_by_default_with_a_safe_allow_list():
    safe = json_safe(
        {
            "headers": {
                "api-key": "sk-SECRET",
                "X-ELS-Insttoken": "INST-SECRET",
                "x-ms-client-request-id": "OPAQUE-BUT-UNKNOWN",
                "Content-Type": "application/json",
                "api-version": "v1",
            },
            "extra_headers": {"Authorization": "Bearer SECRET"},
        }
    )
    headers = safe["headers"]
    assert headers["api-key"] == "<redacted>"
    assert headers["X-ELS-Insttoken"] == "<redacted>"
    # Unknown header names are redacted by policy, not retained by default.
    assert headers["x-ms-client-request-id"] == "<redacted>"
    assert headers["Content-Type"] == "application/json"
    assert headers["api-version"] == "v1"
    assert safe["extra_headers"]["Authorization"] == "<redacted>"
    assert "SECRET" not in json.dumps(safe)


def test_nested_request_bodies_are_redacted_at_depth():
    safe = json_safe(
        {
            "request": {
                "endpoint": "https://example/openai/v1/responses",
                "session": {"auth": {"access_token": "SECRET-TOKEN"}},
                "body": {"input": [{"type": "input_text", "text": "keep me"}]},
                "clients": [{"default_headers": {"api-key": "SECRET-2"}}],
            }
        }
    )
    blob = json.dumps(safe)
    assert "SECRET-TOKEN" not in blob
    assert "SECRET-2" not in blob
    assert "keep me" in blob


def test_unknown_objects_never_fall_back_to_repr():
    class Client:
        def __init__(self):
            self.api_key = "sk-SECRET"

        def __repr__(self):  # the exact leak the repr() fallback allowed
            return f"<Client key={self.api_key}>"

    class Opaque:
        __slots__ = ()

        def __repr__(self):
            return "<Opaque secret=sk-ALSO-SECRET>"

    safe = json_safe({"client": Client(), "opaque": Opaque()})
    assert safe["client"]["api_key"] == "<redacted>"
    assert safe["opaque"] == {
        "unserializable_object": True,
        "type": "Opaque",
        "module": safe["opaque"]["module"],
    }
    assert "SECRET" not in json.dumps(safe)


# ---------------------------------------------------------------------------
# Reasoning capture: a token count is never exposed reasoning.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("tokens", [0, 1024])
def test_reasoning_token_counts_never_claim_exposed_reasoning(tokens: int):
    capture = reasoning_capture(
        {"usage": {"completion_tokens_details": {"reasoning_tokens": tokens}}}
    )
    assert capture["provider_exposed_reasoning_available"] is False
    assert capture["response_paths"] == []
    assert capture["reasoning_token_usage"] == {
        "$.usage.completion_tokens_details.reasoning_tokens": tokens
    }


def test_reasoning_content_is_reported_as_exposed():
    capture = reasoning_capture(
        {
            "reasoning": {"effort": "high", "summary": "Considered the carrier table."},
            "usage": {"reasoning_tokens": 900},
        }
    )
    assert capture["provider_exposed_reasoning_available"] is True
    assert capture["response_paths"] == ["$.reasoning.summary"]
    assert capture["reasoning_token_usage"] == {"$.usage.reasoning_tokens": 900}


def test_effort_only_reasoning_block_is_not_exposed_reasoning():
    capture = reasoning_capture({"reasoning": {"effort": "xhigh"}})
    assert capture["provider_exposed_reasoning_available"] is False


def test_thinking_blocks_are_reported_as_exposed():
    capture = reasoning_capture(
        {"output": [{"type": "reasoning", "summary": [{"text": "step"}]}]}
    )
    assert capture["provider_exposed_reasoning_available"] is True


# ---------------------------------------------------------------------------
# Integrity: pre-manifest tampering, index digests, honest levels.
# ---------------------------------------------------------------------------


def _one_record(trace_root: Path) -> Path:
    with llm_trace_scope(gene="KCNH2", pmid="12345678", stage="paper_curation"):
        capture_llm_call(
            provider="fixture",
            requested_model="model-a",
            resolved_model="model-a",
            request={"input": "exact prompt"},
            call=lambda: SimpleNamespace(output_text="ok"),
        )
    return next(
        path for path in trace_root.rglob("*.json") if path.name != TRACE_MANIFEST_NAME
    )


def test_manifest_detects_tampering_that_happened_before_it_existed(tmp_path: Path):
    """The defect: forging a record pre-manifest was silently blessed by re-hash."""
    trace_root = configure_llm_tracing(tmp_path / "traces", run_id="pre")
    record_path = _one_record(trace_root)

    payload = json.loads(record_path.read_text(encoding="utf-8"))
    payload["request"]["payload"]["input"] = "FORGED prompt"
    record_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    manifest = build_trace_manifest(trace_root, run_id="pre")

    assert manifest["verification"]["level"] == "generated_now"
    assert any(
        error.startswith("write_time_digest_mismatch")
        for error in manifest["verification"]["errors"]
    )
    assert any(
        error.startswith("write_time_digest_mismatch")
        for error in validate_trace_manifest(trace_root, manifest)
    )
    level, _ = integrity_level(manifest, validate_trace_manifest(trace_root, manifest))
    assert level == "generated_now"


def test_untampered_manifest_reaches_write_time_verified(tmp_path: Path):
    trace_root = configure_llm_tracing(tmp_path / "traces", run_id="clean")
    _one_record(trace_root)

    manifest = build_trace_manifest(trace_root, run_id="clean")

    assert manifest["verification"]["level"] == "write_time_verified"
    assert manifest["verification"]["errors"] == []
    assert manifest["index"]["present"] is True
    assert manifest["index"]["sha256"]
    assert validate_trace_manifest(trace_root, manifest) == []
    assert integrity_level(manifest, [])[0] == "write_time_verified"
    assert integrity_level(manifest, [], locked=True)[0] == "locked"


def test_index_swap_is_detected(tmp_path: Path):
    trace_root = configure_llm_tracing(tmp_path / "traces", run_id="idx")
    _one_record(trace_root)
    manifest = build_trace_manifest(trace_root, run_id="idx")

    (trace_root / TRACE_INDEX_NAME).write_text("", encoding="utf-8")

    errors = validate_trace_manifest(trace_root, manifest)
    assert any("changed after manifest" in error for error in errors)


def test_lock_targets_cover_the_write_time_index(tmp_path: Path):
    trace_root = configure_llm_tracing(tmp_path / "traces", run_id="lock")
    _one_record(trace_root)
    build_trace_manifest(trace_root, run_id="lock")

    names = {path.name for path in trace_lock_targets(trace_root)}
    assert TRACE_INDEX_NAME in names
    assert TRACE_MANIFEST_NAME in names


def test_manifest_refuses_to_relabel_another_runs_records(tmp_path: Path):
    """The skip-extract defect: rebuilding a prior run's manifest under a new id."""
    trace_root = configure_llm_tracing(tmp_path / "traces", run_id="run-A")
    _one_record(trace_root)

    with pytest.raises(TraceRunMismatchError, match="run-A"):
        build_trace_manifest(trace_root, run_id="run-B")

    mixed = build_trace_manifest(trace_root, run_id="run-B", allow_mixed_runs=True)
    assert "run-A" in mixed["record_run_ids"]
    assert any(
        error.startswith("mixed_run_trace_records")
        for error in mixed["verification"]["errors"]
    )


# ---------------------------------------------------------------------------
# Accepted-response selection. Every one of these fails against the previous
# implementation, which let the LAST parsing call overwrite the accepted link.
# ---------------------------------------------------------------------------


def _fake_call(text: str):
    return SimpleNamespace(
        choices=[
            SimpleNamespace(message=SimpleNamespace(content=text), finish_reason="stop")
        ]
    )


def test_attempt_numbering_is_monotonic_across_retry_reentry(tmp_path: Path):
    """@llm_retry re-enters call_llm_json; a method-local counter restarted at 1."""
    from utils.llm_utils import BaseLLMCaller

    configure_llm_tracing(tmp_path / "traces", run_id="attempts")
    caller = BaseLLMCaller(model="fixture", max_tokens=64)
    calls = {"n": 0}

    def flaky(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise TimeoutError("transient")
        return _fake_call('{"ok": true}')

    with llm_attempt_ledger() as ledger:
        # Simulate tenacity re-entry: the same logical call, twice through the body.
        with pytest.raises(TimeoutError):
            caller._traced_call("primary", lambda: flaky())
        response, trace_id = caller._traced_call("primary", lambda: flaky())
        note_llm_outcome(trace_id, "parsed")

    numbers = [entry["attempt"] for entry in ledger["attempts"]]
    assert numbers == [1, 2], f"attempt numbers restarted: {numbers}"
    assert response is not None


def test_multi_model_selection_links_the_retained_call_not_the_last(tmp_path: Path):
    """Two successful models; the FIRST wins on yield. The link must follow yield."""
    configure_llm_tracing(tmp_path / "traces", run_id="selection")
    ids: dict[str, str] = {}

    with llm_attempt_ledger() as ledger:
        marks = {}
        for model in ("model-a", "model-b"):
            marks[model] = ledger_mark()
            _, summary = capture_llm_call(
                provider="fixture",
                requested_model=model,
                resolved_model=model,
                request={"model": model, "input": "extract"},
                call=lambda: SimpleNamespace(output_text="{}"),
            )
            trace_id = note_llm_attempt(summary, role="primary")
            note_llm_outcome(trace_id, "parsed")
            ids[model] = trace_id
        # model-a produced 12 variants, model-b only 3 -> keep model-a.
        finalize_attempt_selection(ledger_trace_ids_since(marks["model-a"])[:1])
        links = attempt_link_summary()

    assert links["accepted_response_trace_id"] == ids["model-a"]
    assert links["accepted_response_trace_ids"] == [ids["model-a"]]
    # The later successful call must be marked discarded, not left as accepted.
    assert links["discarded_trace_ids"] == [ids["model-b"]]
    outcomes = {
        link["trace_id"]: link["outcome"] for link in links["attempt_trace_links"]
    }
    assert outcomes[ids["model-a"]] == "accepted"
    assert outcomes[ids["model-b"]] == "discarded"


def test_repair_that_produced_the_data_is_the_accepted_link(tmp_path: Path):
    configure_llm_tracing(tmp_path / "traces", run_id="repair")
    with llm_attempt_ledger():
        _, primary = capture_llm_call(
            provider="fixture",
            requested_model="m",
            resolved_model="m",
            request={"input": "a"},
            call=lambda: SimpleNamespace(output_text="not json"),
        )
        primary_id = note_llm_attempt(primary, role="primary")
        note_llm_outcome(primary_id, "parse_failed")
        _, repair = capture_llm_call(
            provider="fixture",
            requested_model="m",
            resolved_model="m",
            request={"input": "repair"},
            call=lambda: SimpleNamespace(output_text="{}"),
        )
        repair_id = note_llm_attempt(repair, role="json_repair")
        note_llm_outcome(repair_id, "parsed", repaired=True)
        links = attempt_link_summary()

    assert links["accepted_response_trace_id"] == repair_id
    assert links["repaired"] is True
    # A parse failure is a FAILURE, not a discarded-but-usable response.
    assert links["failed_trace_ids"] == [primary_id]
    assert links["discarded_trace_ids"] == []
    assert links["parse_failures"] == 1


def test_single_call_stage_still_reports_an_accepted_link(tmp_path: Path):
    """Filters/triage have no selection step; the one parsed call is accepted."""
    configure_llm_tracing(tmp_path / "traces", run_id="single")
    with llm_attempt_ledger():
        _, summary = capture_llm_call(
            provider="fixture",
            requested_model="m",
            resolved_model="m",
            request={"input": "classify"},
            call=lambda: SimpleNamespace(output_text="{}"),
        )
        trace_id = note_llm_attempt(summary, role="primary")
        note_llm_outcome(trace_id, "parsed")
        links = attempt_link_summary()
    assert links["accepted_response_trace_id"] == trace_id
    assert links["discarded_trace_ids"] == []


def test_all_calls_failing_reports_no_accepted_response(tmp_path: Path):
    configure_llm_tracing(tmp_path / "traces", run_id="allfail")
    with llm_attempt_ledger():
        with pytest.raises(RuntimeError):
            capture_llm_call(
                provider="fixture",
                requested_model="m",
                resolved_model="m",
                request={"input": "x"},
                call=lambda: (_ for _ in ()).throw(RuntimeError("down")),
            )
        trace_id = note_llm_attempt(last_llm_trace(), role="primary")
        note_llm_outcome(trace_id, "error")
        links = attempt_link_summary()
    assert links["accepted_response_trace_id"] is None
    assert links["accepted_response_trace_ids"] == []
    assert links["failed_trace_ids"] == [trace_id]


# ---------------------------------------------------------------------------
# GVF_LLM_TRACE_DIR is a storage BASE, not one run's directory.
# ---------------------------------------------------------------------------


def test_trace_dir_override_gets_a_per_run_child(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("GVF_LLM_TRACE_DIR", str(tmp_path / "vol"))
    monkeypatch.delenv("GVF_LLM_TRACE_RUN_ID", raising=False)

    first = resolve_trace_location("run/one", default_root=tmp_path / "ignored")
    second = resolve_trace_location("run:two", default_root=tmp_path / "ignored")

    assert first.root == tmp_path / "vol" / "run_one"
    assert second.root == tmp_path / "vol" / "run_two"
    assert first.root != second.root, "sequential runs would share a trace tree"
    assert first.newly_selected and second.newly_selected


def test_an_already_selected_child_is_used_verbatim(tmp_path: Path, monkeypatch):
    """No double-append: a nested/subprocess stage must not add another level."""
    child = tmp_path / "vol" / "gvfrun-1"
    monkeypatch.setenv("GVF_LLM_TRACE_DIR", str(child))
    monkeypatch.setenv("GVF_LLM_TRACE_RUN_ID", "gvfrun-1")

    location = resolve_trace_location("some-other-id", default_root=tmp_path / "d")

    assert location.root == child
    assert location.run_id == "gvfrun-1"
    assert location.newly_selected is False
    # Resolving again (a further nesting level) is still the same directory.
    assert resolve_trace_location("x", default_root=tmp_path / "d").root == child


def test_exported_trace_identity_restores_previous_values(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("GVF_LLM_TRACE_DIR", "/original")
    monkeypatch.delenv("GVF_LLM_TRACE_RUN_ID", raising=False)

    with exported_trace_identity(tmp_path / "child", "run-x"):
        assert os.environ["GVF_LLM_TRACE_DIR"] == str(tmp_path / "child")
        assert os.environ["GVF_LLM_TRACE_RUN_ID"] == "run-x"

    assert os.environ["GVF_LLM_TRACE_DIR"] == "/original"
    assert "GVF_LLM_TRACE_RUN_ID" not in os.environ


def test_two_sequential_runs_under_one_base_do_not_trip_run_mismatch(
    tmp_path: Path, monkeypatch
):
    """The end-to-end symptom: mixed records made every later rebuild raise."""
    base = tmp_path / "vol"
    monkeypatch.setenv("GVF_LLM_TRACE_DIR", str(base))
    for run in ("run-1", "run-2"):
        monkeypatch.delenv("GVF_LLM_TRACE_RUN_ID", raising=False)
        location = resolve_trace_location(run, default_root=tmp_path / "unused")
        configure_llm_tracing(location.root, run_id=location.run_id)
        with llm_trace_scope(gene="KCNH2", pmid="111", stage="count_recovery"):
            record_trace_event("count_recovery_decision", {"counts_accepted": 1})
        # Must not raise: each run owns its own child.
        manifest = build_trace_manifest(location.root, run_id=location.run_id)
        assert manifest["record_run_ids"] == [run]
    assert sorted(p.name for p in base.iterdir()) == ["run-1", "run-2"]
