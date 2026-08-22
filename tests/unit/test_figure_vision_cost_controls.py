"""Vision-lane cost controls: dedup, caption triage, token-cap retry, usage.

Evidence base (docs/evidence/strategy_consult_20260822.md): 47/244 gold-120
figure-vision calls re-read byte-identical images; 3 burned the full
4096-token reasoning cap and returned empty; and the raw Responses-API path
recorded usage under counter names cost tooling does not sum, so the vision
lane billed as zero. No live calls — requests.post and the vision read are
mocked throughout.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

import harvesting.figure_text_extractor as fte
import harvesting.figure_variant_reader as fvr
from utils.llm_trace import (
    TRACE_MANIFEST_NAME,
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


def _events(trace_root: Path, event_type: str) -> list[dict]:
    return [
        record["event"]["data"]
        for record in _records(trace_root)
        if (record.get("event") or {}).get("type") == event_type
    ]


def _png(tmp_path: Path, name: str, payload: bytes = b"same-bytes") -> Path:
    path = tmp_path / name
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + payload)
    return path


def _patched_read_one(monkeypatch) -> list[str]:
    """Replace the vision call, returning the list of image names it saw."""
    calls: list[str] = []

    def fake_read_one(img, gene, model, *, pmid=""):
        calls.append(img.name)
        return fvr.FigureReadResult(
            image_path=str(img), variants=[{"protein": "R176W"}]
        )

    monkeypatch.setattr(fvr, "_read_one", fake_read_one)
    return calls


# ---------------------------------------------------------------------------
# SHA-256 image dedup
# ---------------------------------------------------------------------------


class TestImageDedup:
    def test_duplicate_image_is_read_once(self, tmp_path: Path, monkeypatch):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="dedup")
        calls = _patched_read_one(monkeypatch)
        first = _png(tmp_path, "fig1.png")
        duplicate = _png(tmp_path, "fig1_via_other_url.png")

        report = fvr.read_images(
            [first, duplicate], "KCNH2", pmid="12345678", model="anthropic/x"
        )

        assert calls == ["fig1.png"]
        assert report.per_figure[0].skipped is None
        assert report.per_figure[1].skipped == "duplicate_image:fig1.png"
        assert report.per_figure[1].variants == []
        # The paper still keeps the identity from the read that happened.
        assert report.distinct_variants == [{"protein": "R176W"}]
        (skip,) = _events(trace_root, "figure_vision_skip")
        assert skip["skip"] == "duplicate_image"
        assert skip["image"] == "fig1_via_other_url.png"
        assert skip["duplicate_of"] == "fig1.png"

    def test_distinct_images_are_both_read(self, tmp_path: Path, monkeypatch):
        calls = _patched_read_one(monkeypatch)
        a = _png(tmp_path, "a.png", payload=b"bytes-a")
        b = _png(tmp_path, "b.png", payload=b"bytes-b")

        report = fvr.read_images([a, b], "KCNH2", pmid="1", model="anthropic/x")

        assert calls == ["a.png", "b.png"]
        assert [r.skipped for r in report.per_figure] == [None, None]


# ---------------------------------------------------------------------------
# Caption triage gate
# ---------------------------------------------------------------------------


class TestCaptionTriage:
    def test_reads_image_with_no_caption(self, tmp_path: Path, monkeypatch):
        """Fail-open: an unmapped caption is unknown, not skippable."""
        calls = _patched_read_one(monkeypatch)
        img = _png(tmp_path, "fig_pmc_1.jpg")

        report = fvr.read_images(
            [img], "KCNH2", pmid="1", model="anthropic/x", captions={}
        )

        assert calls == ["fig_pmc_1.jpg"]
        assert report.per_figure[0].skipped is None

    def test_reads_image_whose_caption_matches_no_list(
        self, tmp_path: Path, monkeypatch
    ):
        calls = _patched_read_one(monkeypatch)
        img = _png(tmp_path, "fig2.png")

        fvr.read_images(
            [img],
            "KCNH2",
            pmid="1",
            model="anthropic/x",
            captions={"fig2.png": "Sequencing chromatograms for Fig. 2."},
        )

        assert calls == ["fig2.png"]

    def test_skips_functional_only_caption(self, tmp_path: Path, monkeypatch):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="triage")
        calls = _patched_read_one(monkeypatch)
        img = _png(tmp_path, "fig3.png")

        report = fvr.read_images(
            [img],
            "KCNH2",
            pmid="12345678",
            model="anthropic/x",
            captions={"fig3.png": "Representative current traces and Boltzmann fits."},
        )

        assert calls == []
        assert report.per_figure[0].skipped == (
            "caption_triage:functional:current trace"
        )
        (skip,) = _events(trace_root, "figure_vision_skip")
        assert skip["skip"] == "caption_triage"
        assert skip["image"] == "fig3.png"
        assert skip["matched_keyword"] == "functional:current trace"

    def test_keeps_pedigree_caption_despite_functional_words(
        self, tmp_path: Path, monkeypatch
    ):
        calls = _patched_read_one(monkeypatch)
        img = _png(tmp_path, "fig4.png")

        report = fvr.read_images(
            [img],
            "KCNH2",
            pmid="1",
            model="anthropic/x",
            captions={
                "fig4.png": (
                    "Pedigree of kindred K1 and voltage clamp recordings "
                    "of the mutant channel."
                )
            },
        )

        assert calls == ["fig4.png"]
        assert report.per_figure[0].skipped is None

    def test_off_switch_reads_functional_caption(self, tmp_path: Path, monkeypatch):
        monkeypatch.setenv(fvr.FIGURE_TRIAGE_ENV, "off")
        calls = _patched_read_one(monkeypatch)
        img = _png(tmp_path, "fig5.png")

        fvr.read_images(
            [img],
            "KCNH2",
            pmid="1",
            model="anthropic/x",
            captions={"fig5.png": "Western blot and topology schematic."},
        )

        assert calls == ["fig5.png"]


# ---------------------------------------------------------------------------
# Caption lookup
# ---------------------------------------------------------------------------


FULL_CONTEXT = """\
## FIGURE CAPTIONS

### Figure 1.

Pedigree of the proband's family showing affected carriers.

_image_: https://cdn.example.org/content/fig1_small.jpg?width=200
_image_: https://cdn.example.org/content/FIG1_LARGE.jpg

### Figure 2.

Representative current traces from voltage-clamp experiments.

_image_: mmc2.png

## TABLE CAPTIONS
"""


class TestCaptionLookup:
    def test_parser_maps_every_image_line_to_the_block_caption(self):
        captions = fvr.parse_full_context_captions(FULL_CONTEXT)
        # Both URLs of the same figure bind to one caption; basenames are
        # lowercased and stripped of query strings.
        assert "pedigree" in captions["fig1_small.jpg"].lower()
        assert captions["fig1_small.jpg"] == captions["fig1_large.jpg"]
        assert "current traces" in captions["mmc2.png"]

    def test_lookup_prefers_the_on_disk_sidecar_binding(self, tmp_path: Path):
        # fig_pmc_1.jpg exists only in captions.json — the FULL_CONTEXT
        # _image_ basename never matches PMC-renamed files on disk.
        fig_dir = tmp_path / "999_figures"
        fig_dir.mkdir()
        (fig_dir / "captions.json").write_text(
            json.dumps(
                [
                    {
                        "filename": "fig_pmc_1.jpg",
                        "label": "Figure 1",
                        "title": "Pedigree of the family.",
                        "text": "Filled symbols denote affected carriers.",
                    }
                ]
            )
        )
        (tmp_path / "999_FULL_CONTEXT.md").write_text(FULL_CONTEXT)

        captions = fvr.load_figure_captions(tmp_path, "999")

        assert "pedigree" in captions["fig_pmc_1.jpg"].lower()
        assert "mmc2.png" in captions  # FULL_CONTEXT names still contribute
        assert "fig1_small.jpg" in captions


# ---------------------------------------------------------------------------
# Token-cap burn: incomplete envelope retry
# ---------------------------------------------------------------------------


RESPONSES_USAGE = {
    "input_tokens": 528,
    "output_tokens": 4096,
    "output_tokens_details": {"reasoning_tokens": 4027},
}


def _envelope(
    text: "str | None",
    *,
    status: str = "completed",
    reason: "str | None" = None,
    usage: "dict | None" = None,
) -> dict:
    output: list[dict] = [{"type": "reasoning", "summary": []}]
    if text is not None:
        output.append(
            {
                "type": "message",
                "content": [{"type": "output_text", "text": text}],
            }
        )
    data = {"id": "resp_x", "object": "response", "status": status, "output": output}
    if reason is not None:
        data["incomplete_details"] = {"reason": reason}
    if usage is not None:
        data["usage"] = dict(usage)
    return data


def _fake_post(*envelopes: dict):
    served = iter(envelopes)
    calls: list[dict] = []

    def post(url, headers, json, timeout):  # noqa: A002 - matches call site
        calls.append({"url": url, "headers": headers, "body": json})
        body = next(served)
        return SimpleNamespace(
            status_code=200,
            text="{}",
            json=lambda: body,
        )

    return post, calls


@pytest.fixture()
def _azure_env(monkeypatch):
    monkeypatch.setenv("AZURE_AI_API_BASE", "https://example.services.ai.azure.com")
    monkeypatch.setenv("AZURE_AI_API_KEY", "test-key")


def _stub_effort(monkeypatch, effort: "str | None"):
    monkeypatch.setattr(
        fte, "get_settings", lambda: SimpleNamespace(vision_reasoning_effort=effort)
    )


class TestIncompleteEnvelopeRetry:
    def test_retries_once_at_low_effort_then_accepts_empty(
        self, tmp_path: Path, monkeypatch, _azure_env
    ):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="cap")
        _stub_effort(monkeypatch, None)
        burned = _envelope(
            None, status="incomplete", reason="max_output_tokens", usage=RESPONSES_USAGE
        )
        post, calls = _fake_post(burned, burned)

        with patch.object(fte.requests, "post", post):
            out = fte.call_responses_api_vision(
                "read this",
                "data:image/png;base64,AA==",
                "azure_ai/gpt-5.3-codex-1",
                attempt_role="figure_variant_read",
            )

        assert out == ""
        assert len(calls) == 2
        assert "reasoning" not in calls[0]["body"]  # provider-default effort
        assert calls[1]["body"]["reasoning"] == {"effort": "low"}
        outcomes = [
            e["outcome"] for e in _events(trace_root, "figure_vision_token_cap")
        ]
        assert outcomes == ["retry_low_effort", "accepted_empty_after_retry"]

    def test_retry_recovers_text(self, tmp_path: Path, monkeypatch, _azure_env):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="cap2")
        _stub_effort(monkeypatch, "medium")
        post, calls = _fake_post(
            _envelope(None, status="incomplete", reason="max_output_tokens"),
            _envelope('{"variants": []}'),
        )

        with patch.object(fte.requests, "post", post):
            out = fte.call_responses_api_vision(
                "read this",
                "data:image/png;base64,AA==",
                "azure_ai/gpt-5.3-codex-1",
                attempt_role="figure_variant_read",
            )

        assert out == '{"variants": []}'
        assert len(calls) == 2
        assert calls[0]["body"]["reasoning"] == {"effort": "medium"}
        assert calls[1]["body"]["reasoning"] == {"effort": "low"}
        outcomes = [
            e["outcome"] for e in _events(trace_root, "figure_vision_token_cap")
        ]
        assert outcomes == ["retry_low_effort"]

    def test_already_low_effort_accepts_empty_without_retry(
        self, tmp_path: Path, monkeypatch, _azure_env
    ):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="cap3")
        _stub_effort(monkeypatch, "low")
        post, calls = _fake_post(
            _envelope(None, status="incomplete", reason="max_output_tokens")
        )

        with patch.object(fte.requests, "post", post):
            out = fte.call_responses_api_vision(
                "read this",
                "data:image/png;base64,AA==",
                "azure_ai/gpt-5.3-codex-1",
                attempt_role="figure_ocr",
            )

        assert out == ""
        assert len(calls) == 1
        outcomes = [
            e["outcome"] for e in _events(trace_root, "figure_vision_token_cap")
        ]
        assert outcomes == ["accepted_empty_no_retry"]

    def test_clean_complete_envelope_never_retries(
        self, tmp_path: Path, monkeypatch, _azure_env
    ):
        _stub_effort(monkeypatch, None)
        post, calls = _fake_post(_envelope("some text"))

        with patch.object(fte.requests, "post", post):
            out = fte.call_responses_api_vision(
                "read this",
                "data:image/png;base64,AA==",
                "azure_ai/gpt-5.3-codex-1",
                attempt_role="figure_ocr",
            )

        assert out == "some text"
        assert len(calls) == 1


# ---------------------------------------------------------------------------
# Usage capture in the trace record
# ---------------------------------------------------------------------------


class TestResponsesUsageCapture:
    def test_usage_aliased_into_trace_record(
        self, tmp_path: Path, monkeypatch, _azure_env
    ):
        """End-to-end through read_images: the llm_call record carries both the
        provider's counters and the litellm names cost tooling sums."""
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="usage")
        _stub_effort(monkeypatch, None)
        usage = {
            "input_tokens": 528,
            "output_tokens": 100,
            "output_tokens_details": {"reasoning_tokens": 41},
        }
        post, _calls = _fake_post(
            _envelope('{"variants": [{"protein": "R176W"}]}', usage=usage)
        )
        img = _png(tmp_path, "fig1.png")

        with patch.object(fte.requests, "post", post):
            report = fvr.read_images(
                [img], "KCNH2", pmid="12345678", model="azure_ai/gpt-5.3-codex-1"
            )

        assert report.per_figure[0].variants
        (call,) = [
            r for r in _records(trace_root) if r.get("record_type") == "llm_call"
        ]
        recorded = call["response"]["usage"]
        assert recorded["prompt_tokens"] == 528
        assert recorded["completion_tokens"] == 100
        assert recorded["total_tokens"] == 628
        assert recorded["completion_tokens_details"]["reasoning_tokens"] == 41
        # Provider names stay alongside — nothing was renamed away.
        assert recorded["input_tokens"] == 528
        assert recorded["output_tokens"] == 100

    def test_missing_usage_stays_absent(self, tmp_path: Path, monkeypatch, _azure_env):
        trace_root = configure_llm_tracing(tmp_path / "traces", run_id="usage2")
        _stub_effort(monkeypatch, None)
        post, _calls = _fake_post(_envelope('{"variants": []}'))

        with patch.object(fte.requests, "post", post):
            fvr.read_images(
                [_png(tmp_path, "fig1.png")],
                "KCNH2",
                pmid="1",
                model="azure_ai/gpt-5.3-codex-1",
            )

        (call,) = [
            r for r in _records(trace_root) if r.get("record_type") == "llm_call"
        ]
        assert call["response"]["usage"] is None

    def test_normalize_never_fabricates_counters(self):
        data = {"usage": {"output_tokens_details": {"reasoning_tokens": 7}}}
        fte.normalize_responses_usage(data)
        usage = data["usage"]
        assert "prompt_tokens" not in usage
        assert "completion_tokens" not in usage
        assert "total_tokens" not in usage
        assert usage["completion_tokens_details"] == {"reasoning_tokens": 7}
