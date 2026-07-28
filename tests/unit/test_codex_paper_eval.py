"""Regression coverage for the extraction-blinded Codex paper evaluation."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

from benchmarks.codex_paper_eval.build_report_artifact import build_payload
from pipeline.prompts import TABLE_ATTRIBUTION_GUIDANCE
from utils.llm_trace import reset_llm_tracing
from benchmarks.codex_paper_eval.run_eval import (
    EXTRACTION_INSTRUCTIONS,
    choose_source,
    command_extract,
    command_lock,
    digest,
    effective_effort,
    matches,
    material_digest_errors,
    reasoning_params,
    looks_truncated_json,
    selection_metadata,
    supports_images,
    usable_sources,
    write_json,
    write_markdown_report,
)


COUNT_FIELDS = ("carriers", "affected", "unaffected")
GENES = ("SCN5A", "KCNH2", "KCNQ1", "RYR2")


@pytest.fixture(autouse=True)
def _reset_trace_configuration():
    reset_llm_tracing()
    yield
    reset_llm_tracing()


def _count_metric(asserted: int = 1, predicted: int = 1) -> dict:
    return {
        "gold_asserted": asserted,
        "predicted": predicted,
        "recall": predicted / asserted if asserted else None,
        "mae": 0.0 if predicted else None,
        "rmse": 0.0 if predicted else None,
    }


def _report_fixture() -> dict:
    papers = []
    by_gene = {}
    for index, gene in enumerate(GENES, 1):
        count = {field: _count_metric() for field in COUNT_FIELDS}
        paper = {
            "gene": gene,
            "pmid": str(index),
            "tool": "text",
            "tool_rationale": "Running text contained the evidence.",
            "source_completeness": "full_text",
            "elapsed_seconds": 1.0,
            "token_usage": {
                "input_tokens": 20,
                "output_tokens": 5,
                "total_tokens": 25,
            },
            "tp": 1,
            "fp": 0,
            "fn": 0,
            "precision": 1.0,
            "recall": 1.0,
            "f1": 1.0,
            "count": count,
            "matched_variants": [{"predicted": "A1V", "gold": "A1V"}],
            "missed_gold": [],
            "extra_predictions": [],
            "count_errors": [],
        }
        papers.append(paper)
        by_gene[gene] = {
            "papers": 1,
            "tp": 1,
            "fp": 0,
            "fn": 0,
            "precision": 1.0,
            "recall": 1.0,
            "f1": 1.0,
            "elapsed_seconds": 1.0,
            "token_usage": {
                "input_tokens": 20,
                "output_tokens": 5,
                "total_tokens": 25,
            },
            "count": count,
        }
    overall = {
        "papers": 4,
        "tp": 4,
        "fp": 0,
        "fn": 0,
        "precision": 1.0,
        "recall": 1.0,
        "f1": 1.0,
        "elapsed_seconds": 4.0,
        "token_usage": {
            "input_tokens": 80,
            "output_tokens": 20,
            "total_tokens": 100,
        },
        "count": {field: _count_metric(4, 4) for field in COUNT_FIELDS},
    }
    return {
        "run_id": "fixture",
        "seed": 7,
        "locked_at": "2026-07-24T00:00:00+00:00",
        "scored_at": "2026-07-24T00:01:00+00:00",
        "overall": overall,
        "by_gene": by_gene,
        "papers": papers,
        "selection": {
            "mode": "manifest",
            "population": "fixed manifest `fixture.tsv` (4 papers)",
            "description": (
                "Paper selection used the fixed manifest `fixture.tsv` (4 papers). "
                "Routing and extraction were gold-value-blind."
            ),
        },
        "tools_used": {"text": 4},
        "token_usage": {
            "input_tokens": 100,
            "output_tokens": 23,
            "total_tokens": 123,
        },
        "timing": {
            "wall_seconds": 60.0,
            "summed_paper_seconds": 4.0,
            "started_at": "2026-07-24T00:00:00+00:00",
            "completed_at": "2026-07-24T00:01:00+00:00",
        },
    }


def test_deletion_range_does_not_match_single_residue_deletion():
    assert not matches("K1505del", "K1505_Q1507del", "SCN5A")
    assert not matches("K1505_Q1507del", "K1505del", "SCN5A")
    assert matches("K1505_Q1507del", "K1505_Q1507del", "SCN5A")


def test_material_digests_cover_source_artifact_pdf_and_figure(tmp_path: Path):
    source = tmp_path / "source.md"
    artifact = tmp_path / "artifact.json"
    pdf = tmp_path / "paper.pdf"
    figure = tmp_path / "figure.png"
    source.write_text("source")
    artifact.write_text("{}")
    pdf.write_bytes(b"%PDF fixture")
    figure.write_bytes(b"PNG fixture")
    paper = {
        "gene": "SCN5A",
        "pmid": "1",
        "source": str(source),
        "source_sha256": digest(source),
        "artifacts": str(artifact),
        "artifacts_sha256": digest(artifact),
        "pdfs": [str(pdf)],
        "pdf_sha256": {str(pdf): digest(pdf)},
        "figures": [str(figure)],
        "figure_sha256": {str(figure): digest(figure)},
    }

    assert material_digest_errors(paper) == []

    figure.write_bytes(b"changed")
    assert any(
        "figure changed after selection" in error
        for error in material_digest_errors(paper)
    )


def test_api_usage_is_checkpointed_before_response_parsing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    source = tmp_path / "source.md"
    source.write_text("SCN5A A1V carrier evidence")
    selection = {
        "papers": [
            {
                "gene": "SCN5A",
                "pmid": "1",
                "source": str(source),
                "source_sha256": digest(source),
                "artifacts": None,
                "artifacts_sha256": None,
                "pdfs": [],
                "pdf_sha256": {},
                "figures": [],
                "figure_sha256": {},
            }
        ]
    }
    predictions = {
        "token_usage": None,
        "papers": [
            {
                "gene": "SCN5A",
                "pmid": "1",
                "tool": None,
                "token_usage": None,
                "variants": [],
            }
        ],
    }
    write_json(tmp_path / "selection.json", selection)
    write_json(tmp_path / "predictions.json", predictions)

    responses = iter(
        [
            SimpleNamespace(
                output_text=json.dumps(
                    {
                        "tool": "text",
                        "tool_rationale": "Only text was available.",
                        "source_completeness": "full_text",
                    }
                ),
                usage=SimpleNamespace(input_tokens=5, output_tokens=1),
            ),
            SimpleNamespace(
                output_text="not valid JSON",
                usage=SimpleNamespace(input_tokens=7, output_tokens=2),
            ),
        ]
    )
    fake_openai = ModuleType("openai")
    fake_openai.OpenAI = lambda **_kwargs: SimpleNamespace(
        responses=SimpleNamespace(create=lambda **_kwargs: next(responses))
    )
    monkeypatch.setitem(sys.modules, "openai", fake_openai)
    monkeypatch.setenv("AZURE_AI_API_BASE", "https://example.invalid")
    monkeypatch.setenv("AZURE_AI_API_KEY", "test-key")
    args = SimpleNamespace(
        run_dir=tmp_path,
        timeout=1,
        model="test-model",
        force=False,
        max_artifact_chars=100,
        max_source_chars=1000,
        route_preview_chars=500,
        max_ocr_images=1,
        route_reasoning_effort="medium",
        reasoning_effort="high",
        max_output_tokens=100,
        route_max_output_tokens=1600,
        max_output_tokens_ceiling=100000,
        legacy_table_material=False,
    )

    with pytest.raises(json.JSONDecodeError):
        command_extract(args)

    checkpoint = json.loads((tmp_path / "predictions.json").read_text())
    assert checkpoint["token_usage"]["total_tokens"] == 15
    assert checkpoint["papers"][0]["token_usage"]["total_tokens"] == 15
    refs = checkpoint["papers"][0]["llm_trace_refs"]
    assert [(ref["context"]["stage"], ref["success"]) for ref in refs] == [
        ("representation_route", True),
        ("representation_route_decision", None),
        ("paper_curation", True),
    ]
    trace_files = list((tmp_path / "llm_traces" / "SCN5A" / "1").glob("*.json"))
    assert len(trace_files) == 3


def test_selection_metadata_describes_manifest_and_random_modes(tmp_path: Path):
    manifest = tmp_path / "comparison.tsv"
    manifest.write_text("SCN5A\t1\n")
    fixed = selection_metadata(
        {
            "paper_manifest": str(manifest),
            "papers": [{"gene": "SCN5A", "pmid": "1"}],
            "seed": 7,
            "per_gene": 5,
        }
    )
    random = selection_metadata(
        {
            "paper_manifest": None,
            "papers": [{"gene": "SCN5A", "pmid": "1"}] * 4,
            "seed": 11,
            "per_gene": 1,
        }
    )

    assert fixed["mode"] == "manifest"
    assert "comparison.tsv" in fixed["description"]
    assert random["mode"] == "random"
    assert "seed 11" in random["description"]
    assert "high-carrier" not in fixed["description"] + random["description"]


def test_markdown_and_artifact_narratives_are_derived_from_report(tmp_path: Path):
    report = _report_fixture()
    write_json(tmp_path / "report.json", report)

    markdown_path = tmp_path / "report.md"
    write_markdown_report(report, markdown_path)
    markdown = markdown_path.read_text()
    assert "fixture.tsv" in markdown
    assert "gold-ranked" not in markdown
    assert "high-carrier" not in markdown

    payload = build_payload(tmp_path)
    manifest = payload["manifest"]
    executive = next(
        block["body"]
        for block in manifest["blocks"]
        if block["id"] == "executive_summary"
    )
    assert manifest["title"] == "Codex 4-Paper Extraction Evaluation"
    assert "100.0% precision" in executive
    assert "4 gold variant rows" in executive
    assert "123 exact API tokens" in executive
    assert "93.1%" not in json.dumps(payload)


def _rendering(path: Path, rows: int, variants: int, padding: int = 0) -> Path:
    """Write a candidate rendering with a known number of table rows/variants."""
    lines = [f"| {chr(65 + i % 26)}{i + 1}V | {i} | carrier |" for i in range(rows)]
    lines += [f"prose mentioning p.Ala{i + 1}Val here" for i in range(variants)]
    lines.append("x" * padding)
    path.write_text("\n".join(lines))
    return path


def test_choose_source_prefers_a_strictly_richer_rendering(tmp_path: Path):
    full_context = _rendering(tmp_path / "a_FULL_CONTEXT.md", rows=0, variants=1)
    cleaned = _rendering(tmp_path / "a_CLEANED.md", rows=12, variants=9)

    assert choose_source([full_context, cleaned]) == cleaned


def test_choose_source_keeps_priority_when_neither_rendering_dominates(tmp_path: Path):
    """More table rows but less prose is a trade, not an improvement."""
    full_context = _rendering(
        tmp_path / "b_FULL_CONTEXT.md", rows=8, variants=9, padding=5000
    )
    cleaned = _rendering(tmp_path / "b_CLEANED.md", rows=12, variants=9)

    assert choose_source([full_context, cleaned]) == full_context


def test_usable_sources_selects_the_richer_rendering_on_disk(tmp_path: Path):
    paper_dir = tmp_path / "KCNQ1" / "17470695"
    paper_dir.mkdir(parents=True)
    _rendering(paper_dir / "17470695_FULL_CONTEXT.md", rows=0, variants=1, padding=3000)
    cleaned = _rendering(
        paper_dir / "17470695_CLEANED.md", rows=40, variants=30, padding=3000
    )

    papers = usable_sources(tmp_path, "KCNQ1", minimum_chars=100)

    assert [Path(p["source"]) for p in papers] == [cleaned.resolve()]
    audit = papers[0]["source_selection"]
    assert audit["policy"] == "pareto_richness"
    assert audit["selected"] == str(cleaned.resolve())
    assert len(audit["candidates"]) == 2
    assert "Pareto-dominated" in audit["rationale"]


def _extract_captured_prompts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    source_text: str,
    artifact: str,
    route_tool: str = "table",
) -> tuple[dict, list]:
    source = tmp_path / "source.md"
    source.write_text(source_text)
    artifacts = tmp_path / "artifacts.json"
    artifacts.write_text(artifact)
    write_json(
        tmp_path / "selection.json",
        {
            "papers": [
                {
                    "gene": "KCNQ1",
                    "pmid": "17470695",
                    "source": str(source),
                    "source_sha256": digest(source),
                    "artifacts": str(artifacts),
                    "artifacts_sha256": digest(artifacts),
                    "pdfs": [],
                    "pdf_sha256": {},
                    "figures": [],
                    "figure_sha256": {},
                }
            ]
        },
    )
    write_json(
        tmp_path / "predictions.json",
        {
            "schema_version": 2,
            "token_usage": None,
            "papers": [
                {
                    "gene": "KCNQ1",
                    "pmid": "17470695",
                    "tool": None,
                    "token_usage": None,
                    "variants": [],
                }
            ],
        },
    )

    prompts: list = []
    responses = iter(
        [
            SimpleNamespace(
                output_text=json.dumps(
                    {
                        "tool": route_tool,
                        "tool_rationale": "Captions imply variant tables.",
                        "source_completeness": "partial_text",
                    }
                ),
                usage=SimpleNamespace(input_tokens=5, output_tokens=1),
            ),
            SimpleNamespace(
                output_text=json.dumps(
                    {
                        "notes": "No in-scope rows.",
                        "curation_rationale": (
                            "The supplied material had no qualifying human "
                            "variant/count row."
                        ),
                        "variants": [],
                    }
                ),
                usage=SimpleNamespace(input_tokens=7, output_tokens=2),
            ),
        ]
    )

    def create(**kwargs):
        prompts.append(kwargs["input"])
        return next(responses)

    fake_openai = ModuleType("openai")
    fake_openai.OpenAI = lambda **_kwargs: SimpleNamespace(
        responses=SimpleNamespace(create=create)
    )
    monkeypatch.setitem(sys.modules, "openai", fake_openai)
    monkeypatch.setenv("AZURE_AI_API_BASE", "https://example.invalid")
    monkeypatch.setenv("AZURE_AI_API_KEY", "test-key")
    command_extract(
        SimpleNamespace(
            run_dir=tmp_path,
            timeout=1,
            model="test-model",
            force=False,
            max_artifact_chars=1000,
            max_source_chars=10000,
            route_preview_chars=5000,
            max_ocr_images=1,
            route_reasoning_effort="medium",
            reasoning_effort="high",
            max_output_tokens=100,
            route_max_output_tokens=1600,
            max_output_tokens_ceiling=100000,
            legacy_table_material=False,
        )
    )
    return json.loads((tmp_path / "predictions.json").read_text()), prompts


def test_successful_eval_has_four_hash_lockable_trace_stages(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    checkpoint, _ = _extract_captured_prompts(
        tmp_path,
        monkeypatch,
        source_text=(
            "| Variant | carriers |\n|---|---|\n| A1V | 2 |\nKCNQ1 patient evidence."
        ),
        artifact="{}",
        route_tool="table",
    )

    refs = checkpoint["papers"][0]["llm_trace_refs"]
    assert [ref["context"]["stage"] for ref in refs] == [
        "representation_route",
        "representation_route_decision",
        "paper_curation",
        "paper_curation_decision",
    ]
    command_lock(SimpleNamespace(run_dir=tmp_path))
    lock = json.loads((tmp_path / "LOCK.json").read_text())
    assert lock["llm_trace_manifest_sha256"]
    assert lock["llm_trace_report_sha256"]
    manifest = json.loads((tmp_path / "llm_traces" / "trace_manifest.json").read_text())
    assert manifest["llm_call_count"] == 2
    assert manifest["decision_event_count"] == 2
    report = tmp_path / "llm_trace_report.html"
    assert report.is_file()
    assert "Sent to model" in report.read_text(encoding="utf-8")
    assert report.stat().st_mode & 0o222 == 0


def test_table_route_is_not_offered_without_real_table_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """An artifact full of captions must not advertise a table representation.

    Regression for KCNQ1 17470695: the artifact wrapper made the table string
    non-empty, the router picked ``table``, and the payload carried captions plus a
    keyword preview instead of any table rows.
    """
    checkpoint, prompts = _extract_captured_prompts(
        tmp_path,
        monkeypatch,
        source_text="KCNQ1 running text with no table rows at all.",
        artifact='{"table_captions_count": 5, "main_text": "TABLE 2. Mutations"}',
    )
    paper = checkpoint["papers"][0]

    assert paper["representations_available"] == ["text"]
    assert paper["tool"] == "text"
    assert "which was unavailable" in paper["tool_rationale"]
    assert "## TABLE PREVIEW\n\n" in prompts[0]
    assert "### Full/partial running text" in prompts[1]
    assert "no table rows at all" in prompts[1]


def test_parsed_artifact_reaches_the_model_on_a_non_table_route(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """The artifact used to ride along only on the table route."""
    _, prompts = _extract_captured_prompts(
        tmp_path,
        monkeypatch,
        source_text="KCNQ1 running text with no table rows at all.",
        artifact='{"supplements": ["Supplemental_Material.docx"]}',
    )

    assert "### Parsed artifact data" in prompts[1]
    assert "Supplemental_Material.docx" in prompts[1]
    assert "## PARSED ARTIFACT (attached to every route)" in prompts[0]


def test_grok_gets_no_reasoning_block_and_is_labelled_honestly():
    """Grok-4-class deployments 400 on reasoning.effort; they reason by default."""
    assert reasoning_params("gpt-5.6-sol", "high") == {"reasoning": {"effort": "high"}}
    assert effective_effort("gpt-5.6-sol", "high") == "high"

    assert reasoning_params("grok-4.3", "high") == {}
    assert effective_effort("grok-4.3", "high") == "model_default"


def test_ocr_route_is_withheld_from_a_model_without_vision():
    """grok-4.3 rejects image input; offering ocr would fail the paper outright."""
    assert supports_images("gpt-5.6-sol")
    assert not supports_images("grok-4.3")


def test_truncated_json_is_distinguished_from_garbage():
    """Azure reports status=completed even when output is cut off mid-string."""
    assert looks_truncated_json('{"variants": [{"variant": "A1V", "evidence": "row 1')
    assert looks_truncated_json('```json\n{"variants": [{"variant": "A1V"')
    assert not looks_truncated_json('{"variants": []}')
    assert not looks_truncated_json("not valid JSON")


def test_legacy_source_selection_ablation_restores_first_candidate(tmp_path: Path):
    """The ablation must reproduce the old behaviour exactly, or it proves nothing."""
    paper_dir = tmp_path / "KCNQ1" / "17470695"
    paper_dir.mkdir(parents=True)
    full = _rendering(
        paper_dir / "17470695_FULL_CONTEXT.md", rows=0, variants=1, padding=3000
    )
    _rendering(paper_dir / "17470695_CLEANED.md", rows=40, variants=30, padding=3000)

    legacy = usable_sources(
        tmp_path, "KCNQ1", minimum_chars=100, legacy_source_selection=True
    )
    fixed = usable_sources(tmp_path, "KCNQ1", minimum_chars=100)

    assert Path(legacy[0]["source"]) == full.resolve()
    assert Path(fixed[0]["source"]) != full.resolve()


def test_extraction_prompt_carries_table_attribution_guidance():
    """Compilation tables caused the 10973849 over-attribution regression."""
    # Pinned to the canonical copy so the harness cannot silently drift from it,
    # and so the guidance survives this harness being replaced.
    assert TABLE_ATTRIBUTION_GUIDANCE in EXTRACTION_INSTRUCTIONS
    assert "compilation citing other" in TABLE_ATTRIBUTION_GUIDANCE
    assert (
        "families, individuals, alleles, probands, cases" in TABLE_ATTRIBUTION_GUIDANCE
    )
    assert "Count only what this study observed" in TABLE_ATTRIBUTION_GUIDANCE
    assert "{" not in TABLE_ATTRIBUTION_GUIDANCE, "must be str.format-safe"
    # The measured dominant failure was omitting uncountable variants entirely.
    assert (
        "always\nemit the variant even when all three counts are null"
        in EXTRACTION_INSTRUCTIONS
    )


def test_table_route_carries_full_text_not_a_keyword_preview(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Provenance needs the caption/footnotes; a 6k excerpt strips them.

    Regression for KCNH2 10973849, where a compilation table read as first-party data.
    """
    marker = "SENTINEL_PROSE_FAR_FROM_ANY_VARIANT_MENTION"
    source = (
        "| Mutation | n | Type |\n| A1V | 2 | Missense |\n| G2R | 3 | Missense |\n"
        # Long enough to fall outside targeted_preview's reach, short enough to
        # survive max_source_chars, so the assertion tests the route not the cap.
        + ("filler prose. " * 400)
        + marker
    )
    _, prompts = _extract_captured_prompts(
        tmp_path, monkeypatch, source_text=source, artifact="{}", route_tool="table"
    )

    assert "### Structured table rows" in prompts[1]
    assert marker in prompts[1], "full running text must reach the table route"
