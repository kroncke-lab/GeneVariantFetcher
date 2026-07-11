"""Tests for the optional vision-LLM figure-text extraction path.

Covers two layers:
  1. ``harvesting.figure_text_extractor`` — the standalone helper module.
  2. The vision pass inside ``enrich_paywall_full_context`` — injectable
     extractor + GVF_EXTRACT_FIGURE_TEXT env gate.

No real API calls are made: litellm.completion is never imported from here.
"""

from __future__ import annotations

from pathlib import Path
from typing import List
from unittest.mock import MagicMock, patch

import pytest

from harvesting.figure_text_extractor import (
    _IMAGE_SUFFIXES,
    extract_images_to_markdown,
    is_image_path,
)
from harvesting.paywall_context_enrichment import enrich_paywall_full_context


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_dummy_image(path: Path) -> Path:
    """Write a 1×1 white PNG (valid enough to pass read_bytes)."""
    # Minimal valid PNG bytes (1x1 white pixel).
    PNG_1X1 = (
        b"\x89PNG\r\n\x1a\n"
        b"\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
        b"\x08\x02\x00\x00\x00\x90wS\xde\x00\x00\x00\x0cIDATx"
        b"\x9cc\xf8\x0f\x00\x00\x01\x01\x00\x05\x18\xd8N\x00\x00"
        b"\x00\x00IEND\xaeB`\x82"
    )
    path.write_bytes(PNG_1X1)
    return path


# ---------------------------------------------------------------------------
# figure_text_extractor module
# ---------------------------------------------------------------------------


class TestIsImagePath:
    def test_recognises_png(self):
        assert is_image_path(Path("fig1.png"))

    def test_recognises_jpg(self):
        assert is_image_path("figure.JPG")  # case-insensitive

    def test_rejects_pdf(self):
        assert not is_image_path("supplement.pdf")

    def test_rejects_xlsx(self):
        assert not is_image_path(Path("table.xlsx"))

    def test_all_registered_suffixes_pass(self):
        for suf in _IMAGE_SUFFIXES:
            assert is_image_path(f"img{suf}"), f"Expected {suf} to be recognised"


class TestExtractImagesToMarkdown:
    def test_empty_list_returns_empty_string_no_api_call(self):
        # Should return "" immediately — no litellm call.
        result = extract_images_to_markdown([], model="anthropic/claude-sonnet-4-6")
        assert result == ""

    def test_successful_extraction_returns_section(self, tmp_path: Path):
        img = _write_dummy_image(tmp_path / "fig1.png")

        fake_response = MagicMock()
        fake_response.choices[
            0
        ].message.content = "| Variant | Count |\n|---|---|\n| p.Arg176Trp | 4 |"

        with patch(
            "harvesting.figure_text_extractor.litellm_completion",
            return_value=fake_response,
        ):
            result = extract_images_to_markdown(
                [img], model="anthropic/claude-sonnet-4-6"
            )

        assert "## FIGURE IMAGE TEXT" in result
        assert "fig1.png" in result
        assert "p.Arg176Trp" in result

    def test_api_failure_skipped_gracefully(self, tmp_path: Path):
        img = _write_dummy_image(tmp_path / "fig2.png")

        with patch(
            "harvesting.figure_text_extractor.litellm_completion",
            side_effect=RuntimeError("timeout"),
        ):
            result = extract_images_to_markdown(
                [img], model="anthropic/claude-sonnet-4-6"
            )

        # No crash; empty result because the only image failed.
        assert result == ""

    def test_empty_model_response_omitted(self, tmp_path: Path):
        img = _write_dummy_image(tmp_path / "blank.png")

        fake_response = MagicMock()
        fake_response.choices[0].message.content = "   "  # whitespace only

        with patch(
            "harvesting.figure_text_extractor.litellm_completion",
            return_value=fake_response,
        ):
            result = extract_images_to_markdown(
                [img], model="anthropic/claude-sonnet-4-6"
            )

        assert result == ""

    def test_multiple_images_indexed(self, tmp_path: Path):
        imgs = [
            _write_dummy_image(tmp_path / "a.png"),
            _write_dummy_image(tmp_path / "b.jpg"),
        ]

        def _fake_completion(**kwargs):
            url = kwargs["messages"][0]["content"][1]["image_url"]["url"]
            r = MagicMock()
            r.choices[0].message.content = f"Text from {url[:20]}"
            return r

        with patch(
            "harvesting.figure_text_extractor.litellm_completion",
            side_effect=_fake_completion,
        ):
            result = extract_images_to_markdown(
                imgs, model="anthropic/claude-sonnet-4-6"
            )

        assert "Figure Image 1: a.png" in result
        assert "Figure Image 2: b.jpg" in result

    def test_azure_gpt5_uses_responses_api(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        img = _write_dummy_image(tmp_path / "azure.png")
        monkeypatch.setenv("AZURE_AI_API_BASE", "https://example.services.ai.azure.com")
        monkeypatch.setenv("AZURE_AI_API_KEY", "test-key")

        class _Resp:
            status_code = 200
            text = '{"ok": true}'

            def json(self):
                return {
                    "output": [
                        {
                            "type": "message",
                            "content": [
                                {
                                    "type": "output_text",
                                    "text": "| Variant | Notes |\n|---|---|\n| p.Glu95Gly | visible |",
                                }
                            ],
                        }
                    ]
                }

        captured = {}

        def _fake_post(url, headers, json, timeout):
            captured["url"] = url
            captured["headers"] = headers
            captured["body"] = json
            captured["timeout"] = timeout
            return _Resp()

        with (
            patch("harvesting.figure_text_extractor.requests.post", _fake_post),
            patch(
                "harvesting.figure_text_extractor.litellm_completion",
                side_effect=AssertionError("chat completion should not be used"),
            ),
        ):
            result = extract_images_to_markdown([img], model="azure_ai/gpt-5.3-codex-1")

        assert "p.Glu95Gly" in result
        assert captured["url"].endswith("/openai/v1/responses?api-version=v1")
        assert captured["headers"]["api-key"] == "test-key"
        assert captured["body"]["model"] == "gpt-5.3-codex-1"


# ---------------------------------------------------------------------------
# enrich_paywall_full_context — vision gate
# ---------------------------------------------------------------------------


class TestEnrichVisionGate:
    """The vision pass in enrich_paywall_full_context."""

    def _make_supplement_result_with_images(self, tmp_path: Path) -> "list":
        """Return a SupplementFileResult whose nested_files include a real image."""
        from harvesting.supplement_processing_service import SupplementFileResult

        img_path = _write_dummy_image(tmp_path / "extracted_fig1.png")
        sfr = SupplementFileResult(
            filename="supplement.pdf",
            path=str(tmp_path / "supplement.pdf"),
            downloaded=True,
            nested_files=[str(img_path)],
        )
        return [sfr]

    def test_injected_extractor_appends_text(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """injectable image_text_extractor result lands in unified_markdown."""
        import harvesting.paywall_context_enrichment as _mod

        sfr_list = self._make_supplement_result_with_images(tmp_path)

        captured_paths: List[List[Path]] = []

        def _fake_extractor(paths: List[Path]) -> str:
            captured_paths.append(list(paths))
            return "\n\n## FIGURE IMAGE TEXT\n\n### Figure Image 1: extracted_fig1.png\n\np.Gly604Ser variant table\n"

        # Monkeypatch process_supplement_files to return our pre-built results.
        def _fake_process(**kwargs):
            from harvesting.supplement_processing_service import (
                SupplementProcessingResult,
            )

            return SupplementProcessingResult(
                supplement_markdown="",
                downloaded_count=1,
                total_figures_extracted=1,
                file_results=sfr_list,
            )

        monkeypatch.setattr(_mod, "process_supplement_files", _fake_process)

        result = enrich_paywall_full_context(
            body_markdown="# BODY\n\nSome text.\n",
            html=None,
            supp_files=[{"url": "https://example.org/s.pdf", "name": "s.pdf"}],
            pmid="11111111",
            output_dir=tmp_path,
            converter=MagicMock(),
            session=MagicMock(
                **{
                    "get.return_value.__enter__": lambda s: s,
                    "get.return_value.__exit__": lambda s, *a: False,
                    "get.return_value.status_code": 200,
                    "get.return_value.iter_content.return_value": [b"data"],
                }
            ),
            image_text_extractor=_fake_extractor,
        )

        assert "p.Gly604Ser variant table" in result.unified_markdown
        assert "## FIGURE IMAGE TEXT" in result.unified_markdown
        # Extractor was called with the image path from nested_files.
        assert len(captured_paths) == 1
        assert any(p.name == "extracted_fig1.png" for p in captured_paths[0])

    def test_env_disabled_by_default_no_extractor_called(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """Without env var or injected extractor, vision is never invoked."""
        import harvesting.paywall_context_enrichment as _mod

        monkeypatch.delenv("GVF_EXTRACT_FIGURE_TEXT", raising=False)

        called: list = []

        def _extractor_should_not_run(paths):
            called.append(paths)
            return "SHOULD NOT APPEAR"

        # No injected extractor, env is unset — vision must not run.
        result = enrich_paywall_full_context(
            body_markdown="# BODY\n",
            html=None,
            supp_files=None,
            pmid="22222222",
            output_dir=tmp_path,
            converter=MagicMock(),
            session=None,
            download_supplements=False,
            # image_text_extractor intentionally omitted
        )

        assert "SHOULD NOT APPEAR" not in result.unified_markdown
        assert not called

    def test_env_var_enables_extraction(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """GVF_EXTRACT_FIGURE_TEXT=1 triggers the default extractor path."""
        import harvesting.paywall_context_enrichment as _mod

        monkeypatch.setenv("GVF_EXTRACT_FIGURE_TEXT", "1")

        img = _write_dummy_image(tmp_path / "env_fig.png")
        sfr_list = self._make_supplement_result_with_images(tmp_path)

        def _fake_process(**kwargs):
            from harvesting.supplement_processing_service import (
                SupplementProcessingResult,
            )

            return SupplementProcessingResult(
                supplement_markdown="",
                downloaded_count=1,
                total_figures_extracted=1,
                file_results=sfr_list,
            )

        monkeypatch.setattr(_mod, "process_supplement_files", _fake_process)

        captured: list = []

        def _fake_extract(paths, model):
            captured.append((paths, model))
            return "\n\n## FIGURE IMAGE TEXT\n\n### Figure Image 1: env_fig.png\n\nEnv-triggered extraction.\n"

        monkeypatch.setattr(_mod, "extract_images_to_markdown", _fake_extract)

        result = enrich_paywall_full_context(
            body_markdown="# BODY\n",
            html=None,
            supp_files=[{"url": "https://example.org/s.pdf", "name": "s.pdf"}],
            pmid="33333333",
            output_dir=tmp_path,
            converter=MagicMock(),
            session=MagicMock(
                **{
                    "get.return_value.__enter__": lambda s: s,
                    "get.return_value.__exit__": lambda s, *a: False,
                    "get.return_value.status_code": 200,
                    "get.return_value.iter_content.return_value": [b"data"],
                }
            ),
        )

        assert "Env-triggered extraction" in result.unified_markdown
        assert len(captured) == 1
        _, model_used = captured[0]
        # Model comes from settings.get_vision_model() — just verify it's a string.
        assert isinstance(model_used, str) and model_used

    def test_non_image_nested_files_ignored(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """PDF/xlsx nested_files must not be passed to the vision extractor."""
        import harvesting.paywall_context_enrichment as _mod
        from harvesting.supplement_processing_service import SupplementFileResult

        sfr = SupplementFileResult(
            filename="s.zip",
            path=str(tmp_path / "s.zip"),
            downloaded=True,
            nested_files=[
                str(tmp_path / "inner_table.xlsx"),
                str(tmp_path / "readme.txt"),
            ],
        )

        def _fake_process(**kwargs):
            from harvesting.supplement_processing_service import (
                SupplementProcessingResult,
            )

            return SupplementProcessingResult(
                supplement_markdown="",
                downloaded_count=1,
                total_figures_extracted=0,
                file_results=[sfr],
            )

        monkeypatch.setattr(_mod, "process_supplement_files", _fake_process)

        received_paths: list = []

        def _extractor(paths: List[Path]) -> str:
            received_paths.extend(paths)
            return ""

        enrich_paywall_full_context(
            body_markdown="# BODY\n",
            html=None,
            supp_files=[{"url": "https://example.org/s.zip", "name": "s.zip"}],
            pmid="44444444",
            output_dir=tmp_path,
            converter=MagicMock(),
            session=MagicMock(
                **{
                    "get.return_value.__enter__": lambda s: s,
                    "get.return_value.__exit__": lambda s, *a: False,
                    "get.return_value.status_code": 200,
                    "get.return_value.iter_content.return_value": [b"data"],
                }
            ),
            image_text_extractor=_extractor,
        )

        # No image-suffix files → extractor called with empty list → no paths forwarded.
        assert received_paths == []

    def test_vision_exception_does_not_abort_enrichment(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """A crash inside the extractor must not propagate — body is still returned."""
        import harvesting.paywall_context_enrichment as _mod

        img = _write_dummy_image(tmp_path / "crash_fig.png")
        sfr_list = self._make_supplement_result_with_images(tmp_path)

        def _fake_process(**kwargs):
            from harvesting.supplement_processing_service import (
                SupplementProcessingResult,
            )

            return SupplementProcessingResult(
                supplement_markdown="",
                downloaded_count=1,
                total_figures_extracted=1,
                file_results=sfr_list,
            )

        monkeypatch.setattr(_mod, "process_supplement_files", _fake_process)

        def _crashing_extractor(paths):
            raise RuntimeError("vision model unavailable")

        result = enrich_paywall_full_context(
            body_markdown="# BODY\n\nSafe content.\n",
            html=None,
            supp_files=[{"url": "https://example.org/s.pdf", "name": "s.pdf"}],
            pmid="55555555",
            output_dir=tmp_path,
            converter=MagicMock(),
            session=MagicMock(
                **{
                    "get.return_value.__enter__": lambda s: s,
                    "get.return_value.__exit__": lambda s, *a: False,
                    "get.return_value.status_code": 200,
                    "get.return_value.iter_content.return_value": [b"data"],
                }
            ),
            image_text_extractor=_crashing_extractor,
        )

        assert "Safe content" in result.unified_markdown
        assert "## FIGURE IMAGE TEXT" not in result.unified_markdown
