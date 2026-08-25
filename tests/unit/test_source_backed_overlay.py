from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).parents[2]
    / "benchmarks"
    / "codex_paper_eval"
    / "apply_source_backed_overlay.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("source_backed_overlay", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_overlay_sets_counts_and_adds_a_new_variant():
    module = load_module()
    predictions = {
        "papers": [
            {
                "gene": "SCN5A",
                "pmid": "1",
                "variants": [
                    {
                        "variant": "p.Arg1His",
                        "carriers": 1,
                        "affected": None,
                        "unaffected": None,
                    }
                ],
            }
        ]
    }
    overlay = {
        "operations": [
            {
                "operation": "set_counts",
                "gene": "SCN5A",
                "pmid": "1",
                "variant": "p.Arg1His",
                "carriers": 3,
                "affected": 2,
                "unaffected": 1,
                "source": "fixture",
            },
            {
                "operation": "add_variant",
                "gene": "SCN5A",
                "pmid": "1",
                "variant": "deletion of exon 3",
                "carriers": 2,
                "source": "fixture",
            },
        ]
    }

    module.apply_overlay(predictions, overlay)

    variants = predictions["papers"][0]["variants"]
    assert variants[0]["carriers"] == 3
    assert variants[0]["affected"] == 2
    assert variants[0]["unaffected"] == 1
    assert variants[1]["variant"] == "deletion of exon 3"
    assert variants[1]["source_layer"] == "source_backed_diagnostic"


def test_overlay_fails_closed_on_ambiguous_or_duplicate_identity():
    module = load_module()
    predictions = {
        "papers": [
            {
                "gene": "RYR2",
                "pmid": "1",
                "variants": [{"variant": "p.Arg1His"}, {"variant": "p.Arg1His"}],
            }
        ]
    }

    with pytest.raises(ValueError, match="requires one exact variant match"):
        module.apply_overlay(
            predictions,
            {
                "operations": [
                    {
                        "operation": "set_counts",
                        "gene": "RYR2",
                        "pmid": "1",
                        "variant": "p.Arg1His",
                        "carriers": 1,
                    }
                ]
            },
        )

    predictions["papers"][0]["variants"] = [{"variant": "p.Arg1His"}]
    with pytest.raises(ValueError, match="would duplicate"):
        module.apply_overlay(
            predictions,
            {
                "operations": [
                    {
                        "operation": "add_variant",
                        "gene": "RYR2",
                        "pmid": "1",
                        "variant": "p.Arg1His",
                    }
                ]
            },
        )


def test_checked_in_overlay_is_bound_to_current_baseline():
    module = load_module()
    repo = Path(__file__).parents[2]
    predictions = (
        repo
        / "benchmarks/codex_paper_eval/runs/20260824_postfix_gold118/predictions.json"
    )
    overlay_path = (
        repo
        / "benchmarks/codex_paper_eval/runs/20260824_postfix_gold118/diagnostics"
        / "source_backed_overlay_20260824.json"
    )
    overlay = json.loads(overlay_path.read_text())

    assert module.digest(predictions) == overlay["baseline_predictions_sha256"]
