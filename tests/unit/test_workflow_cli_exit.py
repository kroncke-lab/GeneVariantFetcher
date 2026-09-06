"""Workflow failures returned as data must fail both public CLI entry points."""

import sys

import pytest
from typer.testing import CliRunner

import cli
import cli.automated_workflow as workflow


@pytest.mark.parametrize("entrypoint", ["typer", "module"])
@pytest.mark.parametrize("success", [True, False])
def test_workflow_result_controls_cli_exit(monkeypatch, tmp_path, entrypoint, success):
    result = {"success": success, "error": "abstract retrieval failed"}
    monkeypatch.setenv("MODEL_PROVIDER", "azure")
    monkeypatch.setenv("AZURE_AI_API_KEY", "offline-test-key")
    monkeypatch.setattr(cli, "initialize_runtime", lambda: None)
    monkeypatch.setattr(workflow, "initialize_runtime", lambda: None)
    monkeypatch.setattr(cli, "has_llm_provider_key", lambda: True)
    monkeypatch.setattr(workflow, "has_llm_provider_key", lambda: True)
    monkeypatch.setattr(
        cli, "automated_variant_extraction_workflow", lambda **kw: result
    )
    monkeypatch.setattr(
        workflow, "automated_variant_extraction_workflow", lambda **kw: result
    )
    args = ["SCN5A", "--email", "ci@ncbi.test", "--output", str(tmp_path)]
    if entrypoint == "typer":
        outcome = CliRunner().invoke(cli.app, ["extract", *args])
        assert outcome.exit_code == (0 if success else 1), outcome.output
        if not success:
            assert "abstract retrieval failed" in outcome.output
    else:
        monkeypatch.setattr(
            sys, "argv", ["automated_workflow", *args, "--tier-threshold", "1"]
        )
        with pytest.raises(SystemExit) as exc:
            workflow.main()
        assert exc.value.code == (0 if success else 1)
