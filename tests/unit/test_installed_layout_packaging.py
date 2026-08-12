"""Isolated installed-layout smoke test for the shipped import graph.

``pyproject.toml`` ships the runtime packages plus ``data*``. A repository
checkout hides every ``scripts.*`` import
from those packages because the repo root sits on ``sys.path``; an installed
wheel does not contain ``scripts`` at all. Two whole features died there
silently:

* ``cli/gvf_run.py`` and ``cli/automated_workflow.py`` imported the HTML report
  builder from ``scripts.build_llm_trace_html`` inside a ``try``, so no
  installed run ever produced ``llm_trace_report.html``.
* ``cli/compare_variants.py`` had three UNGUARDED ``from
  scripts.ingest_review_adjudications import ...`` calls, so the scorer's
  adjudication-overlay path raised ``ImportError``.

This test builds a minimal installed layout (the shipped packages, no
``scripts/``) in a temp directory and drives the real entry points in a
subprocess whose ``sys.path`` cannot reach the repo root.
"""

from __future__ import annotations

import ast
import os
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SHIPPED_PACKAGES = (
    "pipeline",
    "harvesting",
    "config",
    "utils",
    "gene_literature",
    "cli",
    "data",
)


def _shipped_python_files() -> list[Path]:
    files: list[Path] = []
    for package in SHIPPED_PACKAGES:
        files.extend(sorted((REPO / package).rglob("*.py")))
    return files


def _scripts_imports(tree: ast.Module) -> list[tuple[int, str]]:
    """Every ``scripts.*`` import in a module, with its line number."""
    found: list[tuple[int, str]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] == "scripts":
                found.append((node.lineno, node.module or "scripts"))
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] == "scripts":
                    found.append((node.lineno, alias.name))
    return found


def _guarded_line_ranges(tree: ast.Module) -> list[range]:
    """Line ranges of ``try`` bodies that catch ImportError/Exception."""
    ranges: list[range] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Try):
            continue
        catches = False
        for handler in node.handlers:
            if handler.type is None:
                catches = True
            else:
                names = (
                    [handler.type]
                    if not isinstance(handler.type, ast.Tuple)
                    else list(handler.type.elts)
                )
                for name in names:
                    label = getattr(name, "id", None) or getattr(name, "attr", None)
                    if label in {"ImportError", "ModuleNotFoundError", "Exception"}:
                        catches = True
        if not catches:
            continue
        start = min(child.lineno for child in node.body)
        end = max(getattr(child, "end_lineno", child.lineno) for child in node.body)
        ranges.append(range(start, end + 1))
    return ranges


def test_no_shipped_module_has_an_unguarded_scripts_import():
    """A hard ``scripts.*`` import in a shipped package is broken in a wheel."""
    offenders: list[str] = []
    for path in _shipped_python_files():
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        guarded = _guarded_line_ranges(tree)
        for lineno, module in _scripts_imports(tree):
            if not any(lineno in span for span in guarded):
                offenders.append(
                    f"{path.relative_to(REPO)}:{lineno} imports {module} "
                    "outside any ImportError guard"
                )
    assert offenders == [], "unguarded scripts.* imports:\n" + "\n".join(offenders)


@pytest.fixture(scope="module")
def installed_layout(tmp_path_factory) -> Path:
    """A directory holding ONLY the shipped packages — no ``scripts``."""
    root = tmp_path_factory.mktemp("installed_layout")
    for package in SHIPPED_PACKAGES:
        shutil.copytree(
            REPO / package,
            root / package,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
    assert not (root / "scripts").exists()
    return root


def _run_in_layout(layout: Path, code: str) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    # PYTHONPATH is the layout only; cwd is the layout so sys.path[0] cannot
    # reach the repo root and re-expose `scripts`.
    env["PYTHONPATH"] = str(layout)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["GVF_LLM_TRACE_ENABLED"] = "false"
    env["GVF_DISABLE_LOCAL_DATA"] = "1"
    env.pop("GVF_LLM_TRACE_DIR", None)
    return subprocess.run(
        [sys.executable, "-B", "-c", textwrap.dedent(code)],
        cwd=layout,
        env=env,
        capture_output=True,
        text=True,
        timeout=180,
    )


def test_scripts_package_is_absent_from_the_installed_layout(installed_layout: Path):
    """Guard the guard: if `scripts` were importable the rest proves nothing."""
    result = _run_in_layout(
        installed_layout,
        """
        import importlib.util
        assert importlib.util.find_spec("scripts") is None, "scripts leaked in"
        print("OK")
        """,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_runtime_data_is_present_in_the_installed_layout(installed_layout: Path):
    """Aliases and reference gates must not silently disappear from wheels."""
    result = _run_in_layout(
        installed_layout,
        """
        from pipeline.reference_validation import load_reference_protein
        from utils.variant_normalizer import _load_gene_aliases

        aliases = _load_gene_aliases("KCNH2")
        assert len(aliases) > 4_000
        assert aliases["P.ALA1058GLU"] == "A1058E"
        sequence = load_reference_protein("KCNH2")
        assert sequence is not None and len(sequence) == 1_159
        print("OK", len(aliases), len(sequence))
        """,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_institutional_preflight_uses_shipped_session_helpers(installed_layout: Path):
    result = _run_in_layout(
        installed_layout,
        """
        from cli.institutional_preflight import probe_institutional_access
        from harvesting.paywall_session import (
            hydrate_session_with_browser_cookies, make_session,
        )

        assert callable(probe_institutional_access)
        assert callable(make_session)
        assert callable(hydrate_session_with_browser_cookies)
        print("OK")
        """,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_trace_html_builder_imports_and_runs_in_installed_layout(
    installed_layout: Path,
):
    result = _run_in_layout(
        installed_layout,
        """
        import json, tempfile
        from pathlib import Path
        from utils.llm_trace import (
            build_trace_manifest, configure_llm_tracing, llm_trace_scope,
            record_trace_event, reset_llm_tracing,
        )
        from utils.llm_trace_html import TRACE_REPORT_NAME, build_trace_html_report

        tmp = Path(tempfile.mkdtemp())
        root = configure_llm_tracing(tmp / "llm_traces", run_id="wheel", enabled=True)
        with llm_trace_scope(gene="KCNH2", pmid="12345678", stage="paper_curation"):
            record_trace_event("paper_curation_decision", {"variant_count": 1})
        build_trace_manifest(root, run_id="wheel")
        out = tmp / TRACE_REPORT_NAME
        data = build_trace_html_report(root, output_path=out, run_dir=tmp,
                                       run_id="wheel")
        reset_llm_tracing()
        assert out.is_file() and out.stat().st_size > 1000
        assert data["summary"]["decision_event_count"] == 1
        print("OK", data["integrity"]["level"])
        """,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_count_recovery_source_resolver_imports_in_installed_layout(
    installed_layout: Path,
):
    result = _run_in_layout(
        installed_layout,
        """
        import tempfile
        from pathlib import Path
        from pipeline.count_recovery import (
            corpus_root, default_source_roots, make_source_resolver,
        )
        from cli.gvf_run import step_count_recovery  # Step 3.55 import graph

        tmp = Path(tempfile.mkdtemp())
        (tmp / "pmc_fulltext").mkdir()
        (tmp / "pmc_fulltext" / "12345678_FULL_CONTEXT.md").write_text("body")
        roots = default_source_roots(tmp / "KCNH2.db")
        resolve = make_source_resolver("KCNH2", roots)
        assert resolve("12345678") == "body"
        assert resolve("99999999") is None
        assert corpus_root()
        assert callable(step_count_recovery)
        print("OK")
        """,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_adjudication_overlay_scoring_imports_in_installed_layout(
    installed_layout: Path,
):
    result = _run_in_layout(
        installed_layout,
        """
        from pathlib import Path
        from cli.compare_variants import (
            load_adjudication_overlay, load_adjudication_overlay_db,
            _adjudication_variant_key, _overlay_action,
        )
        from pipeline.adjudication_contract import VERDICT_TO_ACTION, variant_key

        # These three sites were unguarded `from scripts...` imports.
        assert load_adjudication_overlay(Path("/nonexistent.csv")) == {}
        assert load_adjudication_overlay_db(Path("/nonexistent.sqlite3"), "KCNH2") == {}
        assert _adjudication_variant_key("p.Leu552Ser")
        assert _overlay_action({"verdict": "correct_counts"}) == "count_override"
        assert variant_key("p.Leu552Ser") == _adjudication_variant_key("p.Leu552Ser")
        assert VERDICT_TO_ACTION["confirm"] == "gold_confirmed"
        print("OK")
        """,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


def test_gvf_run_help_works_in_installed_layout(installed_layout: Path):
    """The turnkey CLI must at least load and render help without `scripts`."""
    result = _run_in_layout(
        installed_layout,
        """
        import cli
        from typer.testing import CliRunner
        outcome = CliRunner().invoke(cli.app, ["gvf-run", "--help"])
        assert outcome.exit_code == 0, outcome.output
        assert "gvf-run" in outcome.output or "GENE" in outcome.output
        print("OK")
        """,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout
