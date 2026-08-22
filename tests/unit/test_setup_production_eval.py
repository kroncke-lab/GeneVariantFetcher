"""Regression coverage for the registry-pinned production gold evaluation setup."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks.codex_paper_eval import setup_production_eval as setup


def test_direct_cli_entrypoint_resolves_repository_imports(tmp_path: Path):
    result = subprocess.run(
        [sys.executable, str(Path(setup.__file__).resolve()), "--help"],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "create" in result.stdout
    assert "check" in result.stdout


def test_setup_preserves_virtual_environment_interpreter_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    base = tmp_path / "base-python"
    base.write_text("")
    environment_python = tmp_path / "venv" / "bin" / "python"
    environment_python.parent.mkdir(parents=True)
    environment_python.symlink_to(base)
    monkeypatch.setattr(setup.sys, "executable", str(environment_python))

    assert setup.runtime_python() == environment_python.absolute()
    assert setup.runtime_python() != environment_python.resolve()


def test_live_gold120_contract_is_the_corrected_118_attempt_cohort():
    contract = setup.cohort_contract(setup.DEFAULT_MANIFEST, setup.DEFAULT_REGISTRY)

    assert contract["id"] == "gold_120"
    assert contract["attempt_count"] == 118
    assert contract["unique_pmid_count"] == 114
    assert contract["gene_attempt_counts"] == {
        "KCNH2": 28,
        "KCNQ1": 30,
        "RYR2": 30,
        "SCN5A": 30,
    }
    assert ("KCNH2", "10086972") not in contract["rows"]
    assert ("KCNH2", "14642689") not in contract["rows"]


def test_generated_extraction_is_calibrated_blinded_and_nonpublishing(tmp_path: Path):
    script = setup.make_extraction_script(
        run_dir=tmp_path / "run",
        genes=["KCNH2", "KCNQ1", "RYR2", "SCN5A"],
        email="test@example.org",
        python=tmp_path / "python",
    )

    assert script.count("python -m cli gvf-run") == 0  # interpreter is a shell variable
    assert script.count('"$PYTHON" -m cli gvf-run') == 4
    assert script.count("--pmid-file") == 4
    assert script.count("--no-source-recovery") == 4
    assert script.count("--no-corpus-sync") == 4
    assert script.count("--no-publish-review") == 4
    assert script.count("--gold-free-run") == 4
    assert script.count(") &") == 4
    assert script.count("pids+=($!)") == 4
    assert 'for pid in "${pids[@]}"' in script
    assert "At least one gene extraction failed" in script
    assert script.count("operator_logs/") == 4
    assert "--full-coverage" not in script
    assert "--gold-root" not in script
    assert "run_eval.py score" not in script


def test_finalize_projects_and_locks_before_opening_gold(tmp_path: Path):
    script = setup.make_finalize_script(
        run_dir=tmp_path / "run", python=tmp_path / "python"
    )

    assert "benchmarks/codex_paper_eval/db_to_predictions.py" in script
    assert "runs/20260726_fixed48_production/db_to_predictions.py" not in script
    assert script.index("rebind_production_sources.py") < script.index(
        "db_to_predictions.py"
    )
    assert script.index("db_to_predictions.py") < script.index(" lock --run-dir")
    assert script.index(" lock --run-dir") < script.index(" score --run-dir")
    assert "--trust-mode trusted --identity-mode trusted" in script


def test_runbook_labels_trusted_projection_as_primary(tmp_path: Path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    contract = {
        "manifest": str(tmp_path / "tier.tsv"),
        "sha256": "fixture",
        "attempt_count": 1,
        "unique_pmid_count": 1,
        "gene_attempt_counts": {"KCNH2": 1},
    }

    setup.write_runbook(run_dir, contract, ["KCNH2"])

    runbook = " ".join((run_dir / "RUNBOOK.md").read_text().split())
    assert "trusted count and identity projection" in runbook
    assert "must not replace the locked primary" in runbook
    assert "raw (`--trust-mode all`) projection is the blinded primary" not in runbook


def test_check_run_detects_pmid_drift_without_reading_gold(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    contract = setup.cohort_contract(setup.DEFAULT_MANIFEST, setup.DEFAULT_REGISTRY)
    run_dir = tmp_path / "run"
    pmid_dir = run_dir / "pmids"
    pmid_dir.mkdir(parents=True)
    papers = [{"gene": gene, "pmid": pmid} for gene, pmid in contract["rows"]]
    (run_dir / "selection.json").write_text(json.dumps({"papers": papers}))
    fingerprint = {"sha256": "fixture", "file_count": 1}
    (run_dir / "setup.json").write_text(
        json.dumps(
            {
                "cohort": {
                    "manifest": str(setup.DEFAULT_MANIFEST),
                    "registry": str(setup.DEFAULT_REGISTRY),
                },
                "repository": {"runtime_source": fingerprint},
            }
        )
    )
    by_gene: dict[str, list[str]] = {}
    for gene, pmid in contract["rows"]:
        by_gene.setdefault(gene, []).append(pmid)
    for gene, pmids in by_gene.items():
        (pmid_dir / f"{gene}.txt").write_text("\n".join(pmids) + "\n")
    extraction = setup.make_extraction_script(
        run_dir=run_dir,
        genes=sorted(by_gene),
        email="test@example.org",
        python=tmp_path / "python",
    )
    (run_dir / "run_extraction.sh").write_text(extraction)
    monkeypatch.setattr(setup, "runtime_fingerprint", lambda: fingerprint)

    setup.check_run(run_dir)

    kcnq1 = pmid_dir / "KCNQ1.txt"
    kcnq1.write_text(kcnq1.read_text().replace("\n", "\n99999999\n", 1))
    with pytest.raises(setup.SetupError, match="KCNQ1 PMID file drifted"):
        setup.check_run(run_dir)
