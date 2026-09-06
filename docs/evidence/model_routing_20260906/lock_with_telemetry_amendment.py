"""Pre-score, experiment-only audit of a failed-usage measurement amendment.

Never alters the original setup fingerprint. Reconstructs its exact 251-file
fingerprint with only the benchmark run_eval.py at its pre-amendment bytes,
then uses the ordinary source rebind, trusted projection and native lock tools.
This is not permission for drift in any extraction or source-processing code.
"""

import ast
import difflib
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from benchmarks.codex_paper_eval.setup_production_eval import (
    runtime_fingerprint,
    runtime_source_files,
)

HERE = Path(__file__).resolve().parent
HARNESS = ROOT / "benchmarks/codex_paper_eval"
FOLDER = HARNESS / "runs/20260906_model12_astra_medium_verified"
PYTHON = ROOT / ".venv/bin/python"
BASE = "97b59b17"
TARGET = HARNESS / "run_eval.py"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def functions(source):
    return {
        n.name: ast.dump(n, include_attributes=False)
        for n in ast.parse(source).body
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def main():
    receipt = HERE / "runtime_telemetry_amendment.json"
    assert not receipt.exists(), "Preserve the original audit receipt"
    assert not (FOLDER / "LOCK.json").exists()
    final_names = ["grok43", "astra_medium_verified", "grok43_astra_clinical_verified"]
    assert not any(
        (FOLDER.parent / ("20260906_model12_" + n) / "report.json").exists()
        for n in final_names
    ), "Amendment must precede all new scores"
    budget = json.loads((HERE / "budget.json").read_text())
    assert not any(c["status"] == "reserved" for c in budget["calls"])
    current = TARGET.read_text()
    base = subprocess.check_output(
        ["git", "show", BASE + ":benchmarks/codex_paper_eval/run_eval.py"],
        cwd=ROOT,
        text=True,
    )
    # These two model-compatibility helpers were already changed before the
    # verified arm was prepared. All other amendments came later.
    old_tree, new_tree = ast.parse(base), ast.parse(current)
    old_lines, new_lines = (
        base.splitlines(keepends=True),
        current.splitlines(keepends=True),
    )
    for name in ["effective_effort", "reasoning_params"]:
        old = next(n for n in old_tree.body if getattr(n, "name", None) == name)
        new = next(n for n in new_tree.body if getattr(n, "name", None) == name)
        old_lines[old.lineno - 1 : old.end_lineno] = new_lines[
            new.lineno - 1 : new.end_lineno
        ]
    before = "".join(old_lines)
    value = hashlib.sha256()
    inventory = []
    for path in runtime_source_files():
        rel = path.relative_to(ROOT).as_posix()
        old_sha = sha(before.encode() if path == TARGET else path.read_bytes())
        new_sha = sha(path.read_bytes())
        value.update(rel.encode() + b"\0" + old_sha.encode() + b"\n")
        inventory.append(
            {"path": rel, "before_sha256": old_sha, "after_sha256": new_sha}
        )
    reconstructed = {"sha256": value.hexdigest(), "file_count": len(inventory)}
    setup = json.loads((FOLDER / "analysis_setup.json").read_text())
    assert reconstructed == setup["runtime"], "Cannot prove the exact prepared runtime"
    changed_files = [
        r["path"] for r in inventory if r["before_sha256"] != r["after_sha256"]
    ]
    assert changed_files == ["benchmarks/codex_paper_eval/run_eval.py"]
    old_fn, new_fn = functions(before), functions(current)
    added = sorted(new_fn.keys() - old_fn.keys())
    changed = sorted(k for k in old_fn if old_fn[k] != new_fn.get(k))
    assert added == ["documented_unknown_call_usage", "indexed_paper_call_usage"]
    assert changed == [
        "aggregate",
        "command_lock",
        "validate_predictions",
        "write_markdown_report",
    ]
    assert old_fn["score_one"] == new_fn["score_one"]
    assert sha((HERE / "budget_guard.py").read_bytes()) == setup["campaign_hook_sha256"]
    assert (
        sha((FOLDER / "frozen_corpus/source_snapshot.json").read_bytes())
        == setup["source_snapshot_sha256"]
    )
    statuses = sorted((FOLDER / "production_runs").glob("*/*/RUN_STATUS.json"))
    assert len(statuses) == 3
    for path in statuses:
        status = json.loads(path.read_text())
        assert status["status"] == "completed" and not status["stage_failures"]
        assert status["gold_access"]["disabled"]
    snapshot = HERE / "pre_telemetry_run_eval_snapshot.json"
    snapshot.write_text(
        json.dumps(
            {
                "path": str(TARGET.relative_to(ROOT)),
                "sha256": sha(before.encode()),
                "source": before,
            },
            indent=2,
        )
        + "\n"
    )
    diff = HERE / "runtime_telemetry_amendment.diff"
    diff.write_text(
        "".join(
            difflib.unified_diff(
                before.splitlines(keepends=True),
                current.splitlines(keepends=True),
                fromfile="prepared/run_eval.py",
                tofile="locked/run_eval.py",
            )
        )
    )
    payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "classification": "Pre-score measurement amendment on opened calibration; original setup preserved",
        "reason": "Authentic timeout usage is unknown. Retain failed papers and reserves rather than omit them or misstate zero cost.",
        "prepared_runtime": reconstructed,
        "locking_runtime": runtime_fingerprint(),
        "changed_files": changed_files,
        "new_functions": added,
        "changed_functions": changed,
        "unchanged_function_count": sum(old_fn[k] == new_fn.get(k) for k in old_fn),
        "scientific_scope": "No extraction, prompts, models, source processing, variant matching or count formulas changed. Aggregate changes are restricted to token accounting and incomplete-usage labels; score_one is AST-identical.",
        "validation": {
            "offline_unit_passed": 2938,
            "failed_usage_and_non_dispatch_focused_passed": 23,
        },
        "new_scores_absent_on_all_final_arms": True,
        "calls_in_flight": 0,
        "budget_sha256": sha((HERE / "budget.json").read_bytes()),
        "original_setup_sha256": sha((FOLDER / "analysis_setup.json").read_bytes()),
        "snapshot_sha256": sha(snapshot.read_bytes()),
        "diff_sha256": sha(diff.read_bytes()),
        "driver_sha256": sha(Path(__file__).read_bytes()),
        "runtime_inventory": inventory,
    }
    receipt.write_text(json.dumps(payload, indent=2) + "\n")
    commands = [
        [
            HARNESS / "rebind_production_sources.py",
            "--run-dir",
            FOLDER,
            "--production-root",
            FOLDER / "production_runs",
        ],
        [
            HARNESS / "db_to_predictions.py",
            "--run-dir",
            FOLDER,
            "--production-root",
            FOLDER / "production_runs",
            "--trust-mode",
            "trusted",
            "--identity-mode",
            "trusted",
            "--paper-primary",
            "--out",
            FOLDER / "predictions.json",
        ],
        [HARNESS / "run_eval.py", "lock", "--run-dir", FOLDER],
    ]
    for action, command in zip(["rebind", "projection", "lock"], commands):
        assert runtime_fingerprint() == payload["locking_runtime"], (
            "Further runtime drift"
        )
        with (HERE / f"astra_verified_amended_{action}.log").open("w") as log:
            subprocess.run(
                [str(PYTHON)] + list(map(str, command)),
                cwd=ROOT,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        print("COMPLETED", action, flush=True)


if __name__ == "__main__":
    main()
