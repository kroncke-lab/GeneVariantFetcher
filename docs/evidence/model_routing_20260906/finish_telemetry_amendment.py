"""Finish the explicit post-extraction, pre-score accounting amendment.

The first driver was stopped before lock when a second measurement bug was
found: production projection treated absent failed-call usage as zero. Preserve
both audit receipts, re-project with corrected usage, and prove scientific
prediction arrays are identical before invoking the ordinary native lock.
"""

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
RUN = HARNESS / "runs/20260906_model12_astra_medium_verified"
PYTHON = ROOT / ".venv/bin/python"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def scientific(data):
    return json.dumps(
        [
            {
                k: p.get(k)
                for k in [
                    "gene",
                    "pmid",
                    "variants",
                    "comparison_variants",
                    "external_linkage_variants",
                    "unattributed_variants",
                ]
            }
            for p in data["papers"]
        ],
        sort_keys=True,
    )


def main():
    receipt_path = HERE / "production_usage_amendment.json"
    assert not receipt_path.exists()
    assert not (RUN / "LOCK.json").exists()
    names = ["grok43", "astra_medium_verified", "grok43_astra_clinical_verified"]
    assert not any(
        (RUN.parent / ("20260906_model12_" + n) / "report.json").exists() for n in names
    )
    first = json.loads((HERE / "runtime_telemetry_amendment.json").read_text())
    assert sha((HERE / "budget.json").read_bytes()) == first["budget_sha256"]
    db_file = HARNESS / "db_to_predictions.py"
    db_before = subprocess.check_output(
        ["git", "show", "97b59b17:benchmarks/codex_paper_eval/db_to_predictions.py"],
        cwd=ROOT,
    )
    before = json.loads((HERE / "pre_telemetry_run_eval_snapshot.json").read_text())[
        "source"
    ].encode()
    value = hashlib.sha256()
    inventory = []
    for path in runtime_source_files():
        rel = path.relative_to(ROOT).as_posix()
        old = (
            before
            if path == HARNESS / "run_eval.py"
            else db_before
            if path == db_file
            else path.read_bytes()
        )
        value.update(rel.encode() + b"\0" + sha(old).encode() + b"\n")
        inventory.append(
            {
                "path": rel,
                "before_sha256": sha(old),
                "after_sha256": sha(path.read_bytes()),
            }
        )
    reconstructed = {"sha256": value.hexdigest(), "file_count": len(inventory)}
    assert reconstructed == first["prepared_runtime"]
    changed = [r["path"] for r in inventory if r["before_sha256"] != r["after_sha256"]]
    assert changed == [
        "benchmarks/codex_paper_eval/db_to_predictions.py",
        "benchmarks/codex_paper_eval/run_eval.py",
    ]
    snapshot = HERE / "pre_telemetry_projection_snapshot.json"
    snapshot.write_text(
        json.dumps(
            {
                "path": str(db_file.relative_to(ROOT)),
                "sha256": sha(db_before),
                "source": db_before.decode(),
            },
            indent=2,
        )
        + "\n"
    )
    prediction_file = RUN / "predictions.json"
    initial = json.loads(prediction_file.read_text())
    assert len(initial["papers"]) == 12
    science_before = scientific(initial)
    current_runtime = runtime_fingerprint()
    receipt = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "prepared_runtime": reconstructed,
        "locking_runtime": current_runtime,
        "changed_files": changed,
        "original_setup_preserved": True,
        "first_driver_disposition": "Stopped before lock; production export accounting omission found. Earlier audit remains preserved.",
        "new_scores_absent_on_all_final_arms": True,
        "budget_sha256": first["budget_sha256"],
        "driver_sha256": sha(Path(__file__).read_bytes()),
        "before_predictions_sha256": sha(prediction_file.read_bytes()),
        "scientific_predictions_sha256_before": sha(science_before.encode()),
        "projection_prior_snapshot_sha256": sha(snapshot.read_bytes()),
        "runtime_inventory": inventory,
    }
    # Save pre-score attestation before the replacement export.
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    command = [
        PYTHON,
        db_file,
        "--run-dir",
        RUN,
        "--production-root",
        RUN / "production_runs",
        "--trust-mode",
        "trusted",
        "--identity-mode",
        "trusted",
        "--paper-primary",
        "--out",
        prediction_file,
    ]
    with (HERE / "astra_final_projection.log").open("w") as log:
        subprocess.run(
            list(map(str, command)),
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    assert runtime_fingerprint() == current_runtime
    final = json.loads(prediction_file.read_text())
    assert scientific(final) == science_before, "Scientific content changed"
    receipt.update(
        after_predictions_sha256=sha(prediction_file.read_bytes()),
        scientific_predictions_unchanged=True,
    )
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    with (HERE / "astra_final_native_lock.log").open("w") as log:
        subprocess.run(
            list(map(str, [PYTHON, HARNESS / "run_eval.py", "lock", "--run-dir", RUN])),
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    print("Astra locked; scientific predictions unchanged", flush=True)


if __name__ == "__main__":
    main()
