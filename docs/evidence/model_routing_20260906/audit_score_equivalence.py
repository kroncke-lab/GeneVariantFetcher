"""Verify the accounting amendment leaves aggregate scientific scores unchanged."""

import ast
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from benchmarks.codex_paper_eval import run_eval

HERE = Path(__file__).resolve().parent


def main():
    snapshot = json.loads((HERE / "pre_telemetry_run_eval_snapshot.json").read_text())
    assert hashlib.sha256(snapshot["source"].encode()).hexdigest() == snapshot["sha256"]
    old_tree = ast.parse(snapshot["source"])
    node = next(n for n in old_tree.body if getattr(n, "name", None) == "aggregate")
    namespace = vars(run_eval).copy()
    exec(
        compile(
            ast.Module(body=[node], type_ignores=[]), "<prepared aggregate>", "exec"
        ),
        namespace,
    )
    names = ["grok43", "astra_medium_verified", "grok43_astra_clinical_verified"]
    audit = {}
    for name in names:
        path = (
            ROOT
            / "benchmarks/codex_paper_eval/runs"
            / ("20260906_model12_" + name)
            / "report.json"
        )
        report = json.loads(path.read_text())
        old, new = (
            namespace["aggregate"](report["papers"]),
            run_eval.aggregate(report["papers"]),
        )
        old.pop("token_usage")
        new.pop("token_usage")
        assert old == new
        audit[name] = {
            "all_non_usage_aggregate_fields_identical": True,
            "report_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    (HERE / "score_equivalence_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    print("All scientific aggregates identical under prepared and amended harnesses")


if __name__ == "__main__":
    main()
