"""Offline integrity and grade replay; does not dispatch APIs or rewrite evidence."""

import ast
import hashlib
import json
from pathlib import Path
from client import HERE, sha
from score import grade


def normalized_ast(text):
    """Ruff splits `import a, b`; preserve import order while comparing syntax."""

    class SplitImports(ast.NodeTransformer):
        def visit_Import(self, node):
            return [ast.Import(names=[alias]) for alias in node.names]

    return ast.dump(SplitImports().visit(ast.parse(text)))


def main():
    snapshots = json.loads((HERE / "executed_code_snapshots.json").read_text())["files"]

    def bound(name, digest):
        if sha(HERE / name) == digest:
            return
        saved = snapshots[name]
        assert saved["sha256"] == digest
        assert hashlib.sha256(saved["text"].encode()).hexdigest() == digest
        assert normalized_ast(saved["text"]) == normalized_ast(
            (HERE / name).read_text()
        ), name

    prepared = json.loads((HERE / "prepared.json").read_text())
    assert sha(HERE / "source_reference.json") == prepared["reference_sha256"]
    bound("prepare.py", prepared["source_preparation_sha256"])
    for name, digest in prepared["packets"].items():
        assert sha(HERE / "packets" / (name + ".json")) == digest
    for fn in ["execution_frozen.json", "grok_execution_frozen.json"]:
        e = json.loads((HERE / fn).read_text())
        assert sha(HERE / "contract.md") == e["contract_sha256"]
        assert sha(HERE / "prepared.json") == e["prepared_sha256"]
        bound("client.py", e["client_sha256"])
        bound(
            "run.py" if fn == "execution_frozen.json" else "run_grok.py",
            e["run_sha256"],
        )
    refs = json.loads((HERE / "source_reference.json").read_text())["rows"]
    scores = json.loads((HERE / "source_scores.json").read_text())
    bound("score.py", scores["scorer_sha256"])
    count = 0
    for lockname in ["outputs_locked.json", "grok_outputs_locked.json"]:
        assert sha(HERE / lockname) == scores["lock_hashes"][lockname]
        lock = json.loads((HERE / lockname).read_text())
        assert lock["locked_unix"] < scores["graded_unix"]
        execution = (
            "execution_frozen.json"
            if lockname == "outputs_locked.json"
            else "grok_execution_frozen.json"
        )
        assert sha(HERE / execution) == lock["execution_sha256"]
        for c in lock["cells"]:
            prefix = (
                "grok_low_"
                if c["name"].startswith("grok_")
                else (
                    "astra_low_"
                    if c["name"].startswith("astra_low_")
                    else "astra_medium_"
                )
            )
            name = c["name"][len(prefix) :]
            p = json.loads((HERE / "packets" / (name + ".json")).read_text())
            r = grade(p, refs[name], c)
            r["arm"] = prefix[:-1]
            expected = next(
                x
                for x in scores["results"]
                if x["arm"] == r["arm"] and x["packet"] == name
            )
            assert r == expected, c["name"]
            count += 1
    d = json.loads((HERE / "budget.json").read_text())
    summary = json.loads((HERE / "budget_summary.json").read_text())
    assert sha(HERE / "budget.json") == summary["budget_sha256"]
    assert (
        sha(HERE.parent / "model_routing_20260906/budget.json")
        == d["prior_budget_sha256"]
    )
    assert sum(c["accounted_usd"] for c in d["calls"]) <= d["limit_usd"]
    assert summary["combined_accounted_usd"] <= 150
    for c in d["calls"]:
        assert c["status"] != "reserved"
        if "response_sha256" in c:
            assert (
                sha(HERE / "responses" / (c["name"] + ".json")) == c["response_sha256"]
            )
    assert sum(len(refs[k]) for k in refs if k != "myb204_structured") == 93
    print(
        json.dumps(
            dict(
                primary_replayed=count,
                immutable_source_packets=len(prepared["packets"]),
                unique_source_people=93,
                old_and_new_budget_integrity="passed",
                code_snapshot_AST_equivalence="passed",
            )
        )
    )


if __name__ == "__main__":
    main()
