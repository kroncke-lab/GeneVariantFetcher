"""Offline transport and budget audit; response access requires both output locks.

Default mode verifies existing evidence without writing. --write-summary creates
the final budget summary only after both frozen panels and all ledger rows are
settled. This module uses only the standard library and cannot dispatch APIs.
"""

from __future__ import annotations

import argparse
import ast
import fcntl
import hashlib
import json
import math
import os
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
RATES = {"gpt-6-astra": (10, 50), "grok-4.6": (2, 6)}
EXPECTED_PRIOR_HASH = "dbeabed9d309cb03264606ea0b521188a872d20f47d23ab4f97bf02a8c02221e"
EXPECTED_ORIGINAL_HASH = (
    "03f21ea91625de7933ba71b20cf9cbf9e2b60a5c06d5092826c71ecb2cbd2fc8"
)
PANELS = (
    ("primary", "plan.json", "outputs_locked.json", "concurrency_amendment.json"),
    (
        "remediation",
        "remediation_plan.json",
        "remediation_locked.json",
        "remediation_concurrency_amendment.json",
    ),
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(actual, expected, label):
    require(isinstance(actual, (int, float)) and math.isfinite(actual), label)
    require(math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-10), label)


def normalized_ast(source):
    """Allow formatting and Ruff's split imports, but no executable changes."""

    class SplitImports(ast.NodeTransformer):
        def visit_Import(self, node):
            return [ast.Import(names=[alias]) for alias in node.names]

    return ast.dump(SplitImports().visit(ast.parse(source)))


def snapshot_text(saved):
    """Read either historical snapshot spelling and always validate its hash."""
    values = [saved[key] for key in ("text", "content") if key in saved]
    require(
        bool(values) and all(isinstance(value, str) for value in values),
        "Snapshot has no source text",
    )
    require(
        all(value == values[0] for value in values), "Snapshot text/content disagree"
    )
    require(
        hashlib.sha256(values[0].encode()).hexdigest() == saved["sha256"],
        "Snapshot source digest mismatch",
    )
    return values[0]


def without_main_ast(source):
    """Permit deletion of exactly one top-level `if __name__ == '__main__'`."""
    tree = ast.parse(source)
    expected = ast.dump(ast.parse("__name__ == '__main__'", mode="eval").body)
    matches = [
        node
        for node in tree.body
        if isinstance(node, ast.If) and ast.dump(node.test) == expected
    ]
    require(
        len(matches) == 1 and not matches[0].orelse,
        "Expected one main-only guard without else",
    )
    removed = matches[0]
    tree.body.remove(removed)
    return normalized_ast(ast.unparse(tree)), hashlib.sha256(
        ast.dump(removed).encode()
    ).hexdigest()


def verify_client_cleanup(saved_source, digest, current_source):
    """Check the single recorded import-only cleanup; no other AST edits pass."""
    proof = read(HERE / "client_import_only_cleanup.json")
    require(
        proof["file"] == "client.py" and proof["original_sha256"] == digest,
        "Client cleanup identity mismatch",
    )
    require(
        hashlib.sha256(proof["after_source"].encode()).hexdigest()
        == proof["after_sha256"],
        "Client cleanup source digest",
    )
    expected, removed_digest = without_main_ast(saved_source)
    require(
        proof["removed_main_ast_sha256"] == removed_digest,
        "Client cleanup removed-node digest",
    )
    require(
        expected == normalized_ast(proof["after_source"]),
        "Client cleanup changed code outside __main__",
    )
    require(
        expected == normalized_ast(current_source),
        "Current client differs from recorded import-only cleanup",
    )
    locks = [read(HERE / lock_name)["locked_unix"] for _, _, lock_name, _ in PANELS]
    require(
        proof["recorded_unix"] >= max(locks),
        "Client cleanup predates final output locks",
    )


def verify_code(name, digest, snapshots):
    path = HERE / name
    matches = [
        files[name]
        for files in snapshots
        if name in files and files[name]["sha256"] == digest
    ]
    # Validate all supplied source receipts even if the current file hash is
    # exact; the shortcut must not conceal a damaged or differently keyed receipt.
    sources = [snapshot_text(saved) for saved in matches]
    if sha(path) == digest:
        return "exact_hash"
    require(bool(matches), f"No executed source snapshot binds {name}")
    saved_source = sources[0]
    current_source = path.read_text()
    if name == "client.py" and (HERE / "client_import_only_cleanup.json").exists():
        verify_client_cleanup(saved_source, digest, current_source)
        return "recorded_import_only_cleanup"
    require(
        normalized_ast(saved_source) == normalized_ast(current_source),
        f"Executable code differs from recorded snapshot: {name}",
    )
    return "recorded_source_AST_equivalent"


def usage_cost(body, result):
    """Reproduce the retry-free client's conservative charge calculation."""
    model = body["model"]
    require(model in RATES, f"Unknown model: {model}")
    input_rate, output_rate = RATES[model]
    cap = next(
        body[k]
        for k in ("max_completion_tokens", "max_output_tokens", "max_tokens")
        if k in body
    )
    require(type(cap) is int and cap > 0, "Invalid output cap")
    input_bound = len(json.dumps(body).encode()) + 2048
    require(input_bound < 272000, "Input no longer fits the recorded price tier")
    reserve = (input_bound * input_rate * 1.25 + cap * output_rate) / 1e6
    usage = result.get("response", {}).get("usage")
    if not usage:
        require(not result.get("usage"), "Top-level usage disagrees with response")
        return {
            "cap": cap,
            "reserve": reserve,
            "known": False,
            "accounted": reserve,
            "proxy": None,
            "input_tokens": None,
            "output_tokens": None,
            "input_bound": input_bound,
        }
    require(result.get("usage") == usage, "Top-level and response usage differ")
    ni = usage.get("prompt_tokens", usage.get("input_tokens", 0))
    no = max(
        usage.get("completion_tokens", usage.get("output_tokens", 0)),
        usage.get("total_tokens", 0) - ni,
    )
    require(type(ni) is int and ni >= 0, "Invalid input token count")
    require(type(no) is int and no >= 0, "Invalid output token count")
    require(ni <= input_bound, "Actual input exceeded conservative reservation bound")
    proxy = (ni * input_rate + no * output_rate) / 1e6
    accounted = (ni * input_rate * 1.25 + no * output_rate) / 1e6
    return {
        "cap": cap,
        "reserve": reserve,
        "known": True,
        "accounted": accounted,
        "proxy": proxy,
        "input_tokens": ni,
        "output_tokens": no,
        "input_bound": input_bound,
    }


def budget_chain(ledger):
    prior_path = HERE.parent / "model_followup_20260906/budget.json"
    original_path = HERE.parent / "model_routing_20260906/budget.json"
    require(sha(prior_path) == EXPECTED_PRIOR_HASH, "Prior ledger changed")
    require(sha(original_path) == EXPECTED_ORIGINAL_HASH, "Original ledger changed")
    require(ledger["prior_budget_sha256"] == EXPECTED_PRIOR_HASH, "Prior hash mismatch")
    require(
        ledger["original_budget_sha256"] == EXPECTED_ORIGINAL_HASH,
        "Original hash mismatch",
    )
    prior, original = read(prior_path), read(original_path)
    require(
        prior["prior_budget_sha256"] == EXPECTED_ORIGINAL_HASH, "Broken ledger chain"
    )
    original_accounted = original["smoke_uncertainty_reserve_usd"] + sum(
        c.get("accounted_usd", c["reserved_usd"]) for c in original["calls"]
    )
    close(prior["prior_accounted_usd"], original_accounted, "Original total mismatch")
    prior_accounted = original_accounted + sum(
        c["accounted_usd"] for c in prior["calls"]
    )
    close(
        ledger["prior_accounted_usd"], prior_accounted, "Prior accounted total mismatch"
    )
    close(original["api_ceiling_usd"], 150, "Original envelope changed")
    close(ledger["limit_usd"], 2.17, "Current campaign envelope changed")
    require(
        prior_accounted + ledger["limit_usd"] <= 150 + 1e-10, "Overallocated envelope"
    )
    return prior_accounted


def inspect_locked_panels():
    """Load no response until all required locks are present and valid."""
    for _, _, lock_name, _ in PANELS:
        require(
            (HERE / lock_name).exists(), f"Required output lock is absent: {lock_name}"
        )
    primary = read(HERE / "plan.json")
    require(
        sha(HERE / "reference_values.json") == primary["reference_sha256"],
        "Reference changed after freeze",
    )
    for name, digest in primary["packet_sha256"].items():
        require(
            sha(HERE / "packets" / (name + ".json")) == digest,
            f"Packet changed: {name}",
        )
    require(len(primary["packet_sha256"]) == 8, "Unexpected primary paper count")
    panels = []
    for label, plan_name, lock_name, amendment_name in PANELS:
        plan, lock = read(HERE / plan_name), read(HERE / lock_name)
        require(lock["plan_sha256"] == sha(HERE / plan_name), f"{label} plan hash")
        require(plan["prepared_unix"] <= lock["locked_unix"], f"{label} chronology")
        amendment = read(HERE / amendment_name)
        require(
            lock["concurrency_amendment_sha256"] == sha(HERE / amendment_name),
            f"{label} scheduling hash",
        )
        require(
            amendment["plan_sha256"] == lock["plan_sha256"], f"{label} amendment plan"
        )
        require(
            plan["prepared_unix"] <= amendment["recorded_unix"] <= lock["locked_unix"],
            f"{label} amendment chronology",
        )
        names = [p["name"] for p in plan["plans"]]
        locked_names = [c["name"] for c in lock["cells"]]
        require(len(set(names)) == len(names), f"{label} duplicate planned calls")
        require(
            len(set(locked_names)) == len(locked_names),
            f"{label} duplicate locked calls",
        )
        require(set(names) == set(locked_names), f"{label} planned/locked calls differ")
        panels.append((label, plan, lock, amendment))
    remediation = panels[1][1]
    require(
        remediation["original_plan_sha256"] == sha(HERE / "plan.json"),
        "Original plan changed",
    )
    require(
        remediation["primary_lock_sha256"] == sha(HERE / "outputs_locked.json"),
        "Primary lock changed",
    )
    require(
        remediation["primary_score_sha256"] == sha(HERE / "source_scores.json"),
        "Primary scores changed",
    )
    require(
        panels[0][2]["locked_unix"] <= remediation["prepared_unix"],
        "Remediation chronology",
    )
    require(
        len(primary["plans"]) == 24 and len(remediation["plans"]) == 8,
        "Unexpected panel sizes",
    )
    stop = HERE / "remediation_boundary_stop.json"
    require(
        panels[1][3]["stop_receipt_sha256"] == sha(stop),
        "Remediation boundary receipt changed",
    )
    return panels


def build_summary():
    panels = inspect_locked_panels()
    ledger = read(HERE / "budget.json")
    names = [c["name"] for c in ledger["calls"]]
    require(len(set(names)) == len(names), "Duplicate budget call")
    require(
        all(c["status"] != "reserved" for c in ledger["calls"]),
        "Unsettled budget reservation",
    )
    prior_accounted = budget_chain(ledger)
    rows = {c["name"]: c for c in ledger["calls"]}
    verified, reports = set(), []
    for label, plan, lock, _ in panels:
        planned = {p["name"]: p for p in plan["plans"]}
        for cell in lock["cells"]:
            name = cell["name"]
            require(type(cell["dispatched"]) is bool, f"Invalid dispatch flag: {name}")
            target = HERE / "responses" / (name + ".json")
            if not cell["dispatched"]:
                require(
                    name not in rows and not target.exists(),
                    f"Undispatched call has artifacts: {name}",
                )
                continue
            require(name in rows, f"Response missing budget row: {name}")
            row = rows[name]
            require(
                sha(target) == cell["sha256"] == row["response_sha256"],
                f"Response hash: {name}",
            )
            # Both panels have been locked and the ledger is settled before this read.
            result = read(target)
            body = planned[name]["body"]
            require(result["name"] == name, f"Response name: {name}")
            require(result["request"] == body, f"Request changed: {name}")
            require(
                result["timeout_seconds"] == plan["timeout_seconds"],
                f"Deadline: {name}",
            )
            require(result["retries"] == plan["retries"] == 0, f"Retry setting: {name}")
            require(result["status"] == row["status"], f"Status mismatch: {name}")
            require(row["model"] == body["model"], f"Model mismatch: {name}")
            require(
                plan["prepared_unix"] <= row["started_unix"] <= lock["locked_unix"],
                f"Call chronology: {name}",
            )
            calc = usage_cost(body, result)
            require(row["cap"] == calc["cap"], f"Cap mismatch: {name}")
            require(row["usage_known"] is calc["known"], f"Usage flag mismatch: {name}")
            for actual in (row["reserved_usd"], result["reserved_usd"]):
                close(actual, calc["reserve"], f"Reservation math: {name}")
            for actual in (row["accounted_usd"], result["accounted_usd"]):
                close(actual, calc["accounted"], f"Accounting math: {name}")
            close(row["seconds"], result["seconds"], f"Elapsed time mismatch: {name}")
            if calc["known"]:
                for actual in (row["api_proxy_usd"], result["api_proxy_usd"]):
                    close(actual, calc["proxy"], f"Proxy math: {name}")
            else:
                require(
                    "api_proxy_usd" not in row and "api_proxy_usd" not in result,
                    f"Unknown charge relabeled as proxy: {name}",
                )
            reports.append(
                {"name": name, "panel": label, "model": row["model"], **calc}
            )
            verified.add(name)
    require(verified == set(rows), "Budget contains an unplanned or unlocked API call")
    actual_files = {p.stem for p in (HERE / "responses").glob("*.json")}
    require(
        actual_files == verified,
        "Response directory contains an unlocked/unaccounted call",
    )
    known = [r for r in reports if r["known"]]
    unknown = [r for r in reports if not r["known"]]
    accounted = sum(r["accounted"] for r in reports)
    close(
        accounted,
        sum(c["accounted_usd"] for c in ledger["calls"]),
        "Ledger sum mismatch",
    )
    require(accounted <= ledger["limit_usd"] + 1e-10, "Campaign exceeded $2.17")
    combined = prior_accounted + accounted
    require(combined <= 150 + 1e-10, "Campaigns exceeded $150")
    snapshots = [read(HERE / "executed_code_snapshots.json")["files"]]
    for filename in (
        "remediation_code_snapshots.json",
        "audit_transport_snapshot.json",
    ):
        extra = HERE / filename
        if extra.exists():
            snapshots.append(read(extra)["files"])
    code = {}
    for files in snapshots:
        for name, saved in files.items():
            snapshot_text(saved)
            verify_code(name, saved["sha256"], snapshots)
            code[name] = {"recorded_sha256": saved["sha256"], "verified": True}
    for label, _, _, amendment in panels:
        name = (
            "run_parallel.py" if label == "primary" else "run_remediation_parallel.py"
        )
        verify_code(name, amendment["runner_sha256"], snapshots)
        code[name] = {"recorded_sha256": amendment["runner_sha256"], "verified": True}
    primary_integrity = read(HERE / "source_scores.json")["integrity"]
    require(
        sha(HERE / "score.py") == primary_integrity["score_code_sha256"],
        "Primary scorer changed after scoring",
    )
    require(
        sha(HERE.parents[2] / "pipeline/count_recovery.py")
        == primary_integrity["literal_validator_code_sha256"],
        "Literal validator changed after scoring",
    )
    code["score.py"] = {
        "recorded_sha256": primary_integrity["score_code_sha256"],
        "verified": True,
    }
    code["pipeline/count_recovery.py"] = {
        "recorded_sha256": primary_integrity["literal_validator_code_sha256"],
        "verified": True,
    }
    if (HERE / "client_import_only_cleanup.json").exists():
        code["client.py"]["allowed_exception"] = (
            "Only the unused __main__ health-probe block was removed after both output locks; imported transport AST is unchanged."
        )
        code["client.py"]["exception_receipt_sha256"] = sha(
            HERE / "client_import_only_cleanup.json"
        )
    return {
        "actual_api_requests": len(reports),
        "returned_usage_requests": len(known),
        "unknown_usage_requests": len(unknown),
        "returned_api_proxy_usd": sum(r["proxy"] for r in known),
        "unknown_reservations_usd": sum(r["reserve"] for r in unknown),
        "returned_accounting_with_input_premium_usd": sum(
            r["accounted"] for r in known
        ),
        "new_accounted_usd": accounted,
        "prior_accounted_usd": prior_accounted,
        "new_limit_usd": ledger["limit_usd"],
        "combined_accounted_usd": combined,
        "remaining_original_150_envelope_usd": 150 - combined,
        "remaining_current_campaign_usd": ledger["limit_usd"] - accounted,
        "budget_sha256": sha(HERE / "budget.json"),
        "prior_budget_sha256": EXPECTED_PRIOR_HASH,
        "original_budget_sha256": EXPECTED_ORIGINAL_HASH,
        "lock_hashes": {name: sha(HERE / name) for _, _, name, _ in PANELS},
        "panels": {
            label: {
                "planned": len(plan["plans"]),
                "dispatched": sum(c["dispatched"] for c in lock["cells"]),
                "undispatched": sum(not c["dispatched"] for c in lock["cells"]),
            }
            for label, plan, lock, _ in panels
        },
        "code_integrity": code,
        "verified_response_hashes": len(verified),
        "price_basis": "Uncached proxy rates plus 25% input premium in conservative accounting; reasoning uses max(completion, total-input). Unknown charges retain full reservations. Azure invoice not reconciled; CLI costs excluded.",
    }


def record_client_cleanup():
    """Record an already-made main-guard deletion; never edit the client itself."""
    inspect_locked_panels()
    ledger = read(HERE / "budget.json")
    require(
        all(c["status"] != "reserved" for c in ledger["calls"]),
        "Client cleanup requires settled calls",
    )
    target = HERE / "client_import_only_cleanup.json"
    require(not target.exists(), "Client cleanup proof already exists")
    saved = read(HERE / "executed_code_snapshots.json")["files"]["client.py"]
    before = snapshot_text(saved)
    after = (HERE / "client.py").read_text()
    expected, removed_digest = without_main_ast(before)
    require(
        expected == normalized_ast(after),
        "Client changed beyond deleting its main-only block",
    )
    proof = {
        "recorded_unix": time.time(),
        "file": "client.py",
        "original_sha256": saved["sha256"],
        "after_sha256": sha(HERE / "client.py"),
        "after_source": after,
        "removed_main_ast_sha256": removed_digest,
        "note": "After both output panels locked and all API reservations settled, remove only the unused historical direct-run health probes. Frozen runners import this module; all imported transport code is AST-identical. Original executed source and snapshot timing notes remain unchanged. No paid request is dispatched under this cleanup.",
    }
    target.write_text(json.dumps(proof, indent=2) + "\n")
    verify_client_cleanup(before, saved["sha256"], after)
    return {"cleanup_receipt": str(target), "sha256": sha(target)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-summary", action="store_true")
    parser.add_argument("--record-client-cleanup", action="store_true")
    args = parser.parse_args()
    # A shared lock prevents a concurrent client from altering ledger accounting
    # while the audit checks it and (optionally) atomically writes its summary.
    with (HERE / "budget.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_SH)
        if args.record_client_cleanup:
            require(
                not args.write_summary,
                "Record cleanup separately before the final summary",
            )
            print(json.dumps(record_client_cleanup(), indent=2))
            return
        summary = build_summary()
        path = HERE / "budget_summary.json"
        if args.write_summary:
            if path.exists():
                require(
                    read(path) == summary,
                    "Existing summary differs; refusing to overwrite final evidence",
                )
            else:
                temporary = path.with_suffix(".json.tmp")
                temporary.write_text(json.dumps(summary, indent=2) + "\n")
                os.replace(temporary, path)
        elif path.exists():
            require(read(path) == summary, "Saved summary fails offline replay")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
