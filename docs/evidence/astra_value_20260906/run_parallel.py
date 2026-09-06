"""Resume frozen cells concurrently in isolated retry-free worker processes."""

import concurrent.futures
import json
import subprocess
import sys
import time
from pathlib import Path

from client import HERE, call, sha, write


def main():
    plan_path = HERE / "plan.json"
    plan = json.loads(plan_path.read_text())
    assert not (HERE / "outputs_locked.json").exists()
    assert plan["reference_sha256"] == sha(HERE / "reference_values.json")
    for name, digest in plan["packet_sha256"].items():
        assert digest == sha(HERE / "packets" / (name + ".json"))
    if len(sys.argv) == 2:
        p = next(p for p in plan["plans"] if p["name"] == sys.argv[1])
        try:
            call(p["name"], p["body"], timeout=plan["timeout_seconds"])
        except RuntimeError as exc:
            if "budget exhausted" not in str(exc):
                raise
            raise SystemExit(3) from exc
        return
    budget = json.loads((HERE / "budget.json").read_text())
    assert all(c["status"] != "reserved" for c in budget["calls"])
    amendment = HERE / "concurrency_amendment.json"
    assert not amendment.exists()
    write(
        amendment,
        dict(
            recorded_unix=time.time(),
            plan_sha256=sha(plan_path),
            completed_calls=len(budget["calls"]),
            workers=3,
            reason="Serial Grok latency observed from transport metadata. Original runner stopped between completed calls with no pending reservation; no clinical answers inspected. Remaining unchanged requests run in separate processes, preserving per-call budget locks. This changes scheduling only, not prompt, model, cap, deadline, reference, or scoring. No API request canceled or retried.",
            runner_sha256=sha(Path(__file__)),
        ),
    )
    pending = [
        p
        for p in plan["plans"]
        if not (HERE / "responses" / (p["name"] + ".json")).exists()
    ]

    def worker(p):
        return p["name"], subprocess.run(
            [sys.executable, __file__, p["name"]], check=False
        ).returncode

    deferred = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for name, code in pool.map(worker, pending):
            if code == 3:
                deferred.append(name)
            elif code:
                raise RuntimeError(
                    f"Worker failed without a completed response: {name}, {code}"
                )
    # A reservation can temporarily occupy allowance needed by another worker.
    # Re-evaluate only undispatched cells after all initial reservations settle.
    for name in deferred:
        code = subprocess.run([sys.executable, __file__, name], check=False).returncode
        if code not in (0, 3):
            raise RuntimeError(f"Deferred worker failed: {name}, {code}")
    cells = []
    for p in plan["plans"]:
        target = HERE / "responses" / (p["name"] + ".json")
        if target.exists():
            cells.append(dict(name=p["name"], dispatched=True, sha256=sha(target)))
        else:
            cells.append(
                dict(
                    name=p["name"],
                    dispatched=False,
                    reason="API budget exhausted before dispatch",
                )
            )
    write(
        HERE / "outputs_locked.json",
        dict(
            locked_unix=time.time(),
            plan_sha256=sha(plan_path),
            concurrency_amendment_sha256=sha(amendment),
            cells=cells,
        ),
    )


if __name__ == "__main__":
    main()
