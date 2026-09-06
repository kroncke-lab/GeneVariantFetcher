"""Resume unchanged remediation requests with two budget-isolated workers."""

import concurrent.futures
import json
import subprocess
import sys
import threading
import time
from pathlib import Path

from client import HERE, call, sha, write


def main():
    plan_path = HERE / "remediation_plan.json"
    plan = json.loads(plan_path.read_text())
    assert not (HERE / "remediation_locked.json").exists()
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
    receipt = HERE / "remediation_boundary_stop.json"
    assert receipt.exists()
    amendment = HERE / "remediation_concurrency_amendment.json"
    assert not amendment.exists()
    write(
        amendment,
        dict(
            recorded_unix=time.time(),
            plan_sha256=sha(plan_path),
            stop_receipt_sha256=sha(receipt),
            workers=2,
            minimum_start_interval_seconds=22,
            runner_sha256=sha(Path(__file__)),
            reason="Completed serial calls take minutes. Resume remaining identical requests with two independently locked worker processes and globally spaced subprocess launches. Azure key retrieval precedes HTTP dispatch, so launch spacing is not a guaranteed exact dispatch interval. No clinical outputs inspected, no in-flight requests canceled, no prompt/cap/deadline/reference changes.",
        ),
    )
    pending = [
        p
        for p in plan["plans"]
        if not (HERE / "responses" / (p["name"] + ".json")).exists()
    ]
    timing_lock = threading.Lock()
    last_start = [0.0]

    def worker(p):
        with timing_lock:
            pause = max(0, 22 - (time.monotonic() - last_start[0]))
            if pause:
                time.sleep(pause)
            last_start[0] = time.monotonic()
        return p["name"], subprocess.run(
            [sys.executable, __file__, p["name"]], check=False
        ).returncode

    deferred = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for name, code in pool.map(worker, pending):
            if code == 3:
                deferred.append(name)
            elif code:
                raise RuntimeError(f"Worker failed: {name}, {code}")
    for name in deferred:
        _, code = worker(next(p for p in pending if p["name"] == name))
        if code not in (0, 3):
            raise RuntimeError(f"Deferred worker failed: {name}, {code}")
    cells = []
    for p in plan["plans"]:
        target = HERE / "responses" / (p["name"] + ".json")
        cells.append(
            dict(
                name=p["name"],
                dispatched=target.exists(),
                **(
                    {"sha256": sha(target)}
                    if target.exists()
                    else {"reason": "Budget exhausted before dispatch"}
                ),
            )
        )
    write(
        HERE / "remediation_locked.json",
        dict(
            locked_unix=time.time(),
            plan_sha256=sha(plan_path),
            concurrency_amendment_sha256=sha(amendment),
            cells=cells,
        ),
    )


if __name__ == "__main__":
    main()
