"""Execute frozen initial cells, then lock every returned output before grading."""

import json
import time

from client import HERE, call, sha, write


def main():
    plan_path = HERE / "plan.json"
    plan = json.loads(plan_path.read_text())
    assert not (HERE / "outputs_locked.json").exists()
    assert plan["reference_sha256"] == sha(HERE / "reference_values.json")
    for name, digest in plan["packet_sha256"].items():
        assert digest == sha(HERE / "packets" / (name + ".json"))
    cells = []
    for p in plan["plans"]:
        target = HERE / "responses" / (p["name"] + ".json")
        if not target.exists():
            try:
                call(p["name"], p["body"], timeout=plan["timeout_seconds"])
            except RuntimeError as exc:
                if "budget exhausted" not in str(exc):
                    raise
                cells.append(dict(name=p["name"], dispatched=False, reason=str(exc)))
                continue
        cells.append(dict(name=p["name"], dispatched=True, sha256=sha(target)))
    write(
        HERE / "outputs_locked.json",
        dict(locked_unix=time.time(), plan_sha256=sha(plan_path), cells=cells),
    )


if __name__ == "__main__":
    main()
