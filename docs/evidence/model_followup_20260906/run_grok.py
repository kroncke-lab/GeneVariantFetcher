"""Grok's successful health gate permits the same frozen bounded roster task."""

import json
import time
from client import HERE, call, write, sha

prepared = json.loads((HERE / "prepared.json").read_text())
assert not (HERE / "grok_outputs_locked.json").exists()
write(
    HERE / "grok_execution_frozen.json",
    dict(
        time_unix=time.time(),
        arm="grok_low",
        contract_sha256=sha(HERE / "contract.md"),
        client_sha256=sha(HERE / "client.py"),
        run_sha256=sha(__import__("pathlib").Path(__file__)),
        prepared_sha256=sha(HERE / "prepared.json"),
        note="Same eight frozen packets/schema/caps/deadlines; low effort. Added after live health passes, before any Astra roster output inspection or score. Original Astra effort ablation unchanged.",
    ),
)
results = []
for name in prepared["packets"]:
    path = HERE / "packets" / (name + ".json")
    assert sha(path) == prepared["packets"][name]
    p = json.loads(path.read_text())
    tag = "grok_low_" + name
    try:
        call(
            tag,
            dict(
                model="grok-4.6",
                messages=[dict(role="user", content=p["prompt"])],
                reasoning_effort="low",
                max_completion_tokens=p["output_cap"],
                response_format={"type": "json_object"},
            ),
            timeout=p["timeout"],
        )
        results.append(
            dict(
                name=tag,
                path="responses/" + tag + ".json",
                sha256=sha(HERE / "responses" / (tag + ".json")),
            )
        )
    except Exception as e:
        results.append(dict(name=tag, status="not_dispatched", error=str(e)))
write(
    HERE / "grok_outputs_locked.json",
    dict(
        locked_unix=time.time(),
        cells=results,
        execution_sha256=sha(HERE / "grok_execution_frozen.json"),
    ),
)
