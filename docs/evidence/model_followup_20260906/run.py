"""Run frozen paired cells without reading the reference, then lock every output."""

import concurrent.futures
import json
import random
import time
from client import HERE, call, sha, write

prepared = json.loads((HERE / "prepared.json").read_text())
assert not (HERE / "outputs_locked.json").exists()
assert prepared["reference_sha256"] == sha(HERE / "source_reference.json")
write(
    HERE / "execution_frozen.json",
    dict(
        time_unix=time.time(),
        contract_sha256=sha(HERE / "contract.md"),
        client_sha256=sha(HERE / "client.py"),
        run_sha256=sha(__import__("pathlib").Path(__file__)),
        prepared_sha256=sha(HERE / "prepared.json"),
    ),
)
cells = [(name, effort) for name in prepared["packets"] for effort in ["low", "medium"]]
random.Random(prepared["order_seed"]).shuffle(cells)
write(HERE / "dispatch_order.json", cells)


def run(cell):
    name, effort = cell
    path = HERE / "packets" / (name + ".json")
    assert sha(path) == prepared["packets"][name]
    packet = json.loads(path.read_text())
    name = "astra_" + effort + "_" + name
    try:
        call(
            name,
            dict(
                model="gpt-6-astra",
                messages=[dict(role="user", content=packet["prompt"])],
                reasoning_effort=effort,
                max_completion_tokens=packet["output_cap"],
                response_format={"type": "json_object"},
            ),
            timeout=packet["timeout"],
        )
        return dict(
            name=name,
            path="responses/" + name + ".json",
            sha256=sha(HERE / "responses" / (name + ".json")),
        )
    except Exception as e:
        return dict(name=name, status="not_dispatched", error=str(e))


with concurrent.futures.ThreadPoolExecutor(2) as pool:
    results = list(pool.map(run, cells))
write(
    HERE / "outputs_locked.json",
    dict(
        locked_unix=time.time(),
        classification="All planned outputs locked before grading; nondispatches retained",
        cells=results,
        execution_sha256=sha(HERE / "execution_frozen.json"),
    ),
)
