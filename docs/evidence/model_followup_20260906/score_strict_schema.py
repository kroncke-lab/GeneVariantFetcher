"""Post-hoc strict-schema probe grade; adapt objects to the unchanged row comparator."""

import json
import time
from client import HERE, write, sha
from score import COLS, grade

lock = json.loads((HERE / "strict_schema_locked.json").read_text())
refs = json.loads((HERE / "source_reference.json").read_text())["rows"]
out = []
for c in lock["cells"]:
    path = HERE / c["path"]
    assert sha(path) == c["sha256"]
    d = json.loads(path.read_text())
    p = json.loads((HERE / "packets" / (c["packet"] + ".json")).read_text())
    if d["status"] == "returned":
        parsed = json.loads(d["response"]["choices"][0]["message"]["content"])
        converted = dict(
            columns=COLS,
            rows=[[r[k] for k in COLS] for r in parsed["rows"]],
            limitations=parsed["limitations"],
        )
        d["response"]["choices"][0]["message"]["content"] = json.dumps(converted)
    material = HERE / "strict_schema_material"
    material.mkdir(exist_ok=True)
    dest = material / (c["name"] + ".json")
    write(dest, d)
    r = grade(
        p, refs[c["packet"]], dict(path=str(dest.relative_to(HERE)), sha256=sha(dest))
    )
    r["model_probe"] = c["name"]
    out.append(r)
write(
    HERE / "strict_schema_scores.json",
    dict(
        graded_unix=time.time(),
        classification="Post-hoc schema and prompt refinement, not primary randomized effort comparison",
        lock_sha256=sha(HERE / "strict_schema_locked.json"),
        results=out,
    ),
)
for r in out:
    print(
        json.dumps(
            {
                k: r.get(k)
                for k in [
                    "model_probe",
                    "status",
                    "exact_people",
                    "reference_people",
                    "field_differences",
                    "source_binding_errors",
                    "complete_source_bound_roster",
                ]
            }
        )
    )
