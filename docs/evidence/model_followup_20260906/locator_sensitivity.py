"""Post-hoc format-only sensitivity; primary strict scores remain immutable."""

import copy
import json
import re
import time
from client import HERE, sha, write
from score import grade


def expand(value, units):
    m = re.fullmatch(r"(?:Table \d+, )?([LHTNA])(\d+)(?:-\1?(\d+))?", value)
    if not m:
        return [value]
    prefix, first, last = m.groups()
    a, b = int(first), int(last or first)
    if not 0 <= b - a <= 32:
        return [value]
    ids = [prefix + str(i) for i in range(a, b + 1)]
    return ids if all(i in units for i in ids) else [value]


scores = json.loads((HERE / "source_scores.json").read_text())
refs = json.loads((HERE / "source_reference.json").read_text())["rows"]
results = []
changes = []
for rr in scores["results"]:
    name = rr["arm"] + "_" + rr["packet"]
    p = HERE / "responses" / (name + ".json")
    d = json.loads(p.read_text())
    packet = json.loads((HERE / "packets" / (rr["packet"] + ".json")).read_text())
    if rr["status"] != "scored":
        continue
    parsed = json.loads(d["response"]["choices"][0]["message"]["content"])
    for row in parsed["rows"]:
        before = copy.deepcopy(row[0])
        row[0] = list(
            dict.fromkeys(i for v in before for i in expand(v, packet["source_units"]))
        )
        if before != row[0]:
            changes.append(
                dict(
                    arm=rr["arm"],
                    packet=rr["packet"],
                    person=row[3],
                    before=before,
                    after=row[0],
                )
            )
    d["response"]["choices"][0]["message"]["content"] = json.dumps(parsed)
    # In-memory diagnostic material lives outside original responses and locks.
    temp = HERE / "locator_sensitivity_material"
    temp.mkdir(exist_ok=True)
    dst = temp / (name + ".json")
    write(dst, d)
    out = grade(
        packet,
        refs[rr["packet"]],
        dict(path=str(dst.relative_to(HERE)), sha256=sha(dst)),
    )
    out["arm"] = rr["arm"]
    results.append(out)
summary = {
    a: {
        "source_bound_people": sum(r["exact_people"] for r in results if r["arm"] == a),
        "complete_positive_packets": sum(
            r.get("complete_source_bound_roster", False)
            for r in results
            if r["arm"] == a and r["reference_people"]
        ),
        "qualified_variants": sum(
            len(r["source_adjudicated_complete_variant_counts"])
            for r in results
            if r["arm"] == a
        ),
    }
    for a in ["astra_low", "astra_medium", "grok_low"]
}
write(
    HERE / "locator_sensitivity.json",
    dict(
        created_unix=time.time(),
        classification="Post-hoc normalization diagnostic, not the primary registered score or production acceptance",
        primary_scores_sha256=sha(HERE / "source_scores.json"),
        summary=summary,
        changes=changes,
        results=results,
    ),
)
print(json.dumps(summary, indent=2))
