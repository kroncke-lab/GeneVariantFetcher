"""Post-hoc source-only output-contract probes; separate from frozen main arms."""

import json
import time
from client import HERE, call, write, sha

cols = [
    "sources",
    "variants",
    "family",
    "person",
    "genotype",
    "phenotype",
    "index_case",
    "endpoint_status",
]
plans = []
for model, name in [
    ("gpt-6-astra", "ryr304_multivariant"),
    ("grok-4.6", "myb213_patient_ids"),
]:
    p = json.loads((HERE / "packets" / (name + ".json")).read_text())
    props = dict(
        sources={
            "type": "array",
            "items": {"type": "string", "enum": list(p["source_units"])},
        },
        variants={"type": "array", "items": {"type": "string"}},
        family={"type": ["string", "null"]},
        person={"type": "string"},
        genotype={"type": "string", "enum": ["Y", "N", "?"]},
        phenotype={"type": "string"},
        index_case={"type": ["boolean", "null"]},
        endpoint_status={"type": "string", "enum": ["positive", "negative", "unknown"]},
    )
    schema = {
        "type": "object",
        "additionalProperties": False,
        "required": ["rows", "limitations"],
        "properties": {
            "rows": {
                "type": "array",
                "items": {
                    "type": "object",
                    "additionalProperties": False,
                    "required": cols,
                    "properties": props,
                },
            },
            "limitations": {"type": "array", "items": {"type": "string"}},
        },
    }
    before = p["prompt"].index("Return JSON only:")
    after = p["prompt"].index("\nTARGET:")
    prompt = (
        p["prompt"][:before]
        + "Return JSON matching the supplied schema. Each row is a record object. sources must use exact source IDs from the supplied enum, and must include both the person row and any cohort/footnote/join evidence needed to support phenotype or genotype.\n"
        + p["prompt"][after:]
    )
    body = dict(
        model=model,
        messages=[dict(role="user", content=prompt)],
        reasoning_effort="low",
        max_completion_tokens=4096,
        response_format={
            "type": "json_schema",
            "json_schema": {"name": "patient_roster", "strict": True, "schema": schema},
        },
    )
    plans.append(
        dict(
            name=("astra" if model == "gpt-6-astra" else "grok") + "_strict_" + name,
            packet=name,
            packet_sha256=sha(HERE / "packets" / (name + ".json")),
            body=body,
        )
    )
write(
    HERE / "strict_schema_plan.json",
    dict(
        prepared_unix=time.time(),
        classification="Post-hoc output schema/prompt refinement, after primary results and gold overlap were inspected. No gold in prompts. Not part of the registered effort comparison.",
        plans=plans,
    ),
)
results = []
for p in plans:
    call(p["name"], p["body"], timeout=180)
    path = HERE / "responses" / (p["name"] + ".json")
    results.append(
        dict(
            name=p["name"],
            packet=p["packet"],
            path=str(path.relative_to(HERE)),
            sha256=sha(path),
        )
    )
write(
    HERE / "strict_schema_locked.json",
    dict(
        locked_unix=time.time(),
        plan_sha256=sha(HERE / "strict_schema_plan.json"),
        cells=results,
    ),
)
