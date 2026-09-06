"""Post-hoc cheaper reader package; never replaces the locked primary test."""

import copy
import json
import time

from client import HERE, call, sha, write


def main():
    original = json.loads((HERE / "plan.json").read_text())
    assert (HERE / "source_scores.json").exists()
    assert (HERE / "outputs_locked.json").exists()
    assert not (HERE / "remediation_plan.json").exists()
    plans = []
    for p in original["plans"]:
        if p["arm"] != "grok_first":
            continue
        body = copy.deepcopy(p["body"])
        schema = body["response_format"]["json_schema"]["schema"]
        body["messages"][0]["content"] += (
            "\n\nOUTPUT CONTRACT (validated locally after your response):\n"
            + json.dumps(schema, ensure_ascii=False)
        )
        body["response_format"] = {"type": "json_object"}
        body["max_completion_tokens"] = 8192
        plans.append(
            dict(
                name=p["packet"] + "__grok_remediated",
                packet=p["packet"],
                arm="grok_remediated",
                body=body,
                schema=schema,
            )
        )
    plan = dict(
        prepared_unix=time.time(),
        original_plan_sha256=sha(HERE / "plan.json"),
        primary_lock_sha256=sha(HERE / "outputs_locked.json"),
        primary_score_sha256=sha(HERE / "source_scores.json"),
        classification="Post-hoc production-compatible cheaper-reader package, selected after primary scoring. Full eight-paper panel; same source, queries, references and local schema. No reference values or earlier model answers enter requests. Changes inline schema/prompt adapter, provider JSON-object format, cap and deadline together. This is a package comparison, not a causal ablation of any one setting, and not unchanged production transport.",
        rationale="Current recovery uses an inline JSON contract, 8192 output cap and the ordinary 1200-second timeout. Test whether a usable cheaper reader eliminates apparent Astra gains caused by the bounded strict-schema configuration. json_object adds syntactic JSON enforcement; semantic schema remains local.",
        timeout_seconds=1200,
        minimum_start_interval_seconds=22,
        retries=0,
        budget="Shares the existing 2.17 campaign ceiling; reserves every call before dispatch, retains unknown charges, and reports any nondispatch. No new allowance.",
        plans=plans,
    )
    write(HERE / "remediation_plan.json", plan)
    cells = []
    last_start = 0
    for p in plans:
        pause = max(
            0, plan["minimum_start_interval_seconds"] - (time.monotonic() - last_start)
        )
        if pause:
            time.sleep(pause)
        last_start = time.monotonic()
        try:
            call(p["name"], p["body"], timeout=plan["timeout_seconds"])
        except RuntimeError as exc:
            if "budget exhausted" not in str(exc):
                raise
            cells.append(dict(name=p["name"], dispatched=False, reason=str(exc)))
            continue
        target = HERE / "responses" / (p["name"] + ".json")
        cells.append(dict(name=p["name"], dispatched=True, sha256=sha(target)))
    write(
        HERE / "remediation_locked.json",
        dict(
            locked_unix=time.time(),
            plan_sha256=sha(HERE / "remediation_plan.json"),
            cells=cells,
        ),
    )


if __name__ == "__main__":
    main()
