"""Score the separately locked cheaper-reader package without changing primary evidence.

The provider accepts JSON-object output; the original strict schema is still
validated locally. Package changes include prompt adapter, format, cap, deadline
and scheduling. This is not a causal format ablation or unchanged production run.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "astra_value_primary_score", HERE / "score.py"
)
primary = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(primary)
SENSITIVE_QUERY = ("ryr193_abstract", "q1_affected")


def key(row):
    return row["packet"], row["id"]


def verify_concurrency_amendment(plan, lock):
    if "concurrency_amendment_sha256" not in lock:
        return {}
    path = HERE / "remediation_concurrency_amendment.json"
    assert primary.sha(path) == lock["concurrency_amendment_sha256"]
    amendment = primary.read(path)
    assert amendment["plan_sha256"] == lock["plan_sha256"]
    assert plan["prepared_unix"] <= amendment["recorded_unix"] <= lock["locked_unix"]
    receipt_path = HERE / "remediation_boundary_stop.json"
    assert primary.sha(receipt_path) == amendment["stop_receipt_sha256"]
    return {
        "concurrency_amendment_sha256": lock["concurrency_amendment_sha256"],
        "boundary_stop_receipt_sha256": amendment["stop_receipt_sha256"],
    }


def verify_remediation():
    """Require the new lock before opening any remediation response."""
    lock_path = HERE / "remediation_locked.json"
    if not lock_path.exists():
        raise RuntimeError("Remediation outputs must be locked before scoring")
    lock = primary.read(lock_path)
    plan_path = HERE / "remediation_plan.json"
    assert primary.sha(plan_path) == lock["plan_sha256"]
    plan = primary.read(plan_path)
    verify_concurrency_amendment(plan, lock)
    assert plan["original_plan_sha256"] == primary.sha(HERE / "plan.json")
    assert plan["primary_lock_sha256"] == primary.sha(HERE / "outputs_locked.json")
    assert plan["primary_score_sha256"] == primary.sha(HERE / "source_scores.json")
    original, original_lock, packets, references, _ = primary.verify_inputs()
    saved_primary = primary.read(HERE / "source_scores.json")
    assert primary.build_scores() == saved_primary, (
        "Primary scores fail exact offline replay"
    )
    assert original_lock["locked_unix"] <= plan["prepared_unix"] <= lock["locked_unix"]
    originals = {p["packet"]: p for p in original["plans"] if p["arm"] == "grok_first"}
    names = [p["name"] for p in plan["plans"]]
    assert len(set(names)) == len(names)
    locked_names = [c["name"] for c in lock["cells"]]
    assert len(set(locked_names)) == len(locked_names)
    assert set(names) == set(locked_names)
    assert [p["packet"] for p in plan["plans"]] == list(originals)
    locked = {c["name"]: c for c in lock["cells"]}
    responses = {}
    for p in plan["plans"]:
        assert p["arm"] == "grok_remediated"
        assert p["name"] == p["packet"] + "__grok_remediated"
        original_body = originals[p["packet"]]["body"]
        schema = original_body["response_format"]["json_schema"]["schema"]
        assert p["schema"] == schema
        # Reconstruct the entire permitted adapter from answer-free frozen inputs.
        expected_body = copy.deepcopy(original_body)
        expected_body["messages"][0]["content"] += (
            "\n\nOUTPUT CONTRACT (validated locally after your response):\n"
            + json.dumps(schema, ensure_ascii=False)
        )
        expected_body["response_format"] = {"type": "json_object"}
        expected_body["max_completion_tokens"] = 8192
        assert p["body"] == expected_body, (
            "Remediation body changed more than its frozen adapter"
        )
        cell = locked[p["name"]]
        assert type(cell["dispatched"]) is bool
        path = HERE / "responses" / (p["name"] + ".json")
        if not cell["dispatched"]:
            assert not path.exists(), "Undispatched remediation cell has a response"
            responses[p["name"]] = None
            continue
        assert primary.sha(path) == cell["sha256"]
        response = primary.read(path)
        assert response["name"] == p["name"]
        assert response["request"] == p["body"]
        assert response["timeout_seconds"] == plan["timeout_seconds"]
        assert response["retries"] == plan["retries"] == 0
        responses[p["name"]] = response
    return plan, lock, packets, references, saved_primary, responses


def local_contract_cell(plan_cell):
    """Use the frozen schema for local grading, not as a claimed wire format."""
    return {
        **plan_cell,
        "body": {"response_format": {"json_schema": {"schema": plan_cell["schema"]}}},
    }


def subset_comparison(name, packets, remediation_cells, astra_cells):
    """Compare the same planned paper set; keep failed/nondispatched cells."""
    wanted = set(packets)
    grok = [c for c in remediation_cells if c["packet"] in wanted]
    astra = [c for c in astra_cells if c["packet"] in wanted]
    assert {c["packet"] for c in grok} == {c["packet"] for c in astra} == wanted
    assert len(grok) == len(astra) == len(wanted)
    grok_cost, astra_cost = primary.cost_summary(grok), primary.cost_summary(astra)
    both_known = all(
        c["dispatched"] and c["api_proxy_usd"] is not None for c in [*grok, *astra]
    )
    return {
        "name": name,
        "packets": sorted(wanted),
        "grok_remediated": {
            "metrics": primary.metrics([r for c in grok for r in c["queries"]]),
            "operations": grok_cost,
        },
        "astra_primary_180s_4096_strict": {
            "metrics": primary.metrics([r for c in astra for r in c["queries"]]),
            "operations": astra_cost,
        },
        "all_planned_paired_requests_dispatched_with_known_usage": both_known,
        "grok_to_astra_returned_proxy_ratio": (
            grok_cost["returned_api_proxy_usd"] / astra_cost["returned_api_proxy_usd"]
            if both_known and astra_cost["returned_api_proxy_usd"] > 0
            else None
        ),
        "latency_comparability": (
            "Observed package latency only: cap, timeout, prompt, format and scheduling differ. "
            "Primary 180-second failures are not evidence of intrinsic reasoning weakness."
        ),
    }


def scoped_metrics(all_arm_rows):
    full = {arm: primary.metrics(rows) for arm, rows in all_arm_rows.items()}
    sensitivity = {
        arm: primary.metrics([row for row in rows if key(row) != SENSITIVE_QUERY])
        for arm, rows in all_arm_rows.items()
    }
    return full, sensitivity


def build_scores():
    plan, lock, packets, references, saved_primary, responses = verify_remediation()
    cells = [
        primary.grade_cell(
            local_contract_cell(p),
            packets[p["packet"]],
            references[p["packet"]],
            responses[p["name"]],
        )
        for p in plan["plans"]
    ]
    rows = [r for c in cells for r in c["queries"]]
    by_key = {key(r): r for r in rows}
    assert len(by_key) == len(rows)
    original_cells = saved_primary["cells"]
    astra_cells = [c for c in original_cells if c["arm"] == "astra_low"]
    astra_rows = [r for c in astra_cells for r in c["queries"]]
    astra_index = {key(r): r for r in astra_rows}
    previous_unique = saved_primary["astra_unique_numeric_exact_beyond_both_groks"]
    recovery = []
    for previous in previous_unique:
        row = by_key[key(previous)]
        recovery.append(
            {
                "packet": row["packet"],
                "pmid": row["pmid"],
                "id": row["id"],
                "variant": row["variant"],
                "field": row["field"],
                "expected": row["expected"],
                "previous_grok_first_outcome": previous["grok_first_outcome"],
                "previous_grok_repeat_outcome": previous["grok_repeat_outcome"],
                "astra_primary_claim": previous["astra_claim"],
                "grok_remediated_claim": row["claim"],
                "grok_remediated_outcome": row["numeric_outcome"],
                "grok_recovered_numeric_exact": row["numeric_outcome"]
                == "exact_nonnull",
                "grok_mechanical_source_quote_ok": row["mechanical_source_quote_ok"],
                "grok_current_validator_accepted_correct": row["accepted_correct"],
                "astra_current_validator_accepted_correct": astra_index[key(row)][
                    "accepted_correct"
                ],
                "included_in_42_query_sensitivity": key(row) != SENSITIVE_QUERY,
            }
        )
    all_arm_rows = {
        a: [r for c in original_cells if c["arm"] == a for r in c["queries"]]
        for a in primary.ARMS
    }
    all_arm_rows["grok_remediated"] = rows
    full, sensitivity = scoped_metrics(all_arm_rows)
    overlays = {}
    for lane in ("raw_candidates", "current_validator"):
        out = primary.compare_overlay(rows, astra_rows, lane)
        cost = primary.cost_summary(astra_cells)
        gain = out["new_correct_nonnull"]
        net = out["net_additional_correct_minus_wrong"]
        out["cost_assumption"] = (
            "Replay of Astra on every paper; this is not a measured conditional escalation policy."
        )
        out["incremental_astra_returned_proxy_usd"] = cost["returned_api_proxy_usd"]
        out["incremental_astra_accounted_usd"] = cost[
            "accounted_usd_including_unknown_reserves"
        ]
        out["proxy_usd_per_added_correct_nonnull"] = (
            cost["returned_api_proxy_usd"] / gain if gain else None
        )
        out["proxy_usd_per_net_added_correct_minus_wrong"] = (
            cost["returned_api_proxy_usd"] / net if net > 0 else None
        )
        overlays[lane] = out
    unique_packets = {r["packet"] for r in previous_unique}
    matched_dispatched = {c["packet"] for c in cells if c["dispatched"]} & {
        c["packet"] for c in astra_cells if c["dispatched"]
    }
    matched_returned = {
        c["packet"]
        for c in cells
        if c["request_successful_stop"] and c["api_proxy_usd"] is not None
    } & {
        c["packet"]
        for c in astra_cells
        if c["request_successful_stop"] and c["api_proxy_usd"] is not None
    }
    paired_costs = [
        subset_comparison("all_planned_papers", packets, cells, astra_cells),
        subset_comparison(
            "papers_with_previous_Astra_only_numeric_values",
            unique_packets,
            cells,
            astra_cells,
        ),
        subset_comparison(
            "both_packages_dispatched_secondary_subset",
            matched_dispatched,
            cells,
            astra_cells,
        ),
        subset_comparison(
            "both_packages_successful_with_known_usage_secondary_subset",
            matched_returned,
            cells,
            astra_cells,
        ),
    ]
    remaining = [r for r in recovery if not r["grok_recovered_numeric_exact"]]
    return {
        "classification": plan["classification"],
        "integrity": {
            "remediation_plan_sha256": primary.sha(HERE / "remediation_plan.json"),
            "remediation_lock_sha256": primary.sha(HERE / "remediation_locked.json"),
            "primary_plan_sha256": plan["original_plan_sha256"],
            "primary_lock_sha256": plan["primary_lock_sha256"],
            "primary_score_sha256": plan["primary_score_sha256"],
            "reference_sha256": primary.sha(HERE / "reference_values.json"),
            "score_code_sha256": primary.sha(Path(__file__)),
            "reused_primary_score_code_sha256": primary.sha(HERE / "score.py"),
            "prepared_unix": plan["prepared_unix"],
            "locked_unix": lock["locked_unix"],
            **verify_concurrency_amendment(plan, lock),
        },
        "limitations": [
            "This post-hoc package preserves frozen queries, sources, references and local validation; it changes prompt adapter, provider format, cap, deadline and scheduling together.",
            "It is neither a causal strict-schema ablation nor an execution of unchanged production transport.",
            "The regular runtime's 1200-second default differs from the primary test's 180-second deadline; primary timeout advantages cannot alone establish regular-protocol value or intrinsic intelligence.",
            "Numeric agreement is with investigator packet references. Mechanical quote and source-ID checks do not establish clinical semantic support or whole-paper scope.",
            "Current-validator acceptance is replayed without source-reference correctness as a gate. Derived claims remain unaccepted candidates.",
            "Every planned query, including nondispatched and failed cells, remains in the full-panel denominator. Failure is never correct abstention.",
            "The sensitivity omits one disputed RYR2 19398417 symptom-complement inference without altering frozen references or original primary scores.",
            "Matched-paper cost compares corresponding packets. API proxies are not Azure invoices; unknown reserves remain separately accounted.",
            "The hypothetical Astra overlay preserves existing nonnull values and is not a tested trigger-selected production escalation policy.",
        ],
        "full_43_query_metrics": full,
        "sensitivity_42_queries": {
            "excluded_query": {"packet": SENSITIVE_QUERY[0], "id": SENSITIVE_QUERY[1]},
            "reason": "The abstract's 'including two without symptoms' wording may not independently prove that every other carrier was symptomatic. The frozen primary reference remains unchanged.",
            "metrics": sensitivity,
        },
        "grok_remediated_operations": primary.cost_summary(cells),
        "previous_astra_only_values": {
            "number": len(recovery),
            "papers": sorted(unique_packets),
            "recovered_numeric_exact_by_grok": sum(
                r["grok_recovered_numeric_exact"] for r in recovery
            ),
            "recovered_numeric_exact_and_mechanical_pass": sum(
                r["grok_recovered_numeric_exact"]
                and r["grok_mechanical_source_quote_ok"]
                for r in recovery
            ),
            "remaining_numeric_advantage_values": len(remaining),
            "remaining_numeric_advantage_papers": sorted(
                {r["packet"] for r in remaining}
            ),
            "queries": recovery,
        },
        "astra_missing_only_overlay_on_remediated_grok": overlays,
        "matched_paper_quality_and_cost": paired_costs,
        "cells": cells,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify offline replay without rewriting evidence",
    )
    args = parser.parse_args()
    result = build_scores()
    target = HERE / "remediated_scores.json"
    if args.check or target.exists():
        assert primary.read(target) == result, (
            "Existing remediation scores differ; do not overwrite evidence"
        )
    else:
        target.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(
        json.dumps(
            {
                "grok_remediated": result["full_43_query_metrics"]["grok_remediated"],
                "operations": result["grok_remediated_operations"],
                "previous_unique_recovered": result["previous_astra_only_values"][
                    "recovered_numeric_exact_by_grok"
                ],
                "previous_unique_total": result["previous_astra_only_values"]["number"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
