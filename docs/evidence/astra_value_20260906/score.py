"""Offline paired scoring; no API calls, database writes or gold-standard access.

The investigator's packet reference measures numeric agreement, never acceptance.
Mechanical source checks prove locator membership and verbatim quotation only;
they do not prove clinical meaning, completeness, genotype ownership or timepoint.
The existing literal validator is replayed unchanged. Derived claims are retained
as candidates and never passed off as production-accepted literal counts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

from jsonschema import Draft202012Validator

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))

from pipeline.count_recovery import (  # noqa: E402
    PER_VARIANT_ROLES,
    PaperGap,
    VariantGap,
    validate_paper_response,
)

ARMS = ("grok_first", "grok_repeat", "astra_low")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def normalize(text):
    """Whitespace only: case, punctuation and scientific spelling stay exact."""
    return re.sub(r"\s+", " ", text).strip()


def verify_inputs():
    """Read response files only after a complete, hash-bound lock exists."""
    lock_path = HERE / "outputs_locked.json"
    if not lock_path.exists():
        raise RuntimeError("All initial outputs must be locked before scoring")
    lock = read(lock_path)
    plan_path = HERE / "plan.json"
    assert sha(plan_path) == lock["plan_sha256"], "Frozen plan changed"
    plan = read(plan_path)
    assert lock["locked_unix"] >= plan["prepared_unix"]
    if "concurrency_amendment_sha256" in lock:
        amendment_path = HERE / "concurrency_amendment.json"
        assert sha(amendment_path) == lock["concurrency_amendment_sha256"]
        amendment = read(amendment_path)
        assert amendment["plan_sha256"] == lock["plan_sha256"]
        assert (
            plan["prepared_unix"] <= amendment["recorded_unix"] <= lock["locked_unix"]
        )
    ref_path = HERE / "reference_values.json"
    assert sha(ref_path) == plan["reference_sha256"], "Source reference changed"
    references = read(ref_path)
    packets = {}
    for name, digest in plan["packet_sha256"].items():
        path = HERE / "packets" / (name + ".json")
        assert sha(path) == digest, f"Frozen packet changed: {name}"
        packet = read(path)
        assert packet["name"] == name
        assert sha(ROOT / packet["source_path"]) == packet["source_sha256"]
        qids = [q["id"] for q in packet["queries"]]
        assert len(set(qids)) == len(qids)
        assert set(qids) == set(references[name])
        assert all(
            v is None or (type(v) is int and v >= 0) for v in references[name].values()
        )
        packets[name] = packet
    assert set(packets) == set(references)
    planned_names = [p["name"] for p in plan["plans"]]
    locked_names = [p["name"] for p in lock["cells"]]
    assert len(set(planned_names)) == len(planned_names)
    assert len(set(locked_names)) == len(locked_names)
    assert set(planned_names) == set(locked_names), "Incomplete output lock"
    locked = {c["name"]: c for c in lock["cells"]}
    responses = {}
    for p in plan["plans"]:
        cell = locked[p["name"]]
        assert p["arm"] in ARMS
        assert p["body"]["messages"] == [
            {"role": "user", "content": packets[p["packet"]]["prompt"]}
        ]
        assert type(cell["dispatched"]) is bool
        target = HERE / "responses" / (p["name"] + ".json")
        if not cell["dispatched"]:
            assert not target.exists(), "Undispatched cell has a response"
            responses[p["name"]] = None
            continue
        assert sha(target) == cell["sha256"], f"Locked response changed: {target}"
        response = read(target)
        assert response["name"] == p["name"]
        assert response["request"] == p["body"], "Actual request differs from plan"
        responses[p["name"]] = response
    for name in packets:
        assert Counter(
            p["arm"] for p in plan["plans"] if p["packet"] == name
        ) == Counter(ARMS)
    return plan, lock, packets, references, responses


def validator_replay(packet, query, claim):
    """Translate an explicit claim, without consulting the reference value."""
    if claim.get("basis") != "explicit" or type(claim.get("value")) is not int:
        return {
            "attempted": False,
            "accepted": False,
            "reason": "Only model-declared explicit nonnull integers are replayed",
        }
    gap = PaperGap(
        packet["gene"],
        packet["pmid"],
        [VariantGap(1, query["variant"], [query["field"]], paper_derived=True)],
    )
    entry = {
        "variant": query["variant"],
        query["field"]: claim["value"],
        "quote": claim.get("quote"),
        "count_role": PER_VARIANT_ROLES[query["field"]],
        "evidence_locator": {"source_ids": ",".join(claim.get("sources", []))},
    }
    source_text = "\n".join(packet["source_units"].values())
    result = validate_paper_response(gap, [entry], source_text)
    return {
        "attempted": True,
        "accepted": bool(result.accepted),
        "rejections": result.rejected,
        "accepted_evidence": [
            {
                "field": v.field,
                "value": v.value,
                "quote": v.quote,
                "count_role": v.count_role,
                "evidence_locator": v.evidence_locator,
            }
            for v in result.accepted
        ],
    }


def grade_cell(plan_cell, packet, reference, response):
    dispatched = response is not None
    errors = []
    payload = None
    finish = None
    schema_valid = False
    if dispatched:
        choices = response.get("response", {}).get("choices", [])
        if len(choices) == 1:
            finish = choices[0].get("finish_reason")
            content = choices[0].get("message", {}).get("content")
            if isinstance(content, str):
                try:
                    payload = json.loads(content)
                except (ValueError, TypeError) as exc:
                    errors.append("invalid_json: " + str(exc))
            else:
                errors.append("missing_string_content")
        else:
            errors.append("expected_one_choice")
    else:
        errors.append("not_dispatched")
    if payload is not None:
        schema = plan_cell["body"]["response_format"]["json_schema"]["schema"]
        schema_errors = list(Draft202012Validator(schema).iter_errors(payload))
        errors.extend(
            "schema: " + str(e.json_path) + ": " + e.message for e in schema_errors
        )
        schema_valid = not schema_errors
    request_ok = (
        dispatched
        and response.get("status") == "returned"
        and response.get("http_status") == 200
        and finish == "stop"
    )
    if dispatched and not request_ok:
        errors.append("response_not_successful_stop")
    rows = payload.get("claims", []) if isinstance(payload, dict) else []
    if not isinstance(rows, list):
        rows = []
    counts = Counter(
        c.get("id")
        for c in rows
        if isinstance(c, dict) and isinstance(c.get("id"), str)
    )
    qids = {q["id"] for q in packet["queries"]}
    duplicate_ids = sorted(k for k, n in counts.items() if n > 1)
    missing_ids = sorted(qids - set(counts))
    extra_ids = sorted(set(counts) - qids)
    allowed_sources = set(packet["source_units"])
    source_text = normalize("\n".join(packet["source_units"].values()))
    grades = []
    for query in packet["queries"]:
        qid = query["id"]
        candidates = [c for c in rows if isinstance(c, dict) and c.get("id") == qid]
        claim = candidates[0] if len(candidates) == 1 else None
        usable = request_ok and schema_valid and claim is not None
        value = claim.get("value") if usable else None
        if usable and value is not None and (type(value) is not int or value < 0):
            usable = False
            value = None
        expected = reference[qid]
        numeric_exact = usable and value == expected
        if not usable:
            outcome = "response_failure_or_missing_claim"
        elif value is None:
            outcome = "correct_abstention" if expected is None else "null_miss"
        elif numeric_exact:
            outcome = "exact_nonnull"
        else:
            outcome = "wrong_nonnull"
        mechanical_errors = []
        replay = {"attempted": False, "accepted": False}
        if usable:
            sources = claim["sources"]
            if any(s not in allowed_sources for s in sources):
                mechanical_errors.append("invalid_source_id")
            if len(set(sources)) != len(sources):
                mechanical_errors.append("duplicate_source_id")
            if value is None:
                if claim["basis"] != "unknown":
                    mechanical_errors.append("null_value_basis_not_unknown")
                if claim["quote"] is not None:
                    mechanical_errors.append("null_value_quote_not_null")
            else:
                if claim["basis"] not in ("explicit", "derived"):
                    mechanical_errors.append("nonnull_value_basis_unknown")
                if not sources:
                    mechanical_errors.append("nonnull_value_has_no_source_ids")
                quote = claim["quote"]
                if not isinstance(quote, str) or not normalize(quote):
                    mechanical_errors.append("missing_nonempty_quote")
                else:
                    if len(quote) > 600:
                        mechanical_errors.append("quote_exceeds_600_characters")
                    if normalize(quote) not in source_text:
                        mechanical_errors.append("quote_not_verbatim_in_packet")
                replay = validator_replay(packet, query, claim)
        mechanical_ok = usable and not mechanical_errors
        # Acceptance depends on response/evidence rules, never expected correctness.
        accepted = mechanical_ok and replay["accepted"]
        grades.append(
            {
                "packet": packet["name"],
                "gene": packet["gene"],
                "pmid": packet["pmid"],
                **query,
                "expected": expected,
                "claim": claim,
                "response_usable": usable,
                "value": value,
                "numeric_outcome": outcome,
                "numeric_exact": numeric_exact,
                "model_declared_basis": claim.get("basis") if claim else None,
                "mechanical_source_quote_ok": mechanical_ok,
                "mechanical_errors": mechanical_errors,
                "literal_validator": replay,
                "current_validator_lane_accepted": accepted,
                "accepted_correct": accepted and value == expected,
                "accepted_wrong": accepted and value != expected,
            }
        )
    return {
        "name": plan_cell["name"],
        "packet": packet["name"],
        "arm": plan_cell["arm"],
        "dispatched": dispatched,
        "status": response.get("status") if dispatched else "not_dispatched",
        "finish_reason": finish,
        "schema_valid": schema_valid,
        "request_successful_stop": request_ok,
        "errors": errors,
        "duplicate_query_ids": duplicate_ids,
        "missing_query_ids": missing_ids,
        "extra_query_ids": extra_ids,
        "complete_query_set": not (duplicate_ids or missing_ids or extra_ids),
        "seconds": response.get("seconds") if dispatched else None,
        "api_proxy_usd": response.get("api_proxy_usd") if dispatched else None,
        "accounted_usd": response.get("accounted_usd") if dispatched else 0,
        "queries": grades,
    }


def metrics(rows):
    outcomes = Counter(r["numeric_outcome"] for r in rows)
    basis = {}
    for label in ("explicit", "derived", "unknown", None):
        selected = [r for r in rows if r["model_declared_basis"] == label]
        basis[str(label)] = {
            "claims": len(selected),
            "outcomes": dict(Counter(r["numeric_outcome"] for r in selected)),
            "accepted_correct": sum(r["accepted_correct"] for r in selected),
            "accepted_wrong": sum(r["accepted_wrong"] for r in selected),
        }
    return {
        "planned_queries": len(rows),
        "reference_nonnull": sum(r["expected"] is not None for r in rows),
        "reference_positive": sum(
            r["expected"] is not None and r["expected"] > 0 for r in rows
        ),
        "reference_zero": sum(r["expected"] == 0 for r in rows),
        "reference_unknown": sum(r["expected"] is None for r in rows),
        "outcomes": {
            k: outcomes[k]
            for k in (
                "exact_nonnull",
                "wrong_nonnull",
                "null_miss",
                "correct_abstention",
                "response_failure_or_missing_claim",
            )
        },
        "numeric_exact_including_valid_abstentions": sum(
            r["numeric_exact"] for r in rows
        ),
        "mechanical_source_quote_pass": sum(
            r["mechanical_source_quote_ok"] for r in rows
        ),
        "exact_nonnull_with_mechanical_source_quote_pass": sum(
            r["numeric_outcome"] == "exact_nonnull" and r["mechanical_source_quote_ok"]
            for r in rows
        ),
        "current_validator_accepted_correct": sum(r["accepted_correct"] for r in rows),
        "current_validator_accepted_wrong": sum(r["accepted_wrong"] for r in rows),
        "current_validator_raw_accepts_before_mechanical_gate": sum(
            r["literal_validator"]["accepted"] for r in rows
        ),
        "by_model_declared_basis_not_adjudicated_basis": basis,
    }


def cost_summary(cells):
    seconds = [c["seconds"] for c in cells if c["seconds"] is not None]
    return {
        "planned_requests": len(cells),
        "dispatched": sum(c["dispatched"] for c in cells),
        "successful_stop": sum(c["request_successful_stop"] for c in cells),
        "known_usage_requests": sum(c["api_proxy_usd"] is not None for c in cells),
        "unknown_usage_requests": sum(
            c["dispatched"] and c["api_proxy_usd"] is None for c in cells
        ),
        "returned_api_proxy_usd": sum(c["api_proxy_usd"] or 0 for c in cells),
        "accounted_usd_including_unknown_reserves": sum(
            c["accounted_usd"] or 0 for c in cells
        ),
        "elapsed_request_seconds_sum": sum(seconds),
        "elapsed_request_seconds_median": statistics.median(seconds)
        if seconds
        else None,
    }


def chosen_value(row, lane):
    if lane == "raw_candidates":
        return row["value"] if row["response_usable"] else None
    assert lane == "current_validator"
    return row["value"] if row["current_validator_lane_accepted"] else None


def compare_overlay(first_rows, added_rows, lane):
    by_key = {(r["packet"], r["id"]): r for r in added_rows}
    rows = []
    for first in first_rows:
        added = by_key[(first["packet"], first["id"])]
        old, new = chosen_value(first, lane), chosen_value(added, lane)
        value = old if old is not None else new
        expected = first["expected"]
        rows.append(
            {
                "packet": first["packet"],
                "pmid": first["pmid"],
                "id": first["id"],
                "variant": first["variant"],
                "field": first["field"],
                "expected": expected,
                "first_value": old,
                "overlay_value": new,
                "selected_value": value,
                "selected_from": "first"
                if old is not None
                else "overlay"
                if new is not None
                else "none",
                "conflict_preserving_first": old is not None
                and new is not None
                and old != new,
                "new_correct_nonnull": old is None
                and new is not None
                and new == expected,
                "new_wrong_nonnull": old is None
                and new is not None
                and new != expected,
                "final_correct_nonnull": value is not None and value == expected,
                "final_wrong_nonnull": value is not None and value != expected,
            }
        )
    sums = {
        k: sum(r[k] for r in rows)
        for k in (
            "conflict_preserving_first",
            "new_correct_nonnull",
            "new_wrong_nonnull",
            "final_correct_nonnull",
            "final_wrong_nonnull",
        )
    }
    sums["net_additional_correct_minus_wrong"] = (
        sums["new_correct_nonnull"] - sums["new_wrong_nonnull"]
    )
    return {"lane": lane, "planned_queries": len(rows), **sums, "queries": rows}


def build_scores():
    plan, lock, packets, references, responses = verify_inputs()
    cells = [
        grade_cell(
            p, packets[p["packet"]], references[p["packet"]], responses[p["name"]]
        )
        for p in plan["plans"]
    ]
    by_arm = {a: [c for c in cells if c["arm"] == a] for a in ARMS}
    arm_rows = {a: [r for c in by_arm[a] for r in c["queries"]] for a in ARMS}
    summaries = {
        a: {**metrics(arm_rows[a]), "operations": cost_summary(by_arm[a])} for a in ARMS
    }
    overlays = {}
    for arm in ("grok_repeat", "astra_low"):
        overlays[arm] = {}
        for lane in ("raw_candidates", "current_validator"):
            comparison = compare_overlay(arm_rows["grok_first"], arm_rows[arm], lane)
            cost = summaries[arm]["operations"]["returned_api_proxy_usd"]
            gain = comparison["new_correct_nonnull"]
            net = comparison["net_additional_correct_minus_wrong"]
            comparison["incremental_returned_api_proxy_usd"] = cost
            comparison["incremental_accounted_usd"] = summaries[arm]["operations"][
                "accounted_usd_including_unknown_reserves"
            ]
            comparison["proxy_usd_per_added_correct_nonnull"] = (
                cost / gain if gain else None
            )
            comparison["proxy_usd_per_net_added_correct_minus_wrong"] = (
                cost / net if net > 0 else None
            )
            overlays[arm][lane] = comparison
    indexes = {a: {(r["packet"], r["id"]): r for r in arm_rows[a]} for a in ARMS}
    unique = []
    for key, astra in indexes["astra_low"].items():
        groks = [indexes[a][key] for a in ("grok_first", "grok_repeat")]
        if astra["numeric_outcome"] == "exact_nonnull" and all(
            g["numeric_outcome"] != "exact_nonnull" for g in groks
        ):
            unique.append(
                {
                    "packet": astra["packet"],
                    "pmid": astra["pmid"],
                    "id": astra["id"],
                    "variant": astra["variant"],
                    "field": astra["field"],
                    "expected": astra["expected"],
                    "astra_claim": astra["claim"],
                    "grok_first_claim": groks[0]["claim"],
                    "grok_repeat_claim": groks[1]["claim"],
                    "grok_first_outcome": groks[0]["numeric_outcome"],
                    "grok_repeat_outcome": groks[1]["numeric_outcome"],
                    "mechanical_source_quote_ok": astra["mechanical_source_quote_ok"],
                    "current_validator_accepted": astra[
                        "current_validator_lane_accepted"
                    ],
                    "requires_semantic_source_adjudication": True,
                }
            )
    return {
        "classification": plan["classification"],
        "limitations": [
            "Reference numbers are investigator readings of selected opened packets, not an independent human benchmark.",
            "Clinical endpoints and populations are those explicitly dictated by each packet.",
            "Mechanical source checks verify identifier membership and verbatim quotation, not semantic support or completeness.",
            "Numeric reference equality is never used to accept a claim; accepted correctness is assessed after acceptance.",
            "Model-declared explicit/derived labels are not independent adjudications of derivation.",
            "Only explicit nonnull claims are replayed through the unchanged literal validator; no database writes occur.",
            "Missing, failed, malformed or undispatched outputs never earn correct-abstention credit.",
            "Raw-candidate overlays are diagnostics, not validated production merges; existing nonnull values are preserved.",
            "Cost-per-gain uses returned-usage proxies; unknown reserves and full accounted totals are shown separately.",
            "Astra-only numeric candidates require source-semantic adjudication and replication before a selective-pilot claim.",
        ],
        "integrity": {
            "plan_sha256": sha(HERE / "plan.json"),
            "output_lock_sha256": sha(HERE / "outputs_locked.json"),
            "reference_sha256": sha(HERE / "reference_values.json"),
            "score_code_sha256": sha(Path(__file__)),
            "literal_validator_code_sha256": sha(ROOT / "pipeline/count_recovery.py"),
            "locked_unix": lock["locked_unix"],
            "concurrency_amendment_sha256": lock.get("concurrency_amendment_sha256"),
            "source_packet_hashes": plan["packet_sha256"],
        },
        "panel": {
            "papers": len(packets),
            "queries_per_arm": sum(len(p["queries"]) for p in packets.values()),
            "reference_nonnull": sum(
                v is not None for r in references.values() for v in r.values()
            ),
            "reference_positive": sum(
                v is not None and v > 0 for r in references.values() for v in r.values()
            ),
            "reference_zero": sum(
                v == 0 for r in references.values() for v in r.values()
            ),
            "reference_unknown": sum(
                v is None for r in references.values() for v in r.values()
            ),
        },
        "arms": summaries,
        "overlays": overlays,
        "astra_unique_numeric_exact_beyond_both_groks": unique,
        "astra_unique_numeric_exact_papers": sorted({r["packet"] for r in unique}),
        "cells": cells,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify exact offline replay without writing evidence",
    )
    args = parser.parse_args()
    result = build_scores()
    target = HERE / "source_scores.json"
    if args.check or target.exists():
        assert read(target) == result, (
            "Existing scores differ; do not overwrite locked evidence"
        )
    else:
        target.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(
        json.dumps(
            {
                "panel": result["panel"],
                "arms": result["arms"],
                "astra_unique_numeric_exact_papers": result[
                    "astra_unique_numeric_exact_papers"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
