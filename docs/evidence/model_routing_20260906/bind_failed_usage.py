"""Bind existing failed API traces to a finished clinical overlay before locking."""

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from benchmarks.codex_paper_eval.run_eval import indexed_paper_call_usage


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def documented_budget_non_dispatch(paper, audit, folder, budget):
    """Identify an unissued guarded call; never turn a timeout into zero usage."""
    error = str(audit.get("error") or "")
    if not error.startswith(
        "RuntimeError: Campaign API budget prevents call to azure_ai/gpt-6-astra; reservation $"
    ):
        return None
    gene, pmid = paper["gene"], str(paper["pmid"])
    if any(
        row["arm"] == "grok43_astra_clinical_verified"
        and row.get("gene") == gene
        and str(row.get("pmid")) == pmid
        for row in budget["calls"]
    ):
        return None
    trace_root = (folder / "llm_traces").resolve()
    matches = []
    try:
        for line in (trace_root / "trace_index.jsonl").read_text().splitlines():
            entry = json.loads(line)
            context = entry.get("context") or {}
            path = (trace_root / entry["path"]).resolve()
            if not path.is_relative_to(trace_root) or sha(path) != entry["sha256"]:
                return None
            record = json.loads(path.read_text())
            actual = record.get("context") or {}
            if (
                (actual.get("gene"), str(actual.get("pmid")))
                != (context.get("gene"), str(context.get("pmid")))
                or record.get("trace_id") != entry.get("trace_id")
                or record.get("record_type") != entry.get("record_type")
            ):
                return None
            if (actual.get("gene"), str(actual.get("pmid"))) != (gene, pmid):
                continue
            if record.get("record_type") == "llm_call":
                return None
            if (
                record.get("record_type") == "decision_event"
                and actual.get("stage") == "paper_curation"
                and record.get("event")
                == {
                    "type": "paper_curation",
                    "data": {"status": "failed", "error": error},
                }
            ):
                if not any(
                    ref.get("trace_id") == entry["trace_id"]
                    and ref.get("path") == entry["path"]
                    and ref.get("sha256") == entry["sha256"]
                    for ref in paper.get("llm_trace_refs", [])
                ):
                    return None
                matches.append(entry)
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return matches[0] if len(matches) == 1 else None


def main():
    assert not (HERE / "failed_usage_binding.json").exists(), (
        "Preserve the original binding receipt"
    )
    folder = (
        ROOT
        / "benchmarks/codex_paper_eval/runs/20260906_model12_grok43_astra_clinical_verified"
    )
    assert not (folder / "LOCK.json").exists(), "Never mutate a locked result"
    final_runs = ["grok43", "astra_medium_verified", "grok43_astra_clinical_verified"]
    assert not any(
        (folder.parent / ("20260906_model12_" + name) / "report.json").exists()
        for name in final_runs
    ), "All arms must lock before scoring"
    path = folder / "predictions.json"
    data = json.loads(path.read_text())
    assert data.get("completed_at") and len(data["papers"]) == 12
    before_sha = sha(path)
    scientific = json.dumps(
        [
            (p["gene"], p["pmid"], p["variants"], p["comparison_variants"])
            for p in data["papers"]
        ],
        sort_keys=True,
    )
    budget_path = HERE / "budget.json"
    budget_sha = sha(budget_path)
    budget = json.loads(budget_path.read_text())
    assert not any(row["status"] == "reserved" for row in budget["calls"]), (
        "Finish every campaign call before binding usage"
    )
    changed = []
    not_dispatched = []
    audits = {(a["gene"], str(a["pmid"])): a for a in data["clinical_reader_audits"]}
    for paper in data["papers"]:
        usage = paper["token_usage"]
        if usage["telemetry_available"]:
            continue
        blocked = documented_budget_non_dispatch(
            paper, audits[(paper["gene"], str(paper["pmid"]))], folder, budget
        )
        if blocked is not None:
            usage.update(
                telemetry_available=True,
                input_tokens=0,
                output_tokens=0,
                total_tokens=0,
                status="budget_not_dispatched",
                note="The campaign guard refused before dispatch. A hash-bound decision records the refusal; no API call or ledger reservation exists for this paper's clinical arm.",
            )
            not_dispatched.append(
                {
                    "gene": paper["gene"],
                    "pmid": paper["pmid"],
                    "budget_decision_trace_id": blocked["trace_id"],
                }
            )
            continue
        indexed = indexed_paper_call_usage(paper, folder / "llm_traces")
        assert indexed is not None
        matches = list(indexed["unknown_calls"].values())
        assert matches, "Unknown usage requires an authentic API failure"
        assert any(
            row["arm"] == "grok43_astra_clinical_verified"
            and row.get("gene") == paper["gene"]
            and str(row.get("pmid")) == str(paper["pmid"])
            and row["status"] == "failed_usage_unknown_reservation_retained"
            and row["reserved_usd"] > 0
            for row in budget["calls"]
        ), "Retain the failed-call reservation"
        for entry in matches:
            assert sha(folder / "llm_traces" / entry["path"]) == entry["sha256"]
            record = json.loads((folder / "llm_traces" / entry["path"]).read_text())
            assert (
                record["response"]["success"] is False
                and record["response"]["usage"] is None
            )
            if entry["trace_id"] not in {
                r.get("trace_id") for r in paper["llm_trace_refs"]
            }:
                paper["llm_trace_refs"].append(entry)
        usage.update(
            status="unknown_failed_call",
            unknown_failed_call_trace_ids=[e["trace_id"] for e in matches],
            known_usage=indexed["known_usage"],
        )
        changed.append(
            {
                "gene": paper["gene"],
                "pmid": paper["pmid"],
                "failed_trace_ids": usage["unknown_failed_call_trace_ids"],
            }
        )
    assert scientific == json.dumps(
        [
            (p["gene"], p["pmid"], p["variants"], p["comparison_variants"])
            for p in data["papers"]
        ],
        sort_keys=True,
    )
    data["baseline_token_usage"] = data["token_usage"]
    known = {
        f: sum(
            (
                p["token_usage"]
                if p["token_usage"]["telemetry_available"]
                else p["token_usage"]["known_usage"]
            ).get(f)
            or 0
            for p in data["papers"]
        )
        for f in ["input_tokens", "output_tokens", "total_tokens"]
    }
    data["token_usage"] = {
        "telemetry_available": not bool(changed),
        **{f: None if changed else v for f, v in known.items()},
        "known_incremental_usage": known,
        "known_usage": known,
        "unknown_usage_papers": len(changed),
        "note": "Incremental clinical usage; aggregate totals are unknown when a failed API call returned no usage. Baseline usage remains separate. Campaign ledger retains cost reservations.",
    }
    assert sha(budget_path) == budget_sha, "Budget changed during metadata binding"
    path.write_text(json.dumps(data, indent=2) + "\n")
    receipt = {
        "budget_sha256": budget_sha,
        "campaign_calls_in_flight": 0,
        "binder_sha256": sha(Path(__file__)),
        "before_predictions_sha256": before_sha,
        "after_predictions_sha256": sha(path),
        "scientific_predictions_sha256_unchanged": hashlib.sha256(
            scientific.encode()
        ).hexdigest(),
        "attached_failure_references": changed,
        "budget_not_dispatched": not_dispatched,
        "new_score_artifacts_absent_on_all_final_arms": True,
        "operator_attestation": "No new gold scores were inspected before binding; this is opened-paper calibration, not investigator-blind research.",
    }
    (HERE / "failed_usage_binding.json").write_text(
        json.dumps(receipt, indent=2) + "\n"
    )
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
