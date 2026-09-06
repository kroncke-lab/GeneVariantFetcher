"""Summarize call outcomes without treating API success as extraction success."""

import collections
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
RUNS = ROOT / "benchmarks/codex_paper_eval/runs"
ARMS = ["grok43", "astra_medium_verified", "grok43_astra_clinical_verified"]


def main():
    budget = json.loads((HERE / "budget.json").read_text())
    assert not any(r["status"] == "reserved" for r in budget["calls"])
    output = {}
    for arm in ARMS:
        run = RUNS / ("20260906_model12_" + arm)
        assert (run / "LOCK.json").exists()
        # Read original, source-bound production traces once. The clinical arm
        # copies baseline traces under another folder; keep those out of its
        # incremental readout by requiring the PaperCurator component.
        clinical = arm.endswith("clinical_verified")
        paths = (
            (run / "llm_traces").rglob("*.json")
            if clinical
            else (run / "production_runs").glob("*/*/llm_traces/*/*/*.json")
        )
        rows = []
        for path in paths:
            record = json.loads(path.read_text())
            if record.get("record_type") != "llm_call":
                continue
            context = record.get("context", {})
            if clinical and context.get("stage") != "paper_curation":
                continue
            if not clinical and context.get("stage") != "paper_variant_extraction":
                continue
            response = record["response"]
            usage = response.get("usage")
            envelope = response.get("envelope") or {}
            choices = envelope.get("choices") or []
            call = next(
                (r for r in budget["calls"] if r.get("trace_id") == record["trace_id"]),
                None,
            )
            rows.append(
                {
                    "gene": context.get("gene"),
                    "pmid": context.get("pmid"),
                    "model": context.get("model"),
                    "operation": context.get("operation"),
                    "attempt": context.get("attempt"),
                    "trace_id": record["trace_id"],
                    "trace_path": str(path.relative_to(ROOT)),
                    "api_success": response.get("success"),
                    "seconds": response.get("duration_seconds"),
                    "finish_reason": choices[0].get("finish_reason")
                    if choices
                    else None,
                    "visible_output_chars": len(response.get("output_text") or ""),
                    "reasoning_tokens": (
                        usage.get("completion_tokens_details") or {}
                    ).get("reasoning_tokens")
                    if usage
                    else None,
                    "usage": usage,
                    "error": response.get("error"),
                    "sdk_retry_count": call.get("sdk_retry_count") if call else None,
                    "returned_proxy_usd": call.get("api_proxy_usd") if call else None,
                    "unknown_retry_reserve_usd": call.get("unknown_retry_reserve_usd")
                    if call
                    else None,
                }
            )
        rows.sort(key=lambda r: (r["gene"], r["pmid"], r["trace_id"]))
        output[arm] = {
            "scope": "Incremental clinical calls"
            if clinical
            else "Primary extraction, its empty-output retries and JSON repair calls",
            "calls": len(rows),
            "papers_reaching_stage": sorted(
                {r["gene"] + ":" + r["pmid"] for r in rows}
            ),
            "finish_reasons": dict(
                collections.Counter(str(r["finish_reason"]) for r in rows)
            ),
            "empty_visible_returns": sum(
                r["api_success"] and r["visible_output_chars"] == 0 for r in rows
            ),
            "known_reasoning_tokens": sum(r["reasoning_tokens"] or 0 for r in rows),
            "call_records": rows,
        }
    (HERE / "operational_summary.json").write_text(json.dumps(output, indent=2) + "\n")
    print(
        json.dumps(
            {
                a: {k: v for k, v in d.items() if k != "call_records"}
                for a, d in output.items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
