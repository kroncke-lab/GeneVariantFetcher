#!/usr/bin/env python3
"""Shadow a Tier-2 relevance model on a stratified historical workload.

This evaluator never mutates the source run and never changes production model
routing. It replays title/abstract inputs through ``InternFilter`` so the model,
effort, structured-output validity, deterministic fail-open behavior, tokens,
and latency can be compared with the source run's exact traces.

The default cohort is intentionally diagnostic rather than a gold standard:

* ``productive_positive`` papers come from a pinned review cohort whose source
  run extracted at least one target-gene variant;
* ``high_confidence_negative`` papers were rejected by the historical model at
  >=0.95 confidence; and
* ``fail_open_boundary`` papers were rejected by the historical model but kept
  by the deterministic target-gene/cohort guard.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pipeline.filters import InternFilter
from utils.llm_trace import configure_llm_tracing
from utils.models import Paper


DEFAULT_SOURCE_RUN = REPO / "results/BMPR2/20260807_163246"
DEFAULT_REVIEW_MANIFEST = (
    REPO
    / "benchmarks/curated_extraction_eval/review_pmids_50_bmpr2_20260808"
    / "BMPR2_manifest.csv"
)
DEFAULT_OUTDIR = REPO / "validation_runs/20260810_bmpr2_tier2_luna_max"
DEFAULT_DISEASE = (
    "pulmonary arterial hypertension (aliases: heritable pulmonary arterial "
    "hypertension, familial pulmonary arterial hypertension, idiopathic "
    "pulmonary arterial hypertension, primary pulmonary hypertension, PAH, "
    "HPAH, FPAH, IPAH, PPH, pulmonary hypertension, pulmonary veno-occlusive "
    "disease, PVOD)"
)


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def _write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _stable_sample(
    rows: Iterable[dict[str, Any]], count: int, *, seed: str
) -> list[dict[str, Any]]:
    ranked = sorted(
        rows,
        key=lambda row: hashlib.sha256(
            f"{seed}:{row['pmid']}".encode("utf-8")
        ).hexdigest(),
    )
    if len(ranked) < count:
        raise ValueError(f"requested {count} cases but only {len(ranked)} qualify")
    return ranked[:count]


def _load_review_pmids(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [
            str(row["pmid"]).strip()
            for row in csv.DictReader(handle)
            if str(row.get("pmid") or "").strip()
        ]


def build_cohort(
    *,
    source_run: Path,
    review_manifest: Path,
    positive_count: int,
    negative_count: int,
    boundary_count: int,
    seed: str,
) -> list[dict[str, Any]]:
    """Build a deterministic, disjoint diagnostic cohort."""
    progress_rows = _read_jsonl(source_run / "pmid_status/filter_progress.jsonl")
    progress = {str(row["pmid"]): row for row in progress_rows}
    abstract_dir = source_run / "abstract_json"

    def has_input(pmid: str) -> bool:
        path = abstract_dir / f"{pmid}.json"
        if not path.exists():
            return False
        record = json.loads(path.read_text(encoding="utf-8"))
        metadata = record.get("metadata") or {}
        return bool(metadata.get("title") and record.get("abstract"))

    positive_candidates = [
        {"pmid": pmid, "group": "productive_positive"}
        for pmid in _load_review_pmids(review_manifest)
        if pmid in progress and has_input(pmid)
    ]
    positives = _stable_sample(positive_candidates, positive_count, seed=seed)
    positive_pmids = {row["pmid"] for row in positives}

    negative_candidates = [
        {"pmid": pmid, "group": "high_confidence_negative"}
        for pmid, row in progress.items()
        if pmid not in positive_pmids
        and row.get("final_decision") == "FAIL"
        and row.get("stage") == "tier2"
        and float((row.get("tier2") or {}).get("confidence") or 0) >= 0.95
        and has_input(pmid)
    ]
    negatives = _stable_sample(negative_candidates, negative_count, seed=seed)
    excluded = positive_pmids | {row["pmid"] for row in negatives}

    boundary_candidates = [
        {"pmid": pmid, "group": "fail_open_boundary"}
        for pmid, row in progress.items()
        if pmid not in excluded
        and str((row.get("tier2") or {}).get("reason") or "").startswith(
            "Fail-open override:"
        )
        and has_input(pmid)
    ]
    boundaries = _stable_sample(boundary_candidates, boundary_count, seed=seed)

    cohort: list[dict[str, Any]] = []
    for selected in positives + negatives + boundaries:
        pmid = selected["pmid"]
        historical = progress[pmid]
        record = json.loads((abstract_dir / f"{pmid}.json").read_text(encoding="utf-8"))
        metadata = record.get("metadata") or {}
        tier2 = historical.get("tier2") or {}
        historical_reason = str(tier2.get("reason") or historical.get("reason") or "")
        cohort.append(
            {
                **selected,
                "title": metadata.get("title") or "",
                "abstract": record.get("abstract") or "",
                "historical_final_decision": str(
                    tier2.get("decision") or historical.get("final_decision") or ""
                ).lower(),
                "historical_raw_model_decision": (
                    "fail"
                    if historical_reason.startswith("Fail-open override:")
                    else str(
                        tier2.get("decision") or historical.get("final_decision") or ""
                    ).lower()
                ),
                "historical_confidence": tier2.get("confidence"),
                "historical_reason": historical_reason,
            }
        )
    return sorted(cohort, key=lambda row: (row["group"], row["pmid"]))


def _trace_records(trace_root: Path) -> list[dict[str, Any]]:
    index_path = trace_root / "trace_index.jsonl"
    if not index_path.exists():
        return []
    records: list[dict[str, Any]] = []
    for entry in _read_jsonl(index_path):
        path = trace_root / str(entry["path"])
        if path.exists():
            records.append(json.loads(path.read_text(encoding="utf-8")))
    return records


def _telemetry_by_pmid(
    records: Iterable[dict[str, Any]], selected_pmids: set[str]
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    calls: dict[str, dict[str, Any]] = {}
    decisions: dict[str, dict[str, Any]] = {}
    for record in records:
        context = record.get("context") or {}
        pmid = str(context.get("pmid") or "")
        if pmid not in selected_pmids:
            continue
        if record.get("record_type") == "llm_call" and context.get("stage") == (
            "tier2_relevance_filter"
        ):
            response = record.get("response") or {}
            usage = (
                response.get("usage")
                or (response.get("envelope") or {}).get("usage")
                or {}
            )
            entry = calls.setdefault(
                pmid,
                {
                    "calls": 0,
                    "input_tokens": 0,
                    "output_tokens": 0,
                    "reasoning_tokens": 0,
                    "total_tokens": 0,
                    "duration_seconds": 0.0,
                    "requested_models": [],
                    "response_models": [],
                    "reasoning_efforts": [],
                    "failures": 0,
                },
            )
            entry["calls"] += 1
            entry["input_tokens"] += int(
                usage.get("prompt_tokens") or usage.get("input_tokens") or 0
            )
            entry["output_tokens"] += int(
                usage.get("completion_tokens") or usage.get("output_tokens") or 0
            )
            entry["total_tokens"] += int(usage.get("total_tokens") or 0)
            completion_details = usage.get("completion_tokens_details") or {}
            entry["reasoning_tokens"] += int(
                completion_details.get("reasoning_tokens") or 0
            )
            entry["duration_seconds"] += float(response.get("duration_seconds") or 0)
            requested = str((record.get("request") or {}).get("requested_model") or "")
            response_model = str(((response.get("envelope") or {}).get("model") or ""))
            request_payload = (record.get("request") or {}).get("payload") or {}
            reasoning_effort = str(
                request_payload.get("reasoning_effort")
                or (request_payload.get("reasoning") or {}).get("effort")
                or ""
            )
            if requested and requested not in entry["requested_models"]:
                entry["requested_models"].append(requested)
            if response_model and response_model not in entry["response_models"]:
                entry["response_models"].append(response_model)
            if reasoning_effort and reasoning_effort not in entry["reasoning_efforts"]:
                entry["reasoning_efforts"].append(reasoning_effort)
            if not response.get("success"):
                entry["failures"] += 1
        elif record.get("record_type") == "decision_event":
            event = record.get("event") or {}
            if event.get("type") == "tier2_relevance_decision":
                decisions[pmid] = event.get("data") or {}
    return calls, decisions


def _percentile(values: list[float], quantile: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, math.ceil(quantile * len(ordered)) - 1))
    return ordered[index]


def _aggregate_telemetry(per_pmid: dict[str, dict[str, Any]]) -> dict[str, Any]:
    durations = [
        float(row["duration_seconds"]) for row in per_pmid.values() if row.get("calls")
    ]
    return {
        "papers": len(per_pmid),
        "calls": sum(int(row["calls"]) for row in per_pmid.values()),
        "input_tokens": sum(int(row["input_tokens"]) for row in per_pmid.values()),
        "output_tokens": sum(int(row["output_tokens"]) for row in per_pmid.values()),
        "reasoning_tokens": sum(
            int(row["reasoning_tokens"]) for row in per_pmid.values()
        ),
        "total_tokens": sum(int(row["total_tokens"]) for row in per_pmid.values()),
        "summed_call_seconds": round(sum(durations), 6),
        "median_call_seconds": (
            round(statistics.median(durations), 6) if durations else None
        ),
        "p95_call_seconds": (
            round(float(_percentile(durations, 0.95)), 6) if durations else None
        ),
        "provider_failures": sum(int(row["failures"]) for row in per_pmid.values()),
        "requested_models": sorted(
            {model for row in per_pmid.values() for model in row["requested_models"]}
        ),
        "response_models": sorted(
            {model for row in per_pmid.values() for model in row["response_models"]}
        ),
        "reasoning_efforts": sorted(
            {effort for row in per_pmid.values() for effort in row["reasoning_efforts"]}
        ),
    }


def _quality_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    def rate(group_rows: list[dict[str, Any]], key: str, expected: str) -> float | None:
        if not group_rows:
            return None
        return sum(row[key] == expected for row in group_rows) / len(group_rows)

    by_group: dict[str, Any] = {}
    for group in sorted({row["group"] for row in rows}):
        group_rows = [row for row in rows if row["group"] == group]
        expected = "fail" if group == "high_confidence_negative" else "pass"
        by_group[group] = {
            "n": len(group_rows),
            "expected_final_decision": expected,
            "luna_expected_rate": rate(group_rows, "luna_final_decision", expected),
            "historical_expected_rate": rate(
                group_rows, "historical_final_decision", expected
            ),
            "luna_raw_passes": sum(
                row["luna_raw_model_decision"] == "pass" for row in group_rows
            ),
            "luna_raw_fails": sum(
                row["luna_raw_model_decision"] == "fail" for row in group_rows
            ),
            "luna_fail_open_overrides": sum(
                bool(row.get("luna_fail_open_target_gene_signal")) for row in group_rows
            ),
        }

    return {
        "n": len(rows),
        "final_decision_agreement": sum(
            row["luna_final_decision"] == row["historical_final_decision"]
            for row in rows
        )
        / len(rows),
        "raw_model_decision_agreement": sum(
            row["luna_raw_model_decision"] == row["historical_raw_model_decision"]
            for row in rows
        )
        / len(rows),
        "fail_open_call_errors": sum(bool(row.get("luna_call_error")) for row in rows),
        "groups": by_group,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", type=Path, default=DEFAULT_SOURCE_RUN)
    parser.add_argument("--review-manifest", type=Path, default=DEFAULT_REVIEW_MANIFEST)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--model", default="azure_ai/gpt-5.6-luna")
    parser.add_argument("--reasoning-effort", default="max")
    parser.add_argument("--positive-count", type=int, default=50)
    parser.add_argument("--negative-count", type=int, default=50)
    parser.add_argument("--boundary-count", type=int, default=50)
    parser.add_argument("--seed", default="bmpr2-tier2-shadow-v1")
    parser.add_argument("--workers", type=int, default=8)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    source_run = args.source_run.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cohort = build_cohort(
        source_run=source_run,
        review_manifest=args.review_manifest.resolve(),
        positive_count=args.positive_count,
        negative_count=args.negative_count,
        boundary_count=args.boundary_count,
        seed=args.seed,
    )
    _write_json(outdir / "cohort.json", cohort)
    selected_pmids = {row["pmid"] for row in cohort}

    historical_trace_root = source_run / "llm_traces"
    historical_calls, historical_decisions = _telemetry_by_pmid(
        _trace_records(historical_trace_root), selected_pmids
    )

    trace_root = configure_llm_tracing(
        outdir / "llm_traces",
        run_id=f"{args.model}@{args.reasoning_effort}",
    )
    classifier = InternFilter(
        model=args.model,
        reasoning_effort=args.reasoning_effort,
        confidence_threshold=0.3,
        disease=DEFAULT_DISEASE,
    )

    def classify(case: dict[str, Any]) -> dict[str, Any]:
        started = time.monotonic()
        result = classifier.filter(
            Paper(
                pmid=case["pmid"],
                title=case["title"],
                abstract=case["abstract"],
                gene_symbol="BMPR2",
            )
        )
        return {
            "pmid": case["pmid"],
            "luna_final_decision": result.decision.value,
            "luna_confidence": result.confidence,
            "luna_reason": result.reason,
            "luna_fail_open_target_gene_signal": bool(
                (result.metadata or {}).get("fail_open_target_gene_signal")
            ),
            "luna_call_error": (result.metadata or {}).get("error"),
            "worker_wall_seconds": round(time.monotonic() - started, 6),
        }

    shadow_results: dict[str, dict[str, Any]] = {}
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(classify, case): case["pmid"] for case in cohort}
        for completed, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            shadow_results[result["pmid"]] = result
            print(f"[{completed}/{len(cohort)}] PMID {result['pmid']}", flush=True)
    wall_seconds = time.monotonic() - started

    new_records = _trace_records(trace_root)
    luna_calls, luna_decisions = _telemetry_by_pmid(new_records, selected_pmids)

    output_rows: list[dict[str, Any]] = []
    for case in cohort:
        pmid = case["pmid"]
        historical_decision = historical_decisions.get(pmid) or {}
        luna_decision = luna_decisions.get(pmid) or {}
        output_rows.append(
            {
                **case,
                **shadow_results[pmid],
                "historical_raw_model_decision": str(
                    historical_decision.get("model_decision")
                    or case["historical_raw_model_decision"]
                ).lower(),
                "luna_raw_model_decision": str(
                    luna_decision.get("model_decision")
                    or shadow_results[pmid]["luna_final_decision"]
                ).lower(),
                "historical_telemetry": historical_calls.get(pmid) or {},
                "luna_telemetry": luna_calls.get(pmid) or {},
            }
        )

    results_path = outdir / "results.jsonl"
    with results_path.open("w", encoding="utf-8") as handle:
        for row in output_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    summary = {
        "source_run": str(source_run),
        "model": args.model,
        "reasoning_effort": args.reasoning_effort,
        "workers": args.workers,
        "wall_seconds": round(wall_seconds, 6),
        "cohort_sha256": hashlib.sha256(
            (outdir / "cohort.json").read_bytes()
        ).hexdigest(),
        "quality": _quality_summary(output_rows),
        "historical_sol": _aggregate_telemetry(historical_calls),
        "luna": _aggregate_telemetry(luna_calls),
        "limitations": [
            "This is a shadow diagnostic, not a manual relevance gold standard.",
            "Productive positives have extracted variants; negative labels inherit the historical Sol decision.",
            "The source Sol run omitted reasoning effort; Luna telemetry records the deployment-effective effort.",
            "Token and latency telemetry are exact, but Azure contract pricing is not encoded in traces.",
        ],
    }
    _write_json(outdir / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
