#!/usr/bin/env python3
"""Per-model API cost proxy for one production evaluation run, from its LLM traces.

Reads every ``llm_traces/trace_index.jsonl`` under ``<run>/production_runs/<GENE>/*``,
opens each ``llm_call`` record, and sums prompt/completion tokens per model.
Prices are the repository's dated public list-price proxy
(``RATES`` in docs/evidence/model_routing_20260906/budget_guard.py, USD per
million input/output tokens). This is a planning proxy, not an invoice; calls
whose record carries no usage are counted separately and cost nothing here.

Usage: cost_ledger.py <run_dir> <budget_ceiling_usd> <out.json>
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

RATES = {  # USD per million (input, output); mirrors budget_guard.py
    "gpt-6-astra": (10.0, 50.0),
    "grok-4.6": (2.0, 6.0),
    "grok-4.3": (1.25, 2.5),
    "gpt-5.6-sol": (5.0, 30.0),
    "gpt-5.6-luna": (0.20, 1.20),
    "kimi-k2.6": (0.95, 4.0),
}


def rate_for(model: str):
    name = model.lower()
    for key, rate in RATES.items():
        if key in name:
            return rate
    return None


def main(run_dir: Path, ceiling: float, out: Path) -> None:
    per_model: dict[str, dict] = defaultdict(
        lambda: {
            "calls": 0,
            "input_tokens": 0,
            "output_tokens": 0,
            "calls_without_usage": 0,
            "failed_calls": 0,
        }
    )
    unpriced = set()
    for index in sorted(
        run_dir.glob("production_runs/*/*/llm_traces/trace_index.jsonl")
    ):
        root = index.parent
        for line in index.read_text().splitlines():
            if not line.strip():
                continue
            entry = json.loads(line)
            if entry.get("record_type") != "llm_call":
                continue
            path = root / entry["path"]
            if not path.suffix:
                matches = list(path.parent.glob(path.name + "*"))
                if not matches:
                    continue
                path = matches[0]
            record = json.loads(path.read_text())
            model = str((record.get("context") or {}).get("model") or "unknown")
            response = record.get("response") or {}
            usage = response.get("usage") or {}
            bucket = per_model[model]
            bucket["calls"] += 1
            if response.get("success") is False or response.get("error"):
                bucket["failed_calls"] += 1
            if not usage:
                bucket["calls_without_usage"] += 1
                continue
            bucket["input_tokens"] += int(usage.get("prompt_tokens") or 0)
            bucket["output_tokens"] += int(usage.get("completion_tokens") or 0)
    rows = []
    total = 0.0
    for model, bucket in sorted(per_model.items()):
        rate = rate_for(model)
        usd = None
        if rate:
            usd = round(
                bucket["input_tokens"] / 1e6 * rate[0]
                + bucket["output_tokens"] / 1e6 * rate[1],
                4,
            )
            total += usd
        else:
            unpriced.add(model)
        rows.append({"model": model, **bucket, "usd": usd})
    payload = {
        "run_id": run_dir.name,
        "pricing_basis": (
            "repository dated API list-price proxy (RATES in "
            "docs/evidence/model_routing_20260906/budget_guard.py) applied to per-call "
            "trace usage; actual invoice may differ; CLI reviewer costs excluded"
        ),
        "per_model": rows,
        "unpriced_models": sorted(unpriced),
        "api_proxy_usd": round(total, 4),
        "budget_ceiling_usd": ceiling,
        "remaining_usd": round(ceiling - total, 4),
    }
    out.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in payload.items() if k != "pricing_basis"}, indent=2)
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]), float(sys.argv[2]), Path(sys.argv[3]))
