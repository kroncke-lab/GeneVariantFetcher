#!/usr/bin/env python3
"""Export usage-only calibration evidence for tests on a fresh checkout.

Requires the original ignored operator traces. Does not change the SHA-pinned
cost profile or registries, and never exports prompts or response text.
"""

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export(profile_path: Path, destination: Path) -> dict:
    profile = json.loads(profile_path.read_text())
    result = {
        "schema_version": 1,
        "purpose": "Portable usage-only receipts from local calibration traces; no prompts, responses, credentials or paper text. Original calibration/profile bytes remain unchanged.",
        "calibration_sha256": sha(profile_path),
        "calibrations": {},
    }
    for gene, calibration in profile["calibrations"].items():
        source = ROOT / calibration["source"]
        if source.name == "predictions.json":
            continue  # These complete locked predictions are already tracked.
        if sha(source) != calibration["source_sha256"]:
            raise ValueError(f"Calibration source changed: {source}")
        calls = []
        for path in sorted(source.parent.rglob("*.json")):
            if path == source:
                continue
            try:
                record = json.loads(path.read_text())
            except (OSError, ValueError):
                continue
            if not isinstance(record, dict) or record.get("record_type") != "llm_call":
                continue
            response = record.get("response")
            context = record.get("context")
            response = {} if response is None else response
            context = {} if context is None else context
            if not isinstance(response, dict) or not isinstance(context, dict):
                raise ValueError(f"Malformed llm_call response/context: {path}")
            usage = response.get("usage")
            usage = {} if usage is None else usage
            if not isinstance(usage, dict):
                raise ValueError(f"Malformed llm_call usage: {path}")
            input_tokens = int(
                usage.get("prompt_tokens") or usage.get("input_tokens") or 0
            )
            total_tokens = int(usage.get("total_tokens") or 0)
            calls.append(
                {
                    "trace_sha256": sha(path),
                    "model": str(context.get("model") or "unknown"),
                    "input_tokens": input_tokens,
                    "output_tokens": max(total_tokens - input_tokens, 0),
                }
            )
        actual = {}
        for call in calls:
            counts = actual.setdefault(
                call["model"], dict(calls=0, input_tokens=0, output_tokens=0)
            )
            counts["calls"] += 1
            counts["input_tokens"] += call["input_tokens"]
            counts["output_tokens"] += call["output_tokens"]
        if actual != calibration["models"]:
            raise ValueError(f"Trace usage disagrees with calibration: {gene}")
        result["calibrations"][gene] = {
            "source": calibration["source"],
            "source_sha256": calibration["source_sha256"],
            "attempts": calibration["attempts"],
            "calls": calls,
        }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        type=Path,
        default=Path(__file__).with_name("cost_calibration.json"),
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = export(args.profile, args.out)
    print(
        f"Exported {len(result['calibrations'])} usage receipts; SHA-256 {sha(args.out)}"
    )


if __name__ == "__main__":
    main()
