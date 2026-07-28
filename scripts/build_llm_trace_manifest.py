#!/usr/bin/env python3
"""Build or refresh the SHA-256 manifest for an LLM trace directory."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils.llm_trace import TRACE_MANIFEST_NAME, build_trace_manifest  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("trace_dir", type=Path)
    parser.add_argument("--run-id")
    args = parser.parse_args()
    output = args.trace_dir / TRACE_MANIFEST_NAME
    manifest = build_trace_manifest(
        args.trace_dir,
        output_path=output,
        run_id=args.run_id,
    )
    print(
        json.dumps(
            {
                "manifest": str(output),
                "llm_calls": manifest["llm_call_count"],
                "decision_events": manifest["decision_event_count"],
                "records": manifest["trace_count"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
