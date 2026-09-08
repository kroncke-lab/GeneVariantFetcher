#!/usr/bin/env python3
"""Replay the table phenotype coverage audit without extraction or API calls."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pipeline.table_phenotype_coverage import audit_table_phenotype_coverage  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extraction", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve() in {args.source.resolve(), args.extraction.resolve()}:
        parser.error("output must not overwrite a source or extraction")
    extraction_bytes = args.extraction.read_bytes()
    source_bytes = args.source.read_bytes()
    report = {
        "extraction": str(args.extraction.resolve()),
        "extraction_sha256": hashlib.sha256(extraction_bytes).hexdigest(),
        "source": str(args.source.resolve()),
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "api_calls": 0,
        "coverage": audit_table_phenotype_coverage(
            json.loads(extraction_bytes), source_bytes.decode("utf-8")
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["coverage"], indent=2))


if __name__ == "__main__":
    main()
