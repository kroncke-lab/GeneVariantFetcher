#!/usr/bin/env python3
"""Apply an explicit source-backed diagnostic overlay to locked predictions.

This utility never reads gold. It exists to make post-lock source diagnostics
reproducible without mutating the immutable baseline or production databases.
The resulting predictions must be locked in a separate run directory before
the ordinary scorer can read gold.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def apply_overlay(predictions: dict[str, Any], overlay: dict[str, Any]) -> None:
    papers = {
        (str(paper.get("gene")), str(paper.get("pmid"))): paper
        for paper in predictions.get("papers", [])
    }
    for operation in overlay.get("operations", []):
        key = (str(operation.get("gene")), str(operation.get("pmid")))
        if key not in papers:
            raise ValueError(f"overlay paper is absent from predictions: {key}")
        variants = papers[key].setdefault("variants", [])
        variant_text = str(operation.get("variant") or "")
        matches = [row for row in variants if row.get("variant") == variant_text]
        kind = operation.get("operation")
        if kind == "set_counts":
            if len(matches) != 1:
                raise ValueError(
                    f"set_counts requires one exact variant match: {key} "
                    f"{variant_text!r}; found {len(matches)}"
                )
            for field in ("carriers", "affected", "unaffected"):
                if field in operation:
                    matches[0][field] = operation[field]
            matches[0]["diagnostic_source"] = operation.get("source")
        elif kind == "add_variant":
            if matches:
                raise ValueError(
                    f"add_variant would duplicate an exact variant: {key} "
                    f"{variant_text!r}"
                )
            variants.append(
                {
                    "variant": variant_text,
                    "carriers": operation.get("carriers"),
                    "affected": operation.get("affected"),
                    "unaffected": operation.get("unaffected"),
                    "evidence": operation.get("evidence"),
                    "source_location": operation.get("source_location"),
                    "source_layer": "source_backed_diagnostic",
                    "diagnostic_source": operation.get("source"),
                }
            )
        else:
            raise ValueError(f"unsupported overlay operation: {kind!r}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--overlay", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    overlay = load_json(args.overlay)
    expected = str(overlay.get("baseline_predictions_sha256") or "")
    actual = digest(args.predictions)
    if not expected or actual != expected:
        raise SystemExit(
            "baseline prediction digest mismatch: "
            f"expected={expected or '<missing>'} actual={actual}"
        )
    predictions = load_json(args.predictions)
    apply_overlay(predictions, overlay)
    predictions["diagnostic_overlay"] = {
        "overlay": str(args.overlay),
        "overlay_sha256": digest(args.overlay),
        "baseline_predictions_sha256": actual,
        "gold_read": False,
    }
    args.out.write_text(
        json.dumps(predictions, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
