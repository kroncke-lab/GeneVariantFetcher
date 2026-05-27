#!/usr/bin/env python3
"""Apply the count-outlier guard to existing extraction JSONs.

Walks an extractions directory, runs `pipeline.count_outlier_guard` on each
PMID, and (in "flag" mode) annotates per-variant `count_outlier_flags`, or
(in "clear" mode) zeros the flagged count fields while preserving the raw
value in `count_outlier_flags`. Writes a per-run summary report.

This is the retroactive entry point used to lift rows-mode MAE on existing
runs without re-extracting. The detector is gold-free; this script does not
read or write gold standards.

Typical usage:

    .venv/bin/python scripts/apply_count_outlier_guard.py \\
        --extractions results/KCNH2/20260506_102238/extractions \\
        --policy flag \\
        --report-out /tmp/kcnh2_outlier_report.json

When you are ready to actually drop suspected study-wide-N values:

    .venv/bin/python scripts/apply_count_outlier_guard.py \\
        --extractions results/KCNH2/20260506_102238/extractions \\
        --policy clear \\
        --backup-dir results/KCNH2/20260506_102238/outlier_guard_backup_20260527 \\
        --report-out /tmp/kcnh2_outlier_report.json

After --policy clear, rebuild the DB with harvesting.migrate_to_sqlite and
re-score with scripts.run_recall_suite to measure MAE delta.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pipeline.count_outlier_guard import (  # noqa: E402
    DEFAULT_ABSOLUTE_THRESHOLD,
    DEFAULT_MIN_VARIANTS,
    DEFAULT_MULTIPLIER_THRESHOLD,
    apply_outlier_policy,
    detect_count_outliers,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--extractions",
        required=True,
        type=Path,
        help="Path to extractions directory containing <GENE>_PMID_<PMID>.json files.",
    )
    p.add_argument(
        "--policy",
        choices=["off", "flag", "clear"],
        default="flag",
        help=(
            "off: detect and report only; do not modify JSONs. "
            "flag (default): annotate variants with count_outlier_flags metadata. "
            "clear: also zero the flagged count fields (preserves raw under flags)."
        ),
    )
    p.add_argument(
        "--backup-dir",
        type=Path,
        default=None,
        help=(
            "Required when --policy is flag or clear: directory to copy "
            "originals into before modifying."
        ),
    )
    p.add_argument(
        "--report-out",
        type=Path,
        required=True,
        help="Where to write the per-PMID outlier report (JSON).",
    )
    p.add_argument(
        "--min-variants",
        type=int,
        default=DEFAULT_MIN_VARIANTS,
        help=(
            "Minimum variants per PMID before the detector runs. "
            "Defaults to {}.".format(DEFAULT_MIN_VARIANTS)
        ),
    )
    p.add_argument(
        "--multiplier-threshold",
        type=float,
        default=DEFAULT_MULTIPLIER_THRESHOLD,
        help=(
            "Flag values that are more than this multiple of the per-paper "
            "median. Defaults to {}.".format(DEFAULT_MULTIPLIER_THRESHOLD)
        ),
    )
    p.add_argument(
        "--absolute-threshold",
        type=int,
        default=DEFAULT_ABSOLUTE_THRESHOLD,
        help=(
            "Flag values must also exceed this absolute threshold. "
            "Defaults to {}.".format(DEFAULT_ABSOLUTE_THRESHOLD)
        ),
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Detect and report only; never write changes (overrides --policy).",
    )
    return p


def main() -> int:
    args = build_parser().parse_args()
    extractions = args.extractions.expanduser().resolve()
    if not extractions.is_dir():
        sys.exit(f"extractions dir not found: {extractions}")

    policy = args.policy
    if args.dry_run:
        policy = "off"

    if policy in {"flag", "clear"} and args.backup_dir is None:
        sys.exit("--backup-dir is required when --policy is flag or clear")
    backup_dir = (
        args.backup_dir.expanduser().resolve() if args.backup_dir is not None else None
    )
    if backup_dir is not None:
        backup_dir.mkdir(parents=True, exist_ok=True)

    report_rows: list[dict[str, object]] = []
    totals = {
        "pmids_processed": 0,
        "pmids_with_outliers": 0,
        "total_flagged_values": 0,
        "total_cleared_values": 0,
    }

    for json_path in sorted(extractions.glob("*.json")):
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
        except Exception as exc:  # noqa: BLE001
            report_rows.append(
                {
                    "file": json_path.name,
                    "error": f"parse_failure: {exc}",
                }
            )
            continue
        if not isinstance(payload, dict):
            continue
        variants = payload.get("variants")
        if not isinstance(variants, list):
            continue

        totals["pmids_processed"] += 1
        annotations = detect_count_outliers(
            variants,
            min_variants=args.min_variants,
            multiplier_threshold=args.multiplier_threshold,
            absolute_threshold=args.absolute_threshold,
        )
        if not annotations:
            continue

        if policy in {"flag", "clear"} and backup_dir is not None:
            shutil.copy2(json_path, backup_dir / json_path.name)

        result = apply_outlier_policy(variants, annotations, policy=policy)

        if policy in {"flag", "clear"}:
            payload["count_outlier_guard"] = {
                "policy": policy,
                "flagged": result.flagged,
                "cleared": result.cleared,
                "annotations": result.as_dict()["annotations"],
            }
            json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        totals["pmids_with_outliers"] += 1
        totals["total_flagged_values"] += result.flagged
        totals["total_cleared_values"] += result.cleared
        report_rows.append(
            {
                "file": json_path.name,
                "policy": policy,
                **result.as_dict(),
            }
        )

    summary = {
        "extractions_dir": str(extractions),
        "policy": policy,
        "thresholds": {
            "min_variants": args.min_variants,
            "multiplier_threshold": args.multiplier_threshold,
            "absolute_threshold": args.absolute_threshold,
        },
        "totals": totals,
        "rows": report_rows,
    }
    args.report_out.parent.mkdir(parents=True, exist_ok=True)
    args.report_out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(totals, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
