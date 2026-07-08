#!/usr/bin/env python3
"""Apply the count-classifier (#5) to existing extraction JSONs.

Walks an extractions directory, reads each variant's ``count_provenance``
(populated by the LLM per the prompts in pipeline/prompts.py), and either
flags or clears count assignments whose declared count_type is not
``per_variant_carrier``.

Sibling of ``scripts/apply_count_outlier_guard.py``. The two are
complementary: the outlier guard fires on statistical anomalies regardless
of provenance; the classifier fires on declared provenance regardless of
value magnitude. Running both ("clear then clear") is safe — operations
are idempotent on a null count.

Typical workflow on existing runs:

    .venv/bin/python scripts/apply_count_classifier.py \\
        --extractions <run>/extractions \\
        --policy clear \\
        --backup-dir <run>/count_classifier_backup_<date> \\
        --report-out /tmp/classifier_report.json

Then rebuild the DB with harvesting.migrate_to_sqlite and re-score with
scripts.run_recall_suite to measure MAE delta.

The classifier silently skips variants without a ``count_provenance``
block, so it is safe to run on older extraction JSONs from before #4
shipped — they simply pass through unchanged.
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

from pipeline.count_classifier import (  # noqa: E402
    detect_misclassified_counts,
    enforce_per_variant_policy,
)

COUNT_CLASSIFIER_FIELDS = {"carriers", "affected", "unaffected"}


def parse_fields(value: str) -> list[str] | None:
    """Return requested classifier fields; None means all fields."""
    raw = str(value or "").strip().lower()
    if raw in {"", "all", "*"}:
        return None
    fields: list[str] = []
    for item in raw.split(","):
        field = item.strip().lower()
        if not field:
            continue
        if field not in COUNT_CLASSIFIER_FIELDS:
            raise argparse.ArgumentTypeError(
                f"--fields entries must be one of "
                f"{sorted(COUNT_CLASSIFIER_FIELDS)} or 'all'; got {field!r}"
            )
        if field not in fields:
            fields.append(field)
    return fields or None


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
            "flag (default): annotate variants with count_classifier_flags. "
            "clear: also null the offending count fields (raw preserved)."
        ),
    )
    p.add_argument(
        "--backup-dir",
        type=Path,
        default=None,
        help="Required when --policy is flag or clear: directory to copy originals into.",
    )
    p.add_argument(
        "--report-out",
        type=Path,
        required=True,
        help="Where to write the per-PMID classifier report (JSON).",
    )
    p.add_argument(
        "--fields",
        type=parse_fields,
        default=None,
        help=(
            "Comma-separated logical fields to classify: all, carriers, affected, "
            "unaffected. Default: all."
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

    policy = "off" if args.dry_run else args.policy
    if policy in {"flag", "clear"} and args.backup_dir is None:
        sys.exit("--backup-dir is required when --policy is flag or clear")
    backup_dir = (
        args.backup_dir.expanduser().resolve() if args.backup_dir is not None else None
    )
    if backup_dir is not None:
        backup_dir.mkdir(parents=True, exist_ok=True)
    classifier_fields = args.fields

    report_rows: list[dict[str, object]] = []
    totals = {
        "pmids_processed": 0,
        "pmids_with_misclassified": 0,
        "pmids_with_provenance": 0,
        "total_flagged_values": 0,
        "total_cleared_values": 0,
    }

    for json_path in sorted(extractions.glob("*.json")):
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
        except Exception as exc:  # noqa: BLE001
            report_rows.append(
                {"file": json_path.name, "error": f"parse_failure: {exc}"}
            )
            continue
        if not isinstance(payload, dict):
            continue
        variants = payload.get("variants")
        if not isinstance(variants, list):
            continue

        totals["pmids_processed"] += 1
        if any(
            isinstance(v, dict) and isinstance(v.get("count_provenance"), dict)
            for v in variants
        ):
            totals["pmids_with_provenance"] += 1

        annotations = detect_misclassified_counts(variants, fields=classifier_fields)
        if not annotations:
            continue

        if policy in {"flag", "clear"} and backup_dir is not None:
            shutil.copy2(json_path, backup_dir / json_path.name)

        result = enforce_per_variant_policy(variants, annotations, policy=policy)
        if policy in {"flag", "clear"}:
            payload["count_classifier"] = {
                "policy": policy,
                "fields": classifier_fields or "all",
                "flagged": result.flagged,
                "cleared": result.cleared,
                "annotations": result.as_dict()["annotations"],
            }
            json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        totals["pmids_with_misclassified"] += 1
        totals["total_flagged_values"] += result.flagged
        totals["total_cleared_values"] += result.cleared
        report_rows.append(
            {"file": json_path.name, "policy": policy, **result.as_dict()}
        )

    summary = {
        "extractions_dir": str(extractions),
        "policy": policy,
        "fields": classifier_fields or "all",
        "totals": totals,
        "rows": report_rows,
    }
    args.report_out.parent.mkdir(parents=True, exist_ok=True)
    args.report_out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(totals, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
