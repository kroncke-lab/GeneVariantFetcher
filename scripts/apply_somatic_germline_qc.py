#!/usr/bin/env python3
"""Measure (and optionally gate) somatic contamination in a run's extractions.

Gold-free, no-network. Runs ``pipeline.somatic_germline_qc`` over every
``*.json`` in an extractions dir and reports how many extracted variant records
look germline vs somatic vs ambiguous vs unknown — the contamination estimate a
cold-start CANCER-gene run (e.g. BRCA2) needs before its carrier counts can be
trusted. Default is report-only (dry run); ``--write`` persists the per-variant
``somatic_germline_flag`` annotations, and ``--policy drop_somatic`` additionally
removes somatic-labelled records (ambiguous/unknown are always kept).

Examples:
    # Measure contamination in the BRCA2 pilot (no mutation):
    python scripts/apply_somatic_germline_qc.py \
        --extractions results/BRCA2/<run>/extractions --report-out /tmp/brca2_qc.json

    # Persist flags (still non-destructive — annotation only):
    python scripts/apply_somatic_germline_qc.py --extractions <dir> --write
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pipeline.somatic_germline_qc import (  # noqa: E402
    AMBIGUOUS,
    GERMLINE,
    SOMATIC,
    UNKNOWN,
    annotate_variants,
)


def _paper_context(payload: dict) -> str:
    """Title/abstract give strong paper-level somatic vs germline signal."""
    return "\n".join(
        str(payload.get(k, "")) for k in ("title", "abstract") if payload.get(k)
    )


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--extractions", required=True, type=Path)
    ap.add_argument(
        "--policy",
        choices=("flag", "drop_somatic"),
        default="flag",
        help="flag (default): annotate only. drop_somatic: also remove somatic "
        "records (implies --write).",
    )
    ap.add_argument(
        "--write",
        action="store_true",
        help="Persist annotations (and drops) back into the JSON files.",
    )
    ap.add_argument("--report-out", type=Path, default=None)
    return ap


def main() -> int:
    args = build_parser().parse_args()
    extractions = args.extractions.expanduser().resolve()
    if not extractions.is_dir():
        print(f"error: not a directory: {extractions}", file=sys.stderr)
        return 2

    write = args.write or args.policy == "drop_somatic"
    totals = {GERMLINE: 0, SOMATIC: 0, AMBIGUOUS: 0, UNKNOWN: 0}
    per_pmid: list[dict] = []
    files = 0
    dropped_total = 0

    for json_path in sorted(extractions.glob("*.json")):
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        variants = payload.get("variants")
        if not isinstance(variants, list) or not variants:
            continue
        files += 1
        summary = annotate_variants(
            variants, paper_context=_paper_context(payload), policy=args.policy
        )
        for label, n in summary["counts"].items():
            totals[label] += n
        dropped_total += summary["dropped_somatic"]
        if summary["somatic_fraction"] > 0:
            per_pmid.append(
                {
                    "file": json_path.name,
                    "pmid": payload.get("pmid"),
                    "somatic_fraction": round(summary["somatic_fraction"], 3),
                    "counts": summary["counts"],
                }
            )
        if write:
            json_path.write_text(
                json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
            )

    total = sum(totals.values())
    contaminated = totals[SOMATIC] + totals[AMBIGUOUS]
    report = {
        "extractions": str(extractions),
        "policy": args.policy,
        "written": write,
        "files_with_variants": files,
        "total_variants": total,
        "counts": totals,
        "somatic_fraction": round((contaminated / total) if total else 0.0, 4),
        "dropped_somatic": dropped_total,
        "top_contaminated_pmids": sorted(
            per_pmid, key=lambda r: r["somatic_fraction"], reverse=True
        )[:25],
    }

    print(
        json.dumps(
            {k: v for k, v in report.items() if k != "top_contaminated_pmids"}, indent=2
        )
    )
    print(
        f"\nSomatic-contamination estimate: {contaminated}/{total} "
        f"({report['somatic_fraction']:.1%}) of extracted variant records are "
        f"somatic+ambiguous (germline={totals[GERMLINE]}, unknown={totals[UNKNOWN]})."
    )
    if args.report_out:
        args.report_out.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"Report written: {args.report_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
