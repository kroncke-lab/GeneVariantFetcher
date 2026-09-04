#!/usr/bin/env python3
"""Combine already-opened candidate-tranche phenotype agreement figures.

This reads only locked runs already present in the mixed-gold consumption log.
It never opens an unscored answer key and never substitutes historical outputs
from a different protocol. Each source run must already have the per-run CSV
and JSON written by ``build_phenotype_count_recovery.py``.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
FIGURE_SCRIPT = REPO / "scripts" / "build_phenotype_count_recovery.py"


def load_figure_module():
    spec = importlib.util.spec_from_file_location("gvf_phenotype_figure", FIGURE_SCRIPT)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise SystemExit(f"cannot import {FIGURE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


figure = load_figure_module()


DEFAULT_LOG = (
    REPO / "benchmarks" / "evaluation_tiers" / "mixed_gold" / "consumption_log.jsonl"
)
DEFAULT_LOGS = (
    DEFAULT_LOG,
    REPO
    / "benchmarks"
    / "evaluation_tiers"
    / "mixed_gold_continuation_120"
    / "consumption_log.jsonl",
)
DEFAULT_RUNS_ROOT = REPO / "benchmarks" / "codex_paper_eval" / "runs"
DEFAULT_OUT = REPO / "docs" / "figures" / "mixed_gold_opened_candidates"

INTEGER_COLUMNS = {
    "manual_count",
    "ai_count_evaluation",
    "evaluation_zero_filled",
    "residual_ai_minus_manual",
    "absolute_difference",
    "exact_match",
    "zero_nonzero_disagreement",
}

# Hand-spaced examples from the two opened candidate tranches. They include
# exact multi-carrier recoveries, a small disagreement, and conspicuous nulls.
POOLED_ANNOTATIONS = {
    ("affected", "SCN5A", "30059973", "E1784K"): (-166.0, -324.0),
    ("affected", "MYBPC3", "33673806", "c.3181C>T"): (37.0, -134.0),
    ("affected", "RYR2", "15466642", "I4848V"): (80.0, -100.0),
    ("affected", "RYR2", "15466642", "I419F"): (-164.0, -231.0),
    ("unaffected", "KCNQ1", "21138517", "G589D"): (75.0, -125.0),
    ("unaffected", "KCNQ1", "33498651", "P631fsX"): (-155.0, -270.0),
    ("unaffected", "SCN5A", "26803770", "R1644H"): (82.0, -115.0),
    ("unaffected", "SCN5A", "11463728", "V1777M"): (-150.0, -195.0),
}

POOLED_DOCKS = {
    ("affected", "SCN5A", "30059973", "E1784K"): "right",
    ("affected", "RYR2", "15466642", "I419F"): "right",
    ("unaffected", "KCNQ1", "33498651", "P631fsX"): "right",
    ("unaffected", "SCN5A", "11463728", "V1777M"): "right",
}


def candidate_runs(
    consumption_log: Path, runs_root: Path, *, require: bool = True
) -> list[tuple[str, Path]]:
    """Return scored candidate runs in registry order without opening new gold."""
    found: list[tuple[str, Path]] = []
    for line in consumption_log.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        if record.get("comparison_arm") != "candidate":
            continue
        run_dir = runs_root / str(record["run_id"])
        found.append((str(record["tier_id"]), run_dir))
    if not found and require:
        raise SystemExit("no scored candidate arms found in the consumption log")
    return found


def candidate_runs_from_logs(
    consumption_logs: list[Path] | tuple[Path, ...], runs_root: Path
) -> list[tuple[str, Path]]:
    """Collect candidate runs from versioned suites in ledger order."""

    found: list[tuple[str, Path]] = []
    for path in consumption_logs:
        if path.is_file():
            found.extend(candidate_runs(path, runs_root, require=False))
    if not found:
        raise SystemExit("no scored candidate arms found in the consumption logs")
    if len({tier_id for tier_id, _ in found}) != len(found):
        raise SystemExit("duplicate candidate tier across consumption logs")
    return found


def read_rows(path: Path, tier_id: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for source in csv.DictReader(handle):
            row: dict[str, Any] = dict(source)
            for column in INTEGER_COLUMNS:
                row[column] = int(row[column])
            row["ai_count_raw"] = (
                int(row["ai_count_raw"]) if row["ai_count_raw"] != "" else None
            )
            row["analysis_tranche"] = tier_id
            rows.append(row)
    return rows


def build(run_entries: list[tuple[str, Path]], out_dir: Path) -> dict[str, Path]:
    all_rows: list[dict[str, Any]] = []
    source_summaries: list[dict[str, Any]] = []
    integrity: list[dict[str, Any]] = []
    for tier_id, run_dir in run_entries:
        lock = figure.verify_lock(run_dir)
        data_dir = run_dir / "figures" / "data"
        source_summary = json.loads(
            (data_dir / "phenotype_count_recovery.json").read_text(encoding="utf-8")
        )
        if source_summary.get("integrity") != lock["verified_hashes"]:
            raise SystemExit(f"figure integrity does not match lock for {run_dir.name}")
        all_rows.extend(read_rows(data_dir / "phenotype_count_recovery.csv", tier_id))
        source_summaries.append(source_summary)
        integrity.append(
            {
                "tier_id": tier_id,
                "run_id": run_dir.name,
                **lock["verified_hashes"],
            }
        )

    all_rows.sort(
        key=lambda row: (
            row["measure"],
            row["gene"],
            int(row["pmid"]),
            row["manual_variant"],
        )
    )
    keys = [
        (row["measure"], row["gene"], row["pmid"], row["manual_variant"])
        for row in all_rows
    ]
    if len(keys) != len(set(keys)):
        raise SystemExit("candidate tranches are not disjoint at plotted grain")

    missingness: dict[str, dict[str, int]] = {}
    report = {"overall": {"count": {}}}
    for field in figure.FIELDS:
        field_rows = [row for row in all_rows if row["measure"] == field]
        missingness[field] = {
            "matched": len(field_rows),
            "count_extracted": sum(
                row["ai_count_raw"] is not None for row in field_rows
            ),
            "zero_filled": sum(row["ai_count_raw"] is None for row in field_rows),
            "identity_miss": sum(
                int(item["metrics"][field]["identity_miss"])
                for item in source_summaries
            ),
        }
        report["overall"]["count"][field] = {
            "gold_asserted": sum(
                int(item["metrics"][field]["gold_asserted"])
                for item in source_summaries
            ),
            "predicted": missingness[field]["count_extracted"],
        }
    summary = figure.summarize(all_rows, missingness, report)

    paired_keys = {
        field: {
            (row["gene"], row["pmid"], row["manual_variant"])
            for row in all_rows
            if row["measure"] == field
        }
        for field in figure.FIELDS
    }
    paired_overlap = len(paired_keys["affected"] & paired_keys["unaffected"])
    tranche_numbers = [tier_id.rsplit("_", 1)[-1] for tier_id, _ in run_entries]
    title = (
        "Automated versus reference phenotype counts in opened mixed-gold candidates"
    )
    analysis_run = "mixed_gold_opened_candidates_" + "_".join(tranche_numbers)
    locked_at = max(str(item["locked_at"]) for item in source_summaries)

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "svg": out_dir / "phenotype_count_recovery.svg",
        "png": out_dir / "phenotype_count_recovery.png",
        "pdf": out_dir / "phenotype_count_recovery.pdf",
        "csv": out_dir / "phenotype_count_recovery.csv",
        "json": out_dir / "phenotype_count_recovery.json",
    }
    figure.write_csv(all_rows, paths["csv"])
    paths["svg"].write_text(
        figure.render_zero_filled_svg(
            all_rows,
            summary,
            paired_overlap,
            analysis_run,
            locked_at,
            title=title,
            example_annotations=POOLED_ANNOTATIONS,
            example_docks=POOLED_DOCKS,
        ),
        encoding="utf-8",
    )
    paths["json"].write_text(
        json.dumps(
            {
                "analysis_run": analysis_run,
                "source_runs": integrity,
                "locked_at": locked_at,
                "row_grain": "identity-matched variant-paper row and count measure",
                "paired_rows_in_both_measures": paired_overlap,
                "metrics": summary,
                "figure_contract": {
                    "title": title,
                    "evaluation_missing_count": 0,
                    "storage_missing_count": None,
                    "included_arms": "scored candidate arms only",
                    "included_tranches": [tier_id for tier_id, _ in run_entries],
                    "unopened_tranches": "excluded",
                    "marker_size_key_pairs": list(figure.MARKER_SIZE_KEY_COUNTS),
                    "labeled_variant_paper_pairs": [
                        {
                            "measure": measure,
                            "gene": gene,
                            "pmid": pmid,
                            "variant": variant,
                        }
                        for measure, gene, pmid, variant in POOLED_ANNOTATIONS
                        if any(
                            row["measure"] == measure
                            and row["gene"] == gene
                            and row["pmid"] == pmid
                            and row["manual_variant"] == variant
                            for row in all_rows
                        )
                    ],
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    figure.export_svg_copies(paths["svg"], paths["png"], paths["pdf"])
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--consumption-log",
        type=Path,
        action="append",
        help="Candidate-arm burn ledger; repeat to combine versioned suites",
    )
    parser.add_argument("--runs-root", type=Path, default=DEFAULT_RUNS_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--run-dir",
        type=Path,
        action="append",
        help="Explicit locked candidate run; repeat to override log discovery",
    )
    args = parser.parse_args()
    if args.run_dir:
        entries = [
            (f"candidate_{index:02d}", path)
            for index, path in enumerate(args.run_dir, 1)
        ]
    else:
        entries = candidate_runs_from_logs(
            args.consumption_log or DEFAULT_LOGS, args.runs_root
        )
    paths = build(entries, args.out_dir)
    for label, path in paths.items():
        print(f"wrote {label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
