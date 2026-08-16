#!/usr/bin/env python3
"""Exact-match precision of emitted carrier/affected/unaffected integers.

This is the third gold-120 denominator described in
``docs/AFFECTED_UNAFFECTED_PRECISION.md``: on a matched gold row, when *both*
gold and the prediction supply a value for a field,

    exact-match precision = exact / supplied

It is not Gate 2 counted-extra identity precision and not carrier MAE. A null
prediction against a gold ``0`` is not an error here (that convention gap is
count recall).

The script rescores an existing prediction set with the current
``run_eval.py`` matcher and the live gold snapshot. It performs no extraction
and makes no LLM calls, so a candidate guard rule can be measured for free
before anyone spends a paid re-extract arm.

Examples::

    scripts/phenotype_value_precision.py --run-dir <run> --guard none
    scripts/phenotype_value_precision.py --run-dir <run> --guard current \
        --errors-out /tmp/errors.tsv
    scripts/phenotype_value_precision.py --run-dir <run> --compare

``--compare`` prints the baseline and guarded tables side by side, which is the
form quoted in the briefing doc.
"""

from __future__ import annotations

import argparse
import copy
import csv
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Optional

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

RUN_EVAL_PATH = REPO / "benchmarks" / "codex_paper_eval" / "run_eval.py"
DEFAULT_GOLD = REPO / "gene_variant_fetcher_gold_standard" / "normalized"

FIELDS = ("carriers", "affected", "unaffected")


def load_run_eval():
    """Import the benchmark harness by path (it is not an installed module)."""
    spec = importlib.util.spec_from_file_location("gvf_run_eval", RUN_EVAL_PATH)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise SystemExit(f"cannot import run_eval from {RUN_EVAL_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def apply_guard(papers: list[dict[str, Any]], guard: str) -> dict[str, int]:
    """Apply a phenotype guard in memory. Returns per-reason clear counts."""
    if guard == "none":
        return {}
    if guard != "current":
        raise SystemExit(f"unknown guard: {guard}")

    from pipeline.phenotype_count_guard import apply_phenotype_count_guard

    reasons: dict[str, int] = {}
    for paper in papers:
        summary = apply_phenotype_count_guard(paper.get("variants") or [])
        for note in summary.annotations:
            key = f"{note['field']}:{note['reason']}"
            reasons[key] = reasons.get(key, 0) + 1
    return reasons


def score(run_eval, papers: list[dict[str, Any]], selection: dict, gold_root: Path):
    pred_map = {(p["gene"], str(p["pmid"])): p for p in papers}
    scores = []
    for paper in selection["papers"]:
        key = (paper["gene"], paper["pmid"])
        predicted = pred_map.get(key)
        if predicted is None:
            continue
        gold = run_eval.load_gold(gold_root, *key)
        scores.append(run_eval.score_one(*key, predicted, gold))
    return scores


def field_table(scores: list[dict]) -> dict[str, dict[str, Any]]:
    table: dict[str, dict[str, Any]] = {}
    for field in FIELDS:
        supplied = sum(paper["count"][field]["predicted"] for paper in scores)
        errors = [
            error
            for paper in scores
            for error in paper["count_errors"]
            if error["field"] == field
        ]
        exact = supplied - len(errors)
        table[field] = {
            "supplied": supplied,
            "exact": exact,
            "errors": len(errors),
            "precision": exact / supplied if supplied else None,
            # MAE is over every supplied pair, not only the erroring ones.
            "mae": (
                sum(abs(error["error"]) for error in errors) / supplied
                if supplied
                else None
            ),
        }
    return table


def identity_totals(scores: list[dict]) -> dict[str, Any]:
    tp = sum(paper["tp"] for paper in scores)
    fp = sum(paper["fp"] for paper in scores)
    fn = sum(paper["fn"] for paper in scores)
    counted_extra = sum(
        paper["counted_precision"]["counted_extra_rows"] for paper in scores
    )
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "counted_extra_rows": counted_extra,
        "raw_precision": tp / (tp + fp) if tp + fp else None,
        "counted_extra_precision": (
            tp / (tp + counted_extra) if tp + counted_extra else None
        ),
        "recall": tp / (tp + fn) if tp + fn else None,
    }


def _pct(value: Optional[float]) -> str:
    return "n/a" if value is None else f"{value * 100:.2f}%"


def _num(value: Optional[float]) -> str:
    return "n/a" if value is None else f"{value:.4f}"


def print_report(label: str, table: dict, identity: dict) -> None:
    print(f"\n== {label} ==")
    header = f"{'field':<12}{'exact/supplied':>16}{'precision':>12}"
    print(header + f"{'errors':>9}{'MAE':>10}")
    for field in FIELDS:
        row = table[field]
        pair = f"{row['exact']}/{row['supplied']}"
        print(
            f"{field:<12}{pair:>16}{_pct(row['precision']):>12}"
            f"{row['errors']:>9}{_num(row['mae']):>10}"
        )
    print(
        f"identity: TP {identity['tp']} FP {identity['fp']} FN {identity['fn']} "
        f"| recall {_pct(identity['recall'])} "
        f"| raw precision {_pct(identity['raw_precision'])} "
        f"| counted-extra precision {_pct(identity['counted_extra_precision'])} "
        f"({identity['counted_extra_rows']} counted extras)"
    )


def error_rows(scores: list[dict], papers: list[dict]) -> list[dict[str, Any]]:
    """Join each value error back to the prediction row that produced it."""
    by_paper = {(p["gene"], str(p["pmid"])): p for p in papers}
    rows: list[dict[str, Any]] = []
    for paper in scores:
        source = by_paper.get((paper["gene"], str(paper["pmid"]))) or {}
        variants = {
            str(v.get("variant")): v for v in (source.get("variants") or []) if v
        }
        for error in paper["count_errors"]:
            variant = variants.get(str(error["variant"]), {})
            rows.append(
                {
                    "gene": paper["gene"],
                    "pmid": paper["pmid"],
                    "variant": error["variant"],
                    "field": error["field"],
                    "predicted": error["predicted"],
                    "gold": error["gold"],
                    "error": error["error"],
                    "source_layer": variant.get("source_layer") or "",
                    "source_location": variant.get("source_location") or "",
                    "carriers": variant.get("carriers"),
                    "affected": variant.get("affected"),
                    "unaffected": variant.get("unaffected"),
                    "evidence": " ".join(str(variant.get("evidence") or "").split())[
                        :400
                    ],
                }
            )
    rows.sort(key=lambda row: (-abs(row["error"]), row["gene"], str(row["pmid"])))
    return rows


def write_errors(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {len(rows)} value errors to {path}")


def run(args) -> dict[str, Any]:
    run_eval = load_run_eval()
    selection = json.loads((args.run_dir / "selection.json").read_text())
    predictions = json.loads((args.run_dir / "predictions.json").read_text())

    results: dict[str, Any] = {"run_dir": str(args.run_dir), "arms": {}}
    guards = ["none", "current"] if args.compare else [args.guard]
    scored: dict[str, Any] = {}
    for guard in guards:
        papers = copy.deepcopy(predictions["papers"])
        reasons = apply_guard(papers, guard)
        scores = score(run_eval, papers, selection, args.gold_root)
        table = field_table(scores)
        identity = identity_totals(scores)
        label = "baseline (no guard)" if guard == "none" else f"guard={guard}"
        print_report(label, table, identity)
        if reasons:
            joined = ", ".join(f"{k}={v}" for k, v in sorted(reasons.items()))
            print("cleared: " + joined)
        results["arms"][guard] = {
            "fields": table,
            "identity": identity,
            "cleared": reasons,
        }
        scored[guard] = (scores, papers)

    if args.errors_out:
        guard = guards[-1]
        scores, papers = scored[guard]
        rows = error_rows(scores, papers)
        if rows:
            write_errors(rows, args.errors_out)

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(results, indent=2, sort_keys=True))
        print(f"wrote {args.json_out}")
    return results


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--gold-root", type=Path, default=DEFAULT_GOLD)
    parser.add_argument("--guard", default="current", choices=["none", "current"])
    parser.add_argument(
        "--compare",
        action="store_true",
        help="score both the unguarded and guarded arms",
    )
    parser.add_argument("--errors-out", type=Path)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
