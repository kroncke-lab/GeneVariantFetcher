#!/usr/bin/env python3
"""Zero-LLM replay of the table-cohort phenotype projection on a locked run.

The projection in ``pipeline/table_cohort_phenotype.py`` is a deterministic
post-processor over already-extracted table rows, so its effect can be
measured on a finished evaluation run without spending anything: re-run the
module over the archived per-PMID extraction JSON and the run's own frozen
source text, re-apply the always-on phenotype guard, fill the resulting values
into a copy of the locked ``predictions.json`` (fill-null-only, never
overwriting a stored value), and rescore both arms with the harness matcher.

Two modes:

* ``derive`` (default) -- the archived extraction predates the module; the
  ``on`` arm is the locked predictions plus the derived values.
* ``strip`` -- the archived extraction already carries the module's stamps; the
  ``off`` arm nulls every field whose provenance source is the module's stamp.

Either way the two arms differ only by this step, which removes provider
run-to-run variance from the comparison. The script never writes into the
locked run directory.

Example::

    scripts/replay_table_cohort_phenotype.py \
        --run-dir benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_02_candidate \
        --run-dir benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_03_candidate \
        --out-dir docs/evidence/table_cohort_phenotype_20260907/replay
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
os.environ.setdefault("GVF_DISABLE_LOCAL_DATA", "1")

from benchmarks.codex_paper_eval.production_run import (  # noqa: E402
    ProductionRunError,
    resolve_active_gene_run,
)
from pipeline.count_provenance import TABLE_COHORT_PHENOTYPE_SOURCE  # noqa: E402
from pipeline.phenotype_count_guard import apply_phenotype_count_guard  # noqa: E402
from pipeline.table_cohort_phenotype import (  # noqa: E402
    METHOD,
    derive_table_cohort_phenotype_counts,
)

RUN_EVAL_PATH = REPO / "benchmarks" / "codex_paper_eval" / "run_eval.py"
FIELDS = ("carriers", "affected", "unaffected")
CARDIAC = ("KCNH2", "KCNQ1", "SCN5A", "RYR2")


def load_run_eval():
    spec = importlib.util.spec_from_file_location("gvf_run_eval", RUN_EVAL_PATH)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise SystemExit(f"cannot import run_eval from {RUN_EVAL_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _norm_token(value: Any) -> str:
    text = str(value or "").strip().lower()
    text = re.sub(r"^(?:p\.|c\.)", "", text)
    return re.sub(r"[\s()\[\]*]", "", text)


def record_tokens(record: dict[str, Any]) -> set[str]:
    tokens = set()
    for key in (
        "protein_notation",
        "cdna_notation",
        "legacy_notation",
        "source_notation",
    ):
        token = _norm_token(record.get(key))
        if token:
            tokens.add(token)
    return tokens


def prediction_tokens(row: dict[str, Any]) -> set[str]:
    tokens = set()
    for text in [row.get("variant"), *(row.get("merged_notations") or [])]:
        for part in str(text or "").split():
            token = _norm_token(part)
            if token:
                tokens.add(token)
    return tokens


def _table_key(text: Any) -> str:
    text = re.sub(r",\s*row\s+\d+.*$", "", str(text or ""), flags=re.IGNORECASE)
    text = re.sub(r"\s*\((?:regex|router)[^)]*\)\s*$", "", text)
    return re.sub(r"[^a-z0-9]+", "", text.lower())[:40]


def source_text_for(extraction: dict[str, Any], gene_run: Path, pmid: str) -> str:
    metadata = extraction.get("extraction_metadata") or {}
    candidates = [metadata.get("source_file")]
    for name in (f"{pmid}_FULL_CONTEXT.md", f"{pmid}_CLEANED.md"):
        candidates.append(gene_run / "pmc_fulltext" / name)
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate)
        if path.is_file():
            return path.read_text(errors="replace")
    return ""


def disease_for(gene_run: Path) -> Optional[str]:
    path = gene_run / "gene_disease_context.json"
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text())
    except json.JSONDecodeError:
        return None
    for key in ("disease", "disease_name", "phenotype"):
        value = payload.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return None


def derived_records(
    extraction: dict[str, Any],
    source_text: str,
    *,
    gene: str,
    disease: Optional[str],
    mode: str,
    paper_tier: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Return per-record target values plus the module's metadata block."""
    title = (extraction.get("paper_metadata") or {}).get("title")
    if mode == "derive":
        data = derive_table_cohort_phenotype_counts(
            extraction,
            source_text,
            gene_symbol=gene,
            disease=disease,
            title=title,
            enabled=True,
            allow_paper_tier=paper_tier,
        )
        # Mimic the persist-site guard so only values the guard keeps count.
        apply_phenotype_count_guard(data.get("variants") or [])
        metadata = (data.get("extraction_metadata") or {}).get(
            "table_cohort_phenotype_derivation", {}
        )
    else:
        data = copy.deepcopy(extraction)
        metadata = (data.get("extraction_metadata") or {}).get(
            "table_cohort_phenotype_derivation", {}
        )
    records: list[dict[str, Any]] = []
    for variant in data.get("variants") or []:
        if not isinstance(variant, dict):
            continue
        provenance = variant.get("count_provenance") or {}
        stamped = {
            field
            for field in ("affected", "unaffected")
            if str(provenance.get(f"{field}_source") or "").strip().lower()
            == TABLE_COHORT_PHENOTYPE_SOURCE
        }
        derivation = variant.get("phenotype_derivation") or {}
        if derivation.get("method") != METHOD and not stamped:
            continue
        penetrance = variant.get("penetrance_data") or {}
        targets: dict[str, Optional[int]] = {}
        for field in stamped:
            key = "affected_count" if field == "affected" else "unaffected_count"
            targets[field] = None if mode == "strip" else penetrance.get(key)
        if not targets:
            continue
        records.append(
            {
                "tokens": record_tokens(variant),
                "table_key": _table_key(
                    variant.get("source_location")
                    or variant.get("source_table")
                    or (variant.get("patients") or {}).get("source_ref")
                ),
                "targets": targets,
                "role": derivation.get("cohort_role"),
                "tier": derivation.get("tier"),
                "caption": derivation.get("source_table"),
                "count_column": derivation.get("count_column"),
                "identity": variant.get("protein_notation")
                or variant.get("cdna_notation")
                or variant.get("source_notation"),
            }
        )
    return records, metadata


def patch_paper(
    paper: dict[str, Any], records: list[dict[str, Any]], mode: str
) -> list[dict[str, Any]]:
    """Fill (or null) prediction rows from derived records; return the log."""
    log: list[dict[str, Any]] = []
    for row in paper.get("variants") or []:
        tokens = prediction_tokens(row)
        if not tokens:
            continue
        matches = [r for r in records if r["tokens"] & tokens]
        if not matches:
            continue
        if len(matches) > 1:
            row_key = _table_key(row.get("source_location"))
            located = [
                r
                for r in matches
                if r["table_key"]
                and row_key
                and (
                    r["table_key"].startswith(row_key)
                    or row_key.startswith(r["table_key"])
                )
            ]
            if len(located) == 1:
                matches = located
            else:
                values = {json.dumps(r["targets"], sort_keys=True) for r in matches}
                if len(values) != 1:
                    log.append(
                        {
                            "variant": row.get("variant"),
                            "status": "ambiguous_records",
                            "candidates": len(matches),
                        }
                    )
                    continue
        record = matches[0]
        for field, value in record["targets"].items():
            before = row.get(field)
            if mode == "derive":
                if before is not None or value is None:
                    continue
                row[field] = value
            else:
                if before is None:
                    continue
                row[field] = None
            log.append(
                {
                    "variant": row.get("variant"),
                    "field": field,
                    "before": before,
                    "after": row[field],
                    "role": record["role"],
                    "tier": record["tier"],
                    "caption": record["caption"],
                    "count_column": record["count_column"],
                    "source_location": row.get("source_location"),
                }
            )
    return log


def gold_for_variant(
    run_eval, gene: str, variant: str, gold: list[dict]
) -> Optional[dict]:
    for row in gold:
        if run_eval.matches(variant, row["variant"], gene):
            return row
    return None


def summarize_arm(
    scores: list[dict], genes: Optional[set[str]] = None
) -> dict[str, Any]:
    rows = [s for s in scores if genes is None or s["gene"] in genes]
    identity = {key: sum(int(s[key]) for s in rows) for key in ("tp", "fp", "fn")}
    identity["counted_extra_rows"] = sum(
        int(s["counted_precision"]["counted_extra_rows"]) for s in rows
    )
    counts: dict[str, Any] = {}
    for field in FIELDS:
        errors = [e for s in rows for e in s["count_errors"] if e["field"] == field]
        positive = sum(s["count"][field]["gold_asserted_nonzero"] for s in rows)
        supplied = sum(s["count"][field]["predicted"] for s in rows)
        positive_supplied = sum(
            s["count"][field]["predicted_on_nonzero_gold"] for s in rows
        )
        positive_exact = positive_supplied - sum(1 for e in errors if e["gold"] > 0)
        counts[field] = {
            "gold_assertions": sum(s["count"][field]["gold_asserted"] for s in rows),
            "positive_gold": positive,
            "supplied": supplied,
            "positive_supplied": positive_supplied,
            "positive_exact": positive_exact,
            "wrong_supplied": len(errors),
            "exact_supplied": supplied - len(errors),
            "positive_exact_recovery_pct": (
                round(100 * positive_exact / positive, 2) if positive else None
            ),
            "conditional_exact_pct": (
                round(100 * (supplied - len(errors)) / supplied, 2)
                if supplied
                else None
            ),
        }
    return {"attempts": len(rows), "identity": identity, "counts": counts}


def run(args: argparse.Namespace) -> dict[str, Any]:
    run_eval = load_run_eval()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    pooled: dict[str, list[dict]] = {"off": [], "on": []}
    pooled_new_values: list[dict[str, Any]] = []
    report: dict[str, Any] = {
        "mode": args.mode,
        "paper_tier": bool(args.paper_tier),
        "runs": {},
    }

    for run_dir in args.run_dir:
        run_dir = run_dir.resolve()
        run_id = run_dir.name
        selection = json.loads((run_dir / "selection.json").read_text())
        predictions = json.loads((run_dir / "predictions.json").read_text())
        setup = json.loads((run_dir / "setup.json").read_text())
        gold_root = Path(
            args.gold_root or (setup.get("cohort") or {}).get("gold_root") or ""
        )
        if not gold_root.is_dir():
            raise SystemExit(f"{run_id}: gold root not found: {gold_root}")
        production_root = run_dir / "production_runs"
        gene_runs: dict[str, Path] = {}
        for gene in sorted({p["gene"] for p in selection["papers"]}):
            try:
                gene_run, _db, _status = resolve_active_gene_run(production_root, gene)
            except ProductionRunError as exc:
                print(f"[{run_id}] {gene}: {exc}", file=sys.stderr)
                continue
            gene_runs[gene] = gene_run

        base_papers = {(p["gene"], str(p["pmid"])): p for p in predictions["papers"]}
        patched_papers = copy.deepcopy(base_papers)
        paper_logs: dict[str, Any] = {}
        module_tables: dict[str, Any] = {}
        for paper in selection["papers"]:
            key = (paper["gene"], str(paper["pmid"]))
            gene, pmid = key
            gene_run = gene_runs.get(gene)
            if gene_run is None or key not in patched_papers:
                continue
            path = gene_run / "extractions" / f"{gene}_PMID_{pmid}.json"
            if not path.is_file():
                continue
            try:
                extraction = json.loads(path.read_text())
            except json.JSONDecodeError:
                continue
            records, metadata = derived_records(
                extraction,
                source_text_for(extraction, gene_run, pmid),
                gene=gene,
                disease=disease_for(gene_run),
                mode=args.mode,
                paper_tier=bool(args.paper_tier),
            )
            if metadata.get("tables"):
                module_tables[f"{gene}:{pmid}"] = metadata["tables"]
            if not records:
                continue
            log = patch_paper(patched_papers[key], records, args.mode)
            if log:
                paper_logs[f"{gene}:{pmid}"] = log

        # Score both arms.
        arms = {"off": base_papers, "on": patched_papers}
        if args.mode == "strip":
            arms = {"off": patched_papers, "on": base_papers}
        scores: dict[str, dict[tuple[str, str], dict]] = {"off": {}, "on": {}}
        gold_cache: dict[tuple[str, str], list[dict]] = {}
        for arm, papers in arms.items():
            for paper in selection["papers"]:
                key = (paper["gene"], str(paper["pmid"]))
                predicted = papers.get(key)
                if predicted is None:
                    continue
                if key not in gold_cache:
                    gold_cache[key] = run_eval.load_gold(gold_root, *key)
                scores[arm][key] = run_eval.score_one(*key, predicted, gold_cache[key])
                pooled[arm].append(scores[arm][key])

        # Newly supplied values: exact vs wrong, with gold row context.
        new_values: list[dict[str, Any]] = []
        for key, log in paper_logs.items():
            gene, pmid = key.split(":")
            off_score = scores["off"][(gene, pmid)]
            on_score = scores["on"][(gene, pmid)]
            off_errors = {
                (e["variant"], e["field"]): e for e in off_score["count_errors"]
            }
            on_errors = {
                (e["variant"], e["field"]): e for e in on_score["count_errors"]
            }
            gold = gold_cache[(gene, pmid)]
            for entry in log:
                if entry.get("status"):
                    continue
                gold_row = gold_for_variant(run_eval, gene, entry["variant"], gold)
                error_key = (entry["variant"], entry["field"])
                supplied_value = (
                    entry["after"] if args.mode == "derive" else entry["before"]
                )
                if gold_row is None or gold_row.get(entry["field"]) is None:
                    verdict = "unmatched_or_gold_null"
                elif error_key in on_errors and error_key not in off_errors:
                    verdict = "wrong"
                elif error_key in off_errors and error_key not in on_errors:
                    verdict = "wrong_removed"
                elif gold_row.get(entry["field"]) == supplied_value:
                    verdict = "exact"
                else:
                    verdict = "unscored"
                new_values.append(
                    {
                        "run_id": run_id,
                        "gene": gene,
                        "pmid": pmid,
                        "variant": entry["variant"],
                        "field": entry["field"],
                        "value": supplied_value,
                        "gold_carriers": gold_row.get("carriers") if gold_row else None,
                        "gold_affected": gold_row.get("affected") if gold_row else None,
                        "gold_unaffected": gold_row.get("unaffected")
                        if gold_row
                        else None,
                        "verdict": verdict,
                        "gold_real_split": (
                            gold_row is not None
                            and gold_row.get("carriers") is not None
                            and gold_row.get("affected") is not None
                            and gold_row["affected"] != gold_row["carriers"]
                        ),
                        "role": entry["role"],
                        "tier": entry["tier"],
                        "caption": entry["caption"],
                        "count_column": entry["count_column"],
                        "source_location": entry["source_location"],
                    }
                )
        pooled_new_values.extend(new_values)

        per_paper = []
        for key in sorted(scores["on"]):
            off_s, on_s = scores["off"][key], scores["on"][key]
            changed = any(
                off_s["count"][f]["predicted"] != on_s["count"][f]["predicted"]
                or len([e for e in off_s["count_errors"] if e["field"] == f])
                != len([e for e in on_s["count_errors"] if e["field"] == f])
                for f in FIELDS
            )
            if not changed:
                continue
            per_paper.append(
                {
                    "gene": key[0],
                    "pmid": key[1],
                    "identity_tp_fp_fn": [on_s["tp"], on_s["fp"], on_s["fn"]],
                    **{
                        f"{f}_supplied_off_on": [
                            off_s["count"][f]["predicted"],
                            on_s["count"][f]["predicted"],
                        ]
                        for f in FIELDS
                    },
                    **{
                        f"{f}_wrong_off_on": [
                            len([e for e in off_s["count_errors"] if e["field"] == f]),
                            len([e for e in on_s["count_errors"] if e["field"] == f]),
                        ]
                        for f in FIELDS
                    },
                    "tables": module_tables.get(f"{key[0]}:{key[1]}", {}),
                }
            )

        verdicts = Counter((v["field"], v["verdict"]) for v in new_values)
        real_split_hits = [
            v for v in new_values if v["field"] == "affected" and v["gold_real_split"]
        ]
        report["runs"][run_id] = {
            "inputs": {
                "predictions_sha256": sha256(run_dir / "predictions.json"),
                "report_sha256": sha256(run_dir / "report.json")
                if (run_dir / "report.json").is_file()
                else None,
                "gold_root": str(gold_root),
            },
            "arms": {
                arm: {
                    "all_genes": summarize_arm(list(scores[arm].values())),
                    "cardiac_four": summarize_arm(
                        list(scores[arm].values()), set(CARDIAC)
                    ),
                }
                for arm in ("off", "on")
            },
            "new_value_verdicts": {
                f"{f}:{v}": n for (f, v), n in sorted(verdicts.items())
            },
            "new_affected_on_gold_real_split_rows": len(real_split_hits),
            "papers_changed": per_paper,
        }
        (out_dir / f"new_values_{run_id}.csv").write_text("")
        if new_values:
            with (out_dir / f"new_values_{run_id}.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(new_values[0]))
                writer.writeheader()
                writer.writerows(new_values)

    report["pooled"] = {
        arm: {
            "all_genes": summarize_arm(pooled[arm]),
            "cardiac_four": summarize_arm(pooled[arm], set(CARDIAC)),
        }
        for arm in ("off", "on")
    }
    delta_set = Counter((v["field"], v["verdict"]) for v in pooled_new_values)
    report["pooled"]["new_value_verdicts"] = {
        f"{f}:{v}": n for (f, v), n in sorted(delta_set.items())
    }
    report["pooled"]["new_affected_on_gold_real_split_rows"] = sum(
        1
        for v in pooled_new_values
        if v["field"] == "affected" and v["gold_real_split"]
    )
    (out_dir / "replay_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print_report(report)
    return report


def _fmt(value: Any) -> str:
    return "n/a" if value is None else str(value)


def print_report(report: dict[str, Any]) -> None:
    for scope in ("cardiac_four", "all_genes"):
        print(
            f"\n== pooled {scope} ({report['mode']}, paper_tier={report['paper_tier']}) =="
        )
        off = report["pooled"]["off"][scope]
        on = report["pooled"]["on"][scope]
        print(
            f"identity off/on: TP {off['identity']['tp']}/{on['identity']['tp']} "
            f"FP {off['identity']['fp']}/{on['identity']['fp']} "
            f"FN {off['identity']['fn']}/{on['identity']['fn']} "
            f"counted extras {off['identity']['counted_extra_rows']}/"
            f"{on['identity']['counted_extra_rows']}"
        )
        header = (
            f"{'field':<11}{'pos gold':>9}{'supplied off/on':>18}{'pos exact off/on':>19}"
            f"{'exact rec% off/on':>20}{'wrong off/on':>14}{'cond% off/on':>16}"
        )
        print(header)
        for field in FIELDS:
            a, b = off["counts"][field], on["counts"][field]
            print(
                f"{field:<11}{a['positive_gold']:>9}"
                f"{str(a['supplied']) + '/' + str(b['supplied']):>18}"
                f"{str(a['positive_exact']) + '/' + str(b['positive_exact']):>19}"
                f"{_fmt(a['positive_exact_recovery_pct']) + '/' + _fmt(b['positive_exact_recovery_pct']):>20}"
                f"{str(a['wrong_supplied']) + '/' + str(b['wrong_supplied']):>14}"
                f"{_fmt(a['conditional_exact_pct']) + '/' + _fmt(b['conditional_exact_pct']):>16}"
            )
    print("\nnew value verdicts (pooled):", report["pooled"]["new_value_verdicts"])
    print(
        "new affected on gold real-split rows:",
        report["pooled"]["new_affected_on_gold_real_split_rows"],
    )


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, action="append", required=True)
    parser.add_argument("--gold-root", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--mode", choices=("derive", "strip"), default="derive")
    parser.add_argument(
        "--paper-tier",
        action="store_true",
        help="also enable the paper-ascertainment tier (diagnostic only)",
    )
    args = parser.parse_args(argv)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
