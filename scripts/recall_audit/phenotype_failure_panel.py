#!/usr/bin/env python3
"""Read-only, source-audited phenotype failure panel from explicitly opened locks.

The panel is diagnostic/calibration evidence, never a new held-out benchmark.
Count errors are per field; identity extras are exported separately because
the gold does not assert a count for an unmatched prediction. No cohort pooling,
gold edits, extraction, source mutation, or conversion of stored NULL to zero.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts import build_gold_difference_figure as difference  # noqa: E402


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summarize_paper(
    rows: list[dict[str, Any]], paper: dict[str, Any], run_id: str, lane: str
) -> dict[str, Any]:
    """Separate absent identities, absent counts, and supplied discrepancies."""
    result: dict[str, Any] = {
        "run_id": run_id,
        "score_lane": lane,
        "gene": paper["gene"],
        "pmid": paper["pmid"],
        **{key: paper[key] for key in ("tp", "fp", "fn")},
    }
    phenotype = [r for r in rows if r["measure"] in ("affected", "unaffected")]
    supplied = [r for r in phenotype if r["status"] == "supplied"]
    result["phenotype_capture"] = (
        "none_asserted"
        if not phenotype
        else "all_na"
        if not supplied
        else "partial"
        if len(supplied) < len(phenotype)
        else "all_supplied"
    )
    result["phenotype_abs_error"] = sum(r["absolute_difference"] for r in phenotype)
    result["phenotype_supplied_fields"] = len(supplied)
    result["phenotype_positive_supplied_fields"] = sum(
        r["automated_count_raw"] > 0 for r in supplied
    )
    result["phenotype_wrong_supplied_fields"] = sum(
        r["absolute_difference"] > 0 for r in supplied
    )
    for field in difference.FIELDS:
        items = [r for r in rows if r["measure"] == field]
        values = [r for r in items if r["status"] == "supplied"]
        result.update(
            {
                f"{field}_asserted": len(items),
                f"{field}_supplied": len(values),
                f"{field}_supplied_exact": sum(
                    r["absolute_difference"] == 0 for r in values
                ),
                f"{field}_wrong_supplied": sum(
                    r["absolute_difference"] > 0 for r in values
                ),
                f"{field}_abs_error": sum(r["absolute_difference"] for r in items),
                f"{field}_missing_positive_fields": sum(
                    r["status"] != "supplied" and r["gold_count"] > 0 for r in items
                ),
                f"{field}_missing_gold_zero_fields": sum(
                    r["status"] != "supplied" and r["gold_count"] == 0 for r in items
                ),
                # Magnitude, NOT numbers of people proved misclassified. Gold can
                # encode a different phenotype, ascertainment group, or timepoint.
                f"{field}_overcount_units": sum(
                    max(0, r["difference"]) for r in values
                ),
                f"{field}_supplied_undercount_units": sum(
                    max(0, -r["difference"]) for r in values
                ),
            }
        )
        for status in difference.STATUSES:
            result[f"{field}_{status}_error"] = sum(
                r["absolute_difference"] for r in items if r["status"] == status
            )
    return result


def source_inventory(run: Path, gene: str, pmid: str) -> list[dict[str, Any]]:
    """Inventory observed bytes; length/table density is not completeness proof."""
    rows = []
    directories = [("current_corpus", REPO / "corpus" / gene / pmid)]
    for production in sorted((run / "production_runs" / gene).glob("*/")):
        directories.append(("run_artifact", production / "pmc_fulltext"))
    for origin, directory in directories:
        paths = set(directory.glob(f"{pmid}*"))
        paths |= {
            p for folder in list(paths) if folder.is_dir() for p in folder.rglob("*")
        }
        for path in sorted(paths):
            if not path.is_file():
                continue
            row = {
                "run_id": run.name,
                "gene": gene,
                "pmid": pmid,
                "origin": origin,
                "path": str(path.relative_to(REPO)),
                "bytes": path.stat().st_size,
                "sha256": difference.sha256(path),
                "historical_lock_binding": "not independently verified; inventory observed at audit time",
            }
            if path.suffix.lower() in {".md", ".txt"}:
                text = path.read_text(errors="replace")
                row.update(
                    {
                        "characters": len(text),
                        "pipe_table_lines": sum(
                            line.count("|") >= 3 for line in text.splitlines()
                        ),
                        "explicit_abstract_fallback": "ABSTRACT-ONLY FALLBACK"
                        in text[:1500],
                        "explicit_supplement_only": "supplement-only recovery"
                        in text[:1500],
                    }
                )
            rows.append(row)
    return rows


def probe_stages(
    run: Path, gene: str, pmid: str, probe: dict[str, Any]
) -> dict[str, Any]:
    """Phrase/field audit only; a phrase hit does not prove clinical support."""
    from pipeline.table_router import enumerate_markdown_tables
    from scripts.recall_audit.fn_root_cause import RunPaper
    from scripts.recall_audit.tier_source_reachability import _squash

    paper = RunPaper(run, gene, pmid)
    result: dict[str, Any] = {
        "run_id": run.name,
        "gene": gene,
        "pmid": pmid,
        "interpretation": "observed trace/artifact probe; not causal proof or an endpoint adjudication",
        "stages": {},
        "extractions": [],
        "trace_inputs": [],
    }
    if paper.prod is None:
        result["status"] = "production_artifacts_missing"
        return result
    phrases = probe["phrases"]
    for folder, suffix in (
        ("pmc_fulltext", "FULL_CONTEXT"),
        ("pmc_fulltext", "CLEANED"),
        ("scout_output", "DATA_ZONES"),
    ):
        path = paper.prod / folder / f"{pmid}_{suffix}.md"
        if path.exists():
            text = path.read_text(errors="replace")
            result["stages"][suffix] = {
                "path": str(path.relative_to(REPO)),
                "sha256": difference.sha256(path),
                "characters": len(text),
                "markdown_tables": len(enumerate_markdown_tables(text)),
                "phrase_present": {p: _squash(p) in _squash(text) for p in phrases},
            }
    for stage, text in (
        ("llm_request", paper.llm_request()),
        ("llm_response", paper.llm_response()),
    ):
        result["stages"][stage] = {
            "characters_squashed": len(text),
            "phrase_present": {p: _squash(p) in text for p in phrases},
        }
    for path in sorted((paper.prod / "llm_traces" / gene / pmid).glob("*.json")):
        record = difference.read_json(path)
        result["trace_inputs"].append(
            {
                "path": str(path.relative_to(REPO)),
                "sha256": difference.sha256(path),
                "record_type": record.get("record_type"),
            }
        )
    for path in sorted((paper.prod / "extractions").glob(f"*{pmid}*.json")):
        data = difference.read_json(path)
        metadata = data.get("extraction_metadata") or {}
        selected = []
        for variant in data.get("variants", []):
            identity = " ".join(
                str(variant.get(k) or "")
                for k in (
                    "protein_notation",
                    "cdna_notation",
                    "source_notation",
                    "legacy_notation",
                )
            )
            if not any(
                _squash(n) in _squash(identity) for n in probe["identity_needles"]
            ):
                continue
            selected.append(
                {
                    k: variant.get(k)
                    for k in (
                        "protein_notation",
                        "cdna_notation",
                        "patients",
                        "penetrance_data",
                        "count_provenance",
                        "phenotype_count_flags",
                        "claim_verification",
                        "fact_provenance",
                    )
                }
            )
        result["extractions"].append(
            {
                "path": str(path.relative_to(REPO)),
                "sha256": difference.sha256(path),
                "metadata_notes": metadata.get("notes"),
                "variants": selected,
            }
        )
    return result


def build(panel_path: Path, out: Path) -> dict[str, Any]:
    config = difference.read_json(panel_path)
    selected = {str(p["pmid"]): p for p in config["papers"]}
    if len(selected) != len(config["papers"]):
        raise ValueError("panel must select unique PMIDs")
    if not (REPO / "corpus").is_symlink() or not (REPO / "corpus").is_dir():
        raise ValueError("mount Ezekers: corpus must be the existing reachable symlink")
    all_papers, selected_papers, selected_rows, extras, inventory = [], [], [], [], []
    runs, seen, stage_probes = [], set(), []
    hashes = {str(panel_path.relative_to(REPO)): difference.sha256(panel_path)}
    for spec in config["runs"]:
        run = REPO / spec["path"]
        lock = difference.verify_lock(run)
        report = difference.read_json(run / "report.json")
        predictions = difference.read_json(run / "predictions.json")
        gold = difference.resolve_gold_sources(report, None)
        rows = difference.build_rows("diagnostic", report, predictions, gold)
        metrics = {
            field: difference.summarize_field(
                [r for r in rows if r["measure"] == field],
                report["overall"]["count"][field],
            )
            for field in difference.FIELDS
        }
        runs.append(
            {
                "run_id": run.name,
                "score_lane": spec["lane"],
                "locked_at": lock["locked_at"],
                "metrics": metrics,
            }
        )
        for path in [
            run / n
            for n in ("LOCK.json", "report.json", "predictions.json", "selection.json")
        ] + list(gold.values()):
            hashes[str(path.relative_to(REPO))] = difference.sha256(path)
        grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            grouped[(row["gene"], row["pmid"])].append(row)
        pred_map = {(p["gene"], p["pmid"]): p for p in predictions["papers"]}
        for paper in report["papers"]:
            gene, pmid = paper["gene"], paper["pmid"]
            summary = summarize_paper(
                grouped[(gene, pmid)], paper, run.name, spec["lane"]
            )
            all_papers.append(summary)
            if pmid not in selected:
                continue
            seen.add(pmid)
            selected_papers.append(
                {**summary, "selection_role": selected[pmid]["role"]}
            )
            selected_rows.extend(grouped[(gene, pmid)])
            inventory.extend(source_inventory(run, gene, pmid))
            if probe := selected[pmid].get("probe"):
                stage_probes.append(probe_stages(run, gene, pmid, probe))
            predicted = difference.predicted_count_rows(pred_map[(gene, pmid)], gene)
            for name in paper["extra_predictions"]:
                value = predicted.get(name, {})
                raw = value.get("row", {})
                extras.append(
                    {
                        "run_id": run.name,
                        "score_lane": spec["lane"],
                        "gene": gene,
                        "pmid": pmid,
                        "variant": name,
                        **{f: value.get(f) for f in difference.FIELDS},
                        "source_layer": raw.get("source_layer"),
                        "source_location": raw.get("source_location"),
                        "interpretation": "unmatched vs fixture; not adjudicated biological false positive",
                    }
                )
    missing = set(selected) - seen
    if missing:
        raise ValueError(
            f"selected PMIDs absent from opened reports: {sorted(missing)}"
        )
    all_papers.sort(
        key=lambda p: (-p["phenotype_abs_error"], p["run_id"], p["gene"], p["pmid"])
    )
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "unique_panel_pmids": len(selected),
        "panel_gene_run_attempts": len(selected_papers),
        "available_gene_run_attempts": len(all_papers),
        "purpose": "post-score, deliberately enriched diagnostic panel; not confirmation or population estimate",
        "count_rule": "abs(predicted-gold); missing identity/count evaluates to zero only; stored NULL preserved",
        "fp_rule": "identity extras separate from count overestimation; no gold count exists for an extra",
        "gold_rule": "scorer-authoritative gold-v2 values; zero may mean not reported, not verified healthy",
        "source_rule": "current corpus and run source files inventoried separately; not assumed identical or lock-bound",
        "runs": runs,
        "input_sha256": hashes,
    }
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (
        ("all_papers.csv", all_papers),
        ("panel.csv", selected_papers),
        ("panel_rows.csv", selected_rows),
        ("extra_predictions.csv", extras),
        ("source_inventory.csv", inventory),
    ):
        write_csv(out / name, data)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "stage_probes.json").write_text(json.dumps(stage_probes, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    manifest = build(args.panel.resolve(), args.out_dir.resolve())
    print(
        json.dumps(
            {
                k: manifest[k]
                for k in (
                    "unique_panel_pmids",
                    "panel_gene_run_attempts",
                    "available_gene_run_attempts",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
