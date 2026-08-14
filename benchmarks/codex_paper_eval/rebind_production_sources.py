#!/usr/bin/env python3
"""Bind an evaluation selection to the exact material a production run read.

``gvf-run`` stages source into a run-local ``pmc_fulltext`` directory and may
refold supplements or upgrade a weak cached rendering before extraction.  A
production projection must therefore lock those run-local inputs, not merely
the corpus files that made a PMID eligible during ``prepare``.

This is a source-only, pre-score operation.  It does not open a gold file or
change cohort membership.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def material_files(pmc_dir: Path, pmid: str) -> list[Path]:
    files: list[Path] = []
    for candidate in sorted(pmc_dir.glob(f"{pmid}_*")):
        if candidate.is_file():
            files.append(candidate.resolve())
        elif candidate.is_dir():
            files.extend(
                sorted(
                    path.resolve() for path in candidate.rglob("*") if path.is_file()
                )
            )
    return files


def latest_gene_run(production_root: Path, gene: str) -> Path:
    candidates = [path for path in (production_root / gene).iterdir() if path.is_dir()]
    if not candidates:
        raise SystemExit(f"no production run for {gene} under {production_root}")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--production-root", type=Path, required=True)
    args = parser.parse_args()

    selection_path = args.run_dir / "selection.json"
    backup_path = args.run_dir / "selection_eligibility.json"
    if backup_path.exists():
        raise SystemExit(f"refusing to overwrite existing backup: {backup_path}")
    original_sha = digest(selection_path)
    selection = json.loads(selection_path.read_text())
    shutil.copy2(selection_path, backup_path)

    gene_runs: dict[str, Path] = {}
    changed_sources = 0
    for paper in selection["papers"]:
        gene, pmid = paper["gene"], str(paper["pmid"])
        if gene not in gene_runs:
            gene_runs[gene] = latest_gene_run(args.production_root, gene)
        gene_run = gene_runs[gene]
        pmc_dir = gene_run / "pmc_fulltext"
        source = (pmc_dir / f"{pmid}_FULL_CONTEXT.md").resolve()
        if not source.is_file():
            raise SystemExit(f"production source missing for {gene}:{pmid}: {source}")
        if digest(source) != paper.get("source_sha256"):
            changed_sources += 1

        files = material_files(pmc_dir, pmid)
        artifact = (pmc_dir / f"{pmid}_artifacts.json").resolve()
        pdfs = [path for path in files if path.suffix.lower() == ".pdf"]
        figures = [path for path in files if path.suffix.lower() in IMAGE_SUFFIXES]
        source_metrics = {
            "path": str(source),
            "sha256": digest(source),
            "characters": len(source.read_text(errors="replace")),
        }
        paper.update(
            {
                "source": str(source),
                "source_sha256": source_metrics["sha256"],
                "source_bytes": source.stat().st_size,
                "source_selection": {
                    "policy": "production_input_snapshot",
                    "selected": str(source),
                    "selected_metrics": source_metrics,
                    "candidates": [source_metrics],
                    "rationale": (
                        "Bound after extraction, before gold scoring, to the exact "
                        "run-local FULL_CONTEXT input staged by gvf-run. Cohort "
                        "membership and gold values were unchanged."
                    ),
                },
                "artifacts": str(artifact) if artifact.is_file() else None,
                "artifacts_sha256": digest(artifact) if artifact.is_file() else None,
                "pdfs": [str(path) for path in pdfs],
                "pdf_sha256": {str(path): digest(path) for path in pdfs},
                "figures": [str(path) for path in figures],
                "figure_sha256": {str(path): digest(path) for path in figures},
                "production_run_dir": str(gene_run.resolve()),
            }
        )

    selection["production_source_rebinding"] = {
        "rebound_at": datetime.now(timezone.utc).isoformat(),
        "original_selection": str(backup_path.resolve()),
        "original_selection_sha256": original_sha,
        "changed_source_count": changed_sources,
        "paper_count": len(selection["papers"]),
        "gold_values_read": False,
        "cohort_membership_changed": False,
    }
    selection_path.write_text(json.dumps(selection, indent=2, sort_keys=True) + "\n")
    print(
        f"{selection_path}: rebound {len(selection['papers'])} papers; "
        f"{changed_sources} source digests differed from eligibility snapshot"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
