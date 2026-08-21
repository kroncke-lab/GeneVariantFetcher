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
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff"}
TEXT_REPRESENTATION_SUFFIXES = (
    "_FULL_CONTEXT.md",
    "_CLEANED.md",
    "_DATA_ZONES.md",
)
REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from benchmarks.codex_paper_eval.production_run import (  # noqa: E402
    ProductionRunError,
    resolve_active_gene_run,
)


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


def extraction_source(
    gene_run: Path, gene: str, pmid: str, fallback: Path
) -> tuple[Path, Path | None, str]:
    """Resolve and verify the exact primary source recorded by extraction.

    Successful extraction JSON stores the source path, byte size, and SHA-256
    that ``pipeline.steps`` computed immediately before reading it. Papers that
    never produced an extraction artifact retain an explicit fallback binding;
    their empty prediction is still trace/DB-bound, but must not be described as
    having a proven extraction input.
    """
    record = gene_run / "extractions" / f"{gene}_PMID_{pmid}.json"
    if not record.is_file():
        return fallback.resolve(), None, "no_extraction_record_fallback"
    try:
        payload = json.loads(record.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(
            f"invalid extraction record for {gene}:{pmid}: {record}: {exc}"
        )
    metadata = payload.get("extraction_metadata") or {}
    raw_source = str(metadata.get("source_file") or "").strip()
    expected_sha = str(metadata.get("source_sha256") or "").strip()
    expected_size = metadata.get("source_size_bytes")
    if not raw_source or not expected_sha:
        raise SystemExit(
            f"extraction record lacks source path/SHA for {gene}:{pmid}: {record}"
        )
    source = Path(raw_source)
    if not source.is_absolute():
        candidates = (gene_run / source, REPO / source)
        source = next((path for path in candidates if path.is_file()), candidates[-1])
    source = source.resolve()
    if not source.is_file():
        raise SystemExit(
            f"recorded extraction source missing for {gene}:{pmid}: {source}"
        )
    actual_sha = digest(source)
    if actual_sha != expected_sha:
        raise SystemExit(
            f"recorded extraction source SHA mismatch for {gene}:{pmid}: {source}"
        )
    if expected_size is not None and source.stat().st_size != int(expected_size):
        raise SystemExit(
            f"recorded extraction source size mismatch for {gene}:{pmid}: {source}"
        )
    return source, record.resolve(), "extraction_metadata_verified"


def latest_gene_run(production_root: Path, gene: str) -> Path:
    """Compatibility name for exact active-run resolution without mtimes."""
    try:
        run_dir, _db, _status = resolve_active_gene_run(production_root, gene)
    except ProductionRunError as exc:
        raise SystemExit(str(exc)) from exc
    return run_dir


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
    binding_statuses: Counter[str] = Counter()
    source_representations: Counter[str] = Counter()
    for paper in selection["papers"]:
        gene, pmid = paper["gene"], str(paper["pmid"])
        if gene not in gene_runs:
            gene_runs[gene] = latest_gene_run(args.production_root, gene)
        gene_run = gene_runs[gene]
        pmc_dir = gene_run / "pmc_fulltext"
        full_context = (pmc_dir / f"{pmid}_FULL_CONTEXT.md").resolve()
        if not full_context.is_file():
            raise SystemExit(
                f"production FULL_CONTEXT missing for {gene}:{pmid}: {full_context}"
            )
        source, extraction_record, binding_status = extraction_source(
            gene_run, gene, pmid, full_context
        )
        binding_statuses[binding_status] += 1
        source_representations[
            next(
                (
                    suffix.removeprefix("_")
                    for suffix in TEXT_REPRESENTATION_SUFFIXES
                    if source.name.endswith(suffix)
                ),
                source.suffix.lower() or "other",
            )
        ] += 1
        if digest(source) != paper.get("source_sha256"):
            changed_sources += 1

        files = material_files(pmc_dir, pmid)
        artifact = (pmc_dir / f"{pmid}_artifacts.json").resolve()
        pdfs = [path for path in files if path.suffix.lower() == ".pdf"]
        figures = [path for path in files if path.suffix.lower() in IMAGE_SUFFIXES]
        representations = [
            path
            for path in files
            if any(
                path.name.endswith(suffix) for suffix in TEXT_REPRESENTATION_SUFFIXES
            )
        ]
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
                    "policy": binding_status,
                    "selected": str(source),
                    "selected_metrics": source_metrics,
                    "candidates": [
                        {
                            "path": str(path),
                            "sha256": digest(path),
                            "characters": len(path.read_text(errors="replace")),
                        }
                        for path in representations
                    ],
                    "rationale": (
                        "Bound after extraction and before gold scoring to the "
                        "source path/SHA persisted by gvf-run extraction metadata. "
                        "All run-local text representations available to the "
                        "multi-stage extractor are hashed separately. A missing "
                        "extraction record is labeled as a fallback, never as a "
                        "verified read. Cohort membership and gold values were "
                        "unchanged."
                    ),
                },
                "production_extraction_record": (
                    str(extraction_record) if extraction_record else None
                ),
                "production_extraction_record_sha256": (
                    digest(extraction_record) if extraction_record else None
                ),
                "representations": [str(path) for path in representations],
                "representation_sha256": {
                    str(path): digest(path) for path in representations
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
        "binding_status_counts": dict(sorted(binding_statuses.items())),
        "selected_representation_counts": dict(sorted(source_representations.items())),
        "gold_values_read": False,
        "cohort_membership_changed": False,
    }
    selection_path.write_text(json.dumps(selection, indent=2, sort_keys=True) + "\n")
    print(
        f"{selection_path}: rebound {len(selection['papers'])} papers; "
        f"{changed_sources} source digests differed from eligibility snapshot; "
        f"bindings={dict(sorted(binding_statuses.items()))}; "
        f"sources={dict(sorted(source_representations.items()))}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
