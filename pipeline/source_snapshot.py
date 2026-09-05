"""Explicit source-only snapshots for repeatable, acquisition-free evaluations.

This mode preserves incomplete sources as incomplete. It does not promote them
to successful acquisition or allow ordinary production cache reuse to skip a
needed supplement retry.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def freeze_sources(papers: list[dict], output_dir: Path) -> Path:
    """Copy selected text and its source assets; never copy predictions or DBs."""
    output_dir.mkdir(parents=True, exist_ok=False)
    records = []
    for paper in papers:
        gene, pmid = paper["gene"], str(paper["pmid"])
        source = Path(paper["source"])
        expected = paper.get("source_sha256")
        if expected and file_sha256(source) != expected:
            raise ValueError(f"source changed before freezing: {gene}/{pmid}")
        destination = output_dir / gene / pmid
        destination.mkdir(parents=True)
        # The selected rendering is the extraction input, regardless of its
        # old suffix. Preprocessing in each arm still runs normally.
        shutil.copy2(source, destination / f"{pmid}_FULL_CONTEXT.md")
        parent = source.parent
        artifact = parent / f"{pmid}_artifacts.json"
        if artifact.is_file():
            shutil.copy2(artifact, destination / artifact.name)
        for suffix in ("_figures", "_supplements"):
            folder = parent / f"{pmid}{suffix}"
            if folder.is_dir():
                shutil.copytree(folder, destination / folder.name)
        # Raw PDFs next to the article can be read by figure/table routes.
        for pdf in parent.glob(f"{pmid}*.pdf"):
            shutil.copy2(pdf, destination / pdf.name)
        files = {
            path.relative_to(destination).as_posix(): file_sha256(path)
            for path in sorted(destination.rglob("*"))
            if path.is_file()
        }
        records.append({"gene": gene, "pmid": pmid, "files": files})
    manifest = output_dir / "source_snapshot.json"
    manifest.write_text(
        json.dumps({"schema_version": 1, "papers": records}, indent=2) + "\n"
    )
    return manifest


def stage_frozen_sources(
    manifest: Path, gene: str, pmids: list[str], harvest_dir: Path
) -> set[str]:
    """Verify and copy only requested source bytes; fail rather than fetch."""
    data = json.loads(manifest.read_text())
    if data.get("schema_version") != 1:
        raise ValueError("unsupported source snapshot schema")
    selected = {(row["gene"], str(row["pmid"])): row for row in data["papers"]}
    copied = set()
    for pmid in map(str, pmids):
        row = selected.get((gene.upper(), pmid))
        if row is None:
            raise ValueError(f"missing frozen source: {gene}/{pmid}")
        root = manifest.parent / gene.upper() / pmid
        files = row["files"]
        if f"{pmid}_FULL_CONTEXT.md" not in files:
            raise ValueError(f"snapshot lacks selected text: {gene}/{pmid}")
        for relative, expected in files.items():
            source = root / relative
            if Path(relative).is_absolute() or ".." in Path(relative).parts:
                raise ValueError("invalid snapshot member path")
            if not source.is_file() or file_sha256(source) != expected:
                raise ValueError(f"frozen source changed: {gene}/{pmid}/{relative}")
            destination = harvest_dir / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
        copied.add(pmid)
    return copied
