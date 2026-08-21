"""Fail-closed resolution of one production run used by the paper eval."""

from __future__ import annotations

import json
from pathlib import Path


class ProductionRunError(RuntimeError):
    """A production root cannot be resolved to one successful active run."""


def resolve_active_gene_run(
    production_root: Path, gene: str
) -> tuple[Path, Path, dict]:
    """Return ``(run_dir, active_db, status)`` without using mtimes.

    A fresh evaluation should create exactly one completed ``gvf-run`` per gene.
    If a root contains history, callers must narrow the root rather than letting
    a timestamp or backup filename silently decide which science is scored.
    """
    gene_root = production_root / gene
    if not gene_root.is_dir():
        raise ProductionRunError(f"no production root for {gene}: {gene_root}")
    candidates: dict[Path, tuple[Path, dict]] = {}
    for status_path in sorted(gene_root.glob("**/RUN_STATUS.json")):
        try:
            status = json.loads(status_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if (
            status.get("status") != "completed"
            or status.get("exit_code") != 0
            or status.get("stage_failures")
        ):
            continue
        active_raw = str(status.get("active_db") or "").strip()
        active_db = Path(active_raw).expanduser()
        if not active_db.is_absolute():
            active_db = status_path.parent / active_db
        active_db = active_db.resolve()
        if active_db.name != f"{gene}.db" or not active_db.is_file():
            continue
        candidates[status_path.parent.resolve()] = (active_db, status)
    if not candidates:
        raise ProductionRunError(
            f"no completed {gene} run declares an existing active {gene}.db under {gene_root}"
        )
    if len(candidates) != 1:
        raise ProductionRunError(
            f"multiple completed {gene} runs under {gene_root}; narrow the production "
            f"root or remove ambiguity: {', '.join(str(path) for path in candidates)}"
        )
    run_dir, (active_db, status) = next(iter(candidates.items()))
    return run_dir, active_db, status
