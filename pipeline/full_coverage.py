"""Walk-to-taper extraction for full literature coverage.

Drives ``scripts/run_priority_extraction.py`` in fixed offset batches until the
per-batch variant yield tapers (DB ``COUNT(DISTINCT variant_id)`` growth below a
threshold for ``patience`` consecutive batches), instead of a single bounded
top-N priority pass. Extraction is idempotent (already-extracted PMIDs are
reused), so re-covering early ranks is cheap and closes gaps.

Used by ``gvf-run --full-coverage``. Self-contained: no edits to the extraction
pipeline — it composes the existing priority-extraction entry point and reads
the run DB for the taper signal. Sets ``TIER3_MODELS``/``MAX_WORKERS`` for the
batch subprocess so full coverage can use a high-TPM model (default gpt-5.4)
without changing global config.
"""

from __future__ import annotations

import os
import json
import sqlite3
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

ROOT = Path(__file__).resolve().parent.parent
PRIORITY = ROOT / "scripts" / "run_priority_extraction.py"


def _distinct_variants(db: Path) -> int:
    if not db.exists():
        return 0
    con = sqlite3.connect(str(db))
    try:
        return int(
            con.execute(
                "SELECT COUNT(DISTINCT variant_id) FROM variant_papers"
            ).fetchone()[0]
            or 0
        )
    except Exception:
        return 0
    finally:
        con.close()


def run_walk_to_taper(
    gene: str,
    run_dir: Path,
    *,
    model: str = "azure_ai/gpt-5.4",
    max_workers: int = 10,
    step: int = 1000,
    start_offset: int = 0,
    min_new_variants: int = 8,
    patience: int = 2,
    triage_mode: str = "deterministic",
    logger: Optional[Any] = None,
) -> dict:
    run_dir = Path(run_dir)
    db = run_dir / f"{gene}.db"
    cand_tsv = run_dir / "extraction_priority" / "priority_candidates.tsv"
    try:
        max_cand = sum(1 for _ in cand_tsv.open()) - 1
    except Exception:
        max_cand = 0
    if max_cand <= 0:
        if logger:
            logger.warning(
                "full-coverage walk: no candidate ranking at %s; skipping",
                cand_tsv,
            )
        return {"walked": False, "reason": "no_candidates"}

    env = dict(
        os.environ,
        TIER3_MODELS=model,
        MAX_WORKERS=str(max_workers),
        AZURE_MAX_WORKERS=str(max_workers),
    )
    prev = _distinct_variants(db)
    seen_yield = False
    low = 0
    offset = max(0, start_offset)
    batches = 0
    batch_results: list[dict[str, Any]] = []
    while offset < max_cand:
        cmd = [
            sys.executable,
            str(PRIORITY),
            "--gene",
            gene,
            "--run-dir",
            str(run_dir),
            "--top-n",
            str(step),
            "--priority-offset",
            str(offset),
            "--max-workers",
            str(max_workers),
            "--skip-preprocess",
            "--triage-mode",
            triage_mode,
        ]
        rc = subprocess.run(cmd, cwd=str(ROOT), env=env).returncode
        batches += 1
        now = _distinct_variants(db)
        dv = now - prev
        batch_results.append(
            {
                "offset": offset,
                "top_n": step,
                "returncode": rc,
                "new_distinct_variants": dv,
                "distinct_variants": now,
            }
        )
        if logger:
            logger.info(
                "full-coverage walk: offset=%d rc=%d dVariants=%d cum=%d",
                offset,
                rc,
                dv,
                now,
            )
        if rc != 0:
            if logger:
                logger.warning(
                    "full-coverage walk: batch rc=%d at offset %d; stopping",
                    rc,
                    offset,
                )
            summary = {
                "walked": True,
                "batches": batches,
                "variants": now,
                "stopped": "error",
                "start_offset": start_offset,
                "min_new_variants": min_new_variants,
                "batch_results": batch_results,
            }
            _write_summary(run_dir, gene, summary)
            return summary
        if dv >= min_new_variants:
            seen_yield = True
            low = 0
        elif seen_yield:
            low += 1
        if seen_yield and low >= patience:
            if logger:
                logger.info(
                    "full-coverage walk: taper reached at offset=%d (%d batches)",
                    offset,
                    batches,
                )
            summary = {
                "walked": True,
                "batches": batches,
                "variants": now,
                "stopped": "taper",
                "start_offset": start_offset,
                "min_new_variants": min_new_variants,
                "batch_results": batch_results,
            }
            _write_summary(run_dir, gene, summary)
            return summary
        prev = now
        offset += step
    summary = {
        "walked": True,
        "batches": batches,
        "variants": _distinct_variants(db),
        "stopped": "candidate_cap",
        "start_offset": start_offset,
        "min_new_variants": min_new_variants,
        "batch_results": batch_results,
    }
    _write_summary(run_dir, gene, summary)
    return summary


def _write_summary(run_dir: Path, gene: str, summary: dict[str, Any]) -> None:
    outdir = run_dir / "extraction_priority"
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "full_coverage_walk_summary.json").write_text(
        json.dumps({"gene": gene, **summary}, indent=2),
        encoding="utf-8",
    )
