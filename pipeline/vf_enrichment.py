"""variantFeatures enrichment + false-positive quarantine as a gvf-run step.

Composes two existing scripts (no pipeline edits):
  * ``scripts/enrich_from_variantfeatures.py`` — match GVF variants to the
    variantFeatures warehouse → canonical names + in silico scores + an fp_class
    flag, written to a ``vf_enrichment`` table.
  * ``scripts/quarantine_fp.py`` — move the high-confidence wrong-gene FPs
    (``misparse_out_of_range``: residue beyond the gene's protein length) into a
    ``quarantined_variants`` table and out of the active tables.

Turnkey-safe: a gene absent from variantFeatures yields an all-unmatched
enrichment with max protein length 0, so nothing is ever quarantined. Needs no
API keys (reads the local variantFeatures SQLite; pure-SQLite quarantine).
"""

from __future__ import annotations

import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

ROOT = Path(__file__).resolve().parent.parent


def enrich_and_quarantine(
    gene: str,
    db,
    *,
    vf_db: Optional[Path] = None,
    quarantine_classes: str = "misparse_out_of_range",
    stamp: Optional[str] = None,
    logger: Optional[Any] = None,
) -> dict:
    db = Path(db)
    enrich = ROOT / "scripts" / "enrich_from_variantfeatures.py"
    quar = ROOT / "scripts" / "quarantine_fp.py"
    if not enrich.exists():
        if logger:
            logger.warning("vf-enrich: %s missing; skipping", enrich)
        return {"enriched": False, "reason": "script_missing"}

    cmd = [sys.executable, str(enrich), "--gene", gene, "--db", str(db)]
    if vf_db:
        cmd += ["--vf", str(vf_db)]
    rc1 = subprocess.run(cmd, cwd=str(ROOT)).returncode
    if rc1 != 0:
        if logger:
            logger.warning("vf-enrich: enrichment rc=%d; skipping quarantine", rc1)
        return {"enriched": False, "rc": rc1}

    stamp = stamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    rc2 = 0
    if quar.exists():
        rc2 = subprocess.run(
            [
                sys.executable,
                str(quar),
                "--db",
                str(db),
                "--classes",
                quarantine_classes,
                "--stamp",
                stamp,
            ],
            cwd=str(ROOT),
        ).returncode
    if logger:
        logger.info(
            "vf-enrich: enrichment ok; quarantine rc=%d (classes=%s)",
            rc2,
            quarantine_classes,
        )
    return {"enriched": True, "quarantined": rc2 == 0}
