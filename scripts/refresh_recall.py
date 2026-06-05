#!/usr/bin/env python3
"""Idempotent recall-refresh orchestrator — re-run as papers/permissions grow.

Composes the already-idempotent pieces into one command so the recall pipeline
can simply be re-executed over time. Nothing here mutates anything unless there
is genuinely new source to fold in; an unchanged corpus is a no-op.

Steps (per gene):
  1. Fetch publisher supplements that the full-text-API route drops. Elsevier
     ``mmc`` today (``scripts/fetch_elsevier_supplements.py``); Springer/Wiley
     activate automatically once their keys/access are restored (the script
     no-ops cleanly when a key is absent). Skips supplements already on disk.
  2. Bridge the nested corpus to the flat harvest layout
     (``scripts/corpus_to_harvest.py``) the refresh path expects.
  3. ``scripts/refresh_run_db.py`` folds on-disk supplements into FULL_CONTEXT and
     acceptance-gated-re-extracts any paper whose folded source now yields more
     variants (candidate selection by deterministic lift — no hand-picked PMIDs).
     Re-running after a successful pass selects nothing new.
  4. Score the rebuilt DB vs gold and print the delta.
  5. ``--land``: surgically promote the re-extracted papers into the CANONICAL DB
     (delete each changed paper's extraction-origin rows, preserving
     clinvar/pubtator/figure layer rows, then re-migrate), back up first, and
     only keep it if nothing regressed.

Why this is safe to re-run:
  * supplement fetch skips files already present;
  * the fold is sentinel-delimited + idempotent;
  * refresh re-extracts only papers with a real deterministic lift, behind the
    regression + explosion acceptance gates;
  * ``--land`` backs up the canonical DB and reverts on any regression.

  # measure only (no canonical mutation):
  python scripts/refresh_recall.py --gene SCN5A \
      --run-dir validation_runs/.../results/SCN5A/...
  # measure AND promote into the canonical DB if it improves:
  python scripts/refresh_recall.py --gene SCN5A --run-dir <run> --land
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sqlite3
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable

# variant_papers rows written by the recovery layers (NOT extraction) — preserved
# across a surgical land so re-extraction never drops ClinVar/PubTator/figure rows.
_LAYER = (
    "(LOWER(COALESCE(source_location,'')) LIKE '%clinvar%' "
    "OR LOWER(COALESCE(source_location,'')) LIKE '%pubtator%' "
    "OR LOWER(COALESCE(source_location,'')) LIKE '%figure%')"
)


def _run(cmd: list[str], label: str) -> None:
    print(f"\n=== {label} ===\n$ {' '.join(str(c) for c in cmd)}", flush=True)
    subprocess.run(cmd, cwd=str(REPO), check=True)


def _score(gene: str, db: Path, gold_dir: Path, outdir: Path) -> dict | None:
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        PY,
        str(REPO / "scripts" / "run_recall_suite.py"),
        "--score",
        "--genes",
        gene,
        "--db",
        f"{gene}={db}",
        "--outdir",
        str(outdir),
    ]
    subprocess.run(cmd, cwd=str(REPO), check=True)
    summ = outdir / gene / "summary.json"
    try:
        return json.loads(summ.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _uniqv(summary: dict | None) -> tuple[int, int]:
    """(matched, gold) for unique variants from a run_recall_suite summary."""
    if not summary:
        return (0, 0)
    for key in ("unique_variants", "unique_variant_recall", "uniqueVariants"):
        v = summary.get(key)
        if isinstance(v, dict):
            return int(v.get("matched", 0)), int(v.get("gold", v.get("total", 0)))
    return (0, 0)


def _land(
    gene: str, canonical_db: Path, staged_dir: Path, changed_pmids: list[str]
) -> Path:
    """Surgically inject re-extracted papers into a copy of the canonical DB.

    Preserves layer (clinvar/pubtator/figure) variant_papers rows; replaces each
    paper's extraction-origin penetrance/individual/variant_papers with the fresh
    (supplement-folded) re-extraction. Returns the path to the injected copy.
    """
    from harvesting.migrate_to_sqlite import (
        migrate_extraction_file,
        create_database_schema,
    )

    out = canonical_db.parent / f"{gene}.supp_landed.db"
    shutil.copy2(canonical_db, out)
    con = create_database_schema(str(out))  # idempotent + ALTERs any missing columns
    cur = con.cursor()
    for pm in changed_pmids:
        js = staged_dir / f"{gene}_PMID_{pm}.json"
        if not js.exists():
            continue
        cur.execute(
            "DELETE FROM age_dependent_penetrance WHERE penetrance_id IN "
            "(SELECT penetrance_id FROM penetrance_data WHERE pmid=?)",
            (pm,),
        )
        cur.execute("DELETE FROM penetrance_data WHERE pmid=?", (pm,))
        cur.execute("DELETE FROM individual_records WHERE pmid=?", (pm,))
        cur.execute(f"DELETE FROM variant_papers WHERE pmid=? AND NOT {_LAYER}", (pm,))
        migrate_extraction_file(cur, js, replace_existing_paper=True)
    con.commit()
    con.close()
    return out


def _changed_pmids(refresh_summary: Path) -> list[str]:
    """Successful replay PMIDs from a refresh_run_db refresh_summary.json."""
    try:
        data = json.loads(refresh_summary.read_text())
    except (OSError, json.JSONDecodeError):
        return []
    replay = data.get("replay") or data.get("replay_stats") or {}
    return list(replay.get("successful_pmids") or [])


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--gene", required=True)
    ap.add_argument(
        "--run-dir",
        required=True,
        type=Path,
        help="Canonical run dir (has extractions/ + <GENE>.db).",
    )
    ap.add_argument("--corpus", default=str(REPO / "corpus"))
    ap.add_argument(
        "--gold-dir", default=str(REPO / "gene_variant_fetcher_gold_standard")
    )
    ap.add_argument(
        "--workdir",
        default="/tmp/gvf_refresh_recall",
        help="Scratch dir for bridge/scoring.",
    )
    ap.add_argument(
        "--skip-fetch", action="store_true", help="Skip the supplement-fetch step."
    )
    ap.add_argument(
        "--land",
        action="store_true",
        help="Promote into the canonical DB if it improves (backs up first).",
    )
    args = ap.parse_args()

    gene = args.gene.upper()
    run_dir = args.run_dir.expanduser().resolve()
    canonical_db = run_dir / f"{gene}.db"
    work = Path(args.workdir) / gene
    work.mkdir(parents=True, exist_ok=True)
    gold_dir = Path(args.gold_dir)

    # 1. supplement fetch (idempotent; Elsevier today, Springer/Wiley when keyed)
    if not args.skip_fetch and os.getenv("ELSEVIER_API_KEY"):
        _run(
            [
                PY,
                str(REPO / "scripts" / "fetch_elsevier_supplements.py"),
                "--gene",
                gene,
            ],
            "1. fetch Elsevier supplements (idempotent)",
        )
    elif not args.skip_fetch:
        print("\n=== 1. supplement fetch skipped: ELSEVIER_API_KEY not set ===")

    # 2. bridge corpus -> flat harvest
    harvest = work / "harvest"
    _run(
        [
            PY,
            str(REPO / "scripts" / "corpus_to_harvest.py"),
            "--gene",
            gene,
            "--out",
            str(harvest),
        ],
        "2. bridge corpus -> flat harvest",
    )

    # 3. refresh: fold + acceptance-gated re-extract of grown papers (no forced PMIDs)
    out_db = work / f"{gene}.refreshed.db"
    if out_db.exists():
        out_db.unlink()
    _run(
        [
            PY,
            str(REPO / "scripts" / "refresh_run_db.py"),
            "--gene",
            gene,
            "--run-dir",
            str(run_dir),
            "--harvest-dir",
            str(harvest),
            "--stage-extractions",
            "--output-db",
            str(out_db),
            "--skip-recovery",
        ],
        "3. refresh_run_db (fold + gated re-extract of grown papers)",
    )

    # locate the refresh dir this produced (newest)
    refresh_dirs = sorted(run_dir.glob("refresh_*"), key=lambda p: p.stat().st_mtime)
    refresh_dir = refresh_dirs[-1] if refresh_dirs else None
    changed = (
        _changed_pmids(refresh_dir / "refresh_summary.json") if refresh_dir else []
    )
    print(f"\nRe-extracted (accepted) papers this pass: {len(changed)} {changed}")

    # 4. score the refreshed DB
    before = _score(gene, canonical_db, gold_dir, work / "score_canonical")
    after = _score(gene, out_db, gold_dir, work / "score_refreshed")
    bm, bg = _uniqv(before)
    am, ag = _uniqv(after)
    print(
        f"\nunique_variants: canonical {bm}/{bg} -> refreshed {am}/{ag}  (Δ={am - bm})"
    )

    if not changed:
        print(
            "\nNo new re-extractions — corpus unchanged since last refresh (idempotent no-op)."
        )
        return 0

    # 5. optional gated land into the canonical DB
    if args.land:
        if refresh_dir is None:
            print("Cannot land: no refresh dir found.")
            return 1
        landed = _land(gene, canonical_db, refresh_dir / "staged_extractions", changed)
        lm, lg = _uniqv(_score(gene, landed, gold_dir, work / "score_landed"))
        print(f"\nlanded DB unique_variants: {lm}/{lg} (canonical was {bm}/{bg})")
        if lm >= bm:
            backup = canonical_db.parent / f"{gene}.db.before_refresh_recall.db"
            shutil.copy2(canonical_db, backup)
            shutil.copy2(landed, canonical_db)
            print(
                f"PROMOTED -> {canonical_db} (backup: {backup.name}); Δ uniqV = {lm - bm}"
            )
        else:
            print(
                f"NOT promoted: landed {lm} < canonical {bm} (regression). Canonical untouched."
            )
    else:
        print("\n(measure-only; pass --land to promote into the canonical DB.)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
