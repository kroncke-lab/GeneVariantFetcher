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
     only keep it if neither unique-variant recall NOR rows-mode count MAE
     (carriers/affected/unaffected) regressed.

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
    "(LOWER(COALESCE(source_layer,'')) IN ('clinvar', 'pubtator', 'figure') "
    "OR LOWER(COALESCE(source_location,'')) LIKE '%clinvar%' "
    "OR LOWER(COALESCE(source_location,'')) LIKE '%pubtator%' "
    "OR LOWER(COALESCE(source_location,'')) LIKE '%figure%')"
)


def _run(cmd: list[str], label: str) -> None:
    print(f"\n=== {label} ===\n$ {' '.join(str(c) for c in cmd)}", flush=True)
    subprocess.run(cmd, cwd=str(REPO), check=True)


def _supplement_fetch_cmd(gene: str, corpus: Path) -> list[str]:
    return [
        PY,
        str(REPO / "scripts" / "fetch_elsevier_supplements.py"),
        "--gene",
        gene,
        "--corpus",
        str(corpus),
    ]


def _corpus_bridge_cmd(gene: str, corpus: Path, harvest: Path) -> list[str]:
    return [
        PY,
        str(REPO / "scripts" / "corpus_to_harvest.py"),
        "--gene",
        gene,
        "--corpus",
        str(corpus),
        "--out",
        str(harvest),
    ]


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
    """(matched, gold) for unique variants from a run_recall_suite summary.

    run_recall_suite writes ``summary["recall"]["unique_variants"] = {matched,
    gold, recall}``; tolerate a flattened layout too.
    """
    if not summary:
        return (0, 0)
    rec = summary.get("recall") or {}
    v = rec.get("unique_variants") or summary.get("unique_variants")
    if isinstance(v, dict):
        return int(v.get("matched", 0)), int(v.get("gold", v.get("total", 0)))
    return (0, 0)


def _rows(summary: dict | None) -> tuple[int, int]:
    """(matched, gold) for variant rows from a run_recall_suite summary.

    Mirrors ``_uniqv`` for the ``variant_rows`` recall block. Used by the land
    gate so a promotion can't hold unique-variant recall while quietly dropping
    matched variant ROWS (the per-(pmid, variant) coverage the grant scores).
    """
    if not summary:
        return (0, 0)
    rec = summary.get("recall") or {}
    v = rec.get("variant_rows") or summary.get("variant_rows")
    if isinstance(v, dict):
        return int(v.get("matched", 0)), int(v.get("gold", v.get("total", 0)))
    return (0, 0)


# --- MAE non-regression guard ----------------------------------------------
# Promotion requires recall to hold AND rows-mode count accuracy not to regress.
# A re-extraction can lift unique-variant recall while silently inflating carrier
# counts (the 2026-06-06 KCNQ1 carriers MAE 0.882->1.274 spike: one variant
# repeated across supplement tables, summed per (pmid, variant)). Recall is
# necessary but not sufficient. Two failure modes are guarded per count field
# (carriers/affected/unaffected):
#   * MAE worsens beyond max(GVF_LAND_MAE_ABS_TOL, baseline*GVF_LAND_MAE_REL_TOL)
#     — tolerance absorbs re-extraction noise; and
#   * count COVERAGE collapses — a re-extraction must not "improve" MAE by
#     dropping the very counts it is scored on (n_matched falling away by more
#     than GVF_LAND_MAE_COV_TOL when the baseline had counts).
_MAE_FIELDS = ("carriers", "affected", "unaffected")


def _env_float(name: str, default: float) -> float:
    """Parse a float env override at import; fall back on absence/garbage.

    A malformed value (e.g. ``GVF_LAND_MAE_ABS_TOL=foo``) must not crash the
    module before argparse/--help can run.
    """
    raw = os.getenv(name)
    if not raw or not raw.strip():
        return default
    try:
        return float(raw)
    except ValueError:
        print(
            f"warning: {name}={raw!r} is not a float; using {default}",
            file=sys.stderr,
        )
        return default


_MAE_ABS_TOL = _env_float("GVF_LAND_MAE_ABS_TOL", 0.05)
_MAE_REL_TOL = _env_float("GVF_LAND_MAE_REL_TOL", 0.05)
_MAE_COV_TOL = _env_float("GVF_LAND_MAE_COV_TOL", 0.10)
# Per-field absolute MAE tolerance. Carriers (total observed) is the count most
# perturbed when a real supplement-table row legitimately adds carriers to a
# previously-undercounted variant, so it gets a modestly more generous absolute
# tolerance than affected/unaffected — recovering real rows should not be blocked
# by a small carrier-MAE wiggle. The coverage-collapse guard and the relative
# tolerance still apply, so a genuine count blow-up is still caught. Override per
# field via GVF_LAND_MAE_ABS_TOL_<CARRIERS|AFFECTED|UNAFFECTED>.
_MAE_ABS_TOL_BY_FIELD = {
    "carriers": _env_float("GVF_LAND_MAE_ABS_TOL_CARRIERS", 0.10),
    "affected": _env_float("GVF_LAND_MAE_ABS_TOL_AFFECTED", _MAE_ABS_TOL),
    "unaffected": _env_float("GVF_LAND_MAE_ABS_TOL_UNAFFECTED", _MAE_ABS_TOL),
}


def _mae(summary: dict | None, field: str = "carriers") -> float | None:
    """Rows-mode MAE for a count field from a run_recall_suite summary.

    run_recall_suite writes ``summary["mae"][field] = {sum_abs_error,
    n_matched, mae}``. Returns the MAE float, or None when unavailable (no
    matched rows, or an older summary without the block).
    """
    if not summary:
        return None
    block = (summary.get("mae") or {}).get(field) or {}
    v = block.get("mae")
    return float(v) if isinstance(v, (int, float)) else None


def _mae_n(summary: dict | None, field: str = "carriers") -> int:
    """Matched-row count contributing to a field's MAE (0 if unavailable)."""
    if not summary:
        return 0
    block = (summary.get("mae") or {}).get(field) or {}
    try:
        return int(block.get("n_matched") or 0)
    except (TypeError, ValueError):
        return 0


def _mae_regressions(before: dict | None, after: dict | None) -> list[str]:
    """Count fields whose MAE or count-coverage regresses before->after.

    Returns human-readable strings like ``"carriers 0.882->1.274"`` or
    ``"carriers count-coverage 120->0"``; an empty list means the check passes.
    The coverage check fires first so a re-extraction cannot evade the MAE bound
    by dropping the very counts it would be scored on.
    """
    out: list[str] = []
    for field in _MAE_FIELDS:
        bn, an = _mae_n(before, field), _mae_n(after, field)
        if bn > 0 and an < bn * (1 - _MAE_COV_TOL):
            out.append(f"{field} count-coverage {bn}->{an}")
            continue
        b, a = _mae(before, field), _mae(after, field)
        if b is None or a is None:
            continue
        abs_tol = _MAE_ABS_TOL_BY_FIELD.get(field, _MAE_ABS_TOL)
        if a > b + max(abs_tol, b * _MAE_REL_TOL):
            out.append(f"{field} {b:.3f}->{a:.3f}")
    return out


def _promotion_decision(
    lm: int, bm: int, mae_regressions: list[str], lr: int = 0, br: int = 0
) -> tuple[bool, str]:
    """Whether a landed DB may replace the canonical one, with the reason if not.

    Promote only if unique-variant recall holds (``lm >= bm``) AND variant-row
    recall holds (``lr >= br``) AND no count field regressed (``mae_regressions``
    empty). The row-recall check (Codex review guardrail) stops a promotion that
    lifts unique-variant recall while dropping matched variant ROWS; it is
    skipped when ``br`` is 0 (caller did not supply row recall). Returns
    ``(promote, message)``; ``message`` is the NOT-promoted explanation (empty
    when promoting). Pure so the land decision is unit-testable without the
    DB/scoring machinery.
    """
    if lm < bm:
        return False, (
            f"NOT promoted: landed {lm} < canonical {bm} "
            "(unique-variant regression). Canonical untouched."
        )
    if br and lr < br:
        return False, (
            f"NOT promoted: landed variant-row recall {lr} < canonical {br} "
            "(variant-row regression). Canonical untouched."
        )
    if mae_regressions:
        return False, (
            "NOT promoted: rows-mode count accuracy regressed ("
            + "; ".join(mae_regressions)
            + "). Recall held but per-row counts worsened — canonical untouched. "
            "Refine the re-extraction, or relax the bound via GVF_LAND_MAE_ABS_TOL"
            " / GVF_LAND_MAE_REL_TOL / GVF_LAND_MAE_COV_TOL."
        )
    return True, ""


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
            _supplement_fetch_cmd(gene, Path(args.corpus)),
            "1. fetch Elsevier supplements (idempotent)",
        )
    elif not args.skip_fetch:
        print("\n=== 1. supplement fetch skipped: ELSEVIER_API_KEY not set ===")

    # 2. bridge corpus -> flat harvest
    harvest = work / "harvest"
    _run(
        _corpus_bridge_cmd(gene, Path(args.corpus), harvest),
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
        landed_summary = _score(gene, landed, gold_dir, work / "score_landed")
        lm, lg = _uniqv(landed_summary)
        lr, _lrg = _rows(landed_summary)
        br, _brg = _rows(before)
        print(
            f"\nlanded DB unique_variants: {lm}/{lg} (canonical was {bm}/{bg}); "
            f"variant_rows: {lr} (canonical was {br})"
        )
        promote, msg = _promotion_decision(
            lm, bm, _mae_regressions(before, landed_summary), lr=lr, br=br
        )
        if promote:
            backup = canonical_db.parent / f"{gene}.db.before_refresh_recall.db"
            shutil.copy2(canonical_db, backup)
            shutil.copy2(landed, canonical_db)
            print(
                f"PROMOTED -> {canonical_db} (backup: {backup.name}); Δ uniqV = {lm - bm}"
            )
        else:
            print(msg)
    else:
        print("\n(measure-only; pass --land to promote into the canonical DB.)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
