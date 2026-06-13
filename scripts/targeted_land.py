#!/usr/bin/env python3
"""Fast targeted land — re-extract only the quality-lift candidate papers and
surgically inject them into the canonical DB, gold-gated.

``refresh_recall.py --land`` folds the whole (~2,400-paper) corpus on every run
(~1 hour) even when only a handful of papers changed. This lands the same
recovery in minutes by:

  1. SCAN (no LLM): for each paper in the run, compare the gene-scoped
     deterministic TABLE parse of its corpus source to the current extraction;
     pick papers where re-extraction would add cDNA+protein PAIRING (the
     over-counted-cDNA-only class). Same intrinsic, gold-free signal the
     ``refresh_run_db`` selector uses — no PMID/gene/gold literals.
  2. BRIDGE only those papers' sources to a small flat harvest.
  3. RE-EXTRACT + GATE only those papers (``refresh_run_db --only-forced-pmids
     --no-supplement-fold --stage-extractions``): the deterministic-baseline
     acceptance gate decides keep/revert per paper.
  4. SURGICALLY INJECT the accepted papers into a copy of the canonical DB
     (``refresh_recall._land``; preserves clinvar/pubtator/figure layer rows),
     score vs gold before/after, and (with ``--land``) promote only if unique
     recall, variant-row recall, AND rows-mode MAE all hold (backup first;
     revert otherwise).

Gold is used ONLY to gate the final promotion (step 4) and is optional; the
selection (steps 1-3) is gene-agnostic and gold-free, so this is the same
turnkey path a brand-new gene would use.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pipeline.extraction import ExpertExtractor  # noqa: E402
from scripts import refresh_recall as RR  # noqa: E402
from scripts import refresh_run_db as RDB  # noqa: E402
from scripts.corpus_to_harvest import build as bridge_corpus  # noqa: E402

PY = sys.executable

# A current extraction is only worth the (cheap, no-LLM) deterministic re-parse
# when it looks over-counted-but-unpaired: many variants, low pairing ratio.
# This bounds the scan to suspicious papers instead of every paper in the run.
_SUSPICIOUS_MIN_VARIANTS = 15
_SUSPICIOUS_MAX_PAIRED_RATIO = 0.4


def _current_variants(extraction_json: Path) -> list[dict]:
    try:
        data = json.loads(extraction_json.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return []
    variants = data.get("variants")
    return variants if isinstance(variants, list) else []


def scan_quality_lift_candidates(
    gene: str, run_dir: Path, corpus: Path, min_deterministic: int = 5
) -> list[str]:
    """PMIDs where the deterministic table parse adds pairing over the current
    extraction. Pre-filters on the current pairing ratio so the deterministic
    parse only runs on suspicious (over-counted-cDNA-only) papers."""
    extraction_dir = run_dir / "extractions"
    extractor = ExpertExtractor(models=["noop"], tier_threshold=0)
    candidates: list[str] = []
    for ej in sorted(extraction_dir.glob(f"{gene}_PMID_*.json")):
        pmid = ej.name.removeprefix(f"{gene}_PMID_").removesuffix(".json")
        cur = _current_variants(ej)
        cur_paired = RDB._paired_count(cur)
        # Pre-filter: skip papers that are already small or already well-paired.
        if len(cur) < _SUSPICIOUS_MIN_VARIANTS:
            continue
        if cur and cur_paired / len(cur) > _SUSPICIOUS_MAX_PAIRED_RATIO:
            continue
        src = corpus / gene / pmid / f"{pmid}_FULL_CONTEXT.md"
        if not src.exists():
            continue
        det = RDB.deterministic_variant_list(extractor, src, gene)
        if len(det) >= min_deterministic and RDB._paired_count(det) > cur_paired:
            candidates.append(pmid)
    return candidates


def _score_uniqv_rows(gene: str, db: Path, gold_dir: Path, outdir: Path):
    summary = RR._score(gene, db, gold_dir, outdir)
    return summary, RR._uniqv(summary), RR._rows(summary)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--gene", required=True)
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--corpus", type=Path, default=REPO_ROOT / "corpus")
    ap.add_argument("--gold-dir", type=Path, default=None)
    ap.add_argument(
        "--pmids-file",
        type=Path,
        default=None,
        help="Skip the scan; use these PMIDs as the candidates.",
    )
    ap.add_argument("--workdir", type=Path, default=None)
    ap.add_argument(
        "--land",
        action="store_true",
        help="Promote into the canonical DB if it improves (backs up first).",
    )
    args = ap.parse_args()

    gene = args.gene
    run_dir = args.run_dir.resolve()
    canonical_db = run_dir / f"{gene}.db"
    if not canonical_db.exists():
        print(f"No canonical DB at {canonical_db}")
        return 1
    gold_dir = args.gold_dir or (
        REPO_ROOT / "gene_variant_fetcher_gold_standard" / "normalized"
    )
    work = Path(args.workdir or tempfile.mkdtemp(prefix=f"gvf_targeted_{gene}_"))
    work.mkdir(parents=True, exist_ok=True)

    # 1. candidates
    if args.pmids_file:
        candidates = [
            ln.strip()
            for ln in args.pmids_file.read_text().splitlines()
            if ln.strip() and not ln.startswith("#")
        ]
        print(f"[targeted] using {len(candidates)} provided candidate PMIDs")
    else:
        print("[targeted] scanning for quality-lift candidates (no LLM)...")
        candidates = scan_quality_lift_candidates(gene, run_dir, args.corpus)
        print(f"[targeted] found {len(candidates)} quality-lift candidate(s)")
    if not candidates:
        print("[targeted] nothing to do.")
        return 0

    pmids_file = work / "candidates.txt"
    pmids_file.write_text("\n".join(candidates) + "\n", encoding="utf-8")

    # 2. bridge ONLY the candidates
    harvest = work / "harvest"
    harvest.mkdir(exist_ok=True)
    n_fc, _ = bridge_corpus(gene, args.corpus, harvest, set(candidates))
    print(f"[targeted] bridged {n_fc} candidate source(s) to {harvest}")

    # 3. re-extract + gate ONLY the candidates (no fold, staged)
    out_db = work / f"{gene}.targeted.db"
    cmd = [
        PY,
        str(REPO_ROOT / "scripts" / "refresh_run_db.py"),
        "--gene",
        gene,
        "--run-dir",
        str(run_dir),
        "--harvest-dir",
        str(harvest),
        "--pmids-file",
        str(pmids_file),
        "--only-forced-pmids",
        "--no-supplement-fold",
        "--stage-extractions",
        "--skip-recovery",
        "--output-db",
        str(out_db),
    ]
    print(f"[targeted] re-extract + gate {len(candidates)} paper(s)...")
    subprocess.run(cmd, cwd=str(REPO_ROOT), check=True)

    # locate the refresh dir this produced and its accepted PMIDs
    refresh_dirs = sorted(run_dir.glob("refresh_*"), key=lambda p: p.stat().st_mtime)
    if not refresh_dirs:
        print("[targeted] no refresh dir produced; aborting.")
        return 1
    refresh_dir = refresh_dirs[-1]
    accepted = RR._changed_pmids(refresh_dir / "refresh_summary.json")
    print(f"[targeted] gate accepted {len(accepted)} paper(s): {accepted}")
    if not accepted:
        print("[targeted] nothing passed the gate; canonical untouched.")
        return 0

    # 4. surgical inject + score + gated promote
    staged = refresh_dir / "staged_extractions"
    landed = RR._land(gene, canonical_db, staged, accepted)
    before, (bm, _bg), (br, _brg) = _score_uniqv_rows(
        gene, canonical_db, gold_dir, work / "score_before"
    )
    after, (lm, _lg), (lr, _lrg) = _score_uniqv_rows(
        gene, landed, gold_dir, work / "score_after"
    )
    print(
        f"\n[targeted] unique_variants {bm} -> {lm} | variant_rows {br} -> {lr} "
        f"| candidates accepted={len(accepted)}"
    )
    promote, msg = RR._promotion_decision(
        lm, bm, RR._mae_regressions(before, after), lr=lr, br=br
    )
    if not args.land:
        print("\n[targeted] (measure-only; pass --land to promote.)")
        print("PROMOTE" if promote else msg)
        return 0
    if promote:
        backup = canonical_db.parent / f"{gene}.db.before_targeted_land.db"
        import shutil

        shutil.copy2(canonical_db, backup)
        shutil.copy2(landed, canonical_db)
        print(f"\nPROMOTED -> {canonical_db} (backup: {backup.name})")
    else:
        print(f"\n{msg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
