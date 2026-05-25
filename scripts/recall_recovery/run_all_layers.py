#!/usr/bin/env python3
"""One-shot driver: apply every gene-agnostic recall recovery layer in order.

Layers run from DB-observed PMIDs by default. If ``--gold`` is supplied, each
layer is scored so per-layer yield is visible; without gold, recovery still
runs and writes progression/QC artifacts.

  0. baseline    — score the DB as-is (after main extraction + migration)
  1. clinvar     — ingest ClinVar variants citing PMIDs already in the DB
  2. pubtator    — ingest NCBI PubTator3 mutations for PMIDs already in the DB
  3. figures     — vision-LLM figure reader over every PMID with images on disk
  4. v12_merge   — OPT-IN ONLY (``--with-v12 PATH``): copy variants from an
                   older manually-curated DB. Refuses by default because the
                   v12 baseline was only ever produced for KCNH2 and is NOT
                   cold-start capability.

Usage::

    # KCNH2 — opt in to the v12 merge because we have one on disk
    python scripts/recall_recovery/run_all_layers.py \\
        --gene KCNH2 \\
        --db results/KCNH2/20260517_074737/KCNH2.db \\
        --gold gene_variant_fetcher_gold_standard/normalized/KCNH2_recall_input.csv \\
        --pmc-dir results/KCNH2/20260517_074737/pmc_fulltext \\
        --with-v12 results/KCNH2/20260506_102238/end_to_end_20260515_manual_recovery/KCNH2_v12_manual_recovery_20260515.db

    # Cold-start/no-gold gene — no v12 path, just runs gene-agnostic layers
    python scripts/recall_recovery/run_all_layers.py \\
        --gene MYGENE \\
        --db results/MYGENE/<timestamp>/MYGENE.db \\
        --pmc-dir results/MYGENE/<timestamp>/pmc_fulltext
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]

logger = logging.getLogger("run_all_layers")


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def score_layer(
    gene: str,
    db: Path,
    gold: Path,
    outdir: Path,
) -> dict:
    """Run the recall suite and return the aggregate-recall dict."""
    outdir.mkdir(parents=True, exist_ok=True)
    # The recall suite's --gold-dir is the PACKAGE ROOT (parent of normalized/),
    # not the normalized dir itself.
    if gold.parent.name == "normalized":
        gold_dir = gold.parent.parent
    else:
        gold_dir = gold.parent
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "run_recall_suite.py"),
        "--score",
        "--genes",
        gene,
        "--db",
        f"{gene}={db}",
        "--gold-dir",
        str(gold_dir),
        "--outdir",
        str(outdir),
    ]
    logger.info("score → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("run_recall_suite.py failed: %s", result.stderr[-400:])
        return {}
    summary_path = outdir / "summary.json"
    if not summary_path.exists():
        logger.error("no summary.json under %s", outdir)
        return {}
    data = json.loads(summary_path.read_text())
    agg = data.get("aggregate_recall", {})
    return {k: v.get("recall") for k, v in agg.items() if isinstance(v, dict)}


# ---------------------------------------------------------------------------
# Layer runners
# ---------------------------------------------------------------------------


def run_clinvar(gene: str, db: Path, gold: Optional[Path], pmid_source: str) -> dict:
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "recall_recovery" / "ingest_clinvar.py"),
        "--gene",
        gene,
        "--db",
        str(db),
        "--pmid-source",
        pmid_source,
    ]
    if gold:
        cmd.extend(["--gold", str(gold)])
    logger.info("clinvar → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    logger.debug("stdout: %s", result.stdout)
    if result.returncode != 0:
        logger.error("ingest_clinvar failed: %s", result.stderr[-400:])
    try:
        return json.loads(result.stdout.strip().splitlines()[-1])
    except Exception:
        return {"added": None}


def run_pubtator(gene: str, db: Path, gold: Optional[Path], pmid_source: str) -> dict:
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "recall_recovery" / "ingest_pubtator.py"),
        "--gene",
        gene,
        "--db",
        str(db),
        "--pmid-source",
        pmid_source,
    ]
    if gold:
        cmd.extend(["--gold", str(gold)])
    logger.info("pubtator → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    logger.debug("stdout: %s", result.stdout)
    if result.returncode != 0:
        logger.error("ingest_pubtator failed: %s", result.stderr[-400:])
    try:
        return json.loads(result.stdout.strip().splitlines()[-1])
    except Exception:
        return {"added": None}


def run_figures(gene: str, db: Path, pmc_dir: Optional[Path], out: Path) -> dict:
    if not pmc_dir or not pmc_dir.is_dir():
        logger.warning("Skipping figures layer — pmc_dir missing: %s", pmc_dir)
        return {"added": None, "skipped": True}
    out.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "extract_figure_variants.py"),
        "--gene",
        gene,
        "--pmc-dir",
        str(pmc_dir),
        "--auto-pmids",
        "--out",
        str(out),
        "--db",
        str(db),
    ]
    logger.info("figures → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    logger.debug("stdout (tail): %s", result.stdout[-600:])
    if result.returncode != 0:
        logger.error("extract_figure_variants failed: %s", result.stderr[-400:])
    summary_path = out / "summary.json"
    if summary_path.exists():
        return json.loads(summary_path.read_text())
    return {"added": None}


def run_v12_merge(gene: str, db: Path, gold: Path, v12_db: Path) -> dict:
    if gene.upper() != "KCNH2":
        logger.error(
            "v12 merge refused for %s — the v12 manually-curated baseline only "
            "exists for KCNH2. Counting it as recall for a different gene would "
            "be mis-leading.",
            gene,
        )
        return {"refused": True}
    if not v12_db.exists():
        logger.error("v12 db not found: %s", v12_db)
        return {"refused": True}
    # Build the gold pmids file the merger expects
    gold_pmids = sorted(
        {r["pmid"] for r in csv.DictReader(open(gold)) if r.get("pmid")}
    )
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tmp:
        tmp.write("\n".join(gold_pmids))
        gold_pmids_path = tmp.name
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "recall_recovery" / "merge_v12_db.py"),
        "--v12-db",
        str(v12_db),
        "--target-db",
        str(db),
        "--gold-pmids",
        gold_pmids_path,
    ]
    logger.info("v12_merge → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    Path(gold_pmids_path).unlink(missing_ok=True)
    if result.returncode != 0:
        logger.error("merge_v12_db failed: %s", result.stderr[-400:])
    return {"stdout_tail": result.stdout[-400:]}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--gene", required=True, help="Gene symbol (e.g. KCNH2, KCNQ1)")
    p.add_argument("--db", required=True, help="Path to the recall SQLite DB")
    p.add_argument(
        "--gold",
        default=None,
        help=(
            "Path to the gold-standard recall CSV. Optional in DB-PMID mode; "
            "when omitted, recovery layers run but per-layer recall scoring is skipped."
        ),
    )
    p.add_argument(
        "--pmc-dir",
        default=None,
        help="pmc_fulltext directory for the figure-reader layer. If omitted, "
        "the figures layer is skipped.",
    )
    p.add_argument(
        "--with-v12",
        default=None,
        help="Path to a v12-style manually-curated DB to merge in. KCNH2-only "
        "and rejected for any other gene. Default: not applied.",
    )
    p.add_argument(
        "--backup",
        action="store_true",
        help="Snapshot the input DB to <DB>.before_layers.db before mutating it.",
    )
    p.add_argument(
        "--outdir",
        default=None,
        help="Where to put per-layer recall reports + summary. Default: "
        "recall_metrics/run_all_layers_{gene}_{timestamp}/",
    )
    p.add_argument(
        "--skip",
        action="append",
        default=[],
        help="Layer to skip: clinvar, pubtator, figures",
    )
    p.add_argument(
        "--allow-gold-pmid-enrichment",
        action="store_true",
        help=(
            "Use gold PMIDs as the ClinVar/PubTator query set. This is "
            "evaluation-aided and should not be used for cold-start recall claims."
        ),
    )
    p.add_argument(
        "--no-score",
        action="store_true",
        help="Skip per-layer recall scoring even when --gold is supplied.",
    )
    p.add_argument("--verbose", "-v", action="store_true")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    gene = args.gene.upper()
    db = Path(args.db)
    gold = Path(args.gold) if args.gold else None
    pmc_dir = Path(args.pmc_dir) if args.pmc_dir else None
    v12 = Path(args.with_v12) if args.with_v12 else None
    skip = {s.lower() for s in args.skip}
    pmid_source = "gold" if args.allow_gold_pmid_enrichment else "db"

    if not db.exists():
        sys.exit(f"DB not found: {db}")
    if args.allow_gold_pmid_enrichment and not gold:
        sys.exit("--gold is required with --allow-gold-pmid-enrichment")
    if gold and not gold.exists():
        sys.exit(f"Gold CSV not found: {gold}")
    scoring_enabled = bool(gold and not args.no_score)

    if args.backup:
        backup = db.with_suffix(db.suffix + ".before_layers.db")
        shutil.copy2(db, backup)
        logger.info("Backed up DB → %s", backup)
    else:
        logger.warning(
            "Mutating %s in place without a backup. Pass --backup for turnkey runs.",
            db,
        )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = (
        Path(args.outdir)
        if args.outdir
        else (REPO_ROOT / "recall_metrics" / f"run_all_layers_{gene}_{timestamp}")
    )
    outdir.mkdir(parents=True, exist_ok=True)

    progression: list[dict] = []

    def record(label: str, extra: Optional[dict] = None) -> None:
        recall = (
            score_layer(gene, db, gold, outdir / f"score_{label}")
            if scoring_enabled and gold
            else {}
        )
        row = {
            "layer": label,
            **{
                k: round(v * 100, 2) if v is not None else None
                for k, v in recall.items()
            },
        }
        if extra:
            row.update(extra)
        progression.append(row)
        logger.info(
            "[%s] %s",
            label,
            {
                k: row[k]
                for k in ("pmids", "variant_rows", "unique_variants")
                if k in row
            },
        )

    record("0_baseline")

    if "clinvar" not in skip:
        info = run_clinvar(gene, db, gold, pmid_source)
        record("1_clinvar", extra={"clinvar_added": info.get("added")})
    if "pubtator" not in skip:
        info = run_pubtator(gene, db, gold, pmid_source)
        record("2_pubtator", extra={"pubtator_added": info.get("added")})
    if "figures" not in skip:
        info = run_figures(gene, db, pmc_dir, outdir / "figure_reads")
        record("3_figures", extra={"figures_added": info.get("total_db_added")})
    if v12:
        if not gold:
            logger.error("Skipping v12 merge because --gold is required for gold PMIDs")
        else:
            info = run_v12_merge(gene, db, gold, v12)
            record("4_v12_merge", extra={"v12": str(v12)})

    summary = {
        "gene": gene,
        "db": str(db),
        "gold": str(gold) if gold else None,
        "pmc_dir": str(pmc_dir) if pmc_dir else None,
        "v12_db": str(v12) if v12 else None,
        "skipped": sorted(skip),
        "enrichment_pmid_source": pmid_source,
        "scoring_enabled": scoring_enabled,
        "progression": progression,
    }
    (outdir / "progression.json").write_text(json.dumps(summary, indent=2))

    # Also dump a CSV
    with (outdir / "progression.csv").open("w", newline="") as f:
        fieldnames = [
            "layer",
            "pmids",
            "variant_rows",
            "unique_variants",
            "patients",
            "affected",
            "unaffected",
            "clinvar_added",
            "pubtator_added",
            "figures_added",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in progression:
            w.writerow(row)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
