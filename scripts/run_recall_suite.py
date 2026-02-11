#!/usr/bin/env python3
"""
GVF Recall Test Suite Runner

Runs all recall tests and optionally logs metrics for tracking progress.

Usage:
    python scripts/run_recall_suite.py           # Run tests only
    python scripts/run_recall_suite.py --log     # Run and log metrics
    python scripts/run_recall_suite.py --score   # Show current score only
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# Paths
BASE_DIR = Path(__file__).parent.parent
GOLD_STANDARD = (
    BASE_DIR
    / "comparison_results"
    / "KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
)
DOWNLOAD_DIR = Path("/mnt/temp2/kronckbm/gvf_output/verified_downloads_20260208")
METRICS_DIR = BASE_DIR / "recall_metrics"
HISTORY_FILE = METRICS_DIR / "history.jsonl"


def calculate_metrics():
    """Calculate all recall metrics."""
    import glob
    import re

    # Load gold standard
    df = pd.read_excel(GOLD_STANDARD)
    df["PMID_int"] = df["PMID"].apply(lambda x: int(float(x)) if pd.notna(x) else None)
    df["is_LQT_affected"] = df["LQT"].apply(lambda x: x >= 1 if pd.notna(x) else False)

    gold_pmids = set(df["PMID_int"].dropna().astype(int))
    gold_variants = set(df["Variant"].dropna().apply(lambda v: str(v).upper().strip()))

    # Get downloaded PMIDs
    pdfs = glob.glob(str(DOWNLOAD_DIR / "*.pdf"))
    downloaded_pmids = set()
    for pdf in pdfs:
        match = re.match(r"^(\d+)", Path(pdf).name)
        if match:
            downloaded_pmids.add(int(match.group(1)))

    # Paper recall
    paper_overlap = gold_pmids & downloaded_pmids
    paper_recall = len(paper_overlap) / len(gold_pmids) if gold_pmids else 0

    # Variant coverage (by paper - which variants are in downloaded papers)
    variants_in_downloaded = set()
    for _, row in df.iterrows():
        pmid = row["PMID_int"]
        variant = row["Variant"]
        if pd.notna(variant) and pmid in downloaded_pmids:
            variants_in_downloaded.add(str(variant).upper().strip())

    variant_recall = (
        len(variants_in_downloaded) / len(gold_variants) if gold_variants else 0
    )

    # Carrier recall
    total_carriers = len(df)
    carriers_captured = len(df[df["PMID_int"].isin(downloaded_pmids)])
    carrier_recall = carriers_captured / total_carriers if total_carriers else 0

    # LQT recall
    total_lqt = df["is_LQT_affected"].sum()
    lqt_captured = df[df["PMID_int"].isin(downloaded_pmids)]["is_LQT_affected"].sum()
    lqt_recall = lqt_captured / total_lqt if total_lqt else 0

    return {
        "paper_recall": paper_recall,
        "variant_recall": variant_recall,
        "carrier_recall": carrier_recall,
        "lqt_recall": lqt_recall,
        "paper_count": len(paper_overlap),
        "paper_total": len(gold_pmids),
        "variant_count": len(variants_in_downloaded),
        "variant_total": len(gold_variants),
        "carrier_count": carriers_captured,
        "carrier_total": total_carriers,
        "lqt_count": int(lqt_captured),
        "lqt_total": int(total_lqt),
    }


def calculate_score(metrics: dict) -> float:
    """Calculate composite GVF recall score (0-100)."""
    weights = {
        "variant_recall": 0.40,
        "carrier_recall": 0.25,
        "lqt_recall": 0.25,
        "paper_recall": 0.10,
    }

    score = sum(metrics.get(k, 0) * weights[k] * 100 for k in weights)

    return round(score, 1)


def get_grade(score: float) -> str:
    """Convert score to letter grade."""
    if score >= 80:
        return "A"
    elif score >= 60:
        return "B"
    elif score >= 40:
        return "C"
    elif score >= 20:
        return "D"
    else:
        return "F"


def log_metrics(metrics: dict, notes: str = ""):
    """Append metrics to history."""
    METRICS_DIR.mkdir(exist_ok=True)

    entry = {
        "timestamp": datetime.now().isoformat(),
        "metrics": metrics,
        "score": calculate_score(metrics),
        "notes": notes,
    }

    with HISTORY_FILE.open("a") as f:
        f.write(json.dumps(entry) + "\n")

    print(f"Metrics logged to {HISTORY_FILE}")


def show_scorecard(metrics: dict):
    """Display formatted scorecard."""
    score = calculate_score(metrics)
    grade = get_grade(score)

    def bar(pct, width=20):
        filled = int(pct * width)
        return "█" * filled + "░" * (width - filled)

    print()
    print("╔" + "═" * 63 + "╗")
    print("║" + "GVF RECALL SCORECARD".center(63) + "║")
    print("║" + f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}".center(63) + "║")
    print("╠" + "═" * 63 + "╣")
    print(
        f"║  COMPOSITE SCORE: {score:.1f} / 100  [{bar(score / 100)}] Grade {grade}  ║"
    )
    print("╠" + "═" * 63 + "╣")
    print(
        f"║  Paper Recall:    {metrics['paper_recall'] * 100:5.1f}%  [{bar(metrics['paper_recall'])}] {metrics['paper_count']}/{metrics['paper_total']}".ljust(
            64
        )
        + "║"
    )
    print(
        f"║  Variant Recall:  {metrics['variant_recall'] * 100:5.1f}%  [{bar(metrics['variant_recall'])}] {metrics['variant_count']}/{metrics['variant_total']}".ljust(
            64
        )
        + "║"
    )
    print(
        f"║  Carrier Recall:  {metrics['carrier_recall'] * 100:5.1f}%  [{bar(metrics['carrier_recall'])}] {metrics['carrier_count']}/{metrics['carrier_total']}".ljust(
            64
        )
        + "║"
    )
    print(
        f"║  LQT Recall:      {metrics['lqt_recall'] * 100:5.1f}%  [{bar(metrics['lqt_recall'])}] {metrics['lqt_count']}/{metrics['lqt_total']}".ljust(
            64
        )
        + "║"
    )
    print("╠" + "═" * 63 + "╣")
    print(f"║  Target: 80.0  |  Gap: {80.0 - score:.1f}".ljust(64) + "║")
    print("╚" + "═" * 63 + "╝")
    print()


def run_tests():
    """Run pytest recall tests."""
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/recall/", "-v", "--tb=short"],
        cwd=str(BASE_DIR),
    )
    return result.returncode


def main():
    parser = argparse.ArgumentParser(description="GVF Recall Test Suite")
    parser.add_argument("--log", action="store_true", help="Log metrics to history")
    parser.add_argument(
        "--score", action="store_true", help="Show score only (skip tests)"
    )
    parser.add_argument("--notes", default="", help="Notes for log entry")
    args = parser.parse_args()

    # Calculate metrics
    print("Calculating metrics...")
    metrics = calculate_metrics()

    # Show scorecard
    show_scorecard(metrics)

    if args.score:
        return 0

    # Run tests
    print("Running recall tests...\n")
    exit_code = run_tests()

    # Log if requested
    if args.log:
        log_metrics(metrics, args.notes)

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
