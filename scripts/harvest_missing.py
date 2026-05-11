#!/usr/bin/env python3
"""
Harvest-only script for missing gold-standard PMIDs.

Runs ONLY the download/harvesting phase — no LLM calls, no token spend.
Then audits download coverage against the gold standard.

Usage:
    python scripts/harvest_missing.py [--run-dir results/KCNH2/TIMESTAMP]

To do a dry-run audit without downloading:
    python scripts/harvest_missing.py --audit-only
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))
os.chdir(project_root)

# Load .env before anything else
try:
    from dotenv import load_dotenv

    load_dotenv(project_root / ".env")
except ImportError:
    print("Warning: python-dotenv not installed, relying on shell env vars")


def load_gold_pmids() -> set:
    """Load gold standard PMIDs."""
    gold_file = project_root / "comparison_results" / "gold_standard_pmids.txt"
    pmids = set()
    with open(gold_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                pmids.add(line)
    return pmids


def find_latest_run_dir(gene: str = "KCNH2") -> Path:
    """Find the latest run directory for a gene."""
    results_dir = project_root / "results" / gene
    if not results_dir.exists():
        raise FileNotFoundError(f"No results directory: {results_dir}")
    runs = sorted(results_dir.iterdir(), reverse=True)
    for run in runs:
        if run.is_dir() and re.match(r"\d{8}_\d{6}", run.name):
            return run
    raise FileNotFoundError(f"No timestamped run dirs in {results_dir}")


def get_fulltext_pmids(harvest_dir: Path) -> set:
    """Get PMIDs that have FULL_CONTEXT.md files."""
    pmids = set()
    if not harvest_dir.exists():
        return pmids
    for f in harvest_dir.iterdir():
        if f.name.endswith("_FULL_CONTEXT.md"):
            m = re.match(r"(\d+)_FULL_CONTEXT\.md", f.name)
            if m:
                pmids.add(m.group(1))
    return pmids


def get_filter_decisions(run_dir: Path) -> dict:
    """Load filter decisions from the filter progress log."""
    filter_file = run_dir / "pmid_status" / "filter_progress.jsonl"
    decisions = {}
    if not filter_file.exists():
        return decisions
    with open(filter_file) as f:
        for line in f:
            try:
                rec = json.loads(line)
                pmid = str(rec.get("pmid", ""))
                decisions[pmid] = rec
            except Exception:
                pass
    return decisions


def audit(run_dir: Path, gold_pmids: set) -> dict:
    """Audit download coverage against gold standard."""
    harvest_dir = run_dir / "pmc_fulltext"
    fulltext_pmids = get_fulltext_pmids(harvest_dir)
    filter_decisions = get_filter_decisions(run_dir)

    # Categorize each gold PMID
    not_discovered = []
    failed_filter = []
    passed_no_fulltext = []
    has_fulltext = []

    for pmid in sorted(gold_pmids):
        if pmid not in filter_decisions:
            not_discovered.append(pmid)
        elif filter_decisions[pmid].get("final_decision") != "PASS":
            failed_filter.append(pmid)
        elif pmid not in fulltext_pmids:
            passed_no_fulltext.append(pmid)
        else:
            has_fulltext.append(pmid)

    # Check paywalled log for failure reasons
    paywalled_reasons = {}
    paywall_file = harvest_dir / "paywalled_missing.csv"
    if paywall_file.exists():
        import csv

        with open(paywall_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                paywalled_reasons[str(row.get("PMID", ""))] = row.get(
                    "Reason", "unknown"
                )

    results = {
        "total_gold": len(gold_pmids),
        "not_discovered": not_discovered,
        "failed_filter": failed_filter,
        "passed_no_fulltext": passed_no_fulltext,
        "has_fulltext": has_fulltext,
        "paywalled_reasons": paywalled_reasons,
        "total_fulltext_files": len(fulltext_pmids),
    }
    return results


def print_audit(results: dict):
    """Print audit results."""
    total = results["total_gold"]
    nd = len(results["not_discovered"])
    ff = len(results["failed_filter"])
    pnf = len(results["passed_no_fulltext"])
    hf = len(results["has_fulltext"])

    print("\n" + "=" * 65)
    print("GOLD STANDARD PIPELINE FUNNEL")
    print("=" * 65)
    print(f"  Total gold PMIDs:                       {total:4d}")
    print(f"  1. Not discovered (not in search):       {nd:4d}  ({nd/total*100:.1f}%)")
    print(f"  2. Failed Tier 1 filter:                 {ff:4d}  ({ff/total*100:.1f}%)")
    print(
        f"  3. Passed filter, NO full text:          {pnf:4d}  ({pnf/total*100:.1f}%)"
    )
    print(f"  4. Has full text (download success):     {hf:4d}  ({hf/total*100:.1f}%)")
    print(f"  ─────────────────────────────────────────────")
    print(f"  Download rate:  {hf}/{total} = {hf/total*100:.1f}%")
    print(
        f"  Total fulltext files (all PMIDs):        {results['total_fulltext_files']}"
    )

    if results["not_discovered"]:
        print(f"\n  NOT DISCOVERED ({nd} PMIDs — need manual addition):")
        for p in results["not_discovered"]:
            print(f"    PMID {p}")

    if results["failed_filter"]:
        print(f"\n  FAILED FILTER ({ff} PMIDs — wrongly rejected):")
        for p in results["failed_filter"]:
            print(f"    PMID {p}")

    # Categorize missing-fulltext by reason
    paywalled = []
    no_reason = []
    for p in results["passed_no_fulltext"]:
        if p in results["paywalled_reasons"]:
            paywalled.append((p, results["paywalled_reasons"][p]))
        else:
            no_reason.append(p)

    if paywalled:
        print(f"\n  PAYWALLED / BLOCKED ({len(paywalled)} PMIDs):")
        # Group by reason
        by_reason = {}
        for pmid, reason in paywalled:
            short = reason[:80]
            by_reason.setdefault(short, []).append(pmid)
        for reason, pmids in sorted(by_reason.items(), key=lambda x: -len(x[1])):
            print(f"    [{len(pmids)} PMIDs] {reason}")
            if len(pmids) <= 5:
                for p in pmids:
                    print(f"      PMID {p}")

    if no_reason:
        print(f"\n  NEVER ATTEMPTED ({len(no_reason)} PMIDs — download cap victims):")
        print(f"    (These should download on next uncapped run)")
        for p in no_reason[:20]:
            print(f"    PMID {p}")
        if len(no_reason) > 20:
            print(f"    ... and {len(no_reason)-20} more")

    return results


def harvest_missing(run_dir: Path, missing_pmids: list, gene: str = "KCNH2"):
    """Download full text for missing PMIDs using the existing harvester."""
    from harvesting import PMCHarvester

    harvest_dir = run_dir / "pmc_fulltext"
    harvest_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'=' * 65}")
    print(f"HARVESTING {len(missing_pmids)} MISSING PMIDs")
    print(f"{'=' * 65}")
    print(f"Output: {harvest_dir}")

    harvester = PMCHarvester(output_dir=str(harvest_dir), gene_symbol=gene)

    # Harvest in batches to avoid overwhelming APIs
    batch_size = 50
    total_success = 0
    total_failed = 0

    for i in range(0, len(missing_pmids), batch_size):
        batch = missing_pmids[i : i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(missing_pmids) + batch_size - 1) // batch_size
        print(f"\n  Batch {batch_num}/{total_batches}: {len(batch)} PMIDs...")

        try:
            harvester.harvest(batch, delay=2.0)
        except Exception as e:
            print(f"  ERROR in batch {batch_num}: {e}")

    # Count results
    new_fulltext = get_fulltext_pmids(harvest_dir)
    newly_downloaded = set(missing_pmids) & new_fulltext
    still_missing = set(missing_pmids) - new_fulltext

    print(f"\n  Results:")
    print(f"    Successfully downloaded: {len(newly_downloaded)}")
    print(f"    Still missing:           {len(still_missing)}")

    return newly_downloaded, still_missing


def main():
    parser = argparse.ArgumentParser(description="Harvest missing gold-standard PMIDs")
    parser.add_argument(
        "--run-dir", type=str, help="Run directory (auto-detected if omitted)"
    )
    parser.add_argument("--gene", default="KCNH2", help="Gene symbol (default: KCNH2)")
    parser.add_argument(
        "--audit-only", action="store_true", help="Only audit, don't download"
    )
    parser.add_argument(
        "--include-undiscovered",
        action="store_true",
        help="Also try to harvest PMIDs not found by search",
    )
    args = parser.parse_args()

    # Find run directory
    if args.run_dir:
        run_dir = Path(args.run_dir)
    else:
        run_dir = find_latest_run_dir(args.gene)
    print(f"Using run directory: {run_dir}")

    # Load gold standard
    gold_pmids = load_gold_pmids()
    print(f"Gold standard PMIDs: {len(gold_pmids)}")

    # Initial audit
    results = audit(run_dir, gold_pmids)
    print_audit(results)

    if args.audit_only:
        # Save audit report
        report_path = run_dir / "harvest_audit.json"
        report = {
            "total_gold": results["total_gold"],
            "not_discovered": results["not_discovered"],
            "failed_filter": results["failed_filter"],
            "passed_no_fulltext": results["passed_no_fulltext"],
            "has_fulltext_count": len(results["has_fulltext"]),
            "download_rate": f"{len(results['has_fulltext'])/results['total_gold']*100:.1f}%",
        }
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)
        print(f"\nAudit report saved to: {report_path}")
        return

    # Determine what to harvest
    to_harvest = list(results["passed_no_fulltext"])  # Cap victims
    if args.include_undiscovered:
        to_harvest.extend(results["not_discovered"])

    if not to_harvest:
        print("\nAll gold standard PMIDs already have full text!")
        return

    print(f"\nWill attempt to harvest {len(to_harvest)} PMIDs...")

    # Harvest
    newly_downloaded, still_missing = harvest_missing(run_dir, to_harvest, args.gene)

    # Post-harvest audit
    print("\n\n" + "=" * 65)
    print("POST-HARVEST AUDIT")
    print("=" * 65)
    results2 = audit(run_dir, gold_pmids)
    print_audit(results2)

    before = len(results["has_fulltext"])
    after = len(results2["has_fulltext"])
    print(f"\n  Improvement: {before} → {after} gold PMIDs with full text")
    print(
        f"  Download rate: {before/len(gold_pmids)*100:.1f}% → {after/len(gold_pmids)*100:.1f}%"
    )


if __name__ == "__main__":
    main()
