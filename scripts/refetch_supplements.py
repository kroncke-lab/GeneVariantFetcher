#!/usr/bin/env python3
"""
Re-fetch supplements for papers with empty or missing supplement directories.

Runs the full supplement cascade (Unified[PMC+Elsevier+CrossRef] -> PMC HTML
-> DOI scraping) on papers that were originally processed before the pipeline
improvements were in place.

Usage:
    python scripts/refetch_supplements.py                      # Full run
    python scripts/refetch_supplements.py --dry-run             # Preview candidates
    python scripts/refetch_supplements.py --limit 10            # Small batch
    python scripts/refetch_supplements.py --priority-sort       # High-likelihood first
    python scripts/refetch_supplements.py --status              # Show manifest status
    python scripts/refetch_supplements.py --gene SCN5A          # Different gene
"""

import argparse
import csv
import logging
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Add project root
sys.path.insert(0, str(Path(__file__).parent.parent))

from dotenv import load_dotenv

load_dotenv()

from utils.logging_utils import setup_logging
from utils.manifest import Manifest, ManifestEntry, Stage, Status

setup_logging(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
DEFAULT_OUTPUT_DIR = "/mnt/temp2/kronckbm/gvf_output"
DEFAULT_GENE = "KCNH2"
REFETCH_HEADER = "# ===== SUPPLEMENTS (RE-FETCHED) ====="
MANIFEST_FILENAME = "refetch_manifest.json"
INTER_PAPER_DELAY = 2.0  # seconds between papers


# =============================================================================
# Data classes
# =============================================================================


@dataclass
class RefetchCandidate:
    """A paper that may need supplement re-fetching."""

    pmid: str
    full_context_path: Path
    supplements_dir: Path
    pmcid: str = ""
    doi: str = ""
    supplement_ref_count: int = 0
    priority_score: float = 0.0
    has_empty_dir: bool = False


# =============================================================================
# Phase 1: Discovery
# =============================================================================


def find_latest_pmc_dir(output_dir: str, gene: str) -> Optional[Path]:
    """Auto-detect the most recent run directory with a pmc_fulltext folder."""
    gene_dir = Path(output_dir) / gene
    if not gene_dir.exists():
        return None

    run_dirs = [d for d in gene_dir.iterdir() if d.is_dir() and d.name[0].isdigit()]
    if not run_dirs:
        return None

    # Sort by timestamp (YYYYMMDD_HHMMSS), newest first
    run_dirs.sort(key=lambda d: d.name, reverse=True)

    for run_dir in run_dirs:
        pmc_dir = run_dir / "pmc_fulltext"
        if pmc_dir.exists():
            return pmc_dir

    return None


def scan_for_candidates(pmc_dir: Path) -> List[RefetchCandidate]:
    """Find papers with empty or missing supplement directories."""
    from harvesting.supplement_reference_parser import parse_supplement_references

    candidates = []

    for fc_path in sorted(pmc_dir.glob("*_FULL_CONTEXT.md")):
        # Extract PMID from filename: {PMID}_FULL_CONTEXT.md
        pmid = fc_path.stem.replace("_FULL_CONTEXT", "")
        if not pmid.isdigit():
            continue

        supp_dir = pmc_dir / f"{pmid}_supplements"
        has_empty_dir = supp_dir.exists() and not any(supp_dir.iterdir())
        has_no_dir = not supp_dir.exists()

        # Also check if supplements already present (non-empty dir)
        if supp_dir.exists() and any(supp_dir.iterdir()):
            continue  # Already has supplements, skip

        # Check for idempotency: already re-fetched
        try:
            text = fc_path.read_text(encoding="utf-8", errors="ignore")
            if REFETCH_HEADER in text:
                continue  # Already processed
        except Exception:
            continue

        # Estimate supplement likelihood from main text
        ref_count = 0
        try:
            refs = parse_supplement_references(text)
            ref_count = (
                refs["expected_tables"]
                + refs["expected_figures"]
                + len(refs["supplement_urls"])
            )
        except Exception:
            pass

        # Priority: papers that reference supplements are higher value
        priority = ref_count * 10.0
        if has_empty_dir:
            priority += 5.0  # Slightly boost: we know pipeline tried and failed

        candidates.append(
            RefetchCandidate(
                pmid=pmid,
                full_context_path=fc_path,
                supplements_dir=supp_dir,
                supplement_ref_count=ref_count,
                priority_score=priority,
                has_empty_dir=has_empty_dir,
            )
        )

    return candidates


# =============================================================================
# Phase 2: ID Resolution
# =============================================================================


def load_pmcid_map(pmc_dir: Path) -> Dict[str, str]:
    """Load PMID->PMCID mapping from successful_downloads.csv."""
    csv_path = pmc_dir / "successful_downloads.csv"
    mapping: Dict[str, str] = {}
    if not csv_path.exists():
        return mapping

    try:
        with open(csv_path, newline="") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            for row in reader:
                if len(row) >= 2:
                    pmid_val, pmcid_val = row[0].strip(), row[1].strip()
                    if pmid_val and pmcid_val:
                        mapping[pmid_val] = pmcid_val
    except Exception as e:
        logger.warning("Failed to load successful_downloads.csv: %s", e)

    return mapping


def resolve_ids(
    candidates: List[RefetchCandidate],
    pmcid_map: Dict[str, str],
    pmc_api,
) -> None:
    """Resolve PMCID and DOI for each candidate (in-place)."""
    for cand in candidates:
        # PMCID from CSV first
        if cand.pmid in pmcid_map:
            cand.pmcid = pmcid_map[cand.pmid]
        else:
            try:
                pmcid = pmc_api.pmid_to_pmcid(cand.pmid)
                if pmcid:
                    cand.pmcid = pmcid
            except Exception as e:
                logger.debug("PMCID lookup failed for %s: %s", cand.pmid, e)

        # DOI
        try:
            doi = pmc_api.get_doi_from_pmid(cand.pmid)
            if doi:
                cand.doi = doi
        except Exception as e:
            logger.debug("DOI lookup failed for %s: %s", cand.pmid, e)


# =============================================================================
# Phase 3: Re-fetch
# =============================================================================


def refetch_one_paper(
    cand: RefetchCandidate,
    harvester,
    converter,
) -> Tuple[Status, str, int]:
    """Re-fetch supplements for a single paper.

    Returns:
        (status, message, file_count)
    """
    pmid = cand.pmid
    supp_dir = cand.supplements_dir

    # 1. Pre-fetch cleanup: remove empty supplement dir
    if supp_dir.exists() and not any(supp_dir.iterdir()):
        shutil.rmtree(supp_dir)
        logger.info("Removed empty dir: %s", supp_dir)

    # 2. Run full supplement cascade
    try:
        supp_files = harvester.get_supplemental_files(
            cand.pmcid, pmid, cand.doi
        )
    except Exception as e:
        return Status.FAILED, f"get_supplemental_files error: {e}", 0

    # 3. No supplements found — don't create dir
    if not supp_files:
        return Status.SKIPPED, "No supplements found in cascade", 0

    # 4. Download + convert via shared service
    from harvesting.supplement_processing_service import process_supplement_files

    result = process_supplement_files(
        supp_files=supp_files,
        supplements_dir=supp_dir,
        pmid=pmid,
        converter=converter,
        download_callback=lambda url, file_path, pmid_arg, filename, supp: (
            harvester.download_supplement(
                url,
                file_path,
                pmid_arg,
                filename,
                supp.get("base_url"),
                supp.get("original_url"),
            )
        ),
        logger=logger,
        sleep_seconds=0.5,
    )

    # 5. If all downloads failed, remove the dir that mkdir created
    if result.downloaded_count == 0:
        if supp_dir.exists() and not any(supp_dir.iterdir()):
            shutil.rmtree(supp_dir)
            logger.info("Removed empty dir after failed downloads: %s", supp_dir)
        return Status.FAILED, "All downloads failed", 0

    # 6. Append supplement markdown to FULL_CONTEXT.md
    if result.supplement_markdown.strip():
        try:
            with open(cand.full_context_path, "a", encoding="utf-8") as f:
                f.write(f"\n\n{REFETCH_HEADER}\n\n")
                f.write(result.supplement_markdown)
            logger.info(
                "Appended %d chars of supplement text to %s",
                len(result.supplement_markdown),
                cand.full_context_path.name,
            )
        except Exception as e:
            logger.error("Failed to append to FULL_CONTEXT: %s", e)

    return Status.SUCCESS, f"Downloaded {result.downloaded_count} files", result.downloaded_count


# =============================================================================
# Manifest helpers
# =============================================================================


def load_or_create_manifest(pmc_dir: Path, gene: str) -> Manifest:
    """Load existing manifest or create a new one."""
    manifest_path = pmc_dir / MANIFEST_FILENAME
    if manifest_path.exists():
        try:
            return Manifest.load(manifest_path)
        except Exception as e:
            logger.warning("Failed to load manifest, creating new: %s", e)

    return Manifest(stage=Stage.REFETCH, gene=gene)


def get_already_processed(manifest: Manifest) -> set:
    """Get PMIDs already processed (SUCCESS or SKIPPED)."""
    return {
        e.pmid
        for e in manifest.entries
        if e.status in (Status.SUCCESS, Status.SKIPPED)
    }


# =============================================================================
# Main
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Re-fetch supplements for papers with empty supplement dirs"
    )
    parser.add_argument("--gene", default=DEFAULT_GENE, help="Gene symbol")
    parser.add_argument(
        "--output-dir", default=DEFAULT_OUTPUT_DIR, help="Base output directory"
    )
    parser.add_argument(
        "--pmc-dir", default=None, help="Explicit pmc_fulltext directory"
    )
    parser.add_argument("--dry-run", action="store_true", help="Preview only")
    parser.add_argument("--limit", type=int, default=0, help="Max papers to process")
    parser.add_argument(
        "--priority-sort",
        action="store_true",
        help="Sort by supplement-reference likelihood",
    )
    parser.add_argument(
        "--status", action="store_true", help="Show manifest status and exit"
    )
    args = parser.parse_args()

    # Resolve pmc_dir
    if args.pmc_dir:
        pmc_dir = Path(args.pmc_dir)
    else:
        pmc_dir = find_latest_pmc_dir(args.output_dir, args.gene)

    if not pmc_dir or not pmc_dir.exists():
        print(f"ERROR: pmc_fulltext directory not found for {args.gene}")
        print(f"  Searched: {args.output_dir}/{args.gene}/*/pmc_fulltext/")
        sys.exit(1)

    print(f"Using pmc_dir: {pmc_dir}")

    # --status: show manifest and exit
    if args.status:
        manifest_path = pmc_dir / MANIFEST_FILENAME
        if not manifest_path.exists():
            print("No refetch manifest found yet.")
            sys.exit(0)
        manifest = Manifest.load(manifest_path)
        summary = manifest.summary()
        print(f"\nRefetch manifest ({len(manifest)} entries):")
        for status_name, count in sorted(summary.items()):
            print(f"  {status_name}: {count}")
        sys.exit(0)

    # Phase 1: Discovery
    print("\n--- Phase 1: Discovery ---")
    candidates = scan_for_candidates(pmc_dir)
    print(f"Found {len(candidates)} candidates (empty/missing supplement dirs)")

    if not candidates:
        print("Nothing to do!")
        sys.exit(0)

    # Filter out already-processed
    manifest = load_or_create_manifest(pmc_dir, args.gene)
    already_done = get_already_processed(manifest)
    candidates = [c for c in candidates if c.pmid not in already_done]
    print(f"After filtering already-processed: {len(candidates)} remaining")

    if args.priority_sort:
        candidates.sort(key=lambda c: c.priority_score, reverse=True)

    if args.limit > 0:
        candidates = candidates[: args.limit]
        print(f"Limited to {len(candidates)} papers")

    # Dry run: print candidates and exit
    if args.dry_run:
        print(f"\n{'PMID':<12} {'PMCID':<14} {'Refs':<6} {'Priority':<10} {'Status'}")
        print("-" * 65)
        for c in candidates:
            status = "empty_dir" if c.has_empty_dir else "no_dir"
            print(
                f"{c.pmid:<12} {c.pmcid or '?':<14} {c.supplement_ref_count:<6} "
                f"{c.priority_score:<10.1f} {status}"
            )
        print(f"\nTotal: {len(candidates)} candidates")
        total_refs = sum(c.supplement_ref_count for c in candidates)
        with_refs = sum(1 for c in candidates if c.supplement_ref_count > 0)
        print(f"  {with_refs} reference supplements in main text ({total_refs} total refs)")
        sys.exit(0)

    # Phase 2: ID Resolution
    print("\n--- Phase 2: ID Resolution ---")
    from harvesting.pmc_api import PMCAPIClient

    pmc_api = PMCAPIClient()
    pmcid_map = load_pmcid_map(pmc_dir)
    print(f"Loaded {len(pmcid_map)} PMCID mappings from CSV")
    resolve_ids(candidates, pmcid_map, pmc_api)

    with_doi = sum(1 for c in candidates if c.doi)
    with_pmcid = sum(1 for c in candidates if c.pmcid)
    print(f"Resolved: {with_pmcid} PMCIDs, {with_doi} DOIs")

    # Phase 3: Re-fetch
    print("\n--- Phase 3: Re-fetch ---")
    from harvesting.format_converters import FormatConverter
    from harvesting.orchestrator import PMCHarvester

    harvester = PMCHarvester(output_dir=str(pmc_dir), gene_symbol=args.gene)
    converter = FormatConverter()

    success_count = 0
    skip_count = 0
    fail_count = 0
    total_files = 0

    for i, cand in enumerate(candidates, 1):
        print(f"\n[{i}/{len(candidates)}] PMID {cand.pmid} (PMCID={cand.pmcid or '?'}, DOI={cand.doi or '?'})")

        status, msg, file_count = refetch_one_paper(cand, harvester, converter)

        # Record in manifest
        entry = ManifestEntry(
            pmid=cand.pmid,
            status=status,
            error_message=msg if status != Status.SUCCESS else None,
            files_created=[str(cand.supplements_dir)] if status == Status.SUCCESS else [],
        )
        manifest.add_entry(entry)

        # Save manifest after each paper (atomic, enables resume)
        manifest.save(pmc_dir / MANIFEST_FILENAME)

        if status == Status.SUCCESS:
            success_count += 1
            total_files += file_count
            print(f"  -> SUCCESS: {msg}")
        elif status == Status.SKIPPED:
            skip_count += 1
            print(f"  -> SKIPPED: {msg}")
        else:
            fail_count += 1
            print(f"  -> FAILED: {msg}")

        # Inter-paper delay
        if i < len(candidates):
            time.sleep(INTER_PAPER_DELAY)

    # Phase 4: Summary
    print("\n--- Summary ---")
    print(f"  Processed: {len(candidates)}")
    print(f"  Success:   {success_count} ({total_files} files downloaded)")
    print(f"  Skipped:   {skip_count}")
    print(f"  Failed:    {fail_count}")
    print(f"  Manifest:  {pmc_dir / MANIFEST_FILENAME}")


if __name__ == "__main__":
    main()
