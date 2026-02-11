#!/usr/bin/env python3
"""
Repair script for malformed extraction JSON files.

Fixes common issues:
1. Missing PMID in paper_metadata (extracts from filename)
2. Missing gene_symbol in variants (extracts from filename or paper_metadata)
3. Missing required fields that cause NOT NULL constraint failures

Usage:
    python scripts/repair_extraction_jsons.py /path/to/extractions
    python scripts/repair_extraction_jsons.py /path/to/extractions --dry-run
    python scripts/repair_extraction_jsons.py /path/to/extractions --backup
"""

import argparse
import json
import shutil
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.pmid_utils import extract_gene_from_filename, extract_pmid_from_filename


def validate_extraction_json(
    data: Dict[str, Any], filename: str
) -> Dict[str, List[str]]:
    """
    Validate an extraction JSON and return issues found.

    Returns:
        Dict with 'errors' (blocking issues) and 'warnings' (non-blocking)
    """
    issues = {
        "errors": [],
        "warnings": [],
    }

    # Check paper_metadata
    paper_meta = data.get("paper_metadata", {})
    if not paper_meta:
        issues["errors"].append("Missing paper_metadata section")
    else:
        if not paper_meta.get("pmid"):
            issues["errors"].append("Missing paper_metadata.pmid")
        elif paper_meta.get("pmid") == "UNKNOWN":
            issues["warnings"].append("paper_metadata.pmid is 'UNKNOWN'")

        if not paper_meta.get("title"):
            issues["warnings"].append("Missing paper_metadata.title")

    # Check variants
    variants = data.get("variants", [])
    if not variants:
        issues["warnings"].append("No variants found in extraction")
    else:
        for i, variant in enumerate(variants):
            if not variant.get("gene_symbol"):
                issues["errors"].append(f"Variant {i}: missing gene_symbol")

            # Check for at least one notation
            has_notation = any(
                [
                    variant.get("cdna_notation"),
                    variant.get("protein_notation"),
                    variant.get("genomic_position"),
                ]
            )
            if not has_notation:
                issues["warnings"].append(
                    f"Variant {i}: no notation (cdna/protein/genomic)"
                )

            # Check individuals
            individuals = variant.get("individuals", [])
            for j, ind in enumerate(individuals):
                if not ind.get("affected_status"):
                    issues["warnings"].append(
                        f"Variant {i}, Individual {j}: missing affected_status"
                    )

    # Check extraction_metadata
    if not data.get("extraction_metadata"):
        issues["warnings"].append("Missing extraction_metadata section")

    return issues


def repair_extraction_json(
    data: Dict[str, Any], filename: str, gene_symbol_override: Optional[str] = None
) -> Tuple[Dict[str, Any], List[str]]:
    """
    Repair common issues in an extraction JSON.

    Args:
        data: The JSON data to repair
        filename: Original filename (for PMID/gene extraction)
        gene_symbol_override: Force this gene symbol if provided

    Returns:
        Tuple of (repaired_data, list of repairs made)
    """
    repairs = []
    data = json.loads(json.dumps(data))  # Deep copy

    # Extract info from filename
    filename_pmid = extract_pmid_from_filename(filename)
    filename_gene = extract_gene_from_filename(filename)

    # Use override or filename-derived gene
    default_gene = gene_symbol_override or filename_gene

    # Repair paper_metadata
    if "paper_metadata" not in data:
        data["paper_metadata"] = {}
        repairs.append("Created missing paper_metadata section")

    paper_meta = data["paper_metadata"]

    # Fix missing PMID
    if not paper_meta.get("pmid") or paper_meta.get("pmid") == "UNKNOWN":
        if filename_pmid:
            old_value = paper_meta.get("pmid", "MISSING")
            paper_meta["pmid"] = filename_pmid
            repairs.append(
                f"Set paper_metadata.pmid from filename: {old_value} -> {filename_pmid}"
            )
        else:
            repairs.append("WARNING: Could not determine PMID from filename")

    # Fix missing title
    if not paper_meta.get("title"):
        paper_meta["title"] = f"Paper {paper_meta.get('pmid', 'Unknown')}"
        repairs.append("Set default paper_metadata.title")

    # Ensure gene_symbol in paper_metadata
    if not paper_meta.get("gene_symbol") and default_gene:
        paper_meta["gene_symbol"] = default_gene
        repairs.append(f"Set paper_metadata.gene_symbol: {default_gene}")

    # Repair variants
    variants = data.get("variants", [])
    for i, variant in enumerate(variants):
        # Fix missing gene_symbol in variant
        if not variant.get("gene_symbol"):
            # Try to get from paper_metadata, then filename
            gene = paper_meta.get("gene_symbol") or default_gene or "UNKNOWN_GENE"
            variant["gene_symbol"] = gene
            repairs.append(f"Variant {i}: set gene_symbol to '{gene}'")

        # Ensure individuals have affected_status
        individuals = variant.get("individuals", [])
        for j, ind in enumerate(individuals):
            if not ind.get("affected_status"):
                ind["affected_status"] = "Ambiguous"
                repairs.append(
                    f"Variant {i}, Individual {j}: set affected_status to 'Ambiguous'"
                )

    # Ensure extraction_metadata exists
    if "extraction_metadata" not in data:
        data["extraction_metadata"] = {
            "total_variants_found": len(variants),
            "extraction_confidence": "repaired",
            "repair_timestamp": datetime.now().isoformat(),
        }
        repairs.append("Created missing extraction_metadata section")
    else:
        data["extraction_metadata"]["repair_timestamp"] = datetime.now().isoformat()
        if repairs:
            data["extraction_metadata"]["was_repaired"] = True
            data["extraction_metadata"]["repairs_applied"] = len(repairs)

    return data, repairs


def process_directory(
    directory: Path,
    dry_run: bool = False,
    backup: bool = False,
    gene_symbol: Optional[str] = None,
    validate_only: bool = False,
) -> Dict[str, Any]:
    """
    Process all extraction JSON files in a directory.

    Args:
        directory: Path to extractions directory
        dry_run: If True, don't modify files
        backup: If True, create .bak files before modifying
        gene_symbol: Override gene symbol for all files
        validate_only: Only validate, don't repair

    Returns:
        Summary statistics
    """
    results = {
        "timestamp": datetime.now().isoformat(),
        "directory": str(directory),
        "dry_run": dry_run,
        "total_files": 0,
        "files_with_errors": 0,
        "files_repaired": 0,
        "files_skipped": 0,
        "total_repairs": 0,
        "by_issue": Counter(),
        "details": [],
    }

    json_files = list(directory.glob("*_PMID_*.json"))
    results["total_files"] = len(json_files)

    print(f"Found {len(json_files)} extraction JSON files")

    for filepath in sorted(json_files):
        file_result = {
            "filename": filepath.name,
            "status": "unknown",
            "issues": {},
            "repairs": [],
        }

        try:
            with open(filepath, "r", encoding="utf-8") as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            file_result["status"] = "invalid_json"
            file_result["error"] = str(e)
            results["files_skipped"] += 1
            results["details"].append(file_result)
            print(f"  SKIP {filepath.name}: Invalid JSON - {e}")
            continue

        # Validate
        issues = validate_extraction_json(data, filepath.name)
        file_result["issues"] = issues

        # Count issues
        for error in issues["errors"]:
            results["by_issue"][error.split(":")[0]] += 1
        for warning in issues["warnings"]:
            results["by_issue"][warning.split(":")[0]] += 1

        if issues["errors"]:
            results["files_with_errors"] += 1

        if validate_only:
            if issues["errors"] or issues["warnings"]:
                file_result["status"] = "has_issues"
                print(
                    f"  ISSUES {filepath.name}: {len(issues['errors'])} errors, {len(issues['warnings'])} warnings"
                )
            else:
                file_result["status"] = "valid"
            results["details"].append(file_result)
            continue

        # Repair if there are errors
        if issues["errors"]:
            repaired_data, repairs = repair_extraction_json(
                data, filepath.name, gene_symbol
            )
            file_result["repairs"] = repairs
            results["total_repairs"] += len(repairs)

            if repairs:
                if dry_run:
                    file_result["status"] = "would_repair"
                    print(f"  WOULD REPAIR {filepath.name}: {len(repairs)} repairs")
                    for r in repairs:
                        print(f"    - {r}")
                else:
                    # Create backup if requested
                    if backup:
                        backup_path = filepath.with_suffix(".json.bak")
                        shutil.copy2(filepath, backup_path)

                    # Write repaired file
                    with open(filepath, "w", encoding="utf-8") as f:
                        json.dump(repaired_data, f, indent=2)

                    file_result["status"] = "repaired"
                    results["files_repaired"] += 1
                    print(f"  REPAIRED {filepath.name}: {len(repairs)} repairs")
                    for r in repairs:
                        print(f"    - {r}")
        else:
            file_result["status"] = "ok"

        results["details"].append(file_result)

    # Convert Counter to dict for JSON serialization
    results["by_issue"] = dict(results["by_issue"])

    return results


def print_summary(results: Dict[str, Any]):
    """Print a summary of processing results."""
    print("\n" + "=" * 60)
    print("REPAIR SUMMARY")
    print("=" * 60)
    print(f"Directory: {results['directory']}")
    print(f"Dry run: {results['dry_run']}")
    print()
    print(f"Total files:        {results['total_files']}")
    print(f"Files with errors:  {results['files_with_errors']}")
    print(f"Files repaired:     {results['files_repaired']}")
    print(f"Files skipped:      {results['files_skipped']}")
    print(f"Total repairs made: {results['total_repairs']}")
    print()

    if results["by_issue"]:
        print("ISSUES BY TYPE:")
        for issue, count in sorted(results["by_issue"].items(), key=lambda x: -x[1]):
            print(f"  {issue}: {count}")


def main():
    parser = argparse.ArgumentParser(
        description="Repair malformed extraction JSON files"
    )
    parser.add_argument(
        "directory", type=Path, help="Directory containing extraction JSON files"
    )
    parser.add_argument(
        "--dry-run",
        "-n",
        action="store_true",
        help="Show what would be repaired without modifying files",
    )
    parser.add_argument(
        "--backup", "-b", action="store_true", help="Create .bak files before modifying"
    )
    parser.add_argument(
        "--validate-only",
        "-v",
        action="store_true",
        help="Only validate files, don't repair",
    )
    parser.add_argument(
        "--gene-symbol", "-g", type=str, help="Override gene symbol for all files"
    )
    parser.add_argument(
        "--output", "-o", type=Path, help="Output JSON file for detailed results"
    )

    args = parser.parse_args()

    if not args.directory.exists():
        print(f"Error: Directory not found: {args.directory}")
        return 1

    print(f"Processing: {args.directory}")
    if args.dry_run:
        print("DRY RUN - no files will be modified")
    if args.validate_only:
        print("VALIDATE ONLY - checking for issues")

    results = process_directory(
        args.directory,
        dry_run=args.dry_run,
        backup=args.backup,
        gene_symbol=args.gene_symbol,
        validate_only=args.validate_only,
    )

    print_summary(results)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
        print(f"\nDetailed results written to: {args.output}")

    return 0


if __name__ == "__main__":
    exit(main())
