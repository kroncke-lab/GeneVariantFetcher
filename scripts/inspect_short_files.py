#!/usr/bin/env python3
"""
Batch inspection script for short/failed DATA_ZONES.md and FULL_CONTEXT.md files.

This script analyzes files that are suspiciously small (< 500 chars) and identifies
common failure patterns like CAPTCHA, 403 errors, paywalls, etc.

Usage:
    python scripts/inspect_short_files.py /path/to/pmc_fulltext
    python scripts/inspect_short_files.py /path/to/pmc_fulltext --threshold 1000
    python scripts/inspect_short_files.py /path/to/pmc_fulltext --output report.json
"""

import argparse
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

# Common error patterns that indicate scraper/download failures
ERROR_PATTERNS = {
    "captcha": [
        r"captcha",
        r"verify you are human",
        r"robot check",
        r"are you a robot",
        r"security check",
        r"challenge-platform",
        r"cf-turnstile",
        r"hcaptcha",
        r"recaptcha",
        r"i('|')m not a robot",
        r"prove you('|')re human",
        r"human verification",
    ],
    "cloudflare": [
        r"cloudflare",
        r"cf-ray",
        r"cf[-_]?browser",
        r"checking your browser",
        r"please wait\.\.\.",
        r"enable javascript and cookies",
        r"ray id:",
        r"performance & security by cloudflare",
        r"ddos protection",
        r"just a moment",
        r"attention required",
    ],
    "javascript_required": [
        r"enable javascript",
        r"javascript is required",
        r"javascript is disabled",
        r"requires javascript",
        r"please enable javascript",
        r"javascript must be enabled",
        r"this site requires javascript",
        r"you need to enable javascript",
        r"browser does not support javascript",
        r"noscript",
    ],
    "403_forbidden": [
        r"403 forbidden",
        r"access denied",
        r"permission denied",
        r"not authorized",
        r"unauthorized access",
        r"you don('|')t have permission",
        r"forbidden",
        r"access to .* denied",
        r"request blocked",
    ],
    "404_not_found": [
        r"404 not found",
        r"page not found",
        r"resource not found",
        r"the page you requested",
        r"could not be found",
        r"no longer available",
        r"has been removed",
        r"does not exist",
    ],
    "paywall": [
        r"paywall",
        r"subscribe to",
        r"subscription required",
        r"purchase this article",
        r"buy this article",
        r"institutional access",
        r"sign in to access",
        r"login required",
        r"full text not available",
        r"purchase access",
        r"rent this article",
        r"get access",
        r"buy instant access",
        r"members only",
        r"premium content",
        r"pay per view",
    ],
    "rate_limit": [
        r"rate limit",
        r"too many requests",
        r"slow down",
        r"try again later",
        r"request limit exceeded",
        r"throttl",
        r"quota exceeded",
        r"exceeded.*limit",
    ],
    "timeout": [
        r"timeout",
        r"timed out",
        r"connection timed out",
        r"request timed out",
        r"server.*not respond",
        r"gateway timeout",
        r"504 gateway",
    ],
    "empty_zones": [
        r"no high-value data zones identified",
        r"no relevant data zones",
        r"no data zones found",
    ],
    "conversion_error": [
        r"error converting",
        r"failed to convert",
        r"conversion failed",
        r"could not parse",
        r"invalid pdf",
        r"corrupted file",
        r"unable to extract",
        r"extraction failed",
        r"malformed",
        r"encoding error",
    ],
    "network_error": [
        r"network error",
        r"connection refused",
        r"connection reset",
        r"dns resolution failed",
        r"host not found",
        r"unreachable",
        r"socket error",
        r"ssl error",
        r"certificate error",
    ],
    "bot_detection": [
        r"automated access",
        r"bot detection",
        r"suspicious activity",
        r"unusual traffic",
        r"automated queries",
        r"scraping",
        r"crawler",
        r"selenium",
        r"webdriver",
    ],
    "geo_restriction": [
        r"not available in your",
        r"region restricted",
        r"country not supported",
        r"geographic restriction",
        r"unavailable in your location",
    ],
    "session_expired": [
        r"session expired",
        r"session timeout",
        r"please log in again",
        r"authentication required",
        r"session invalid",
    ],
}


def extract_pmid_from_filename(filename: str) -> str:
    """Extract PMID from filename like 'PMID_12345678_FULL_CONTEXT.md'"""
    match = re.search(r"(\d{7,8})", filename)
    return match.group(1) if match else "UNKNOWN"


def detect_error_patterns(content: str) -> Dict[str, List[str]]:
    """Detect error patterns in file content."""
    content_lower = content.lower()
    detected = {}

    for category, patterns in ERROR_PATTERNS.items():
        matches = []
        for pattern in patterns:
            if re.search(pattern, content_lower):
                # Find the actual match in original content for context
                match = re.search(pattern, content_lower)
                if match:
                    # Get surrounding context (50 chars before and after)
                    start = max(0, match.start() - 50)
                    end = min(len(content), match.end() + 50)
                    context = content[start:end].strip()
                    matches.append(context)
        if matches:
            detected[category] = matches[:3]  # Limit to 3 examples per category

    return detected


def analyze_content_structure(content: str) -> Dict[str, Any]:
    """Analyze the structure of the content."""
    lines = content.split("\n")
    non_empty_lines = [line for line in lines if line.strip()]

    # Check for markdown structure
    has_headers = bool(re.search(r"^#+\s", content, re.MULTILINE))
    has_tables = "|" in content and "-|-" in content or "| --- |" in content
    has_lists = bool(re.search(r"^[\*\-]\s", content, re.MULTILINE))

    # Word count
    words = len(content.split())

    # Check for actual paper content indicators
    has_abstract = "abstract" in content.lower()
    has_introduction = "introduction" in content.lower()
    has_methods = any(
        m in content.lower() for m in ["method", "materials", "procedure"]
    )
    has_results = "results" in content.lower()
    has_references = any(
        r in content.lower() for r in ["references", "bibliography", "cited"]
    )

    content_score = sum(
        [
            has_abstract,
            has_introduction,
            has_methods,
            has_results,
            has_references,
            has_tables,
        ]
    )

    return {
        "total_chars": len(content),
        "total_lines": len(lines),
        "non_empty_lines": len(non_empty_lines),
        "word_count": words,
        "has_headers": has_headers,
        "has_tables": has_tables,
        "has_lists": has_lists,
        "content_indicators": {
            "abstract": has_abstract,
            "introduction": has_introduction,
            "methods": has_methods,
            "results": has_results,
            "references": has_references,
        },
        "content_score": content_score,  # 0-6, higher is better
    }


def get_companion_file(filepath: Path) -> Optional[Path]:
    """Get the companion file (DATA_ZONES <-> FULL_CONTEXT)."""
    name = filepath.name
    parent = filepath.parent

    if "_DATA_ZONES.md" in name:
        companion_name = name.replace("_DATA_ZONES.md", "_FULL_CONTEXT.md")
    elif "_FULL_CONTEXT.md" in name:
        companion_name = name.replace("_FULL_CONTEXT.md", "_DATA_ZONES.md")
    else:
        return None

    companion = parent / companion_name
    return companion if companion.exists() else None


def analyze_companion_file(filepath: Path) -> Optional[Dict[str, Any]]:
    """Analyze the companion file to understand the root cause."""
    companion = get_companion_file(filepath)
    if not companion:
        return None

    try:
        content = companion.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return None

    errors = detect_error_patterns(content)
    structure = analyze_content_structure(content)

    return {
        "companion_file": companion.name,
        "companion_chars": len(content),
        "companion_errors": errors,
        "companion_content_score": structure["content_score"],
        "companion_has_tables": structure["has_tables"],
    }


def classify_file(
    content: str, filename: str, threshold: int, filepath: Optional[Path] = None
) -> Dict[str, Any]:
    """Classify a file based on its content and companion file analysis."""
    pmid = extract_pmid_from_filename(filename)
    char_count = len(content)

    # Detect error patterns
    errors = detect_error_patterns(content)

    # Analyze structure
    structure = analyze_content_structure(content)

    # Analyze companion file if this is a DATA_ZONES file
    companion_analysis = None
    root_cause = None

    if filepath and "_DATA_ZONES.md" in filename:
        companion_analysis = analyze_companion_file(filepath)

        if companion_analysis:
            comp_chars = companion_analysis["companion_chars"]
            comp_errors = companion_analysis["companion_errors"]
            comp_score = companion_analysis["companion_content_score"]

            # Determine root cause
            if comp_chars < threshold and comp_errors:
                # FULL_CONTEXT also failed - download/scraping issue
                primary_error = list(comp_errors.keys())[0]
                root_cause = f"download_failed_{primary_error}"
            elif comp_chars < threshold:
                # FULL_CONTEXT is small but no errors - possibly empty PDF
                root_cause = "source_empty_or_corrupt"
            elif comp_score >= 3 and not comp_errors:
                # FULL_CONTEXT looks valid - scout genuinely found nothing
                root_cause = "legitimate_no_zones"
            elif comp_errors:
                # FULL_CONTEXT has content but also errors
                primary_error = list(comp_errors.keys())[0]
                root_cause = f"partial_content_{primary_error}"
            else:
                # FULL_CONTEXT has content, low score, no errors
                root_cause = "low_relevance_content"

    # Determine classification
    if char_count < threshold:
        if root_cause:
            classification = root_cause
        elif errors:
            primary_error = list(errors.keys())[0]
            classification = f"error_{primary_error}"
        elif structure["content_score"] == 0:
            classification = "empty_or_minimal"
        else:
            classification = "truncated_content"
    else:
        if errors:
            classification = "has_errors_but_large"
        elif structure["content_score"] >= 3:
            classification = "likely_valid"
        else:
            classification = "questionable_content"

    result = {
        "pmid": pmid,
        "filename": filename,
        "char_count": char_count,
        "classification": classification,
        "error_patterns": errors,
        "structure": structure,
        "first_500_chars": content[:500] if char_count > 0 else "",
    }

    if companion_analysis:
        result["companion_analysis"] = companion_analysis
        result["root_cause"] = root_cause

    return result


def scan_directory(
    directory: Path, threshold: int = 500, file_patterns: List[str] = None
) -> Dict[str, Any]:
    """Scan a directory for short/problematic files."""

    if file_patterns is None:
        file_patterns = ["*_DATA_ZONES.md", "*_FULL_CONTEXT.md"]

    results = {
        "scan_timestamp": datetime.now().isoformat(),
        "directory": str(directory),
        "threshold": threshold,
        "summary": {
            "total_files": 0,
            "short_files": 0,
            "by_classification": Counter(),
            "by_error_type": Counter(),
        },
        "short_files": [],
        "all_files": [],
    }

    for pattern in file_patterns:
        for filepath in directory.glob(pattern):
            try:
                content = filepath.read_text(encoding="utf-8", errors="replace")
            except Exception as e:
                content = f"[ERROR READING FILE: {e}]"

            analysis = classify_file(content, filepath.name, threshold, filepath)
            analysis["filepath"] = str(filepath)
            analysis["file_type"] = (
                "DATA_ZONES" if "DATA_ZONES" in filepath.name else "FULL_CONTEXT"
            )

            results["summary"]["total_files"] += 1
            results["all_files"].append(
                {
                    "pmid": analysis["pmid"],
                    "filename": analysis["filename"],
                    "char_count": analysis["char_count"],
                    "classification": analysis["classification"],
                }
            )

            if analysis["char_count"] < threshold:
                results["summary"]["short_files"] += 1
                results["short_files"].append(analysis)
                results["summary"]["by_classification"][analysis["classification"]] += 1

                for error_type in analysis["error_patterns"].keys():
                    results["summary"]["by_error_type"][error_type] += 1

    # Convert Counters to dicts for JSON serialization
    results["summary"]["by_classification"] = dict(
        results["summary"]["by_classification"]
    )
    results["summary"]["by_error_type"] = dict(results["summary"]["by_error_type"])

    return results


def print_report(results: Dict[str, Any], verbose: bool = False):
    """Print a human-readable report."""
    summary = results["summary"]

    print("\n" + "=" * 70)
    print("FILE INSPECTION REPORT")
    print("=" * 70)
    print(f"Directory: {results['directory']}")
    print(f"Threshold: {results['threshold']} characters")
    print(f"Scan time: {results['scan_timestamp']}")
    print()

    print("SUMMARY")
    print("-" * 40)
    print(f"Total files scanned:  {summary['total_files']}")
    print(
        f"Short/problem files:  {summary['short_files']} ({100 * summary['short_files'] / max(1, summary['total_files']):.1f}%)"
    )
    print()

    if summary["by_classification"]:
        print("BY CLASSIFICATION:")
        # Group by root cause category for better readability
        download_issues = []
        content_issues = []
        other_issues = []

        for classification, count in sorted(
            summary["by_classification"].items(), key=lambda x: -x[1]
        ):
            if any(
                x in classification
                for x in [
                    "download_failed",
                    "error_captcha",
                    "error_cloudflare",
                    "error_403",
                    "error_javascript",
                    "error_rate_limit",
                    "error_timeout",
                    "error_network",
                    "error_bot",
                ]
            ):
                download_issues.append((classification, count))
            elif any(
                x in classification
                for x in ["legitimate_no_zones", "low_relevance", "source_empty"]
            ):
                content_issues.append((classification, count))
            else:
                other_issues.append((classification, count))

        if download_issues:
            print("\n  [DOWNLOAD/ACCESS FAILURES - Retry candidates]")
            for classification, count in download_issues:
                print(f"    {classification}: {count}")

        if content_issues:
            print("\n  [CONTENT ISSUES - May be legitimate]")
            for classification, count in content_issues:
                print(f"    {classification}: {count}")

        if other_issues:
            print("\n  [OTHER]")
            for classification, count in other_issues:
                print(f"    {classification}: {count}")
        print()

    if summary["by_error_type"]:
        print("BY ERROR TYPE (in file content):")
        for error_type, count in sorted(
            summary["by_error_type"].items(), key=lambda x: -x[1]
        ):
            print(f"  {error_type}: {count}")
        print()

    # Calculate actionable summary
    actionable = sum(
        1
        for f in results["short_files"]
        if "download_failed" in f.get("classification", "")
        or f.get("classification", "").startswith("error_")
    )
    legitimate = sum(
        1
        for f in results["short_files"]
        if f.get("classification") == "legitimate_no_zones"
    )

    if actionable or legitimate:
        print("ACTIONABLE SUMMARY")
        print("-" * 40)
        print(f"  Files to RETRY (download failures): {actionable}")
        print(f"  Files that are LEGITIMATE (no zones): {legitimate}")
        print(
            f"  Files needing INVESTIGATION: {summary['short_files'] - actionable - legitimate}"
        )
        print()

    if verbose and results["short_files"]:
        print("\nDETAILED FILE ANALYSIS")
        print("-" * 70)

        # Group by classification
        by_class = defaultdict(list)
        for f in results["short_files"]:
            by_class[f["classification"]].append(f)

        for classification, files in sorted(by_class.items()):
            print(f"\n### {classification.upper()} ({len(files)} files) ###\n")

            for f in files[:5]:  # Show first 5 per category
                print(f"PMID: {f['pmid']}")
                print(f"File: {f['filename']}")
                print(f"Size: {f['char_count']} chars")

                if f.get("root_cause"):
                    print(f"Root Cause: {f['root_cause']}")

                if f.get("companion_analysis"):
                    ca = f["companion_analysis"]
                    print(
                        f"FULL_CONTEXT: {ca['companion_chars']} chars, score={ca['companion_content_score']}"
                    )
                    if ca["companion_errors"]:
                        print(
                            f"FULL_CONTEXT errors: {', '.join(ca['companion_errors'].keys())}"
                        )

                if f["error_patterns"]:
                    print(f"This file errors: {', '.join(f['error_patterns'].keys())}")

                print(f"Preview: {f['first_500_chars'][:200]}...")
                print()

            if len(files) > 5:
                remaining_pmids = [f["pmid"] for f in files[5:]]
                print(
                    f"... and {len(files) - 5} more: {', '.join(remaining_pmids[:10])}"
                )
                if len(remaining_pmids) > 10:
                    print(f"    (plus {len(remaining_pmids) - 10} more)")
                print()


def generate_remediation_csv(results: Dict[str, Any], output_path: Path):
    """Generate a CSV for remediation tracking."""
    import csv

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "PMID",
                "Filename",
                "CharCount",
                "Classification",
                "RootCause",
                "ErrorTypes",
                "ContentScore",
                "CompanionChars",
                "CompanionErrors",
                "Action",
            ]
        )

        for file_info in results["short_files"]:
            classification = file_info["classification"]
            root_cause = file_info.get("root_cause", "")

            # Determine action based on classification
            if "download_failed" in classification or classification.startswith(
                "error_"
            ):
                action = "RETRY"
            elif classification == "legitimate_no_zones":
                action = "SKIP"  # Legitimate - no action needed
            elif classification == "source_empty_or_corrupt":
                action = "CHECK_SOURCE"
            else:
                action = "INVESTIGATE"

            # Get companion info
            companion = file_info.get("companion_analysis", {})
            companion_chars = companion.get("companion_chars", "")
            companion_errors = "|".join(companion.get("companion_errors", {}).keys())

            writer.writerow(
                [
                    file_info["pmid"],
                    file_info["filename"],
                    file_info["char_count"],
                    classification,
                    root_cause,
                    "|".join(file_info["error_patterns"].keys()),
                    file_info["structure"]["content_score"],
                    companion_chars,
                    companion_errors,
                    action,
                ]
            )

    print(f"Remediation CSV written to: {output_path}")


def export_retry_pmids(results: Dict[str, Any], output_path: Path):
    """Export PMIDs that need retry to a text file (one per line)."""
    retry_pmids = []

    for file_info in results["short_files"]:
        classification = file_info["classification"]
        # Include files that failed due to download/access issues
        if any(
            x in classification
            for x in [
                "download_failed",
                "error_captcha",
                "error_cloudflare",
                "error_403",
                "error_javascript",
                "error_rate_limit",
                "error_timeout",
                "error_network",
                "error_bot",
                "error_paywall",
            ]
        ):
            pmid = file_info["pmid"]
            if pmid and pmid != "UNKNOWN" and pmid not in retry_pmids:
                retry_pmids.append(pmid)

    with open(output_path, "w", encoding="utf-8") as f:
        for pmid in sorted(retry_pmids):
            f.write(f"{pmid}\n")

    print(f"Exported {len(retry_pmids)} PMIDs for retry to: {output_path}")
    return retry_pmids


def main():
    parser = argparse.ArgumentParser(
        description="Inspect short/failed markdown files for error patterns"
    )
    parser.add_argument(
        "directory",
        type=Path,
        help="Directory containing *_DATA_ZONES.md and *_FULL_CONTEXT.md files",
    )
    parser.add_argument(
        "--threshold",
        "-t",
        type=int,
        default=500,
        help="Character threshold for 'short' files (default: 500)",
    )
    parser.add_argument(
        "--output", "-o", type=Path, help="Output JSON file for detailed results"
    )
    parser.add_argument(
        "--csv", type=Path, help="Output CSV file for remediation tracking"
    )
    parser.add_argument(
        "--retry-list",
        type=Path,
        help="Output text file with PMIDs to retry (one per line)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Show detailed file-by-file analysis",
    )
    parser.add_argument(
        "--data-zones-only", action="store_true", help="Only scan DATA_ZONES.md files"
    )
    parser.add_argument(
        "--full-context-only",
        action="store_true",
        help="Only scan FULL_CONTEXT.md files",
    )

    args = parser.parse_args()

    if not args.directory.exists():
        print(f"Error: Directory not found: {args.directory}")
        return 1

    # Determine file patterns
    patterns = None
    if args.data_zones_only:
        patterns = ["*_DATA_ZONES.md"]
    elif args.full_context_only:
        patterns = ["*_FULL_CONTEXT.md"]

    print(f"Scanning directory: {args.directory}")
    print(f"Threshold: {args.threshold} characters")

    results = scan_directory(args.directory, args.threshold, patterns)

    print_report(results, verbose=args.verbose)

    if args.output:
        # Remove first_500_chars from detailed output to keep file size reasonable
        for f in results["short_files"]:
            f.pop("first_500_chars", None)

        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
        print(f"\nDetailed results written to: {args.output}")

    if args.csv:
        generate_remediation_csv(results, args.csv)

    if args.retry_list:
        export_retry_pmids(results, args.retry_list)

    return 0


if __name__ == "__main__":
    exit(main())
