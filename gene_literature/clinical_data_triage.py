#!/usr/bin/env python3
"""
Clinical Data Triage Tool

A standalone tool to triage scientific papers for original clinical data.
Determines if a paper contains NEW patient-level clinical data vs. reviews,
animal studies, or other non-clinical research.

Usage:
    # Single paper triage
    python clinical_data_triage.py --title "..." --abstract "..." --gene SCN5A

    # Batch triage from CSV
    python clinical_data_triage.py --input papers.csv --gene BRCA1 --output results.json

    # Interactive mode
    python clinical_data_triage.py --interactive
"""

import argparse
import json
import csv
import sys
import logging
from pathlib import Path
from typing import List, Dict, Optional

from pipeline.filters import ClinicalDataTriageFilter
from models import Paper

# Configure logging using centralized utility
from utils.logging_utils import setup_logging, get_logger
setup_logging(level=logging.INFO)
logger = get_logger(__name__)


def triage_single_paper(
    title: str,
    abstract: str,
    gene: str,
    pmid: Optional[str] = None,
    model: str = "gpt-4o-mini"
) -> dict:
    """
    Triage a single paper.

    Args:
        title: Paper title.
        abstract: Paper abstract or introduction.
        gene: Gene symbol.
        pmid: Optional PMID.
        model: LLM model to use.

    Returns:
        Triage result dictionary.
    """
    triage_filter = ClinicalDataTriageFilter(model=model)
    result = triage_filter.triage(
        title=title,
        abstract=abstract,
        gene=gene,
        pmid=pmid
    )
    return result


def triage_from_csv(
    input_file: Path,
    gene: str,
    output_file: Optional[Path] = None,
    model: str = "gpt-4o-mini"
) -> List[dict]:
    """
    Triage multiple papers from a CSV file.

    Expected CSV columns: pmid, title, abstract
    Optional columns: gene (overrides command-line gene)

    Args:
        input_file: Path to input CSV file.
        gene: Default gene symbol (used if not in CSV).
        output_file: Optional path to save results as JSON.
        model: LLM model to use.

    Returns:
        List of triage results.
    """
    logger.info(f"Reading papers from {input_file}")

    triage_filter = ClinicalDataTriageFilter(model=model)
    results = []

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)

        for i, row in enumerate(reader, 1):
            pmid = row.get('pmid', f'row_{i}')
            title = row.get('title', '')
            abstract = row.get('abstract', '')
            paper_gene = row.get('gene', gene)

            if not title or not abstract:
                logger.warning(f"Skipping row {i} (PMID: {pmid}): missing title or abstract")
                continue

            logger.info(f"Triaging paper {i} (PMID: {pmid})")

            result = triage_filter.triage(
                title=title,
                abstract=abstract,
                gene=paper_gene,
                pmid=pmid
            )

            results.append(result)

    logger.info(f"Triaged {len(results)} papers")

    # Count decisions
    keep_count = sum(1 for r in results if r['decision'] == 'KEEP')
    drop_count = sum(1 for r in results if r['decision'] == 'DROP')

    logger.info(f"Results: {keep_count} KEEP, {drop_count} DROP")

    # Save to file if specified
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Results saved to {output_file}")

    return results


def interactive_mode(model: str = "gpt-4o-mini"):
    """
    Interactive triage mode - prompt user for input.

    Args:
        model: LLM model to use.
    """
    print("\n" + "="*80)
    print("Clinical Data Triage Tool - Interactive Mode")
    print("="*80)
    print("\nEnter paper details below (or 'quit' to exit):\n")

    triage_filter = ClinicalDataTriageFilter(model=model)

    while True:
        try:
            # Get input
            print("-" * 80)
            gene = input("Gene symbol (e.g., SCN5A): ").strip()

            if gene.lower() == 'quit':
                break

            pmid = input("PMID (optional, press Enter to skip): ").strip() or None
            title = input("Paper title: ").strip()

            if not title:
                print("Title is required. Please try again.\n")
                continue

            print("\nPaste abstract (press Enter twice when done):")
            abstract_lines = []
            while True:
                line = input()
                if not line:
                    break
                abstract_lines.append(line)

            abstract = " ".join(abstract_lines).strip()

            if not abstract:
                print("Abstract is required. Please try again.\n")
                continue

            # Triage
            print("\n⏳ Triaging paper...")
            result = triage_filter.triage(
                title=title,
                abstract=abstract,
                gene=gene,
                pmid=pmid
            )

            # Display result
            print("\n" + "="*80)
            print("TRIAGE RESULT")
            print("="*80)
            print(json.dumps(result, indent=2))
            print("="*80 + "\n")

            # Ask to continue
            continue_choice = input("Triage another paper? (y/n): ").strip().lower()
            if continue_choice != 'y':
                break

        except KeyboardInterrupt:
            print("\n\nExiting...")
            break
        except Exception as e:
            logger.error(f"Error: {e}")
            continue

    print("\nGoodbye!")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Clinical Data Triage Tool - Determine if papers contain original clinical data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single paper triage
  python clinical_data_triage.py \\
    --title "A novel SCN5A mutation in Brugada syndrome" \\
    --abstract "We report a family with..." \\
    --gene SCN5A

  # Batch triage from CSV
  python clinical_data_triage.py \\
    --input papers.csv \\
    --gene BRCA1 \\
    --output triage_results.json

  # Interactive mode
  python clinical_data_triage.py --interactive
        """
    )

    # Mode selection
    mode_group = parser.add_mutually_exclusive_group(required=False)
    mode_group.add_argument(
        '--interactive', '-i',
        action='store_true',
        help='Run in interactive mode'
    )
    mode_group.add_argument(
        '--input',
        type=Path,
        help='Input CSV file with papers (columns: pmid, title, abstract)'
    )

    # Single paper arguments
    parser.add_argument('--title', type=str, help='Paper title')
    parser.add_argument('--abstract', type=str, help='Paper abstract')
    parser.add_argument('--pmid', type=str, help='PubMed ID (optional)')

    # Common arguments
    parser.add_argument(
        '--gene', '-g',
        type=str,
        default='the gene of interest',
        help='Gene symbol (e.g., SCN5A, BRCA1)'
    )
    parser.add_argument(
        '--output', '-o',
        type=Path,
        help='Output file for results (JSON format)'
    )
    parser.add_argument(
        '--model', '-m',
        type=str,
        default='gpt-4o-mini',
        help='LLM model to use (default: gpt-4o-mini)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Determine mode and execute
    if args.interactive:
        # Interactive mode
        interactive_mode(model=args.model)

    elif args.input:
        # Batch mode from CSV
        if not args.input.exists():
            logger.error(f"Input file not found: {args.input}")
            sys.exit(1)

        results = triage_from_csv(
            input_file=args.input,
            gene=args.gene,
            output_file=args.output,
            model=args.model
        )

        # Print summary
        if not args.output:
            print("\nResults:")
            print(json.dumps(results, indent=2))

    elif args.title and args.abstract:
        # Single paper mode
        result = triage_single_paper(
            title=args.title,
            abstract=args.abstract,
            gene=args.gene,
            pmid=args.pmid,
            model=args.model
        )

        # Save or print
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(result, f, indent=2)
            logger.info(f"Result saved to {args.output}")
        else:
            print("\nTriage Result:")
            print(json.dumps(result, indent=2))

    else:
        # No valid mode selected
        parser.print_help()
        print("\nError: Please specify either --interactive, --input, or --title/--abstract")
        sys.exit(1)


if __name__ == '__main__':
    main()
