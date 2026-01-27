#!/usr/bin/env python3
"""
Test script for pedigree extraction functionality.

Usage:
    # Test on existing PDF file
    python scripts/test_pedigree_extraction.py --pdf /path/to/supplement.pdf

    # Test on existing figures directory
    python scripts/test_pedigree_extraction.py --figures /path/to/pmid_figures/

    # Harvest a specific PMID and test extraction
    python scripts/test_pedigree_extraction.py --pmid 12345678 --gene KCNH2
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from config.settings import get_settings
from harvesting.format_converters import FormatConverter
from pipeline.pedigree_extractor import PedigreeExtractor


def test_pdf_extraction(pdf_path: Path, output_dir: Path):
    """Test image extraction from a PDF."""
    print(f"\n{'='*60}")
    print(f"Testing PDF image extraction")
    print(f"PDF: {pdf_path}")
    print(f"Output: {output_dir}")
    print(f"{'='*60}\n")

    converter = FormatConverter()

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract with images
    markdown, images = converter.pdf_to_markdown_with_images(
        pdf_path,
        output_dir=output_dir,
        extract_images=True,
    )

    print(f"Extracted text: {len(markdown)} characters")
    print(f"Extracted images: {len(images)}")

    if images:
        print("\nImages extracted:")
        for img in images:
            print(
                f"  - {img['filename']}: {img['width']}x{img['height']} ({img['size_bytes']} bytes)"
            )
            if "path" in img:
                print(f"    Saved to: {img['path']}")

    return images


def test_pedigree_detection(figures_dir: Path, detect_only: bool = False):
    """Test pedigree detection on a figures directory."""
    print(f"\n{'='*60}")
    print(f"Testing pedigree detection")
    print(f"Figures directory: {figures_dir}")
    print(f"Mode: {'detect only' if detect_only else 'full extraction'}")
    print(f"{'='*60}\n")

    settings = get_settings()
    print(f"Vision model: {settings.vision_model}")
    print(f"Confidence threshold: {settings.pedigree_confidence_threshold}")
    print()

    extractor = PedigreeExtractor(
        model=settings.vision_model,
        detection_confidence_threshold=settings.pedigree_confidence_threshold,
    )

    results = extractor.process_figures_directory(figures_dir, detect_only=detect_only)

    print(f"\nResults: {len(results)} pedigree(s) found")

    if results:
        for i, result in enumerate(results, 1):
            print(f"\n--- Pedigree {i}: {result.get('image')} ---")
            print(f"Confidence: {result.get('confidence', 0):.0%}")

            if not detect_only and "individuals" in result:
                individuals = result["individuals"]
                print(f"Individuals: {len(individuals)}")
                print(f"Inheritance: {result.get('inheritance_pattern', 'unknown')}")

                affected = sum(
                    1
                    for ind in individuals
                    if ind.get("affection_status") == "affected"
                )
                carriers = sum(1 for ind in individuals if ind.get("is_carrier"))
                print(f"Affected: {affected}, Carriers: {carriers}")

                # Print individual details
                print("\nIndividuals:")
                for ind in individuals:
                    status = ind.get("affection_status", "?")
                    carrier = " (carrier)" if ind.get("is_carrier") else ""
                    proband = " *PROBAND*" if ind.get("is_proband") else ""
                    print(
                        f"  {ind.get('id', '?')}: {ind.get('sex', '?')}, {status}{carrier}{proband}"
                    )

        # Generate summary
        print("\n" + "=" * 60)
        print("Generated summary for extraction:")
        print("=" * 60)
        summary = extractor.summarize_for_extraction(results)
        print(summary)

    return results


def test_full_pipeline(pmid: str, gene: str, output_dir: Path):
    """Test the full pipeline on a specific PMID."""
    print(f"\n{'='*60}")
    print(f"Testing full pipeline")
    print(f"PMID: {pmid}")
    print(f"Gene: {gene}")
    print(f"Output: {output_dir}")
    print(f"{'='*60}\n")

    from harvesting.orchestrator import PMCHarvester

    harvester = PMCHarvester(output_dir=str(output_dir), gene_symbol=gene)
    success, result = harvester.process_pmid(pmid)

    if success:
        print(f"\n✅ Harvest successful: {result}")

        # Check for extracted figures
        figures_dir = output_dir / f"{pmid}_figures"
        if figures_dir.exists():
            images = list(figures_dir.glob("*"))
            print(f"Figures extracted: {len(images)}")

            # Check for pedigree results
            pedigree_json = output_dir / f"{pmid}_PEDIGREES.json"
            if pedigree_json.exists():
                import json

                with open(pedigree_json) as f:
                    pedigrees = json.load(f)
                print(f"Pedigrees found: {len(pedigrees)}")
            else:
                print("No pedigrees detected")
        else:
            print("No figures were extracted")
    else:
        print(f"\n❌ Harvest failed: {result}")

    return success


def main():
    parser = argparse.ArgumentParser(
        description="Test pedigree extraction functionality"
    )
    parser.add_argument("--pdf", type=Path, help="Test image extraction on a PDF file")
    parser.add_argument(
        "--figures", type=Path, help="Test pedigree detection on a figures directory"
    )
    parser.add_argument("--pmid", type=str, help="Harvest and test a specific PMID")
    parser.add_argument(
        "--gene", type=str, default="KCNH2", help="Gene symbol (default: KCNH2)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("test_pedigree_output"),
        help="Output directory",
    )
    parser.add_argument(
        "--detect-only",
        action="store_true",
        help="Only detect pedigrees, don't extract",
    )

    args = parser.parse_args()

    if args.pdf:
        if not args.pdf.exists():
            print(f"Error: PDF file not found: {args.pdf}")
            return 1

        figures_dir = args.output / "figures"
        images = test_pdf_extraction(args.pdf, figures_dir)

        if images:
            print("\n" + "-" * 60)
            test_pedigree_detection(figures_dir, detect_only=args.detect_only)

    elif args.figures:
        if not args.figures.exists():
            print(f"Error: Figures directory not found: {args.figures}")
            return 1

        test_pedigree_detection(args.figures, detect_only=args.detect_only)

    elif args.pmid:
        test_full_pipeline(args.pmid, args.gene, args.output)

    else:
        parser.print_help()
        print("\nExamples:")
        print("  python scripts/test_pedigree_extraction.py --pdf supplement.pdf")
        print(
            "  python scripts/test_pedigree_extraction.py --figures /path/to/figures/"
        )
        print(
            "  python scripts/test_pedigree_extraction.py --pmid 12345678 --gene KCNH2"
        )
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
