"""
Example usage of the Tiered Biomedical Extraction Pipeline.

This script demonstrates how to:
1. Run the full pipeline for a gene symbol
2. Customize filters and models
3. Access and export results
"""

import json
import os
from dotenv import load_dotenv
from pipeline import BiomedicalExtractionPipeline, run_pipeline_for_gene
from pipeline.filters import KeywordFilter, InternFilter
from pipeline.extraction import ExpertExtractor

# Load environment variables from .env file
load_dotenv()

# Set your API keys as environment variables before running:
# export OPENAI_API_KEY="your-api-key-here"
# or
# export ANTHROPIC_API_KEY="your-api-key-here"


def example_1_basic_usage():
    """Example 1: Basic usage with default settings."""
    print("\n" + "="*80)
    print("EXAMPLE 1: Basic Pipeline Usage")
    print("="*80 + "\n")

    # Simple one-liner to run the pipeline
    results = run_pipeline_for_gene(
        gene_symbol="BRCA1",
        email="your_email@example.com",  # Replace with your email
        max_papers=10  # Limit to 10 papers for testing
    )

    print(f"\nProcessed {results['stats'].total_papers_sourced} papers")
    print(f"Successful extractions: {results['stats'].total_extraction_successes}")
    print(f"Total variants found: {results['stats'].total_variants_extracted}")


def example_2_custom_configuration():
    """Example 2: Customized pipeline with specific models and filters."""
    print("\n" + "="*80)
    print("EXAMPLE 2: Custom Configuration")
    print("="*80 + "\n")

    # Create custom keyword filter with stricter requirements
    keyword_filter = KeywordFilter(
        min_keyword_matches=3  # Require at least 3 keyword matches
    )

    # Use Claude instead of GPT-4 for extraction (optional)
    # expert_extractor = ExpertExtractor(model="claude-3-sonnet-20240229")

    # Or use GPT-4o for extraction
    expert_extractor = ExpertExtractor(model="gpt-4o")

    # Create pipeline with custom components
    pipeline = BiomedicalExtractionPipeline(
        email="your_email@example.com",
        keyword_filter=keyword_filter,
        expert_extractor=expert_extractor
    )

    # Run for a specific gene
    results = pipeline.run(
        gene_symbol="TP53",
        max_papers=5
    )

    print(f"\nCost savings from tiered filtering: ${results['stats'].estimated_cost_savings:.2f}")


def example_3_export_results():
    """Example 3: Export results to JSON file."""
    print("\n" + "="*80)
    print("EXAMPLE 3: Export Results")
    print("="*80 + "\n")

    # Run pipeline
    results = run_pipeline_for_gene(
        gene_symbol="SCN5A",
        email="your_email@example.com",
        max_papers=3
    )

    # Extract successful extractions
    successful_extractions = results['successful_extractions']

    # Create output directory
    os.makedirs("output", exist_ok=True)

    # Export each extraction to JSON
    for extraction in successful_extractions:
        if extraction.extracted_data:
            filename = f"output/SCN5A_PMID_{extraction.pmid}.json"

            with open(filename, 'w') as f:
                json.dump(extraction.extracted_data, f, indent=2)

            print(f"✓ Exported: {filename}")

    print(f"\nExported {len(successful_extractions)} extractions to output/ directory")


def example_4_detailed_filtering_logs():
    """Example 4: Examine detailed filtering decisions."""
    print("\n" + "="*80)
    print("EXAMPLE 4: Detailed Filtering Analysis")
    print("="*80 + "\n")

    pipeline = BiomedicalExtractionPipeline(email="your_email@example.com")

    results = pipeline.run(
        gene_symbol="CFTR",
        max_papers=5
    )

    # Analyze what happened at each tier
    print("\nDetailed Filtering Analysis:")
    print("-" * 80)

    for pipeline_result in results['results']:
        print(f"\nPMID: {pipeline_result.pmid}")
        print(f"Final tier reached: {pipeline_result.final_tier_reached.value}")

        for filter_result in pipeline_result.filter_results:
            decision_emoji = "✓" if filter_result.decision.value == "pass" else "✗"
            print(f"  {decision_emoji} {filter_result.tier.value}: {filter_result.reason}")

        if pipeline_result.extraction_result:
            if pipeline_result.extraction_result.success:
                variants = (
                    pipeline_result.extraction_result.extracted_data
                    .get("extraction_metadata", {})
                    .get("total_variants_found", 0)
                )
                print(f"  ✓ Extraction: {variants} variants found")
            else:
                print(f"  ✗ Extraction failed: {pipeline_result.extraction_result.error}")


def example_5_disable_expensive_tier():
    """Example 5: Run only the filtering tiers (skip expensive extraction)."""
    print("\n" + "="*80)
    print("EXAMPLE 5: Filtering Only (No Extraction)")
    print("="*80 + "\n")

    # Disable Tier 3 to only run filtering (useful for testing)
    pipeline = BiomedicalExtractionPipeline(
        email="your_email@example.com",
        enable_tier3=False  # Skip expensive extraction
    )

    results = pipeline.run(
        gene_symbol="KCNQ1",
        max_papers=20  # Can process more papers since we're not extracting
    )

    print("\nFiltering completed without extraction:")
    print(f"Papers that would proceed to extraction: {results['stats'].passed_tier2_intern}")
    print(f"Papers filtered out (cost saved): {results['stats'].failed_tier1 + results['stats'].failed_tier2}")


if __name__ == "__main__":
    # Check if API keys are set
    if not (os.getenv("OPENAI_API_KEY") or os.getenv("ANTHROPIC_API_KEY")):
        print("⚠️  WARNING: No API keys found!")
        print("Please set OPENAI_API_KEY or ANTHROPIC_API_KEY environment variable.")
        print("\nExample:")
        print("  export OPENAI_API_KEY='your-key-here'")
        print("\nYou can still run the examples, but LLM-based tiers will fail.")
        print()

    # Run examples (uncomment the ones you want to try)

    # example_1_basic_usage()
    # example_2_custom_configuration()
    # example_3_export_results()
    # example_4_detailed_filtering_logs()
    # example_5_disable_expensive_tier()

    print("\n" + "="*80)
    print("To run examples, uncomment the function calls in __main__")
    print("="*80 + "\n")
