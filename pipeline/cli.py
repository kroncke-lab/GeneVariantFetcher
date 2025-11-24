"""
CLI Entry Point for the Gene Variant Fetcher Pipeline.

Provides Typer-based command-line interface for all pipeline operations.
"""

import typer
from pathlib import Path
from typing import Optional

from .sourcing import PaperSourcer
from .harvesting import PMCHarvester
from .aggregation import DataAggregator
from config.settings import get_settings

app = typer.Typer(help="Gene Variant Fetcher pipeline CLI")


@app.command()
def source(
    query: str = typer.Argument(..., help="Gene symbol or search query"),
    limit: int = typer.Option(50, help="Max records to fetch"),
    output: Optional[Path] = typer.Option(None, help="Output file for PMIDs (default: {query}_pmids.txt)")
):
    """Run the sourcing stage and persist candidate papers."""
    typer.echo(f"Sourcing papers for query: {query}")

    sourcer = PaperSourcer()
    results = sourcer.fetch_papers(query, max_results_per_source=limit)

    # Save to file
    if not output:
        output = Path(f"{query}_pmids.txt")

    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w') as f:
        for pmid in results:
            f.write(f"{pmid}\n")

    typer.echo(f"Sourced {len(results)} papers for query '{query}'")
    typer.echo(f"PMIDs saved to: {output}")


@app.command()
def harvest(
    input_path: str = typer.Argument(..., help="Path to file containing PMIDs"),
    output_path: str = typer.Option("data/harvests", help="Destination directory")
):
    """Harvest full text for previously sourced papers."""
    typer.echo(f"Harvesting papers from: {input_path}")

    harvester = PMCHarvester()
    count = harvester.run(input_path=input_path, output_path=output_path)

    typer.echo(f"Harvested {count} articles into {output_path}")


@app.command()
def aggregate(
    extraction_dir: Path = typer.Argument(..., help="Directory containing extraction JSON files"),
    gene_symbol: str = typer.Argument(..., help="Gene symbol (e.g., BRCA1, TTR)"),
    output: Optional[Path] = typer.Option(None, help="Output file path")
):
    """Aggregate penetrance data from extraction results."""
    if not output:
        output = extraction_dir.parent / f"{gene_symbol}_penetrance_summary.json"

    typer.echo(f"Aggregating penetrance data for {gene_symbol} from {extraction_dir}")

    aggregator = DataAggregator()
    summary = aggregator.aggregate_from_directory(extraction_dir, gene_symbol, output)

    typer.echo(f"Aggregation complete!")
    typer.echo(f"Total variants: {summary['total_variants']}")
    typer.echo(f"Validation errors: {summary['validation']['error_count']}")
    typer.echo(f"Validation warnings: {summary['validation']['warning_count']}")
    typer.echo(f"Output saved to: {output}")


@app.command()
def run(
    gene_symbol: str = typer.Argument(..., help="Gene symbol to process"),
    max_pmids: int = typer.Option(50, help="Maximum PMIDs to source"),
    max_downloads: int = typer.Option(20, help="Maximum papers to download"),
    skip_harvest: bool = typer.Option(False, help="Skip harvesting step"),
    skip_aggregate: bool = typer.Option(False, help="Skip aggregation step")
):
    """Run the complete pipeline for a gene symbol."""
    typer.echo(f"🚀 Starting complete pipeline for gene: {gene_symbol}")
    typer.echo(f"Max PMIDs: {max_pmids}, Max downloads: {max_downloads}")

    # Step 1: Source PMIDs
    typer.echo("\n📚 Step 1: Sourcing papers...")
    sourcer = PaperSourcer()
    pmids = sourcer.fetch_papers(gene_symbol, max_results_per_source=max_pmids)

    if not pmids:
        typer.echo("❌ No PMIDs found, exiting")
        return

    # Save PMIDs
    pmids_file = Path(f"{gene_symbol}_pmids.txt")
    with open(pmids_file, 'w') as f:
        for pmid in pmids[:max_downloads]:  # Limit to max_downloads
            f.write(f"{pmid}\n")

    typer.echo(f"✓ Sourced {len(pmids)} PMIDs, processing {min(len(pmids), max_downloads)}")

    if skip_harvest:
        typer.echo("⏭️ Skipping harvest step")
        return

    # Step 2: Harvest full text
    typer.echo("\n📥 Step 2: Harvesting full text...")
    harvester = PMCHarvester()
    harvested_count = harvester.run(str(pmids_file))

    typer.echo(f"✓ Harvested {harvested_count} papers")

    if skip_aggregate:
        typer.echo("⏭️ Skipping aggregation step")
        return

    # Step 3: Aggregate (if extraction results exist)
    extraction_dir = Path("data/harvests")
    if extraction_dir.exists():
        typer.echo("\n📊 Step 3: Aggregating results...")
        aggregate(extraction_dir, gene_symbol)
    else:
        typer.echo("⏭️ No extraction directory found, skipping aggregation")

    typer.echo(f"\n✅ Pipeline complete for {gene_symbol}!")


@app.command()
def config():
    """Show current configuration settings."""
    settings = get_settings()

    typer.echo("🔧 Current Configuration:")
    typer.echo(f"  OpenAI API Key: {'✓ Set' if settings.openai_api_key else '✗ Not set'}")
    typer.echo(f"  Anthropic API Key: {'✓ Set' if settings.anthropic_api_key else '✗ Not set'}")
    typer.echo(f"  NCBI Email: {settings.ncbi_email}")
    typer.echo(f"  NCBI API Key: {'✓ Set' if settings.ncbi_api_key else '✗ Not set'}")
    typer.echo(f"  Intern Model: {settings.intern_model}")
    typer.echo(f"  Extractor Model: {settings.extractor_model}")
    typer.echo(f"  Rate Limit: {settings.rate_limit_per_minute} per minute")


if __name__ == "__main__":
    app()