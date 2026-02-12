"""
CLI tools for GeneVariantFetcher.

This package contains standalone command-line utilities:
- automated_workflow: Full pipeline CLI
- fetch_manager: Semi-manual download helper for paywalled papers
- browser_fetch: Browser-automated download with Cloudflare/CAPTCHA handling
- scout: Standalone Data Scout for identifying high-value data zones
"""

import os
from pathlib import Path

import typer
from typing_extensions import Annotated

from cli.automated_workflow import automated_variant_extraction_workflow
from cli.audit_paywalls import run_paywall_audit
from cli.scout import run_scout

__all__ = ["automated_variant_extraction_workflow", "run_scout", "run_paywall_audit", "app"]

app = typer.Typer(
    help="GeneVariantFetcher CLI - Tools for extracting genetic variant data from literature"
)


@app.command("extract")
def extract_command(
    gene: Annotated[str, typer.Argument(help="Gene symbol (e.g., BRCA1, SCN5A, TP53)")],
    email: Annotated[
        str, typer.Option("--email", "-e", help="Your email for NCBI E-utilities")
    ],
    output: Annotated[
        str,
        typer.Option(
            "--output", "-o", help="Output directory for all data and analyses"
        ),
    ],
    max_pmids: Annotated[
        int, typer.Option("--max-pmids", help="Maximum PMIDs to fetch")
    ] = 100,
    max_downloads: Annotated[
        int, typer.Option("--max-downloads", help="Maximum papers to download")
    ] = 50,
    tier_threshold: Annotated[
        int,
        typer.Option(
            "--tier-threshold",
            help="Threshold for trying next model if first finds fewer variants",
        ),
    ] = 1,
    clinical_triage: Annotated[
        bool,
        typer.Option(
            "--clinical-triage",
            help="Use ClinicalDataTriageFilter for Tier 2 filtering",
        ),
    ] = False,
    auto_synonyms: Annotated[
        bool,
        typer.Option(
            "--auto-synonyms/--no-auto-synonyms",
            help="Automatically discover and use gene synonyms (default: on)",
        ),
    ] = True,
    synonyms: Annotated[
        list[str],
        typer.Option(
            "--synonym",
            help="Manually specify gene synonym (can be used multiple times)",
        ),
    ] = None,
    scout_first: Annotated[
        bool,
        typer.Option(
            "--scout-first",
            help="Run Data Scout before extraction to identify high-value data zones for better context",
        ),
    ] = False,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable verbose logging")
    ] = False,
):
    """Run the complete automated workflow from gene symbol to extracted variant data."""
    import logging

    from utils.logging_utils import setup_logging

    if verbose:
        setup_logging(level=logging.DEBUG)
    else:
        setup_logging(level=logging.INFO)

    # Check for API keys
    if not os.getenv("OPENAI_API_KEY"):
        typer.echo("⚠️  ERROR: OPENAI_API_KEY not found in environment!", err=True)
        typer.echo("Please set OPENAI_API_KEY in your .env file", err=True)
        raise typer.Exit(1)

    try:
        automated_variant_extraction_workflow(
            gene_symbol=gene,
            email=email,
            output_dir=output,
            max_pmids=max_pmids,
            max_papers_to_download=max_downloads,
            tier_threshold=tier_threshold,
            use_clinical_triage=clinical_triage,
            auto_synonyms=auto_synonyms,
            synonyms=synonyms or [],
            scout_first=scout_first,
        )
    except KeyboardInterrupt:
        typer.echo("\n⚠️  Workflow interrupted by user", err=True)
        raise typer.Exit(1)
    except Exception as e:
        typer.echo(f"\n❌ Workflow failed with error: {e}", err=True)
        raise typer.Exit(1)


@app.command("audit-paywalls")
def audit_paywalls_command(
    harvest_dir: Annotated[
        str,
        typer.Option(
            "--harvest-dir",
            help="Harvest directory containing paywalled_missing.csv (typically .../pmc_fulltext)",
        ),
    ],
    out_dir: Annotated[
        str,
        typer.Option(
            "--out-dir",
            help="Output directory for audit report/json",
        ),
    ],
    limit: Annotated[
        int,
        typer.Option("--limit", help="Optional limit for quick tests"),
    ] = 0,
):
    """Audit paywalled_missing.csv to distinguish true paywalls from captcha/blocked or missed OA."""
    paywalled_csv = Path(harvest_dir) / "paywalled_missing.csv"
    out_base = Path(out_dir)
    result = run_paywall_audit(
        paywalled_csv=paywalled_csv,
        out_json=out_base / "paywall_audit.json",
        out_md=out_base / "paywall_audit.md",
        limit=limit or None,
    )
    typer.echo(f"Wrote: {result['out_md']}")
    typer.echo(f"Wrote: {result['out_json']}")


@app.command("scout")
def scout_command(
    input_path: Annotated[
        str,
        typer.Argument(
            help="Input directory with downloaded papers or manifest.json file"
        ),
    ],
    output_dir: Annotated[
        str, typer.Argument(help="Output directory for scout results")
    ],
    gene: Annotated[str, typer.Argument(help="Gene symbol (e.g., BRCA1, SCN5A)")],
    manifest_out: Annotated[
        str, typer.Option("--manifest-out", help="Custom path for output manifest")
    ] = None,
    min_relevance: Annotated[
        float,
        typer.Option(
            "--min-relevance", help="Minimum relevance score for TEXT zones (0.0-1.0)"
        ),
    ] = 0.1,
    max_zones: Annotated[
        int, typer.Option("--max-zones", help="Maximum zones per paper")
    ] = 30,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable debug logging")
    ] = False,
):
    """Run the Data Scout stage to identify high-value data zones in papers."""
    import logging

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)

    try:
        manifest = run_scout(
            input_path=input_path,
            output_dir=output_dir,
            gene=gene,
            manifest_out=manifest_out,
            min_relevance=min_relevance,
            max_zones=max_zones,
        )

        # Exit with error if all failed
        if hasattr(manifest, "entries"):
            from utils.manifest import Status

            if all(e.status != Status.SUCCESS for e in manifest.entries):
                raise typer.Exit(1)

    except Exception as e:
        typer.echo(f"Fatal error: {e}", err=True)
        raise typer.Exit(1)
