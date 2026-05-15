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
from typing import Optional

import typer
from typing_extensions import Annotated

from cli.automated_workflow import automated_variant_extraction_workflow
from cli.audit_paywalls import run_paywall_audit
from cli.discover import discover_command, run_discover
from cli.extract import run_extraction
from cli.reharvest import reharvest_command
from cli.scout import run_scout

__all__ = [
    "automated_variant_extraction_workflow",
    "run_scout",
    "run_paywall_audit",
    "run_extraction",
    "run_discover",
    "reharvest_command",
    "app",
]

app = typer.Typer(
    help="GeneVariantFetcher CLI - Tools for extracting genetic variant data from literature"
)


def _read_pmid_file(path: Path) -> list[str]:
    """Read PMIDs from a file, one per line. '#' comments and blank lines ignored."""
    if not path.exists():
        raise typer.BadParameter(f"PMID list file not found: {path}")
    tokens: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Allow trailing comments after the PMID (e.g. "12345 # paper title")
        token = line.split("#", 1)[0].strip()
        if token:
            tokens.append(token)
    return tokens


def _parse_pmid_arg(arg: str) -> list[str]:
    """Parse a --pmids argument: '@path/to/file.txt' or 'pmid1,pmid2,...'.

    Returns a deduplicated list preserving first-seen order. PMIDs are validated
    as all-digit tokens; non-digit tokens are dropped with no warning so that
    file headers and stray whitespace pass through cleanly.
    """
    if arg.startswith("@"):
        raw_tokens = _read_pmid_file(Path(arg[1:]).expanduser())
    else:
        raw_tokens = [tok.strip() for tok in arg.split(",")]

    seen: set[str] = set()
    pmids: list[str] = []
    for tok in raw_tokens:
        if tok.isdigit() and tok not in seen:
            seen.add(tok)
            pmids.append(tok)
    return pmids


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
    ] = 1500,
    max_downloads: Annotated[
        Optional[int],
        typer.Option(
            "--max-downloads",
            help="Maximum papers to download (default: no limit)",
        ),
    ] = None,
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
    disease: Annotated[
        str,
        typer.Option(
            "--disease",
            "-d",
            help=(
                "Optional disease term (e.g. 'atrial fibrillation'). When set, "
                "disease clause is appended to PubMed queries and Tier-2 filter "
                "prompt prioritizes original patient/functional data. When omitted, "
                "behavior is unchanged."
            ),
        ),
    ] = None,
    pmids: Annotated[
        str,
        typer.Option(
            "--pmids",
            help=(
                "Skip PubMed/PubMind/EuropePMC discovery AND Tier 1/2 filtering, "
                "running the pipeline on an explicit PMID list. Accepts a "
                "comma-separated list (e.g. 12345,67890) or a file path "
                "prefixed with '@' (e.g. @comparison_results/gold_standard_pmids.txt). "
                "Lines starting with '#' and blank lines are ignored in files."
            ),
        ),
    ] = None,
    pmid_file: Annotated[
        Path,
        typer.Option(
            "--pmid-file",
            help=(
                "Same effect as '--pmids @file' but as a dedicated flag for "
                "clarity. One PMID per line; lines starting with '#' and blank "
                "lines are ignored. Trailing '# comment' is allowed on PMID "
                "lines. May be combined with --pmids; the union (deduped) is used."
            ),
        ),
    ] = None,
    browser_html_fallback: Annotated[
        Optional[bool],
        typer.Option(
            "--browser-html-fallback/--no-browser-html-fallback",
            help=(
                "Toggle Tier 3.5 browser-based HTML fallback (Playwright) for"
                " papers that fail the publisher-API path but are freely"
                " readable post-embargo. ON by default; pass"
                " --no-browser-html-fallback to disable. When omitted, the"
                " ENABLE_BROWSER_HTML_FALLBACK env var / settings default wins."
            ),
        ),
    ] = None,
    model_provider: Annotated[
        str,
        typer.Option(
            "--model-provider",
            help=(
                "LLM provider for all tiers: 'anthropic' (default), 'azure', or"
                " 'openai'. Pass 'azure' to fall back to Azure AI Foundry"
                " deployments. Per-tier env vars (TIER2_MODEL, TIER3_MODELS,"
                " TABLE_ROUTER_MODEL, VISION_MODEL) still win."
            ),
        ),
    ] = None,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable verbose logging")
    ] = False,
):
    """Run the complete automated workflow from gene symbol to extracted variant data."""
    import logging

    from utils.logging_utils import setup_logging

    if browser_html_fallback is not None:
        os.environ["ENABLE_BROWSER_HTML_FALLBACK"] = (
            "true" if browser_html_fallback else "false"
        )
        # Clear cached settings so the new env var is picked up.
        try:
            from config.settings import get_settings as _get_settings

            _get_settings.cache_clear()
        except Exception:
            pass

    if verbose:
        setup_logging(level=logging.DEBUG)
    else:
        setup_logging(level=logging.INFO)

    # If --model-provider was supplied, push it into the env BEFORE the
    # workflow runs so settings.get_*_model() picks the right tier defaults.
    # Reset any cached Settings so the next get_settings() re-reads env.
    if model_provider:
        provider_normalized = model_provider.strip().lower()
        valid_providers = {"azure", "anthropic", "openai"}
        if provider_normalized not in valid_providers:
            typer.echo(
                f"⚠️  ERROR: --model-provider must be one of {sorted(valid_providers)},"
                f" got '{model_provider}'",
                err=True,
            )
            raise typer.Exit(1)
        os.environ["MODEL_PROVIDER"] = provider_normalized
        from config.settings import reset_settings_cache

        reset_settings_cache()

    # Require at least one LLM provider key. ANTHROPIC_API_KEY, OPENAI_API_KEY,
    # or AZURE_AI_API_KEY all satisfy the pipeline since LiteLLM resolves
    # credentials per-call based on the configured TIER*_MODEL strings.
    selected_provider = (
        (model_provider.strip().lower() if model_provider else None)
        or os.getenv("MODEL_PROVIDER", "").strip().lower()
        or "anthropic"
    )
    if selected_provider == "anthropic":
        if not os.getenv("ANTHROPIC_API_KEY"):
            typer.echo(
                "⚠️  ERROR: --model-provider=anthropic requires ANTHROPIC_API_KEY in your .env file.",
                err=True,
            )
            raise typer.Exit(1)
    elif not (
        os.getenv("OPENAI_API_KEY")
        or os.getenv("AZURE_AI_API_KEY")
        or os.getenv("ANTHROPIC_API_KEY")
    ):
        typer.echo("⚠️  ERROR: No LLM provider API key found in environment!", err=True)
        typer.echo(
            "Set OPENAI_API_KEY, AZURE_AI_API_KEY, or ANTHROPIC_API_KEY in your .env file.",
            err=True,
        )
        raise typer.Exit(1)

    explicit_pmids: list[str] | None = None
    if pmids or pmid_file:
        merged: list[str] = []
        seen: set[str] = set()
        if pmids:
            for p in _parse_pmid_arg(pmids):
                if p not in seen:
                    seen.add(p)
                    merged.append(p)
        if pmid_file:
            for tok in _read_pmid_file(pmid_file.expanduser()):
                if tok.isdigit() and tok not in seen:
                    seen.add(tok)
                    merged.append(tok)
        explicit_pmids = merged
        if not explicit_pmids:
            sources = []
            if pmids:
                sources.append(f"--pmids={pmids!r}")
            if pmid_file:
                sources.append(f"--pmid-file={str(pmid_file)!r}")
            typer.echo(
                f"⚠️  ERROR: explicit PMID list is empty after parsing "
                f"{', '.join(sources)}",
                err=True,
            )
            raise typer.Exit(1)
        typer.echo(
            f"📋 Using {len(explicit_pmids)} explicit PMIDs; discovery and "
            f"Tier 1/2 filtering will be skipped."
        )

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
            disease=disease,
            pmids=explicit_pmids,
        )
    except KeyboardInterrupt:
        typer.echo("\n⚠️  Workflow interrupted by user", err=True)
        raise typer.Exit(1)
    except Exception as e:
        typer.echo(f"\n❌ Workflow failed with error: {e}", err=True)
        raise typer.Exit(1)


app.command("discover")(discover_command)
app.command("reharvest")(reharvest_command)


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


@app.command("extract-folder")
def extract_folder_command(
    input_path: Annotated[
        str,
        typer.Argument(help="Input directory with scout outputs or manifest.json file"),
    ],
    output_dir: Annotated[
        str, typer.Argument(help="Output directory for extraction JSON files")
    ],
    gene: Annotated[str, typer.Argument(help="Gene symbol (e.g., KCNH2, BRCA1)")],
    manifest_out: Annotated[
        str,
        typer.Option("--manifest-out", help="Custom path for output manifest"),
    ] = None,
    model: Annotated[
        list[str],
        typer.Option(
            "--model",
            help="Model for extraction (repeatable for tiered fallback)",
        ),
    ] = None,
    tier_threshold: Annotated[
        int,
        typer.Option(
            "--tier-threshold",
            help="If first model finds fewer variants, try next model",
        ),
    ] = 1,
    full_text: Annotated[
        bool,
        typer.Option(
            "--full-text",
            help="Use FULL_CONTEXT.md files instead of DATA_ZONES.md",
        ),
    ] = False,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable debug logging")
    ] = False,
):
    """Extract variants from pre-downloaded papers (DATA_ZONES or FULL_CONTEXT files)."""
    import logging

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)

    try:
        from cli.extract import ValidationError

        manifest = run_extraction(
            input_path=Path(input_path),
            output_dir=Path(output_dir),
            gene=gene,
            manifest_out=Path(manifest_out) if manifest_out else None,
            models=model or None,
            tier_threshold=tier_threshold,
            use_full_text=full_text,
        )

        # Exit with error if all failed
        if manifest.entries:
            from utils.manifest import Status

            if all(e.status != Status.SUCCESS for e in manifest.entries):
                raise typer.Exit(1)

    except ValidationError as e:
        typer.echo(f"Input validation failed: {e}", err=True)
        raise typer.Exit(2)
    except Exception as e:
        typer.echo(f"Fatal error: {e}", err=True)
        raise typer.Exit(1)
