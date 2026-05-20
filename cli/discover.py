"""
Discover-only command for GeneVariantFetcher.

Runs the upstream half of the pipeline:
  1. (Optional) Synonym discovery (NCBI Gene)
  2. PubMed/PubMind/EuropePMC PMID search (with optional --disease clause)
  3. Abstract + metadata fetch
  4. Tier-1 keyword filter + Tier-2 LLM filter

Outputs a single CSV per gene with one row per PMID:
  pmid, title, abstract_snippet, relevance_score, data_type, reason, year, journal

Does NOT download full text or extract variants. Reuses the existing
pipeline/steps.py building blocks so behavior matches the extract command's
upstream stages exactly.
"""

from __future__ import annotations

import csv
import json
import logging
import os
from pathlib import Path
from typing import Optional

import typer
from typing_extensions import Annotated

from utils.bootstrap import has_llm_provider_key, initialize_runtime
from utils.logging_utils import get_logger, setup_logging

logger = get_logger(__name__)


def _classify_data_type(reason: str) -> str:
    """Best-effort label for the kind of data the Tier-2 LLM said the paper has.

    Pure heuristic over the LLM's reason string — used only for the discovery
    CSV. Not load-bearing.
    """
    if not reason:
        return "unknown"
    r = reason.lower()
    has_patient = any(
        token in r
        for token in (
            "patient",
            "case report",
            "case series",
            "cohort",
            "proband",
            "family",
            "clinical",
        )
    )
    has_functional = any(
        token in r
        for token in (
            "functional",
            "in vitro",
            "patch clamp",
            "electrophysiolog",
            "expression",
            "biophysical",
        )
    )
    if has_patient and has_functional:
        return "both"
    if has_patient:
        return "patient"
    if has_functional:
        return "functional"
    return "other"


def _abstract_snippet(abstract: Optional[str], max_chars: int = 240) -> str:
    if not abstract:
        return ""
    text = " ".join(abstract.split())
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 1].rstrip() + "…"


def run_discover(
    gene_symbol: str,
    email: str,
    output_dir: Path,
    max_pmids: int = 200,
    auto_synonyms: bool = True,
    synonyms: Optional[list[str]] = None,
    disease: Optional[str] = None,
    use_clinical_triage: bool = False,
) -> Path:
    """Run discover-only pipeline and write a per-gene CSV.

    Returns the path to the written CSV.
    """
    initialize_runtime()

    from config.settings import get_settings
    from pipeline.steps import (
        discover_synonyms,
        fetch_abstracts,
        fetch_pmids,
        filter_papers,
    )

    settings = get_settings()

    output_dir = Path(output_dir)
    gene_dir = output_dir / gene_symbol
    gene_dir.mkdir(parents=True, exist_ok=True)

    log_file = gene_dir / f"{gene_symbol}_discover.log"
    setup_logging(level=logging.INFO, log_file=log_file)
    logger.info("Logging to %s", log_file)

    logger.info("=" * 80)
    logger.info(
        "DISCOVER ONLY — gene=%s disease=%s max_pmids=%d",
        gene_symbol,
        disease,
        max_pmids,
    )
    logger.info("=" * 80)

    all_synonyms: list[str] = list(synonyms) if synonyms else []
    if auto_synonyms:
        logger.info("Discovering synonyms from NCBI Gene…")
        syn_result = discover_synonyms(
            gene_symbol=gene_symbol,
            email=email,
            existing_synonyms=all_synonyms,
            api_key=os.getenv("NCBI_API_KEY"),
        )
        all_synonyms = syn_result.data.get("synonyms", all_synonyms)
        logger.info("Synonyms in use: %s", all_synonyms or "(none)")

    pmid_result = fetch_pmids(
        gene_symbol=gene_symbol,
        email=email,
        output_path=gene_dir,
        max_results=max_pmids,
        synonyms=all_synonyms or None,
        use_pubmind=settings.use_pubmind,
        use_pubmed=settings.use_pubmed and not settings.pubmind_only,
        use_europepmc=settings.use_europepmc,
        api_key=os.getenv("NCBI_API_KEY"),
        disease=disease,
    )
    if not pmid_result.success:
        raise RuntimeError(f"PMID discovery failed: {pmid_result.error}")

    pmids: list[str] = list(pmid_result.data.get("pmids", []))
    logger.info("Discovered %d unique PMIDs", len(pmids))
    if not pmids:
        csv_path = gene_dir / f"{gene_symbol}_discover.csv"
        csv_path.write_text(
            "pmid,title,abstract_snippet,relevance_score,data_type,reason,year,journal\n",
            encoding="utf-8",
        )
        logger.warning("No PMIDs found — wrote empty CSV at %s", csv_path)
        return csv_path

    abs_result = fetch_abstracts(pmids=pmids, output_path=gene_dir, email=email)
    if not abs_result.success:
        raise RuntimeError(f"Abstract fetch failed: {abs_result.error}")
    abstract_records = abs_result.data.get("abstract_records", {})

    filter_result = filter_papers(
        pmids=pmids,
        abstract_records=abstract_records,
        gene_symbol=gene_symbol,
        output_path=gene_dir,
        enable_tier1=settings.enable_tier1,
        enable_tier2=settings.enable_tier2,
        use_clinical_triage=use_clinical_triage,
        tier1_min_keywords=settings.tier1_min_keywords,
        tier2_confidence_threshold=settings.tier2_confidence_threshold,
        disease=disease,
    )
    if not filter_result.success:
        raise RuntimeError(f"Filtering failed: {filter_result.error}")

    filtered_pmids = set(filter_result.data.get("filtered_pmids", []))
    dropped = dict(filter_result.data.get("dropped_pmids", []))

    # Pull per-PMID Tier-2 reason/confidence from the progress jsonl that
    # filter_papers writes alongside its other outputs.
    progress_file = gene_dir / "pmid_status" / "filter_progress.jsonl"
    per_pmid_progress: dict[str, dict] = {}
    if progress_file.exists():
        with progress_file.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                pmid = str(rec.get("pmid", ""))
                if pmid:
                    per_pmid_progress[pmid] = rec

    csv_path = gene_dir / f"{gene_symbol}_discover.csv"
    rows = 0
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "pmid",
                "passed_filters",
                "title",
                "abstract_snippet",
                "relevance_score",
                "data_type",
                "reason",
                "year",
                "journal",
            ]
        )
        for pmid in pmids:
            record_path = abstract_records.get(pmid)
            title = ""
            abstract = ""
            year = ""
            journal = ""
            if record_path and Path(record_path).exists():
                try:
                    with open(record_path, "r", encoding="utf-8") as f:
                        rec = json.load(f)
                    md = rec.get("metadata", {})
                    title = md.get("title") or ""
                    abstract = rec.get("abstract") or ""
                    year = md.get("year") or ""
                    journal = md.get("journal") or ""
                except Exception as exc:  # noqa: BLE001
                    logger.debug("Could not read abstract record for %s: %s", pmid, exc)

            progress = per_pmid_progress.get(pmid, {})
            tier2 = progress.get("tier2") or {}
            confidence = tier2.get("confidence")
            reason = tier2.get("reason") or progress.get("reason") or ""
            if not reason and pmid in dropped:
                reason = dropped[pmid]

            writer.writerow(
                [
                    pmid,
                    "yes" if pmid in filtered_pmids else "no",
                    title,
                    _abstract_snippet(abstract),
                    f"{confidence:.2f}" if isinstance(confidence, (int, float)) else "",
                    _classify_data_type(reason),
                    reason,
                    year,
                    journal,
                ]
            )
            rows += 1

    logger.info(
        "Wrote %d rows to %s (%d passed filters)",
        rows,
        csv_path,
        len(filtered_pmids),
    )
    return csv_path


# Wired into cli/__init__.py
def discover_command(
    gene: Annotated[str, typer.Argument(help="Gene symbol (e.g. NPPA, KCNH2)")],
    email: Annotated[
        str, typer.Option("--email", "-e", help="Your email for NCBI E-utilities")
    ],
    output: Annotated[
        str,
        typer.Option(
            "--output",
            "-o",
            help="Output directory; CSV is written to <output>/<gene>/<gene>_discover.csv",
        ),
    ],
    max_pmids: Annotated[
        int, typer.Option("--max-pmids", help="Maximum PMIDs to fetch")
    ] = 200,
    disease: Annotated[
        str,
        typer.Option(
            "--disease",
            "-d",
            help="Optional disease term to scope discovery (e.g. 'atrial fibrillation')",
        ),
    ] = None,
    auto_synonyms: Annotated[
        bool,
        typer.Option(
            "--auto-synonyms/--no-auto-synonyms",
            help="Discover gene synonyms from NCBI Gene (default: on)",
        ),
    ] = True,
    synonyms: Annotated[
        list[str],
        typer.Option(
            "--synonym",
            help="Manually specify gene synonym (repeatable)",
        ),
    ] = None,
    clinical_triage: Annotated[
        bool,
        typer.Option(
            "--clinical-triage",
            help="Use ClinicalDataTriageFilter for Tier 2 (instead of InternFilter)",
        ),
    ] = False,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable verbose logging")
    ] = False,
):
    """Run the upstream discovery half of the pipeline (no full-text, no extraction)."""
    initialize_runtime()

    if verbose:
        setup_logging(level=logging.DEBUG)
    else:
        setup_logging(level=logging.INFO)

    if not has_llm_provider_key():
        typer.echo("⚠️  ERROR: No LLM provider API key found in environment!", err=True)
        typer.echo(
            "Set OPENAI_API_KEY, AZURE_AI_API_KEY, or ANTHROPIC_API_KEY in your .env file.",
            err=True,
        )
        raise typer.Exit(1)

    try:
        csv_path = run_discover(
            gene_symbol=gene,
            email=email,
            output_dir=Path(output),
            max_pmids=max_pmids,
            auto_synonyms=auto_synonyms,
            synonyms=synonyms or None,
            disease=disease,
            use_clinical_triage=clinical_triage,
        )
        typer.echo(f"✓ Wrote discovery CSV: {csv_path}")
    except KeyboardInterrupt:
        typer.echo("\n⚠️  Discover interrupted by user", err=True)
        raise typer.Exit(1)
    except Exception as exc:
        typer.echo(f"\n❌ Discover failed: {exc}", err=True)
        raise typer.Exit(1)
