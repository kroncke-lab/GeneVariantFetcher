"""
Shared workflow step implementations.

This module contains the core logic for each pipeline step, used by both
the CLI (automated_workflow.py) and GUI (worker.py) implementations.

Each step function takes explicit parameters and returns a result dict,
allowing the callers to handle their own context (checkpoints, logging, etc.).
"""

import csv
import json
import logging
import os
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# =============================================================================
# Step Result Types
# =============================================================================


@dataclass
class StepResult:
    """Result from a workflow step."""

    success: bool
    stats: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None
    data: Dict[str, Any] = field(default_factory=dict)


# =============================================================================
# Step 0: Discover Synonyms
# =============================================================================


def discover_synonyms(
    gene_symbol: str,
    email: str,
    existing_synonyms: Optional[List[str]] = None,
    api_key: Optional[str] = None,
) -> StepResult:
    """
    Discover gene synonyms from NCBI Gene database.

    Args:
        gene_symbol: Primary gene symbol
        email: Email for NCBI API
        existing_synonyms: Already known synonyms to merge with
        api_key: Optional NCBI API key for higher rate limits

    Returns:
        StepResult with synonyms in data["synonyms"]
    """
    from gene_literature.synonym_finder import (
        SynonymFinder,
        automatic_synonym_selection,
    )

    all_synonyms = list(existing_synonyms) if existing_synonyms else []

    try:
        synonym_finder = SynonymFinder(
            email=email,
            api_key=api_key or os.getenv("NCBI_API_KEY"),
        )

        found_synonyms = synonym_finder.find_gene_synonyms(
            gene_symbol,
            include_other_designations=False,
        )

        auto_selected = automatic_synonym_selection(
            gene_symbol,
            found_synonyms,
            include_official=True,
            include_aliases=True,
            include_other_designations=False,
            only_relevant=False,
        )

        # Merge with existing synonyms (avoid duplicates)
        existing_set = set(s.lower() for s in all_synonyms)
        for syn in auto_selected:
            if syn.lower() not in existing_set and syn.lower() != gene_symbol.lower():
                all_synonyms.append(syn)
                existing_set.add(syn.lower())

        return StepResult(
            success=True,
            stats={
                "synonyms_found": len(auto_selected),
                "total_synonyms": len(all_synonyms),
            },
            data={"synonyms": all_synonyms},
        )

    except Exception as e:
        logger.warning(f"Failed to discover synonyms: {e}")
        return StepResult(
            success=False,
            error=str(e),
            data={"synonyms": all_synonyms},
        )


# =============================================================================
# Step 1: Fetch PMIDs
# =============================================================================


def fetch_pmids(
    gene_symbol: str,
    email: str,
    output_path: Path,
    max_results: int = 100,
    synonyms: Optional[List[str]] = None,
    use_pubmind: bool = True,
    use_pubmed: bool = True,
    use_europepmc: bool = False,
    api_key: Optional[str] = None,
) -> StepResult:
    """
    Fetch PMIDs from literature sources.

    Args:
        gene_symbol: Gene to search for
        email: Email for NCBI API
        output_path: Directory to save PMID files
        max_results: Maximum PMIDs to fetch
        synonyms: Gene synonyms to include in search
        use_pubmind: Enable PubMind source
        use_pubmed: Enable PubMed source
        use_europepmc: Enable Europe PMC source
        api_key: Optional NCBI API key

    Returns:
        StepResult with PMIDs in data["pmids"]
    """
    from config.settings import get_settings
    from gene_literature.discovery import discover_pmids_for_gene

    settings = get_settings()
    settings.use_pubmind = use_pubmind
    settings.use_pubmed = use_pubmed
    settings.use_europepmc = use_europepmc

    pubmind_file = output_path / f"{gene_symbol}_pmids_pubmind.txt"
    pubmed_file = output_path / f"{gene_symbol}_pmids_pubmed.txt"
    combined_file = output_path / f"{gene_symbol}_pmids.txt"

    try:
        pmid_discovery = discover_pmids_for_gene(
            gene_symbol=gene_symbol,
            email=email,
            max_results=max_results,
            pubmind_output=pubmind_file,
            pubmed_output=pubmed_file,
            combined_output=combined_file,
            api_key=api_key or os.getenv("NCBI_API_KEY"),
            settings=settings,
            synonyms=synonyms,
        )

        return StepResult(
            success=True,
            stats={
                "pubmind_count": len(pmid_discovery.pubmind_pmids),
                "pubmed_count": len(pmid_discovery.pubmed_pmids),
                "europepmc_count": len(pmid_discovery.europepmc_pmids),
                "total_unique": len(pmid_discovery.combined_pmids),
            },
            data={"pmids": list(pmid_discovery.combined_pmids)},
        )

    except Exception as e:
        logger.error(f"PMID discovery failed: {e}")
        return StepResult(success=False, error=str(e), data={"pmids": []})


# =============================================================================
# Step 2: Fetch Abstracts
# =============================================================================


def fetch_abstracts(
    pmids: List[str],
    output_path: Path,
    email: str,
) -> StepResult:
    """
    Fetch and save abstracts for PMIDs.

    Args:
        pmids: List of PMIDs to fetch
        output_path: Base output directory
        email: Email for NCBI API

    Returns:
        StepResult with abstract records in data["abstract_records"]
    """
    from harvesting.abstracts import fetch_and_save_abstracts

    abstract_dir = output_path / "abstract_json"

    try:
        abstract_records = fetch_and_save_abstracts(
            pmids=pmids,
            output_dir=str(abstract_dir),
            email=email,
        )

        return StepResult(
            success=True,
            stats={"abstracts_fetched": len(abstract_records)},
            data={"abstract_records": abstract_records, "abstract_dir": abstract_dir},
        )

    except Exception as e:
        logger.error(f"Abstract fetch failed: {e}")
        return StepResult(success=False, error=str(e))


# =============================================================================
# Step 3: Filter Papers
# =============================================================================


def filter_papers(
    pmids: List[str],
    abstract_records: Dict[str, str],
    gene_symbol: str,
    output_path: Path,
    enable_tier1: bool = True,
    enable_tier2: bool = True,
    use_clinical_triage: bool = False,
    tier1_min_keywords: int = 1,
    tier2_confidence_threshold: float = 0.5,
) -> StepResult:
    """
    Filter papers using tiered filtering.

    Args:
        pmids: List of PMIDs to filter
        abstract_records: Dict mapping PMID to abstract JSON path
        gene_symbol: Gene being searched
        output_path: Base output directory
        enable_tier1: Enable keyword filtering
        enable_tier2: Enable LLM filtering
        use_clinical_triage: Use ClinicalDataTriageFilter for Tier 2
        tier1_min_keywords: Minimum keywords for Tier 1
        tier2_confidence_threshold: Confidence threshold for Tier 2

    Returns:
        StepResult with filtered PMIDs in data["filtered_pmids"]
    """
    from pipeline.filters import ClinicalDataTriageFilter, InternFilter, KeywordFilter
    from utils.models import FilterDecision, FilterResult, FilterTier, Paper

    keyword_filter = KeywordFilter(min_keyword_matches=tier1_min_keywords)
    tier2_filter = (
        ClinicalDataTriageFilter()
        if use_clinical_triage
        else InternFilter(confidence_threshold=tier2_confidence_threshold)
    )

    # Persistent, resume-safe filtering output
    pmid_status_dir = output_path / "pmid_status"
    pmid_status_dir.mkdir(parents=True, exist_ok=True)

    progress_file = pmid_status_dir / "filter_progress.jsonl"
    filtered_pmids_file = pmid_status_dir / "filtered_pmids.txt"
    filtered_out_file = pmid_status_dir / "filtered_out.csv"

    filtered_pmids: List[str] = []
    dropped_pmids: List[Tuple[str, str]] = []

    processed_pmids: set[str] = set()

    # Resume: load prior progress if present
    if progress_file.exists():
        try:
            with open(progress_file, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except Exception:
                        continue

                    pmid = str(rec.get("pmid", ""))
                    if not pmid:
                        continue

                    processed_pmids.add(pmid)
                    decision = str(rec.get("final_decision", "")).upper()
                    reason = str(rec.get("reason", ""))
                    if decision == "PASS":
                        filtered_pmids.append(pmid)
                    elif decision == "FAIL":
                        dropped_pmids.append((pmid, reason))
        except Exception as e:
            logger.warning(f"Failed to load filter progress ({progress_file}): {e}")

    # Append-only progress writer (flush each record so SIGKILL/OOM doesn't lose work)
    progress_fh = open(progress_file, "a", encoding="utf-8")

    def _write_progress(record: Dict[str, Any]) -> None:
        progress_fh.write(json.dumps(record, ensure_ascii=False) + "\n")
        progress_fh.flush()

    for pmid in pmids:
        if pmid in processed_pmids:
            continue
        record_path = abstract_records.get(pmid)

        if not record_path or not Path(record_path).exists():
            reason = "Missing abstract JSON"
            dropped_pmids.append((pmid, reason))
            _write_progress({
                "pmid": pmid,
                "final_decision": "FAIL",
                "reason": reason,
                "stage": "abstract_load",
            })
            continue

        try:
            with open(record_path, "r", encoding="utf-8") as f:
                record = json.load(f)
        except Exception as e:
            reason = f"Abstract JSON read error: {e}"
            dropped_pmids.append((pmid, reason))
            _write_progress({
                "pmid": pmid,
                "final_decision": "FAIL",
                "reason": reason,
                "stage": "abstract_load",
            })
            continue

        metadata = record.get("metadata", {})
        paper = Paper(
            pmid=pmid,
            title=metadata.get("title"),
            abstract=record.get("abstract"),
            authors=metadata.get("authors"),
            journal=metadata.get("journal"),
            publication_date=metadata.get("year"),
            gene_symbol=gene_symbol,
            source="PubMed",
        )

        # Tier 1: Keyword filter
        tier1_result = None
        if enable_tier1:
            tier1_result = keyword_filter.filter(paper)
            if tier1_result.decision is not FilterDecision.PASS:
                reason = tier1_result.reason
                dropped_pmids.append((pmid, reason))
                _write_progress({
                    "pmid": pmid,
                    "final_decision": "FAIL",
                    "reason": reason,
                    "stage": "tier1",
                    "tier1": {
                        "decision": tier1_result.decision.value,
                        "reason": tier1_result.reason,
                        "confidence": tier1_result.confidence,
                    },
                })
                continue

        # Tier 2: LLM filter
        if enable_tier2:
            if use_clinical_triage:
                triage_result = tier2_filter.triage_paper(paper, gene_symbol)
                decision = (
                    FilterDecision.PASS
                    if triage_result.get("decision") == "KEEP"
                    else FilterDecision.FAIL
                )
                confidence = triage_result.get("confidence")
                reason = triage_result.get("reason", "No reason provided")

                if (
                    decision is FilterDecision.PASS
                    and confidence is not None
                    and confidence < tier2_confidence_threshold
                ):
                    decision = FilterDecision.FAIL
                    reason = f"Low confidence ({confidence:.2f}): {reason}"

                tier2_result = FilterResult(
                    decision=decision,
                    tier=FilterTier.TIER_2_INTERN,
                    reason=reason,
                    pmid=pmid,
                    confidence=confidence,
                )
            else:
                tier2_result = tier2_filter.filter(paper)

            if tier2_result.decision is not FilterDecision.PASS:
                reason = tier2_result.reason
                dropped_pmids.append((pmid, reason))
                _write_progress({
                    "pmid": pmid,
                    "final_decision": "FAIL",
                    "reason": reason,
                    "stage": "tier2",
                    "tier1": (
                        {
                            "decision": tier1_result.decision.value,
                            "reason": tier1_result.reason,
                            "confidence": tier1_result.confidence,
                        }
                        if tier1_result is not None
                        else None
                    ),
                    "tier2": {
                        "decision": tier2_result.decision.value,
                        "reason": tier2_result.reason,
                        "confidence": tier2_result.confidence,
                    },
                })
                continue

        filtered_pmids.append(pmid)
        _write_progress({
            "pmid": pmid,
            "final_decision": "PASS",
            "reason": "passed",
            "stage": "done",
            "tier1": (
                {
                    "decision": tier1_result.decision.value,
                    "reason": tier1_result.reason,
                    "confidence": tier1_result.confidence,
                }
                if tier1_result is not None
                else None
            ),
            "tier2": (
                {
                    "decision": tier2_result.decision.value,
                    "reason": tier2_result.reason,
                    "confidence": tier2_result.confidence,
                }
                if enable_tier2
                else None
            ),
        })

    # Close progress file handle (important when running long jobs)
    try:
        progress_fh.close()
    except Exception:
        pass

    # Save final dropped PMIDs + passed PMIDs
    with open(filtered_out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["PMID", "Reason"])
        writer.writerows(dropped_pmids)

    with open(filtered_pmids_file, "w", encoding="utf-8") as f:
        for pmid in filtered_pmids:
            f.write(str(pmid) + "\n")

    return StepResult(
        success=True,
        stats={
            "passed_filter": len(filtered_pmids),
            "dropped": len(dropped_pmids),
            "resumed_processed": len(processed_pmids),
        },
        data={
            "filtered_pmids": filtered_pmids,
            "dropped_pmids": dropped_pmids,
            "progress_file": str(progress_file),
            "filtered_pmids_file": str(filtered_pmids_file),
            "filtered_out_file": str(filtered_out_file),
        },
    )


# =============================================================================
# Step 4: Download Full-Text
# =============================================================================


def _consolidate_prior_downloads(
    pmids: List[str],
    harvest_dir: Path,
    output_base: Path,
) -> set:
    """
    Pre-download step: scan prior run directories for already-acquired
    FULL_CONTEXT.md files and copy them into the current harvest directory.

    This avoids re-downloading papers that were successfully fetched in
    previous runs (PMC, publisher APIs, browser downloads, etc.).

    Args:
        pmids: List of PMIDs we need content for
        harvest_dir: Current run's harvest directory (pmc_fulltext/)
        output_base: The gvf_output base directory (parent of gene dirs)

    Returns:
        Set of PMIDs that were recovered from prior runs
    """
    import shutil

    pmid_set = set(str(p) for p in pmids)
    recovered = set()

    # Already present in current harvest dir
    already_present = set()
    for f in harvest_dir.glob("*_FULL_CONTEXT.md"):
        already_present.add(f.name.replace("_FULL_CONTEXT.md", ""))

    needed = pmid_set - already_present
    if not needed:
        return recovered

    # Scan all known prior output directories for FULL_CONTEXT.md files
    # Walk the entire output base to find any prior downloads
    prior_files = {}  # pmid -> (path, size)
    for root, dirs, files in os.walk(str(output_base)):
        # Skip the current harvest dir to avoid self-references
        if str(harvest_dir) in root:
            continue
        for f in files:
            if f.endswith("_FULL_CONTEXT.md"):
                pmid = f.replace("_FULL_CONTEXT.md", "")
                if pmid in needed:
                    full_path = Path(root) / f
                    size = full_path.stat().st_size
                    # Keep the largest version (most complete)
                    if pmid not in prior_files or size > prior_files[pmid][1]:
                        prior_files[pmid] = (full_path, size)

    # Copy recovered files and their associated directories (supplements, figures)
    for pmid, (src_path, size) in prior_files.items():
        if size < 500:
            continue  # Skip near-empty files

        dst_path = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
        try:
            shutil.copy2(str(src_path), str(dst_path))
            recovered.add(pmid)

            # Also copy supplements and figures directories if they exist
            src_dir = src_path.parent
            for suffix in ["_supplements", "_figures"]:
                src_assoc = src_dir / f"{pmid}{suffix}"
                dst_assoc = harvest_dir / f"{pmid}{suffix}"
                if src_assoc.is_dir() and not dst_assoc.exists():
                    shutil.copytree(str(src_assoc), str(dst_assoc))

        except Exception as e:
            logger.warning(f"Failed to recover PMID {pmid} from {src_path}: {e}")

    if recovered:
        logger.info(
            f"✓ Recovered {len(recovered)} papers from prior runs "
            f"(skipping re-download)"
        )
        print(
            f"\n📂 Pre-download consolidation: recovered {len(recovered)} papers "
            f"from prior runs"
        )

    return recovered


def download_fulltext(
    pmids: List[str],
    output_path: Path,
    gene_symbol: str,
    max_papers: Optional[int] = None,
    delay: float = 2.0,
    prior_output_base: Optional[Path] = None,
) -> StepResult:
    """
    Download full-text papers from PMC, with prior-run consolidation.

    Before downloading, scans prior run directories for already-acquired
    FULL_CONTEXT.md files and copies them in. Only attempts fresh downloads
    for truly missing PMIDs.

    Args:
        pmids: List of PMIDs to download
        output_path: Base output directory
        gene_symbol: Gene symbol for naming
        max_papers: Maximum papers to download (None = all)
        delay: Delay between downloads in seconds
        prior_output_base: Base directory to scan for prior downloads.
            If None, uses output_path's grandparent (typically gvf_output/).

    Returns:
        StepResult with download stats
    """
    import pandas as pd

    from harvesting import PMCHarvester

    harvest_dir = output_path / "pmc_fulltext"
    harvester = PMCHarvester(output_dir=str(harvest_dir), gene_symbol=gene_symbol)

    pmids_to_download = pmids[:max_papers] if max_papers else pmids

    # --- Pre-download consolidation: recover papers from prior runs ---
    if prior_output_base is None:
        # Default: walk up to the gvf_output directory
        # output_path is typically gvf_output/GENE/TIMESTAMP/
        # We want gvf_output/ as the scan root
        prior_output_base = output_path.parent.parent
        # Sanity check - make sure we're not scanning from /
        if len(str(prior_output_base)) < 10:
            prior_output_base = output_path

    harvest_dir.mkdir(parents=True, exist_ok=True)
    recovered_pmids = _consolidate_prior_downloads(
        pmids_to_download, harvest_dir, prior_output_base
    )

    # Only send un-recovered PMIDs to the harvester
    remaining_pmids = [
        p for p in pmids_to_download if str(p) not in recovered_pmids
    ]
    logger.info(
        f"Download plan: {len(recovered_pmids)} recovered, "
        f"{len(remaining_pmids)} to download fresh"
    )

    if remaining_pmids:
        harvester.harvest(remaining_pmids, delay=delay)

    # Check results — include both freshly downloaded and recovered from prior runs
    success_log = harvest_dir / "successful_downloads.csv"
    successfully_downloaded = set(recovered_pmids)  # Start with recovered

    if success_log.exists():
        try:
            df = pd.read_csv(success_log)
            successfully_downloaded.update(str(p) for p in df["PMID"].tolist())
        except Exception as e:
            logger.warning(f"Could not read success log: {e}")

    # Get paywalled/failed PMIDs
    paywalled_log = harvest_dir / "paywalled_missing.csv"
    paywalled_pmids = {}

    if paywalled_log.exists():
        try:
            df = pd.read_csv(paywalled_log)
            for _, row in df.iterrows():
                paywalled_pmids[str(row["PMID"])] = row.get("Reason", "Unknown")
        except Exception:
            pass

    # Identify abstract-only PMIDs
    abstract_only = [
        p for p in pmids_to_download if str(p) not in successfully_downloaded
    ]

    return StepResult(
        success=True,
        stats={
            "attempted": len(pmids_to_download),
            "downloaded": len(successfully_downloaded),
            "recovered_from_prior": len(recovered_pmids),
            "freshly_downloaded": len(successfully_downloaded) - len(recovered_pmids),
            "paywalled": len(paywalled_pmids),
            "abstract_only": len(abstract_only),
        },
        data={
            "harvest_dir": harvest_dir,
            "downloaded_pmids": list(successfully_downloaded),
            "recovered_pmids": list(recovered_pmids),
            "abstract_only_pmids": abstract_only,
            "paywalled_pmids": paywalled_pmids,
        },
    )


# =============================================================================
# Step 5: Run Data Scout
# =============================================================================


def preprocess_papers(
    harvest_dir: Path,
    gene_symbol: str,
    abstracts_dir: Optional[Path] = None,
) -> StepResult:
    """
    Deterministic pre-processing step. Runs BEFORE Data Scout or LLM extraction.
    
    - Strips XML/HTML noise from FULL_CONTEXT.md files
    - Removes reference sections, copyright boilerplate
    - Injects PubMed abstracts as guaranteed baseline
    - Classifies files (full_text, abstract_only, empty, etc.)
    
    Modifies files IN-PLACE to minimize token usage in downstream LLM steps.
    No API calls — pure regex/string operations.
    """
    from gvf.preprocessor import PaperPreprocessor
    
    if not harvest_dir.exists():
        return StepResult(
            success=True,
            stats={"skipped": True, "reason": "no_harvest_dir"},
        )
    
    # Auto-detect abstracts directory
    if abstracts_dir is None:
        # Look for pubmed_abstracts in parent/sibling dirs
        candidates = [
            harvest_dir.parent / "pubmed_abstracts",
            harvest_dir.parent.parent / "pubmed_abstracts",
            harvest_dir.parent.parent / gene_symbol / "pubmed_abstracts",
        ]
        for c in candidates:
            if c.exists():
                abstracts_dir = c
                break
    
    preprocessor = PaperPreprocessor(
        abstracts_dir=str(abstracts_dir) if abstracts_dir else None
    )
    
    files = list(harvest_dir.glob("*_FULL_CONTEXT.md"))
    if not files:
        return StepResult(
            success=True,
            stats={"skipped": True, "reason": "no_files"},
        )
    
    processed = 0
    classifications = {}
    total_original = 0
    total_cleaned = 0
    
    for f in files:
        try:
            pmid = f.name.replace("_FULL_CONTEXT.md", "").replace(f"KCNH2_PMID_", "").replace(f"{gene_symbol}_PMID_", "")
            text = f.read_text(encoding="utf-8", errors="replace")
            total_original += len(text)
            
            classification = preprocessor.classify(text)
            classifications[classification] = classifications.get(classification, 0) + 1
            
            cleaned = preprocessor.clean(text)
            cleaned = preprocessor.inject_abstract(cleaned, pmid)
            total_cleaned += len(cleaned)
            
            # Write back in-place
            f.write_text(cleaned, encoding="utf-8")
            processed += 1
        except Exception as e:
            logger.warning(f"Preprocess error for {f.name}: {e}")
    
    token_savings_pct = round((1 - total_cleaned / max(1, total_original)) * 100, 1) if total_original > total_cleaned else 0
    
    logger.info(f"✓ Preprocessed {processed} papers")
    logger.info(f"  Classifications: {classifications}")
    if token_savings_pct > 0:
        logger.info(f"  Token savings: ~{token_savings_pct}% reduction in text size")
    
    return StepResult(
        success=True,
        stats={
            "processed": processed,
            "classifications": classifications,
            "original_bytes": total_original,
            "cleaned_bytes": total_cleaned,
            "token_savings_pct": token_savings_pct,
        },
    )


def run_data_scout(
    harvest_dir: Path,
    gene_symbol: str,
    min_relevance: float = 0.3,
) -> StepResult:
    """
    Run Data Scout to create DATA_ZONES files.

    Args:
        harvest_dir: Directory with FULL_CONTEXT.md files
        gene_symbol: Gene symbol for context
        min_relevance: Minimum relevance threshold

    Returns:
        StepResult with scouting stats
    """
    from config.settings import get_settings
    from pipeline.data_scout import GeneticDataScout

    if not harvest_dir.exists():
        return StepResult(
            success=True,
            stats={"skipped": True, "reason": "no_fulltext_dir"},
        )

    settings = get_settings()

    # Find files needing scouting
    full_context_files = list(harvest_dir.glob("*_FULL_CONTEXT.md"))
    existing_zones = set(
        f.name.replace("_DATA_ZONES.md", "")
        for f in harvest_dir.glob("*_DATA_ZONES.md")
    )

    files_to_scout = [
        f
        for f in full_context_files
        if f.name.replace("_FULL_CONTEXT.md", "") not in existing_zones
    ]

    if not files_to_scout:
        return StepResult(
            success=True,
            stats={"skipped": True, "already_scouted": len(existing_zones)},
        )

    scout = GeneticDataScout(
        gene_symbol=gene_symbol,
        min_relevance_score=min_relevance,
        max_zones=settings.scout_max_zones if settings else 30,
    )

    scouted = 0
    errors = 0

    for md_file in files_to_scout:
        try:
            content = md_file.read_text(encoding="utf-8")
            # Extract PMID from filename (e.g., "12345678_FULL_CONTEXT.md" -> "12345678")
            pmid = md_file.name.replace("_FULL_CONTEXT.md", "")
            report = scout.scan(content, pmid=pmid)

            if report.zones_kept > 0:
                condensed = scout.format_markdown(report, content)
                output_file = md_file.with_name(
                    md_file.name.replace("_FULL_CONTEXT.md", "_DATA_ZONES.md")
                )
                output_file.write_text(condensed, encoding="utf-8")
                scouted += 1
        except Exception as e:
            logger.warning(f"Error scouting {md_file.name}: {e}")
            errors += 1

    return StepResult(
        success=True,
        stats={"scouted": scouted, "errors": errors, "skipped": len(existing_zones)},
    )


# =============================================================================
# Step 6: Extract Variants
# =============================================================================


def extract_variants(
    harvest_dir: Path,
    extraction_dir: Path,
    gene_symbol: str,
    abstract_records: Optional[Dict[str, str]] = None,
    abstract_only_pmids: Optional[List[str]] = None,
    tier_threshold: int = 1,
    max_workers: int = 8,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> StepResult:
    """
    Extract variant data from papers.

    Args:
        harvest_dir: Directory with markdown files
        extraction_dir: Directory to save extractions
        gene_symbol: Gene symbol
        abstract_records: Dict mapping PMID to abstract JSON path
        abstract_only_pmids: PMIDs that need abstract-only extraction
        tier_threshold: Model cascade threshold
        max_workers: Max parallel workers
        progress_callback: Optional callback(completed, total)

    Returns:
        StepResult with extraction stats
    """
    from pipeline.extraction import ExpertExtractor
    from utils.models import Paper
    from utils.pmid_utils import extract_pmid_from_filename

    extraction_dir.mkdir(exist_ok=True)

    # Find markdown files (prefer DATA_ZONES over FULL_CONTEXT)
    data_zones = {
        f.name.replace("_DATA_ZONES.md", ""): f
        for f in harvest_dir.glob("*_DATA_ZONES.md")
    }
    full_context = {
        f.name.replace("_FULL_CONTEXT.md", ""): f
        for f in harvest_dir.glob("*_FULL_CONTEXT.md")
    }

    markdown_files = []
    for pmid in set(data_zones.keys()) | set(full_context.keys()):
        if pmid in data_zones:
            markdown_files.append(data_zones[pmid])
        elif pmid in full_context:
            markdown_files.append(full_context[pmid])

    # Prepare abstract-only papers
    abstract_papers = []
    if abstract_only_pmids and abstract_records:
        for pmid in abstract_only_pmids:
            record_path = abstract_records.get(pmid)
            if record_path and Path(record_path).exists():
                abstract_papers.append((pmid, record_path))

    extractor_local = threading.local()

    def get_extractor() -> ExpertExtractor:
        """Get a thread-local extractor to avoid shared mutable model state."""
        extractor = getattr(extractor_local, "instance", None)
        if extractor is None:
            extractor = ExpertExtractor(
                tier_threshold=tier_threshold, fulltext_dir=str(harvest_dir)
            )
            extractor_local.instance = extractor
        return extractor

    extractions = []
    failures = []
    total = len(markdown_files) + len(abstract_papers)
    completed = 0

    def process_fulltext(md_file: Path) -> Tuple[Any, Optional[Tuple[str, str]]]:
        pmid = extract_pmid_from_filename(md_file)
        if not pmid:
            return (None, (md_file.name, "Could not extract PMID"))

        try:
            content = md_file.read_text(encoding="utf-8")
            paper = Paper(
                pmid=pmid,
                title=f"Paper {pmid}",
                full_text=content,
                gene_symbol=gene_symbol,
            )

            result = get_extractor().extract(paper)

            if result.success:
                output_file = extraction_dir / f"{gene_symbol}_PMID_{pmid}.json"
                with open(output_file, "w") as f:
                    json.dump(result.extracted_data, f, indent=2)
                return (result, None)
            else:
                return (None, (pmid, result.error))

        except Exception as e:
            return (None, (pmid, str(e)))

    def process_abstract(
        pmid: str, record_path: str
    ) -> Tuple[Any, Optional[Tuple[str, str]]]:
        try:
            with open(record_path, "r", encoding="utf-8") as f:
                record = json.load(f)

            abstract_text = record.get("abstract")
            if not abstract_text:
                return (None, (pmid, "No abstract available"))

            metadata = record.get("metadata", {})
            paper = Paper(
                pmid=pmid,
                title=metadata.get("title", f"Paper {pmid}"),
                abstract=abstract_text,
                gene_symbol=gene_symbol,
            )

            result = get_extractor().extract(paper)

            if result.success:
                if result.extracted_data:
                    if "extraction_metadata" not in result.extracted_data:
                        result.extracted_data["extraction_metadata"] = {}
                    result.extracted_data["extraction_metadata"]["abstract_only"] = True

                output_file = extraction_dir / f"{gene_symbol}_PMID_{pmid}.json"
                with open(output_file, "w") as f:
                    json.dump(result.extracted_data, f, indent=2)
                return (result, None)
            else:
                return (None, (pmid, result.error))

        except Exception as e:
            return (None, (pmid, str(e)))

    # Process in parallel
    workers = min(max_workers, total) if total > 0 else 1

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {}

        for md_file in markdown_files:
            futures[executor.submit(process_fulltext, md_file)] = ("fulltext", md_file)

        for pmid, record_path in abstract_papers:
            futures[executor.submit(process_abstract, pmid, record_path)] = (
                "abstract",
                pmid,
            )

        for future in as_completed(futures):
            completed += 1
            if progress_callback:
                progress_callback(completed, total)

            try:
                result, failure = future.result()
                if result:
                    extractions.append(result)
                if failure:
                    failures.append(failure)
            except Exception as e:
                source_type, source_info = futures[future]
                pmid = (
                    extract_pmid_from_filename(source_info)
                    if source_type == "fulltext"
                    else source_info
                )
                failures.append((pmid, str(e)))

    # Save failures
    if failures:
        failures_file = (
            extraction_dir.parent / "pmid_status" / "extraction_failures.csv"
        )
        failures_file.parent.mkdir(parents=True, exist_ok=True)
        with open(failures_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["PMID", "Error"])
            writer.writerows(failures)

    # Count variants
    total_variants = sum(
        e.extracted_data.get("extraction_metadata", {}).get("total_variants_found", 0)
        for e in extractions
        if e.extracted_data
    )

    return StepResult(
        success=True,
        stats={
            "papers_extracted": len(extractions),
            "failures": len(failures),
            "total_variants": total_variants,
        },
        data={
            "extractions": extractions,
            "failures": failures,
        },
    )


# =============================================================================
# Step 7: Aggregate Data
# =============================================================================


def aggregate_data(
    extraction_dir: Path,
    gene_symbol: str,
    output_path: Path,
) -> StepResult:
    """
    Aggregate penetrance data across extractions.

    Args:
        extraction_dir: Directory with extraction JSONs
        gene_symbol: Gene symbol
        output_path: Base output directory for summary file

    Returns:
        StepResult with aggregation stats
    """
    from pipeline.aggregation import aggregate_penetrance

    summary_file = output_path / f"{gene_symbol}_penetrance_summary.json"

    try:
        penetrance_summary = aggregate_penetrance(
            extraction_dir=extraction_dir,
            gene_symbol=gene_symbol,
            output_file=summary_file,
        )

        return StepResult(
            success=True,
            stats={
                "variants_aggregated": penetrance_summary.get("total_variants", 0),
            },
            data={"summary": penetrance_summary},
        )

    except Exception as e:
        logger.error(f"Aggregation failed: {e}")
        return StepResult(success=False, error=str(e))


# =============================================================================
# Step 8: Migrate to SQLite
# =============================================================================


def migrate_to_sqlite(
    extraction_dir: Path,
    db_path: Path,
) -> StepResult:
    """
    Migrate extraction JSONs to SQLite database.

    Args:
        extraction_dir: Directory with extraction JSONs
        db_path: Path for output database

    Returns:
        StepResult with migration stats
    """
    from harvesting.migrate_to_sqlite import (
        create_database_schema,
        migrate_extraction_directory,
    )

    try:
        conn = create_database_schema(str(db_path))
        stats = migrate_extraction_directory(conn, extraction_dir)
        conn.close()

        return StepResult(
            success=True,
            stats={
                "total_files": stats["total_files"],
                "successful": stats["successful"],
                "failed": stats["failed"],
            },
        )

    except Exception as e:
        logger.error(f"Migration failed: {e}")
        return StepResult(success=False, error=str(e))
