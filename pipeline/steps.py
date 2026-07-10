"""
Shared workflow step implementations.

This module contains the core logic for each pipeline step used by CLI
orchestration.

Each step function takes explicit parameters and returns a result dict,
allowing the callers to handle their own context (checkpoints, logging, etc.).
"""

import csv
import hashlib
import json
import logging
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from config.constants import DEFAULT_MAX_WORKERS
from pipeline.source_quality import is_usable_fulltext_source

logger = logging.getLogger(__name__)


_NON_RETRYABLE_EXTRACTION_FAILURES = (
    "could not extract pmid",
    "no abstract available",
    "source is abstract-only fallback",
)


def _resolve_default_workers() -> int:
    """Pick the provider-aware worker default (Anthropic 10, Azure 3).

    Falls back to DEFAULT_MAX_WORKERS if Settings can't load (e.g. tests that
    construct minimal env). Wrapped in a function so each call reflects the
    current MODEL_PROVIDER — important when the CLI flips provider at runtime.
    """
    try:
        from config.settings import get_settings

        return get_settings().get_max_workers()
    except Exception:
        return DEFAULT_MAX_WORKERS


def _is_retryable_extraction_failure(error: Any) -> bool:
    """Return whether a failed extraction is worth one more model pass.

    Most extraction failures are provider/rate-limit/model exceptions and are
    cheap to retry after the main batch drains. Source-availability failures are
    deterministic and should not burn another LLM call.
    """
    if error is None:
        return True
    message = str(error).strip().lower()
    if not message:
        return True
    return not any(token in message for token in _NON_RETRYABLE_EXTRACTION_FAILURES)


def _artifact_supplement_summary(harvest_dir: Path, pmid: str) -> dict[str, Any]:
    """Return supplement/source counters from a run-local artifact log."""
    artifact_path = harvest_dir / f"{pmid}_artifacts.json"
    if not artifact_path.exists():
        return {}
    try:
        artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"artifact_log": str(artifact_path), "artifact_read_error": str(exc)}

    summary = artifact.get("summary") or {}
    main_text = artifact.get("main_text") or {}
    supplements = artifact.get("supplements") or []
    return {
        "artifact_log": str(artifact_path),
        "main_text_source": main_text.get("source"),
        "main_text_chars": summary.get("main_text_chars") or main_text.get("chars"),
        "supplement_refs": main_text.get("supplement_descriptions_count", 0),
        "supplements_downloaded": summary.get("supplement_count", len(supplements)),
        "supplements_converted": summary.get("supplements_converted"),
        "supplements_total_chars": summary.get("supplements_total_chars"),
    }


def _cap_qc_record_from_extraction(
    extraction: Any, harvest_dir: Path
) -> Optional[dict[str, Any]]:
    data = extraction.extracted_data or {}
    metadata = data.get("extraction_metadata") or {}
    scanner_skip = metadata.get("scanner_merge_skipped") or {}
    table_overflow = metadata.get("table_merge_overflow") or {}
    if not scanner_skip and not table_overflow:
        return None

    pmid = str(
        extraction.pmid or (data.get("paper_metadata") or {}).get("pmid") or ""
    ).strip()
    variants = data.get("variants") or []
    record: dict[str, Any] = {
        "pmid": pmid,
        "model_used": metadata.get("model_used") or extraction.model_used,
        "source_type": metadata.get("source_type"),
        "source_file": metadata.get("source_file"),
        "source_size_bytes": metadata.get("source_size_bytes"),
        "final_variant_count": len(variants),
        "scanner_cap_tripped": bool(scanner_skip),
        "scanner_candidate_count": scanner_skip.get("candidate_count"),
        "scanner_safety_cap": scanner_skip.get("safety_cap"),
        "table_merge_cap_tripped": bool(table_overflow),
        "table_candidate_count": table_overflow.get("candidate_count"),
        "table_deduped_count": table_overflow.get("deduped_count"),
        "table_normal_safety_cap": table_overflow.get("normal_safety_cap"),
        "table_overflow_merge_cap": table_overflow.get("overflow_merge_cap"),
        "table_selected_for_merge": table_overflow.get("selected_for_merge"),
        "table_omitted_after_dedupe": table_overflow.get("omitted_after_dedupe"),
        "table_structured_candidate_count": table_overflow.get(
            "structured_candidate_count"
        ),
        "table_regex_only_candidate_count": table_overflow.get(
            "regex_only_candidate_count"
        ),
        "table_merged_added": table_overflow.get("merged_added"),
    }
    if pmid:
        record.update(_artifact_supplement_summary(harvest_dir, pmid))
    return record


def _write_dense_table_overflow_report(
    *,
    extractions: list[Any],
    harvest_dir: Path,
    output_dir: Path,
) -> tuple[Optional[Path], Optional[Path], dict[str, Any]]:
    """Write run-level QC for scanner/table cap trips."""
    records = [
        record
        for extraction in extractions
        if (record := _cap_qc_record_from_extraction(extraction, harvest_dir))
    ]
    summary = {
        "records": len(records),
        "scanner_cap_trips": sum(1 for row in records if row["scanner_cap_tripped"]),
        "table_merge_cap_trips": sum(
            1 for row in records if row["table_merge_cap_tripped"]
        ),
        "table_candidates_omitted_after_dedupe": sum(
            int(row.get("table_omitted_after_dedupe") or 0) for row in records
        ),
        "missing_supplement_ref_pmids": sum(
            1
            for row in records
            if int(row.get("supplement_refs") or 0)
            > int(row.get("supplements_downloaded") or 0)
        ),
    }

    pmid_status_dir = output_dir / "pmid_status"
    pmid_status_dir.mkdir(parents=True, exist_ok=True)
    json_path = pmid_status_dir / "dense_table_overflow.json"
    tsv_path = pmid_status_dir / "dense_table_overflow.tsv"
    json_path.write_text(
        json.dumps({"summary": summary, "records": records}, indent=2),
        encoding="utf-8",
    )

    fieldnames = [
        "pmid",
        "scanner_cap_tripped",
        "scanner_candidate_count",
        "scanner_safety_cap",
        "table_merge_cap_tripped",
        "table_candidate_count",
        "table_deduped_count",
        "table_selected_for_merge",
        "table_omitted_after_dedupe",
        "table_merged_added",
        "supplement_refs",
        "supplements_downloaded",
        "supplements_converted",
        "main_text_source",
        "final_variant_count",
        "source_type",
        "source_file",
    ]
    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in records:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    return json_path, tsv_path, summary


def _source_has_folded_supplements(source_file: Path) -> bool:
    """Return true when a markdown source carries folded supplement content."""
    try:
        from harvesting.supplement_fold import FOLD_BEGIN

        return FOLD_BEGIN in source_file.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False


def _resolve_filter_workers() -> int:
    """Pick the worker count for the Tier-2 LLM filter stage.

    Reads Settings.filter_max_workers directly so the filter stage gets a
    higher default (typically 20) independent of MAX_WORKERS, which users
    legitimately set low (e.g. MAX_WORKERS=1) for extraction reasons. Filter
    calls are small per-request and Anthropic comfortably handles 20+
    concurrent abstract-classification prompts.
    """
    try:
        from config.settings import get_settings

        return get_settings().filter_max_workers
    except Exception:
        return DEFAULT_MAX_WORKERS


COUNT_CLASSIFIER_FIELDS = {"carriers", "affected", "unaffected"}


def _parse_count_classifier_fields(value: Any) -> tuple[str, ...] | None:
    """Normalize COUNT_CLASSIFIER_FIELDS; None means all classifier fields."""
    if value is None:
        return None
    raw = str(value).strip().lower()
    if raw in {"", "all", "*"}:
        return None
    fields: list[str] = []
    for item in raw.split(","):
        field = item.strip().lower()
        if not field:
            continue
        if field not in COUNT_CLASSIFIER_FIELDS:
            raise ValueError(
                f"count classifier field must be one of "
                f"{sorted(COUNT_CLASSIFIER_FIELDS)} or 'all'; got {field!r}"
            )
        if field not in fields:
            fields.append(field)
    return tuple(fields) or None


def _resolve_count_policies() -> Tuple[str, str, tuple[str, ...] | None]:
    """Return (guard_policy, classifier_policy, classifier_fields) from Settings.

    Both default to "off" (strict no-op) so behavior is unchanged unless a
    stage is explicitly opted in via COUNT_GUARD_POLICY / COUNT_CLASSIFIER_POLICY.
    COUNT_CLASSIFIER_FIELDS narrows which fields the classifier may flag/clear;
    None means all fields. Falls back to ("off", "off", None) if Settings can't
    load (e.g. minimal-env tests).
    """
    try:
        from config.settings import get_settings

        settings = get_settings()
        return (
            getattr(settings, "count_guard_policy", "off") or "off",
            getattr(settings, "count_classifier_policy", "off") or "off",
            _parse_count_classifier_fields(
                getattr(settings, "count_classifier_fields", "all")
            ),
        )
    except Exception:
        return ("off", "off", None)


def _nonhuman_veterinary_source_reason(source_text: str) -> Optional[str]:
    """Return a reason when a source is clearly non-human/veterinary.

    These papers can still be useful for variant identities, but their carrier
    counts should not populate human clinical evidence rows. Require a strong
    species signal so background mentions of animal models do not trigger it.
    """
    if not source_text:
        return None
    head = source_text[:20_000]
    if re.search(
        r"(?im)^#{1,3}\s+.*\b(?:cats?\s+and\s+humans?|humans?\s+and\s+cats?)\b",
        head,
    ) or re.search(r"(?im)^#{1,4}\s+HUMAN\b", head):
        return None
    strong_patterns = (
        r"\bMaine\s+Coon\s+cats?\b",
        r"\bfeline\s+hypertrophic\s+cardiomyopathy\b",
        r"\bveterinary\b",
        r"\bAnimals?\s*:\s*[^\n.]{0,240}\b(?:cats?|dogs?|mice|rats?|horses?|"
        r"zebrafish|canine|feline|equine|murine|porcine|bovine)\b",
    )
    for pattern in strong_patterns:
        if re.search(pattern, head, re.IGNORECASE):
            return "source is explicitly non-human/veterinary"

    animal_terms = re.findall(
        r"\b(?:cats?|feline|dogs?|canine|mice|mouse|murine|rats?|zebrafish|"
        r"horses?|equine|porcine|bovine)\b",
        head,
        flags=re.IGNORECASE,
    )
    title_lines = "\n".join(
        line for line in head.splitlines()[:40] if line.lstrip().startswith("#")
    )
    title_animal_signal = re.search(
        r"\b(?:cats?|feline|dogs?|canine|mice|mouse|murine|rats?|zebrafish|"
        r"horses?|equine)\b",
        title_lines,
        flags=re.IGNORECASE,
    )
    if title_animal_signal and len(animal_terms) >= 4:
        return "title and body indicate a non-human source"
    return None


def _apply_nonhuman_clinical_count_guard(
    extracted_data: Optional[Dict[str, Any]], source_text: str
) -> None:
    """Clear human clinical count fields for clearly non-human sources."""
    reason = _nonhuman_veterinary_source_reason(source_text)
    if not reason or not isinstance(extracted_data, dict):
        return
    variants = extracted_data.get("variants")
    if not isinstance(variants, list) or not variants:
        return

    from pipeline.count_outlier_guard import COUNT_FIELDS, _clear_field, _read_field

    flagged = 0
    cleared = 0
    for variant in variants:
        if not isinstance(variant, dict):
            continue
        for field_name, paths in COUNT_FIELDS.items():
            value = _read_field(variant, paths)
            if value is None:
                continue
            flags = variant.setdefault("nonhuman_source_flags", {})
            flags[field_name] = {
                "raw": value,
                "reason": reason,
                "policy": "clear",
            }
            flagged += 1
            _clear_field(variant, paths)
            cleared += 1

    if flagged:
        md = extracted_data.setdefault("extraction_metadata", {})
        md["nonhuman_count_guard"] = {
            "policy": "clear",
            "flagged": flagged,
            "cleared": cleared,
            "reason": reason,
        }


def _apply_count_hygiene(extracted_data: Optional[Dict[str, Any]]) -> None:
    """Run the count-outlier guard then the count classifier in flag/clear mode.

    Applied in-place to the freshly extracted variants right before they are
    persisted to per-PMID JSON. Mirrors the calling convention of
    scripts/apply_count_outlier_guard.py and scripts/apply_count_classifier.py:
    detect first, then apply the configured policy (guard before classifier;
    both are idempotent on a null count).

    Strict no-op when both policies are "off" (the default) — it returns before
    importing the guard/classifier modules so the default extraction path is
    byte-identical to before. In "flag" mode it only annotates variants with
    count_outlier_flags / count_classifier_flags metadata and preserves raw
    counts; "clear" additionally zeros the flagged counts (raw preserved under
    the flag).
    """
    guard_policy, classifier_policy, classifier_fields = _resolve_count_policies()
    if guard_policy == "off" and classifier_policy == "off":
        return
    if not isinstance(extracted_data, dict):
        return
    variants = extracted_data.get("variants")
    if not isinstance(variants, list) or not variants:
        return

    if guard_policy != "off":
        from pipeline.count_outlier_guard import (
            apply_outlier_policy,
            detect_count_outliers,
        )

        guard_annotations = detect_count_outliers(variants)
        guard_result = apply_outlier_policy(
            variants, guard_annotations, policy=guard_policy
        )
        if guard_result.flagged:
            md = extracted_data.setdefault("extraction_metadata", {})
            md["count_outlier_guard"] = {
                "policy": guard_policy,
                "flagged": guard_result.flagged,
                "cleared": guard_result.cleared,
            }

    if classifier_policy != "off":
        from pipeline.count_classifier import (
            detect_misclassified_counts,
            enforce_per_variant_policy,
        )

        classifier_annotations = detect_misclassified_counts(
            variants,
            fields=classifier_fields,
        )
        classifier_result = enforce_per_variant_policy(
            variants, classifier_annotations, policy=classifier_policy
        )
        if classifier_result.flagged:
            md = extracted_data.setdefault("extraction_metadata", {})
            md["count_classifier"] = {
                "policy": classifier_policy,
                "fields": list(classifier_fields) if classifier_fields else "all",
                "flagged": classifier_result.flagged,
                "cleared": classifier_result.cleared,
            }


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
    from utils.gene_metadata import get_gene_aliases

    all_synonyms = list(existing_synonyms) if existing_synonyms else []
    existing_set = {syn.lower() for syn in all_synonyms}
    for alias in get_gene_aliases(gene_symbol, include_query_aliases=False):
        if alias.lower() != gene_symbol.lower() and alias.lower() not in existing_set:
            all_synonyms.append(alias)
            existing_set.add(alias.lower())

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
    max_results: int = 1500,
    synonyms: Optional[List[str]] = None,
    use_pubmind: bool = True,
    use_pubmed: bool = True,
    use_europepmc: bool = False,
    api_key: Optional[str] = None,
    disease: Optional[str] = None,
    disease_terms: Optional[List[str]] = None,
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
        disease: Optional disease term to scope PubMed queries (e.g.
            "atrial fibrillation"). When None (default), behavior is identical
            to the gene-only baseline.
        disease_terms: Optional disease aliases for PubMed query expansion.

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
            disease=disease,
            disease_terms=disease_terms,
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
    tier2_confidence_threshold: float = 0.3,
    filter_max_workers: Optional[int] = None,
    disease: Optional[str] = None,
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
        ClinicalDataTriageFilter(disease=disease)
        if use_clinical_triage
        else InternFilter(
            confidence_threshold=tier2_confidence_threshold, disease=disease
        )
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
    progress_lock = threading.Lock()

    def _write_progress(record: Dict[str, Any]) -> None:
        # Single writer thread/lock — flush per-record for crash safety
        line = json.dumps(record, ensure_ascii=False) + "\n"
        with progress_lock:
            progress_fh.write(line)
            progress_fh.flush()

    def _classify_pmid(pmid: str) -> Tuple[str, Optional[str], Dict[str, Any]]:
        """
        Run tier1/tier2 filters for a single PMID.

        Returns (pmid, fail_reason_or_None_if_pass, progress_record).
        Side-effect-free except for the LLM call inside tier2_filter.filter().
        """
        record_path = abstract_records.get(pmid)
        if not record_path or not Path(record_path).exists():
            reason = "Missing abstract JSON"
            return (
                pmid,
                reason,
                {
                    "pmid": pmid,
                    "final_decision": "FAIL",
                    "reason": reason,
                    "stage": "abstract_load",
                },
            )

        try:
            with open(record_path, "r", encoding="utf-8") as f:
                record = json.load(f)
        except Exception as e:
            reason = f"Abstract JSON read error: {e}"
            return (
                pmid,
                reason,
                {
                    "pmid": pmid,
                    "final_decision": "FAIL",
                    "reason": reason,
                    "stage": "abstract_load",
                },
            )

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

        tier1_result = None
        if enable_tier1:
            tier1_result = keyword_filter.filter(paper)
            if tier1_result.decision is not FilterDecision.PASS:
                reason = tier1_result.reason
                return (
                    pmid,
                    reason,
                    {
                        "pmid": pmid,
                        "final_decision": "FAIL",
                        "reason": reason,
                        "stage": "tier1",
                        "tier1": {
                            "decision": tier1_result.decision.value,
                            "reason": tier1_result.reason,
                            "confidence": tier1_result.confidence,
                        },
                    },
                )

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
                return (
                    pmid,
                    reason,
                    {
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
                    },
                )

        return (
            pmid,
            None,
            {
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
            },
        )

    pending_pmids = [p for p in pmids if p not in processed_pmids]

    # Tier 2 is LLM-bound; parallelize across PMIDs to cut filter wallclock
    # from hours to minutes on large gene queries. Per-record progress writes
    # remain atomic (lock-protected, flushed per record), so kill -9 mid-batch
    # is safe and resume picks up where we left off.
    # FILTER_MAX_WORKERS env var lets you ratchet down concurrency when the
    # backing Azure deployment has tight per-minute quota. When neither the
    # env var nor an explicit kwarg is supplied, we fall through to the
    # filter-specific Settings default (Settings.filter_max_workers, typically
    # 20 for Anthropic). This is independent of MAX_WORKERS, which extraction
    # uses and which users often pin to 1 — at MAX_WORKERS=1 the filter would
    # otherwise serialize to ~21 RPM despite Anthropic having spare capacity.
    env_workers = os.environ.get("FILTER_MAX_WORKERS")
    if env_workers:
        try:
            max_workers = max(1, int(env_workers))
        except ValueError:
            max_workers = (
                filter_max_workers
                if filter_max_workers and filter_max_workers > 0
                else _resolve_filter_workers()
            )
    else:
        max_workers = (
            filter_max_workers
            if filter_max_workers and filter_max_workers > 0
            else _resolve_filter_workers()
        )

    if pending_pmids and enable_tier2:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(_classify_pmid, p): p for p in pending_pmids}
            for fut in as_completed(futures):
                pmid = futures[fut]
                try:
                    out_pmid, fail_reason, record = fut.result()
                except Exception as e:
                    reason = f"Filter exception: {e}"
                    dropped_pmids.append((pmid, reason))
                    _write_progress(
                        {
                            "pmid": pmid,
                            "final_decision": "FAIL",
                            "reason": reason,
                            "stage": "exception",
                        }
                    )
                    continue

                if fail_reason is None:
                    filtered_pmids.append(out_pmid)
                else:
                    dropped_pmids.append((out_pmid, fail_reason))
                _write_progress(record)
    else:
        # Tier 1 only or empty pending list — keep the simple sequential path
        for pmid in pending_pmids:
            out_pmid, fail_reason, record = _classify_pmid(pmid)
            if fail_reason is None:
                filtered_pmids.append(out_pmid)
            else:
                dropped_pmids.append((out_pmid, fail_reason))
            _write_progress(record)

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


def _resolve_corpus_dir() -> Optional[Path]:
    """Resolve the consolidated source corpus dir (GVF_CORPUS_DIR or <repo>/corpus)."""
    env = os.environ.get("GVF_CORPUS_DIR")
    corpus = Path(env) if env else Path(__file__).resolve().parents[1] / "corpus"
    return corpus if corpus.is_dir() else None


def _consolidate_from_corpus(
    pmids: List[str],
    harvest_dir: Path,
    gene_symbol: str,
    corpus_dir: Optional[Path],
) -> set:
    """Reuse already-fetched source from the consolidated corpus cache.

    Looks up ``corpus/<GENE>/<PMID>/`` and, when the cached full text is usable
    (``is_usable_fulltext_source`` — i.e. not a paywall/abstract stub), copies
    it plus its ``_figures``/``_supplements`` into the run so the harvester
    skips re-fetching it. A stub/compromised cached copy is deliberately NOT
    reused, so a fresh run (e.g. after adding a publisher key) will re-attempt
    it. This is the cross-run idempotency cache.

    Returns the set of PMIDs satisfied from the corpus.
    """
    import shutil

    from pipeline.source_quality import is_usable_fulltext_source

    recovered: set = set()
    if not corpus_dir:
        return recovered
    gene_dir = corpus_dir / gene_symbol.upper()
    if not gene_dir.is_dir():
        return recovered

    already_present = {
        f.name.replace("_FULL_CONTEXT.md", "")
        for f in harvest_dir.glob("*_FULL_CONTEXT.md")
    }
    for pmid in (str(p) for p in pmids):
        if pmid in already_present:
            continue
        src_ft = gene_dir / pmid / f"{pmid}_FULL_CONTEXT.md"
        if not src_ft.is_file() or not is_usable_fulltext_source(src_ft):
            continue  # absent or stub/compromised -> let the harvester try
        try:
            shutil.copy2(str(src_ft), str(harvest_dir / f"{pmid}_FULL_CONTEXT.md"))
            for extra in (f"{pmid}_CLEANED.md", f"{pmid}_artifacts.json"):
                s = gene_dir / pmid / extra
                if s.is_file():
                    shutil.copy2(str(s), str(harvest_dir / extra))
            for suffix in ("_figures", "_supplements"):
                s = gene_dir / pmid / f"{pmid}{suffix}"
                d = harvest_dir / f"{pmid}{suffix}"
                if s.is_dir() and not d.exists():
                    shutil.copytree(str(s), str(d))
            recovered.add(pmid)
        except Exception as e:  # noqa: BLE001 - best-effort cache reuse
            logger.warning(f"corpus cache: failed to reuse PMID {pmid}: {e}")

    if recovered:
        logger.info(
            f"✓ corpus cache: reused {len(recovered)} papers from {gene_dir} "
            f"(skipping re-download)"
        )
        print(
            f"\n📦 Corpus cache: reused {len(recovered)} papers from corpus/"
            f"{gene_symbol.upper()} (no re-fetch needed)"
        )
    return recovered


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

    # Copy recovered files and their associated directories (supplements, figures).
    # Uses config.constants.REUSE_FULL_CONTEXT_BYTES (5 KB) — same threshold the
    # orchestrator's is_thin_full_context() check applies, so consolidation
    # and the thin-context gate stay in lockstep. Pre-2026-05-18 this was 500
    # bytes; that mismatch caused abstract stubs to be copied across runs and
    # then immediately re-fetched, defeating the cache.
    from config.constants import REUSE_FULL_CONTEXT_BYTES

    reuse_threshold = int(
        os.environ.get("GVF_REUSE_FULL_CONTEXT_BYTES", REUSE_FULL_CONTEXT_BYTES)
    )
    for pmid, (src_path, size) in prior_files.items():
        if size < reuse_threshold:
            continue  # Skip thin / abstract-only files

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
    from harvesting.content_validation import validate_content_quality

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

    # Corpus cache first (authoritative, quality-gated cross-run cache), then
    # fall back to the legacy prior-run walk for anything not yet in corpus.
    recovered_pmids = _consolidate_from_corpus(
        pmids_to_download, harvest_dir, gene_symbol, _resolve_corpus_dir()
    )
    recovered_pmids |= _consolidate_prior_downloads(
        pmids_to_download, harvest_dir, prior_output_base
    )

    # Only send un-recovered PMIDs to the harvester
    remaining_pmids = [p for p in pmids_to_download if str(p) not in recovered_pmids]
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

    # Post-download content validation across all sources.
    # This catches junk pages, partial content, and malformed outputs before extraction.
    invalid_pmids = {}
    for pmid in list(successfully_downloaded):
        full_context_path = harvest_dir / f"{pmid}_FULL_CONTEXT.md"
        if not full_context_path.exists():
            invalid_pmids[str(pmid)] = "Missing FULL_CONTEXT.md after download"
            successfully_downloaded.discard(str(pmid))
            continue

        try:
            content = full_context_path.read_text(encoding="utf-8", errors="replace")
            is_valid, reason = validate_content_quality(content)
            if not is_valid:
                invalid_pmids[str(pmid)] = reason
                successfully_downloaded.discard(str(pmid))
                logger.warning(
                    "PMID %s failed post-download content validation: %s",
                    pmid,
                    reason,
                )
        except Exception as e:
            invalid_pmids[str(pmid)] = f"Validation read error: {e}"
            successfully_downloaded.discard(str(pmid))
            logger.warning("PMID %s failed content validation read: %s", pmid, e)

    if invalid_pmids:
        paywalled_pmids.update(invalid_pmids)

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
            "post_validation_failed": len(invalid_pmids),
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

    Writes cleaned content to *_CLEANED.md files to avoid mutating source files.
    No API calls — pure regex/string operations.
    """
    from pipeline.preprocessor import PaperPreprocessor

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

    written = 0
    for f in files:
        try:
            pmid = (
                f.name.replace("_FULL_CONTEXT.md", "")
                .replace(f"KCNH2_PMID_", "")
                .replace(f"{gene_symbol}_PMID_", "")
            )
            text = f.read_text(encoding="utf-8", errors="replace")
            total_original += len(text)

            classification = preprocessor.classify(text)
            classifications[classification] = classifications.get(classification, 0) + 1

            cleaned = preprocessor.clean(text)
            cleaned = preprocessor.inject_abstract(cleaned, pmid)
            total_cleaned += len(cleaned)

            cleaned_path = f.with_name(
                f.name.replace("_FULL_CONTEXT.md", "_CLEANED.md")
            )
            cleaned_path.write_text(cleaned, encoding="utf-8")
            written += 1
            processed += 1
        except Exception as e:
            logger.warning(f"Preprocess error for {f.name}: {e}")

    token_savings_pct = (
        round((1 - total_cleaned / max(1, total_original)) * 100, 1)
        if total_original > total_cleaned
        else 0
    )

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
            "cleaned_files_written": written,
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
    from pipeline.data_scout import GeneticDataScout, select_scout_source_path

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
            source_file = select_scout_source_path(md_file)
            content = source_file.read_text(encoding="utf-8")
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
    disease: Optional[str] = None,
    abstract_records: Optional[Dict[str, str]] = None,
    abstract_only_pmids: Optional[List[str]] = None,
    candidate_pmids: Optional[List[str]] = None,
    force_pmids: Optional[List[str]] = None,
    tier_threshold: int = 1,
    max_workers: Optional[int] = None,
    priority_top_n: Optional[int] = None,
    priority_offset: int = 0,
    priority_report_dir: Optional[Path] = None,
    priority_disease_terms: Optional[List[str]] = None,
    triage_mode: Optional[str] = None,
    triage_model: Optional[str] = None,
    triage_report_dir: Optional[Path] = None,
    triage_include_defer: bool = False,
    triage_max_llm_candidates: Optional[int] = None,
    retry_failed_extractions: bool = True,
    extraction_retry_attempts: int = 1,
    extraction_retry_max_workers: Optional[int] = 1,
    extraction_retry_backoff_seconds: float = 30.0,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> StepResult:
    """
    Extract variant data from papers.

    Args:
        harvest_dir: Directory with markdown files
        extraction_dir: Directory to save extractions
        gene_symbol: Gene symbol
        disease: Optional disease term used to interpret affected/unaffected counts
        abstract_records: Dict mapping PMID to abstract JSON path
        abstract_only_pmids: PMIDs that need abstract-only extraction
        candidate_pmids: Optional PMID surface to extract from. When provided,
            full-text files outside this set are ignored. This is important for
            resume/targeted retry runs whose harvest directory may contain older
            cached papers.
        force_pmids: Optional PMID list that must be re-extracted even when an
            existing JSON has a matching source fingerprint. Use for targeted
            replays after extractor/parser fixes.
        tier_threshold: Model cascade threshold
        max_workers: Max parallel workers (None → provider-aware default)
        priority_top_n: Optional cap for prioritized pre-LLM extraction
        priority_offset: Number of top-ranked candidates to skip before selection
        priority_report_dir: Optional directory for prioritization audit artifacts
        priority_disease_terms: Disease aliases used by the prioritizer
        triage_mode: Optional cheap triage mode ("deterministic", "llm", "hybrid")
        triage_model: Optional model for LLM triage (defaults to Tier-2 model)
        triage_report_dir: Optional directory for triage audit artifacts
        triage_include_defer: Include "defer" papers in extraction queue
        triage_max_llm_candidates: Optional cap on LLM triage calls
        retry_failed_extractions: Run a low-concurrency retry pass for failed
            extraction calls after the first batch completes
        extraction_retry_attempts: Number of retry passes to run
        extraction_retry_max_workers: Max workers for retry passes
        extraction_retry_backoff_seconds: Cooldown before each retry pass
        progress_callback: Optional callback(completed, total)

    Returns:
        StepResult with extraction stats
    """
    from pipeline.extraction import ExpertExtractor
    from utils.models import ExtractionResult, Paper
    from utils.pmid_utils import extract_pmid_from_filename

    extraction_dir.mkdir(exist_ok=True)
    candidate_pmid_order: List[str] = []
    candidate_pmid_set: Optional[set[str]] = None
    if candidate_pmids is not None:
        seen_candidate_pmids: set[str] = set()
        for pmid in candidate_pmids:
            clean = str(pmid).strip()
            if not clean or clean in seen_candidate_pmids:
                continue
            seen_candidate_pmids.add(clean)
            candidate_pmid_order.append(clean)
        candidate_pmid_set = set(candidate_pmid_order)
    force_pmid_set = {
        clean for pmid in (force_pmids or []) if (clean := str(pmid).strip())
    }

    # Find markdown files (prefer DATA_ZONES over CLEANED over FULL_CONTEXT)
    data_zones = {
        f.name.replace("_DATA_ZONES.md", ""): f
        for f in harvest_dir.glob("*_DATA_ZONES.md")
        if is_usable_fulltext_source(f)
    }
    cleaned_context = {
        f.name.replace("_CLEANED.md", ""): f
        for f in harvest_dir.glob("*_CLEANED.md")
        if is_usable_fulltext_source(f)
    }
    full_context = {
        f.name.replace("_FULL_CONTEXT.md", ""): f
        for f in harvest_dir.glob("*_FULL_CONTEXT.md")
        if is_usable_fulltext_source(f)
    }

    markdown_files = []
    for pmid in sorted(
        set(data_zones.keys()) | set(cleaned_context.keys()) | set(full_context.keys())
    ):
        if candidate_pmid_set is not None and pmid not in candidate_pmid_set:
            continue
        if pmid in data_zones:
            chosen_file = data_zones[pmid]
        elif pmid in cleaned_context:
            chosen_file = cleaned_context[pmid]
        elif pmid in full_context:
            chosen_file = full_context[pmid]
        else:
            continue

        full_context_file = full_context.get(pmid)
        if (
            full_context_file is not None
            and chosen_file != full_context_file
            and _source_has_folded_supplements(full_context_file)
            and not _source_has_folded_supplements(chosen_file)
        ):
            chosen_file = full_context_file
        markdown_files.append(chosen_file)

    markdown_pmids = {
        pmid
        for md_file in markdown_files
        if (pmid := extract_pmid_from_filename(md_file))
    }

    # Prepare abstract-only papers. If a PMID has any usable markdown source,
    # never submit an abstract task for the same output JSON; the two workers
    # otherwise race and a stale abstract result can overwrite full-text output.
    abstract_papers = []
    if abstract_records:
        if candidate_pmid_set is not None:
            abstract_source_pmids = candidate_pmid_order
        else:
            abstract_source_pmids = abstract_only_pmids or []
        for pmid in abstract_source_pmids:
            if pmid in markdown_pmids:
                continue
            record_path = abstract_records.get(pmid)
            if record_path and Path(record_path).exists():
                abstract_papers.append((pmid, record_path))

    priority_summary: dict[str, Any] = {}
    priority_result = None
    if priority_top_n and priority_top_n > 0:
        from pipeline.extraction_priority import prioritize_extraction_sources

        report_dir = priority_report_dir or (
            extraction_dir.parent / "extraction_priority"
        )
        priority_terms = priority_disease_terms or ([disease] if disease else [])
        priority_result = prioritize_extraction_sources(
            markdown_files=markdown_files,
            abstract_papers=abstract_papers,
            gene_symbol=gene_symbol,
            disease_terms=priority_terms,
            abstract_records=abstract_records,
            top_n=priority_top_n,
            offset=priority_offset,
            report_dir=report_dir,
        )
        markdown_files = priority_result.selected_markdown_files
        abstract_papers = priority_result.selected_abstract_papers
        selected = priority_result.selected_candidates
        priority_summary = {
            "priority_top_n": priority_top_n,
            "priority_offset": priority_offset,
            "priority_candidates": len(priority_result.candidates),
            "priority_selected": len(selected),
            "priority_fulltext_selected": sum(
                1 for candidate in selected if candidate.source_kind == "fulltext"
            ),
            "priority_abstract_selected": sum(
                1 for candidate in selected if candidate.source_kind == "abstract"
            ),
            "priority_report_dir": str(report_dir),
        }
        logger.info(
            "Prioritized extraction selected %s/%s candidates (full text: %s, abstract-only: %s); audit: %s",
            priority_summary["priority_selected"],
            priority_summary["priority_candidates"],
            priority_summary["priority_fulltext_selected"],
            priority_summary["priority_abstract_selected"],
            report_dir,
        )

        if triage_mode and triage_mode.lower() not in {"off", "none"}:
            from pipeline.extraction_triage import (
                apply_triage_filter,
                triage_priority_result,
            )

            triage_dir = triage_report_dir or (
                extraction_dir.parent / "extraction_triage"
            )
            triage_result = triage_priority_result(
                priority_result,
                gene_symbol=gene_symbol,
                disease=disease or "",
                mode=triage_mode,
                model=triage_model,
                max_llm_candidates=triage_max_llm_candidates,
                report_dir=triage_dir,
            )
            priority_result = apply_triage_filter(
                priority_result,
                triage_result,
                include_defer=triage_include_defer,
            )
            markdown_files = priority_result.selected_markdown_files
            abstract_papers = priority_result.selected_abstract_papers
            decisions = triage_result.decisions
            priority_summary.update(
                {
                    "triage_mode": triage_mode,
                    "triage_model": triage_model,
                    "triage_report_dir": str(triage_dir),
                    "triage_total": len(decisions),
                    "triage_extract_now": sum(
                        1
                        for decision in decisions
                        if decision.decision == "extract_now"
                    ),
                    "triage_defer": sum(
                        1 for decision in decisions if decision.decision == "defer"
                    ),
                    "triage_skip": sum(
                        1 for decision in decisions if decision.decision == "skip"
                    ),
                    "triage_selected_for_extraction": len(markdown_files)
                    + len(abstract_papers),
                    "triage_include_defer": triage_include_defer,
                }
            )
            logger.info(
                "Triage selected %s/%s priority candidates for extraction "
                "(extract_now=%s, defer=%s, skip=%s); audit: %s",
                priority_summary["triage_selected_for_extraction"],
                priority_summary["triage_total"],
                priority_summary["triage_extract_now"],
                priority_summary["triage_defer"],
                priority_summary["triage_skip"],
                triage_dir,
            )

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

    def source_fingerprint(source_file: Path) -> dict[str, Any]:
        """Small source identity block saved into extraction metadata."""
        digest = hashlib.sha256()
        with open(source_file, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                digest.update(chunk)
        stat = source_file.stat()
        return {
            "source_file": str(source_file),
            "source_size_bytes": stat.st_size,
            "source_sha256": digest.hexdigest(),
        }

    def metadata_mentions_abstract_only(metadata: dict[str, Any]) -> bool:
        fields: list[str] = []
        for key in ("source_type", "source_kind", "notes"):
            value = metadata.get(key)
            if value:
                fields.append(str(value))
        challenges = metadata.get("challenges")
        if isinstance(challenges, list):
            fields.extend(str(item) for item in challenges)
        elif challenges:
            fields.append(str(challenges))
        joined = " ".join(fields).lower()
        return (
            bool(metadata.get("abstract_only"))
            or "abstract-only" in joined
            or "abstract only" in joined
            or "full text not available" in joined
            or "full text could not be retrieved" in joined
        )

    def load_existing_extraction(
        output_file: Path,
        pmid: str,
        *,
        source_kind: str,
        fingerprint: Optional[dict[str, Any]] = None,
    ) -> Optional[ExtractionResult]:
        """Return a successful result for an already-written extraction JSON."""
        if not output_file.exists():
            return None
        try:
            with open(output_file, "r", encoding="utf-8") as f:
                extracted_data = json.load(f)
            metadata = (
                extracted_data.get("extraction_metadata", {})
                if isinstance(extracted_data, dict)
                else {}
            )
            if source_kind == "fulltext" and metadata_mentions_abstract_only(metadata):
                logger.info(
                    "Ignoring stale abstract-only extraction for PMID %s; full-text source is available",
                    pmid,
                )
                return None

            if fingerprint:
                old_sha = metadata.get("source_sha256")
                old_size = metadata.get("source_size_bytes")
                if old_sha and old_sha != fingerprint.get("source_sha256"):
                    logger.info(
                        "Ignoring stale extraction for PMID %s; source SHA changed",
                        pmid,
                    )
                    return None
                if old_size and old_size != fingerprint.get("source_size_bytes"):
                    logger.info(
                        "Ignoring stale extraction for PMID %s; source size changed",
                        pmid,
                    )
                    return None

            return ExtractionResult(
                pmid=pmid,
                success=True,
                extracted_data=extracted_data,
                model_used=metadata.get("model_used"),
            )
        except Exception as e:
            logger.warning(
                "Could not reuse existing extraction for PMID %s: %s", pmid, e
            )
            return None

    extractions = []
    failures: list[tuple[str, str]] = []
    total = len(markdown_files) + len(abstract_papers)
    completed = 0
    fulltext_by_pmid: dict[str, Path] = {}
    for md_file in markdown_files:
        pmid = extract_pmid_from_filename(md_file)
        if pmid:
            fulltext_by_pmid[pmid] = md_file
    abstract_by_pmid: dict[str, str] = {
        pmid: record_path for pmid, record_path in abstract_papers
    }

    def process_fulltext(md_file: Path) -> Tuple[Any, Optional[Tuple[str, str]]]:
        pmid = extract_pmid_from_filename(md_file)
        if not pmid:
            return (None, (md_file.name, "Could not extract PMID"))

        try:
            if not is_usable_fulltext_source(md_file):
                return (None, (pmid, "Source is abstract-only fallback"))

            output_file = extraction_dir / f"{gene_symbol}_PMID_{pmid}.json"
            fingerprint = source_fingerprint(md_file)
            if pmid not in force_pmid_set:
                existing = load_existing_extraction(
                    output_file,
                    pmid,
                    source_kind="fulltext",
                    fingerprint=fingerprint,
                )
                if existing is not None:
                    return (existing, None)

            content = md_file.read_text(encoding="utf-8")
            paper = Paper(
                pmid=pmid,
                title=f"Paper {pmid}",
                full_text=content,
                gene_symbol=gene_symbol,
                disease=disease,
            )

            result = get_extractor().extract(paper)

            if result.success:
                # Persist model_used into the saved JSON's extraction_metadata
                # so cascading is auditable after the run finishes.
                if result.extracted_data:
                    md = result.extracted_data.setdefault("extraction_metadata", {})
                    md.update(fingerprint)
                    md["source_type"] = "fulltext"
                    md["abstract_only"] = False
                    md["model_used"] = result.model_used
                from pipeline.target_gene_specificity import (
                    apply_target_gene_specificity,
                )

                apply_target_gene_specificity(
                    result.extracted_data,
                    gene_symbol=gene_symbol,
                    source_text=content,
                )
                _apply_nonhuman_clinical_count_guard(
                    result.extracted_data, source_text=content
                )
                # Count hygiene (flag/clear), no-op when both policies are off.
                _apply_count_hygiene(result.extracted_data)
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
            output_file = extraction_dir / f"{gene_symbol}_PMID_{pmid}.json"
            record = None
            record_file = Path(record_path)
            fingerprint = source_fingerprint(record_file)
            if pmid not in force_pmid_set:
                existing = load_existing_extraction(
                    output_file,
                    pmid,
                    source_kind="abstract",
                    fingerprint=fingerprint,
                )
                if existing is not None:
                    return (existing, None)

            with open(record_file, "r", encoding="utf-8") as f:
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
                disease=disease,
            )

            result = get_extractor().extract(paper)

            if result.success:
                if result.extracted_data:
                    md = result.extracted_data.setdefault("extraction_metadata", {})
                    md.update(fingerprint)
                    md["source_type"] = "abstract_only"
                    md["abstract_only"] = True
                    md["model_used"] = result.model_used
                from pipeline.target_gene_specificity import (
                    apply_target_gene_specificity,
                )

                apply_target_gene_specificity(
                    result.extracted_data,
                    gene_symbol=gene_symbol,
                    source_text=abstract_text,
                )
                _apply_nonhuman_clinical_count_guard(
                    result.extracted_data, source_text=abstract_text
                )

                # Count hygiene (flag/clear), no-op when both policies are off.
                _apply_count_hygiene(result.extracted_data)
                with open(output_file, "w") as f:
                    json.dump(result.extracted_data, f, indent=2)
                return (result, None)
            else:
                return (None, (pmid, result.error))

        except Exception as e:
            return (None, (pmid, str(e)))

    # Process in parallel. None / non-positive max_workers → fall back to the
    # provider-aware default (Anthropic 10, Azure 3) so callers that don't
    # care about tuning still benefit from the right default for their
    # configured LLM provider.
    effective_max = (
        max_workers if max_workers and max_workers > 0 else _resolve_default_workers()
    )
    workers = min(effective_max, total) if total > 0 else 1

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

    initial_failures = [(str(pmid), str(error or "")) for pmid, error in failures]
    retry_attempted_pmids: set[str] = set()
    retry_recovered_pmids: set[str] = set()
    retry_skipped_records: list[dict[str, str]] = []
    retry_attempt_records: list[dict[str, Any]] = []
    retry_report_path: Optional[Path] = None

    if retry_failed_extractions and extraction_retry_attempts > 0 and failures:
        retry_workers_config = (
            extraction_retry_max_workers
            if extraction_retry_max_workers and extraction_retry_max_workers > 0
            else 1
        )
        retry_workers = max(1, retry_workers_config)
        retry_backoff = max(0.0, float(extraction_retry_backoff_seconds or 0.0))
        pending_failures = initial_failures
        skipped_pmids: set[str] = set()

        for attempt_number in range(1, extraction_retry_attempts + 1):
            retry_tasks: list[tuple[str, str, Path | str, str]] = []
            for pmid, error in pending_failures:
                if not _is_retryable_extraction_failure(error):
                    if pmid not in skipped_pmids:
                        skipped_pmids.add(pmid)
                        retry_skipped_records.append(
                            {
                                "pmid": pmid,
                                "error": error,
                                "reason": "non_retryable_failure",
                            }
                        )
                    continue
                if pmid in fulltext_by_pmid:
                    retry_tasks.append(
                        ("fulltext", pmid, fulltext_by_pmid[pmid], error)
                    )
                elif pmid in abstract_by_pmid:
                    retry_tasks.append(
                        ("abstract", pmid, abstract_by_pmid[pmid], error)
                    )
                else:
                    if pmid not in skipped_pmids:
                        skipped_pmids.add(pmid)
                        retry_skipped_records.append(
                            {
                                "pmid": pmid,
                                "error": error,
                                "reason": "source_not_found",
                            }
                        )

            if not retry_tasks:
                pending_failures = []
                break

            if retry_backoff > 0:
                logger.info(
                    "Waiting %.1fs before extraction retry pass %s/%s",
                    retry_backoff,
                    attempt_number,
                    extraction_retry_attempts,
                )
                time.sleep(retry_backoff)

            retry_worker_count = min(retry_workers, len(retry_tasks))
            logger.info(
                "Retrying %s failed extractions (attempt %s/%s, workers=%s)",
                len(retry_tasks),
                attempt_number,
                extraction_retry_attempts,
                retry_worker_count,
            )

            attempt_failures: list[tuple[str, str]] = []
            with ThreadPoolExecutor(max_workers=retry_worker_count) as executor:
                retry_futures = {}
                for source_type, pmid, source, initial_error in retry_tasks:
                    retry_attempted_pmids.add(pmid)
                    if source_type == "fulltext":
                        future = executor.submit(process_fulltext, source)
                    else:
                        future = executor.submit(process_abstract, pmid, str(source))
                    retry_futures[future] = (source_type, pmid, initial_error)

                for future in as_completed(retry_futures):
                    source_type, pmid, initial_error = retry_futures[future]
                    try:
                        result, failure = future.result()
                    except Exception as e:
                        result = None
                        failure = (pmid, str(e))

                    if result:
                        extractions.append(result)
                        retry_recovered_pmids.add(pmid)
                        retry_attempt_records.append(
                            {
                                "attempt": attempt_number,
                                "pmid": pmid,
                                "source_type": source_type,
                                "status": "succeeded",
                                "initial_error": initial_error,
                            }
                        )
                    if failure:
                        failed_pmid, error = failure
                        failed_pmid = str(failed_pmid)
                        error = str(error or "")
                        attempt_failures.append((failed_pmid, error))
                        retry_attempt_records.append(
                            {
                                "attempt": attempt_number,
                                "pmid": failed_pmid,
                                "source_type": source_type,
                                "status": "failed",
                                "initial_error": initial_error,
                                "error": error,
                            }
                        )

            pending_failures = attempt_failures
            logger.info(
                "Extraction retry pass %s/%s recovered %s/%s failures",
                attempt_number,
                extraction_retry_attempts,
                len(retry_recovered_pmids),
                len(retry_attempted_pmids),
            )
            if not pending_failures:
                break

        skipped_failures = [
            (record["pmid"], record["error"])
            for record in retry_skipped_records
            if record["pmid"] not in retry_recovered_pmids
        ]
        failures = pending_failures + skipped_failures

    if initial_failures:
        retry_report_path = (
            extraction_dir.parent / "pmid_status" / "extraction_retry_summary.json"
        )
        retry_report_path.parent.mkdir(parents=True, exist_ok=True)
        retry_report_path.write_text(
            json.dumps(
                {
                    "enabled": retry_failed_extractions,
                    "attempts_configured": extraction_retry_attempts,
                    "retry_workers": extraction_retry_max_workers,
                    "retry_backoff_seconds": extraction_retry_backoff_seconds,
                    "initial_failures": [
                        {"pmid": pmid, "error": error}
                        for pmid, error in initial_failures
                    ],
                    "attempts": retry_attempt_records,
                    "skipped": retry_skipped_records,
                    "final_failures": [
                        {"pmid": pmid, "error": error} for pmid, error in failures
                    ],
                },
                indent=2,
            ),
            encoding="utf-8",
        )

    # Save final failures after retry. If all initial failures recovered, write
    # a header-only CSV so stale failure rows from an earlier pass do not linger.
    if failures or initial_failures:
        failures_file = (
            extraction_dir.parent / "pmid_status" / "extraction_failures.csv"
        )
        failures_file.parent.mkdir(parents=True, exist_ok=True)
        with open(failures_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["PMID", "Error"])
            writer.writerows(failures)

    retry_failed_pmids = {
        pmid for pmid, _error in failures if pmid in retry_attempted_pmids
    }
    dense_table_report_path, dense_table_report_tsv, cap_qc_summary = (
        _write_dense_table_overflow_report(
            extractions=extractions,
            harvest_dir=harvest_dir,
            output_dir=extraction_dir.parent,
        )
    )

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
            "initial_failures": len(initial_failures),
            "total_variants": total_variants,
            "candidate_pmid_count": len(candidate_pmid_order)
            if candidate_pmid_set is not None
            else None,
            "extraction_retry_enabled": retry_failed_extractions,
            "extraction_retry_attempts": extraction_retry_attempts,
            "extraction_retry_attempted": len(retry_attempted_pmids),
            "extraction_retry_succeeded": len(retry_recovered_pmids),
            "extraction_retry_failed": len(retry_failed_pmids),
            "extraction_retry_skipped": len(retry_skipped_records),
            "extraction_retry_report": str(retry_report_path)
            if retry_report_path
            else None,
            "dense_table_overflow_records": cap_qc_summary.get("records", 0),
            "scanner_cap_trips": cap_qc_summary.get("scanner_cap_trips", 0),
            "table_merge_cap_trips": cap_qc_summary.get("table_merge_cap_trips", 0),
            "table_candidates_omitted_after_dedupe": cap_qc_summary.get(
                "table_candidates_omitted_after_dedupe", 0
            ),
            "missing_supplement_ref_pmids": cap_qc_summary.get(
                "missing_supplement_ref_pmids", 0
            ),
            "dense_table_overflow_report": str(dense_table_report_path)
            if dense_table_report_path
            else None,
            **priority_summary,
        },
        data={
            "extractions": extractions,
            "failures": failures,
            "initial_failures": initial_failures,
            "retry_report_path": str(retry_report_path) if retry_report_path else None,
            "dense_table_overflow_report": str(dense_table_report_path)
            if dense_table_report_path
            else None,
            "dense_table_overflow_tsv": str(dense_table_report_tsv)
            if dense_table_report_tsv
            else None,
            "priority_report_dir": (
                str(priority_result.report_dir)
                if priority_result and priority_result.report_dir
                else None
            ),
            "triage_report_dir": priority_summary.get("triage_report_dir"),
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

        failed = stats.get("failed", 0)
        if failed:
            # Partial success. Each failed file was rolled back by its SAVEPOINT
            # (no rows written), so this is missing evidence, not corruption --
            # but surface it loudly so a partial migration is never mistaken for
            # a clean one by the caller or the final status.
            logger.warning(
                "⚠️  SQLite migration: %d/%d extraction file(s) failed to "
                "migrate and were rolled back (no rows written for them)",
                failed,
                stats.get("total_files", 0),
            )
            for msg in stats.get("errors", [])[:10]:
                logger.warning("   - %s", msg)

        return StepResult(
            success=True,
            stats={
                "total_files": stats["total_files"],
                "successful": stats["successful"],
                "failed": failed,
                "errors": stats.get("errors", []),
            },
        )

    except Exception as e:
        logger.error(f"Migration failed: {e}")
        return StepResult(success=False, error=str(e))
