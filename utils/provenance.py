"""Best-effort run provenance for reproducibility.

Captures the git SHA (and dirty flag), a content hash of the extraction prompts
plus extractor code, the fully-resolved model routing, and a dependency-lock
hash, so a run's ``run_manifest.json`` records exactly what produced it. Every
collector is best-effort: a failure yields a ``None`` field, never an exception,
so provenance capture can never break a run.
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]

# Files whose content defines "how variants were extracted". A change to any of
# them should change the prompt/extractor hash.
PROMPT_EXTRACTOR_FILES = [
    "config/settings.py",
    "cli/gvf_run.py",
    "pipeline/prompts.py",
    "pipeline/extraction.py",
    "pipeline/extraction_triage.py",
    "pipeline/table_router.py",
    "pipeline/pedigree_extractor.py",
    "pipeline/claim_verifier.py",
    "pipeline/count_recovery.py",
    "pipeline/count_repair.py",
    "pipeline/carrier_guard.py",
    "pipeline/count_outlier_guard.py",
    "pipeline/paper_final_check.py",
    "pipeline/paper_final_check_gate.py",
    "pipeline/paper_final_check_triage.py",
    "pipeline/reference_validation.py",
    "pipeline/somatic_germline_qc.py",
    "pipeline/target_gene_specificity.py",
    "pipeline/trust_gate.py",
    "pipeline/vf_enrichment.py",
    "pipeline/steps.py",
    # behavior that also moves recall/counts, not just the prompt text:
    "pipeline/filters.py",
    "pipeline/count_classifier.py",
    "pipeline/extraction_priority.py",
    "harvesting/figure_text_extractor.py",
    "harvesting/figure_variant_reader.py",
    "harvesting/migrate_to_sqlite.py",
    "utils/llm_utils.py",
    "utils/provenance.py",
    "utils/protein_notation.py",
    "utils/variant_normalizer.py",
    "utils/variant_scanner.py",
]

# Files whose content defines the dependency set.
DEPENDENCY_FILES = ["requirements.lock", "pyproject.toml"]


def _run_git(args: list[str]) -> Optional[str]:
    try:
        out = subprocess.run(
            ["git", *args],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    return out.stdout.strip()


def git_sha() -> Optional[str]:
    return _run_git(["rev-parse", "HEAD"])


def git_is_dirty() -> Optional[bool]:
    status = _run_git(["status", "--porcelain"])
    if status is None:
        return None
    return bool(status.strip())


def hash_files(rel_paths: Iterable[str], root: Path = REPO_ROOT) -> Optional[str]:
    """Stable sha256 over the given files' path+content.

    A listed file that is absent is folded in as a ``<MISSING>`` sentinel rather
    than skipped, so a rename or removal (e.g. a stale entry in
    ``PROMPT_EXTRACTOR_FILES``) changes the hash and cannot silently vanish while
    the digest still looks valid. Returns None only for an empty input list.
    """
    rel_list = list(rel_paths)
    if not rel_list:
        return None
    hasher = hashlib.sha256()
    for rel in rel_list:
        path = root / rel
        hasher.update(rel.encode("utf-8"))
        hasher.update(b"\0")
        hasher.update(path.read_bytes() if path.is_file() else b"<MISSING>")
        hasher.update(b"\0")
    return hasher.hexdigest()


def prompt_extractor_hash() -> Optional[str]:
    return hash_files(PROMPT_EXTRACTOR_FILES)


def dependency_lock_hash() -> Optional[str]:
    return hash_files(DEPENDENCY_FILES)


def missing_files(rel_paths: Iterable[str], root: Path = REPO_ROOT) -> list[str]:
    """Listed paths that are not present. A stale hash-input list is a
    reproducibility hole, so it is recorded in the manifest, not swallowed."""
    return [rel for rel in rel_paths if not (root / rel).is_file()]


def resolved_model_routing() -> dict[str, Any]:
    """Snapshot the resolved model routing (not the raw tier fields).

    Calls the provider-aware resolvers on ``config.settings`` so provider
    switching (anthropic/azure/openai) is captured as the models that will
    actually be used. Each field is guarded so a missing resolver yields None.
    """
    try:
        from config.settings import get_settings

        settings = get_settings()
    except Exception:
        return {}

    routing: dict[str, Any] = {
        "model_provider": getattr(settings, "model_provider", None)
    }

    def _try(key: str, fn: Optional[Callable[[], Any]]) -> None:
        if fn is None:
            routing[key] = None
            return
        try:
            routing[key] = fn()
        except Exception:
            routing[key] = None

    _try("tier2_model", getattr(settings, "get_tier2_model", None))
    _try("tier3_models", getattr(settings, "get_tier3_models", None))
    _try(
        "tier3_adjudicator_models",
        getattr(settings, "get_tier3_adjudicator_models", None),
    )
    _try("table_router_model", getattr(settings, "get_table_router_model", None))
    _try("vision_model", getattr(settings, "get_vision_model", None))
    _try(
        "final_adjudicator_models",
        getattr(settings, "get_final_adjudicator_models", None),
    )
    _try("final_arbiter_model", getattr(settings, "get_final_arbiter_model", None))
    _try(
        "paper_final_check_model",
        getattr(settings, "get_paper_final_check_model", None),
    )
    # Count recovery WRITES scientific counts, so a run's provenance has to
    # record which model filled them and at what effort.
    _try(
        "count_recovery_model",
        getattr(settings, "get_count_recovery_model", None),
    )
    _try("early_debate_models", getattr(settings, "get_early_debate_models", None))
    for attr in (
        "enable_tier1",
        "enable_tier2",
        "enable_tier3",
        "tier1_min_keywords",
        "tier1_use_llm",
        "tier2_temperature",
        "tier2_max_tokens",
        "tier2_confidence_threshold",
        "tier2_reasoning_effort",
        "tier3_temperature",
        "tier3_max_tokens",
        "tier3_threshold",
        "tier3_reasoning_effort",
        "enable_tier3_ensemble_qa",
        "tier3_adjudicator_reasoning_effort",
        "tier3_adjudication_risk_threshold",
        "tier3_evidence_packet_max_chars",
        "tier3_adjudication_max_tokens",
        "tier3_max_verifier_cards",
        "table_router_max_tokens",
        "table_router_reasoning_effort",
        "vision_reasoning_effort",
        "final_adjudicator_reasoning_effort",
        "final_arbiter_reasoning_effort",
        "count_guard_policy",
        "count_classifier_policy",
        "count_classifier_fields",
        "paper_final_check_reasoning_effort",
        "paper_final_check_enabled",
        "paper_final_check_gate_enabled",
        "paper_final_check_gate_min_severity",
        "paper_final_check_gate_require_source_grounded",
        "paper_summary_source_grounded",
        "paper_summary_max_source_chars",
        "count_recovery_enabled",
        "count_recovery_reasoning_effort",
        "count_recovery_fields",
        "count_recovery_max_variants_per_call",
        "count_recovery_max_source_chars",
        "strict_cohort_labels",
        "reference_validation_policy",
        "enable_table_router",
        "scout_enabled",
        "scout_min_relevance",
        "scout_max_zones",
        "scout_use_condensed",
        "extraction_max_chars",
        "scanner_merge_confidence",
        "scanner_max_hints",
        "extract_figures",
        "extract_pedigrees",
    ):
        routing[attr] = getattr(settings, attr, None)
    return routing


def collect_provenance() -> dict[str, Any]:
    """Gather all provenance fields. Never raises."""
    provenance: dict[str, Any] = {}
    for key, fn in (
        ("git_sha", git_sha),
        ("git_dirty", git_is_dirty),
        ("prompt_extractor_sha256", prompt_extractor_hash),
        ("dependency_lock_sha256", dependency_lock_hash),
        ("model_routing", resolved_model_routing),
    ):
        try:
            provenance[key] = fn()
        except Exception:
            provenance[key] = None
    try:
        provenance["prompt_extractor_files_missing"] = missing_files(
            PROMPT_EXTRACTOR_FILES
        )
    except Exception:
        provenance["prompt_extractor_files_missing"] = None
    return provenance
