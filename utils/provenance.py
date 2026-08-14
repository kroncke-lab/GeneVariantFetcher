"""Best-effort run provenance for reproducibility.

Captures the git SHA (and dirty flag), versioned hashes of executable extractor
tokens and resolved runtime configuration, a legacy raw-byte hash for historical
comparison, and a dependency-lock hash. A run's ``run_manifest.json`` therefore
records what produced it without changing for comments or cosmetic whitespace.
Every collector is best-effort: a failure yields a ``None`` field, never an
exception, so provenance capture can never break a run.
"""

from __future__ import annotations

import hashlib
import io
import json
import subprocess
import sys
import tokenize
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]

# Versioned contract for the extractor fingerprint. Python minor is explicit:
# CPython changed f-string tokenization in 3.12, so cross-minor comparisons must
# not pretend the token stream is the same algorithm.
EXTRACTOR_HASH_ALGORITHM = (
    f"python-token-v1-py{sys.version_info.major}.{sys.version_info.minor}"
)
PROVENANCE_SCHEMA_VERSION = 2

# Historical raw-byte input set retained so old run manifests remain auditable.
LEGACY_PROMPT_EXTRACTOR_FILES = [
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

# Closed set of Python modules whose executable tokens can change scientific
# extraction, recovery, or trust outputs. Runtime configuration is fingerprinted
# separately from resolved values; this list deliberately excludes settings.py
# and this fingerprint implementation itself.
EXTRACTOR_CODE_FILES = [
    path
    for path in LEGACY_PROMPT_EXTRACTOR_FILES
    if path not in {"config/settings.py", "utils/provenance.py"}
] + [
    "config/constants.py",
    "utils/source_layers.py",
]

# Compatibility name for callers that inspect the active extractor file set.
PROMPT_EXTRACTOR_FILES = EXTRACTOR_CODE_FILES

# Files whose content defines the dependency set.
DEPENDENCY_FILES = ["requirements.lock", "pyproject.toml"]

# Explicit resolved-runtime allowlist. Paths, credentials, timestamps, output
# directories, and cohort identifiers are intentionally absent: they identify a
# run but do not define extractor behavior.
EXTRACTOR_MODEL_RESOLVERS = (
    ("tier2_model", "get_tier2_model"),
    ("tier3_models", "get_tier3_models"),
    ("tier3_adjudicator_models", "get_tier3_adjudicator_models"),
    ("table_router_model", "get_table_router_model"),
    ("vision_model", "get_vision_model"),
    ("final_adjudicator_models", "get_final_adjudicator_models"),
    ("final_arbiter_model", "get_final_arbiter_model"),
    ("paper_final_check_model", "get_paper_final_check_model"),
    ("count_recovery_model", "get_count_recovery_model"),
    ("early_debate_models", "get_early_debate_models"),
)
EXTRACTOR_CONFIG_FIELDS = (
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
)


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
    ``EXTRACTOR_CODE_FILES``) changes the hash and cannot silently vanish while
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


def _normalized_python_tokens(source: bytes) -> bytes:
    """Stable executable-token representation for one Python source file.

    Comments, encoding markers, continuation-only newlines, and cosmetic
    whitespace are excluded. Statement boundaries and block structure remain;
    string literals (including prompts and docstrings), names, numbers, and
    operators are retained verbatim. The algorithm identifier includes the
    Python minor version because f-string tokenization changed in Python 3.12.
    """

    kept: list[tuple[int, str]] = []
    ignored = {
        tokenize.ENCODING,
        tokenize.COMMENT,
        tokenize.NL,
        tokenize.ENDMARKER,
    }
    for token in tokenize.tokenize(io.BytesIO(source).readline):
        if token.type in ignored:
            continue
        if token.type in {tokenize.INDENT, tokenize.DEDENT, tokenize.NEWLINE}:
            kept.append((token.type, ""))
        else:
            kept.append((token.type, token.string))
    return json.dumps(kept, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def hash_python_tokens(
    rel_paths: Iterable[str], root: Path = REPO_ROOT
) -> Optional[str]:
    """Hash paths plus normalized Python tokens, folding missing files in."""

    rel_list = list(rel_paths)
    if not rel_list:
        return None
    hasher = hashlib.sha256()
    for rel in rel_list:
        path = root / rel
        hasher.update(rel.encode("utf-8"))
        hasher.update(b"\0")
        if not path.is_file():
            payload = b"<MISSING>"
        elif path.suffix == ".py":
            payload = _normalized_python_tokens(path.read_bytes())
        else:
            payload = path.read_bytes()
        hasher.update(payload)
        hasher.update(b"\0")
    return hasher.hexdigest()


def hash_json(value: Any) -> str:
    """Hash a JSON-compatible value with stable ordering and separators."""

    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def extractor_code_hash() -> Optional[str]:
    return hash_python_tokens(EXTRACTOR_CODE_FILES)


def extractor_config_hash(routing: Optional[dict[str, Any]] = None) -> str:
    """Hash the allowlisted, resolved scientific runtime configuration."""

    return hash_json(resolved_model_routing() if routing is None else routing)


def prompt_extractor_hash(
    routing: Optional[dict[str, Any]] = None,
) -> Optional[str]:
    """Versioned combined digest for extractor code plus resolved config."""

    code_hash = extractor_code_hash()
    if code_hash is None:
        return None
    return hash_json(
        {
            "algorithm": EXTRACTOR_HASH_ALGORITHM,
            "code_sha256": code_hash,
            "config_sha256": extractor_config_hash(routing),
        }
    )


def legacy_prompt_extractor_hash() -> Optional[str]:
    """Pre-v2 raw-byte digest retained solely for historical comparison."""

    return hash_files(LEGACY_PROMPT_EXTRACTOR_FILES)


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

    for key, resolver_name in EXTRACTOR_MODEL_RESOLVERS:
        _try(key, getattr(settings, resolver_name, None))
    for attr in EXTRACTOR_CONFIG_FIELDS:
        routing[attr] = getattr(settings, attr, None)
    return routing


def collect_provenance() -> dict[str, Any]:
    """Gather all provenance fields. Never raises."""
    provenance: dict[str, Any] = {
        "provenance_schema_version": PROVENANCE_SCHEMA_VERSION,
        "extractor_hash_algorithm": EXTRACTOR_HASH_ALGORITHM,
        "extractor_hash_python": f"{sys.version_info.major}.{sys.version_info.minor}",
    }

    def _capture(key: str, fn: Callable[[], Any]) -> None:
        try:
            provenance[key] = fn()
        except Exception:
            provenance[key] = None

    _capture("model_routing", resolved_model_routing)
    routing = provenance["model_routing"]
    if not isinstance(routing, dict):
        routing = {}
        provenance["model_routing"] = routing

    for key, fn in (
        ("git_sha", git_sha),
        ("git_dirty", git_is_dirty),
        ("extractor_code_sha256", extractor_code_hash),
        ("extractor_config_sha256", lambda: extractor_config_hash(routing)),
        ("prompt_extractor_sha256", lambda: prompt_extractor_hash(routing)),
        ("legacy_prompt_extractor_sha256", legacy_prompt_extractor_hash),
        ("dependency_lock_sha256", dependency_lock_hash),
    ):
        _capture(key, fn)
    try:
        missing = missing_files(EXTRACTOR_CODE_FILES)
        provenance["extractor_code_files_missing"] = missing
        provenance["prompt_extractor_files_missing"] = missing
    except Exception:
        provenance["extractor_code_files_missing"] = None
        provenance["prompt_extractor_files_missing"] = None
    return provenance
