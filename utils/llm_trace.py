"""Durable, secret-safe traces for LLM calls and curation decisions.

The recorder is deliberately provider-neutral.  It stores the exact textual
request content, safe request parameters, the provider response envelope,
timing/usage/error metadata, and explicit pipeline decision events.  Provider
credentials are redacted and inline image bytes are represented by a digest
rather than duplicated into every trace file.

Private hidden chain-of-thought is **not** available through normal model APIs
and this module never claims to record it.  Three distinct things are recorded,
and they are deliberately kept apart:

``reasoning_capture.response_paths``
    Paths at which the provider returned actual reasoning *content* — a
    reasoning summary, a ``reasoning_content`` string, or a thinking block.
``reasoning_capture.reasoning_token_usage``
    Billing counters such as ``usage.completion_tokens_details.reasoning_tokens``.
    A token count — zero or positive — is **never** evidence that reasoning
    content was exposed, so it can never set
    ``provider_exposed_reasoning_available``.
Explicit rationales
    Model-authored fields such as ``tool_rationale`` or ``count_rationale``
    live in the response envelope and in ``decision_event`` records.  They are
    ordinary model output, not hidden chain-of-thought.

Schema history
--------------
``TRACE_SCHEMA_VERSION = 2`` adds ``reasoning_capture.reasoning_token_usage``,
``context.attempt``-aware records, per-record ``run_id`` enforcement, and the
``verification`` block on the manifest.  Version 1 records remain readable:
every consumer in this repo treats the new keys as optional.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import logging
import os
import re
import threading
import time
import uuid
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Mapping

logger = logging.getLogger(__name__)

TRACE_SCHEMA_VERSION = 2
TRACE_MANIFEST_NAME = "trace_manifest.json"
TRACE_INDEX_NAME = "trace_index.jsonl"

#: Integrity levels, weakest first.  The HTML report and the CLI must name the
#: level they actually have rather than saying "verified" for all three.
INTEGRITY_GENERATED_NOW = "generated_now"
INTEGRITY_WRITE_TIME_VERIFIED = "write_time_verified"
INTEGRITY_LOCKED = "locked"
INTEGRITY_LEVEL_LABELS = {
    INTEGRITY_GENERATED_NOW: (
        "Manifest generated from the files on disk just now. It proves the "
        "records are internally consistent, not that they are unmodified since "
        "they were written."
    ),
    INTEGRITY_WRITE_TIME_VERIFIED: (
        "Every record still matches the SHA-256 digest recorded when the call "
        "was written, so no record changed after capture."
    ),
    INTEGRITY_LOCKED: (
        "Every record matches its write-time digest and the manifest is covered "
        "by a lock file, so nothing changed after capture or after lock."
    ),
}

#: Redaction key names, normalized (lower-cased, ``-`` folded to ``_``).  Header
#: spellings matter here: Azure sends ``api-key``, Anthropic ``x-api-key``, and
#: Elsevier ``X-ELS-Insttoken``.
_SECRET_KEYS = {
    "access_token",
    "api_key",
    "apikey",
    "auth",
    "auth_token",
    "authorization",
    "bearer",
    "client_secret",
    "cookie",
    "cookies",
    "credential",
    "credentials",
    "id_token",
    "insttoken",
    "password",
    "private_key",
    "proxy_authorization",
    "refresh_token",
    "secret",
    "session_token",
    "set_cookie",
    "signature",
    "subscription_key",
    "token",
    "x_api_key",
    "x_els_insttoken",
}
#: Suffix rules so provider-prefixed spellings (``azure_ai_api_key``,
#: ``openai_api_key``, ``els_insttoken``) redact without being enumerated.
_SECRET_KEY_SUFFIXES = (
    "_access_token",
    "_api_key",
    "_apikey",
    "_auth_token",
    "_authorization",
    "_cookie",
    "_credential",
    "_credentials",
    "_insttoken",
    "_key",
    "_password",
    "_secret",
    "_signature",
    "_token",
)
#: Any mapping stored under one of these keys is treated as an HTTP header map
#: and redacted by default — a caller that passes ``extra_headers=`` must not be
#: able to write a live key into a trace just because the header name is new.
_HEADER_MAP_KEYS = {
    "default_headers",
    "extra_headers",
    "header",
    "headers",
    "http_headers",
    "request_headers",
    "response_headers",
}
#: The only header names kept verbatim. Everything else in a header map becomes
#: ``<redacted>``; the *names* are retained so the request stays reviewable.
_SAFE_HEADER_NAMES = {
    "accept",
    "accept_encoding",
    "accept_language",
    "anthropic_beta",
    "anthropic_version",
    "api_version",
    "cache_control",
    "connection",
    "content_encoding",
    "content_length",
    "content_type",
    "date",
    "host",
    "openai_beta",
    "openai_organization",
    "openai_version",
    "request_id",
    "retry_after",
    "user_agent",
    "x_request_id",
    "x_stainless_lang",
    "x_stainless_package_version",
    "x_stainless_runtime",
}
_FALSE_VALUES = {"0", "false", "no", "off", "disabled"}

#: Keys whose value is reasoning *content* when non-empty.
_REASONING_CONTENT_KEYS = {
    "reasoning_content",
    "reasoning_summary",
    "reasoning_text",
    "thinking",
    "thought",
    "thoughts",
}
#: Nested keys under a ``reasoning``-ish container that carry the content.
_REASONING_CONTAINER_KEYS = {"reasoning", "thinking", "thought", "thoughts"}
_REASONING_CONTENT_SUBKEYS = ("summary", "content", "text", "encrypted_content")
#: ``output``/``content`` block types that are reasoning content.
_REASONING_BLOCK_TYPES = {
    "reasoning",
    "reasoning_summary",
    "redacted_thinking",
    "thinking",
}
#: Token counters. These are billing telemetry and never imply exposed content.
_TOKEN_COUNT_KEY_RE = re.compile(r"(?:^|_)(?:tokens?|token_count)$", re.I)
_REASONING_KEY_RE = re.compile(r"(?:reasoning|thought|thinking)", re.I)
_PMID_RE = re.compile(r"\bPMID(?:\s*[:#=_-]\s*|\s+)(\d{5,10})\b", re.I)
_TARGET_GENE_RE = re.compile(
    r"\b(?:target\s+gene|gene(?:\s+symbol)?)\s*[:=]\s*([A-Za-z][A-Za-z0-9-]{1,15})\b",
    re.I,
)
_SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass(frozen=True)
class TraceConfig:
    root: Path
    run_id: str | None
    enabled: bool


_config: TraceConfig | None = None
_config_lock = threading.RLock()
_write_lock = threading.RLock()
_sequence = itertools.count(1)
# ContextVars are per-thread and per-task by construction: a value set inside a
# ThreadPoolExecutor worker is invisible to its siblings.  That isolation is the
# whole point — shared *instance* state on a caller object cross-attributed one
# thread's PMID onto another's trace.
_trace_context: ContextVar[dict[str, Any]] = ContextVar(
    "gvf_llm_trace_context", default={}
)
_last_trace: ContextVar[dict[str, Any] | None] = ContextVar(
    "gvf_llm_last_trace", default=None
)
_attempt_ledger: ContextVar[dict[str, Any] | None] = ContextVar(
    "gvf_llm_attempt_ledger", default=None
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_llm_tracing(
    root: Path | str,
    *,
    run_id: str | None = None,
    enabled: bool = True,
) -> Path:
    """Configure the process-wide trace root and return its resolved path."""
    resolved = Path(root).expanduser().resolve()
    with _config_lock:
        global _config
        _config = TraceConfig(root=resolved, run_id=run_id, enabled=bool(enabled))
    if enabled:
        resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def reset_llm_tracing() -> None:
    """Clear explicit trace configuration (primarily useful for tests)."""
    with _config_lock:
        global _config
        _config = None


def tracing_enabled_by_environment() -> bool:
    raw = os.environ.get("GVF_LLM_TRACE_ENABLED", "true").strip().lower()
    return raw not in _FALSE_VALUES


def safe_run_id(run_id: str) -> str:
    """Filesystem-safe form of a run id, for use as a directory name."""
    return _safe_name(run_id, "run")


@dataclass(frozen=True)
class TraceLocation:
    """A resolved per-run trace directory plus the id its records carry."""

    root: Path
    run_id: str
    base: Path
    #: True when the caller must export root/run_id for nested + subprocess stages.
    newly_selected: bool


def resolve_trace_location(
    run_id: str,
    *,
    default_root: Path | str,
) -> TraceLocation:
    """Pick this run's trace directory, treating ``GVF_LLM_TRACE_DIR`` as a BASE.

    ``GVF_LLM_TRACE_DIR`` names durable *storage* (an encrypted volume, a large
    disk), not one run's directory. Pointing every run straight at it mixed
    records across runs and then made every manifest rebuild raise
    :class:`TraceRunMismatchError`. So the override selects a per-run child named
    after the run id.

    ``GVF_LLM_TRACE_RUN_ID`` is how an already-selected child is handed down. When
    it is set, the configured directory is a child an outer stage already
    resolved, and it is used **verbatim** — that is what stops a nested stage (the
    extraction workflow inside ``gvf-run``, or a subprocess) from appending a
    second level to a path its parent already picked, and what keeps one
    ``gvf-run``'s extraction and post-extraction stages in ONE tree under ONE id.

    Returns the directory, the run id its records will carry, and whether this
    caller newly selected it (and therefore must export both for nested stages).
    """
    env_root = os.environ.get("GVF_LLM_TRACE_DIR", "").strip()
    env_run_id = os.environ.get("GVF_LLM_TRACE_RUN_ID", "").strip()
    if env_root:
        configured = Path(env_root).expanduser()
        if env_run_id:
            return TraceLocation(
                root=configured,
                run_id=env_run_id,
                base=configured.parent,
                newly_selected=False,
            )
        return TraceLocation(
            root=configured / safe_run_id(run_id),
            run_id=run_id,
            base=configured,
            newly_selected=True,
        )
    default = Path(default_root)
    if env_run_id:
        # An outer stage handed down an id without a directory: keep the id so
        # records stay attributable to one run, but use our own default root.
        return TraceLocation(
            root=default, run_id=env_run_id, base=default.parent, newly_selected=False
        )
    return TraceLocation(
        root=default, run_id=run_id, base=default.parent, newly_selected=True
    )


@contextmanager
def exported_trace_identity(root: Path | str, run_id: str) -> Iterator[None]:
    """Publish the resolved trace dir + run id, then restore the prior values.

    Nested in-process stages and subprocesses read these to join THIS run's tree.
    Restoring on exit matters for sequential runs in one process (and for tests):
    a leaked ``GVF_LLM_TRACE_RUN_ID`` would make the next run adopt the previous
    run's identity and write into its directory.
    """
    previous = {
        "GVF_LLM_TRACE_DIR": os.environ.get("GVF_LLM_TRACE_DIR"),
        "GVF_LLM_TRACE_RUN_ID": os.environ.get("GVF_LLM_TRACE_RUN_ID"),
    }
    os.environ["GVF_LLM_TRACE_DIR"] = str(root)
    os.environ["GVF_LLM_TRACE_RUN_ID"] = run_id
    try:
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _active_config() -> TraceConfig | None:
    with _config_lock:
        if _config is not None:
            return _config if _config.enabled else None
    env_root = os.environ.get("GVF_LLM_TRACE_DIR", "").strip()
    if not env_root or not tracing_enabled_by_environment():
        return None
    return TraceConfig(
        root=Path(env_root).expanduser().resolve(),
        run_id=os.environ.get("GVF_LLM_TRACE_RUN_ID", "").strip() or None,
        enabled=True,
    )


@contextmanager
def llm_trace_scope(**context: Any) -> Iterator[None]:
    """Temporarily add paper/stage/component metadata to nested trace calls."""
    merged = dict(_trace_context.get())
    merged.update({key: value for key, value in context.items() if value is not None})
    token = _trace_context.set(merged)
    try:
        yield
    finally:
        _trace_context.reset(token)


def current_trace_context() -> dict[str, Any]:
    return dict(_trace_context.get())


def last_llm_trace() -> dict[str, Any] | None:
    """The index summary of the most recent provider call *on this thread*."""
    return _last_trace.get()


#: Per-attempt outcomes. ``parsed`` means the provider returned usable content —
#: a *candidate* for acceptance, not a claim that its content became the result.
#: Only a selection step (or, absent one, the last candidate) yields ``accepted``.
OUTCOME_PENDING = "pending"
OUTCOME_PARSED = "parsed"
OUTCOME_ACCEPTED = "accepted"
OUTCOME_DISCARDED = "discarded"
OUTCOME_PARSE_FAILED = "parse_failed"
OUTCOME_ERROR = "error"
_CANDIDATE_OUTCOMES = frozenset({OUTCOME_PARSED, OUTCOME_ACCEPTED})
_FAILURE_OUTCOMES = frozenset({OUTCOME_PARSE_FAILED, OUTCOME_ERROR})


@contextmanager
def llm_attempt_ledger() -> Iterator[dict[str, Any]]:
    """Collect the exact trace ids one *logical* LLM call produced.

    A single "extract this paper" call can hit the provider up to ~16 times:
    ``@llm_retry`` attempts, an empty-content retry, and a JSON-repair call
    whose output can become the accepted data.  Every one of those writes its
    own ``llm_call`` record, and without this ledger a curator cannot tell which
    record produced the stored result.

    Open a ledger around the logical call, then attach
    :func:`attempt_link_summary` to the decision event::

        with llm_attempt_ledger() as ledger:
            data = caller.call_llm_json(prompt)
        record_trace_event("paper_extraction_selection", {..., **attempt_link_summary(ledger)})

    Acceptance is **derived, never mutated in place**. A provider call whose
    content parses is marked ``parsed`` — a candidate. *Which* candidate's
    content became the stored result is decided by
    :func:`finalize_attempt_selection` (for stages that pick among several
    models) or, absent a selection, by taking the last candidate. The earlier
    design let every parsing call overwrite ``accepted_response_trace_id``, so a
    multi-model extractor that kept model A's higher-yield result still reported
    model B's call as accepted and never marked A's discarded.

    :func:`note_llm_attempt` / :func:`note_llm_outcome` are no-ops when no
    ledger is active, so instrumented callers stay cheap and uninstrumented
    ones stay correct.
    """
    ledger: dict[str, Any] = {
        "attempts": [],
        "next_attempt": 1,
        "selection": None,
        "repaired": False,
    }
    token = _attempt_ledger.set(ledger)
    try:
        yield ledger
    finally:
        _attempt_ledger.reset(token)


def current_attempt_ledger() -> dict[str, Any] | None:
    return _attempt_ledger.get()


def next_attempt_number() -> int:
    """Monotonic provider-attempt number within the active logical ledger.

    The counter lives in the LEDGER, not in the calling method, because
    ``@llm_retry`` (tenacity) re-enters ``call_llm_json`` wholesale on a
    retriable provider error. A counter local to that method restarted at 1 on
    every re-entry, so one logical call emitted several records all labelled
    ``attempt 1``. Returns 1 when no ledger is active.
    """
    ledger = _attempt_ledger.get()
    if ledger is None:
        return 1
    number = int(ledger.get("next_attempt") or 1)
    ledger["next_attempt"] = number + 1
    return number


def ledger_mark() -> int:
    """Index into the active ledger's attempt list, for scoping a sub-run.

    Pair with :func:`ledger_trace_ids_since` to learn exactly which provider
    calls one candidate produced — the association a multi-model selection needs
    in order to link the calls that contributed to the result it actually kept.
    """
    ledger = _attempt_ledger.get()
    return len(ledger["attempts"]) if ledger else 0


def ledger_trace_ids_since(mark: int) -> list[str]:
    """Trace ids recorded after ``mark``, chronologically."""
    ledger = _attempt_ledger.get()
    if ledger is None:
        return []
    return [
        str(entry["trace_id"])
        for entry in ledger["attempts"][mark:]
        if entry.get("trace_id")
    ]


def finalize_attempt_selection(accepted_trace_ids: Iterable[Any]) -> None:
    """Declare which calls contributed to the result the caller actually kept.

    Every other call that reached the provider successfully becomes
    ``discarded``. Call this once, AFTER the winner has been chosen.
    """
    ledger = _attempt_ledger.get()
    if ledger is None:
        return
    ledger["selection"] = [
        str(trace_id) for trace_id in (accepted_trace_ids or []) if trace_id
    ]


def note_llm_attempt(
    summary: Mapping[str, Any] | None,
    *,
    attempt: int | None = None,
    role: str,
) -> str | None:
    """Register one provider attempt. ``role`` distinguishes the call's purpose.

    Roles in use: ``primary``, ``empty_content_retry``, ``json_repair``,
    ``continuation``, ``adjudication``, ``fallback``, ``table_router``,
    ``figure_ocr``, ``figure_variant_read``, ``pedigree_detect``,
    ``pedigree_extract``, ``count_recovery``, ``claim_debate``.

    ``attempt=None`` (preferred) takes the next monotonic number from the ledger.
    """
    trace_id = (
        str(summary.get("trace_id"))
        if isinstance(summary, Mapping) and summary.get("trace_id")
        else None
    )
    ledger = _attempt_ledger.get()
    if ledger is not None:
        ledger["attempts"].append(
            {
                "attempt": next_attempt_number() if attempt is None else attempt,
                "role": role,
                "trace_id": trace_id,
                "outcome": OUTCOME_PENDING,
            }
        )
    return trace_id


def note_llm_outcome(
    trace_id: str | None,
    outcome: str,
    *,
    repaired: bool | None = None,
) -> None:
    """Record one attempt's raw outcome. Acceptance is derived, not decided here.

    ``outcome`` is one of ``parsed`` (usable content — a candidate),
    ``accepted`` (explicitly the stored result, for single-call stages),
    ``discarded``, ``parse_failed``, ``error``.
    """
    ledger = _attempt_ledger.get()
    if ledger is None:
        return
    for entry in ledger["attempts"]:
        if entry.get("trace_id") == trace_id:
            entry["outcome"] = outcome
    if repaired is not None:
        ledger["repaired"] = bool(repaired)


def _resolve_accepted(active: Mapping[str, Any]) -> list[str]:
    """Chronologically ordered trace ids whose content became the stored result."""
    attempts = list(active.get("attempts") or [])
    candidates = [
        str(entry["trace_id"])
        for entry in attempts
        if entry.get("trace_id") and entry.get("outcome") in _CANDIDATE_OUTCOMES
    ]
    selection = active.get("selection")
    if selection is not None:
        chosen = set(selection)
        return [trace_id for trace_id in candidates if trace_id in chosen]
    explicit = [
        str(entry["trace_id"])
        for entry in attempts
        if entry.get("trace_id") and entry.get("outcome") == OUTCOME_ACCEPTED
    ]
    if explicit:
        return explicit
    # No selection step: the last usable response is the one that was used.
    return candidates[-1:] if candidates else []


def attempt_link_summary(
    ledger: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Flatten a ledger into decision-event fields, safe when nothing was traced.

    ``accepted_response_trace_id`` is the deterministic PRIMARY link: the LAST
    accepted call, i.e. the one whose parsed content became the stored result (a
    JSON repair, when a repair produced the data).
    ``accepted_response_trace_ids`` lists every contributing call in
    chronological order, so a result assembled from a primary call plus a
    continuation names both. ``discarded_trace_ids`` are calls that reached the
    provider but whose content was not used; ``failed_trace_ids`` are parse
    failures and provider errors — kept separate so the two stay distinguishable.
    """
    empty = {
        "accepted_response_trace_id": None,
        "accepted_response_trace_ids": [],
        "discarded_trace_ids": [],
        "failed_trace_ids": [],
        "provider_attempts": 0,
        "attempt_trace_links": [],
        "repaired": False,
        "parse_failures": 0,
    }
    active = ledger if ledger is not None else _attempt_ledger.get()
    if not active:
        return empty

    attempts = list(active.get("attempts") or [])
    accepted = _resolve_accepted(active)
    accepted_set = set(accepted)
    discarded: list[str] = []
    failed: list[str] = []
    links: list[dict[str, Any]] = []
    for entry in attempts:
        trace_id = str(entry["trace_id"]) if entry.get("trace_id") else None
        outcome = str(entry.get("outcome") or OUTCOME_PENDING)
        if trace_id and trace_id in accepted_set:
            outcome = OUTCOME_ACCEPTED
        elif outcome in _FAILURE_OUTCOMES:
            if trace_id:
                failed.append(trace_id)
        else:
            # Reached the provider; its content was not the stored result.
            outcome = OUTCOME_DISCARDED
            if trace_id:
                discarded.append(trace_id)
        links.append(
            {
                "attempt": entry.get("attempt"),
                "role": entry.get("role"),
                "trace_id": trace_id,
                "outcome": outcome,
            }
        )
    return {
        "accepted_response_trace_id": accepted[-1] if accepted else None,
        "accepted_response_trace_ids": accepted,
        "discarded_trace_ids": discarded,
        "failed_trace_ids": failed,
        "provider_attempts": len(attempts),
        "attempt_trace_links": links,
        "repaired": bool(active.get("repaired")),
        "parse_failures": sum(
            1 for entry in attempts if entry.get("outcome") == OUTCOME_PARSE_FAILED
        ),
    }


def normalize_trace_key(key: Any) -> str:
    """Fold a mapping/header key to its redaction-comparison form.

    Hyphens and underscores are equivalent, so ``api-key``, ``API_Key`` and
    ``Api-Key`` all normalize to ``api_key``. Leading ``:`` (HTTP/2
    pseudo-headers) is dropped.
    """
    return str(key).strip().lstrip(":").replace("-", "_").replace(" ", "_").lower()


def _is_secret_key(key: str) -> bool:
    normalized = normalize_trace_key(key)
    return normalized in _SECRET_KEYS or normalized.endswith(_SECRET_KEY_SUFFIXES)


def _is_header_map_key(key: Any) -> bool:
    return normalize_trace_key(key) in _HEADER_MAP_KEYS


def _is_safe_header_name(key: Any) -> bool:
    normalized = normalize_trace_key(key)
    return normalized in _SAFE_HEADER_NAMES and not _is_secret_key(normalized)


def _inline_data_digest(value: str) -> dict[str, Any] | None:
    if not value.startswith("data:") or ";base64," not in value[:120]:
        return None
    header, encoded = value.split(",", 1)
    return {
        "inline_data_omitted": True,
        "media_type": header[5:].split(";", 1)[0],
        "encoding": "base64",
        "encoded_characters": len(encoded),
        "sha256_of_data_url": sha256_bytes(value.encode("utf-8")),
    }


def json_safe(value: Any, *, _depth: int = 0, _in_header_map: bool = False) -> Any:
    """Convert SDK objects to JSON-safe values while redacting credentials.

    ``_in_header_map`` switches the mapping rule from allow-everything-except-
    known-secrets to deny-everything-except-known-safe, because a header map is
    the one place where an unrecognized key name is overwhelmingly likely to be
    a credential (``x-ms-*``, ``x-goog-*``, a publisher token, a signed URL).
    """
    if _depth > 40:
        return "<maximum serialization depth reached>"
    if value is None or isinstance(value, (bool, int, float)):
        return value
    if isinstance(value, str):
        return _inline_data_digest(value) or value
    if isinstance(value, bytes):
        return {
            "binary_omitted": True,
            "bytes": len(value),
            "sha256": sha256_bytes(value),
        }
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        safe: dict[str, Any] = {}
        for key, item in value.items():
            key_text = str(key)
            if _is_secret_key(key_text):
                safe[key_text] = "<redacted>"
            elif _in_header_map and not _is_safe_header_name(key_text):
                safe[key_text] = "<redacted>"
            elif _is_header_map_key(key_text) and isinstance(item, Mapping):
                safe[key_text] = json_safe(item, _depth=_depth + 1, _in_header_map=True)
            else:
                safe[key_text] = json_safe(item, _depth=_depth + 1)
        return safe
    if isinstance(value, (list, tuple, set)):
        return [
            json_safe(item, _depth=_depth + 1, _in_header_map=_in_header_map)
            for item in value
        ]
    model_dump = getattr(value, "model_dump", None)
    if callable(model_dump):
        try:
            return json_safe(
                model_dump(mode="json"),
                _depth=_depth + 1,
                _in_header_map=_in_header_map,
            )
        except TypeError:
            try:
                return json_safe(
                    model_dump(), _depth=_depth + 1, _in_header_map=_in_header_map
                )
            except Exception:
                pass
        except Exception:
            pass
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        try:
            return json_safe(
                to_dict(), _depth=_depth + 1, _in_header_map=_in_header_map
            )
        except Exception:
            pass
    if hasattr(value, "__dict__"):
        try:
            return json_safe(
                vars(value), _depth=_depth + 1, _in_header_map=_in_header_map
            )
        except Exception:
            pass
    # NEVER fall back to repr(): an SDK client, an httpx transport, or an auth
    # object routinely renders its credential inside __repr__. Record the type
    # instead so the trace stays reviewable without leaking.
    return {
        "unserializable_object": True,
        "type": type(value).__name__,
        "module": getattr(type(value), "__module__", None),
    }


def _walk_text(value: Any) -> Iterator[str]:
    if isinstance(value, str):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _walk_text(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _walk_text(item)


def infer_trace_context(request: Mapping[str, Any]) -> dict[str, str]:
    """Best-effort PMID/gene inference when a caller did not provide a scope."""
    combined = "\n".join(_walk_text(request))
    inferred: dict[str, str] = {}
    if match := _PMID_RE.search(combined):
        inferred["pmid"] = match.group(1)
    if match := _TARGET_GENE_RE.search(combined):
        inferred["gene"] = match.group(1).upper()
    return inferred


def _safe_name(value: Any, fallback: str) -> str:
    cleaned = _SAFE_NAME_RE.sub("_", str(value or "")).strip("._")
    return cleaned[:80] or fallback


def _trace_parent(root: Path, context: Mapping[str, Any]) -> Path:
    gene = context.get("gene")
    pmid = context.get("pmid")
    if gene and pmid:
        return root / _safe_name(gene, "UNKNOWN") / _safe_name(pmid, "unknown")
    if pmid:
        return root / "_by_pmid" / _safe_name(pmid, "unknown")
    return root / "_unscoped"


def _is_present(value: Any) -> bool:
    return value not in (None, "", [], {}) and not (
        isinstance(value, str) and not value.strip()
    )


def _reasoning_paths(
    value: Any,
    path: str = "$",
    *,
    content: list[str] | None = None,
    tokens: dict[str, Any] | None = None,
) -> tuple[list[str], dict[str, Any]]:
    """Split reasoning *content* paths from reasoning *token counters*.

    A token counter (``reasoning_tokens``, ``thoughts_token_count``) is billing
    telemetry.  Treating it as evidence of exposed reasoning — which the prior
    implementation did, including for the literal value ``0`` — overstated the
    single claim this recorder must never overstate.
    """
    content = [] if content is None else content
    tokens = {} if tokens is None else tokens

    if isinstance(value, Mapping):
        for key, item in value.items():
            child = f"{path}.{key}"
            normalized = normalize_trace_key(key)
            if _REASONING_KEY_RE.search(normalized) and _TOKEN_COUNT_KEY_RE.search(
                normalized
            ):
                if isinstance(item, (int, float)) and not isinstance(item, bool):
                    tokens[child] = item
                _reasoning_paths(item, child, content=content, tokens=tokens)
                continue
            if normalized in _REASONING_CONTENT_KEYS and _is_present(item):
                content.append(child)
            elif normalized in _REASONING_CONTAINER_KEYS:
                if isinstance(item, str) and _is_present(item):
                    content.append(child)
                elif isinstance(item, Mapping):
                    for subkey in _REASONING_CONTENT_SUBKEYS:
                        if _is_present(item.get(subkey)):
                            content.append(f"{child}.{subkey}")
                elif isinstance(item, list) and any(_is_present(x) for x in item):
                    content.append(child)
            if (
                normalize_trace_key(value.get("type") or "") in _REASONING_BLOCK_TYPES
                and normalized in ("summary", "content", "text")
                and _is_present(item)
            ):
                content.append(child)
            _reasoning_paths(item, child, content=content, tokens=tokens)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _reasoning_paths(item, f"{path}[{index}]", content=content, tokens=tokens)
    return content, tokens


def reasoning_capture(response_payload: Any) -> dict[str, Any]:
    """Describe what reasoning the provider actually exposed, and what it billed."""
    content_paths, token_usage = _reasoning_paths(response_payload)
    return {
        "provider_exposed_reasoning_available": bool(content_paths),
        "response_paths": sorted(set(content_paths)),
        "reasoning_token_usage": dict(sorted(token_usage.items())),
        "note": (
            "response_paths lists provider-returned reasoning CONTENT only. "
            "reasoning_token_usage is billing telemetry and never implies that "
            "any reasoning content was exposed. Hidden chain-of-thought is not "
            "available through model APIs and is not recorded."
        ),
    }


def _output_text(response: Any) -> str | None:
    if isinstance(response, Mapping):
        direct_value = response.get("output_text")
        if isinstance(direct_value, str):
            return direct_value
        for item in response.get("output", []) or []:
            if not isinstance(item, Mapping) or item.get("type") != "message":
                continue
            for content in item.get("content", []) or []:
                if (
                    isinstance(content, Mapping)
                    and content.get("type") == "output_text"
                    and isinstance(content.get("text"), str)
                ):
                    return content["text"]
    direct = getattr(response, "output_text", None)
    if isinstance(direct, str):
        return direct
    try:
        content = response.choices[0].message.content
    except Exception:
        return None
    return content if isinstance(content, str) else None


def _usage(response_payload: Any) -> Any:
    if isinstance(response_payload, Mapping):
        if "usage" in response_payload:
            return response_payload["usage"]
        response = response_payload.get("response")
        if isinstance(response, Mapping) and "usage" in response:
            return response["usage"]
    return None


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _append_index(root: Path, summary: Mapping[str, Any]) -> None:
    root.mkdir(parents=True, exist_ok=True)
    line = json.dumps(summary, sort_keys=True, ensure_ascii=False) + "\n"
    with _write_lock:
        with (root / TRACE_INDEX_NAME).open("a", encoding="utf-8") as handle:
            handle.write(line)
            handle.flush()


def _write_record(
    config: TraceConfig,
    *,
    record_type: str,
    context: Mapping[str, Any],
    payload: dict[str, Any],
) -> dict[str, Any] | None:
    sequence = next(_sequence)
    trace_id = f"{time.time_ns()}-{sequence:06d}-{uuid.uuid4().hex[:8]}"
    component = _safe_name(context.get("component"), record_type)
    stage = _safe_name(context.get("stage"), record_type)
    filename = f"{trace_id}_{stage}_{component}.json"
    path = _trace_parent(config.root, context) / filename
    record = {
        "schema_version": TRACE_SCHEMA_VERSION,
        "record_type": record_type,
        "trace_id": trace_id,
        "run_id": config.run_id,
        "context": json_safe(dict(context)),
        "reasoning_limit": (
            "Private hidden chain-of-thought is not exposed by model APIs. "
            "Only explicit model rationales and provider-returned reasoning "
            "summaries/content can be recorded."
        ),
        **payload,
    }
    try:
        with _write_lock:
            _atomic_write_json(path, record)
            relative = str(path.relative_to(config.root))
            file_sha256 = sha256_path(path)
            summary = {
                "trace_id": trace_id,
                "record_type": record_type,
                "path": relative,
                "sha256": file_sha256,
                "context": json_safe(dict(context)),
                "success": (record.get("response") or {}).get("success"),
                "started_at": (record.get("request") or {}).get("started_at"),
                "completed_at": (record.get("response") or {}).get("completed_at")
                or record.get("recorded_at"),
            }
            _append_index(config.root, summary)
        return summary
    except Exception as exc:  # tracing must never break extraction
        logger.warning("Could not write LLM trace under %s: %s", config.root, exc)
        return None


def capture_llm_call(
    *,
    provider: str,
    requested_model: str,
    resolved_model: str | None,
    request: Mapping[str, Any],
    call: Callable[[], Any],
) -> tuple[Any, dict[str, Any] | None]:
    """Execute one provider call and, when configured, persist its full trace."""
    config = _active_config()
    if config is None:
        _last_trace.set(None)
        return call(), None

    context = current_trace_context()
    for key, value in infer_trace_context(request).items():
        context.setdefault(key, value)
    context.setdefault("provider", provider)
    context.setdefault("model", requested_model)

    safe_request = json_safe(dict(request))
    started_at = utc_now()
    started = time.monotonic()
    response: Any = None
    error: BaseException | None = None
    try:
        response = call()
        return_value = response
    except BaseException as exc:
        error = exc
        return_value = None

    completed_at = utc_now()
    duration = round(time.monotonic() - started, 6)
    safe_response = json_safe(response) if error is None else None
    summary = _write_record(
        config,
        record_type="llm_call",
        context=context,
        payload={
            "request": {
                "started_at": started_at,
                "provider": provider,
                "requested_model": requested_model,
                "resolved_model": resolved_model or requested_model,
                "payload": safe_request,
            },
            "response": {
                "completed_at": completed_at,
                "duration_seconds": duration,
                "success": error is None,
                "output_text": _output_text(response) if error is None else None,
                "envelope": safe_response,
                "usage": _usage(safe_response),
                "error": (
                    None
                    if error is None
                    else {
                        "type": type(error).__name__,
                        "message": str(error),
                    }
                ),
            },
            "reasoning_capture": reasoning_capture(safe_response),
        },
    )
    # Publish the summary on this thread so a wrapper that cannot change its
    # return type (litellm_completion) can still link the accepted response.
    _last_trace.set(summary)
    if error is not None:
        raise error
    return return_value, summary


def record_trace_event(
    event_type: str,
    data: Mapping[str, Any],
    **context: Any,
) -> dict[str, Any] | None:
    """Persist a non-provider event such as route fallback or model selection."""
    config = _active_config()
    if config is None:
        return None
    merged = current_trace_context()
    merged.update({key: value for key, value in context.items() if value is not None})
    merged.setdefault("stage", event_type)
    return _write_record(
        config,
        record_type="decision_event",
        context=merged,
        payload={
            "recorded_at": utc_now(),
            "event": {
                "type": event_type,
                "data": json_safe(dict(data)),
            },
        },
    )


class TraceRunMismatchError(RuntimeError):
    """A manifest rebuild would relabel records written by a different run."""


def read_trace_index(trace_root: Path | str) -> dict[str, dict[str, Any]]:
    """Return ``{relative_path: last write-time index entry}``.

    ``trace_index.jsonl`` is appended the moment a record lands, so its
    ``sha256`` is the *write-time* digest.  Comparing it against a digest taken
    later is the only way a manifest can detect tampering that happened before
    the manifest existed.
    """
    path = Path(trace_root) / TRACE_INDEX_NAME
    entries: dict[str, dict[str, Any]] = {}
    if not path.is_file():
        return entries
    try:
        raw_lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        logger.warning("could not read %s: %s", path, exc)
        return entries
    for line in raw_lines:
        line = line.strip()
        if not line:
            continue
        try:
            entry = json.loads(line)
        except json.JSONDecodeError:
            continue
        relative = str(entry.get("path") or "")
        if relative:
            entries[relative] = entry
    return entries


def _trace_index_digest(trace_root: Path | str) -> str | None:
    path = Path(trace_root) / TRACE_INDEX_NAME
    return sha256_path(path) if path.is_file() else None


def build_trace_manifest(
    trace_root: Path | str,
    *,
    output_path: Path | str | None = None,
    run_id: str | None = None,
    allow_mixed_runs: bool = False,
) -> dict[str, Any]:
    """Hash every trace record into a review-friendly manifest.

    Two guarantees the prior implementation did not provide:

    * Each record's current digest is compared against the **write-time** digest
      in ``trace_index.jsonl``.  A record forged before the manifest existed now
      surfaces as ``write_time_digest_mismatch`` instead of being silently
      blessed by a fresh re-hash.
    * When ``run_id`` is given and records from a *different* run are present,
      the rebuild raises :class:`TraceRunMismatchError` rather than relabelling
      another run's evidence under this run's id.  Pass ``allow_mixed_runs`` to
      record the mixture explicitly instead.

    The manifest is unbounded by design: one entry per record.  Callers that
    need a size bound should shard the trace root, not truncate the manifest.
    """
    root = Path(trace_root).resolve()
    index = read_trace_index(root)
    records: list[dict[str, Any]] = []
    integrity_errors: list[str] = []
    verified_count = 0
    unverifiable_count = 0
    record_run_ids: set[str] = set()

    if root.is_dir():
        for path in sorted(root.rglob("*.json")):
            if path.name == TRACE_MANIFEST_NAME:
                continue
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                payload = {}
            response = payload.get("response") or {}
            event = payload.get("event") or {}
            relative = str(path.relative_to(root))
            current_sha = sha256_path(path)
            write_time_sha = (index.get(relative) or {}).get("sha256")
            if write_time_sha is None:
                unverifiable_count += 1
                integrity_errors.append(
                    f"missing_write_time_digest: {relative} is not in "
                    f"{TRACE_INDEX_NAME}"
                )
            elif write_time_sha != current_sha:
                integrity_errors.append(
                    f"write_time_digest_mismatch: {relative} changed after it was "
                    f"written (index {str(write_time_sha)[:12]}… vs current "
                    f"{current_sha[:12]}…)"
                )
            else:
                verified_count += 1
            record_run_id = payload.get("run_id")
            if record_run_id:
                record_run_ids.add(str(record_run_id))
            records.append(
                {
                    "path": relative,
                    "sha256": current_sha,
                    "write_time_sha256": write_time_sha,
                    "write_time_verified": bool(
                        write_time_sha and write_time_sha == current_sha
                    ),
                    "trace_id": payload.get("trace_id"),
                    "run_id": record_run_id,
                    "record_type": payload.get("record_type"),
                    "context": payload.get("context") or {},
                    "success": response.get("success"),
                    "event_type": event.get("type"),
                    "started_at": (payload.get("request") or {}).get("started_at"),
                    "completed_at": response.get("completed_at")
                    or payload.get("recorded_at"),
                }
            )

    foreign_runs = sorted(record_run_ids - {run_id}) if run_id else []
    if foreign_runs and not allow_mixed_runs:
        raise TraceRunMismatchError(
            f"{root} holds trace records from run(s) {foreign_runs}; refusing to "
            f"rebuild its manifest under run_id={run_id!r}. Use a per-run trace "
            f"directory, or pass allow_mixed_runs=True to record the mixture."
        )
    if foreign_runs:
        integrity_errors.append(
            f"mixed_run_trace_records: {foreign_runs} present alongside {run_id!r}"
        )

    level = (
        INTEGRITY_WRITE_TIME_VERIFIED
        if records and not integrity_errors
        else INTEGRITY_GENERATED_NOW
    )
    manifest = {
        "schema_version": TRACE_SCHEMA_VERSION,
        "run_id": run_id,
        "generated_at": utc_now(),
        "trace_root": str(root),
        "trace_count": len(records),
        "llm_call_count": sum(
            1 for record in records if record["record_type"] == "llm_call"
        ),
        "decision_event_count": sum(
            1 for record in records if record["record_type"] == "decision_event"
        ),
        "record_run_ids": sorted(record_run_ids),
        "index": {
            "name": TRACE_INDEX_NAME,
            "present": bool(index),
            "sha256": _trace_index_digest(root),
            "entry_count": len(index),
        },
        "verification": {
            "level": level,
            "level_description": INTEGRITY_LEVEL_LABELS[level],
            "write_time_verified_records": verified_count,
            "records_without_write_time_digest": unverifiable_count,
            "errors": integrity_errors,
        },
        "records": records,
        "reasoning_limit": (
            "The manifest covers explicit prompts, outputs, rationales, and "
            "provider-exposed reasoning content; it does not claim hidden "
            "chain-of-thought."
        ),
    }
    destination = (
        Path(output_path) if output_path is not None else root / TRACE_MANIFEST_NAME
    )
    _atomic_write_json(destination, manifest)
    return manifest


def validate_trace_manifest(
    trace_root: Path | str,
    manifest: Mapping[str, Any],
    *,
    check_write_time_digests: bool = True,
) -> list[str]:
    """Return integrity errors for a previously generated trace manifest.

    Checks three independent things: the manifest's own digests still match the
    files, the files still match their *write-time* digests from
    ``trace_index.jsonl``, and the index itself has not been swapped out.
    """
    root = Path(trace_root).resolve()
    errors: list[str] = list((manifest.get("verification") or {}).get("errors") or [])
    index = read_trace_index(root) if check_write_time_digests else {}
    expected_paths: set[str] = set()
    for record in manifest.get("records") or []:
        relative = str(record.get("path") or "")
        if not relative:
            errors.append("trace manifest record is missing path")
            continue
        expected_paths.add(relative)
        path = (root / relative).resolve()
        try:
            path.relative_to(root)
        except ValueError:
            errors.append(f"trace path escapes trace root: {relative}")
            continue
        if not path.is_file():
            errors.append(f"trace file is missing: {relative}")
            continue
        current_sha = sha256_path(path)
        if current_sha != record.get("sha256"):
            errors.append(f"trace file changed after manifest: {relative}")
            continue
        if not check_write_time_digests:
            continue
        write_time_sha = (index.get(relative) or {}).get("sha256")
        if write_time_sha is None:
            errors.append(
                f"missing_write_time_digest: {relative} is not in {TRACE_INDEX_NAME}"
            )
        elif write_time_sha != current_sha:
            errors.append(
                f"write_time_digest_mismatch: {relative} changed after it was written"
            )

    recorded_index = manifest.get("index") or {}
    if recorded_index.get("sha256"):
        actual_index_sha = _trace_index_digest(root)
        if actual_index_sha is None:
            errors.append(f"{TRACE_INDEX_NAME} is missing")
        elif actual_index_sha != recorded_index["sha256"]:
            errors.append(f"{TRACE_INDEX_NAME} changed after manifest")

    actual_paths = (
        {
            str(path.relative_to(root))
            for path in root.rglob("*.json")
            if path.name != TRACE_MANIFEST_NAME
        }
        if root.is_dir()
        else set()
    )
    if actual_paths != expected_paths:
        missing = sorted(expected_paths - actual_paths)
        extra = sorted(actual_paths - expected_paths)
        errors.append(f"trace file set mismatch: missing={missing} extra={extra}")
    # Preserve first-seen order while removing repeats from the manifest's own
    # recorded errors so a caller does not report the same defect twice.
    return list(dict.fromkeys(errors))


def integrity_level(
    manifest: Mapping[str, Any],
    errors: list[str] | tuple[str, ...],
    *,
    locked: bool = False,
) -> tuple[str, str]:
    """Return the honest ``(level, description)`` for a manifest + its errors.

    ``locked`` is the caller's assertion that a lock file covers this manifest
    (``benchmarks/codex_paper_eval`` sets it after ``command_lock``).  Without a
    clean write-time comparison the level never rises above ``generated_now``,
    so a manifest generated moments ago from the files it validates can no
    longer print "verified".
    """
    if errors:
        return INTEGRITY_GENERATED_NOW, (
            "Integrity checks reported problems; see the listed errors."
        )
    verification = manifest.get("verification") or {}
    if not manifest.get("records"):
        return INTEGRITY_GENERATED_NOW, INTEGRITY_LEVEL_LABELS[INTEGRITY_GENERATED_NOW]
    if verification.get("level") != INTEGRITY_WRITE_TIME_VERIFIED:
        return INTEGRITY_GENERATED_NOW, INTEGRITY_LEVEL_LABELS[INTEGRITY_GENERATED_NOW]
    level = INTEGRITY_LOCKED if locked else INTEGRITY_WRITE_TIME_VERIFIED
    return level, INTEGRITY_LEVEL_LABELS[level]


def trace_lock_targets(trace_root: Path | str) -> list[Path]:
    """Every file a lock must cover: records, the write-time index, the manifest.

    ``trace_index.jsonl`` is the write-time digest source, so leaving it out of
    the lock set leaves the one artifact that can prove tampering unprotected.
    """
    root = Path(trace_root)
    if not root.is_dir():
        return []
    targets = sorted(path for path in root.rglob("*.json") if path.is_file())
    index_path = root / TRACE_INDEX_NAME
    if index_path.is_file():
        targets.append(index_path)
    return targets
