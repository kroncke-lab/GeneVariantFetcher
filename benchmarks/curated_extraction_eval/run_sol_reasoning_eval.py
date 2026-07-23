#!/usr/bin/env python3
"""Run the frozen curated fixture through the simplified GPT-5.6 Sol protocol.

This runner is intentionally separate from ``run_benchmark.py``.  It performs
no discovery, fetching, source recovery, or fallback extraction.  Every paper
is read from the benchmark's local ``sources/`` snapshot first and, only when
that exact file is absent, from the manifest's local ``corpus_source`` path.

Each reasoning effort gets isolated extraction, telemetry, database, and score
artifacts.  Successful paper results are committed with a hash receipt; resume
accepts a result only when the source, cleaned context, model, literal effort,
and both output hashes still match.  In particular, ``max`` is never aliased to
``xhigh``.  Subsets and pilots are useful for smoke tests but are never scored.
Only a complete run of the canonical fixture uses the production scorer.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import os
import re
import sqlite3
import sys
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any, Iterable, Mapping, Sequence

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from benchmarks.curated_extraction_eval.run_benchmark import (  # noqa: E402
    LOCAL_SOURCE_CORPUS,
    MANIFEST,
    fixture_sha256,
    genes_from_manifest,
    headline,
    load_manifest,
    per_paper_table,
    run_scorer,
)
from harvesting.migrate_to_sqlite import (  # noqa: E402
    create_database_schema,
    migrate_extraction_directory,
    validate_extraction_data,
)


SCHEMA_VERSION = 1
ALLOWED_EFFORTS = ("none", "low", "medium", "high", "xhigh", "max")
DEFAULT_MODEL = "azure_ai/gpt-5.6-sol"
DEFAULT_OUTDIR = HERE / ".sol_runs"
RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


class SolEvalError(RuntimeError):
    """Base class for a fail-closed experiment-runner error."""


class ResumeIdentityError(SolEvalError):
    """An existing run artifact does not match the requested run identity."""


@dataclass(frozen=True)
class PriceConfig:
    input_per_million: float | None = None
    cached_input_per_million: float | None = None
    output_per_million: float | None = None

    def validate(self) -> None:
        for name, value in asdict(self).items():
            if value is not None and value < 0:
                raise SolEvalError(f"{name} must be non-negative")


@dataclass(frozen=True)
class PaperInput:
    key: str
    gene: str
    pmid: str
    title: str
    strategy: str
    source_path: str
    source_bytes: int
    source_chars: int
    source_sha256: str
    cleaned_chars: int
    cleaned_source_sha256: str

    def public_record(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class TaskSpec:
    ordinal: int
    paper: PaperInput
    effort: str


@dataclass
class TaskResult:
    key: str
    effort: str
    status: str
    reused: bool = False
    error: str | None = None


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    payload = json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ResumeIdentityError(f"invalid JSON artifact {path}: {exc}") from exc


def _load_manifest_at(path: Path) -> dict[str, dict[str, str]]:
    """Use the benchmark loader for the canonical fixture; allow test manifests."""
    if path.resolve() == MANIFEST.resolve():
        return load_manifest()
    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = f"{row.get('gene', '').strip().upper()}:{row.get('pmid', '').strip()}"
            if key in rows:
                raise SolEvalError(f"duplicate manifest paper {key}")
            rows[key] = row
    return rows


def parse_efforts(raw: str | Iterable[str]) -> tuple[str, ...]:
    pieces = raw.split(",") if isinstance(raw, str) else list(raw)
    efforts = tuple(
        str(piece).strip().lower() for piece in pieces if str(piece).strip()
    )
    if not efforts:
        raise SolEvalError("at least one reasoning effort is required")
    unknown = [effort for effort in efforts if effort not in ALLOWED_EFFORTS]
    if unknown:
        raise SolEvalError(
            "unsupported reasoning effort(s): "
            + ", ".join(unknown)
            + "; expected literal values "
            + ",".join(ALLOWED_EFFORTS)
        )
    duplicates = sorted({effort for effort in efforts if efforts.count(effort) > 1})
    if duplicates:
        raise SolEvalError("duplicate reasoning effort(s): " + ", ".join(duplicates))
    return efforts


def _parse_csv_set(raw: str | None, *, upper: bool = False) -> set[str] | None:
    if not raw:
        return None
    values = {piece.strip() for piece in raw.split(",") if piece.strip()}
    return {value.upper() for value in values} if upper else values


def select_manifest_rows(
    manifest: Mapping[str, dict[str, str]],
    *,
    genes: set[str] | None = None,
    pmids: set[str] | None = None,
    pilot: int | None = None,
) -> dict[str, dict[str, str]]:
    """Select a deterministic subset; pilots round-robin gene/strategy strata."""
    all_genes = set(genes_from_manifest(dict(manifest)))
    if genes:
        unknown = sorted(genes - all_genes)
        if unknown:
            raise SolEvalError("unknown --genes: " + ", ".join(unknown))

    selected: dict[str, dict[str, str]] = {}
    unmatched_pmids = set(pmids or ())
    for key, row in manifest.items():
        gene = row.get("gene", "").strip().upper()
        pmid = row.get("pmid", "").strip()
        if genes and gene not in genes:
            continue
        if pmids:
            exact = f"{gene}:{pmid}"
            if pmid not in pmids and exact not in pmids:
                continue
            unmatched_pmids.discard(pmid)
            unmatched_pmids.discard(exact)
        selected[key] = row
    if unmatched_pmids:
        raise SolEvalError(
            "unknown --pmids selectors: " + ", ".join(sorted(unmatched_pmids))
        )
    if not selected:
        raise SolEvalError("selection contains no papers")

    if pilot is not None:
        if pilot <= 0:
            raise SolEvalError("--pilot must be positive")
        if pilot < len(selected):
            strata: dict[tuple[str, str], list[tuple[str, dict[str, str]]]] = {}
            for key, row in selected.items():
                stratum = (
                    row.get("gene", "").strip().upper(),
                    row.get("strategy", "").strip().lower(),
                )
                strata.setdefault(stratum, []).append((key, row))
            for values in strata.values():
                values.sort(key=lambda item: item[0])
            chosen: list[tuple[str, dict[str, str]]] = []
            ordered_strata = sorted(strata)
            while len(chosen) < pilot:
                progressed = False
                for stratum in ordered_strata:
                    if strata[stratum] and len(chosen) < pilot:
                        chosen.append(strata[stratum].pop(0))
                        progressed = True
                if not progressed:
                    break
            selected = dict(chosen)
    return dict(sorted(selected.items()))


def resolve_source(
    row: Mapping[str, str], *, sources_dir: Path, repo: Path = REPO
) -> Path:
    """Resolve only the two approved local source locations, in fixed order."""
    gene = row.get("gene", "").strip().upper()
    pmid = row.get("pmid", "").strip()
    frozen = sources_dir / gene / pmid / f"{pmid}_FULL_CONTEXT.md"
    if frozen.is_file():
        return frozen.resolve()

    corpus_value = row.get("corpus_source", "").strip()
    if corpus_value and not re.match(r"^[A-Za-z][A-Za-z0-9+.-]*://", corpus_value):
        corpus_path = Path(corpus_value)
        if not corpus_path.is_absolute():
            corpus_path = repo / corpus_path
        if corpus_path.is_file():
            return corpus_path.resolve()
    raise SolEvalError(
        f"no local full context for {gene}:{pmid}; checked {frozen} and "
        f"manifest corpus_source={corpus_value!r} (network discovery is disabled)"
    )


def _load_sol_module() -> ModuleType:
    module = importlib.import_module("pipeline.sol_extractor")
    for attr in ("extract_paper", "strip_irrelevant_context", "SOL_PROTOCOL_VERSION"):
        if not hasattr(module, attr):
            raise SolEvalError(f"pipeline.sol_extractor is missing required {attr}")
    return module


def cleaner_identity(module: ModuleType) -> dict[str, Any]:
    module_path = Path(str(module.__file__)).resolve()
    return {
        "protocol_version": str(module.SOL_PROTOCOL_VERSION),
        "cleaner_version": str(getattr(module, "SOL_CLEANER_VERSION", "unknown")),
        "module": str(module_path),
        "module_sha256": _sha256_file(module_path),
        "callable": "strip_irrelevant_context",
    }


def build_source_manifest(
    selected: Mapping[str, dict[str, str]],
    *,
    sources_dir: Path,
    sol_module: ModuleType,
    repo: Path = REPO,
) -> list[PaperInput]:
    papers: list[PaperInput] = []
    cleaner = sol_module.strip_irrelevant_context
    for key, row in sorted(selected.items()):
        path = resolve_source(row, sources_dir=sources_dir, repo=repo)
        raw = path.read_bytes()
        text = raw.decode("utf-8", errors="replace")
        cleaned = cleaner(text)
        if not isinstance(cleaned, str):
            raise SolEvalError(f"cleaner returned {type(cleaned).__name__} for {key}")
        papers.append(
            PaperInput(
                key=key,
                gene=row.get("gene", "").strip().upper(),
                pmid=row.get("pmid", "").strip(),
                title=row.get("title", "").strip(),
                strategy=row.get("strategy", "").strip(),
                source_path=str(path),
                source_bytes=len(raw),
                source_chars=len(text),
                source_sha256=_sha256_bytes(raw),
                cleaned_chars=len(cleaned),
                cleaned_source_sha256=_sha256_bytes(cleaned.encode("utf-8")),
            )
        )
    return papers


def interleaved_schedule(
    papers: Sequence[PaperInput], efforts: Sequence[str]
) -> list[TaskSpec]:
    """Rotate effort order per paper to reduce deterministic time-order bias."""
    schedule: list[TaskSpec] = []
    ordinal = 0
    for paper_index, paper in enumerate(sorted(papers, key=lambda item: item.key)):
        for offset in range(len(efforts)):
            effort = efforts[(paper_index + offset) % len(efforts)]
            schedule.append(TaskSpec(ordinal=ordinal, paper=paper, effort=effort))
            ordinal += 1
    return schedule


def _fixture_hash(selected: Mapping[str, dict[str, str]], manifest_path: Path) -> str:
    """Reuse the canonical fixture fingerprint implementation."""
    genes = genes_from_manifest(dict(selected))
    # fixture_sha256 deliberately hashes missing gold as <MISSING>, which keeps
    # synthetic dry-run manifests deterministic without inventing a second hash.
    return fixture_sha256(genes, dict(selected))


def _prepare_run(
    *,
    run_dir: Path,
    run_id: str,
    config: dict[str, Any],
    source_records: list[dict[str, Any]],
    schedule_records: list[dict[str, Any]],
) -> tuple[str, str]:
    config_sha = _sha256_bytes(_json_bytes(config))
    sources_sha = _sha256_bytes(_json_bytes(source_records))
    schedule_sha = _sha256_bytes(_json_bytes(schedule_records))
    run_record = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "config_sha256": config_sha,
        "source_manifest_sha256": sources_sha,
        "schedule_sha256": schedule_sha,
    }
    paths = {
        "config": run_dir / "config_manifest.json",
        "sources": run_dir / "source_manifest.json",
        "schedule": run_dir / "schedule.json",
        "run": run_dir / "run_manifest.json",
    }
    existing = [name for name, path in paths.items() if path.exists()]
    if existing:
        if len(existing) != len(paths):
            raise ResumeIdentityError(
                f"partial run identity in {run_dir}; found {existing}, expected all manifests"
            )
        comparisons = {
            "config_manifest.json": (_read_json(paths["config"]), config),
            "source_manifest.json": (_read_json(paths["sources"]), source_records),
            "schedule.json": (_read_json(paths["schedule"]), schedule_records),
        }
        for label, (observed, expected) in comparisons.items():
            if observed != expected:
                raise ResumeIdentityError(f"resume identity mismatch in {label}")
        observed_run = _read_json(paths["run"])
        for key, expected in run_record.items():
            if observed_run.get(key) != expected:
                raise ResumeIdentityError(
                    f"resume identity mismatch in run_manifest.{key}"
                )
    else:
        if run_dir.exists() and any(run_dir.iterdir()):
            raise ResumeIdentityError(
                f"non-empty run directory has no identity manifests: {run_dir}"
            )
        run_dir.mkdir(parents=True, exist_ok=True)
        _atomic_write_json(paths["config"], config)
        _atomic_write_json(paths["sources"], source_records)
        _atomic_write_json(paths["schedule"], schedule_records)
        _atomic_write_json(paths["run"], {**run_record, "created_at": _utc_now()})
    return config_sha, sources_sha


def _task_paths(run_dir: Path, task: TaskSpec) -> dict[str, Path]:
    stem = f"{task.paper.gene}_PMID_{task.paper.pmid}"
    root = run_dir / "efforts" / task.effort
    return {
        "extraction": root / "extractions" / task.paper.gene / f"{stem}.json",
        "telemetry": root / "telemetry" / task.paper.gene / f"{stem}.json",
        "receipt": root / "receipts" / task.paper.gene / f"{stem}.json",
    }


def _task_identity(task: TaskSpec, *, model: str, config_sha: str) -> dict[str, str]:
    return {
        "config_sha256": config_sha,
        "key": task.paper.key,
        "gene": task.paper.gene,
        "pmid": task.paper.pmid,
        "model": model,
        "reasoning_effort": task.effort,
        "source_sha256": task.paper.source_sha256,
        "cleaned_source_sha256": task.paper.cleaned_source_sha256,
    }


def _validate_extraction_identity(extraction: Any, task: TaskSpec) -> None:
    if not isinstance(extraction, dict):
        raise SolEvalError("extract_paper extraction result is not an object")
    paper_meta = extraction.get("paper_metadata")
    if not isinstance(paper_meta, dict):
        raise SolEvalError("extraction is missing paper_metadata")
    observed_pmid = str(paper_meta.get("pmid", "")).strip()
    if observed_pmid != task.paper.pmid:
        raise SolEvalError(
            f"extraction PMID mismatch: {observed_pmid!r} != {task.paper.pmid!r}"
        )
    for index, variant in enumerate(extraction.get("variants") or []):
        if not isinstance(variant, dict):
            raise SolEvalError(f"variant {index} is not an object")
        observed_gene = str(variant.get("gene_symbol", "")).strip().upper()
        if observed_gene != task.paper.gene:
            raise SolEvalError(
                f"variant {index} gene mismatch: {observed_gene!r} != {task.paper.gene!r}"
            )
    filename = f"{task.paper.gene}_PMID_{task.paper.pmid}.json"
    valid, errors, _warnings = validate_extraction_data(extraction, filename)
    if not valid:
        raise SolEvalError("migration validation failed: " + "; ".join(errors))


def _validate_telemetry_identity(
    telemetry: Any, task: TaskSpec, *, model: str
) -> dict[str, Any]:
    if not isinstance(telemetry, dict):
        raise SolEvalError("extract_paper telemetry result is not an object")
    checks = {
        "model": model,
        "reasoning_effort": task.effort,
        "requested_reasoning_effort": task.effort,
        "effective_reasoning_effort": task.effort,
        "source_sha256": task.paper.source_sha256,
        "cleaned_source_sha256": task.paper.cleaned_source_sha256,
    }
    for key, expected in checks.items():
        observed = telemetry.get(key)
        if observed is not None and str(observed) != expected:
            raise SolEvalError(
                f"telemetry {key} mismatch: {observed!r} != {expected!r}"
            )
    telemetry = dict(telemetry)
    telemetry.update(
        {
            "gene": task.paper.gene,
            "pmid": task.paper.pmid,
            "key": task.paper.key,
            "model": model,
            "reasoning_effort": task.effort,
            "source_sha256": task.paper.source_sha256,
            "cleaned_source_sha256": task.paper.cleaned_source_sha256,
        }
    )
    return telemetry


def _normalize_result(result: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    if isinstance(result, tuple) and len(result) == 2:
        extraction, telemetry = result
    elif hasattr(result, "extraction") and hasattr(result, "telemetry"):
        extraction, telemetry = result.extraction, result.telemetry
    else:
        raise SolEvalError(
            "extract_paper must return (extraction, telemetry) or an object with "
            ".extraction/.telemetry"
        )
    return extraction, telemetry


def _resume_success(
    task: TaskSpec, *, run_dir: Path, model: str, config_sha: str
) -> bool:
    paths = _task_paths(run_dir, task)
    if not paths["receipt"].exists():
        return False
    receipt = _read_json(paths["receipt"])
    identity = _task_identity(task, model=model, config_sha=config_sha)
    if receipt.get("status") != "success" or receipt.get("identity") != identity:
        raise ResumeIdentityError(
            f"resume receipt identity mismatch for {task.paper.key}/{task.effort}"
        )
    for name in ("extraction", "telemetry"):
        path = paths[name]
        if not path.is_file():
            raise ResumeIdentityError(
                f"receipt exists but {name} is missing for {task.paper.key}/{task.effort}"
            )
        observed_hash = _sha256_file(path)
        if receipt.get(f"{name}_sha256") != observed_hash:
            raise ResumeIdentityError(
                f"receipt hash mismatch for {task.paper.key}/{task.effort}/{name}"
            )
    extraction = _read_json(paths["extraction"])
    telemetry = _read_json(paths["telemetry"])
    _validate_extraction_identity(extraction, task)
    _validate_telemetry_identity(telemetry, task, model=model)
    return True


def _execute_task(
    task: TaskSpec,
    *,
    run_dir: Path,
    model: str,
    config_sha: str,
    sol_module: ModuleType,
) -> TaskResult:
    if _resume_success(task, run_dir=run_dir, model=model, config_sha=config_sha):
        return TaskResult(task.paper.key, task.effort, "success", reused=True)

    paths = _task_paths(run_dir, task)
    started = time.monotonic()
    try:
        source_path = Path(task.paper.source_path)
        raw = source_path.read_bytes()
        if _sha256_bytes(raw) != task.paper.source_sha256:
            raise ResumeIdentityError(f"source changed during run: {source_path}")
        source_text = raw.decode("utf-8", errors="replace")
        cleaned = sol_module.strip_irrelevant_context(source_text)
        cleaned_hash = _sha256_bytes(cleaned.encode("utf-8"))
        if cleaned_hash != task.paper.cleaned_source_sha256:
            raise ResumeIdentityError(
                f"cleaned context changed during run for {task.paper.key}"
            )

        result = sol_module.extract_paper(
            gene=task.paper.gene,
            pmid=task.paper.pmid,
            source_text=source_text,
            source_sha256=task.paper.source_sha256,
            reasoning_effort=task.effort,
            model=model,
            title=task.paper.title or None,
        )
        extraction, telemetry = _normalize_result(result)
        _validate_extraction_identity(extraction, task)
        telemetry = _validate_telemetry_identity(telemetry, task, model=model)
        telemetry.setdefault("elapsed_seconds", time.monotonic() - started)
        telemetry.setdefault("extractor_status", telemetry.get("status"))

        # A migration-ready empty extraction is useful for auditing a model
        # failure, but it is not a successful experimental observation. Never
        # mint a success receipt for provider failures, refusals, incomplete
        # responses, or exhausted tool loops: resume must retry those papers
        # instead of freezing the failure as zero recall.
        response_meta = (extraction.get("extraction_metadata") or {}).get(
            "responses"
        ) or {}
        if (
            telemetry.get("extractor_status") == "model_error"
            or response_meta.get("ok") is False
        ):
            telemetry["status"] = "failed"
            telemetry["error_type"] = "ModelToolLoopError"
            telemetry["error"] = (
                telemetry.get("response_error")
                or response_meta.get("error")
                or f"model tool loop ended with {telemetry.get('model_status')!r}"
            )
            _atomic_write_json(paths["extraction"], extraction)
            _atomic_write_json(paths["telemetry"], telemetry)
            return TaskResult(
                task.paper.key,
                task.effort,
                "failed",
                error=str(telemetry["error"]),
            )

        telemetry["status"] = "success"

        _atomic_write_json(paths["extraction"], extraction)
        _atomic_write_json(paths["telemetry"], telemetry)
        receipt = {
            "schema_version": SCHEMA_VERSION,
            "status": "success",
            "identity": _task_identity(task, model=model, config_sha=config_sha),
            "extraction_sha256": _sha256_file(paths["extraction"]),
            "telemetry_sha256": _sha256_file(paths["telemetry"]),
            "completed_at": _utc_now(),
        }
        _atomic_write_json(paths["receipt"], receipt)
        return TaskResult(task.paper.key, task.effort, "success")
    except ResumeIdentityError:
        # A changed source/cleaner or a corrupted receipt invalidates the run,
        # not merely one model call.  Abort instead of laundering identity drift
        # into an ordinary extraction failure.
        raise
    except Exception as exc:  # Keep every failed paper in the denominator.
        failure = {
            "schema_version": SCHEMA_VERSION,
            "status": "failed",
            "key": task.paper.key,
            "gene": task.paper.gene,
            "pmid": task.paper.pmid,
            "model": model,
            "reasoning_effort": task.effort,
            "source_sha256": task.paper.source_sha256,
            "cleaned_source_sha256": task.paper.cleaned_source_sha256,
            "elapsed_seconds": time.monotonic() - started,
            "error_type": type(exc).__name__,
            "error": str(exc)[:4000],
            "traceback": "".join(traceback.format_exception(exc))[-12000:],
            "failed_at": _utc_now(),
        }
        _atomic_write_json(paths["telemetry"], failure)
        return TaskResult(
            task.paper.key,
            task.effort,
            "failed",
            error=f"{type(exc).__name__}: {exc}",
        )


def _token_usage(telemetry: Mapping[str, Any]) -> dict[str, int]:
    blocks: list[Mapping[str, Any]] = []

    def add_block(candidate: Any) -> None:
        if isinstance(candidate, Mapping) and candidate not in blocks:
            blocks.append(candidate)

    # Prefer a declared aggregate usage block.  The protocol adapter may put
    # the Responses result under ``responses`` or another telemetry envelope,
    # so walk bounded JSON-like children as a compatibility fallback.
    add_block(telemetry.get("usage"))
    add_block(telemetry)
    queue: list[Any] = [telemetry]
    seen: set[int] = set()
    while queue and len(seen) < 100:
        current = queue.pop(0)
        if isinstance(current, Mapping):
            if id(current) in seen:
                continue
            seen.add(id(current))
            add_block(current.get("usage"))
            queue.extend(current.values())
        elif isinstance(current, (list, tuple)):
            queue.extend(current)

    def first(*keys: str) -> int:
        for block in blocks:
            for key in keys:
                value = block.get(key)
                if isinstance(value, (int, float)) and value >= 0:
                    return int(value)
        return 0

    input_tokens = first("input_tokens", "prompt_tokens")
    output_tokens = first("output_tokens", "completion_tokens")
    cached_tokens = first("cached_input_tokens", "cached_tokens")
    cache_write_tokens = first(
        "cache_write_input_tokens",
        "cache_creation_input_tokens",
        "cache_write_tokens",
    )
    reasoning_tokens = first("reasoning_tokens")
    reported_total = first("total_tokens")
    details = None
    for block in blocks:
        candidate = block.get("input_tokens_details") or block.get(
            "prompt_tokens_details"
        )
        if isinstance(candidate, Mapping):
            details = candidate
            break
    if details and not cached_tokens:
        cached_tokens = int(details.get("cached_tokens") or 0)
    if details and not cache_write_tokens:
        cache_write_tokens = int(
            details.get("cache_write_input_tokens")
            or details.get("cache_creation_input_tokens")
            or details.get("cache_write_tokens")
            or 0
        )
    output_details = None
    for block in blocks:
        candidate = block.get("output_tokens_details") or block.get(
            "completion_tokens_details"
        )
        if isinstance(candidate, Mapping):
            output_details = candidate
            break
    if output_details and not reasoning_tokens:
        reasoning_tokens = int(output_details.get("reasoning_tokens") or 0)
    return {
        "input_tokens": input_tokens,
        "cached_input_tokens": min(cached_tokens, input_tokens),
        "cache_write_input_tokens": cache_write_tokens,
        "output_tokens": output_tokens,
        "reasoning_tokens": reasoning_tokens,
        # Reasoning tokens are already a subset of output tokens for billing;
        # retain them as a diagnostic dimension without adding them twice.
        "total_tokens": reported_total or input_tokens + output_tokens,
    }


def _estimated_cost(usage: Mapping[str, int], prices: PriceConfig) -> float | None:
    if prices.input_per_million is None or prices.output_per_million is None:
        return None
    cached_price = (
        prices.cached_input_per_million
        if prices.cached_input_per_million is not None
        else prices.input_per_million
    )
    cached = usage["cached_input_tokens"]
    uncached = max(0, usage["input_tokens"] - cached)
    return (
        uncached * prices.input_per_million
        + cached * cached_price
        + usage["output_tokens"] * prices.output_per_million
    ) / 1_000_000


def _effort_telemetry_rows(
    effort: str, papers: Sequence[PaperInput], run_dir: Path
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for paper in papers:
        task = TaskSpec(0, paper, effort)
        paths = _task_paths(run_dir, task)
        path = paths["telemetry"]
        telemetry = _read_json(path) if path.is_file() else {"status": "missing"}
        extraction = (
            _read_json(paths["extraction"]) if paths["extraction"].is_file() else {}
        )
        usage = _token_usage(telemetry)
        response_meta = (
            (extraction.get("extraction_metadata") or {}).get("responses") or {}
            if isinstance(extraction, Mapping)
            else {}
        )
        nested_response = telemetry.get("responses")
        if isinstance(nested_response, Mapping):
            response_meta = {**response_meta, **nested_response}
        operational = telemetry.get("telemetry")
        if not isinstance(operational, Mapping):
            operational = telemetry.get("response_telemetry")
        if not isinstance(operational, Mapping):
            operational = response_meta.get("telemetry")
        if not isinstance(operational, Mapping):
            operational = {}

        response_ok = response_meta.get("ok")
        explicit_loop_failed = telemetry.get("model_tool_loop_failed")
        if explicit_loop_failed is None:
            explicit_loop_failed = telemetry.get("tool_loop_failed")
        extractor_status = str(telemetry.get("extractor_status") or "").lower()
        model_tool_loop_failed = (
            bool(explicit_loop_failed)
            or response_ok is False
            or extractor_status in {"model_error", "tool_loop_error"}
        )
        model_status = (
            telemetry.get("model_status")
            or response_meta.get("status")
            or operational.get("status")
        )
        rows.append(
            {
                "key": paper.key,
                "gene": paper.gene,
                "pmid": paper.pmid,
                "strategy": paper.strategy,
                "status": telemetry.get("status", "unknown"),
                "elapsed_seconds": telemetry.get("elapsed_seconds"),
                "api_calls": telemetry.get("api_calls", operational.get("api_calls")),
                "rounds": telemetry.get("rounds", operational.get("rounds")),
                "retries": telemetry.get("retries", operational.get("retries")),
                "tool_calls": telemetry.get(
                    "tool_calls",
                    telemetry.get(
                        "tool_calls_executed",
                        operational.get("tool_calls_executed"),
                    ),
                ),
                "tool_errors": telemetry.get(
                    "tool_errors", operational.get("tool_errors")
                ),
                "model_tool_loop_failed": model_tool_loop_failed,
                "model_status": model_status,
                "model_error": telemetry.get("model_error")
                or response_meta.get("error"),
                "error_type": telemetry.get("error_type"),
                "error": telemetry.get("error"),
                **usage,
            }
        )
    return rows


def _receipt_hashes_for_effort(
    effort: str,
    papers: Sequence[PaperInput],
    *,
    run_dir: Path,
    model: str,
    config_sha: str,
) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for paper in papers:
        task = TaskSpec(0, paper, effort)
        if not _resume_success(
            task, run_dir=run_dir, model=model, config_sha=config_sha
        ):
            raise SolEvalError(f"missing successful receipt for {paper.key}/{effort}")
        receipt_path = _task_paths(run_dir, task)["receipt"]
        hashes[paper.key] = _sha256_file(receipt_path)
    return hashes


def _validate_database(path: Path, expected_pmids: set[str]) -> dict[str, Any]:
    conn = sqlite3.connect(path)
    try:
        integrity_rows = conn.execute("PRAGMA integrity_check").fetchall()
        integrity = [str(row[0]) for row in integrity_rows]
        foreign_keys = [list(row) for row in conn.execute("PRAGMA foreign_key_check")]
        observed_pmids = {
            str(row[0]) for row in conn.execute("SELECT pmid FROM papers")
        }
    finally:
        conn.close()
    if integrity != ["ok"]:
        raise SolEvalError(f"SQLite integrity_check failed for {path}: {integrity}")
    if foreign_keys:
        raise SolEvalError(
            f"SQLite foreign_key_check failed for {path}: {foreign_keys[:10]}"
        )
    if observed_pmids != expected_pmids:
        raise SolEvalError(
            f"database PMID completeness failed for {path}: expected "
            f"{sorted(expected_pmids)}, observed {sorted(observed_pmids)}"
        )
    return {
        "integrity_check": integrity,
        "foreign_key_violations": foreign_keys,
        "papers_expected": len(expected_pmids),
        "papers_observed": len(observed_pmids),
    }


def build_gene_databases(
    effort: str,
    papers: Sequence[PaperInput],
    *,
    run_dir: Path,
    model: str,
    config_sha: str,
) -> tuple[dict[str, Path], dict[str, Any]]:
    effort_root = run_dir / "efforts" / effort
    manifest_path = effort_root / "database_manifest.json"
    receipt_hashes = _receipt_hashes_for_effort(
        effort, papers, run_dir=run_dir, model=model, config_sha=config_sha
    )
    identity = {
        "schema_version": SCHEMA_VERSION,
        "config_sha256": config_sha,
        "effort": effort,
        "receipt_hashes": receipt_hashes,
    }
    expected_by_gene: dict[str, set[str]] = {}
    for paper in papers:
        expected_by_gene.setdefault(paper.gene, set()).add(paper.pmid)

    if manifest_path.exists():
        observed = _read_json(manifest_path)
        if observed.get("identity") != identity:
            raise ResumeIdentityError(
                f"database manifest identity mismatch for {effort}"
            )
        dbs: dict[str, Path] = {}
        for gene, expected_pmids in sorted(expected_by_gene.items()):
            path = effort_root / "databases" / f"{gene}.db"
            record = (observed.get("databases") or {}).get(gene, {})
            if not path.is_file() or record.get("sha256") != _sha256_file(path):
                raise ResumeIdentityError(f"database hash mismatch for {effort}/{gene}")
            _validate_database(path, expected_pmids)
            dbs[gene] = path
        return dbs, observed

    db_records: dict[str, Any] = {}
    dbs = {}
    for gene, expected_pmids in sorted(expected_by_gene.items()):
        extraction_dir = effort_root / "extractions" / gene
        expected_files = {f"{gene}_PMID_{pmid}.json" for pmid in expected_pmids}
        observed_files = {path.name for path in extraction_dir.glob("*.json")}
        if observed_files != expected_files:
            raise SolEvalError(
                f"extraction completeness failed for {effort}/{gene}: expected "
                f"{sorted(expected_files)}, observed {sorted(observed_files)}"
            )

        database_dir = effort_root / "databases"
        database_dir.mkdir(parents=True, exist_ok=True)
        final_path = database_dir / f"{gene}.db"
        temporary = database_dir / f".{gene}.{uuid.uuid4().hex}.tmp.db"
        conn = create_database_schema(str(temporary))
        try:
            migration = migrate_extraction_directory(conn, extraction_dir)
        finally:
            conn.close()
        if (
            migration.get("total_files") != len(expected_files)
            or migration.get("successful") != len(expected_files)
            or migration.get("failed") != 0
        ):
            raise SolEvalError(
                f"migration completeness failed for {effort}/{gene}: {migration}"
            )
        integrity = _validate_database(temporary, expected_pmids)
        os.replace(temporary, final_path)
        dbs[gene] = final_path
        db_records[gene] = {
            "path": str(final_path.relative_to(run_dir)),
            "sha256": _sha256_file(final_path),
            "migration": migration,
            **integrity,
        }
    manifest = {"identity": identity, "databases": db_records, "created_at": _utc_now()}
    _atomic_write_json(manifest_path, manifest)
    return dbs, manifest


def score_complete_effort(
    effort: str,
    papers: Sequence[PaperInput],
    manifest_rows: Mapping[str, dict[str, str]],
    *,
    run_dir: Path,
    model: str,
    config_sha: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    dbs, database_manifest = build_gene_databases(
        effort, papers, run_dir=run_dir, model=model, config_sha=config_sha
    )
    effort_root = run_dir / "efforts" / effort
    score_receipt_path = effort_root / "score_receipt.json"
    score_identity = {
        "schema_version": SCHEMA_VERSION,
        "config_sha256": config_sha,
        "effort": effort,
        "database_hashes": {
            gene: record["sha256"]
            for gene, record in sorted(database_manifest["databases"].items())
        },
    }
    if score_receipt_path.exists():
        receipt = _read_json(score_receipt_path)
        if receipt.get("identity") != score_identity:
            raise ResumeIdentityError(f"score receipt identity mismatch for {effort}")
        score_dir = run_dir / receipt["score_dir"]
        summary_path = score_dir / "summary.json"
        if not summary_path.is_file() or receipt.get("summary_sha256") != _sha256_file(
            summary_path
        ):
            raise ResumeIdentityError(f"score summary hash mismatch for {effort}")
        summary = _read_json(summary_path)
        papers_table = per_paper_table(score_dir, dict(manifest_rows))
        return summary, papers_table

    score_dir = effort_root / "scores" / f"score-{uuid.uuid4().hex}"
    summary = run_scorer(dbs, score_dir)
    papers_table = per_paper_table(score_dir, dict(manifest_rows))
    _atomic_write_json(
        score_receipt_path,
        {
            "identity": score_identity,
            "score_dir": str(score_dir.relative_to(run_dir)),
            "summary_sha256": _sha256_file(score_dir / "summary.json"),
            "completed_at": _utc_now(),
        },
    )
    return summary, papers_table


def merge_paper_metrics(
    telemetry_rows: list[dict[str, Any]], per_paper: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    metrics = {(str(row["gene"]), str(row["pmid"])): row for row in per_paper}
    return [
        {**row, "score": metrics.get((row["gene"], row["pmid"]))}
        for row in telemetry_rows
    ]


def _summarize_effort(
    effort: str,
    telemetry_rows: list[dict[str, Any]],
    *,
    prices: PriceConfig,
    score_summary: dict[str, Any] | None,
    scoring_reason: str | None,
) -> dict[str, Any]:
    totals = {
        "input_tokens": sum(row["input_tokens"] for row in telemetry_rows),
        "cached_input_tokens": sum(
            row["cached_input_tokens"] for row in telemetry_rows
        ),
        "cache_write_input_tokens": sum(
            row["cache_write_input_tokens"] for row in telemetry_rows
        ),
        "output_tokens": sum(row["output_tokens"] for row in telemetry_rows),
        "reasoning_tokens": sum(row["reasoning_tokens"] for row in telemetry_rows),
        "total_tokens": sum(row["total_tokens"] for row in telemetry_rows),
    }
    failed = [
        {
            "key": row["key"],
            "error_type": row.get("error_type"),
            "error": row.get("error"),
        }
        for row in telemetry_rows
        if row.get("status") != "success"
    ]
    elapsed_values = [
        float(row["elapsed_seconds"])
        for row in telemetry_rows
        if isinstance(row.get("elapsed_seconds"), (int, float))
    ]
    return {
        "effort": effort,
        "papers_total": len(telemetry_rows),
        "papers_successful": len(telemetry_rows) - len(failed),
        "papers_failed": len(failed),
        "success_rate": (
            (len(telemetry_rows) - len(failed)) / len(telemetry_rows)
            if telemetry_rows
            else None
        ),
        "failures": failed,
        "model_tool_loop_failures": sum(
            bool(row.get("model_tool_loop_failed")) for row in telemetry_rows
        ),
        "usage": totals,
        "elapsed_paper_seconds_sum": sum(elapsed_values),
        "estimated_cost": _estimated_cost(totals, prices),
        "scored": score_summary is not None,
        "scoring_skipped_reason": scoring_reason,
        "headline": headline(score_summary) if score_summary is not None else None,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--sources-dir", type=Path, default=LOCAL_SOURCE_CORPUS)
    parser.add_argument("--genes", help="Comma-separated gene subset (never scored).")
    parser.add_argument(
        "--pmids",
        help="Comma-separated PMID or GENE:PMID selectors (never scored as a subset).",
    )
    parser.add_argument(
        "--pilot",
        type=int,
        help="Deterministic round-robin paper count for a smoke run (never scored).",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument(
        "--efforts",
        default=",".join(ALLOWED_EFFORTS),
        help="Comma-separated literal efforts: " + ",".join(ALLOWED_EFFORTS),
    )
    parser.add_argument(
        "--run-id",
        default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve, clean, hash, and print the plan without writing or calling a model.",
    )
    parser.add_argument("--input-price-per-million", type=float)
    parser.add_argument("--cached-input-price-per-million", type=float)
    parser.add_argument("--output-price-per-million", type=float)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not RUN_ID_RE.fullmatch(args.run_id):
        raise SolEvalError(
            "--run-id must start with an alphanumeric character and contain only "
            "letters, digits, dot, underscore, or hyphen"
        )
    if args.workers <= 0:
        raise SolEvalError("--workers must be positive")
    efforts = parse_efforts(args.efforts)
    prices = PriceConfig(
        input_per_million=args.input_price_per_million,
        cached_input_per_million=args.cached_input_price_per_million,
        output_per_million=args.output_price_per_million,
    )
    prices.validate()

    manifest_path = args.manifest.resolve()
    sources_dir = args.sources_dir.resolve()
    manifest = _load_manifest_at(manifest_path)
    selected = select_manifest_rows(
        manifest,
        genes=_parse_csv_set(args.genes, upper=True),
        pmids=_parse_csv_set(args.pmids),
        pilot=args.pilot,
    )
    sol_module = _load_sol_module()
    papers = build_source_manifest(
        selected, sources_dir=sources_dir, sol_module=sol_module
    )
    schedule = interleaved_schedule(papers, efforts)
    selected_is_canonical_full_fixture = (
        manifest_path == MANIFEST.resolve()
        and set(selected) == set(manifest)
        and len(selected) == len(manifest)
    )
    config = {
        "schema_version": SCHEMA_VERSION,
        "protocol": "simplified_sol_reasoning_eval",
        "manifest_path": str(manifest_path),
        "sources_dir": str(sources_dir),
        "fixture_sha256": _fixture_hash(selected, manifest_path),
        "canonical_full_fixture": selected_is_canonical_full_fixture,
        "paper_keys": [paper.key for paper in papers],
        "model": args.model,
        "reasoning_efforts": list(efforts),
        "reasoning_effort_policy": "literal_no_aliases",
        "workers": args.workers,
        "cleaner": cleaner_identity(sol_module),
        "runtime_code": {
            "driver_sha256": _sha256_file(Path(__file__).resolve()),
            "responses_api_sha256": _sha256_file(
                Path(importlib.import_module("utils.responses_api").__file__).resolve()
            ),
        },
        "prices_per_million_tokens": asdict(prices),
    }
    source_records = [paper.public_record() for paper in papers]
    schedule_records = [
        {
            "ordinal": task.ordinal,
            "key": task.paper.key,
            "effort": task.effort,
        }
        for task in schedule
    ]
    if args.dry_run:
        print(
            json.dumps(
                {
                    "dry_run": True,
                    "run_id": args.run_id,
                    "papers": len(papers),
                    "tasks": len(schedule),
                    "efforts": list(efforts),
                    "model": args.model,
                    "workers": args.workers,
                    "canonical_full_fixture": selected_is_canonical_full_fixture,
                    "fixture_sha256": config["fixture_sha256"],
                    "cleaner": config["cleaner"],
                    "source_manifest_sha256": _sha256_bytes(
                        _json_bytes(source_records)
                    ),
                    "schedule_sha256": _sha256_bytes(_json_bytes(schedule_records)),
                    "scoring": "eligible after complete success"
                    if selected_is_canonical_full_fixture
                    else "disabled for subset/pilot/custom manifest",
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    run_dir = args.outdir.resolve() / args.run_id
    config_sha, _sources_sha = _prepare_run(
        run_dir=run_dir,
        run_id=args.run_id,
        config=config,
        source_records=source_records,
        schedule_records=schedule_records,
    )
    print(
        f"Sol reasoning evaluation: {len(papers)} papers x {len(efforts)} efforts "
        f"= {len(schedule)} tasks; workers={args.workers}; run={run_dir}"
    )
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        results = list(
            executor.map(
                lambda task: _execute_task(
                    task,
                    run_dir=run_dir,
                    model=args.model,
                    config_sha=config_sha,
                    sol_module=sol_module,
                ),
                schedule,
            )
        )

    comparison: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "run_id": args.run_id,
        "model": args.model,
        "fixture_sha256": config["fixture_sha256"],
        "canonical_full_fixture": selected_is_canonical_full_fixture,
        "papers_total": len(papers),
        "efforts": {},
        "completed_at": _utc_now(),
    }
    any_failures = False
    for effort in efforts:
        effort_results = [result for result in results if result.effort == effort]
        successful = sum(result.status == "success" for result in effort_results)
        telemetry_rows = _effort_telemetry_rows(effort, papers, run_dir)
        score_summary: dict[str, Any] | None = None
        per_paper: list[dict[str, Any]] = []
        scoring_reason: str | None = None
        if successful != len(papers):
            any_failures = True
            scoring_reason = (
                f"{len(papers) - successful}/{len(papers)} extraction tasks failed; "
                "failures remain in the denominator"
            )
        elif not selected_is_canonical_full_fixture:
            scoring_reason = (
                "subset/pilot/custom manifest is incomplete; scoring disabled"
            )
        else:
            try:
                score_summary, per_paper = score_complete_effort(
                    effort,
                    papers,
                    selected,
                    run_dir=run_dir,
                    model=args.model,
                    config_sha=config_sha,
                )
            except (Exception, SystemExit) as exc:
                any_failures = True
                scoring_reason = (
                    f"migration/scoring failed: {type(exc).__name__}: {exc}"
                )

        merged = merge_paper_metrics(telemetry_rows, per_paper)
        effort_root = run_dir / "efforts" / effort
        _atomic_write_json(effort_root / "paper_metrics.json", merged)
        effort_summary = _summarize_effort(
            effort,
            telemetry_rows,
            prices=prices,
            score_summary=score_summary,
            scoring_reason=scoring_reason,
        )
        _atomic_write_json(effort_root / "summary.json", effort_summary)
        comparison["efforts"][effort] = effort_summary
        print(
            f"  {effort:6s}: {successful}/{len(papers)} extraction successes; "
            + (
                "scored"
                if score_summary is not None
                else f"not scored ({scoring_reason})"
            )
        )

    _atomic_write_json(run_dir / "comparison.json", comparison)
    print(f"Comparison: {run_dir / 'comparison.json'}")
    return 1 if any_failures else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SolEvalError as exc:
        raise SystemExit(f"error: {exc}") from exc
