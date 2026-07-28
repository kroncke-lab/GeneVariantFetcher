"""Self-contained, per-paper HTML viewer for GVF LLM traces.

Packaged (not in ``scripts/``) because ``cli/gvf_run.py`` and
``cli/automated_workflow.py`` build the report on every run and ``scripts`` is
excluded from the wheel — importing it from there made the tracing feature's
headline artifact silently dead in exactly the deployment an operator uses.
``scripts/build_llm_trace_html.py`` remains as a thin CLI wrapper.

Design constraints, in order:

* **Honest integrity.** The banner names the level it actually has
  (``generated_now`` / ``write_time_verified`` / ``locked``) instead of printing
  "verified" for a manifest generated moments ago from the files it validates.
* **Bounded size.** Trace bodies are full prompts; a 40-paper run measured
  0.49 MB of HTML per paper, i.e. ~400 MB at 800 papers. Embedded bodies are
  capped and long runs shard per paper. Nothing is dropped silently: every
  omission is listed in the report and logged.
* **Followable retries.** Accepted / discarded / repaired attempts are labelled
  and cross-linked, so a curator can tell which of up to ~16 identical-looking
  records produced the stored result.
* **XSS-safe embedding.** ``</`` is escaped in the JSON blob and every value is
  escaped again at render time.
"""

from __future__ import annotations

import html
import json
import logging
import os
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

from utils.llm_trace import (
    TRACE_INDEX_NAME,
    TRACE_MANIFEST_NAME,
    build_trace_manifest,
    integrity_level,
    sha256_bytes,
    validate_trace_manifest,
)

logger = logging.getLogger(__name__)

TRACE_REPORT_NAME = "llm_trace_report.html"
REPORT_SCHEMA_VERSION = 2

#: Per-string embedding cap. Above this a body is truncated with a visible
#: marker naming the on-disk record that holds the full text.
DEFAULT_MAX_FIELD_CHARS = 24_000
#: Above this many trace groups the report shards: an index page plus one file
#: per paper. Keeps any single file openable in a browser.
DEFAULT_MAX_PAPERS_PER_FILE = 60
#: Hard cap on records embedded for one paper. Beyond it the remainder is listed
#: in the omissions panel with their record paths.
DEFAULT_MAX_RECORDS_PER_PAPER = 400

#: Stages that should each produce at least one normalized decision event.
#: Used for the coverage panel: a stage with calls but no decision is a real gap
#: in adjudicability, not a cosmetic one.
EXPECTED_DECISION_EVENTS: dict[str, str] = {
    "tier2_relevance_filter": "tier2_relevance_decision",
    "clinical_data_triage": "clinical_data_triage_decision",
    "extraction_priority_triage": "extraction_priority_triage_decision",
    "paper_variant_extraction": "paper_extraction_selection",
    "clinical_table_routing": "table_routing_decision",
    "variant_claim_verification": "variant_claim_verification_decision",
    "variant_claim_debate": "claim_debate_decision",
    "paper_final_check": "paper_final_check_decision",
    "paper_source_grounded_summary": "paper_source_grounded_summary_decision",
    "count_recovery": "count_recovery_decision",
    "paper_relevance": "paper_relevance_decision",
    "synonym_relevance": "synonym_relevance_decision",
    "figure_text_extraction": "figure_text_extraction_decision",
    "figure_variant_read": "figure_variant_read_decision",
    "pedigree_detection": "pedigree_detection_decision",
    "pedigree_extraction": "pedigree_extraction_decision",
    "representation_route": "representation_route_decision",
    "paper_curation": "paper_curation_decision",
}
#: Stages whose decision event is deliberately folded into a parent stage's
#: event, so a missing event of their own is not a coverage gap.
DECISION_FOLDED_INTO_PARENT: dict[str, str] = {
    "paper_extraction_continuation": "paper_extraction_selection",
    "paper_extraction_adjudication": "paper_extraction_selection",
}
#: Reverse map: which stage a decision event describes. ``record_trace_event``
#: defaults a decision record's ``stage`` to the event type, so without this a
#: decision would be filed under a stage of its own name and never satisfy the
#: stage whose calls it explains.
_EVENT_OWNER_STAGE: dict[str, str] = {
    event: stage for stage, event in EXPECTED_DECISION_EVENTS.items()
}

_SAFE_FILE_RE = re.compile(r"[^A-Za-z0-9._-]+")


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def _paper_key(gene: Any, pmid: Any) -> str:
    gene_text = str(gene or "").strip().upper()
    pmid_text = str(pmid or "").strip()
    if pmid_text:
        return f"{gene_text or 'UNKNOWN'}|{pmid_text}"
    if gene_text:
        return f"{gene_text}|_unscoped"
    return "_unscoped"


def _safe_filename(value: str) -> str:
    cleaned = _SAFE_FILE_RE.sub("_", value).strip("._")
    return (cleaned[:100] or "group") + ".html"


def _paper_metadata(run_dir: Path | None) -> dict[str, dict[str, Any]]:
    """Read optional paper labels and selections from whatever the run wrote.

    Benchmark runs write ``selection.json``/``predictions.json``; production runs
    write ``source_completeness.json`` and ``RUN_STATUS.json``. Reading only the
    benchmark files meant production reports never rendered a title or a
    source-choice chip, so the documented adjudication step did not apply there.
    """
    metadata: dict[str, dict[str, Any]] = {}
    if run_dir is None:
        return metadata

    for filename in ("selection.json", "predictions.json"):
        path = run_dir / filename
        if not path.is_file():
            continue
        try:
            payload = _read_json(path)
        except (OSError, json.JSONDecodeError):
            continue
        for paper in payload.get("papers") or []:
            if not isinstance(paper, Mapping):
                continue
            key = _paper_key(paper.get("gene"), paper.get("pmid"))
            current = metadata.setdefault(key, {})
            for field in (
                "title",
                "tool",
                "tool_rationale",
                "source_completeness",
                "curation_rationale",
                "source_selection",
                "representations_available",
            ):
                if paper.get(field) not in (None, "", [], {}):
                    current[field] = paper[field]
            if isinstance(paper.get("variants"), list):
                current["variant_count"] = len(paper["variants"])

    _merge_production_metadata(run_dir, metadata)
    return metadata


def _merge_production_metadata(
    run_dir: Path, metadata: dict[str, dict[str, Any]]
) -> None:
    """Fold production ``source_completeness.json`` into the per-paper metadata."""
    path = run_dir / "source_completeness.json"
    if not path.is_file():
        return
    try:
        payload = _read_json(path)
    except (OSError, json.JSONDecodeError):
        return
    entries: Iterable[Any]
    if isinstance(payload, Mapping):
        entries = payload.get("papers") or payload.get("pmids") or []
    elif isinstance(payload, list):
        entries = payload
    else:
        return
    for entry in entries:
        if not isinstance(entry, Mapping):
            continue
        key = _paper_key(
            entry.get("gene") or entry.get("gene_symbol"), entry.get("pmid")
        )
        current = metadata.setdefault(key, {})
        for source_field, target_field in (
            ("title", "title"),
            ("status", "source_status"),
            ("source_file", "source_file"),
            ("chars", "source_chars"),
            ("reason", "source_reason"),
        ):
            value = entry.get(source_field)
            if value not in (None, "", [], {}) and target_field not in current:
                current[target_field] = value


def _record_timestamp(record: Mapping[str, Any]) -> str:
    return str(
        (record.get("request") or {}).get("started_at")
        or (record.get("response") or {}).get("completed_at")
        or record.get("recorded_at")
        or ""
    )


def _usage_total(record: Mapping[str, Any]) -> int:
    usage = (record.get("response") or {}).get("usage") or {}
    if not isinstance(usage, Mapping):
        return 0
    for key in ("total_tokens", "total_token_count"):
        value = usage.get(key)
        if isinstance(value, (int, float)):
            return int(value)
    total = 0
    for key in (
        "input_tokens",
        "output_tokens",
        "prompt_tokens",
        "completion_tokens",
    ):
        value = usage.get(key)
        if isinstance(value, (int, float)):
            total += int(value)
    return total


def _load_manifest(
    trace_root: Path, *, run_id: str | None
) -> tuple[dict[str, Any], list[str], bool]:
    """Return ``(manifest, integrity_errors, generated_now)``.

    ``generated_now`` records that this call *created* the manifest, which is
    exactly the case where "verified" would be a tautology.
    """
    path = trace_root / TRACE_MANIFEST_NAME
    if not path.is_file():
        manifest = build_trace_manifest(
            trace_root, output_path=path, run_id=run_id, allow_mixed_runs=True
        )
        return manifest, validate_trace_manifest(trace_root, manifest), True
    try:
        manifest = _read_json(path)
    except (OSError, json.JSONDecodeError) as exc:
        return {}, [f"Could not read {TRACE_MANIFEST_NAME}: {exc}"], False
    return manifest, validate_trace_manifest(trace_root, manifest), False


def _bound_strings(
    value: Any,
    *,
    limit: int,
    record_path: str,
    omissions: list[dict[str, Any]],
    _depth: int = 0,
) -> Any:
    """Cap every embedded string, with a visible marker naming the source record.

    Never silent: the marker states how much was cut and where the full text is,
    and the omission is appended for the report's omissions panel and the log.
    """
    if _depth > 40:
        return value
    if isinstance(value, str):
        if len(value) <= limit:
            return value
        omissions.append(
            {
                "kind": "body_truncated",
                "record": record_path,
                "characters_total": len(value),
                "characters_embedded": limit,
                "sha256": sha256_bytes(value.encode("utf-8")),
            }
        )
        dropped = len(value) - limit
        return (
            value[:limit]
            + f"\n\n[… truncated {dropped:,} of {len(value):,} characters."
            + f" Full text: {record_path} …]"
        )
    if isinstance(value, Mapping):
        return {
            key: _bound_strings(
                item,
                limit=limit,
                record_path=record_path,
                omissions=omissions,
                _depth=_depth + 1,
            )
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [
            _bound_strings(
                item,
                limit=limit,
                record_path=record_path,
                omissions=omissions,
                _depth=_depth + 1,
            )
            for item in value
        ]
    return value


def _attempt_role(record: Mapping[str, Any]) -> str:
    context = record.get("context") or {}
    operation = str(context.get("operation") or "")
    if operation:
        return operation
    return str(record.get("record_type") or "")


def _stage_coverage(records: list[Mapping[str, Any]]) -> dict[str, Any]:
    """Per-stage call/decision/accepted-link counts, plus the expected-event gaps.

    Matching is **within the stage**. The prior version compared a stage's
    expected event against a repo-wide set of every event type seen anywhere, so
    one paper that emitted ``count_recovery_decision`` silently satisfied the
    ``count_recovery`` stage of a different paper that emitted nothing.
    """
    stages: dict[str, dict[str, Any]] = {}
    event_types: set[str] = set()
    for record in records:
        context = record.get("context") or {}
        record_type = record.get("record_type")
        if record_type == "decision_event":
            # A decision is filed under the stage it DESCRIBES. record_trace_event
            # defaults `stage` to the event type, so key on the event type and map
            # it back to its owning stage.
            event = str((record.get("event") or {}).get("type") or "")
            stage = _EVENT_OWNER_STAGE.get(event) or str(
                context.get("stage") or event or "unknown"
            )
        else:
            stage = str(context.get("stage") or record_type or "unknown")
        entry = stages.setdefault(
            stage,
            {
                "stage": stage,
                "llm_calls": 0,
                "decisions": 0,
                "failures": 0,
                "accepted_links": 0,
                "expected_event": EXPECTED_DECISION_EVENTS.get(stage),
                "folded_into": DECISION_FOLDED_INTO_PARENT.get(stage),
                "event_types": [],
            },
        )
        if record_type == "llm_call":
            entry["llm_calls"] += 1
            if (record.get("response") or {}).get("success") is False:
                entry["failures"] += 1
        elif record_type == "decision_event":
            entry["decisions"] += 1
            event = str((record.get("event") or {}).get("type") or "")
            if event:
                event_types.add(event)
                if event not in entry["event_types"]:
                    entry["event_types"].append(event)
            data = (record.get("event") or {}).get("data") or {}
            if isinstance(data, Mapping) and (
                data.get("accepted_response_trace_id")
                or data.get("accepted_response_trace_ids")
            ):
                entry["accepted_links"] += 1

    missing: list[dict[str, Any]] = []
    for stage, entry in stages.items():
        if entry["folded_into"]:
            continue
        expected = entry["expected_event"]
        if not entry["llm_calls"]:
            continue
        if expected and expected not in entry["event_types"]:
            missing.append(
                {
                    "stage": stage,
                    "expected_event": expected,
                    "llm_calls": entry["llm_calls"],
                    "reason": "stage made model calls but emitted no decision event",
                }
            )
        elif not expected and not entry["decisions"]:
            missing.append(
                {
                    "stage": stage,
                    "expected_event": "(unregistered stage)",
                    "llm_calls": entry["llm_calls"],
                    "reason": (
                        "stage made model calls and has no registered decision "
                        "event — its route is not adjudicable"
                    ),
                }
            )
    # A stage whose decisions exist but never link an accepted call is a distinct
    # gap: the "why" is recorded, the "from which exact call" is not.
    for stage, entry in stages.items():
        if entry["decisions"] and not entry["accepted_links"]:
            missing.append(
                {
                    "stage": stage,
                    "expected_event": entry["expected_event"] or "(any)",
                    "llm_calls": entry["llm_calls"],
                    "reason": "decision events carry no accepted_response_trace_id",
                }
            )
    return {
        "stages": sorted(stages.values(), key=lambda item: item["stage"]),
        "event_types": sorted(event_types),
        "missing_decision_links": sorted(
            missing, key=lambda item: (item["stage"], item["reason"])
        ),
    }


def _relative_trace_href(output_path: Path | None, trace_root: Path) -> str | None:
    """POSIX-relative path from the report to the trace root, when one exists.

    Only meaningful when the report sits on the same filesystem tree as the
    records, which is the normal ``<run_dir>/llm_trace_report.html`` +
    ``<run_dir>/llm_traces/`` layout. With ``GVF_LLM_TRACE_DIR`` on a separate
    volume there is no sane relative path, so the viewer shows the path as text
    instead of emitting a broken link.
    """
    if output_path is None:
        return None
    try:
        relative = os.path.relpath(trace_root, start=output_path.parent)
    except (OSError, ValueError):
        return None
    posix = Path(relative).as_posix()
    # Refuse anything that climbs out more than one level: it is almost certainly
    # a cross-volume path, and a long ../../.. chain leaks the layout.
    if posix.count("../") > 1 or posix.startswith("/"):
        return None
    return posix


def collect_trace_report_data(
    trace_root: Path | str,
    *,
    run_dir: Path | str | None = None,
    title: str | None = None,
    run_id: str | None = None,
    max_field_chars: int = DEFAULT_MAX_FIELD_CHARS,
    max_records_per_paper: int = DEFAULT_MAX_RECORDS_PER_PAPER,
    locked: bool = False,
    output_path: Path | str | None = None,
) -> dict[str, Any]:
    """Collect trace JSON into a browser-friendly per-paper data structure."""
    root = Path(trace_root).expanduser().resolve()
    resolved_run_dir = Path(run_dir).expanduser().resolve() if run_dir else None
    manifest, integrity_errors, generated_now = _load_manifest(root, run_id=run_id)
    metadata = _paper_metadata(resolved_run_dir)
    grouped: dict[str, dict[str, Any]] = {}
    omissions: list[dict[str, Any]] = []

    record_paths = (
        sorted(
            path for path in root.rglob("*.json") if path.name != TRACE_MANIFEST_NAME
        )
        if root.is_dir()
        else []
    )
    for path in record_paths:
        relative = str(path.relative_to(root))
        try:
            record = _read_json(path)
        except (OSError, json.JSONDecodeError) as exc:
            record = {
                "schema_version": None,
                "record_type": "unreadable",
                "trace_id": path.stem,
                "context": {},
                "recorded_at": "",
                "read_error": str(exc),
            }
        record = dict(record)
        record["_report_path"] = relative
        context = record.get("context") or {}
        gene = str(context.get("gene") or "").upper()
        pmid = str(context.get("pmid") or "")
        key = _paper_key(gene, pmid)
        paper = grouped.setdefault(
            key,
            {
                "id": key,
                "gene": gene or ("Other" if key == "_unscoped" else "UNKNOWN"),
                "pmid": pmid,
                "title": "",
                "metadata": {},
                "records": [],
            },
        )
        paper["records"].append(record)

    all_records: list[Mapping[str, Any]] = [
        record for paper in grouped.values() for record in paper["records"]
    ]
    coverage = _stage_coverage(all_records)
    total_record_count = len(all_records)

    for key, paper in grouped.items():
        paper["records"].sort(
            key=lambda record: (
                _record_timestamp(record),
                str(record.get("trace_id") or ""),
            )
        )
        paper["metadata"] = metadata.get(key, {})
        paper["title"] = paper["metadata"].get("title", "")
        paper["call_count"] = sum(
            record.get("record_type") == "llm_call" for record in paper["records"]
        )
        paper["decision_count"] = sum(
            record.get("record_type") == "decision_event" for record in paper["records"]
        )
        paper["failure_count"] = sum(
            (record.get("response") or {}).get("success") is False
            for record in paper["records"]
        )
        paper["token_count"] = sum(_usage_total(record) for record in paper["records"])
        paper["models"] = sorted(
            {
                str((record.get("request") or {}).get("requested_model"))
                for record in paper["records"]
                if (record.get("request") or {}).get("requested_model")
            }
        )
        paper["run_ids"] = sorted(
            {
                str(record.get("run_id"))
                for record in paper["records"]
                if record.get("run_id")
            }
        )
        paper["stages"] = sorted(
            {
                str((record.get("context") or {}).get("stage"))
                for record in paper["records"]
                if (record.get("context") or {}).get("stage")
            }
        )
        paper["accepted_trace_ids"] = sorted(
            {
                str(
                    ((record.get("event") or {}).get("data") or {}).get(
                        "accepted_response_trace_id"
                    )
                )
                for record in paper["records"]
                if isinstance((record.get("event") or {}).get("data"), Mapping)
                and ((record.get("event") or {}).get("data") or {}).get(
                    "accepted_response_trace_id"
                )
            }
        )
        paper["discarded_trace_ids"] = sorted(
            {
                str(trace_id)
                for record in paper["records"]
                if isinstance((record.get("event") or {}).get("data"), Mapping)
                for trace_id in (
                    ((record.get("event") or {}).get("data") or {}).get(
                        "discarded_trace_ids"
                    )
                    or []
                )
            }
        )

        if len(paper["records"]) > max_records_per_paper:
            dropped = paper["records"][max_records_per_paper:]
            omissions.append(
                {
                    "kind": "records_omitted",
                    "group": key,
                    "records_total": len(paper["records"]),
                    "records_embedded": max_records_per_paper,
                    "record_paths": [
                        str(record.get("_report_path") or "") for record in dropped
                    ][:200],
                }
            )
            paper["records"] = paper["records"][:max_records_per_paper]

        bounded: list[dict[str, Any]] = []
        for record in paper["records"]:
            relative = str(record.get("_report_path") or "")
            before = len(omissions)
            trimmed = _bound_strings(
                record,
                limit=max_field_chars,
                record_path=relative,
                omissions=omissions,
            )
            # Flag the card so the viewer can point at the complete record.
            trimmed["_truncated_fields"] = len(omissions) - before
            bounded.append(trimmed)
        paper["records"] = bounded
        for record in paper["records"]:
            record["_attempt_role"] = _attempt_role(record)
            record["_is_accepted"] = str(record.get("trace_id") or "") in set(
                paper["accepted_trace_ids"]
            )
            record["_is_discarded"] = str(record.get("trace_id") or "") in set(
                paper["discarded_trace_ids"]
            )

    papers = sorted(
        grouped.values(),
        key=lambda paper: (
            not bool(paper["pmid"]),
            paper["gene"],
            paper["pmid"],
        ),
    )
    level, level_description = integrity_level(
        manifest, list(integrity_errors), locked=locked
    )
    if generated_now and not integrity_errors:
        # The manifest did not exist before this call: it cannot attest to
        # anything beyond internal consistency, whatever the record digests say.
        level, level_description = (
            "generated_now",
            "Manifest was generated by this report from the files on disk. It "
            "proves internal consistency, not that records are unmodified since "
            "they were written.",
        )
    generated_at = datetime.now(timezone.utc).isoformat()
    for omission in omissions:
        logger.info("trace report omission: %s", json.dumps(omission, sort_keys=True))
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "mode": "single",
        "title": title or "LLM Curation Trace Review",
        "generated_at": generated_at,
        "trace_root": str(root),
        "trace_href": _relative_trace_href(
            Path(output_path) if output_path else None, root
        ),
        "run_id": run_id or manifest.get("run_id"),
        "record_run_ids": manifest.get("record_run_ids") or [],
        "integrity": {
            "manifest_present": bool(manifest),
            "manifest_generated_at": manifest.get("generated_at"),
            "manifest_generated_by_this_report": generated_now,
            "level": level,
            "level_description": level_description,
            "valid": not integrity_errors,
            "errors": list(integrity_errors),
            "trace_count": manifest.get("trace_count", total_record_count),
            "index_name": TRACE_INDEX_NAME,
            "index_present": bool((manifest.get("index") or {}).get("present")),
            "write_time_verified_records": (manifest.get("verification") or {}).get(
                "write_time_verified_records"
            ),
        },
        "coverage": coverage,
        "size_policy": {
            "max_field_chars": max_field_chars,
            "max_records_per_paper": max_records_per_paper,
            "note": (
                "Embedded bodies are capped so the report stays openable. Every "
                "truncation states how much was cut and names the on-disk trace "
                "record with the full text."
            ),
        },
        "omissions": omissions,
        "summary": {
            "paper_count": sum(bool(paper["pmid"]) for paper in papers),
            "group_count": len(papers),
            "trace_count": total_record_count,
            "llm_call_count": sum(
                record.get("record_type") == "llm_call" for record in all_records
            ),
            "decision_event_count": sum(
                record.get("record_type") == "decision_event" for record in all_records
            ),
            "failure_count": sum(
                (record.get("response") or {}).get("success") is False
                for record in all_records
            ),
            "token_count": sum(_usage_total(record) for record in all_records),
        },
        "reasoning_limit": manifest.get("reasoning_limit")
        or (
            "This report shows explicit prompts, outputs, rationales, and "
            "provider-exposed reasoning content. Private hidden chain-of-thought "
            "is not exposed by model APIs and is not recorded. Reasoning token "
            "counts are billing telemetry, not exposed reasoning."
        ),
        "papers": papers,
    }


def _embedded_trace_ids(data: Mapping[str, Any]) -> list[str]:
    """Trace ids whose records are embedded in THIS payload.

    Computed in Python, per rendered file, so the viewer never has to guess: a
    jump button is only offered for a target that is actually on the page.
    Sharding and the record cap both remove records, and a decision may name a
    trace id that was never recorded at all (``None``) — all three used to
    produce dead buttons.
    """
    return sorted(
        {
            str(record["trace_id"])
            for paper in (data.get("papers") or [])
            for record in (paper.get("records") or [])
            if record.get("trace_id")
        }
    )


def render_trace_report(data: Mapping[str, Any]) -> str:
    """Render a complete local HTML application with embedded trace data."""
    data = {**data, "embedded_trace_ids": _embedded_trace_ids(data)}
    # "</" is broken so no embedded body can close the <script> element, and the
    # two Unicode line separators are escaped because JS treats them as newlines
    # inside a string literal even though JSON does not.
    data_json = (
        json.dumps(data, ensure_ascii=False, separators=(",", ":"), default=str)
        .replace("</", "<\\/")
        .replace(" ", "\\u2028")
        .replace(" ", "\\u2029")
    )
    title = html.escape(str(data.get("title") or "LLM Curation Trace Review"))
    return _HTML_TEMPLATE.replace("__TITLE__", title).replace("__DATA__", data_json)


def _shard_payload(
    data: Mapping[str, Any], paper: Mapping[str, Any], index_href: str
) -> dict[str, Any]:
    payload = {key: value for key, value in data.items() if key != "papers"}
    payload["mode"] = "shard"
    payload["index_href"] = index_href
    payload["papers"] = [paper]
    payload["title"] = f"{data.get('title')} · {paper.get('gene')} {paper.get('pmid')}"
    # Shards live one directory deeper than the index, so the relative path to
    # the trace root needs one more level up. Drop the link if that would climb
    # out of a same-tree layout.
    href = data.get("trace_href")
    if href:
        candidate = f"../{href}"
        payload["trace_href"] = None if candidate.count("../") > 2 else candidate
    return payload


def _index_payload(
    data: Mapping[str, Any], shard_hrefs: Mapping[str, str]
) -> dict[str, Any]:
    payload = {key: value for key, value in data.items() if key != "papers"}
    payload["mode"] = "index"
    payload["papers"] = [
        {
            **{key: value for key, value in paper.items() if key not in ("records",)},
            "records": [],
            "href": shard_hrefs[paper["id"]],
        }
        for paper in data.get("papers") or []
    ]
    return payload


def build_trace_html_report(
    trace_root: Path | str,
    *,
    output_path: Path | str,
    run_dir: Path | str | None = None,
    title: str | None = None,
    run_id: str | None = None,
    max_field_chars: int = DEFAULT_MAX_FIELD_CHARS,
    max_papers_per_file: int = DEFAULT_MAX_PAPERS_PER_FILE,
    max_records_per_paper: int = DEFAULT_MAX_RECORDS_PER_PAPER,
    locked: bool = False,
) -> dict[str, Any]:
    """Build the report and return its data summary.

    Runs with more than ``max_papers_per_file`` trace groups are written as an
    index page plus one file per group under ``<output>_papers/``, so no single
    HTML file grows past a size a browser can open. The index still carries the
    full summary, integrity level, stage coverage, and omission list.
    """
    destination = Path(output_path)
    data = collect_trace_report_data(
        trace_root,
        run_dir=run_dir,
        title=title,
        run_id=run_id,
        max_field_chars=max_field_chars,
        max_records_per_paper=max_records_per_paper,
        locked=locked,
        output_path=destination,
    )
    papers = data.get("papers") or []
    if len(papers) <= max_papers_per_file:
        _atomic_write_text(destination, render_trace_report(data))
        data["output_files"] = [str(destination)]
        data["sharded"] = False
        return data

    shard_dir = destination.with_name(f"{destination.stem}_papers")
    shard_hrefs: dict[str, str] = {}
    used: set[str] = set()
    for paper in papers:
        name = _safe_filename(str(paper["id"]))
        while name in used:
            name = f"{Path(name).stem}_{len(used)}.html"
        used.add(name)
        shard_hrefs[paper["id"]] = f"{shard_dir.name}/{name}"

    index_href = f"../{destination.name}"
    written: list[str] = []
    for paper in papers:
        shard_path = shard_dir / Path(shard_hrefs[paper["id"]]).name
        _atomic_write_text(
            shard_path, render_trace_report(_shard_payload(data, paper, index_href))
        )
        written.append(str(shard_path))
    index_data = _index_payload(data, shard_hrefs)
    _atomic_write_text(destination, render_trace_report(index_data))
    logger.info(
        "trace report sharded: %d group(s) across %s (index %s)",
        len(papers),
        shard_dir,
        destination,
    )
    data["output_files"] = [str(destination), *written]
    data["sharded"] = True
    data["shard_dir"] = str(shard_dir)
    return data


_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en" data-theme="light">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="light dark">
  <title>__TITLE__</title>
  <style>
    :root {
      --bg: #f3f1eb;
      --surface: #fffefa;
      --surface-2: #f8f6f0;
      --ink: #17212b;
      --muted: #66727d;
      --line: #d9d7cf;
      --navy: #15334a;
      --teal: #087f76;
      --teal-soft: #dff3ef;
      --amber: #aa6812;
      --amber-soft: #fff0d6;
      --red: #a73d3d;
      --red-soft: #fbe5e2;
      --violet: #6556a8;
      --violet-soft: #ece8fb;
      --shadow: 0 18px 50px rgba(31, 42, 50, .08);
      --mono: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
      --sans: Inter, ui-sans-serif, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }
    html[data-theme="dark"] {
      --bg: #11181d;
      --surface: #182229;
      --surface-2: #202c34;
      --ink: #edf2f3;
      --muted: #a5b0b7;
      --line: #34434c;
      --navy: #cde7f5;
      --teal: #56c9bd;
      --teal-soft: #173d3b;
      --amber: #f0b55c;
      --amber-soft: #48351d;
      --red: #f29a94;
      --red-soft: #492927;
      --violet: #b9acef;
      --violet-soft: #302c4c;
      --shadow: 0 20px 55px rgba(0, 0, 0, .26);
    }
    * { box-sizing: border-box; }
    html { scroll-behavior: smooth; }
    body {
      margin: 0;
      background:
        radial-gradient(circle at 8% -10%, rgba(8,127,118,.12), transparent 28rem),
        var(--bg);
      color: var(--ink);
      font: 14px/1.55 var(--sans);
    }
    button, input, select { font: inherit; }
    button { color: inherit; }
    a { color: var(--teal); }
    .topbar {
      min-height: 76px;
      padding: 14px 22px;
      display: flex;
      align-items: center;
      gap: 18px;
      border-bottom: 1px solid var(--line);
      background: color-mix(in srgb, var(--surface) 92%, transparent);
      backdrop-filter: blur(18px);
      position: sticky;
      top: 0;
      z-index: 20;
    }
    .brand { min-width: 240px; }
    .eyebrow {
      color: var(--teal);
      font-size: 10px;
      font-weight: 800;
      letter-spacing: .16em;
      text-transform: uppercase;
    }
    h1 { margin: 2px 0 0; font-size: 20px; line-height: 1.15; letter-spacing: -.025em; }
    .run-id { margin-top: 2px; color: var(--muted); font: 10px var(--mono); }
    .summary-strip { display: flex; gap: 8px; flex: 1; flex-wrap: wrap; }
    .metric {
      min-width: 82px;
      border-left: 1px solid var(--line);
      padding-left: 12px;
    }
    .metric strong { display: block; color: var(--navy); font-size: 16px; }
    .metric span { color: var(--muted); font-size: 11px; }
    .top-actions { display: flex; gap: 7px; }
    .icon-button, .button, .filter-button, .copy-button {
      border: 1px solid var(--line);
      background: var(--surface);
      border-radius: 9px;
      cursor: pointer;
      transition: border-color .15s, transform .15s, background .15s;
    }
    .icon-button { width: 38px; height: 38px; }
    .icon-button:hover, .button:hover, .filter-button:hover, .copy-button:hover {
      border-color: var(--teal);
      transform: translateY(-1px);
    }
    .icon-button:focus-visible, .button:focus-visible,
    .filter-button:focus-visible, .copy-button:focus-visible,
    .paper-button:focus-visible, .jump-button:focus-visible {
      outline: 3px solid color-mix(in srgb, var(--teal) 35%, transparent);
      outline-offset: 2px;
    }
    .layout {
      display: grid;
      grid-template-columns: 306px minmax(0, 1fr);
      min-height: calc(100vh - 77px);
    }
    aside {
      border-right: 1px solid var(--line);
      background: color-mix(in srgb, var(--surface) 76%, transparent);
      padding: 18px 14px;
      position: sticky;
      top: 77px;
      height: calc(100vh - 77px);
      overflow: auto;
    }
    .search {
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 10px;
      background: var(--surface);
      color: var(--ink);
      padding: 10px 12px;
      outline: none;
    }
    .search:focus { border-color: var(--teal); box-shadow: 0 0 0 3px rgba(8,127,118,.12); }
    .nav-label {
      margin: 17px 8px 8px;
      color: var(--muted);
      font-size: 10px;
      font-weight: 800;
      letter-spacing: .13em;
      text-transform: uppercase;
    }
    .paper-nav { display: grid; gap: 6px; }
    .paper-button {
      width: 100%;
      padding: 10px 11px;
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 2px 10px;
      text-align: left;
      border: 1px solid transparent;
      border-radius: 10px;
      background: transparent;
      color: var(--ink);
      cursor: pointer;
      text-decoration: none;
    }
    .paper-button:hover { background: var(--surface-2); }
    .paper-button.active {
      background: var(--surface);
      border-color: var(--line);
      box-shadow: 0 8px 24px rgba(31,42,50,.06);
    }
    .paper-id { font-weight: 750; }
    .paper-title {
      grid-column: 1 / -1;
      color: var(--muted);
      font-size: 11px;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }
    .nav-count {
      align-self: center;
      color: var(--muted);
      font: 10px var(--mono);
    }
    .failure-dot { color: var(--red); }
    main { min-width: 0; padding: 28px clamp(18px, 4vw, 54px) 70px; }
    .paper-header {
      max-width: 1120px;
      margin: 0 auto 20px;
      display: grid;
      gap: 14px;
    }
    .paper-kicker { color: var(--teal); font: 12px var(--mono); }
    .paper-header h2 {
      margin: 0;
      color: var(--navy);
      font-size: clamp(25px, 3vw, 38px);
      line-height: 1.06;
      letter-spacing: -.035em;
    }
    .paper-subtitle { margin: 0; color: var(--muted); max-width: 850px; }
    .source-note {
      max-width: 900px;
      padding: 10px 12px;
      border: 1px solid color-mix(in srgb, var(--amber) 35%, var(--line));
      border-radius: 10px;
      background: var(--amber-soft);
      color: var(--ink);
      font-size: 12px;
    }
    .source-note strong { color: var(--amber); }
    .chips { display: flex; gap: 7px; flex-wrap: wrap; }
    .chip {
      display: inline-flex;
      align-items: center;
      gap: 5px;
      padding: 5px 8px;
      border-radius: 999px;
      background: var(--surface);
      border: 1px solid var(--line);
      color: var(--muted);
      font-size: 11px;
    }
    .chip strong { color: var(--ink); }
    .chip.accepted { background: var(--teal-soft); border-color: var(--teal); color: var(--teal); }
    .chip.discarded { background: var(--amber-soft); border-color: var(--amber); color: var(--amber); }
    .panel-block {
      max-width: 1120px;
      margin: 0 auto 18px;
      padding: 11px 13px;
      border: 1px solid var(--line);
      border-radius: 11px;
      color: var(--muted);
      background: var(--surface);
    }
    .panel-block h3 {
      margin: 0 0 8px;
      color: var(--ink);
      font-size: 11px;
      letter-spacing: .12em;
      text-transform: uppercase;
    }
    .integrity {
      max-width: 1120px;
      margin: 0 auto 18px;
      padding: 11px 13px;
      border: 1px solid var(--line);
      border-radius: 11px;
      color: var(--muted);
      background: var(--surface);
      display: flex;
      align-items: flex-start;
      gap: 9px;
    }
    .integrity.good { border-color: color-mix(in srgb, var(--teal) 45%, var(--line)); }
    .integrity.partial { background: var(--amber-soft); border-color: var(--amber); color: var(--ink); }
    .integrity.bad { background: var(--red-soft); border-color: var(--red); color: var(--ink); }
    .coverage-table { width: 100%; border-collapse: collapse; font-size: 12px; }
    .coverage-table th, .coverage-table td {
      padding: 5px 8px;
      border-bottom: 1px solid var(--line);
      text-align: left;
    }
    .coverage-table th { color: var(--muted); font-size: 10px; text-transform: uppercase; letter-spacing: .1em; }
    .coverage-table td.gap { color: var(--amber); font-weight: 700; }
    .toolbar {
      max-width: 1120px;
      margin: 0 auto 18px;
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
    }
    .toolbar .search { width: min(360px, 100%); }
    .filter-button { padding: 8px 11px; color: var(--muted); }
    .filter-button.active { background: var(--navy); border-color: var(--navy); color: var(--surface); }
    .timeline {
      max-width: 1120px;
      margin: 0 auto;
      position: relative;
      display: grid;
      gap: 15px;
    }
    .timeline::before {
      content: "";
      position: absolute;
      left: 17px;
      top: 8px;
      bottom: 8px;
      width: 1px;
      background: var(--line);
    }
    .trace-card {
      margin-left: 39px;
      position: relative;
      border: 1px solid var(--line);
      border-radius: 14px;
      background: var(--surface);
      box-shadow: var(--shadow);
      overflow: hidden;
    }
    .trace-card::before {
      content: "";
      position: absolute;
      left: -29px;
      top: 22px;
      width: 9px;
      height: 9px;
      border: 4px solid var(--bg);
      border-radius: 50%;
      background: var(--teal);
      box-sizing: content-box;
    }
    .trace-card.decision::before { background: var(--violet); }
    .trace-card.failed::before { background: var(--red); }
    .trace-card.accepted { border-color: color-mix(in srgb, var(--teal) 55%, var(--line)); }
    .trace-card.discarded { opacity: .86; }
    .trace-head {
      padding: 14px 16px 12px;
      display: flex;
      align-items: flex-start;
      gap: 12px;
      border-bottom: 1px solid var(--line);
    }
    .step {
      min-width: 34px;
      height: 24px;
      display: grid;
      place-items: center;
      border-radius: 7px;
      background: var(--teal-soft);
      color: var(--teal);
      font: 700 11px var(--mono);
    }
    .decision .step { background: var(--violet-soft); color: var(--violet); }
    .failed .step { background: var(--red-soft); color: var(--red); }
    .trace-heading { min-width: 0; flex: 1; }
    .trace-heading h3 { margin: 0; font-size: 15px; letter-spacing: -.01em; }
    .trace-meta { margin-top: 3px; color: var(--muted); font-size: 11px; }
    .trace-id { font: 10px var(--mono); color: var(--muted); }
    .trace-body { padding: 15px 16px 16px; display: grid; gap: 12px; }
    .call-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }
    .panel {
      min-width: 0;
      border: 1px solid var(--line);
      border-radius: 10px;
      background: var(--surface-2);
      overflow: hidden;
    }
    .panel-title {
      padding: 8px 11px;
      border-bottom: 1px solid var(--line);
      display: flex;
      justify-content: space-between;
      align-items: center;
      color: var(--muted);
      font-size: 10px;
      font-weight: 800;
      letter-spacing: .11em;
      text-transform: uppercase;
    }
    pre {
      margin: 0;
      padding: 12px;
      max-height: 360px;
      overflow: auto;
      white-space: pre-wrap;
      overflow-wrap: anywhere;
      color: var(--ink);
      font: 12px/1.55 var(--mono);
      tab-size: 2;
    }
    .decision-grid { display: grid; grid-template-columns: minmax(180px, .65fr) 1.35fr; gap: 12px; }
    .decision-facts { display: grid; grid-template-columns: max-content 1fr; gap: 7px 11px; padding: 12px; }
    .decision-facts dt { color: var(--muted); font-size: 11px; }
    .decision-facts dd { margin: 0; overflow-wrap: anywhere; }
    .rationale {
      padding: 14px;
      border-left: 3px solid var(--violet);
      background: var(--violet-soft);
      border-radius: 4px 10px 10px 4px;
    }
    .rationale strong {
      display: block;
      margin-bottom: 5px;
      color: var(--violet);
      font-size: 10px;
      letter-spacing: .12em;
      text-transform: uppercase;
    }
    .rationale p { margin: 0 0 7px; }
    .rationale p:last-child { margin-bottom: 0; }
    details.raw { border-top: 1px solid var(--line); }
    details.raw > summary {
      padding: 10px 16px;
      color: var(--muted);
      cursor: pointer;
      font-size: 11px;
      user-select: none;
    }
    details.raw[open] > summary { border-bottom: 1px solid var(--line); }
    .copy-button { padding: 3px 7px; color: var(--muted); font-size: 10px; }
    .jump-button {
      border: 0;
      padding: 0;
      background: transparent;
      color: var(--teal);
      text-decoration: underline;
      text-underline-offset: 2px;
      cursor: pointer;
    }
    .empty {
      max-width: 1120px;
      margin: 40px auto;
      padding: 34px;
      border: 1px dashed var(--line);
      border-radius: 14px;
      color: var(--muted);
      text-align: center;
    }
    .reasoning-note {
      max-width: 1120px;
      margin: 22px auto 0;
      color: var(--muted);
      font-size: 11px;
    }
    mark { background: var(--amber-soft); color: var(--ink); border-radius: 2px; }
    @media (max-width: 900px) {
      .topbar { position: static; }
      .summary-strip { display: none; }
      .layout { grid-template-columns: 1fr; }
      aside {
        position: static;
        height: auto;
        border-right: 0;
        border-bottom: 1px solid var(--line);
      }
      .paper-nav { grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); max-height: 240px; overflow: auto; }
      main { padding-top: 22px; }
      .call-grid, .decision-grid { grid-template-columns: 1fr; }
    }
    @media (max-width: 560px) {
      .topbar { padding: 12px 14px; }
      .brand { min-width: 0; flex: 1; }
      .top-actions .print { display: none; }
      .trace-card { margin-left: 27px; }
      .timeline::before { left: 8px; }
      .trace-card::before { left: -26px; }
      main { padding-inline: 12px; }
      .trace-head { padding-inline: 12px; }
      .trace-body { padding-inline: 12px; }
    }
    @media print {
      body { background: #fff; color: #111; }
      .topbar, aside, .toolbar, .copy-button { display: none !important; }
      .layout { display: block; }
      main { padding: 0; }
      .trace-card { break-inside: avoid; box-shadow: none; margin-left: 0; }
      .trace-card::before, .timeline::before { display: none; }
      pre { max-height: none; overflow: visible; font-size: 9px; }
      details.raw { display: none; }
    }
  </style>
</head>
<body>
  <header class="topbar">
    <div class="brand">
      <div class="eyebrow">GeneVariantFetcher · audit trail</div>
      <h1 id="reportTitle"></h1>
      <div class="run-id" id="runId"></div>
    </div>
    <div class="summary-strip" id="summaryStrip"></div>
    <div class="top-actions">
      <button class="icon-button print" type="button" onclick="window.print()" title="Print selected paper" aria-label="Print selected paper">⌘P</button>
      <button class="icon-button" type="button" id="themeButton" onclick="toggleTheme()" title="Toggle color theme" aria-label="Toggle color theme">◐</button>
    </div>
  </header>
  <div class="layout">
    <aside>
      <input class="search" id="paperSearch" type="search" placeholder="Search gene, PMID, title…" aria-label="Search papers">
      <div class="nav-label">Papers and trace groups</div>
      <nav class="paper-nav" id="paperNav" aria-label="Paper traces"></nav>
    </aside>
    <main>
      <section class="paper-header" id="paperHeader"></section>
      <div class="integrity" id="integrity"></div>
      <div id="coverage"></div>
      <div id="omissions"></div>
      <div class="toolbar" id="toolbar">
        <input class="search" id="traceSearch" type="search" placeholder="Search this timeline…" aria-label="Search timeline">
        <button class="filter-button active" data-filter="all" type="button">All</button>
        <button class="filter-button" data-filter="llm_call" type="button">Model calls</button>
        <button class="filter-button" data-filter="decision_event" type="button">Decisions</button>
        <button class="filter-button" data-filter="accepted" type="button">Accepted</button>
        <button class="filter-button" data-filter="retry" type="button">Retries &amp; repairs</button>
        <button class="filter-button" data-filter="failed" type="button">Failures</button>
      </div>
      <section class="timeline" id="timeline"></section>
      <p class="reasoning-note" id="reasoningNote"></p>
    </main>
  </div>
  <script>
    const DATA = __DATA__;
    const state = {
      paperId: null,
      paperQuery: "",
      traceQuery: "",
      filter: "all"
    };

    const $ = (selector) => document.querySelector(selector);
    const escapeHTML = (value) => String(value == null ? "" : value).replace(
      /[&<>"']/g,
      char => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#039;"}[char])
    );
    const pretty = (value) => {
      if (typeof value === "string") {
        try { return JSON.stringify(JSON.parse(value), null, 2); } catch (_) { return value; }
      }
      return JSON.stringify(value == null ? null : value, null, 2);
    };
    const formatNumber = (value) => new Intl.NumberFormat().format(Number(value || 0));
    const shortId = (value) => {
      const text = String(value || "");
      return text.length > 22 ? text.slice(0, 10) + "…" + text.slice(-7) : text;
    };
    const formatTime = (value) => {
      if (!value) return "time unavailable";
      const date = new Date(value);
      return Number.isNaN(date.getTime()) ? value : date.toLocaleString();
    };
    const usageTotal = (record) => {
      const usage = record.response?.usage || {};
      if (Number.isFinite(usage.total_tokens)) return usage.total_tokens;
      if (Number.isFinite(usage.total_token_count)) return usage.total_token_count;
      return ["input_tokens","output_tokens","prompt_tokens","completion_tokens"]
        .reduce((sum, key) => sum + (Number(usage[key]) || 0), 0);
    };
    const recordTime = (record) => record.request?.started_at || record.response?.completed_at || record.recorded_at || "";
    const recordStage = (record) => record.context?.stage || record.event?.type || record.record_type || "trace";
    const recordText = (record) => JSON.stringify(record).toLowerCase();
    const RETRY_ROLES = new Set(["json_repair", "empty_content_retry"]);
    const isRetry = (record) => record._is_discarded === true
      || RETRY_ROLES.has(String(record._attempt_role || ""))
      || Number(record.context?.attempt || 0) > 1;

    function init() {
      $("#reportTitle").textContent = DATA.title;
      $("#runId").textContent = DATA.run_id
        ? `run ${DATA.run_id}${(DATA.record_run_ids || []).length > 1 ? " · MIXED RUN DIRECTORY" : ""}`
        : "run id unavailable";
      $("#reasoningNote").textContent = DATA.reasoning_limit;
      renderSummary();
      const requested = decodeURIComponent(location.hash.replace(/^#/, ""));
      state.paperId = DATA.papers.some(p => p.id === requested) ? requested : DATA.papers[0]?.id || null;
      $("#paperSearch").addEventListener("input", event => {
        state.paperQuery = event.target.value.toLowerCase();
        renderPaperNav();
      });
      $("#traceSearch").addEventListener("input", event => {
        state.traceQuery = event.target.value.toLowerCase();
        renderTimeline();
      });
      document.querySelectorAll(".filter-button").forEach(button => {
        button.addEventListener("click", () => {
          state.filter = button.dataset.filter;
          document.querySelectorAll(".filter-button").forEach(item => item.classList.toggle("active", item === button));
          renderTimeline();
        });
      });
      const preferredTheme = safeStorageGet("gvf-trace-theme");
      if (preferredTheme) document.documentElement.dataset.theme = preferredTheme;
      if (DATA.mode === "index") $("#toolbar").style.display = "none";
      render();
    }

    function renderSummary() {
      const items = [
        [DATA.summary.paper_count, "papers"],
        [DATA.summary.llm_call_count, "model calls"],
        [DATA.summary.decision_event_count, "decisions"],
        [DATA.summary.failure_count, "failures"],
        [formatNumber(DATA.summary.token_count), "tokens"]
      ];
      $("#summaryStrip").innerHTML = items.map(([value, label]) =>
        `<div class="metric"><strong>${escapeHTML(value)}</strong><span>${escapeHTML(label)}</span></div>`
      ).join("");
    }

    function render() {
      renderPaperNav();
      renderPaperHeader();
      renderIntegrity();
      renderCoverage();
      renderOmissions();
      renderTimeline();
    }

    function visiblePapers() {
      if (!state.paperQuery) return DATA.papers;
      return DATA.papers.filter(paper => {
        const haystack = [paper.gene, paper.pmid, paper.title, JSON.stringify(paper.metadata)].join(" ").toLowerCase();
        return haystack.includes(state.paperQuery);
      });
    }

    function renderPaperNav() {
      const papers = visiblePapers();
      $("#paperNav").innerHTML = papers.map(paper => {
        const label = paper.pmid ? `${paper.gene} · PMID ${paper.pmid}` : `${paper.gene} · non-paper`;
        const count = `${paper.call_count}C · ${paper.decision_count}D`;
        const inner = `
          <span class="paper-id">${escapeHTML(label)}</span>
          <span class="nav-count ${paper.failure_count ? "failure-dot" : ""}">${escapeHTML(count)}</span>
          <span class="paper-title">${escapeHTML(paper.title || paper.metadata?.tool_rationale || "Trace group")}</span>`;
        return paper.href
          ? `<a class="paper-button" href="${escapeHTML(paper.href)}">${inner}</a>`
          : `<button class="paper-button ${paper.id === state.paperId ? "active" : ""}" type="button" data-paper="${escapeHTML(paper.id)}">${inner}</button>`;
      }).join("") || `<div class="empty">No papers match that search.</div>`;
      document.querySelectorAll("[data-paper]").forEach(button => {
        button.addEventListener("click", () => {
          state.paperId = button.dataset.paper;
          location.hash = encodeURIComponent(state.paperId);
          state.traceQuery = "";
          $("#traceSearch").value = "";
          render();
          window.scrollTo({top: 0, behavior: "smooth"});
        });
      });
    }

    function selectedPaper() {
      return DATA.papers.find(paper => paper.id === state.paperId);
    }

    function renderPaperHeader() {
      if (DATA.mode === "index") {
        $("#paperHeader").innerHTML = `
          <div class="paper-kicker">${escapeHTML(DATA.summary.group_count)} trace groups · sharded report</div>
          <h2>${escapeHTML(DATA.title)}</h2>
          <p class="paper-subtitle">This run has too many trace groups for one HTML file, so each paper has its own page. Pick a paper on the left. Every page carries the same integrity level and coverage summary.</p>`;
        return;
      }
      const paper = selectedPaper();
      if (!paper) {
        $("#paperHeader").innerHTML = `<div class="empty">No trace records were found.</div>`;
        return;
      }
      const label = paper.pmid ? `PMID ${paper.pmid}` : "Non-paper model activity";
      const title = paper.title || paper.metadata?.curation_rationale || "LLM execution and decision trace";
      const models = paper.models.length ? paper.models.join(", ") : "model unavailable";
      const variantChip = Number.isFinite(paper.metadata?.variant_count)
        ? `<span class="chip"><strong>${paper.metadata.variant_count}</strong> variants selected</span>` : "";
      const source = paper.metadata?.source_selection;
      const sourceName = source?.selected ? String(source.selected).split(/[\\/]/).pop() : "";
      const sourceNote = source ? `<div class="source-note">
        <strong>Source choice:</strong> ${escapeHTML(sourceName || "selected rendering")}
        from ${escapeHTML(source.candidates?.length || 0)} candidate(s). ${escapeHTML(source.rationale || "")}
      </div>` : (paper.metadata?.source_file ? `<div class="source-note">
        <strong>Source read:</strong> ${escapeHTML(String(paper.metadata.source_file).split(/[\\/]/).pop())}
        ${paper.metadata.source_chars ? `· ${escapeHTML(formatNumber(paper.metadata.source_chars))} chars` : ""}
        ${paper.metadata.source_status ? `· ${escapeHTML(paper.metadata.source_status)}` : ""}
      </div>` : "");
      const indexLink = DATA.index_href
        ? `<p class="paper-subtitle"><a href="${escapeHTML(DATA.index_href)}">← all papers in this run</a></p>` : "";
      $("#paperHeader").innerHTML = `
        <div class="paper-kicker">${escapeHTML(paper.gene)} · ${escapeHTML(label)}</div>
        <h2>${escapeHTML(title)}</h2>
        ${indexLink}
        <p class="paper-subtitle">${escapeHTML(paper.metadata?.tool_rationale || "Follow the timeline below from the model request through the normalized curation decision.")}</p>
        ${sourceNote}
        <div class="chips">
          <span class="chip"><strong>${paper.call_count}</strong> model calls</span>
          <span class="chip"><strong>${paper.decision_count}</strong> decisions</span>
          <span class="chip accepted"><strong>${(paper.accepted_trace_ids || []).length}</strong> accepted</span>
          <span class="chip discarded"><strong>${(paper.discarded_trace_ids || []).length}</strong> discarded</span>
          <span class="chip"><strong>${formatNumber(paper.token_count)}</strong> tokens</span>
          <span class="chip"><strong>Model</strong> ${escapeHTML(models)}</span>
          <span class="chip"><strong>Stages</strong> ${escapeHTML((paper.stages || []).join(", ") || "unlabelled")}</span>
          ${variantChip}
        </div>`;
    }

    function renderIntegrity() {
      const integrity = DATA.integrity;
      const element = $("#integrity");
      const LEVEL_TITLE = {
        locked: "Locked and verified against write-time digests.",
        write_time_verified: "Verified against write-time digests.",
        generated_now: "Manifest generated now — not a tamper guarantee."
      };
      if (!integrity.valid) {
        element.className = "integrity bad";
        element.innerHTML = `<span>!</span><div><strong>Trace integrity needs attention.</strong>
          ${integrity.errors.map(escapeHTML).join("<br>")}</div>`;
        return;
      }
      element.className = integrity.level === "generated_now" ? "integrity partial" : "integrity good";
      element.innerHTML = `<span>${integrity.level === "generated_now" ? "◐" : "✓"}</span><div>
        <strong>${escapeHTML(LEVEL_TITLE[integrity.level] || integrity.level)}</strong>
        ${formatNumber(integrity.trace_count)} records; ${formatNumber(integrity.write_time_verified_records || 0)}
        match the SHA-256 digest recorded when they were written
        (${escapeHTML(integrity.index_name)} ${integrity.index_present ? "present" : "MISSING"}).
        ${integrity.manifest_generated_at ? `Manifest generated ${escapeHTML(formatTime(integrity.manifest_generated_at))}.` : ""}
        <div style="margin-top:5px">${escapeHTML(integrity.level_description || "")}</div></div>`;
    }

    function renderCoverage() {
      const coverage = DATA.coverage || {};
      const stages = coverage.stages || [];
      if (!stages.length) { $("#coverage").innerHTML = ""; return; }
      const missing = new Map((coverage.missing_decision_links || []).map(item => [item.stage, item]));
      $("#coverage").innerHTML = `<div class="panel-block">
        <h3>Route coverage</h3>
        <table class="coverage-table">
          <thead><tr><th>Stage</th><th>Calls</th><th>Decisions</th><th>Accepted links</th><th>Failures</th><th>Gap</th></tr></thead>
          <tbody>${stages.map(stage => `<tr>
            <td>${escapeHTML(stage.stage)}</td>
            <td>${escapeHTML(stage.llm_calls)}</td>
            <td>${escapeHTML(stage.decisions)}</td>
            <td>${escapeHTML(stage.accepted_links)}</td>
            <td>${escapeHTML(stage.failures)}</td>
            <td class="${missing.has(stage.stage) ? "gap" : ""}">${escapeHTML(missing.get(stage.stage)?.reason || "—")}</td>
          </tr>`).join("")}</tbody>
        </table></div>`;
    }

    function renderOmissions() {
      const omissions = DATA.omissions || [];
      const policy = DATA.size_policy || {};
      if (!omissions.length) {
        $("#omissions").innerHTML = `<div class="panel-block"><h3>Embedding policy</h3>
          Nothing was omitted. Bodies are embedded up to
          ${escapeHTML(formatNumber(policy.max_field_chars))} characters each.</div>`;
        return;
      }
      const rows = omissions.slice(0, 200).map(item => `<tr>
        <td>${escapeHTML(item.kind)}</td>
        <td>${escapeHTML(item.record || item.group || "")}</td>
        <td>${escapeHTML(item.characters_total != null
              ? `${formatNumber(item.characters_embedded)} of ${formatNumber(item.characters_total)} chars embedded`
              : `${formatNumber(item.records_embedded)} of ${formatNumber(item.records_total)} records embedded`)}</td>
      </tr>`).join("");
      $("#omissions").innerHTML = `<div class="panel-block"><h3>Omitted from this page (${escapeHTML(omissions.length)})</h3>
        <div>${escapeHTML(policy.note || "")} The full text is always on disk under the trace root
        <code>${escapeHTML(DATA.trace_root)}</code>.</div>
        <table class="coverage-table" style="margin-top:8px">
          <thead><tr><th>Kind</th><th>Record / group</th><th>Extent</th></tr></thead>
          <tbody>${rows}</tbody></table>
        ${omissions.length > 200 ? `<div style="margin-top:6px">…and ${escapeHTML(omissions.length - 200)} more (see the run log).</div>` : ""}
        </div>`;
    }

    function filteredRecords(paper) {
      return paper.records.filter(record => {
        if (state.filter === "failed" && record.response?.success !== false) return false;
        if (state.filter === "accepted" && record._is_accepted !== true) return false;
        if (state.filter === "retry" && !isRetry(record)) return false;
        if (["all","failed","accepted","retry"].indexOf(state.filter) === -1
            && record.record_type !== state.filter) return false;
        return !state.traceQuery || recordText(record).includes(state.traceQuery);
      });
    }

    function extractPrompt(payload) {
      if (!payload) return "(request payload unavailable)";
      if (typeof payload.input === "string") return payload.input;
      if (payload.input) return pretty(payload.input);
      if (payload.messages) {
        return payload.messages.map(message => {
          const role = String(message.role || "message").toUpperCase();
          return `${role}\n${typeof message.content === "string" ? message.content : pretty(message.content)}`;
        }).join("\n\n");
      }
      return pretty(payload);
    }

    function extractOutput(record) {
      if (record.response?.output_text) return pretty(record.response.output_text);
      const envelope = record.response?.envelope;
      if (envelope?.choices?.[0]?.message?.content) return pretty(envelope.choices[0].message.content);
      if (Array.isArray(envelope?.content)) {
        const texts = envelope.content.map(item => item?.text).filter(Boolean);
        if (texts.length) return texts.join("\n\n");
      }
      return envelope ? pretty(envelope) : "(response unavailable)";
    }

    function copyText(text, button) {
      const done = () => {
        const original = button.textContent;
        button.textContent = "Copied";
        setTimeout(() => button.textContent = original, 1200);
      };
      const copyWithTextarea = () => {
        const area = document.createElement("textarea");
        area.value = text;
        document.body.appendChild(area);
        area.select();
        try {
          if (document.execCommand("copy")) done();
          else button.textContent = "Copy unavailable";
        } catch (_) {
          button.textContent = "Copy unavailable";
        } finally {
          area.remove();
        }
      };
      if (navigator.clipboard?.writeText) {
        navigator.clipboard.writeText(text).then(done).catch(copyWithTextarea);
      } else copyWithTextarea();
    }

    function safeStorageGet(key) {
      try { return localStorage.getItem(key); } catch (_) { return null; }
    }

    function safeStorageSet(key, value) {
      try { localStorage.setItem(key, value); } catch (_) { /* optional */ }
    }

    function rationaleEntries(data) {
      const preferred = [
        "curation_rationale", "inclusion_rationale", "count_rationale",
        "tool_rationale", "rationale", "reasoning", "selection_policy",
        "fallback_reason", "override_reason", "decision_source"
      ];
      const entries = [];
      preferred.forEach(key => {
        const value = data?.[key];
        if (value != null && value !== "" && !entries.some(item => item[0] === key)) entries.push([key, value]);
      });
      return entries;
    }

    // A jump button is only rendered when the target record is actually embedded
    // in THIS page. A null trace id, or an id whose record was sharded to another
    // file or omitted by the size policy, must not produce a dead button. The
    // list is computed by the builder for exactly this payload.
    const embeddedTraceIds = new Set((DATA.embedded_trace_ids || []).map(String));
    const jumpable = (traceId) => Boolean(traceId) && embeddedTraceIds.has(String(traceId));
    const traceRef = (traceId, label) => {
      if (!traceId) return "";
      const text = escapeHTML(label || shortId(traceId));
      return jumpable(traceId)
        ? `<button class="jump-button" type="button" data-jump="${escapeHTML(traceId)}">${text}</button>`
        : `<span class="trace-id" title="not embedded in this page">${text} (not on this page)</span>`;
    };

    function renderFact(key, value) {
      if (key === "accepted_response_trace_id" && value) {
        return `${traceRef(value, "Open accepted model call ↑")}
          <div class="trace-id">${escapeHTML(shortId(value))}</div>`;
      }
      return escapeHTML(value);
    }

    function renderTraceIdList(key, values) {
      if (!Array.isArray(values) || !values.length) return "";
      const refs = values.filter(Boolean).map(value => traceRef(value)).join(" · ");
      if (!refs) return "";
      return `<dt>${escapeHTML(key.replaceAll("_"," "))}</dt><dd>${refs}</dd>`;
    }

    function decisionCard(record, index) {
      const data = record.event?.data || {};
      const rationale = rationaleEntries(data);
      const excluded = new Set(rationale.map(item => item[0]));
      const facts = Object.entries(data).filter(([key, value]) =>
        !excluded.has(key) && value != null && value !== "" && typeof value !== "object"
        && key !== "accepted_response_trace_ids"
      );
      const acceptedIds = Array.isArray(data.accepted_response_trace_ids)
        ? data.accepted_response_trace_ids.filter(Boolean) : [];
      const links =
        (acceptedIds.length > 1
          ? renderTraceIdList("all accepted trace ids", acceptedIds) : "")
        + renderTraceIdList("discarded trace ids", data.discarded_trace_ids)
        + renderTraceIdList("failed trace ids", data.failed_trace_ids)
        + (Array.isArray(data.attempt_trace_links) && data.attempt_trace_links.length
            ? `<dt>attempts</dt><dd>${data.attempt_trace_links.map(item => {
                const label = `#${item.attempt} ${item.role} → ${item.outcome}`;
                return item.trace_id
                  ? traceRef(item.trace_id, label)
                  : `<span class="trace-id">${escapeHTML(label)} (untraced)</span>`;
              }).join("<br>")}</dd>`
            : "");
      return `<article class="trace-card decision" id="trace-${escapeHTML(record.trace_id || index)}">
        ${cardHead(record, index, "Decision")}
        <div class="trace-body">
          <div class="decision-grid">
            <div class="panel">
              <div class="panel-title">Decision facts</div>
              <dl class="decision-facts">
                ${facts.length ? facts.map(([key, value]) => `<dt>${escapeHTML(key.replaceAll("_"," "))}</dt><dd>${renderFact(key, value)}</dd>`).join("") : "<dt>Event</dt><dd>Recorded without scalar fields</dd>"}
                ${links}
              </dl>
            </div>
            <div class="rationale">
              <strong>Why this outcome</strong>
              ${rationale.length ? rationale.map(([key, value]) => `<p><b>${escapeHTML(key.replaceAll("_"," "))}:</b> ${escapeHTML(typeof value === "string" ? value : pretty(value))}</p>`).join("") : `<p>No explicit rationale field was attached to this event. Inspect the raw event below.</p>`}
            </div>
          </div>
        </div>
        ${rawDetails(record)}
      </article>`;
    }

    function callCard(record, index) {
      const failed = record.response?.success === false;
      const request = record.request?.payload || {};
      const output = failed
        ? `${record.response?.error?.type || "Error"}: ${record.response?.error?.message || "Provider call failed"}`
        : extractOutput(record);
      const classes = ["trace-card"];
      if (failed) classes.push("failed");
      if (record._is_accepted) classes.push("accepted");
      if (record._is_discarded) classes.push("discarded");
      let kind = failed ? "Failed call" : "Model call";
      if (record._is_accepted) kind = "Accepted call";
      else if (record._is_discarded) kind = "Discarded call";
      else if (isRetry(record)) kind = "Retry / repair";
      return `<article class="${classes.join(" ")}" id="trace-${escapeHTML(record.trace_id || index)}">
        ${cardHead(record, index, kind)}
        <div class="trace-body">
          <div class="call-grid">
            <section class="panel">
              <div class="panel-title">Sent to model <button class="copy-button" type="button" data-copy="prompt-${index}">Copy</button></div>
              <pre id="prompt-${index}">${escapeHTML(extractPrompt(request))}</pre>
            </section>
            <section class="panel">
              <div class="panel-title">${failed ? "Provider error" : "Model returned"} <button class="copy-button" type="button" data-copy="output-${index}">Copy</button></div>
              <pre id="output-${index}">${escapeHTML(output)}</pre>
            </section>
          </div>
          ${reasoningNote(record)}
        </div>
        ${rawDetails(record)}
      </article>`;
    }

    function reasoningNote(record) {
      const capture = record.reasoning_capture;
      if (!capture) return "";
      const paths = capture.response_paths || [];
      const tokens = capture.reasoning_token_usage || {};
      const tokenText = Object.entries(tokens)
        .map(([key, value]) => `${key} = ${formatNumber(value)}`).join(", ");
      return `<div class="panel-block" style="margin:0">
        <h3>Provider reasoning</h3>
        ${paths.length
          ? `Reasoning CONTENT returned at: ${paths.map(escapeHTML).join(", ")}.`
          : "No reasoning content was returned by the provider."}
        ${tokenText ? `<div style="margin-top:4px">Reasoning token usage (billing telemetry, not exposed reasoning): ${escapeHTML(tokenText)}.</div>` : ""}
      </div>`;
    }

    function cardHead(record, index, kind) {
      const context = record.context || {};
      const model = record.request?.requested_model || context.model || "model unavailable";
      const provider = record.request?.provider || context.provider || "";
      const duration = record.response?.duration_seconds;
      const tokens = usageTotal(record);
      const meta = [
        model,
        provider,
        formatTime(recordTime(record)),
        Number.isFinite(duration) ? `${duration.toFixed(2)}s` : "",
        tokens ? `${formatNumber(tokens)} tokens` : "",
        context.attempt ? `attempt ${context.attempt}` : "",
        record._attempt_role && record._attempt_role !== record.record_type ? `role ${record._attempt_role}` : "",
        context.component || ""
      ].filter(Boolean).join(" · ");
      return `<header class="trace-head">
        <span class="step">${String(index + 1).padStart(2, "0")}</span>
        <div class="trace-heading">
          <h3>${escapeHTML(recordStage(record).replaceAll("_", " "))} <span class="chip">${escapeHTML(kind)}</span></h3>
          <div class="trace-meta">${escapeHTML(meta)}</div>
        </div>
        <span class="trace-id" title="${escapeHTML(record.trace_id || "")}">${escapeHTML(shortId(record.trace_id || ""))}</span>
      </header>`;
    }

    // Link the on-disk record only when the report and the trace root share a
    // filesystem-relative path (DATA.trace_href is set by the builder). Otherwise
    // show the absolute path as TEXT: a file:// href built from an unrelated
    // absolute path is both useless and a needless local-path disclosure.
    function recordRef(record) {
      const relative = record._report_path || "";
      if (!relative) return "(path unavailable)";
      if (!DATA.trace_href) return escapeHTML(relative);
      const href = DATA.trace_href + "/" + relative.split("/").map(encodeURIComponent).join("/");
      return `<a href="${escapeHTML(href)}">${escapeHTML(relative)}</a>`;
    }

    function rawDetails(record) {
      const truncated = record._truncated_fields
        ? `<div style="padding:8px 16px;color:var(--amber)">This card shows a bounded copy.
             The complete request/response is in the trace record:
             ${recordRef(record)}</div>` : "";
      return `<details class="raw"><summary>Exact trace record · ${escapeHTML(record._report_path || "")}</summary>
        ${truncated}<pre>${escapeHTML(pretty(record))}</pre></details>`;
    }

    function renderTimeline() {
      if (DATA.mode === "index") { $("#timeline").innerHTML = ""; return; }
      const paper = selectedPaper();
      if (!paper) {
        $("#timeline").innerHTML = `<div class="empty">No trace records were found.</div>`;
        return;
      }
      const records = filteredRecords(paper);
      $("#timeline").innerHTML = records.length
        ? records.map((record, index) =>
            record.record_type === "decision_event" ? decisionCard(record, index) : callCard(record, index)
          ).join("")
        : `<div class="empty">No events match the current timeline filter.</div>`;
      document.querySelectorAll("[data-copy]").forEach(button => {
        button.addEventListener("click", () => {
          const target = document.getElementById(button.dataset.copy);
          if (target) copyText(target.textContent, button);
        });
      });
      document.querySelectorAll("[data-jump]").forEach(button => {
        button.addEventListener("click", () => {
          state.filter = "all";
          document.querySelectorAll(".filter-button").forEach(item =>
            item.classList.toggle("active", item.dataset.filter === "all")
          );
          renderTimeline();
          requestAnimationFrame(() => {
            document.getElementById(`trace-${button.dataset.jump}`)?.scrollIntoView({
              behavior: "smooth",
              block: "center"
            });
          });
        });
      });
    }

    function toggleTheme() {
      const next = document.documentElement.dataset.theme === "dark" ? "light" : "dark";
      document.documentElement.dataset.theme = next;
      safeStorageSet("gvf-trace-theme", next);
    }

    window.addEventListener("hashchange", () => {
      const requested = decodeURIComponent(location.hash.replace(/^#/, ""));
      if (DATA.papers.some(paper => paper.id === requested)) {
        state.paperId = requested;
        render();
      }
    });
    init();
  </script>
</body>
</html>
"""
