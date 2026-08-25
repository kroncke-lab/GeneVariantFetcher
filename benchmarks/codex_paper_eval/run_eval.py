#!/usr/bin/env python3
"""Blinded, reproducible local-paper evaluation for Codex.

``prepare`` may use gold only to establish PMID eligibility and count-field
presence; it never exports gold values or row counts. ``score`` refuses to run
unless the prediction file's SHA-256 digest matches the immutable lock.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import math
import os
import random
import re
import statistics
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from pipeline.prompts import TABLE_ATTRIBUTION_GUIDANCE  # noqa: E402
from utils.gold_standard import (  # noqa: E402
    authoritative_gold_count,
    gold_row_excluded,
)
from utils.llm_trace import (  # noqa: E402
    TRACE_INDEX_NAME,
    TRACE_MANIFEST_NAME,
    build_trace_manifest,
    capture_llm_call,
    configure_llm_tracing,
    llm_trace_scope,
    record_trace_event,
    trace_lock_targets,
    validate_trace_manifest,
)
from utils.llm_trace_html import (  # noqa: E402
    TRACE_REPORT_NAME,
    build_trace_html_report,
)

# Genes eligible for seeded random sampling: the cardiac four with a manual,
# fully human-curated gold standard (gene_variant_fetcher_gold_standard/).
CARDIAC_GENES = ("SCN5A", "KCNH2", "KCNQ1", "RYR2")
# All genes a fixed manifest may reference. Genes outside CARDIAC_GENES score
# against the curated benchmark's adjudicated gold_overrides answer key (see
# gold_csv_path) — adjudicated, but not the manual gold standard — so report
# their results separately; never fold them into the cardiac headline.
GENES = CARDIAC_GENES + ("BRCA2",)
COUNT_FIELDS = ("carriers", "affected", "unaffected")
DEFAULT_CORPUS = REPO / "corpus"
DEFAULT_GOLD = REPO / "gene_variant_fetcher_gold_standard" / "normalized"
GOLD_OVERRIDES = REPO / "benchmarks" / "curated_extraction_eval" / "gold_overrides"
# Cheap proxy for "does this rendering still carry variant-level evidence". Only
# ever compared between candidate renderings of the same paper, so figure labels
# and other systematic noise land on both sides and cancel out.
VARIANT_TOKEN_RE = re.compile(
    r"\b[A-Z]\d{1,4}[A-Z]\b|\bp\.[A-Za-z]{3}\d{1,4}[A-Za-z]{3}\b"
)
AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text())


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def production_trace_lock_entries(
    predictions: dict,
) -> tuple[list[dict], list[Path], list[str]]:
    """Validate trace manifests referenced by an external production projection.

    Schema-1 projections cannot carry the harness-native per-paper trace refs,
    but they can still bind each complete ``gvf-run`` trace manifest into the
    prediction lock.  The manifest in turn hashes every call/decision record and
    is checked against the write-time trace index here and again at score time.
    """
    supplied = predictions.get("production_trace_manifests") or []
    errors: list[str] = []
    locked: list[dict] = []
    roots: list[Path] = []
    if not supplied:
        if predictions.get("strategy") == "production_gvf_run":
            errors.append(
                "production_gvf_run predictions require production_trace_manifests"
            )
        return locked, roots, errors
    if not isinstance(supplied, list):
        return locked, roots, ["production_trace_manifests must be a list"]
    seen_genes: set[str] = set()
    for index, entry in enumerate(supplied):
        label = f"production_trace_manifests[{index}]"
        if not isinstance(entry, dict):
            errors.append(f"{label} must be an object")
            continue
        gene = str(entry.get("gene") or "").upper()
        raw_path = str(entry.get("manifest") or "")
        expected_sha = str(entry.get("sha256") or "")
        if not gene or not raw_path or not expected_sha:
            errors.append(f"{label} requires gene, manifest, and sha256")
            continue
        if gene in seen_genes:
            errors.append(f"duplicate production trace manifest for {gene}")
            continue
        seen_genes.add(gene)
        path = Path(raw_path)
        if not path.is_absolute():
            path = REPO / path
        path = path.resolve()
        if path.name != TRACE_MANIFEST_NAME or not path.is_file():
            errors.append(f"{label}: trace manifest not found: {path}")
            continue
        actual_sha = digest(path)
        if actual_sha != expected_sha:
            errors.append(f"{label}: trace manifest digest mismatch")
            continue
        try:
            manifest = read_json(path)
            manifest_errors = validate_trace_manifest(path.parent, manifest)
        except Exception as exc:  # noqa: BLE001 - surface malformed audit artifacts
            errors.append(f"{label}: could not validate trace manifest: {exc}")
            continue
        errors.extend(f"{label}: {error}" for error in manifest_errors)
        for field in ("run_id", "llm_call_count", "decision_event_count"):
            if entry.get(field) != manifest.get(field):
                errors.append(f"{label}: {field} does not match trace manifest")
        integrity_level = (manifest.get("verification") or {}).get("level")
        if entry.get("integrity_level") != integrity_level:
            errors.append(f"{label}: integrity_level does not match trace manifest")
        if integrity_level != "write_time_verified":
            errors.append(f"{label}: trace is not write-time verified")
        if int(manifest.get("llm_call_count") or 0) <= 0:
            errors.append(f"{label}: trace contains no model calls")
        locked.append(
            {
                "gene": gene,
                "manifest": raw_path,
                "sha256": actual_sha,
                "run_id": manifest.get("run_id"),
                "llm_call_count": manifest.get("llm_call_count"),
                "decision_event_count": manifest.get("decision_event_count"),
                "integrity_level": integrity_level,
            }
        )
        roots.append(path.parent)
    prediction_genes = {
        str(p.get("gene") or "").upper() for p in predictions.get("papers", [])
    }
    if supplied and seen_genes != prediction_genes:
        errors.append(
            "production trace genes do not match prediction genes: "
            f"traces={sorted(seen_genes)} predictions={sorted(prediction_genes)}"
        )
    return locked, roots, errors


def is_table_line(line: str) -> bool:
    """Whether one line looks like a markdown or HTML table row.

    Deliberately narrow. A corpus census showed tab- and whitespace-aligned
    "tables" are almost entirely reference lists and PDF-justified prose, so
    widening this would inject noise rather than recover rows.
    """
    return line.count("|") >= 2 or bool(
        re.match(r"^\s*</?(?:table|tr|td|th)\b", line, flags=re.I)
    )


def table_row_count(text: str) -> int:
    return sum(1 for line in text.splitlines() if is_table_line(line))


def source_richness(path: Path) -> tuple[int, int, int]:
    """Extraction-relevant content signals for one candidate rendering."""
    text = path.read_text(errors="replace")
    return table_row_count(text), len(set(VARIANT_TOKEN_RE.findall(text))), len(text)


def choose_source(candidates: list[Path]) -> Path:
    """Pick the richest rendering, preserving candidate order on ties.

    Switching away from the first candidate requires Pareto dominance: at least
    as many table rows, distinct variant tokens, and characters, and strictly
    more of one. A rendering that trades prose for tables (or the reverse) leaves
    the priority order alone, so this can only add material, never drop any.
    """
    best = candidates[0]
    best_score = source_richness(best)
    for candidate in candidates[1:]:
        score = source_richness(candidate)
        if score != best_score and all(
            new >= old for new, old in zip(score, best_score)
        ):
            best, best_score = candidate, score
    return best


def usable_sources(
    corpus: Path,
    gene: str,
    minimum_chars: int,
    legacy_source_selection: bool = False,
    include_pmids: set[str] | None = None,
) -> list[dict]:
    papers = []
    for paper_dir in sorted(
        (corpus / gene).iterdir() if (corpus / gene).is_dir() else []
    ):
        if not paper_dir.is_dir() or not paper_dir.name.isdigit():
            continue
        pmid = paper_dir.name
        if include_pmids is not None and pmid not in include_pmids:
            continue
        candidates = [
            paper_dir / f"{pmid}_FULL_CONTEXT.md",
            paper_dir / f"{pmid}_CLEANED.md",
        ]
        usable = [
            p for p in candidates if p.is_file() and p.stat().st_size >= minimum_chars
        ]
        if not usable:
            continue
        # Legacy took the first candidate that cleared the size floor, so
        # _FULL_CONTEXT.md always won regardless of what it actually contained.
        source = usable[0] if legacy_source_selection else choose_source(usable)
        candidate_scores = []
        for candidate in usable:
            table_rows, variant_tokens, characters = source_richness(candidate)
            candidate_scores.append(
                {
                    "path": str(candidate.resolve()),
                    "sha256": digest(candidate),
                    "table_rows": table_rows,
                    "distinct_variant_tokens": variant_tokens,
                    "characters": characters,
                }
            )
        # choose_source() must return one of `usable`; a bare next() would raise
        # StopIteration (which propagates as a confusing RuntimeError inside a
        # generator) if it ever returned something outside the scored set.
        selected_path = str(source.resolve())
        selected_score = next(
            (item for item in candidate_scores if item["path"] == selected_path),
            None,
        )
        if selected_score is None:
            raise SystemExit(
                f"PMID {pmid}: chosen source {selected_path} is not among the "
                f"{len(candidate_scores)} scored candidates; selection.json would "
                "not describe the source actually read."
            )
        artifacts = paper_dir / f"{pmid}_artifacts.json"
        artifact_path = artifacts.resolve() if artifacts.is_file() else None
        pdf_paths = sorted(p.resolve() for p in paper_dir.rglob("*.pdf") if p.is_file())
        figure_paths = sorted(
            p.resolve()
            for p in paper_dir.rglob("*")
            if p.is_file()
            and p.suffix.lower() in {".png", ".jpg", ".jpeg", ".tif", ".tiff"}
        )
        pdfs = [str(path) for path in pdf_paths]
        figures = [str(path) for path in figure_paths]
        papers.append(
            {
                "gene": gene,
                "pmid": pmid,
                "source": str(source.resolve()),
                "source_sha256": digest(source),
                "source_bytes": source.stat().st_size,
                "source_selection": {
                    "policy": (
                        "first_usable_candidate"
                        if legacy_source_selection
                        else "pareto_richness"
                    ),
                    "selected": str(source.resolve()),
                    "selected_metrics": selected_score,
                    "candidates": candidate_scores,
                    "rationale": (
                        "Legacy ablation selected the first candidate clearing the "
                        "size floor."
                        if legacy_source_selection
                        else (
                            "Selected a candidate only when it Pareto-dominated the "
                            "priority candidate on table rows, distinct variant "
                            "tokens, and characters; preserved priority order when "
                            "candidate renderings traded off."
                        )
                    ),
                },
                "artifacts": str(artifact_path) if artifact_path else None,
                "artifacts_sha256": digest(artifact_path) if artifact_path else None,
                "pdfs": pdfs,
                "pdf_sha256": {str(path): digest(path) for path in pdf_paths},
                "figures": figures,
                "figure_sha256": {str(path): digest(path) for path in figure_paths},
            }
        )
    return papers


def material_digest_errors(paper: dict) -> list[str]:
    """Validate every local representation recorded during ``prepare``."""
    label = f"{paper.get('gene')}:{paper.get('pmid')}"
    errors: list[str] = []

    def check_one(kind: str, raw_path: str | None, expected: str | None) -> None:
        if not raw_path:
            if expected:
                errors.append(f"{label}: {kind} digest recorded without a path")
            return
        if not expected:
            errors.append(f"{label}: missing {kind} digest")
            return
        path = Path(raw_path)
        if not path.is_file():
            errors.append(f"{label}: {kind} file is missing: {path}")
        elif digest(path) != expected:
            errors.append(f"{label}: {kind} changed after selection: {path}")

    def check_many(kind: str, paths: list[str], recorded) -> None:
        if not isinstance(recorded, dict):
            errors.append(f"{label}: missing {kind} digest map")
            return
        path_set = set(paths)
        digest_set = set(recorded)
        if path_set != digest_set:
            errors.append(
                f"{label}: {kind} digest paths differ from selected paths: "
                f"missing={sorted(path_set - digest_set)} "
                f"extra={sorted(digest_set - path_set)}"
            )
        for raw_path in paths:
            check_one(kind, raw_path, recorded.get(raw_path))

    check_one("source", paper.get("source"), paper.get("source_sha256"))
    if paper.get("production_extraction_record") or paper.get(
        "production_extraction_record_sha256"
    ):
        check_one(
            "production extraction record",
            paper.get("production_extraction_record"),
            paper.get("production_extraction_record_sha256"),
        )
    if paper.get("representations") or paper.get("representation_sha256"):
        check_many(
            "text representation",
            list(paper.get("representations") or []),
            paper.get("representation_sha256"),
        )
    check_one(
        "artifact",
        paper.get("artifacts"),
        paper.get("artifacts_sha256"),
    )
    check_many("PDF", list(paper.get("pdfs") or []), paper.get("pdf_sha256"))
    check_many(
        "figure",
        list(paper.get("figures") or []),
        paper.get("figure_sha256"),
    )
    return errors


def selection_material_errors(selection: dict) -> list[str]:
    return [
        error
        for paper in selection.get("papers", [])
        for error in material_digest_errors(paper)
    ]


def gold_csv_path(gold_root: Path, gene: str) -> Path:
    """Resolve the gold CSV for one gene: the explicit root first, then overrides.

    The manual cardiac gold standard does not carry every gene a manifest may
    name. A gene absent from ``gold_root`` (BRCA2) resolves to the curated
    benchmark's curator-adjudicated ``gold_overrides`` answer key. Callers record
    the resolved path in run artifacts so a fallback is never silent.
    """
    primary = gold_root / f"{gene}_recall_input.csv"
    if primary.is_file():
        return primary
    fallback = GOLD_OVERRIDES / f"{gene}_recall_input.csv"
    if fallback.is_file():
        return fallback
    raise SystemExit(f"no gold CSV for {gene}: neither {primary} nor {fallback} exists")


def gold_count_eligible_pmids(gold_root: Path, gene: str) -> set[str]:
    """Return PMIDs with gold rows and at least one assertion for every count field.

    Only PMID membership and field presence are used during selection. Gold values
    and gold row counts are never written into the selection or extraction prompt.
    """
    path = gold_csv_path(gold_root, gene)
    coverage: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if gold_row_excluded(row):
                continue
            pmid = str(row.get("pmid", "")).strip()
            if not pmid or not str(row.get("variant", "")).strip():
                continue
            for field in COUNT_FIELDS:
                if authoritative_gold_count(row, field, parser=to_int) is not None:
                    coverage[pmid].add(field)
    return {pmid for pmid, fields in coverage.items() if fields == set(COUNT_FIELDS)}


def read_paper_manifest(path: Path) -> list[tuple[str, str]]:
    papers: list[tuple[str, str]] = []
    with path.open() as fh:
        for line_number, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t, ]+", line)
            if len(parts) != 2 or parts[0] not in GENES or not parts[1].isdigit():
                raise SystemExit(f"{path}:{line_number}: expected '<GENE> <PMID>'")
            papers.append((parts[0], parts[1]))
    if len(papers) != len(set(papers)):
        raise SystemExit(f"{path}: duplicate gene/PMID entries")
    return papers


def selection_metadata(selection: dict) -> dict:
    paper_count = len(selection.get("papers", []))
    manifest = selection.get("paper_manifest")
    if manifest:
        manifest_name = Path(manifest).name
        population = f"fixed manifest `{manifest_name}` ({paper_count} papers)"
        description = (
            f"Paper selection used the {population} from the downloaded-source, "
            "gold-count-eligible pool. Routing, extraction, counts, evidence, and "
            "source locations were gold-value-blind."
        )
        mode = "manifest"
    else:
        per_gene = selection.get("per_gene")
        seed = selection.get("seed")
        population = (
            f"seeded sample of {paper_count} papers "
            f"({per_gene} per cardiac gene; seed {seed})"
        )
        description = (
            f"Paper selection used a {population} from the downloaded-source, "
            "gold-count-eligible pool. Routing, extraction, counts, evidence, and "
            "source locations were gold-value-blind."
        )
        mode = "random"
    return {
        "mode": mode,
        "paper_manifest": manifest,
        "seed": selection.get("seed"),
        "per_gene": selection.get("per_gene"),
        "population": population,
        "description": description,
    }


def command_prepare(args) -> None:
    rng = random.Random(args.seed)
    run_dir = args.runs_dir / args.run_id
    if run_dir.exists():
        raise SystemExit(f"run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)
    selected: list[dict] = []
    eligible: dict[str, int] = {}
    pools: dict[str, dict[str, dict]] = {}
    # Manifest mode builds pools only for the genes the manifest names, so a
    # cardiac-only manifest never requires BRCA2 gold or corpus and vice versa.
    # Random mode samples the cardiac genes only: the extra genes have too few
    # gold papers to sample, and widening would change historical seeded draws.
    requested = (
        read_paper_manifest(args.paper_manifest) if args.paper_manifest else None
    )
    run_genes = (
        tuple(dict.fromkeys(gene for gene, _ in requested))
        if requested is not None
        else CARDIAC_GENES
    )
    for gene in run_genes:
        eligible_pmids = gold_count_eligible_pmids(args.gold_root, gene)
        source_pool = {
            paper["pmid"]: paper
            for paper in usable_sources(
                args.corpus_root,
                gene,
                args.minimum_chars,
                legacy_source_selection=args.legacy_source_selection,
                include_pmids=eligible_pmids,
            )
        }
        pools[gene] = {
            pmid: paper for pmid, paper in source_pool.items() if pmid in eligible_pmids
        }
        eligible[gene] = len(pools[gene])

    if requested is not None:
        for gene, pmid in requested:
            if pmid not in pools[gene]:
                raise SystemExit(
                    f"{gene} {pmid}: missing usable source or complete gold count coverage"
                )
            selected.append(pools[gene][pmid])
    else:
        for gene in run_genes:
            pool = list(pools[gene].values())
            if len(pool) < args.per_gene:
                raise SystemExit(
                    f"{gene}: only {len(pool)} gold-count-eligible usable papers, "
                    f"need {args.per_gene}"
                )
            selected.extend(rng.sample(pool, args.per_gene))
        rng.shuffle(selected)

    now = datetime.now(timezone.utc).isoformat()
    selection = {
        "schema_version": 1,
        "run_id": args.run_id,
        "seed": args.seed,
        "per_gene": args.per_gene,
        "minimum_source_chars": args.minimum_chars,
        "eligible_counts": eligible,
        # Which answer key each gene resolves to (manual cardiac gold vs the
        # adjudicated gold_overrides fallback) — recorded so scoring provenance
        # is explicit in the pre-gold lock.
        "gold_sources": {
            gene: str(gold_csv_path(args.gold_root, gene)) for gene in run_genes
        },
        "paper_manifest": str(args.paper_manifest.resolve())
        if args.paper_manifest
        else None,
        "prepared_at": now,
        "papers": selected,
        "blinding": (
            "Gold was used only to confirm PMID eligibility and the presence of "
            "carrier/affected/unaffected assertions. No gold values or row counts "
            "were written into this run or supplied to extraction."
        ),
    }
    predictions = {
        "schema_version": 2,
        "run_id": args.run_id,
        "started_at": now,
        "extraction_started_at": None,
        "completed_at": None,
        "extraction_elapsed_seconds": None,
        "token_usage": {
            "telemetry_available": False,
            "input_tokens": None,
            "output_tokens": None,
            "total_tokens": None,
            "estimate_method": None,
        },
        "papers": [
            {
                "gene": p["gene"],
                "pmid": p["pmid"],
                "tool": None,
                "tool_rationale": None,
                "elapsed_seconds": None,
                "source_completeness": None,
                "representations_available": None,
                "token_usage": None,
                "llm_trace_refs": [],
                "notes": None,
                "curation_rationale": None,
                "variants": [],
            }
            for p in selected
        ],
    }
    write_json(run_dir / "selection.json", selection)
    write_json(run_dir / "predictions.json", predictions)
    print(run_dir)


def validate_predictions(selection: dict, predictions: dict) -> list[str]:
    errors = []
    expected = {(p["gene"], p["pmid"]) for p in selection["papers"]}
    actual = {
        (p.get("gene"), str(p.get("pmid"))) for p in predictions.get("papers", [])
    }
    if expected != actual:
        errors.append(
            f"paper set mismatch: missing={sorted(expected - actual)} extra={sorted(actual - expected)}"
        )
    # Schema 1 is the external-import contract (e.g. a production gvf-run DB
    # projected by db_to_predictions.py): tool/rationale/variants are required,
    # but per-paper wall time and exact token telemetry don't exist there —
    # gvf-run does not aggregate them. Schema 2 is harness-native extraction,
    # where both are recorded per call and therefore mandatory.
    native = int(predictions.get("schema_version") or 1) >= 2
    for p in predictions.get("papers", []):
        label = f"{p.get('gene')}:{p.get('pmid')}"
        required = ("tool", "tool_rationale", "source_completeness") + (
            ("elapsed_seconds",) if native else ()
        )
        for key in required:
            if p.get(key) in (None, ""):
                errors.append(f"{label}: missing {key}")
        if (
            int(predictions.get("schema_version") or 1) >= 2
            and not str(p.get("curation_rationale") or "").strip()
        ):
            errors.append(f"{label}: missing curation_rationale")
        for i, row in enumerate(p.get("variants", [])):
            if not (row.get("variant") or "").strip():
                errors.append(f"{label} variant[{i}]: missing variant")
            if not (row.get("evidence") or "").strip():
                errors.append(f"{label} variant[{i}]: missing evidence")
            if not (row.get("source_location") or "").strip():
                errors.append(f"{label} variant[{i}]: missing source_location")
            if int(predictions.get("schema_version") or 1) >= 2:
                for rationale in ("inclusion_rationale", "count_rationale"):
                    if not str(row.get(rationale) or "").strip():
                        errors.append(f"{label} variant[{i}]: missing {rationale}")
            for field in COUNT_FIELDS:
                value = row.get(field)
                if value is not None and (
                    not isinstance(value, int) or isinstance(value, bool) or value < 0
                ):
                    errors.append(
                        f"{label} {row.get('variant')}:{field} must be nonnegative int or null"
                    )
        if p.get("tool") not in {"text", "table", "pdf", "ocr"}:
            errors.append(f"{label}: invalid tool {p.get('tool')!r}")
        usage = p.get("token_usage") or {}
        total_tokens = usage.get("total_tokens")
        if native and (
            not usage.get("telemetry_available")
            or not isinstance(total_tokens, int)
            or isinstance(total_tokens, bool)
            or total_tokens < 0
        ):
            errors.append(f"{label}: missing exact token telemetry")
        if int(predictions.get("schema_version") or 1) >= 2:
            refs = p.get("llm_trace_refs") or []
            stages = {
                (ref.get("context") or {}).get("stage")
                for ref in refs
                if isinstance(ref, dict)
            }
            required_stages = {
                "representation_route",
                "representation_route_decision",
                "paper_curation",
                "paper_curation_decision",
            }
            missing_stages = sorted(required_stages - stages)
            if missing_stages:
                errors.append(
                    f"{label}: missing locked LLM trace stages {missing_stages}"
                )
    return errors


def command_lock(args) -> None:
    run_dir = args.run_dir
    lock_path = run_dir / "LOCK.json"
    if lock_path.exists():
        raise SystemExit(f"already locked: {lock_path}")
    selection_path = run_dir / "selection.json"
    prediction_path = run_dir / "predictions.json"
    selection, predictions = read_json(selection_path), read_json(prediction_path)
    trace_root = run_dir / "llm_traces"
    trace_manifest_path = trace_root / TRACE_MANIFEST_NAME
    trace_report_path = run_dir / TRACE_REPORT_NAME
    trace_manifest = None
    trace_index_path = trace_root / TRACE_INDEX_NAME
    production_trace_locks, production_trace_roots, production_trace_errors = (
        production_trace_lock_entries(predictions)
    )
    if int(predictions.get("schema_version") or 1) >= 2:
        # Rebuilding at lock time is deliberate — extraction may have appended
        # after its own manifest — but it now cross-checks each record against
        # the WRITE-TIME digest in trace_index.jsonl, so a record forged between
        # extraction and lock surfaces as an error instead of being re-blessed.
        trace_manifest = build_trace_manifest(
            trace_root,
            output_path=trace_manifest_path,
            run_id=selection.get("run_id"),
        )
        build_trace_html_report(
            trace_root,
            output_path=trace_report_path,
            run_dir=run_dir,
            title=f"{selection.get('run_id') or 'Paper evaluation'} · LLM trace review",
            run_id=selection.get("run_id"),
            locked=True,
        )
    errors = validate_predictions(selection, predictions)
    errors.extend(selection_material_errors(selection))
    errors.extend(production_trace_errors)
    if trace_manifest is not None:
        errors.extend(validate_trace_manifest(trace_root, trace_manifest))
    if errors:
        raise SystemExit("prediction validation failed:\n- " + "\n- ".join(errors))
    prelock_gold_usage = predictions.get("prelock_gold_usage") or {}
    if "read_only_layer_scoring_possible" in prelock_gold_usage:
        production_projection = bool(
            prelock_gold_usage["read_only_layer_scoring_possible"]
        )
    else:
        # Legacy production projections did not record whether gvf-run was
        # allowed to auto-discover gold. Treat only that missing-provenance
        # case conservatively; an explicit false value is affirmative evidence.
        production_projection = int(predictions.get("schema_version") or 1) < 2
    lock_statement = (
        "Prediction content was finalized before external scoring. This external "
        "production projection may come from a workflow that read registered gold "
        "for read-only layer scorecards; those scores do not feed back into "
        "extraction."
        if production_projection
        else (
            "Predictions finalized before gold values or gold row counts were "
            "exposed to extraction; score is the first phase that reads those values."
        )
    )
    lock = {
        "locked_at": datetime.now(timezone.utc).isoformat(),
        "selection_sha256": digest(selection_path),
        "predictions_sha256": digest(prediction_path),
        "llm_trace_manifest_sha256": (
            digest(trace_manifest_path) if trace_manifest is not None else None
        ),
        "llm_trace_report_sha256": (
            digest(trace_report_path) if trace_manifest is not None else None
        ),
        # trace_index.jsonl holds the WRITE-TIME digest of every record. Leaving
        # it out of the lock left the one artifact that can prove tampering
        # unprotected, so "unchanged since lock" did not cover it.
        "llm_trace_index_sha256": (
            digest(trace_index_path)
            if trace_manifest is not None and trace_index_path.is_file()
            else None
        ),
        "llm_trace_integrity_level": (
            trace_manifest["verification"]["level"]
            if trace_manifest is not None
            else None
        ),
        "production_trace_manifests": production_trace_locks,
        "statement": lock_statement,
    }
    write_json(lock_path, lock)
    prediction_path.chmod(0o444)
    if trace_manifest is not None:
        trace_manifest_path.chmod(0o444)
        trace_report_path.chmod(0o444)
        # trace_lock_targets() includes trace_index.jsonl, which the prior
        # chmod set omitted.
        for target in trace_lock_targets(trace_root):
            target.chmod(0o444)
    for production_trace_root in production_trace_roots:
        for target in trace_lock_targets(production_trace_root):
            target.chmod(0o444)
    print(lock_path)


ROUTE_INSTRUCTIONS = """You are routing one biomedical paper for blinded curation.
Use only the representation previews supplied below. Do not use databases, prior
extractions, benchmarks, or gold standards.

Target gene: {gene}. PMID: {pmid}.

Choose exactly one authoritative representation:
- text: running full text is the clearest source.
- table: structured table rows carry the variant-level person counts. Note that a
  table may be a compilation that cites other studies rather than a report of this
  study's own observations; rows alone do not make it first-party data.
- pdf: PDF-layout text is more complete or preserves a table the markdown loses.
- ocr: figure/pedigree images are necessary because textual representations omit
  the genotype/phenotype evidence.

Prefer the representation that best preserves variant-level carrier, affected,
and unaffected evidence. Return JSON only:
{{
  "tool": "text|table|pdf|ocr",
  "tool_rationale": "...",
  "source_completeness": "full_text|partial_text|abstract_only"
}}

REPRESENTATION CATALOG:
{catalog}
"""


EXTRACTION_INSTRUCTIONS = (
    """You are independently curating one biomedical paper.
You are BLINDED: do not use databases, prior extraction results, benchmark files,
or any gold standard. Use only the supplied local paper material.

Target gene: {gene}. PMID: {pmid}.
Authoritative representation selected: {tool}.

Extract every human {gene} variant for which this paper supplies carrier/patient
phenotype evidence. For each variant:
- carriers: distinct genotype-positive people, including probands and genotyped
  relatives; exclude controls/population allele counts and non-carrier relatives.
- affected: carrier people with the relevant cardiac phenotype in this paper.
- unaffected: carrier people explicitly asymptomatic/clinically normal.
Use null when a count cannot be determined, and 0 only when zero is explicit; always
emit the variant even when all three counts are null. Do not turn cohort size,
alleles, families, events, or functional experiments into person counts. Do not
double-count the same people across prose and tables.

"""
    + TABLE_ATTRIBUTION_GUIDANCE
    + """

The selected representation is supplied below. Cite concise evidence and a
specific section/table/page/figure location for every row. Return JSON only:
{{
  "notes": "...",
  "curation_rationale": "Why evidence was included or excluded, including ambiguities.",
  "variants": [{{
    "variant": "...", "carriers": null, "affected": null, "unaffected": null,
    "evidence": "...", "source_location": "...",
    "inclusion_rationale": "Why this is in-scope human variant evidence from this paper.",
    "count_rationale": "How each integer or null count was derived without double counting."
  }}]
}}

LOCAL MATERIAL:
{material}
"""
)


def parse_json_response(text: str) -> dict:
    value = text.strip()
    if value.startswith("```"):
        value = value.split("\n", 1)[1].rsplit("```", 1)[0]
    return json.loads(value)


def looks_truncated_json(text: str) -> bool:
    """Whether ``text`` began as a JSON object but was cut off mid-emission.

    Azure reports ``status="completed"`` even when a reasoning model runs out of
    output budget mid-string, so the response envelope cannot be trusted for this.
    Garbage that never started as JSON is a real failure and must still raise.
    """
    value = text.strip()
    if value.startswith("```"):
        value = value.split("\n", 1)[-1]
    value = value.strip()
    return value.startswith("{") and not value.endswith("}")


def targeted_preview(text: str, gene: str, max_chars: int) -> str:
    if not text:
        return ""
    pattern = re.compile(
        rf"\b{re.escape(gene)}\b|"
        r"\b(?:p\.)?[A-Z][a-z]{0,2}\d{1,5}(?:[A-Z][a-z]{0,2}|\*|X|Ter)\b|"
        r"\b(?:carrier|affected|unaffected|asymptomatic|proband|pedigree|variant)\b",
        flags=re.I,
    )
    pieces = [text[: min(1800, max_chars)]]
    for match in pattern.finditer(text):
        start, end = max(0, match.start() - 350), min(len(text), match.end() + 650)
        pieces.append(text[start:end])
        if sum(len(piece) for piece in pieces) >= max_chars:
            break
    return "\n\n[...]\n\n".join(pieces)[:max_chars]


def markdown_table_material(text: str, max_chars: int) -> str:
    blocks: list[str] = []
    current: list[str] = []
    for line in text.splitlines():
        if is_table_line(line):
            current.append(line)
        elif current:
            if len(current) >= 2:
                blocks.append("\n".join(current))
            current = []
    if len(current) >= 2:
        blocks.append("\n".join(current))
    return "\n\n".join(blocks)[:max_chars]


def extract_pdf_text(pdf_paths: list[str], max_chars: int) -> str:
    parts: list[str] = []
    for raw_path in pdf_paths:
        path = Path(raw_path)
        try:
            completed = subprocess.run(
                ["pdftotext", "-layout", str(path), "-"],
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except (FileNotFoundError, OSError, subprocess.SubprocessError):
            continue
        text = completed.stdout.strip()
        if text:
            parts.append(f"### PDF {path.name}\n\n{text}")
        if sum(len(part) for part in parts) >= max_chars:
            break
    return "\n\n".join(parts)[:max_chars]


def image_data_url(path: Path) -> str:
    mime = {
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".tif": "image/tiff",
        ".tiff": "image/tiff",
    }.get(path.suffix.lower(), "image/png")
    return f"data:{mime};base64,{base64.b64encode(path.read_bytes()).decode()}"


def reasoning_params(model: str, effort: str) -> dict:
    """Responses-API reasoning kwargs the deployment will actually accept.

    Delegates the model allow-list to ``utils.llm_utils`` so the benchmark and
    production agree about which deployments receive an effort setting. They did
    not before: this function sent ``{"reasoning": {"effort": ...}}`` to
    everything except grok, while production gated on
    ``gpt-5|gpt5|o1|o3|o4-mini`` — so for any other deployment the benchmark
    result did not describe the production configuration for the same model
    string. Grok-4-class deployments reason by default and reject
    ``reasoning.effort`` with a 400, so they still get no reasoning block.
    """
    from utils.llm_utils import build_responses_reasoning_param

    if "grok" in model.lower():
        return {}
    return build_responses_reasoning_param(model, effort)


def effective_effort(model: str, effort: str) -> str:
    """The effort label to record, so reports never claim an unsent setting."""
    return effort if reasoning_params(model, effort) else "model_default"


def supports_images(model: str) -> bool:
    """Whether the deployment accepts ``input_image`` content.

    grok-4.3 on Azure rejects image input outright, so offering it an ocr route
    would fail the paper instead of falling back to a textual representation. The
    parsed artifact still carries figure captions on every route.
    """
    return "grok" not in model.lower()


def response_usage(response) -> dict[str, int]:
    usage = getattr(response, "usage", None)
    return {
        "input_tokens": int(getattr(usage, "input_tokens", 0) or 0),
        "output_tokens": int(getattr(usage, "output_tokens", 0) or 0),
        "total_tokens": int(getattr(usage, "total_tokens", 0) or 0),
    }


def add_usage(
    predictions: dict, target: dict, response, model: str, effort: str
) -> None:
    increment = response_usage(response)
    for container in (predictions, target):
        current = container.get("token_usage") or {}
        input_tokens = int(current.get("input_tokens") or 0) + increment["input_tokens"]
        output_tokens = (
            int(current.get("output_tokens") or 0) + increment["output_tokens"]
        )
        container["token_usage"] = {
            "telemetry_available": True,
            "input_tokens": input_tokens,
            "output_tokens": output_tokens,
            "total_tokens": input_tokens + output_tokens,
            "estimate_method": None,
            "model": model,
            "reasoning_efforts": sorted(
                set((current.get("reasoning_efforts") or []) + [effort])
            ),
        }


def add_trace_ref(target: dict, summary: dict[str, Any] | None) -> None:
    if summary is not None:
        target.setdefault("llm_trace_refs", []).append(summary)


def command_extract(args) -> None:
    """Route and run isolated per-paper reads through the Azure OpenAI endpoint."""
    from openai import OpenAI

    run_dir = args.run_dir
    if (run_dir / "LOCK.json").exists():
        raise SystemExit("refusing to extract: run is already locked")
    selection = read_json(run_dir / "selection.json")
    configure_llm_tracing(
        run_dir / "llm_traces",
        run_id=selection.get("run_id"),
        enabled=True,
    )
    material_errors = selection_material_errors(selection)
    if material_errors:
        raise SystemExit(
            "selected material validation failed:\n- " + "\n- ".join(material_errors)
        )
    predictions_path = run_dir / "predictions.json"
    predictions = read_json(predictions_path)
    by_key = {(p["gene"], str(p["pmid"])): p for p in predictions["papers"]}
    base = (os.environ.get("AZURE_AI_API_BASE") or "").strip().rstrip("/")
    key = (os.environ.get("AZURE_AI_API_KEY") or "").strip()
    if not base or not key:
        raise SystemExit("AZURE_AI_API_BASE and AZURE_AI_API_KEY are required")
    client = OpenAI(base_url=base, api_key=key, timeout=args.timeout)
    command_started = time.monotonic()
    predictions["extraction_started_at"] = datetime.now(timezone.utc).isoformat()
    write_json(predictions_path, predictions)
    total_papers = len(selection["papers"])
    for index, paper in enumerate(selection["papers"], 1):
        target = by_key[(paper["gene"], paper["pmid"])]
        if target.get("tool") and not args.force:
            print(
                f"[{index}/{total_papers}] skip {paper['gene']} {paper['pmid']} "
                "(already complete)",
                flush=True,
            )
            continue
        source_text = Path(paper["source"]).read_text(errors="replace")
        artifact_text = ""
        if paper.get("artifacts"):
            artifact_text = Path(paper["artifacts"]).read_text(errors="replace")[
                : args.max_artifact_chars
            ]
        table_text = markdown_table_material(source_text, args.max_source_chars)
        if args.legacy_table_material and artifact_text:
            # Legacy folded the artifact into the table string, which is what made
            # "table" look available on papers with no table rows at all.
            table_text = (
                f"### Parsed artifact data\n\n{artifact_text}\n\n"
                f"### Tables preserved in text\n\n{table_text}"
            )
        pdf_text = extract_pdf_text(paper.get("pdfs") or [], args.max_source_chars)
        figures = [Path(path) for path in paper.get("figures") or []]
        representations = {
            "text": bool(source_text.strip()),
            # Real table rows only. The parsed artifact carries table and figure
            # captions, not table bodies, so folding it in here advertised a table
            # route for papers whose table rows never survived acquisition — and
            # the route then delivered captions plus a keyword preview instead of
            # the running text the model could otherwise have read.
            "table": bool(table_text.strip()),
            "pdf": bool(pdf_text.strip()),
            "ocr": bool(figures) and supports_images(args.model),
        }
        target["representations_available"] = [
            name for name, available in representations.items() if available
        ]
        # Distinguishes "this paper has no figures" from "this model cannot read
        # them", so a text-only route is not mistaken for a routing decision.
        target["figures_withheld_no_vision"] = bool(figures) and not supports_images(
            args.model
        )
        catalog = (
            f"Available={target['representations_available']}\n\n"
            f"## TEXT PREVIEW\n{targeted_preview(source_text, paper['gene'], args.route_preview_chars)}\n\n"
            f"## TABLE PREVIEW\n{targeted_preview(table_text, paper['gene'], args.route_preview_chars)}\n\n"
            f"## PDF PREVIEW\n{targeted_preview(pdf_text, paper['gene'], args.route_preview_chars)}\n\n"
            f"## OCR INVENTORY\n"
            + "\n".join(path.name for path in figures[: args.max_ocr_images])
            + (
                # Attached to whichever route is chosen, so it is context for judging
                # source_completeness, not a fifth option. Captions naming tables that
                # have no rows above are the signal for partial_text.
                f"\n\n## PARSED ARTIFACT (attached to every route)\n"
                f"{targeted_preview(artifact_text, paper['gene'], args.route_preview_chars)}"
                if artifact_text and not args.legacy_table_material
                else ""
            )
        )
        route_prompt = ROUTE_INSTRUCTIONS.format(
            gene=paper["gene"], pmid=paper["pmid"], catalog=catalog
        )
        started = time.monotonic()
        print(
            f"[{index}/{total_papers}] route {paper['gene']} {paper['pmid']}",
            flush=True,
        )
        route_request = {
            "model": args.model,
            "input": route_prompt,
            "max_output_tokens": args.route_max_output_tokens,
            **reasoning_params(args.model, args.route_reasoning_effort),
        }
        with llm_trace_scope(
            gene=paper["gene"],
            pmid=paper["pmid"],
            stage="representation_route",
            component="codex_paper_eval",
            attempt=1,
        ):
            route_response, route_trace = capture_llm_call(
                provider="azure_openai_responses",
                requested_model=args.model,
                resolved_model=args.model,
                request=route_request,
                call=lambda: client.responses.create(**route_request),
            )
        add_trace_ref(target, route_trace)
        add_usage(
            predictions,
            target,
            route_response,
            args.model,
            effective_effort(args.model, args.route_reasoning_effort),
        )
        write_json(predictions_path, predictions)
        route = parse_json_response(route_response.output_text)
        requested_tool = str(route.get("tool", "")).lower()
        tool = requested_tool
        route_fallback = None
        if tool not in representations or not representations[tool]:
            available_fallback = next(
                name
                for name in ("text", "table", "pdf", "ocr")
                if representations[name]
            )
            route["tool_rationale"] = (
                f"{route.get('tool_rationale', '')} Requested {tool or 'no tool'}, "
                f"which was unavailable; used {available_fallback}."
            ).strip()
            route_fallback = {
                "requested_tool": requested_tool or None,
                "selected_tool": available_fallback,
                "reason": "requested representation was unavailable",
            }
            tool = available_fallback
        with llm_trace_scope(
            gene=paper["gene"],
            pmid=paper["pmid"],
            stage="representation_route_decision",
            component="codex_paper_eval",
        ):
            route_decision_trace = record_trace_event(
                "representation_route_decision",
                {
                    "representations_available": target["representations_available"],
                    "model_requested_tool": requested_tool or None,
                    "selected_tool": tool,
                    "fallback": route_fallback,
                    "tool_rationale": route.get("tool_rationale"),
                    "source_completeness": route.get("source_completeness"),
                    "route_call_trace_id": (
                        route_trace.get("trace_id") if route_trace else None
                    ),
                },
            )
        add_trace_ref(target, route_decision_trace)

        if tool == "table":
            # The table route carries the full running text, not a keyword preview.
            # Deciding whether a row is this study's own observation needs the
            # caption, footnotes, symbol definitions, and the prose that introduces
            # the table; a 6k excerpt strips exactly that, which is how a compilation
            # table gets read as first-party data.
            supporting = (
                source_text[: args.max_source_chars]
                if not args.legacy_table_material
                else targeted_preview(
                    source_text, paper["gene"], args.route_preview_chars
                )
            )
            material = (
                f"### Structured table rows\n\n{table_text[: args.max_source_chars]}"
                f"\n\n### Full running text (for table scope and provenance)\n\n"
                f"{supporting}"
            )
        elif tool == "pdf":
            material = f"### PDF layout text\n\n{pdf_text[: args.max_source_chars]}"
        elif tool == "ocr":
            material = "### Supporting targeted text\n\n" + targeted_preview(
                source_text, paper["gene"], args.route_preview_chars
            )
        else:
            material = f"### Full/partial running text\n\n{source_text[: args.max_source_chars]}"
        if artifact_text and not args.legacy_table_material:
            # Captions and figure/supplement metadata are useful on every route, and
            # the table route no longer has to spend its source budget carrying them.
            material = f"### Parsed artifact data\n\n{artifact_text}\n\n{material}"

        prompt = EXTRACTION_INSTRUCTIONS.format(
            gene=paper["gene"], pmid=paper["pmid"], tool=tool, material=material
        )
        response_input: str | list[dict] = prompt
        if tool == "ocr":
            content: list[dict] = [{"type": "input_text", "text": prompt}]
            for image_path in figures[: args.max_ocr_images]:
                try:
                    content.append(
                        {"type": "input_image", "image_url": image_data_url(image_path)}
                    )
                except OSError:
                    continue
            response_input = [{"role": "user", "content": content}]

        print(
            f"[{index}/{total_papers}] read {paper['gene']} {paper['pmid']} via {tool}",
            flush=True,
        )
        budget = args.max_output_tokens
        for attempt in (1, 2):
            extraction_request = {
                "model": args.model,
                "input": response_input,
                "max_output_tokens": budget,
                **reasoning_params(args.model, args.reasoning_effort),
            }
            with llm_trace_scope(
                gene=paper["gene"],
                pmid=paper["pmid"],
                stage="paper_curation",
                component="codex_paper_eval",
                attempt=attempt,
                selected_representation=tool,
            ):
                extraction_response, extraction_trace = capture_llm_call(
                    provider="azure_openai_responses",
                    requested_model=args.model,
                    resolved_model=args.model,
                    request=extraction_request,
                    call=lambda: client.responses.create(**extraction_request),
                )
            add_trace_ref(target, extraction_trace)
            add_usage(
                predictions,
                target,
                extraction_response,
                args.model,
                effective_effort(args.model, args.reasoning_effort),
            )
            write_json(predictions_path, predictions)
            try:
                result = parse_json_response(extraction_response.output_text)
                break
            except json.JSONDecodeError:
                # Reasoning models spend hidden tokens against this budget before
                # emitting JSON, so a long variant table can truncate mid-string.
                # One paid retry beats aborting a 48-paper run on one paper.
                truncated = getattr(
                    extraction_response, "status", None
                ) == "incomplete" or looks_truncated_json(
                    extraction_response.output_text
                )
                if attempt == 2 or not truncated:
                    raise
                budget = min(budget * 4, args.max_output_tokens_ceiling)
                target["output_budget_retry"] = budget
                write_json(predictions_path, predictions)
                print(
                    f"[{index}/{total_papers}] retry {paper['gene']} {paper['pmid']} "
                    f"— output truncated, raising budget to {budget}",
                    flush=True,
                )
        elapsed = time.monotonic() - started
        # Trace references are recorder-owned metadata. Do not let an unexpected
        # model-emitted key replace the audited call history.
        trace_refs = list(target.get("llm_trace_refs") or [])
        target.update(result)
        target["llm_trace_refs"] = trace_refs
        target["tool"] = tool
        target["tool_rationale"] = route.get("tool_rationale")
        target["source_completeness"] = route.get("source_completeness")
        target["elapsed_seconds"] = round(elapsed, 3)
        with llm_trace_scope(
            gene=paper["gene"],
            pmid=paper["pmid"],
            stage="paper_curation_decision",
            component="codex_paper_eval",
        ):
            curation_decision_trace = record_trace_event(
                "paper_curation_decision",
                {
                    "selected_model": args.model,
                    "selected_representation": tool,
                    "variant_count": len(target.get("variants") or []),
                    "curation_rationale": target.get("curation_rationale")
                    or target.get("notes"),
                    "accepted_response_trace_id": (
                        extraction_trace.get("trace_id") if extraction_trace else None
                    ),
                    "selection_policy": (
                        "The parsed response from the final successful extraction "
                        "attempt became the locked paper prediction. A retry occurs "
                        "only after a truncated JSON emission."
                    ),
                },
            )
        add_trace_ref(target, curation_decision_trace)
        write_json(predictions_path, predictions)
        print(
            f"[{index}/{total_papers}] done {paper['gene']} {paper['pmid']} "
            f"variants={len(target.get('variants', []))} seconds={elapsed:.1f}",
            flush=True,
        )
    predictions["completed_at"] = datetime.now(timezone.utc).isoformat()
    predictions["extraction_elapsed_seconds"] = round(
        time.monotonic() - command_started, 3
    )
    # Token usage is checkpointed after every response so interrupted runs retain it.
    write_json(predictions_path, predictions)
    trace_root = run_dir / "llm_traces"
    build_trace_manifest(
        trace_root,
        output_path=trace_root / TRACE_MANIFEST_NAME,
        run_id=selection.get("run_id"),
    )
    build_trace_html_report(
        trace_root,
        output_path=run_dir / TRACE_REPORT_NAME,
        run_dir=run_dir,
        title=f"{selection.get('run_id') or 'Paper evaluation'} · LLM trace review",
    )


def to_int(value):
    if value is None or str(value).strip() == "":
        return None
    try:
        return int(float(str(value)))
    except ValueError:
        return None


def load_gold(gold_root: Path, gene: str, pmid: str) -> list[dict]:
    path = gold_csv_path(gold_root, gene)
    with path.open(newline="") as fh:
        rows = []
        for row in csv.DictReader(fh):
            if str(row.get("pmid", "")).strip() != pmid:
                continue
            if gold_row_excluded(row):
                continue
            rows.append(
                {
                    "variant": str(row.get("variant", "")).strip(),
                    **{
                        field: authoritative_gold_count(row, field, parser=to_int)
                        for field in COUNT_FIELDS
                    },
                }
            )
        return rows


_VARIANT_CANDIDATE_CACHE: dict[tuple[str, str], list[str]] = {}


def variant_candidates(value: str, gene: str) -> list[str]:
    """Return embedded notations without changing the locked prediction."""
    cache_key = (value, gene)
    cached = _VARIANT_CANDIDATE_CACHE.get(cache_key)
    if cached is not None:
        return cached
    candidates = [value]
    patterns = (
        r"\bp\.\(?[A-Z][a-z]{2}\d+(?:[A-Z][a-z]{2}|Ter|fs[^\s,;)]*)\)?",
        (
            r"\bc\.[0-9*?+-]+(?:_[0-9*?+-]+)?"
            r"(?:delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+|[ACGT]>[ACGT])"
        ),
        # A trailing "*" is not a \b boundary, so keep the bare-star
        # alternative unanchored on the right.
        r"\b[A-Z]\d{1,5}(?:[A-Z]\b|X\b|\*)",
    )
    for pattern in patterns:
        candidates.extend(re.findall(pattern, value, flags=re.I))
    candidates.extend(
        re.findall(
            r"(?<![A-Za-z0-9])(?:[A-Z])?\d{1,5}"
            r"(?:_(?:[A-Z])?\d{1,5})?(?:del|ins)[A-Z]*",
            value,
            flags=re.I,
        )
    )
    # Reference-less protein-range duplications occur in older mutation tables
    # (``p.360_361dupKQ``). Keep this spelling narrow: a comma-form token such
    # as ``p.QKQR360,361dup`` is not the same grammar and must not become a
    # source-grounding alias for linked database rows.
    candidates.extend(
        re.findall(
            r"(?<![A-Za-z0-9])p\.(\d{1,5}_\d{1,5}dup[A-Z]+)\b",
            value,
            flags=re.I,
        )
    )
    candidates.extend(
        f"c.{token}"
        for token in re.findall(
            (
                r"(?<![A-Za-z0-9._])([0-9*?+-]+(?:_[0-9*?+-]+)?"
                r"(?:delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+))"
            ),
            value,
            flags=re.I,
        )
    )
    unprefixed_cdna = re.fullmatch(
        (
            r"\s*([0-9*?+-]+(?:_[0-9*?+-]+)?"
            r"(?:delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+|[ACGT]>[ACGT]))\s*"
        ),
        value,
        flags=re.I,
    )
    if unprefixed_cdna:
        candidates.append(f"c.{unprefixed_cdna.group(1)}")
    for original, position, replacement in re.findall(
        r"p\.\(?([A-Z][a-z]{2})(\d{1,5})([A-Z][a-z]{2}|Ter|X|\*)\)?",
        value,
    ):
        left = AA3_TO_1.get(original)
        right = "X" if replacement in {"Ter", "X", "*"} else AA3_TO_1.get(replacement)
        if left and right:
            candidates.append(f"{left}{position}{right}")
    for original, position, replacement in re.findall(
        r"\(([A-Z][a-z]{2})(\d{1,5})([A-Z][a-z]{2}|Ter|X|\*)\)",
        value,
    ):
        left = AA3_TO_1.get(original)
        right = "X" if replacement in {"Ter", "X", "*"} else AA3_TO_1.get(replacement)
        if left and right:
            candidates.append(f"{left}{position}{right}")
    for original, position in re.findall(
        r"p\.\(?([A-Z][a-z]{2})(\d{1,5})del\)?",
        value,
    ):
        if left := AA3_TO_1.get(original):
            candidates.append(f"{left}{position}del")
    for left_aa, left_pos, right_aa, right_pos, operation, inserted in re.findall(
        (
            r"p\.\(?([A-Z][a-z]{2})(\d{1,5})_([A-Z][a-z]{2})(\d{1,5})"
            r"(del|ins)([A-Z][a-z]{2})?\)?"
        ),
        value,
    ):
        left, right = AA3_TO_1.get(left_aa), AA3_TO_1.get(right_aa)
        inserted_aa = AA3_TO_1.get(inserted) if inserted else ""
        if left and right:
            candidates.append(
                f"{left}{left_pos}_{right}{right_pos}{operation}{inserted_aa}"
            )
    for original, position, _new_aa in re.findall(
        (
            r"p\.\(?([A-Z][a-z]{2})(\d{1,5})([A-Z][a-z]{2})"
            r"(?:fs(?:Ter|X)?\d*|X\d+|fs\d+X)"
        ),
        value,
    ):
        if left := AA3_TO_1.get(original):
            candidates.append(f"{left}{position}fsX")
    for residue, position in re.findall(
        r"\b([A-Z])(\d{1,5})fs(?:[/+*]?\d+|X|Ter\d*)?\b", value, flags=re.I
    ):
        candidates.append(f"{residue.upper()}{position}fsX")
    for residue, position in re.findall(
        r"\b([A-Z])(\d{1,5})[A-Z]fs(?:Ter|X)?\d*\b", value, flags=re.I
    ):
        candidates.append(f"{residue.upper()}{position}fsX")
    # Legacy compact frameshifts that carry the stop distance inline:
    # R192CFS91X, P400FS/62X.
    for residue, position in re.findall(
        r"\b([A-Z])(\d{1,5})[A-Z]?fs[/+*]?\d*(?:X|Ter)\b", value, flags=re.I
    ):
        candidates.append(f"{residue.upper()}{position}fsX")
    # 3-letter frameshift with no new-residue letter: p.Arg192fs.
    for original, position in re.findall(r"p\.\(?([A-Z][a-z]{2})(\d{1,5})fs", value):
        if left := AA3_TO_1.get(original):
            candidates.append(f"{left}{position}fsX")
    # 1-letter stop with a bare '*' (p.Q530*): '*' is not a \b boundary, so
    # the generic pattern above misses it; emit the legacy X form directly.
    for residue, position in re.findall(r"\b([A-Z])(\d{1,5})\*", value):
        candidates.append(f"{residue}{position}X")
    for residue, position in re.findall(
        r"\bfs([A-Z])(\d{1,5})(?:[/+*]?\d+|aa)*\b", value, flags=re.I
    ):
        candidates.append(f"{residue.upper()}{position}fsX")
    for residue, position in re.findall(
        r"\bfs([A-Z][a-z]{2})(\d{1,5})(?:[/+*]?\d+|aa)*\b", value
    ):
        if left := AA3_TO_1.get(residue):
            candidates.append(f"{left}{position}fsX")
    # Curated splice shorthands tolerate a separator and a synonymous residue
    # letter before the marker: A344SP, A344/SP, A344A/SPLICE.
    for residue, position in re.findall(
        r"\b([A-Z])(\d{1,5})(?:[A-Z]?/)?(?:sp|splice)\b", value, flags=re.I
    ):
        candidates.append(f"{residue.upper()}{position}SP")
    for position, residue in re.findall(r"\bdel(\d{1,5})([A-Z])\b", value, flags=re.I):
        candidates.append(f"{residue.upper()}{position}del")
    for residue, position in re.findall(r"\bdel([A-Z])(\d{1,5})\b", value, flags=re.I):
        candidates.append(f"{residue.upper()}{position}del")
    from utils.structural_alleles import expand_structural_keys

    candidates.extend(expand_structural_keys(value, gene))
    result = list(dict.fromkeys(c.strip() for c in candidates))
    _VARIANT_CANDIDATE_CACHE[cache_key] = result
    return result


def cdna_indel_position(value: str) -> int | None:
    match = re.search(
        r"(?:^|[^A-Za-z])c?\.?(\d+)(?:_\d+)?(?:del|dup|ins)",
        value,
        flags=re.I,
    )
    return int(match.group(1)) // 3 + 1 if match else None


def cdna_splice_codon(value: str) -> int | None:
    """Codon adjacent to an intronic-offset substitution (splice-site cDNA).

    Curators encode splice-site variants as a terminal protein event at the
    flanking codon (M159X for c.477+1G>A). Only intronic-offset substitutions
    qualify — exonic substitutions must match through translation, never
    through a codon bridge.
    """
    match = re.fullmatch(
        r"(?:C\.)?(\d+)[+-]\d+[ACGT]>[ACGT]",
        re.sub(r"\s+", "", value.upper()),
    )
    return (int(match.group(1)) + 2) // 3 if match else None


def splice_event_position(value: str) -> int | None:
    compact = re.sub(r"^(?:P\.)|[\s/]", "", value.upper())
    match = re.fullmatch(r"[A-Z](\d+)[A-Z]?(?:SP|SPLICE)", compact)
    return int(match.group(1)) if match else None


_CODONS_BY_AA: dict[str, tuple[str, ...]] = {}
for _codon, _aa in {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "TAA": "X",
    "TAG": "X",
    "TGA": "X",
}.items():
    _CODONS_BY_AA.setdefault(_aa, ())
    _CODONS_BY_AA[_aa] = (*_CODONS_BY_AA[_aa], _codon)


def cdna_protein_translation_consistent(cdna: str, protein: str) -> bool:
    """True when an exonic cDNA substitution and a protein substitution are
    the same allele: same codon index, and some codon of the reference
    residue turns into a codon of the alternate residue under exactly that
    nucleotide change at that in-codon offset. Never matches on codon index
    alone (c.1826A>G can be p.D609G, never p.D609N)."""
    cdna_match = re.fullmatch(
        r"(?:C\.)?(\d+)([ACGT])>([ACGT])",
        re.sub(r"\s+", "", cdna.upper()),
    )
    protein_match = re.fullmatch(
        r"(?:P\.)?([A-Z])(\d{1,5})([A-Z]|\*)",
        re.sub(r"\s+", "", protein.upper()),
    )
    if not cdna_match or not protein_match:
        return False
    base, ref_nt, alt_nt = (
        int(cdna_match.group(1)),
        cdna_match.group(2),
        cdna_match.group(3),
    )
    ref_aa, position, alt_aa = (
        protein_match.group(1),
        int(protein_match.group(2)),
        "X" if protein_match.group(3) == "*" else protein_match.group(3),
    )
    if (base + 2) // 3 != position:
        return False
    offset = (base - 1) % 3
    alt_codons = set(_CODONS_BY_AA.get(alt_aa, ()))
    return any(
        codon[offset] == ref_nt
        and codon[:offset] + alt_nt + codon[offset + 1 :] in alt_codons
        for codon in _CODONS_BY_AA.get(ref_aa, ())
    )


def protein_event_position(value: str) -> int | None:
    compact = re.sub(r"^(?:P\.)|\s+", "", value.upper())
    match = re.fullmatch(r"[A-Z](\d+)(?:FS(?:X|TER)?\d*|X|\*|DEL)", compact)
    if match:
        return int(match.group(1))
    stop = re.fullmatch(r"(\d+)STOP", compact)
    return int(stop.group(1)) if stop else None


def terminal_event_position(value: str) -> int | None:
    compact = re.sub(r"^(?:P\.)|\s+", "", value.upper())
    match = re.fullmatch(r"(?:[A-Z](\d+)(?:X|\*)|(\d+)STOP)", compact)
    if not match:
        return None
    return int(match.group(1) or match.group(2))


def structural_event_key(value: str) -> tuple[int, int | None, str, str] | None:
    compact = re.sub(r"^(?:P\.)|[\s()]", "", value.upper())
    match = re.search(
        r"(?:[A-Z])?(\d+)(?:_(?:[A-Z])?(\d+))?(DEL|INS)([A-Z]*)",
        compact,
    )
    if not match:
        return None
    start, end, operation, inserted = match.groups()
    return int(start), int(end) if end else None, operation, inserted


_CDNA_DELETION_RE = re.compile(r"^C?\.?(\d+)(?:_(\d+))?DEL([ACGT]{2,})$")


def cdna_deletion_span(value: str) -> tuple[int, int | None, str] | None:
    """Parse a multi-base cDNA deletion into ``(start, end, sequence)``.

    ``c.692_693delCA`` yields ``(692, 693, "CA")``; the legacy single-coordinate
    spelling ``c.693delCA`` yields ``(693, None, "CA")``. Single-base deletions
    are excluded: ``c.123delA`` is already unambiguous, and bridging it would
    let a one-base event match a range.
    """
    compact = re.sub(r"\s+", "", value.upper())
    match = _CDNA_DELETION_RE.fullmatch(compact)
    if not match:
        return None
    start, end, sequence = match.groups()
    return int(start), int(end) if end else None, sequence


def cdna_deletion_endpoint_match(left: str, right: str) -> bool:
    """True when a one-coordinate deletion names an endpoint of the same span.

    Papers write the same two-base event as ``c.693delCA`` or as
    ``c.692_693delCA``; curators only ever use the span form (0 of 6971 gold
    rows use the single-coordinate spelling). Requiring the deleted bases to be
    identical and the lone coordinate to *be* one of the two endpoints keeps
    this from inventing a position: it never bridges two open spellings, never
    bridges different sequences, and never has to guess which end was written.

    Span width is deliberately not required to equal the sequence length --
    gold itself carries ``c.4066_4068delTT`` and ``c.5355_5354delCT``, so a
    width check would reject curated rows.
    """
    first, second = cdna_deletion_span(left), cdna_deletion_span(right)
    if not first or not second:
        return False
    left_start, left_end, left_seq = first
    right_start, right_end, right_seq = second
    if left_seq != right_seq:
        return False
    if (left_end is None) == (right_end is None):
        return False
    if left_end is None:
        return left_start in (right_start, right_end)
    return right_start in (left_start, left_end)


_ARROW_RE = re.compile(r"\s*(?:-+>|→)\s*")


def matches(a: str, b: str, gene: str) -> bool:
    from utils.variant_normalizer import normalize_variant, variants_match

    if not a or not b:
        return False
    from utils.structural_alleles import structural_keys_match

    if structural_keys_match(a, b, gene):
        return True
    # Legacy arrow spellings (2592+1G->A, G→A) mean '>'.
    a = _ARROW_RE.sub(">", a)
    b = _ARROW_RE.sub(">", b)
    compact_left = re.sub(r"^(?:C\.)|\s+", "", a.upper())
    compact_right = re.sub(r"^(?:C\.)|\s+", "", b.upper())
    cdna_pattern = re.compile(r"^\d+(?:[+-]\d+)?[ACGT]*>[ACGT]+$")
    if (
        compact_left == compact_right
        and cdna_pattern.match(compact_left)
        and cdna_pattern.match(compact_right)
    ):
        return True
    for left in variant_candidates(a, gene):
        for right in variant_candidates(b, gene):
            left_key = re.sub(r"^(?:P\.)|\s+", "", left.upper())
            right_key = re.sub(r"^(?:P\.)|\s+", "", right.upper())
            if left_key == right_key:
                return True
            left_cdna, right_cdna = (
                cdna_indel_position(left),
                cdna_indel_position(right),
            )
            left_protein = protein_event_position(left)
            right_protein = protein_event_position(right)
            if left_cdna is not None and left_cdna == right_protein:
                return True
            if right_cdna is not None and right_cdna == left_protein:
                return True
            left_terminal = terminal_event_position(left)
            right_terminal = terminal_event_position(right)
            if left_terminal is not None and left_terminal == right_terminal:
                return True
            # Splice bridge: an intronic-offset cDNA substitution matches a
            # terminal/splice protein event at the flanking codon (M159X vs
            # c.477+1G>A). Requires exactly one side to be intronic cDNA and
            # the other a protein X/fs/del/SP event — missense never bridges.
            left_splice = cdna_splice_codon(left)
            right_splice = cdna_splice_codon(right)
            if left_splice is not None and right_splice is None:
                target = protein_event_position(right)
                if target is None:
                    target = splice_event_position(right)
                if target is not None and left_splice == target:
                    return True
            if right_splice is not None and left_splice is None:
                target = protein_event_position(left)
                if target is None:
                    target = splice_event_position(left)
                if target is not None and right_splice == target:
                    return True
            # Exonic substitution bridge: nucleotide-verified translation
            # only (c.1127G>A <-> R376H); codon index alone never matches.
            if cdna_protein_translation_consistent(left, right):
                return True
            if cdna_protein_translation_consistent(right, left):
                return True
            if cdna_deletion_endpoint_match(left, right):
                return True
            left_structural = structural_event_key(left)
            right_structural = structural_event_key(right)
            if left_structural and right_structural:
                left_start, left_end, left_op, left_inserted = left_structural
                right_start, right_end, right_op, right_inserted = right_structural
                # A one-residue deletion and a range deletion are distinct events.
                # Legacy insertions may omit the second flanking residue, but
                # deletion endpoints must agree exactly.
                ends_compatible = left_end == right_end or (
                    left_op == right_op == "INS"
                    and (left_end is None or right_end is None)
                )
                if (
                    left_start == right_start
                    and ends_compatible
                    and left_op == right_op
                    and left_inserted == right_inserted
                ):
                    return True
            try:
                if normalize_variant(left, gene) == normalize_variant(right, gene):
                    return True
                if variants_match(left, right, gene):
                    return True
            except Exception:
                continue
    return False


def twin_identical(a: str, b: str, gene: str) -> bool:
    """Strict equivalent-allele identity for prediction de-duplication.

    Deliberately NARROWER than matches(): the scored matcher carries
    position-level bridges (splice codon, terminal position, cDNA-indel
    codon) that are correct scoring fuzz against curated gold shorthand but
    are NOT allele identity — using them as identity fuses distinct alleles
    (c.477+1G>A and c.477+2T>C both splice-bridge to M159X). Identity here
    means: identical normalized notation, identical structural coordinates,
    normalizer equality, a nucleotide-verified cDNA<->protein translation, or
    a multi-base deletion whose one-coordinate spelling names a span endpoint.
    """
    from utils.variant_normalizer import normalize_variant, variants_match

    if not a or not b:
        return False
    from utils.structural_alleles import structural_keys_match

    if structural_keys_match(a, b, gene):
        return True
    a = _ARROW_RE.sub(">", a)
    b = _ARROW_RE.sub(">", b)
    for left in variant_candidates(a, gene):
        for right in variant_candidates(b, gene):
            left_key = re.sub(r"^(?:P\.)|\s+", "", left.upper())
            right_key = re.sub(r"^(?:P\.)|\s+", "", right.upper())
            if left_key == right_key:
                return True
            left_structural = structural_event_key(left)
            right_structural = structural_event_key(right)
            if (
                left_structural
                and right_structural
                and left_structural == right_structural
            ):
                return True
            if cdna_protein_translation_consistent(left, right):
                return True
            if cdna_protein_translation_consistent(right, left):
                return True
            # Same deleted bases, and the one-coordinate spelling names an
            # endpoint of the span: allele identity, not a position bridge.
            # ``c.693delCA`` could in principle also mean ``c.693_694delCA``,
            # so this is only safe because merge_notation_twins refuses a row
            # that is identical to more than one kept row.
            if cdna_deletion_endpoint_match(left, right):
                return True
            try:
                if normalize_variant(left, gene) == normalize_variant(right, gene):
                    return True
                if variants_match(left, right, gene):
                    return True
            except Exception:
                continue
    return False


def merge_notation_twins(rows: list[dict], gene: str) -> tuple[list[dict], int]:
    """Collapse same-paper prediction rows that are the same variant in a
    different notation (vertical tables emit the cDNA and protein halves as
    two rows; the matcher scores one TP and strands the twin as an FP).

    Hard bounds (grok-4.6 final review): identity is twin_identical(), never
    the scored matches() with its position-level bridges; refuse when the row
    matches more than one kept row (ambiguity); count fields transfer ONLY
    when the kept row is all-null across COUNT_FIELDS (the classic vertical
    split: identity-only protein row + count-bearing cDNA row). Complementary
    unions of two partially-counted rows are refused — they can invent a
    count vector no single source row asserted.
    """
    merged: list[dict] = []
    twins = 0
    for row in rows:
        homes = [
            m
            for m in merged
            if twin_identical(row.get("variant"), m.get("variant"), gene)
        ]
        conflict = any(
            m.get(field) is not None
            and row.get(field) is not None
            and m.get(field) != row.get(field)
            for m in homes
            for field in COUNT_FIELDS
        )
        if len(homes) != 1 or conflict:
            merged.append(dict(row))
            continue
        home = homes[0]
        twins += 1
        if all(home.get(field) is None for field in COUNT_FIELDS):
            for field in COUNT_FIELDS:
                if row.get(field) is not None:
                    home[field] = row[field]
    return merged, twins


def score_one(gene: str, pmid: str, predicted: dict, gold: list[dict]) -> dict:
    pred, merged_twin_rows = merge_notation_twins(predicted.get("variants", []), gene)
    used = set()
    pairs, fps = [], []
    for p in pred:
        hit = next(
            (
                i
                for i, g in enumerate(gold)
                if i not in used and matches(p["variant"], g["variant"], gene)
            ),
            None,
        )
        if hit is None:
            fps.append(p)
        else:
            used.add(hit)
            pairs.append((p, gold[hit]))
    misses = [g for i, g in enumerate(gold) if i not in used]
    tp, fp, fn = len(pairs), len(fps), len(misses)
    precision = tp / (tp + fp) if tp + fp else (1.0 if not gold else 0.0)
    recall = tp / (tp + fn) if tp + fn else None
    f1 = (
        2 * precision * recall / (precision + recall)
        if recall is not None and precision + recall
        else 0.0
    )
    count = {}
    disagreements = []
    count_bearing_tp = sum(
        any(p.get(field) is not None for field in COUNT_FIELDS) for p, _g in pairs
    )
    counted_extra_rows = sum(
        any(p.get(field) is not None for field in COUNT_FIELDS) for p in fps
    )
    precision_vs_counted_gold_pmids = (
        tp / (tp + counted_extra_rows) if tp + counted_extra_rows else None
    )
    count_bearing_precision = (
        count_bearing_tp / (count_bearing_tp + counted_extra_rows)
        if count_bearing_tp + counted_extra_rows
        else None
    )
    for field in COUNT_FIELDS:
        gold_asserted = sum(1 for row in gold if row.get(field) is not None)
        gold_asserted_zero = sum(1 for row in gold if row.get(field) == 0)
        observed = [
            (p.get(field), g.get(field), p["variant"])
            for p, g in pairs
            if g.get(field) is not None and p.get(field) is not None
        ]
        errors = [pv - gv for pv, gv, _ in observed]
        # Gold encodes "no <field> reported" as an explicit 0 while the
        # pipeline deliberately abstains with NULL, so pooled count recall
        # mixes a convention gap with real attribution misses. Stratify.
        observed_zero = sum(1 for _pv, gv, _v in observed if gv == 0)
        gold_asserted_nonzero = gold_asserted - gold_asserted_zero
        count[field] = {
            "gold_asserted": gold_asserted,
            "predicted": len(observed),
            "recall": len(observed) / gold_asserted if gold_asserted else None,
            "mae": statistics.fmean(abs(e) for e in errors) if errors else None,
            "rmse": math.sqrt(statistics.fmean(e * e for e in errors))
            if errors
            else None,
            "gold_asserted_zero": gold_asserted_zero,
            "gold_asserted_nonzero": gold_asserted_nonzero,
            "predicted_on_zero_gold": observed_zero,
            "predicted_on_nonzero_gold": len(observed) - observed_zero,
            "recall_zero_gold": (
                observed_zero / gold_asserted_zero if gold_asserted_zero else None
            ),
            "recall_nonzero_gold": (
                (len(observed) - observed_zero) / gold_asserted_nonzero
                if gold_asserted_nonzero
                else None
            ),
        }
        disagreements.extend(
            {
                "variant": v,
                "field": field,
                "predicted": pv,
                "gold": gv,
                "error": pv - gv,
            }
            for pv, gv, v in observed
            if pv != gv
        )
    return {
        "gene": gene,
        "pmid": pmid,
        "tool": predicted.get("tool"),
        "tool_rationale": predicted.get("tool_rationale"),
        "source_completeness": predicted.get("source_completeness"),
        "elapsed_seconds": predicted.get("elapsed_seconds"),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "merged_notation_twin_rows": merged_twin_rows,
        "token_usage": predicted.get("token_usage"),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "counted_precision": {
            # This is the repository's established
            # precision_vs_counted_gold_pmids definition: every matched gold row
            # is signal; only extra predictions that assert a patient count enter
            # the false-positive denominator.
            "matched_rows": tp,
            "counted_extra_rows": counted_extra_rows,
            "precision_vs_counted_gold_pmids": precision_vs_counted_gold_pmids,
            # Keep the stricter, differently-denominated diagnostic explicit so
            # it cannot be compared accidentally with the metric above.
            "count_bearing_matched_rows": count_bearing_tp,
            "precision_among_count_bearing_predictions": count_bearing_precision,
        },
        "count": count,
        "matched_variants": [
            {"predicted": p["variant"], "gold": g["variant"]} for p, g in pairs
        ],
        "missed_gold": [g["variant"] for g in misses],
        "extra_predictions": [p["variant"] for p in fps],
        "count_errors": disagreements,
    }


def aggregate(scores: list[dict]) -> dict:
    tp, fp, fn = (sum(s[k] for s in scores) for k in ("tp", "fp", "fn"))
    precision = tp / (tp + fp) if tp + fp else (1.0 if not fn else 0.0)
    recall = tp / (tp + fn) if tp + fn else None
    result = {
        "papers": len(scores),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "merged_notation_twin_rows": sum(
            int(s.get("merged_notation_twin_rows") or 0) for s in scores
        ),
        "precision": precision,
        "recall": recall,
        "f1": 2 * precision * recall / (precision + recall)
        if recall is not None and precision + recall
        else 0.0,
        "elapsed_seconds": sum(float(s.get("elapsed_seconds") or 0) for s in scores),
        "token_usage": {
            field: sum(
                int((s.get("token_usage") or {}).get(field) or 0) for s in scores
            )
            for field in ("input_tokens", "output_tokens", "total_tokens")
        },
        "count": {},
    }
    counted_extra_rows = sum(
        int((s.get("counted_precision") or {}).get("counted_extra_rows") or 0)
        for s in scores
    )
    count_bearing_matched_rows = sum(
        int((s.get("counted_precision") or {}).get("count_bearing_matched_rows") or 0)
        for s in scores
    )
    result["counted_precision"] = {
        "matched_rows": tp,
        "counted_extra_rows": counted_extra_rows,
        "precision_vs_counted_gold_pmids": (
            tp / (tp + counted_extra_rows) if tp + counted_extra_rows else None
        ),
        "count_bearing_matched_rows": count_bearing_matched_rows,
        "precision_among_count_bearing_predictions": (
            count_bearing_matched_rows
            / (count_bearing_matched_rows + counted_extra_rows)
            if count_bearing_matched_rows + counted_extra_rows
            else None
        ),
    }
    for field in COUNT_FIELDS:
        asserted = sum(s["count"][field]["gold_asserted"] for s in scores)
        observed = sum(s["count"][field]["predicted"] for s in scores)
        # Include zero-error observations, which are not in count_errors.
        abs_sum = sum(
            (s["count"][field]["mae"] or 0) * s["count"][field]["predicted"]
            for s in scores
        )
        sq_sum = sum(
            ((s["count"][field]["rmse"] or 0) ** 2) * s["count"][field]["predicted"]
            for s in scores
        )
        asserted_zero = sum(
            int(s["count"][field].get("gold_asserted_zero") or 0) for s in scores
        )
        observed_zero = sum(
            int(s["count"][field].get("predicted_on_zero_gold") or 0) for s in scores
        )
        asserted_nonzero = asserted - asserted_zero
        observed_nonzero = observed - observed_zero
        result["count"][field] = {
            "gold_asserted": asserted,
            "predicted": observed,
            "recall": observed / asserted if asserted else None,
            "mae": abs_sum / observed if observed else None,
            "rmse": math.sqrt(sq_sum / observed) if observed else None,
            "gold_asserted_zero": asserted_zero,
            "gold_asserted_nonzero": asserted_nonzero,
            "predicted_on_zero_gold": observed_zero,
            "predicted_on_nonzero_gold": observed_nonzero,
            "recall_zero_gold": (
                observed_zero / asserted_zero if asserted_zero else None
            ),
            "recall_nonzero_gold": (
                observed_nonzero / asserted_nonzero if asserted_nonzero else None
            ),
        }
    return result


def format_rate(value) -> str:
    return "n/a" if value is None else f"{100 * value:.1f}%"


def format_number(value, digits: int = 3) -> str:
    return "n/a" if value is None else f"{value:.{digits}f}"


def write_evidence_csv(predictions: dict, path: Path) -> None:
    columns = [
        "gene",
        "pmid",
        "tool",
        "source_completeness",
        "variant",
        *COUNT_FIELDS,
        "evidence",
        "source_location",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for paper in predictions["papers"]:
            for variant in paper.get("variants", []):
                writer.writerow(
                    {
                        "gene": paper["gene"],
                        "pmid": paper["pmid"],
                        "tool": paper.get("tool"),
                        "source_completeness": paper.get("source_completeness"),
                        **{column: variant.get(column) for column in columns[4:]},
                    }
                )


def write_paper_metrics_csv(scores: list[dict], path: Path) -> None:
    columns = [
        "gene",
        "pmid",
        "tool",
        "source_completeness",
        "tp",
        "fp",
        "fn",
        "precision",
        "recall",
        "f1",
        "merged_notation_twin_rows",
        "matched_rows",
        "counted_extra_rows",
        "precision_vs_counted_gold_pmids",
        "count_bearing_matched_rows",
        "precision_among_count_bearing_predictions",
        "elapsed_seconds",
        "input_tokens",
        "output_tokens",
        "total_tokens",
    ]
    for field in COUNT_FIELDS:
        columns.extend(
            [
                f"{field}_gold_asserted",
                f"{field}_predicted",
                f"{field}_recall",
                f"{field}_mae",
                f"{field}_rmse",
                f"{field}_gold_asserted_zero",
                f"{field}_gold_asserted_nonzero",
                f"{field}_predicted_on_zero_gold",
                f"{field}_predicted_on_nonzero_gold",
                f"{field}_recall_zero_gold",
                f"{field}_recall_nonzero_gold",
            ]
        )
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for score in scores:
            usage = score.get("token_usage") or {}
            row = {
                **{key: score.get(key) for key in columns[:10]},
                **(score.get("counted_precision") or {}),
                "merged_notation_twin_rows": score.get("merged_notation_twin_rows"),
                "elapsed_seconds": score.get("elapsed_seconds"),
                "input_tokens": usage.get("input_tokens"),
                "output_tokens": usage.get("output_tokens"),
                "total_tokens": usage.get("total_tokens"),
            }
            for field in COUNT_FIELDS:
                for metric, value in score["count"][field].items():
                    row[f"{field}_{metric}"] = value
            writer.writerow(row)


def matcher_adjudication_basis(predicted: str, gold: str) -> str:
    pair = f"{predicted} {gold}"
    if "(" in predicted and ")" in predicted:
        return "embedded parenthetical protein notation"
    if re.search(r"\bc\.", predicted, flags=re.I) and re.search(
        r"(?:fs|del|X|\*)", gold, flags=re.I
    ):
        return "cDNA indel mapped to its protein codon event"
    if re.search(r"(?:splice|\bsp\b)", pair, flags=re.I):
        return "legacy splice-label equivalence"
    if re.search(r"(?:del|ins)", pair, flags=re.I):
        return "legacy deletion/insertion shorthand"
    if re.search(r"(?:stop|Ter|\*)", pair, flags=re.I):
        return "legacy stop-label equivalence"
    return "embedded or normalized equivalent notation"


def write_matcher_adjudication_csv(
    scores: list[dict], raw_report: dict, path: Path
) -> int:
    raw_by_paper = {
        (paper["gene"], str(paper["pmid"])): paper
        for paper in raw_report.get("papers", [])
    }
    rows = []
    for score in scores:
        raw = raw_by_paper.get((score["gene"], str(score["pmid"])), {})
        raw_extras = Counter(raw.get("extra_predictions", []))
        raw_misses = Counter(raw.get("missed_gold", []))
        for pair in score["matched_variants"]:
            predicted, gold = pair["predicted"], pair["gold"]
            if raw_extras[predicted] and raw_misses[gold]:
                rows.append(
                    {
                        "gene": score["gene"],
                        "pmid": score["pmid"],
                        "predicted": predicted,
                        "gold": gold,
                        "basis": matcher_adjudication_basis(predicted, gold),
                    }
                )
                raw_extras[predicted] -= 1
                raw_misses[gold] -= 1
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["gene", "pmid", "predicted", "gold", "basis"]
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def write_markdown_report(report: dict, path: Path) -> None:
    overall = report["overall"]
    counted_precision = overall.get("counted_precision") or {}
    token_usage = report["token_usage"] or {}
    telemetry_available = bool(report.get("papers")) and all(
        bool((score.get("token_usage") or {}).get("telemetry_available"))
        for score in report["papers"]
    )
    integrity = report.get("integrity") or {}
    native_traces_locked = bool(integrity.get("llm_trace_manifest_sha256"))
    production_traces_locked = bool(integrity.get("production_trace_manifests"))
    timing = report.get("timing") or {}
    wall_timing_available = bool(timing.get("started_at")) and bool(
        timing.get("completed_at")
    )
    papers_per_gene = {
        gene: metric["papers"] for gene, metric in report["by_gene"].items()
    }
    balanced_papers = len(set(papers_per_gene.values())) == 1
    per_gene_label = (
        f"{next(iter(papers_per_gene.values()))} per cardiac gene"
        if balanced_papers
        else f"per-gene counts {papers_per_gene}"
    )
    selection_info = report.get("selection") or {}
    selection_description = selection_info.get(
        "description",
        (
            "Paper selection used the recorded evaluation set. Routing, extraction, "
            "counts, evidence, and source locations were gold-value-blind."
        ),
    )
    selection_population = selection_info.get(
        "population",
        f"recorded evaluation set ({overall['papers']} papers)",
    )
    prelock_gold_usage = report.get("prelock_gold_usage") or {}
    if prelock_gold_usage.get("read_only_layer_scoring_possible"):
        blinding_line = (
            "- Blinding: prediction content was finalized and SHA-256 locked "
            "before this external score. The source production workflow may have "
            "read registered gold for read-only layer scorecards before the "
            "projection lock; those scores did not feed back into extraction. "
            "This is not the stricter native lock-before-any-gold-read protocol."
        )
    else:
        blinding_line = (
            "- Blinding: gold was used only for PMID eligibility and count-field "
            "presence during selection; extraction exported no gold values or row "
            "counts, and predictions were locked before `score` opened gold."
        )
    token_lines = (
        [
            (
                f"- Exact API telemetry: **{token_usage.get('total_tokens', 0):,} "
                f"total tokens** ({token_usage.get('input_tokens', 0):,} input; "
                f"{token_usage.get('output_tokens', 0):,} output)."
            )
        ]
        if telemetry_available
        else [
            "- Exact API token and timing telemetry was not captured for this "
            "legacy production projection; zero placeholders must not be "
            "interpreted as zero cost."
        ]
    )
    timing_lines = (
        [
            (
                f"- Elapsed: **{timing['wall_seconds']:.1f}s wall clock**; "
                f"{overall['elapsed_seconds']:.1f}s summed per-paper route + read time."
            )
        ]
        if wall_timing_available
        else (
            [
                f"- Traced API duration: **{overall['elapsed_seconds']:.1f}s** "
                "summed across papers; end-to-end wall time was not captured."
            ]
            if telemetry_available
            else []
        )
    )
    native_trace_evidence_lines = (
        [
            "- `llm_traces/<GENE>/<PMID>/`: exact textual requests, safe parameters, raw provider response envelopes, parse attempts, and explicit route/final-selection events for each model call.",
            "- `llm_traces/trace_manifest.json`: SHA-256 inventory locked before gold scoring; provider-returned reasoning summaries are retained, while private hidden chain-of-thought is not available.",
            "- `llm_trace_report.html`: self-contained per-paper browser view of the locked trace timeline, prompts, responses, rationales, retries, and integrity state.",
        ]
        if native_traces_locked
        else [
            "- Exact per-call LLM traces are not attached to this evaluation lock; legacy runs require a rerun for a trace-complete audit."
        ]
    )
    trace_evidence_lines = (
        [
            "- The production `gvf-run` trace manifest for every gene, including its exact call/decision records and write-time digest index, is SHA-256-bound in `predictions.json` and `LOCK.json` and revalidated before scoring.",
            "- Each source run retains its own `llm_traces/<GENE>/<PMID>/` records and `llm_trace_report.html`; the evaluation projection does not copy or relabel those run-scoped records.",
        ]
        if production_traces_locked
        else native_trace_evidence_lines
    )
    lines = [
        f"# Codex extraction-blinded paper evaluation — `{report['run_id']}`",
        "",
        "## Technical summary",
        "",
        (
            f"This hash-locked run evaluated **{overall['papers']} papers** "
            f"(**{per_gene_label}**) after selecting only PMIDs with downloaded "
            f"source and gold assertions for carriers, affected, and unaffected. "
            f"Codex predictions were finalized before scoring."
        ),
        "",
        (
            f"- Variant precision **{format_rate(overall['precision'])}**, recall "
            f"**{format_rate(overall['recall'])}**, F1 "
            f"**{format_rate(overall['f1'])}** "
            f"({overall['tp']} TP, {overall['fp']} FP, {overall['fn']} FN)."
        ),
        (
            "- Precision versus counted extras "
            f"**{format_rate(counted_precision.get('precision_vs_counted_gold_pmids'))}** "
            f"({counted_precision.get('matched_rows', 0)} matched rows; "
            f"{counted_precision.get('counted_extra_rows', 0)} extra rows with "
            "patient counts). The stricter count-bearing-only diagnostic is "
            f"**{format_rate(counted_precision.get('precision_among_count_bearing_predictions'))}** "
            "and has a different numerator; it is not comparable to the "
            "repository's counted-extra precision floor."
        ),
        *token_lines,
        *timing_lines,
        (
            "- Notation twins merged before scoring: "
            f"**{overall.get('merged_notation_twin_rows', 0)}** same-paper "
            "prediction rows that were the same variant in another notation "
            "(equivalent-allele identity only; ambiguous or count-conflicting "
            "rows were left separate)."
        ),
        f"- Representation choices: {report['tools_used']}.",
        "",
        "## Blinding and scorer audit",
        "",
        f"- {selection_description}",
        blinding_line,
    ]
    audit = report.get("scoring_audit")
    if audit:
        raw = audit["pre_adjudication"]
        lines.extend(
            [
                (
                    "- The preserved pre-adjudication matcher scored "
                    f"{raw['tp']} TP / {raw['fp']} FP / {raw['fn']} FN "
                    f"(precision {format_rate(raw['precision'])}, recall "
                    f"{format_rate(raw['recall'])}, F1 {format_rate(raw['f1'])}). "
                    f"A post-lock notation audit recovered "
                    f"{audit['recovered_equivalent_matches']} equivalent labels; "
                    "no prediction text or count changed."
                ),
                (
                    "- `matcher_adjudication.csv` lists every recovered pair and "
                    "its equivalence class; raw matcher outputs remain preserved."
                ),
            ]
        )
    lines.extend(
        [
            "",
            "## Count fidelity",
            "",
            (
                "Count recall is the share of all gold count assertions for which the "
                "locked prediction supplied a value; MAE/RMSE are computed only where "
                "both gold and prediction supplied a value."
            ),
            "",
            "| field | supplied / gold assertions | count recall | MAE | RMSE |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for field in COUNT_FIELDS:
        metric = overall["count"][field]
        lines.append(
            f"| {field} | {metric['predicted']} / {metric['gold_asserted']} | "
            f"{format_rate(metric['recall'])} | {format_number(metric['mae'])} | "
            f"{format_number(metric['rmse'])} |"
        )

    lines.extend(
        [
            "",
            (
                'Gold encodes "no such individuals reported" as an explicit 0 '
                "while the pipeline deliberately abstains with NULL, so pooled "
                "count recall mixes that convention gap with real attribution "
                "misses. The stratified view separates them; the non-zero column "
                "is the actionable attribution number."
            ),
            "",
            (
                "| field | non-zero gold: supplied / asserted | non-zero recall "
                "| zero gold: supplied / asserted | zero recall |"
            ),
            "|---|---:|---:|---:|---:|",
        ]
    )
    for field in COUNT_FIELDS:
        metric = overall["count"][field]
        lines.append(
            f"| {field} | {metric.get('predicted_on_nonzero_gold', 'n/a')} / "
            f"{metric.get('gold_asserted_nonzero', 'n/a')} | "
            f"{format_rate(metric.get('recall_nonzero_gold'))} | "
            f"{metric.get('predicted_on_zero_gold', 'n/a')} / "
            f"{metric.get('gold_asserted_zero', 'n/a')} | "
            f"{format_rate(metric.get('recall_zero_gold'))} |"
        )

    lines.extend(
        [
            "",
            "## Per-gene results",
            "",
            "| gene | TP | FP | FN | precision | recall | F1 | precision vs counted extras | count-bearing-only precision | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|",
        ]
    )
    for gene in (g for g in GENES if g in report["by_gene"]):
        metric = report["by_gene"][gene]

        def count_cell(field: str) -> str:
            count = metric["count"][field]
            return (
                f"{format_rate(count['recall'])} / "
                f"{format_number(count['mae'])} / {format_number(count['rmse'])}"
            )

        lines.append(
            f"| {gene} | {metric['tp']} | {metric['fp']} | {metric['fn']} | "
            f"{format_rate(metric['precision'])} | {format_rate(metric['recall'])} | "
            f"{format_rate(metric['f1'])} | "
            f"{format_rate((metric.get('counted_precision') or {}).get('precision_vs_counted_gold_pmids'))} | "
            f"{format_rate((metric.get('counted_precision') or {}).get('precision_among_count_bearing_predictions'))} | "
            f"{count_cell('carriers')} | "
            f"{count_cell('affected')} | {count_cell('unaffected')} |"
        )

    lines.extend(
        [
            "",
            "## Per-paper results",
            "",
            "| gene | PMID | tool | TP | FP | FN | precision | recall | F1 | carrier recall / MAE | affected recall / MAE | unaffected recall / MAE | seconds | tokens |",
            "|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|",
        ]
    )
    for score in report["papers"]:

        def short_count(field: str) -> str:
            count = score["count"][field]
            return f"{format_rate(count['recall'])} / {format_number(count['mae'])}"

        seconds_cell = (
            f"{float(score.get('elapsed_seconds') or 0):.1f}"
            if telemetry_available
            else "n/a"
        )
        tokens_cell = (
            f"{int((score.get('token_usage') or {}).get('total_tokens') or 0):,}"
            if telemetry_available
            else "n/a"
        )
        lines.append(
            f"| {score['gene']} | {score['pmid']} | {score.get('tool') or 'n/a'} | "
            f"{score['tp']} | {score['fp']} | {score['fn']} | "
            f"{format_rate(score['precision'])} | {format_rate(score['recall'])} | "
            f"{format_rate(score['f1'])} | {short_count('carriers')} | "
            f"{short_count('affected')} | {short_count('unaffected')} | "
            f"{seconds_cell} | {tokens_cell} |"
        )

    lines.extend(["", "## Errors and representation choices", ""])
    for score in report["papers"]:
        lines.extend(
            [
                f"### {score['gene']} PMID {score['pmid']}",
                "",
                (
                    f"**{score.get('tool') or 'unspecified'}** — "
                    f"{score.get('tool_rationale') or 'No rationale recorded.'}"
                ),
                "",
            ]
        )
        if score["missed_gold"]:
            lines.append(f"- Missed gold variants: {', '.join(score['missed_gold'])}")
        if score["extra_predictions"]:
            lines.append(
                f"- Extra predictions: {', '.join(score['extra_predictions'])}"
            )
        if score["count_errors"]:
            rendered = "; ".join(
                f"{error['variant']} {error['field']} "
                f"{error['predicted']} vs {error['gold']} "
                f"(error {error['error']:+d})"
                for error in score["count_errors"]
            )
            lines.append(f"- Count disagreements: {rendered}")
        if not (
            score["missed_gold"] or score["extra_predictions"] or score["count_errors"]
        ):
            lines.append("- No scored variant or count disagreement.")
        lines.append("")

    lines.extend(
        [
            "## Scope, method, and limitations",
            "",
            f"- Population: {selection_population}; {per_gene_label}; every PMID has downloaded source and at least one gold assertion in each count field.",
            blinding_line,
            "- Variant metrics are micro-averaged over gold rows. Precision treats unmatched predictions as false positives, although the curated recall packet may omit some real variants.",
            "- Count MAE/RMSE are conditional on a supplied value. Count recall must be read alongside them because abstentions and missed variants are excluded from error magnitude.",
            "- Source acquisition and gold completeness are separate from model reading quality; abstract-only or incomplete source is retained and labeled rather than silently excluded.",
            "- The audited notation score is primary; the preserved raw score bounds sensitivity to post-lock matching adjudication.",
            "",
            "## Reproducibility and evidence",
            "",
            "- `selection.json`: selected PMIDs, source paths, source hashes, and available representations.",
            "- `predictions.json`: immutable per-paper tools, rationales, extracted variants, counts, evidence quotes, source locations, and telemetry when captured.",
            *trace_evidence_lines,
            "- `evidence.csv`: flat evidence ledger for every predicted variant.",
            "- `paper_metrics.csv`: exact per-paper metrics.",
            "- `LOCK.json`: SHA-256 digests proving prediction finalization before scoring.",
            "- `report.json`: complete machine-readable score, errors, timing, and token usage.",
            "- `matcher_adjudication.csv`: post-lock notation-equivalence audit; no extraction was edited.",
            "- `report_raw_matcher.json` and `report_raw_matcher.md`: preserved pre-adjudication score.",
            "- `validation_notes.md`: independent arithmetic, integrity checks, failure concentration, count outliers, and Claude comparison.",
            "- `model_comparison.csv`: compact Codex/Claude comparison with scorer and telemetry caveats.",
            "- `report_queries.sql`: executable DuckDB queries for the bounded analytical report datasets.",
            "",
            "## Recommended next steps",
            "",
            "1. Adjudicate extra predictions against the paper before treating precision as a production false-positive rate.",
            "2. Review count outliers by source location and distinguish model mistakes from gold disagreements.",
            "3. Add automatic fallback routing for data-rich papers that return zero or very few variants, then repeat with the same lock and count-recall definitions.",
            "",
            "## Further questions",
            "",
            "- Does table/PDF/OCR routing improve recall enough to justify its additional routing-call tokens?",
            "- How much of the residual error is source incompleteness versus count-role interpretation?",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def command_score(args) -> None:
    run_dir = args.run_dir
    selection_path, prediction_path, lock_path = (
        run_dir / "selection.json",
        run_dir / "predictions.json",
        run_dir / "LOCK.json",
    )
    if not lock_path.exists():
        raise SystemExit("refusing to score: predictions are not locked")
    lock = read_json(lock_path)
    if (
        digest(selection_path) != lock["selection_sha256"]
        or digest(prediction_path) != lock["predictions_sha256"]
    ):
        raise SystemExit("refusing to score: locked input digest mismatch")
    selection, predictions = read_json(selection_path), read_json(prediction_path)
    trace_manifest_digest = lock.get("llm_trace_manifest_sha256")
    if trace_manifest_digest:
        trace_root = run_dir / "llm_traces"
        trace_manifest_path = trace_root / TRACE_MANIFEST_NAME
        trace_report_path = run_dir / TRACE_REPORT_NAME
        if (
            not trace_manifest_path.is_file()
            or digest(trace_manifest_path) != trace_manifest_digest
        ):
            raise SystemExit("refusing to score: locked LLM trace manifest mismatch")
        trace_errors = validate_trace_manifest(
            trace_root, read_json(trace_manifest_path)
        )
        if trace_errors:
            raise SystemExit(
                "refusing to score: LLM trace integrity failed:\n- "
                + "\n- ".join(trace_errors)
            )
        trace_report_digest = lock.get("llm_trace_report_sha256")
        if trace_report_digest and (
            not trace_report_path.is_file()
            or digest(trace_report_path) != trace_report_digest
        ):
            raise SystemExit("refusing to score: locked LLM trace report mismatch")
    production_trace_locks, _production_trace_roots, production_trace_errors = (
        production_trace_lock_entries(predictions)
    )
    if production_trace_locks != (lock.get("production_trace_manifests") or []):
        production_trace_errors.append(
            "production trace lock entries differ from locked predictions"
        )
    if production_trace_errors:
        raise SystemExit(
            "refusing to score: production LLM trace integrity failed:\n- "
            + "\n- ".join(production_trace_errors)
        )
    pred_map = {(p["gene"], str(p["pmid"])): p for p in predictions["papers"]}
    scores = []
    for paper in selection["papers"]:
        key = (paper["gene"], paper["pmid"])
        scores.append(score_one(*key, pred_map[key], load_gold(args.gold_root, *key)))
    present_genes = [gene for gene in GENES if any(s["gene"] == gene for s in scores)]
    by_gene = {
        gene: aggregate([s for s in scores if s["gene"] == gene])
        for gene in present_genes
    }
    report = {
        "run_id": selection["run_id"],
        "seed": selection["seed"],
        "locked_at": lock["locked_at"],
        "scored_at": datetime.now(timezone.utc).isoformat(),
        "overall": aggregate(scores),
        "by_gene": by_gene,
        "papers": scores,
        "selection": selection_metadata(selection),
        "gold_sources": {
            gene: str(gold_csv_path(args.gold_root, gene)) for gene in present_genes
        },
        "tools_used": dict(Counter(s.get("tool") or "unspecified" for s in scores)),
        "token_usage": predictions.get("token_usage"),
        "timing": {
            "wall_seconds": float(predictions.get("extraction_elapsed_seconds") or 0),
            "summed_paper_seconds": sum(
                float(score.get("elapsed_seconds") or 0) for score in scores
            ),
            "started_at": predictions.get("extraction_started_at"),
            "completed_at": predictions.get("completed_at"),
        },
        "integrity": {
            "selection_sha256": lock["selection_sha256"],
            "predictions_sha256": lock["predictions_sha256"],
            "llm_trace_manifest_sha256": trace_manifest_digest,
            "llm_trace_report_sha256": lock.get("llm_trace_report_sha256"),
            "production_trace_manifests": production_trace_locks,
        },
        "blinding": selection.get("blinding"),
        "prelock_gold_usage": predictions.get("prelock_gold_usage")
        or (
            {
                "read_only_layer_scoring_possible": True,
                "scores_feed_back_into_predictions": False,
            }
            if int(predictions.get("schema_version") or 1) < 2
            else None
        ),
    }
    raw_report_path = run_dir / "report_raw_matcher.json"
    if raw_report_path.exists():
        raw_report = read_json(raw_report_path)
        adjudication_count = write_matcher_adjudication_csv(
            scores, raw_report, run_dir / "matcher_adjudication.csv"
        )
        report["scoring_audit"] = {
            "pre_adjudication": raw_report["overall"],
            "recovered_equivalent_matches": adjudication_count,
            "predictions_changed": False,
        }
    write_json(run_dir / "report.json", report)
    write_evidence_csv(predictions, run_dir / "evidence.csv")
    write_paper_metrics_csv(scores, run_dir / "paper_metrics.csv")
    write_markdown_report(report, run_dir / "report.md")
    print(run_dir / "report.json")


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--per-gene", type=int, default=5)
    p.add_argument("--minimum-chars", type=int, default=2000)
    p.add_argument("--corpus-root", type=Path, default=DEFAULT_CORPUS)
    p.add_argument("--gold-root", type=Path, default=DEFAULT_GOLD)
    p.add_argument(
        "--paper-manifest",
        type=Path,
        help="optional two-column GENE/PMID manifest; gold values remain hidden",
    )
    p.add_argument("--runs-dir", type=Path, default=HERE / "runs")
    p.add_argument(
        "--legacy-source-selection",
        action="store_true",
        help="ablation: take the first candidate rendering instead of the richest",
    )
    p.add_argument("--run-id", default=lambda: None)
    p.set_defaults(func=command_prepare)
    p = sub.add_parser("lock")
    p.add_argument("--run-dir", type=Path, required=True)
    p.set_defaults(func=command_lock)
    p = sub.add_parser("extract")
    p.add_argument("--run-dir", type=Path, required=True)
    p.add_argument("--model", default="gpt-5.6-sol")
    p.add_argument("--reasoning-effort", default="high")
    p.add_argument("--route-reasoning-effort", default="medium")
    p.add_argument("--max-output-tokens", type=int, default=24000)
    # Models that reason by default spend hidden tokens against this budget before
    # emitting the router JSON, so they need more headroom than 1600.
    p.add_argument("--route-max-output-tokens", type=int, default=1600)
    # Upper bound for the one truncation retry. 100000 is verified accepted by the
    # Azure grok-4.3 and gpt-5.6 deployments.
    p.add_argument("--max-output-tokens-ceiling", type=int, default=100000)
    p.add_argument(
        "--legacy-table-material",
        action="store_true",
        help="ablation: fold the artifact into the table string as the old code did",
    )
    p.add_argument("--max-source-chars", type=int, default=120000)
    p.add_argument("--max-artifact-chars", type=int, default=30000)
    p.add_argument("--route-preview-chars", type=int, default=6000)
    p.add_argument("--max-ocr-images", type=int, default=8)
    p.add_argument("--timeout", type=float, default=900)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=command_extract)
    p = sub.add_parser("score")
    p.add_argument("--run-dir", type=Path, required=True)
    p.add_argument("--gold-root", type=Path, default=DEFAULT_GOLD)
    p.set_defaults(func=command_score)
    return ap


def main() -> None:
    ap = parser()
    args = ap.parse_args()
    if args.command == "prepare" and not isinstance(args.run_id, str):
        args.run_id = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    args.func(args)


if __name__ == "__main__":
    main()
