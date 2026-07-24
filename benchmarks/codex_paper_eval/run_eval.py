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

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

GENES = ("SCN5A", "KCNH2", "KCNQ1", "RYR2")
COUNT_FIELDS = ("carriers", "affected", "unaffected")
DEFAULT_CORPUS = REPO / "corpus"
DEFAULT_GOLD = REPO / "gene_variant_fetcher_gold_standard" / "normalized"
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


def usable_sources(corpus: Path, gene: str, minimum_chars: int) -> list[dict]:
    papers = []
    for paper_dir in sorted(
        (corpus / gene).iterdir() if (corpus / gene).is_dir() else []
    ):
        if not paper_dir.is_dir() or not paper_dir.name.isdigit():
            continue
        pmid = paper_dir.name
        candidates = [
            paper_dir / f"{pmid}_FULL_CONTEXT.md",
            paper_dir / f"{pmid}_CLEANED.md",
        ]
        source = next(
            (
                p
                for p in candidates
                if p.is_file() and p.stat().st_size >= minimum_chars
            ),
            None,
        )
        if source is None:
            continue
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


def gold_count_eligible_pmids(gold_root: Path, gene: str) -> set[str]:
    """Return PMIDs with gold rows and at least one assertion for every count field.

    Only PMID membership and field presence are used during selection. Gold values
    and gold row counts are never written into the selection or extraction prompt.
    """
    path = gold_root / f"{gene}_recall_input.csv"
    coverage: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            pmid = str(row.get("pmid", "")).strip()
            if not pmid or not str(row.get("variant", "")).strip():
                continue
            for field in COUNT_FIELDS:
                if str(row.get(field, "")).strip() != "":
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
    for gene in GENES:
        source_pool = {
            paper["pmid"]: paper
            for paper in usable_sources(args.corpus_root, gene, args.minimum_chars)
        }
        eligible_pmids = gold_count_eligible_pmids(args.gold_root, gene)
        pools[gene] = {
            pmid: paper for pmid, paper in source_pool.items() if pmid in eligible_pmids
        }
        eligible[gene] = len(pools[gene])

    if args.paper_manifest:
        requested = read_paper_manifest(args.paper_manifest)
        for gene, pmid in requested:
            if pmid not in pools[gene]:
                raise SystemExit(
                    f"{gene} {pmid}: missing usable source or complete gold count coverage"
                )
            selected.append(pools[gene][pmid])
    else:
        for gene in GENES:
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
        "schema_version": 1,
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
                "notes": None,
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
    for p in predictions.get("papers", []):
        label = f"{p.get('gene')}:{p.get('pmid')}"
        for key in ("tool", "tool_rationale", "elapsed_seconds", "source_completeness"):
            if p.get(key) in (None, ""):
                errors.append(f"{label}: missing {key}")
        for i, row in enumerate(p.get("variants", [])):
            if not (row.get("variant") or "").strip():
                errors.append(f"{label} variant[{i}]: missing variant")
            if not (row.get("evidence") or "").strip():
                errors.append(f"{label} variant[{i}]: missing evidence")
            if not (row.get("source_location") or "").strip():
                errors.append(f"{label} variant[{i}]: missing source_location")
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
        if not usage.get("telemetry_available") or not usage.get("total_tokens"):
            errors.append(f"{label}: missing exact token telemetry")
    return errors


def command_lock(args) -> None:
    run_dir = args.run_dir
    lock_path = run_dir / "LOCK.json"
    if lock_path.exists():
        raise SystemExit(f"already locked: {lock_path}")
    selection_path = run_dir / "selection.json"
    prediction_path = run_dir / "predictions.json"
    selection, predictions = read_json(selection_path), read_json(prediction_path)
    errors = validate_predictions(selection, predictions)
    errors.extend(selection_material_errors(selection))
    if errors:
        raise SystemExit("prediction validation failed:\n- " + "\n- ".join(errors))
    lock = {
        "locked_at": datetime.now(timezone.utc).isoformat(),
        "selection_sha256": digest(selection_path),
        "predictions_sha256": digest(prediction_path),
        "statement": (
            "Predictions finalized before gold values or gold row counts were "
            "exposed to extraction; score is the first phase that reads those values."
        ),
    }
    write_json(lock_path, lock)
    prediction_path.chmod(0o444)
    print(lock_path)


ROUTE_INSTRUCTIONS = """You are routing one biomedical paper for blinded curation.
Use only the representation previews supplied below. Do not use databases, prior
extractions, benchmarks, or gold standards.

Target gene: {gene}. PMID: {pmid}.

Choose exactly one authoritative representation:
- text: running full text is the clearest source.
- table: structured table/artifact rows carry the variant-level person counts.
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


EXTRACTION_INSTRUCTIONS = """You are independently curating one biomedical paper.
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
Use null when a count cannot be determined, and 0 only when zero is explicit.
Do not turn cohort size, alleles, families, events, or functional experiments
into person counts. Do not double-count the same people across prose and tables.

The selected representation is supplied below. Cite concise evidence and a
specific section/table/page/figure location for every row. Return JSON only:
{{
  "notes": "...",
  "variants": [{{
    "variant": "...", "carriers": null, "affected": null, "unaffected": null,
    "evidence": "...", "source_location": "..."
  }}]
}}

LOCAL MATERIAL:
{material}
"""


def parse_json_response(text: str) -> dict:
    value = text.strip()
    if value.startswith("```"):
        value = value.split("\n", 1)[1].rsplit("```", 1)[0]
    return json.loads(value)


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
        is_table_line = line.count("|") >= 2 or bool(
            re.match(r"^\s*</?(?:table|tr|td|th)\b", line, flags=re.I)
        )
        if is_table_line:
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


def command_extract(args) -> None:
    """Route and run isolated per-paper reads through the Azure OpenAI endpoint."""
    from openai import OpenAI

    run_dir = args.run_dir
    if (run_dir / "LOCK.json").exists():
        raise SystemExit("refusing to extract: run is already locked")
    selection = read_json(run_dir / "selection.json")
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
        if artifact_text:
            table_text = (
                f"### Parsed artifact data\n\n{artifact_text}\n\n"
                f"### Tables preserved in text\n\n{table_text}"
            )
        pdf_text = extract_pdf_text(paper.get("pdfs") or [], args.max_source_chars)
        figures = [Path(path) for path in paper.get("figures") or []]
        representations = {
            "text": bool(source_text.strip()),
            "table": bool(table_text.strip()),
            "pdf": bool(pdf_text.strip()),
            "ocr": bool(figures),
        }
        target["representations_available"] = [
            name for name, available in representations.items() if available
        ]
        catalog = (
            f"Available={target['representations_available']}\n\n"
            f"## TEXT PREVIEW\n{targeted_preview(source_text, paper['gene'], args.route_preview_chars)}\n\n"
            f"## TABLE PREVIEW\n{targeted_preview(table_text, paper['gene'], args.route_preview_chars)}\n\n"
            f"## PDF PREVIEW\n{targeted_preview(pdf_text, paper['gene'], args.route_preview_chars)}\n\n"
            f"## OCR INVENTORY\n"
            + "\n".join(path.name for path in figures[: args.max_ocr_images])
        )
        route_prompt = ROUTE_INSTRUCTIONS.format(
            gene=paper["gene"], pmid=paper["pmid"], catalog=catalog
        )
        started = time.monotonic()
        print(
            f"[{index}/{total_papers}] route {paper['gene']} {paper['pmid']}",
            flush=True,
        )
        route_response = client.responses.create(
            model=args.model,
            input=route_prompt,
            reasoning={"effort": args.route_reasoning_effort},
            max_output_tokens=1600,
        )
        add_usage(
            predictions,
            target,
            route_response,
            args.model,
            args.route_reasoning_effort,
        )
        write_json(predictions_path, predictions)
        route = parse_json_response(route_response.output_text)
        tool = str(route.get("tool", "")).lower()
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
            tool = available_fallback

        if tool == "table":
            material = (
                f"### Structured table/artifact material\n\n{table_text[: args.max_source_chars]}"
                f"\n\n### Supporting targeted text\n\n"
                f"{targeted_preview(source_text, paper['gene'], args.route_preview_chars)}"
            )
        elif tool == "pdf":
            material = f"### PDF layout text\n\n{pdf_text[: args.max_source_chars]}"
        elif tool == "ocr":
            material = "### Supporting targeted text\n\n" + targeted_preview(
                source_text, paper["gene"], args.route_preview_chars
            )
        else:
            material = f"### Full/partial running text\n\n{source_text[: args.max_source_chars]}"

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
        extraction_response = client.responses.create(
            model=args.model,
            input=response_input,
            reasoning={"effort": args.reasoning_effort},
            max_output_tokens=args.max_output_tokens,
        )
        add_usage(
            predictions,
            target,
            extraction_response,
            args.model,
            args.reasoning_effort,
        )
        write_json(predictions_path, predictions)
        result = parse_json_response(extraction_response.output_text)
        elapsed = time.monotonic() - started
        target.update(result)
        target["tool"] = tool
        target["tool_rationale"] = route.get("tool_rationale")
        target["source_completeness"] = route.get("source_completeness")
        target["elapsed_seconds"] = round(elapsed, 3)
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


def to_int(value):
    if value is None or str(value).strip() == "":
        return None
    try:
        return int(float(str(value)))
    except ValueError:
        return None


def load_gold(gold_root: Path, gene: str, pmid: str) -> list[dict]:
    path = gold_root / f"{gene}_recall_input.csv"
    with path.open(newline="") as fh:
        rows = []
        for row in csv.DictReader(fh):
            if str(row.get("pmid", "")).strip() != pmid:
                continue
            rows.append(
                {
                    "variant": str(row.get("variant", "")).strip(),
                    **{field: to_int(row.get(field)) for field in COUNT_FIELDS},
                }
            )
        return rows


def variant_candidates(value: str, gene: str) -> list[str]:
    """Return embedded notations without changing the locked prediction."""
    candidates = [value]
    patterns = (
        r"\bp\.\(?[A-Z][a-z]{2}\d+(?:[A-Z][a-z]{2}|Ter|fs[^\s,;)]*)\)?",
        (
            r"\bc\.[0-9*?+-]+(?:_[0-9*?+-]+)?"
            r"(?:delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+|[ACGT]>[ACGT])"
        ),
        r"\b[A-Z]\d{1,5}(?:[A-Z]|\*|X)\b",
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
    candidates.extend(
        f"c.{token}"
        for token in re.findall(
            (
                r"(?<![A-Za-z.])([0-9*?+-]+(?:_[0-9*?+-]+)?"
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
    for residue, position in re.findall(
        r"\bfs([A-Z])(\d{1,5})(?:[/+*]?\d+|aa)*\b", value, flags=re.I
    ):
        candidates.append(f"{residue.upper()}{position}fsX")
    for residue, position in re.findall(
        r"\bfs([A-Z][a-z]{2})(\d{1,5})(?:[/+*]?\d+|aa)*\b", value
    ):
        if left := AA3_TO_1.get(residue):
            candidates.append(f"{left}{position}fsX")
    for residue, position in re.findall(
        r"\b([A-Z])(\d{1,5})(?:sp|splice)\b", value, flags=re.I
    ):
        candidates.append(f"{residue.upper()}{position}SP")
    for position, residue in re.findall(r"\bdel(\d{1,5})([A-Z])\b", value, flags=re.I):
        candidates.append(f"{residue.upper()}{position}del")
    for residue, position in re.findall(r"\bdel([A-Z])(\d{1,5})\b", value, flags=re.I):
        candidates.append(f"{residue.upper()}{position}del")
    structural = re.search(r"exons?\s+(\d+)\s*[–—-]\s*(\d+)", value, flags=re.I)
    if structural and re.search(r"\bdup(?:lication)?\b", value, flags=re.I):
        candidates.append(f"EXON{structural.group(1)}_{structural.group(2)}DUP")
    single_exon_del = re.search(r"\bexon\s+(\d+)\b.*\bdel", value, flags=re.I)
    if single_exon_del:
        candidates.append(f"EXON{single_exon_del.group(1)}DEL")
    if gene == "SCN5A" and re.search(r"(?:Δ|DELTA)\s*KPQ", value, flags=re.I):
        candidates.append("K1505_Q1507del")
    if gene == "KCNH2":
        insertion = re.search(
            r"\bINS\s+[ACGT]+\s+(\d+)(?:\s*[–—-]\s*\d+)?",
            value,
            flags=re.I,
        )
        if insertion:
            candidates.append(f"G{insertion.group(1)}ins")
    return list(dict.fromkeys(c.strip() for c in candidates))


def cdna_indel_position(value: str) -> int | None:
    match = re.search(
        r"(?:^|[^A-Za-z])c?\.?(\d+)(?:_\d+)?(?:del|dup|ins)",
        value,
        flags=re.I,
    )
    return int(match.group(1)) // 3 + 1 if match else None


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


def matches(a: str, b: str, gene: str) -> bool:
    from utils.variant_normalizer import normalize_variant, variants_match

    if not a or not b:
        return False
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


def score_one(gene: str, pmid: str, predicted: dict, gold: list[dict]) -> dict:
    pred = predicted.get("variants", [])
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
            fps.append(p["variant"])
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
    for field in COUNT_FIELDS:
        gold_asserted = sum(1 for row in gold if row.get(field) is not None)
        observed = [
            (p.get(field), g.get(field), p["variant"])
            for p, g in pairs
            if g.get(field) is not None and p.get(field) is not None
        ]
        errors = [pv - gv for pv, gv, _ in observed]
        count[field] = {
            "gold_asserted": gold_asserted,
            "predicted": len(observed),
            "recall": len(observed) / gold_asserted if gold_asserted else None,
            "mae": statistics.fmean(abs(e) for e in errors) if errors else None,
            "rmse": math.sqrt(statistics.fmean(e * e for e in errors))
            if errors
            else None,
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
        "token_usage": predicted.get("token_usage"),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "count": count,
        "matched_variants": [
            {"predicted": p["variant"], "gold": g["variant"]} for p, g in pairs
        ],
        "missed_gold": [g["variant"] for g in misses],
        "extra_predictions": fps,
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
        result["count"][field] = {
            "gold_asserted": asserted,
            "predicted": observed,
            "recall": observed / asserted if asserted else None,
            "mae": abs_sum / observed if observed else None,
            "rmse": math.sqrt(sq_sum / observed) if observed else None,
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
            ]
        )
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for score in scores:
            usage = score.get("token_usage") or {}
            row = {
                **{key: score.get(key) for key in columns[:11]},
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
    token_usage = report["token_usage"] or {}
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
            f"- Exact API telemetry: **{token_usage.get('total_tokens', 0):,} total "
            f"tokens** ({token_usage.get('input_tokens', 0):,} input; "
            f"{token_usage.get('output_tokens', 0):,} output)."
        ),
        (
            f"- Elapsed: **{report['timing']['wall_seconds']:.1f}s wall clock**; "
            f"{overall['elapsed_seconds']:.1f}s summed per-paper route + read time."
        ),
        f"- Representation choices: {report['tools_used']}.",
        "",
        "## Blinding and scorer audit",
        "",
        f"- {selection_description}",
        (
            "- `selection.json` contains source metadata and hashes but no gold "
            "values or gold row counts. `predictions.json` was made read-only and "
            "SHA-256 locked before scoring first opened the gold CSVs."
        ),
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
            "## Per-gene results",
            "",
            "| gene | TP | FP | FN | precision | recall | F1 | carrier count recall / MAE / RMSE | affected count recall / MAE / RMSE | unaffected count recall / MAE / RMSE |",
            "|---|---:|---:|---:|---:|---:|---:|---|---|---|",
        ]
    )
    for gene in GENES:
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
            f"{format_rate(metric['f1'])} | {count_cell('carriers')} | "
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

        lines.append(
            f"| {score['gene']} | {score['pmid']} | {score.get('tool') or 'n/a'} | "
            f"{score['tp']} | {score['fp']} | {score['fn']} | "
            f"{format_rate(score['precision'])} | {format_rate(score['recall'])} | "
            f"{format_rate(score['f1'])} | {short_count('carriers')} | "
            f"{short_count('affected')} | {short_count('unaffected')} | "
            f"{float(score.get('elapsed_seconds') or 0):.1f} | "
            f"{int((score.get('token_usage') or {}).get('total_tokens') or 0):,} |"
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
            "- Blinding: gold was used only for PMID eligibility and count-field presence during selection; extraction exported no gold values or row counts, and predictions were made read-only and SHA-256 locked before `score` opened gold.",
            "- Variant metrics are micro-averaged over gold rows. Precision treats unmatched predictions as false positives, although the curated recall packet may omit some real variants.",
            "- Count MAE/RMSE are conditional on a supplied value. Count recall must be read alongside them because abstentions and missed variants are excluded from error magnitude.",
            "- Source acquisition and gold completeness are separate from model reading quality; abstract-only or incomplete source is retained and labeled rather than silently excluded.",
            "- The audited notation score is primary; the preserved raw score bounds sensitivity to post-lock matching adjudication.",
            "",
            "## Reproducibility and evidence",
            "",
            "- `selection.json`: selected PMIDs, source paths, source hashes, and available representations.",
            "- `predictions.json`: immutable per-paper tools, rationales, extracted variants, counts, evidence quotes, source locations, elapsed time, and token telemetry.",
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
    pred_map = {(p["gene"], str(p["pmid"])): p for p in predictions["papers"]}
    scores = []
    for paper in selection["papers"]:
        key = (paper["gene"], paper["pmid"])
        scores.append(score_one(*key, pred_map[key], load_gold(args.gold_root, *key)))
    by_gene = {
        gene: aggregate([s for s in scores if s["gene"] == gene]) for gene in GENES
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
        },
        "blinding": selection.get("blinding"),
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
