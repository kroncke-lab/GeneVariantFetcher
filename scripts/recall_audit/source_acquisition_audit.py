#!/usr/bin/env python3
"""Build a gold-free source acquisition and refresh worklist for one run.

The recall-disagreement reports are useful when a curated gold standard exists,
but a new gene needs the same operational loop without knowing the answers:

* Which PMIDs have no usable full text and should be fetched?
* Which PMIDs have usable source text but stale/missing/zero-variant extraction
  JSONs and should be replayed from source?
* How many run PMIDs are currently backed by usable full text?

This script reads only run-local artifacts: ``pmc_fulltext/``, extraction JSONs,
``paywalled_missing.csv``, source-completeness JSON, and PMID list files. It
does not read recall outputs or gold standards.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any
from urllib.parse import unquote

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from scripts.recall_audit.common import (  # noqa: E402
        context_status,
        normalize_pmid,
        repo_path,
        source_has_unresolved_variant_supplement_refs,
        write_csv_rows,
    )
except ModuleNotFoundError:  # pragma: no cover
    from common import (  # type: ignore  # noqa: E402
        context_status,
        normalize_pmid,
        repo_path,
        source_has_unresolved_variant_supplement_refs,
        write_csv_rows,
    )

DATA_AVAILABLE_STATUSES = {
    "recovered_pmc",
    "recovered_browser",
    "recovered_supplement_only",
}
MISSING_SOURCE_STATUSES = {"abstract_only", "paywall_stub", "not_attempted"}

FIELDS = [
    "gene",
    "pmid",
    "action",
    "route",
    "priority_score",
    "source_status",
    "source_path",
    "source_bytes",
    "available_context_path",
    "available_context_bytes",
    "extraction_path",
    "extraction_variants",
    "extraction_abstract_only",
    "source_unbound",
    "source_sha_mismatch",
    "missing_variant_supplement",
    "zero_variant_qc",
    "single_carrier_qc",
    "doi",
    "publisher_hint",
    "notes",
]

FETCH_INPUT_FIELDS = ["PMID", "DOI", "route"]
SOURCE_OVERRIDE_FIELDS = [
    "gene",
    "pmid",
    "action",
    "route",
    "available_context_path",
    "available_context_bytes",
    "notes",
]

DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\"'<>]+", re.IGNORECASE)


def _clean_doi(value: str) -> str:
    doi = unquote(value.strip())
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi, flags=re.IGNORECASE)
    return doi.strip().strip(".,;:)]}>")


def _extract_doi(value: Any) -> str:
    if value is None:
        return ""
    match = DOI_RE.search(unquote(str(value)))
    return _clean_doi(match.group(0)) if match else ""


def infer_publisher(doi: str) -> str:
    lower = doi.lower()
    if lower.startswith("10.1016/"):
        return "elsevier"
    if lower.startswith(("10.1002/", "10.1111/")):
        return "wiley"
    if lower.startswith("10.1007/"):
        return "springer"
    if lower.startswith("10.1038/"):
        return "springer_nature"
    if lower.startswith("10.1093/"):
        return "oxford"
    if lower.startswith("10.1089/"):
        return "liebert_sage"
    if lower.startswith("10.1161/"):
        return "aha"
    return "unknown"


def fetch_route_for_doi(doi: str) -> tuple[str, str]:
    publisher = infer_publisher(doi)
    if publisher == "elsevier":
        return "fetch_elsevier_insttoken", publisher
    if publisher == "wiley":
        return "fetch_wiley_tdm", publisher
    if publisher == "springer":
        return "fetch_springer_api", publisher
    if publisher == "springer_nature":
        return "fetch_springer_nature_or_browser", publisher
    if publisher == "oxford":
        return "fetch_oxford_or_browser", publisher
    if publisher == "aha":
        return "fetch_browser_aha", publisher
    if publisher == "liebert_sage":
        return "blocked_liebert_sage_cloudflare", publisher
    return "fetch_paywalled_generic", publisher


def _column(fieldnames: list[str], *candidates: str) -> str | None:
    lower_to_field = {name.strip().lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate in lower_to_field:
            return lower_to_field[candidate]
    return None


def load_pmid_doi_csv(path: Path) -> dict[str, str]:
    """Read PMID -> DOI from heterogeneous fetch/paywall input CSVs."""
    doi_by_pmid: dict[str, str] = {}
    if not path.exists():
        return doi_by_pmid
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return doi_by_pmid
        pmid_col = _column(reader.fieldnames, "pmid", "pubmed_id", "pubmed", "id")
        if pmid_col is None:
            return doi_by_pmid
        doi_col = _column(reader.fieldnames, "doi", "article_doi")
        scan_cols = [
            field
            for field in ("URL", "url", "Reason", "reason", "final_url", "path")
            if field in reader.fieldnames
        ]
        for row in reader:
            pmid = normalize_pmid(row.get(pmid_col))
            if not pmid or pmid in doi_by_pmid:
                continue
            doi = _clean_doi(row.get(doi_col) or "") if doi_col else ""
            if not doi:
                for field in scan_cols:
                    doi = _extract_doi(row.get(field))
                    if doi:
                        break
            doi_by_pmid[pmid] = doi
    return doi_by_pmid


def _doi_from_source_artifacts(pmid: str, source_file: Path | None) -> str:
    if source_file is None:
        return ""
    parent = source_file.parent
    candidates = [
        parent / pmid / "result.json",
        parent / f"{pmid}_artifacts.json",
    ]
    for candidate in candidates:
        if not candidate.exists():
            continue
        payload = _json_load(candidate)
        doi = _clean_doi(str(payload.get("doi") or ""))
        if doi:
            return doi
    return ""


def _doi_from_source_text(source_file: Path | None) -> str:
    """Extract the article DOI from an on-disk source / abstract file.

    Abstract-only stubs carry no structured DOI, but the PubMed abstract they
    do contain prints the article DOI on its NLM citation line
    (``doi: 10.xxxx/...``). Prefer that ``doi:``-labelled match so we never grab
    a DOI from the reference list, then fall back to the first DOI in the text.
    Without this, a stub reaches ``fetch_paywalled`` with an empty DOI and is
    dropped before any publisher/proxy route is tried.
    """
    if source_file is None or not source_file.exists():
        return ""
    try:
        text = source_file.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    labelled = re.search(r"doi:\s*(10\.\d{4,9}/[^\s\"'<>]+)", text, flags=re.IGNORECASE)
    if labelled:
        return _clean_doi(labelled.group(1))
    return _extract_doi(text)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_load(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _variant_count(data: dict[str, Any]) -> int:
    variants = data.get("variants")
    if isinstance(variants, list):
        return len(variants)
    metadata = data.get("extraction_metadata") or {}
    try:
        return int(metadata.get("total_variants_found") or 0)
    except (TypeError, ValueError):
        return 0


def _metadata_mentions_abstract_only(metadata: dict[str, Any]) -> bool:
    fields = [str(metadata.get(key) or "") for key in ("source_type", "notes")]
    challenges = metadata.get("challenges")
    if isinstance(challenges, list):
        fields.extend(str(item) for item in challenges)
    elif challenges:
        fields.append(str(challenges))
    text = " ".join(fields).lower()
    return (
        bool(metadata.get("abstract_only"))
        or "abstract-only" in text
        or "abstract only" in text
        or "full text unavailable" in text
        or "full text not available" in text
        or "full text could not be retrieved" in text
    )


def _metadata_source_unbound(
    metadata: dict[str, Any], *, output_file: Path, extraction_dir: Path
) -> bool:
    source_file = metadata.get("source_file")
    if not source_file:
        return False
    source_path = Path(str(source_file)).expanduser()
    if source_path.suffix.lower() == ".json":
        return True
    try:
        source_path.resolve().relative_to(extraction_dir.resolve())
        return True
    except (OSError, ValueError):
        pass
    try:
        return source_path.resolve() == output_file.resolve()
    except OSError:
        return False


def _read_source_completeness(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _pmids_from_text_file(path: Path) -> set[str]:
    pmids: set[str] = set()
    if not path.exists():
        return pmids
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        value = line.split("#", 1)[0].strip()
        if value:
            pmid = normalize_pmid(value.split()[0])
            if pmid:
                pmids.add(pmid)
    return pmids


def collect_pmids(
    *,
    run_dir: Path,
    gene: str,
    extraction_dir: Path,
    harvest_dir: Path,
    paywalled_csv: Path,
    source_completeness: dict[str, Any],
    include_discovery_pmids: bool = False,
) -> set[str]:
    pmids: set[str] = set()
    if include_discovery_pmids:
        for path in run_dir.glob("*_pmids.txt"):
            pmids.update(_pmids_from_text_file(path))
    for path in extraction_dir.glob(f"{gene}_PMID_*.json"):
        match = re.search(r"_PMID_(\d+)\.json$", path.name)
        if match:
            pmids.add(match.group(1))
    for path in harvest_dir.glob("*_FULL_CONTEXT.md"):
        pmids.add(normalize_pmid(path.name.removesuffix("_FULL_CONTEXT.md")))
    if paywalled_csv.exists():
        pmids.update(load_pmid_doi_csv(paywalled_csv))
    for key in ("abstract_only_pmids", "zero_variant_pmids", "single_carrier_pmids"):
        values = source_completeness.get(key) or []
        if isinstance(values, list):
            pmids.update(
                normalize_pmid(value) for value in values if normalize_pmid(value)
            )
    return {pmid for pmid in pmids if pmid}


def _source_candidates_for_pmid(harvest_dir: Path, pmid: str) -> list[Path]:
    """Return same-PMID markdown sources in extraction priority order."""
    return [
        path
        for suffix in ("_DATA_ZONES.md", "_CLEANED.md", "_FULL_CONTEXT.md")
        if (path := harvest_dir / f"{pmid}{suffix}").exists()
    ]


def _best_source_for_pmid(harvest_dir: Path, pmid: str) -> Path | None:
    """Pick the best run-local source without letting stubs hide real full text."""
    candidates = _source_candidates_for_pmid(harvest_dir, pmid)
    for path in candidates:
        if context_status(path) in DATA_AVAILABLE_STATUSES:
            return path
    if candidates:
        return max(candidates, key=lambda path: (path.stat().st_size, str(path)))
    return None


def _source_sha_mismatch(metadata: dict[str, Any], source_file: Path | None) -> bool:
    if source_file is None or not source_file.exists():
        return False
    existing_sha = metadata.get("source_sha256")
    return bool(existing_sha and existing_sha != _sha256(source_file))


def classify_pmid(
    *,
    gene: str,
    pmid: str,
    harvest_dir: Path,
    extraction_dir: Path,
    doi: str,
    zero_variant_pmids: set[str],
    single_carrier_pmids: set[str],
) -> dict[str, Any]:
    source_file = _best_source_for_pmid(harvest_dir, pmid)
    status = context_status(source_file)
    if not doi:
        doi = _doi_from_source_artifacts(pmid, source_file)
    if not doi:
        doi = _doi_from_source_text(source_file)
    extraction_file = extraction_dir / f"{gene}_PMID_{pmid}.json"
    data = _json_load(extraction_file) if extraction_file.exists() else {}
    metadata = data.get("extraction_metadata") or {}
    variants = _variant_count(data)
    extraction_abstract_only = _metadata_mentions_abstract_only(metadata)
    source_unbound = (
        _metadata_source_unbound(
            metadata, output_file=extraction_file, extraction_dir=extraction_dir
        )
        if metadata
        else False
    )
    sha_mismatch = _source_sha_mismatch(metadata, source_file)
    zero_variant_qc = pmid in zero_variant_pmids or (
        source_file is not None and status in DATA_AVAILABLE_STATUSES and variants == 0
    )
    missing_variant_supplement = (
        source_file is not None
        and status in DATA_AVAILABLE_STATUSES
        and source_has_unresolved_variant_supplement_refs(source_file, gene)
    )
    single_carrier_qc = pmid in single_carrier_pmids
    notes: list[str] = []

    route = ""
    publisher = infer_publisher(doi) if doi else ""
    if status in MISSING_SOURCE_STATUSES:
        if doi:
            route, publisher = fetch_route_for_doi(doi)
            action = "manual_or_blocked" if route.startswith("blocked_") else "fetch"
        else:
            action = "fetch"
            route = "doi_lookup_then_fetch"
        if action == "manual_or_blocked":
            notes.append("no usable run-local full text; known blocked publisher")
        else:
            notes.append("no usable run-local full text")
    elif missing_variant_supplement:
        if doi:
            route, publisher = fetch_route_for_doi(doi)
            action = "manual_or_blocked" if route.startswith("blocked_") else "fetch"
        else:
            action = "fetch"
            route = "doi_lookup_then_fetch"
        if action == "manual_or_blocked":
            notes.append(
                "usable main text references a missing target-gene variant supplement; known blocked publisher"
            )
        else:
            notes.append(
                "usable main text references a missing target-gene variant supplement"
            )
    elif not extraction_file.exists():
        action = "refresh_replay"
        route = "missing_extraction_for_usable_source"
        notes.append("usable source exists but extraction JSON is missing")
    elif extraction_abstract_only:
        action = "refresh_replay"
        route = "stale_abstract_extraction_for_usable_source"
        notes.append("extraction metadata is abstract-only but usable source exists")
    elif source_unbound:
        action = "refresh_replay"
        route = "unbound_extraction_source"
        notes.append("extraction source metadata is not bound to source text")
    elif sha_mismatch:
        action = "refresh_replay"
        route = "source_fingerprint_mismatch"
        notes.append("source fingerprint changed since extraction")
    elif zero_variant_qc and status in DATA_AVAILABLE_STATUSES:
        action = "refresh_replay"
        route = "zero_variant_fulltext_qc"
        notes.append("usable full text produced zero variants")
    elif single_carrier_qc:
        action = "inspect"
        route = "single_carrier_count_qc"
        notes.append("all extracted variants have single-carrier counts")
    else:
        action = "none"
        route = "covered"

    priority = 0
    if action == "fetch":
        priority += 100
    if action == "refresh_replay":
        priority += 80
    if zero_variant_qc:
        priority += 30
    if extraction_abstract_only or source_unbound or sha_mismatch:
        priority += 20
    if source_file and source_file.exists():
        priority += min(source_file.stat().st_size // 5000, 30)
    if single_carrier_qc:
        priority += 10
    if missing_variant_supplement:
        priority += 35

    return {
        "gene": gene,
        "pmid": pmid,
        "action": action,
        "route": route,
        "priority_score": priority,
        "source_status": status,
        "source_path": str(source_file) if source_file else "",
        "source_bytes": source_file.stat().st_size if source_file else 0,
        "available_context_path": str(source_file)
        if source_file and status in DATA_AVAILABLE_STATUSES
        else "",
        "available_context_bytes": source_file.stat().st_size
        if source_file and status in DATA_AVAILABLE_STATUSES
        else 0,
        "extraction_path": str(extraction_file) if extraction_file.exists() else "",
        "extraction_variants": variants,
        "extraction_abstract_only": extraction_abstract_only,
        "source_unbound": source_unbound,
        "source_sha_mismatch": sha_mismatch,
        "missing_variant_supplement": missing_variant_supplement,
        "zero_variant_qc": zero_variant_qc,
        "single_carrier_qc": single_carrier_qc,
        "doi": doi,
        "publisher_hint": publisher,
        "notes": "; ".join(notes),
    }


def build_audit(
    *,
    gene: str,
    run_dir: Path,
    harvest_dir: Path | None = None,
    extraction_dir: Path | None = None,
    paywalled_csv: Path | None = None,
    source_completeness_path: Path | None = None,
    include_discovery_pmids: bool = False,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    gene = gene.upper()
    harvest_dir = harvest_dir or run_dir / "pmc_fulltext"
    extraction_dir = extraction_dir or run_dir / "extractions"
    paywalled_csv = paywalled_csv or harvest_dir / "paywalled_missing.csv"
    source_completeness_path = (
        source_completeness_path or run_dir / "source_completeness.json"
    )
    source_completeness = _read_source_completeness(source_completeness_path)
    doi_by_pmid = load_pmid_doi_csv(paywalled_csv)
    zero_variant_pmids = {
        normalize_pmid(value)
        for value in source_completeness.get("zero_variant_pmids", []) or []
        if normalize_pmid(value)
    }
    single_carrier_pmids = {
        normalize_pmid(value)
        for value in source_completeness.get("single_carrier_pmids", []) or []
        if normalize_pmid(value)
    }
    pmids = collect_pmids(
        run_dir=run_dir,
        gene=gene,
        extraction_dir=extraction_dir,
        harvest_dir=harvest_dir,
        paywalled_csv=paywalled_csv,
        source_completeness=source_completeness,
        include_discovery_pmids=include_discovery_pmids,
    )
    rows = [
        classify_pmid(
            gene=gene,
            pmid=pmid,
            harvest_dir=harvest_dir,
            extraction_dir=extraction_dir,
            doi=doi_by_pmid.get(pmid, ""),
            zero_variant_pmids=zero_variant_pmids,
            single_carrier_pmids=single_carrier_pmids,
        )
        for pmid in sorted(pmids)
    ]
    rows.sort(
        key=lambda row: (
            row["action"] == "none",
            row["action"] == "inspect",
            -int(row["priority_score"] or 0),
            row["pmid"],
        )
    )
    summary = summarize(rows)
    summary["run_dir"] = str(run_dir)
    summary["harvest_dir"] = str(harvest_dir)
    summary["extraction_dir"] = str(extraction_dir)
    summary["paywalled_csv"] = str(paywalled_csv) if paywalled_csv.exists() else None
    summary["include_discovery_pmids"] = include_discovery_pmids
    summary["source_completeness"] = (
        str(source_completeness_path) if source_completeness_path.exists() else None
    )
    return rows, summary


def _coverage_item(count: int, total: int) -> dict[str, Any]:
    return {
        "pmids": count,
        "total_pmids": total,
        "coverage": round(count / total, 6) if total else None,
    }


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    total = len(rows)
    by_action = Counter(str(row["action"]) for row in rows)
    by_route = Counter(str(row["route"]) for row in rows)
    by_source_status = Counter(str(row["source_status"]) for row in rows)
    usable = {
        row["pmid"]
        for row in rows
        if str(row["source_status"]) in DATA_AVAILABLE_STATUSES
    }
    selected_fetch = {row["pmid"] for row in rows if row["action"] == "fetch"}
    selected_refresh = {
        row["pmid"] for row in rows if row["action"] == "refresh_replay"
    }
    selected_manual = {
        row["pmid"] for row in rows if row["action"] == "manual_or_blocked"
    }
    missing_variant_supplement = {
        row["pmid"] for row in rows if row["missing_variant_supplement"]
    }
    zero_variant_usable = {
        row["pmid"]
        for row in rows
        if row["zero_variant_qc"] and row["source_status"] in DATA_AVAILABLE_STATUSES
    }
    return {
        "pmids": total,
        "by_action": dict(sorted(by_action.items())),
        "by_route": dict(sorted(by_route.items())),
        "by_source_status": dict(sorted(by_source_status.items())),
        "pmid_coverage": {
            "usable_fulltext_current": _coverage_item(len(usable), total),
            "selected_for_fetch": _coverage_item(len(selected_fetch), total),
            "selected_for_source_refresh": _coverage_item(len(selected_refresh), total),
            "selected_for_manual_or_blocked": _coverage_item(
                len(selected_manual), total
            ),
            "missing_variant_supplement": _coverage_item(
                len(missing_variant_supplement), total
            ),
            "zero_variant_usable_fulltext": _coverage_item(
                len(zero_variant_usable), total
            ),
        },
        "notes": {
            "selected_for_fetch": (
                "Run PMIDs without usable run-local full text; feed fetch_input.csv "
                "to scripts/fetch_paywalled.py."
            ),
            "selected_for_source_refresh": (
                "Run PMIDs with usable source text that should be replayed through "
                "scripts/refresh_run_db.py --stage-extractions."
            ),
            "selected_for_manual_or_blocked": (
                "Run PMIDs without usable source text that are known blocked "
                "publisher cases rather than normal fetch queue items."
            ),
        },
    }


def write_fetch_input(rows: list[dict[str, Any]], path: Path) -> None:
    fetch_rows = [
        {"PMID": row["pmid"], "DOI": row["doi"], "route": row["route"]}
        for row in rows
        if row["action"] == "fetch"
    ]
    write_csv_rows(fetch_rows, FETCH_INPUT_FIELDS, path)


def write_source_override(rows: list[dict[str, Any]], path: Path) -> None:
    refresh_rows = [
        {
            "gene": row["gene"],
            "pmid": row["pmid"],
            "action": "refresh_replay",
            "route": row["route"],
            "available_context_path": row["available_context_path"],
            "available_context_bytes": row["available_context_bytes"],
            "notes": row["notes"],
        }
        for row in rows
        if row["action"] == "refresh_replay" and row["available_context_path"]
    ]
    write_csv_rows(refresh_rows, SOURCE_OVERRIDE_FIELDS, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True, help="Gene symbol")
    parser.add_argument("--run-dir", required=True, type=Path, help="GVF run directory")
    parser.add_argument("--harvest-dir", type=Path, default=None)
    parser.add_argument("--extraction-dir", type=Path, default=None)
    parser.add_argument("--paywalled-csv", type=Path, default=None)
    parser.add_argument("--source-completeness", type=Path, default=None)
    parser.add_argument(
        "--include-discovery-pmids",
        action="store_true",
        help=(
            "Also include raw *_pmids.txt discovery lists. Off by default so "
            "source recovery queues focus on PMIDs that reached harvest or "
            "extraction signals."
        ),
    )
    parser.add_argument("--out", required=True, type=Path, help="Worklist CSV")
    parser.add_argument("--summary-out", required=True, type=Path)
    parser.add_argument("--fetch-input-out", type=Path, default=None)
    parser.add_argument("--source-override-out", type=Path, default=None)
    args = parser.parse_args()

    run_dir = repo_path(args.run_dir).expanduser()
    rows, summary = build_audit(
        gene=args.gene,
        run_dir=run_dir,
        harvest_dir=repo_path(args.harvest_dir).expanduser()
        if args.harvest_dir
        else None,
        extraction_dir=repo_path(args.extraction_dir).expanduser()
        if args.extraction_dir
        else None,
        paywalled_csv=repo_path(args.paywalled_csv).expanduser()
        if args.paywalled_csv
        else None,
        source_completeness_path=repo_path(args.source_completeness).expanduser()
        if args.source_completeness
        else None,
        include_discovery_pmids=args.include_discovery_pmids,
    )
    write_csv_rows(rows, FIELDS, repo_path(args.out).expanduser())
    summary_path = repo_path(args.summary_out).expanduser()
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if args.fetch_input_out:
        write_fetch_input(rows, repo_path(args.fetch_input_out).expanduser())
    if args.source_override_out:
        write_source_override(rows, repo_path(args.source_override_out).expanduser())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
