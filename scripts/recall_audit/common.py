#!/usr/bin/env python3
"""Shared helpers for recall-audit scripts.

These scripts intentionally use the live scored artifact layout rather than a
separate fixture schema:

  validation_runs/.../<run>/dbs/{GENE}.db
  validation_runs/.../<run>/recall_score/{GENE}/discrepancies.csv
  gene_variant_fetcher_gold_standard/normalized/{GENE}_recall_input.csv
"""

from __future__ import annotations

import csv
import json
import os
import re
import sqlite3
import sys
from pathlib import Path
from typing import Any, Iterable, Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

STATUS_DOC = REPO_ROOT / "docs" / "RECALL_STATUS.md"
RUN_DIR_ENV = "GVF_RECALL_RUN_DIR"
DEFAULT_GOLD_DIR = REPO_ROOT / "gene_variant_fetcher_gold_standard" / "normalized"
DEFAULT_RESULTS_DIR = REPO_ROOT / "results"
DEFAULT_CONTEXT_MIN_BYTES = 5_000
SUPPLEMENT_REF_RE = re.compile(
    r"\b(?:Supplement(?:al|ary)?|Supporting)\s+"
    r"(?:Table|Tables|Data|File)\s*([A-Za-z0-9]+)",
    re.IGNORECASE,
)
SUPPLEMENT_BODY_RE = re.compile(
    r"(?im)^\s*(?:#{1,6}\s*)?(?:\*\*)?\s*"
    r"(?:Supplement(?:al|ary)?|Supporting)\s+"
    r"(?:Table|Data|File)\s*([A-Za-z0-9]+)\b",
)
VARIANT_SUPPLEMENT_TERMS = (
    "mutation",
    "mutations",
    "variant",
    "variants",
    "genetic",
)
SUPPLEMENT_POINTER_TERMS = (
    "listed",
    "shown",
    "provided",
    "available",
    "summarized",
    "summarised",
    "contained",
    "found",
)


def repo_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def current_status_run_dir() -> Path | None:
    if not STATUS_DOC.exists():
        return None
    text = STATUS_DOC.read_text(encoding="utf-8")
    match = re.search(r"`([^`]+/recall_score/summary\.json)`", text)
    if not match:
        return None
    return repo_path(match.group(1)).parent.parent


def normalize_pmid(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text:
        return ""
    try:
        return str(int(float(text)))
    except ValueError:
        return text


def parse_number(value: Any) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_int(value: Any) -> Optional[int]:
    number = parse_number(value)
    if number is None:
        return None
    return int(number)


def parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _supplement_ref_label(value: str) -> str:
    return re.sub(r"\W+", "", value).upper()


def unresolved_variant_supplement_refs(text: str, gene: str) -> list[str]:
    """Return referenced supplement-table labels that are not present as bodies.

    This intentionally looks for explicit prose pointers such as "All SCN5A
    mutations are listed in Supplemental Table 2" and then checks whether a
    same-number supplemental table body appears as a heading/caption. It avoids
    treating ordinary main-text mentions of supplemental figures or statistics
    as missing variant data.
    """
    if not text:
        return []
    body_labels = {
        _supplement_ref_label(match.group(1))
        for match in SUPPLEMENT_BODY_RE.finditer(text)
    }
    unresolved: set[str] = set()
    for match in SUPPLEMENT_REF_RE.finditer(text):
        label = _supplement_ref_label(match.group(1))
        if not label or label in body_labels:
            continue
        start = max(0, match.start() - 180)
        end = min(len(text), match.end() + 180)
        snippet = text[start:end]
        relative_match = match.start() - start
        prefix = snippet[:relative_match]
        suffix = snippet[relative_match:]
        left_boundary = max(prefix.rfind("."), prefix.rfind("\n"))
        right_period = suffix.find(".")
        right_newline = suffix.find("\n")
        right_candidates = [idx for idx in (right_period, right_newline) if idx != -1]
        right_boundary = (
            relative_match + min(right_candidates) if right_candidates else len(snippet)
        )
        sentence = snippet[left_boundary + 1 : right_boundary].lower()
        if (
            gene
            and gene.lower() not in sentence
            and not any(term in sentence for term in VARIANT_SUPPLEMENT_TERMS)
        ):
            continue
        if not any(term in sentence for term in SUPPLEMENT_POINTER_TERMS):
            continue
        unresolved.add(label)
    return sorted(unresolved)


def source_has_unresolved_variant_supplement_refs(path: Path | None, gene: str) -> bool:
    if path is None or not path.exists():
        return False
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    return bool(unresolved_variant_supplement_refs(text, gene))


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(
    rows: Iterable[dict[str, Any]], fieldnames: list[str], out_path: Optional[Path]
) -> None:
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        handle = out_path.open("w", newline="", encoding="utf-8")
        close = True
    else:
        handle = sys.stdout
        close = False
    try:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    finally:
        if close:
            handle.close()


def resolve_run_dir(value: str | Path | None) -> Path:
    if value:
        run_dir = repo_path(value)
    elif os.environ.get(RUN_DIR_ENV):
        run_dir = repo_path(os.environ[RUN_DIR_ENV])
    else:
        run_dir = current_status_run_dir()
        if run_dir is None:
            raise SystemExit(
                "No recall run directory configured. Pass --run-dir, set "
                f"{RUN_DIR_ENV}, or add a recall_score/summary.json path to "
                f"{STATUS_DOC.relative_to(REPO_ROOT)}."
            )

    if not run_dir.exists():
        raise SystemExit(
            f"Recall run directory does not exist: {run_dir}. Pass --run-dir "
            f"or set {RUN_DIR_ENV} to an available scored run."
        )
    return run_dir


def resolve_recall_score(
    recall_score: str | Path | None, run_dir: str | Path | None = None
) -> Path:
    if recall_score:
        return repo_path(recall_score)
    return resolve_run_dir(run_dir) / "recall_score"


def resolve_gene_db(
    gene: str, db: str | Path | None = None, run_dir: str | Path | None = None
) -> Path:
    if db:
        return repo_path(db)
    return resolve_run_dir(run_dir) / "dbs" / f"{gene.upper()}.db"


def resolve_gold_path(gene: str, gold_dir: str | Path | None = None) -> Path:
    root = repo_path(gold_dir) if gold_dir else DEFAULT_GOLD_DIR
    return root / f"{gene.upper()}_recall_input.csv"


def load_gold_rows(
    gene: str, pmid: str | None = None, gold_dir: str | Path | None = None
) -> list[dict[str, str]]:
    path = resolve_gold_path(gene, gold_dir)
    rows = read_csv_rows(path)
    if pmid is None:
        return rows
    wanted = normalize_pmid(pmid)
    return [row for row in rows if normalize_pmid(row.get("pmid")) == wanted]


def display_variant(row: dict[str, Any]) -> str:
    for key in (
        "protein_notation",
        "cdna_notation",
        "genomic_position",
        "variant",
        "sqlite_variant_raw",
        "excel_variant_raw",
    ):
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def canonical_variant(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if not text:
        return ""
    try:
        from cli.compare_variants import normalize_variant, to_canonical_form

        return to_canonical_form(text) or normalize_variant(text) or text
    except Exception:
        return re.sub(r"\s+", "", text).upper()


def query_pipeline_rows(db_path: Path, pmid: str | None = None) -> list[dict[str, Any]]:
    where = ""
    params: tuple[Any, ...] = ()
    if pmid is not None:
        where = "WHERE pd.pmid = ?"
        params = (normalize_pmid(pmid),)
    sql = f"""
        SELECT
            v.variant_id,
            v.gene_symbol,
            v.cdna_notation,
            v.protein_notation,
            v.genomic_position,
            v.clinical_significance,
            v.evidence_level,
            pd.penetrance_id,
            pd.pmid,
            pd.total_carriers_observed,
            pd.affected_count,
            pd.unaffected_count,
            pd.uncertain_count,
            pd.penetrance_percentage,
            vp.source_location,
            vp.additional_notes,
            p.title
        FROM penetrance_data pd
        JOIN variants v ON v.variant_id = pd.variant_id
        LEFT JOIN variant_papers vp
            ON vp.variant_id = pd.variant_id AND vp.pmid = pd.pmid
        LEFT JOIN papers p ON p.pmid = pd.pmid
        {where}
        ORDER BY pd.pmid, v.protein_notation, v.cdna_notation, pd.penetrance_id
    """
    with sqlite3.connect(db_path) as conn:
        conn.row_factory = sqlite3.Row
        return [dict(row) for row in conn.execute(sql, params)]


def format_table(headers: list[str], rows: list[list[Any]], max_width: int = 48) -> str:
    prepared = [[_clip(cell, max_width) for cell in row] for row in rows]
    widths = [
        max(len(str(header)), *(len(str(row[i])) for row in prepared))
        if prepared
        else len(str(header))
        for i, header in enumerate(headers)
    ]
    header_line = " | ".join(
        str(header).ljust(widths[i]) for i, header in enumerate(headers)
    )
    sep = "-+-".join("-" * width for width in widths)
    body = [
        " | ".join(str(cell).ljust(widths[i]) for i, cell in enumerate(row))
        for row in prepared
    ]
    return "\n".join([header_line, sep, *body])


def _clip(value: Any, max_width: int) -> str:
    text = "" if value is None else str(value)
    if len(text) <= max_width:
        return text
    return text[: max_width - 3] + "..."


def find_full_contexts(
    results_dir: Path, gene: str, pmid: str, *, global_search: bool = False
) -> list[Path]:
    wanted = f"{normalize_pmid(pmid)}_FULL_CONTEXT.md"
    gene_root = results_dir / gene.upper()
    candidates: list[Path] = []
    if gene_root.exists():
        candidates.extend(gene_root.glob(f"*/pmc_fulltext/{wanted}"))
        candidates.extend(gene_root.glob(f"**/{wanted}"))
    if global_search and not candidates and results_dir.exists():
        candidates.extend(results_dir.glob(f"**/{wanted}"))
    candidates = [p for p in candidates if p.exists()]
    unique = sorted(
        set(candidates), key=lambda p: (p.stat().st_size, str(p)), reverse=True
    )
    return unique


def context_status(
    path: Optional[Path], min_bytes: int = DEFAULT_CONTEXT_MIN_BYTES
) -> str:
    if path is None:
        return "not_attempted"
    size = path.stat().st_size
    sample = path.read_text(encoding="utf-8", errors="ignore")[:20_000].lower()
    if (
        "# abstract-only fallback" in sample
        or "full text could not be retrieved" in sample
        or "this document contains only the pubmed abstract" in sample
    ):
        return "abstract_only"
    if size < min_bytes:
        return "paywall_stub"
    has_supp = "supplement" in sample or "supplementary" in sample
    has_body = any(
        marker in sample
        for marker in ("abstract", "introduction", "methods", "results")
    )
    if has_supp and not has_body:
        return "recovered_supplement_only"
    if any(
        marker in sample
        for marker in ("browser", "sciencedirect", "science direct", "authenticated")
    ):
        return "recovered_browser"
    return "recovered_pmc"


def json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n",
        encoding="utf-8",
    )
