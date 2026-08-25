#!/usr/bin/env python3
"""Score the pipeline against an assistant's completed curation packet.

Closes the cold-start measurement loop: takes the filled
``curation_template.csv`` from ``build_curation_packet.py`` and the gene's DB,
converts the human answers into a normalized gold recall input, and runs the
existing recall scorer (recall + precision + F2) restricted to the curated
PMIDs. The result is a real, defensible recall/F2 number for a gene that had no
gold standard.

Germline is the heritable-carrier target: ``somatic``-labelled rows are excluded
by default (they are not carriers); ``unknown`` is kept. Explicit ``NONE`` rows
are absent from the gold variant list but their PMIDs remain in the exhaustive
paper scope, so pipeline false positives on those papers are penalized. Blank
variants are unfinished work and fail validation rather than masquerading as
negative papers. Family/unknown count roles are excluded from person-level MAE.

Usage:
    python scripts/score_curation_packet.py \
        --filled-csv curation_template_FILLED.csv \
        --db results/BRCA2/<...>/BRCA2.db --gene BRCA2 \
        --out recall_metrics/brca2_gold_50
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable
VALID_STATUSES = {"complete"}
VALID_SCOPES = {"germline", "somatic", "unknown"}
VALID_COUNT_ROLES = {"individual", "family", "unknown", "not_reported"}


class CurationValidationError(ValueError):
    """The filled answer key is incomplete or internally inconsistent."""


def _norm(s: str | None) -> str:
    return (s or "").strip()


def _count_value(value: str, *, pmid: str, field: str) -> str:
    if not value:
        return ""
    try:
        parsed = int(value)
    except ValueError as exc:
        raise CurationValidationError(
            f"PMID {pmid}: {field} must be a non-negative integer or blank"
        ) from exc
    if parsed < 0:
        raise CurationValidationError(
            f"PMID {pmid}: {field} must be a non-negative integer or blank"
        )
    return str(parsed)


def convert(
    filled_csv: Path,
    include_somatic: bool,
    *,
    expected_pmids: set[str] | None = None,
) -> tuple[list[dict], dict]:
    """Validate a completed curation CSV and produce normalized gold rows."""
    rows: list[dict] = []
    pmids: set[str] = set()
    kinds_by_pmid: dict[str, list[str]] = {}
    seen_variant_rows: set[tuple[str, str]] = set()
    dropped_somatic = none_papers = role_excluded_counts = 0
    with filled_csv.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        required = {
            "pmid",
            "variant",
            "curation_status",
            "germline_or_somatic",
            "carriers",
            "carrier_count_role",
            "affected",
            "unaffected",
            "evidence_note",
        }
        missing_columns = required - set(reader.fieldnames or [])
        if missing_columns:
            raise CurationValidationError(
                "missing required columns: " + ", ".join(sorted(missing_columns))
            )
        for line_number, r in enumerate(reader, start=2):
            pmid = _norm(r.get("pmid"))
            if not pmid or pmid.upper() == "EXAMPLE":
                continue
            if not pmid.isdigit():
                raise CurationValidationError(
                    f"line {line_number}: invalid PMID {pmid!r}"
                )
            pmids.add(pmid)
            variant = _norm(r.get("variant"))
            status = _norm(r.get("curation_status")).lower()
            if status not in VALID_STATUSES:
                raise CurationValidationError(
                    f"PMID {pmid}: curation_status must be complete before scoring"
                )
            evidence_note = _norm(r.get("evidence_note"))
            if not evidence_note:
                raise CurationValidationError(f"PMID {pmid}: evidence_note is required")
            if not variant:
                raise CurationValidationError(
                    f"PMID {pmid}: blank variant is unfinished; use NONE explicitly"
                )
            if variant.upper() == "NONE":
                none_scope = _norm(r.get("germline_or_somatic"))
                none_role = _norm(r.get("carrier_count_role")).lower()
                none_counts = {
                    field: _norm(r.get(field))
                    for field in ("carriers", "affected", "unaffected")
                }
                if (
                    none_scope
                    or any(none_counts.values())
                    or none_role
                    not in {
                        "",
                        "not_reported",
                    }
                ):
                    raise CurationValidationError(
                        f"PMID {pmid}: NONE cannot carry scope or count values"
                    )
                kinds_by_pmid.setdefault(pmid, []).append("none")
                continue

            kinds_by_pmid.setdefault(pmid, []).append("variant")
            variant_key = (pmid, variant.casefold())
            if variant_key in seen_variant_rows:
                raise CurationValidationError(
                    f"PMID {pmid}: duplicate variant row {variant!r}"
                )
            seen_variant_rows.add(variant_key)
            scope = _norm(r.get("germline_or_somatic")).lower()
            if scope not in VALID_SCOPES:
                raise CurationValidationError(
                    f"PMID {pmid}: germline_or_somatic must be germline, somatic, or unknown"
                )
            carriers = _count_value(
                _norm(r.get("carriers")), pmid=pmid, field="carriers"
            )
            affected = _count_value(
                _norm(r.get("affected")), pmid=pmid, field="affected"
            )
            unaffected = _count_value(
                _norm(r.get("unaffected")), pmid=pmid, field="unaffected"
            )
            count_role = _norm(r.get("carrier_count_role")).lower()
            if count_role not in VALID_COUNT_ROLES:
                raise CurationValidationError(
                    f"PMID {pmid}: invalid carrier_count_role {count_role!r}"
                )
            if carriers and count_role == "not_reported":
                raise CurationValidationError(
                    f"PMID {pmid}: a supplied carrier count cannot be not_reported"
                )
            if not carriers and count_role == "individual":
                raise CurationValidationError(
                    f"PMID {pmid}: individual carrier_count_role requires carriers"
                )
            if count_role in {"family", "unknown"}:
                if affected or unaffected:
                    raise CurationValidationError(
                        f"PMID {pmid}: affected/unaffected require an individual count role"
                    )
                if carriers:
                    role_excluded_counts += 1
                carriers = ""

            if scope == "somatic" and not include_somatic:
                dropped_somatic += 1
                continue
            rows.append(
                {
                    "variant": variant,
                    "pmid": pmid,
                    "carriers": carriers,
                    "affected": affected,
                    "unaffected": unaffected,
                }
            )

    for pmid, kinds in kinds_by_pmid.items():
        if "none" in kinds and ("variant" in kinds or len(kinds) > 1):
            raise CurationValidationError(
                f"PMID {pmid}: NONE cannot be mixed with variant rows or repeated"
            )
        if kinds == ["none"]:
            none_papers += 1
    if expected_pmids is not None:
        missing_pmids = expected_pmids - pmids
        unexpected_pmids = pmids - expected_pmids
        if missing_pmids or unexpected_pmids:
            parts = []
            if missing_pmids:
                parts.append("missing=" + ",".join(sorted(missing_pmids)))
            if unexpected_pmids:
                parts.append("unexpected=" + ",".join(sorted(unexpected_pmids)))
            raise CurationValidationError(
                "curated PMID roster mismatch: " + "; ".join(parts)
            )
    summary = {
        "curated_pmids": len(pmids),
        "pmids": sorted(pmids),
        "gold_variant_rows": len(rows),
        "no_variant_papers": none_papers,
        "dropped_somatic": dropped_somatic,
        "nonindividual_counts_excluded_from_mae": role_excluded_counts,
        "include_somatic": include_somatic,
        "paper_exhaustive": True,
    }
    return rows, summary


def write_recall_input(rows: list[dict], gene: str, gold_dir: Path) -> Path:
    norm = gold_dir / "normalized"
    norm.mkdir(parents=True, exist_ok=True)
    path = norm / f"{gene}_recall_input.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f, fieldnames=["variant", "pmid", "carriers", "affected", "unaffected"]
        )
        w.writeheader()
        w.writerows(rows)
    return path


def write_pmid_scope(pmids: set[str], gene: str, gold_dir: Path) -> Path:
    """Persist the exhaustive paper scope, including explicit NONE papers."""
    path = gold_dir / f"{gene}_curated_pmids.txt"
    path.write_text("\n".join(sorted(pmids)) + "\n", encoding="utf-8")
    return path


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--filled-csv", required=True, type=Path)
    ap.add_argument("--db", required=True, type=Path)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument(
        "--packet-meta",
        type=Path,
        help="packet_meta.json; defaults to the filled CSV's directory.",
    )
    ap.add_argument(
        "--include-somatic",
        action="store_true",
        help="Count somatic-labelled variants too (default: germline/unknown only).",
    )
    args = ap.parse_args()

    gene = args.gene.upper()
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    filled_csv = args.filled_csv.expanduser().resolve()
    packet_meta = (
        args.packet_meta.expanduser().resolve()
        if args.packet_meta
        else filled_csv.parent / "packet_meta.json"
    )
    if not packet_meta.exists():
        print(f"error: packet metadata not found: {packet_meta}", file=sys.stderr)
        return 2
    meta = json.loads(packet_meta.read_text(encoding="utf-8"))
    expected_pmids = {str(pmid).strip() for pmid in meta.get("pmids") or []}
    if not expected_pmids:
        print("error: packet metadata has no PMID roster", file=sys.stderr)
        return 2
    try:
        rows, summary = convert(
            filled_csv,
            args.include_somatic,
            expected_pmids=expected_pmids,
        )
    except CurationValidationError as exc:
        print(f"error: curation validation failed: {exc}", file=sys.stderr)
        return 2
    gold_dir = out / "_gold"
    write_recall_input(rows, gene, gold_dir)
    pmid_scope = write_pmid_scope(expected_pmids, gene, gold_dir)
    (out / "curation_coverage.json").write_text(json.dumps(summary, indent=2))
    print("Curation coverage:", json.dumps(summary))

    cmd = [
        PY,
        str(REPO / "scripts" / "run_recall_suite.py"),
        "--score",
        "--genes",
        gene,
        "--db",
        f"{gene}={args.db.expanduser().resolve()}",
        "--gold-dir",
        str(gold_dir),
        "--gold-pmids",
        f"{gene}={pmid_scope}",
        "--gold-paper-exhaustive",
        "--review-gold-sync",
        "off",
        "--review-gold-tier",
        "all",
        "--skip-disagreement-artifacts",
        "--outdir",
        str(out),
    ]
    print("Scoring:", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=str(REPO))
    if proc.returncode != 0:
        return proc.returncode

    summ_path = out / "summary.json"
    if summ_path.exists():
        s = json.loads(summ_path.read_text())
        rec = s.get("aggregate_recall", {})
        uv = rec.get("unique_variants", {})
        f1 = s.get("aggregate_f1") or {}
        f2 = s.get("aggregate_fbeta") or {}
        print("\n=== COLD-START RESULT (vs curated gold) ===")
        print(
            f"  unique-variant recall: {uv.get('matched')}/{uv.get('gold')} = {uv.get('recall')}"
        )
        prec = (s.get("aggregate_precision") or {}).get("precision_vs_gold_pmids")
        print(f"  paper-exhaustive precision: {prec}")
        if f1:
            print(f"  F1: {f1.get('fbeta_vs_gold_pmids')}")
        if f2:
            print(f"  F2: {f2.get('fbeta_vs_gold_pmids')}")
        print(f"  full report: {summ_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
