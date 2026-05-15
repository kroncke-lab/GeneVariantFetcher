#!/usr/bin/env python3
"""Export the RYR2 gold-standard counts from the OneDrive RyR2 Variant Database xlsx.

Unlike scripts/build_gold_standard_from_varbrowser.py (which dumps Variant_Browser
Azure SQL tables for KCNH2/KCNQ1/SCN5A), RYR2 has no Variant_Browser per-PMID
clinical table. The literature-curated counts live in a hand-maintained xlsx on
OneDrive:

    /path/to/RYR2_20241129.xlsx

This script:
  1. Reads the xlsx (one sheet, 6109 rows x 34 cols).
  2. Dumps every row verbatim to raw/RYR2_clinical_counts.csv with a 1-based
     source_row_ordinal (audit key).
  3. Splits rows into literature (digit-bearing PMID) vs gnomAD-only (no PMID);
     stashes the gnomAD subset in raw/RYR2_gnomad_rows.csv.
  4. Normalizes literature rows into the same long-format shape as the other
     three genes -> normalized/RYR2_clinical_counts.csv (with extra RYR2-specific
     phenotype columns appended for completeness).
  5. Emits normalized/RYR2_recall_input.csv (variant, pmid, carriers, affected,
     unaffected) using the same convention as the other genes.
  6. Appends the RYR2 long-format rows to normalized/clinical_counts_long.csv
     (preserving existing KCNH2/KCNQ1/SCN5A rows).
  7. Updates qc/summary.json and manifest.json in-place.

Usage:
    .venv/bin/python scripts/build_ryr2_gold_standard_from_xlsx.py \
        --input "/path/to/RYR2_20241129.xlsx" \
        --out   /path/to/gene_variant_fetcher_gold_standard

Idempotent: re-running overwrites all RYR2 outputs and re-syncs the unified
long-format CSV / manifest / summary.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

GVF_REPO = Path(__file__).resolve().parents[1]
if str(GVF_REPO) not in sys.path:
    sys.path.insert(0, str(GVF_REPO))


def _load_variant_normalizer():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "_gvf_variant_normalizer",
        GVF_REPO / "utils" / "variant_normalizer.py",
    )
    if spec is None or spec.loader is None:
        raise ImportError("could not build importlib spec for variant_normalizer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


try:
    _vn_mod = _load_variant_normalizer()
    normalize_variant = _vn_mod.normalize_variant
except Exception as exc:
    print(
        f"WARNING: could not import variant_normalizer ({exc!r}); "
        "variant_normalized will fall back to variant_raw.",
        file=sys.stderr,
    )

    def normalize_variant(variant: str, gene_symbol: str = "RYR2") -> str:
        return variant or ""


SOURCE_TABLE_LABEL = "RYR2_20241129.xlsx"


# Long-format column contract (matches scripts/build_gold_standard_from_varbrowser.py)
NORMALIZED_COLUMNS = [
    "gene",
    "source_table",
    "source_row_ordinal",
    "source_row_id",
    "source_type",
    "lit_flag",
    "cohort_flag",
    "pmid",
    "cohort_label",
    "year",
    "variant_raw",
    "variant_normalized",
    "variant_normalize_method",
    "affected_lqt",
    "affected_lqt3",
    "affected_brs1",
    "affected_sqts",
    "affected_other",
    "scd",
    "sids",
    "unaffected",
    "ambiguous",
    "unknown_phenotype",
    "homozygous",
    "compound_hets",
    "het_count",
    "total_stated",
    "carriers_total_computed",
    "affected_total",
    "count_balance_delta",
    "other_diag",
    "quality_flags",
    "crosswalk_match_method",
    "browser_variant",
]

# RYR2-only phenotype columns appended to the per-gene normalized CSV. Not
# carried into clinical_counts_long.csv (which keeps the unified 34-col schema).
RYR2_EXTRA_COLUMNS = [
    "affected_cpvt",
    "affected_arvc",
    "cardiac_arrest",
    "syncope",
    "seizure",
    "vf",
    "vt",
    "neural_phenotype",
    "epilepsy",
    "adhd",
    "developmental_delay",
    "num_other_diag",
    "exon",
    "variant_codon",
]


_TRUE_BLANKS = {"", "-", "—", "?", "n/a", "na", "none", "null", "nan"}


def _parse_int(value: Any, *, col_name: str, flags: List[str]) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if value != value:  # NaN
            return None
        return int(value)
    text = str(value).strip()
    if not text or text.lower() in _TRUE_BLANKS:
        return None
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return int(float(text))
    except ValueError:
        flags.append(f"parse_int_failure:{col_name}")
        return None


def _classify_pmid(
    raw_pmid: Any,
) -> Tuple[str, Optional[str], Optional[str], List[str]]:
    """Return (source_type, pmid_digits_or_None, cohort_label_or_None, extra_flags)."""
    flags: List[str] = []
    if raw_pmid is None:
        return "blank", None, None, flags
    if isinstance(raw_pmid, float) and raw_pmid != raw_pmid:
        return "blank", None, None, flags
    text = str(raw_pmid).strip()
    if not text:
        return "blank", None, None, flags
    # pandas may yield ints or floats for numeric cells
    if isinstance(raw_pmid, (int,)):
        return "pubmed", str(raw_pmid), None, flags
    if isinstance(raw_pmid, float):
        return "pubmed", str(int(raw_pmid)), None, flags
    if text.isdigit():
        return "pubmed", text, None, flags
    lowered = text.lower()
    if "personal communication" in lowered or "personal_communication" in lowered:
        return "personal_communication", None, None, flags
    if "cohort" in lowered:
        return "cohort", None, text, flags
    digits = re.findall(r"\d+", text)
    if digits and len(digits[0]) >= 6:
        flags.append("unparseable_pmid")
        return "pubmed", digits[0], None, flags
    return "personal_communication", None, text, flags


def _normalize_variant(variant: str, flags: List[str]) -> Tuple[str, str]:
    raw = (variant or "").strip()
    if not raw:
        return "", "raw_fallback"
    try:
        norm = normalize_variant(raw, "RYR2")
        if norm:
            return norm, "normalizer"
    except Exception:
        pass
    flags.append("variant_normalize_fallback")
    return raw, "raw_fallback"


def _maybe_year_from_pmid(pmid: Optional[str]) -> Optional[int]:
    # The RYR2 xlsx has no Year column. Leave year null; downstream consumers
    # that need it can resolve via PubMed.
    return None


def _row_ryr2(
    raw: Dict[str, Any], ordinal: int
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Return (long_format_row, extras_dict_for_per_gene_csv)."""
    flags: List[str] = []

    source_type, pmid, cohort_label, pmid_flags = _classify_pmid(raw.get("PMID"))
    flags.extend(pmid_flags)

    variant_raw = raw.get("var")
    if variant_raw is None or (
        isinstance(variant_raw, float) and variant_raw != variant_raw
    ):
        variant_raw = ""
    else:
        variant_raw = str(variant_raw).strip()

    variant_norm, method = _normalize_variant(variant_raw, flags)

    affected = _parse_int(
        raw.get("Affected (mutation positive)"), col_name="Affected", flags=flags
    )
    unaffected = _parse_int(
        raw.get("unaffected (mutation positive)"), col_name="Unaffected", flags=flags
    )
    scd = _parse_int(raw.get("Sudden Cardiac Death"), col_name="SCD", flags=flags)
    cardiac_arrest = _parse_int(raw.get("CA"), col_name="CA", flags=flags)
    arvc = _parse_int(raw.get("ARVC"), col_name="ARVC", flags=flags)
    cpvt = _parse_int(raw.get("CPVT"), col_name="CPVT", flags=flags)
    lqts = _parse_int(raw.get("LQTS"), col_name="LQTS", flags=flags)
    syncope = _parse_int(raw.get("syncope"), col_name="syncope", flags=flags)
    seizure = _parse_int(raw.get("seizure"), col_name="seizure", flags=flags)
    vf = _parse_int(raw.get("VF"), col_name="VF", flags=flags)
    vt = _parse_int(raw.get("VT"), col_name="VT", flags=flags)
    neural = _parse_int(raw.get("Neural Phenotype"), col_name="Neural", flags=flags)
    epilepsy = _parse_int(raw.get("epilepsy"), col_name="epilepsy", flags=flags)
    adhd = _parse_int(raw.get("ADHD"), col_name="ADHD", flags=flags)
    devdelay = _parse_int(
        raw.get("developmental delay"), col_name="DevDelay", flags=flags
    )
    num_other = _parse_int(
        raw.get("Number of Other Diagnoses"), col_name="NumOther", flags=flags
    )

    # carriers = affected + unaffected (per user spec; matches the recall_input
    # convention for the other three genes where carriers = sum of phenotype
    # buckets and affected = sum of disease-affected columns).
    carriers_total_computed = (affected or 0) + (unaffected or 0)
    affected_total = affected or 0

    # Map RYR2 phenotypes onto the unified long-format schema where possible:
    #   LQTS -> affected_lqt (LQT1 in KCNQ1, LQT in KCNH2)
    #   Sudden Cardiac Death -> scd
    # The remaining RYR2-specific phenotypes (CPVT, ARVC, syncope, seizure,
    # VF, VT, Neural, epilepsy, ADHD, dev. delay, CA, num_other) live only in
    # the per-gene RYR2_clinical_counts.csv extras block.

    long_row = {
        "gene": "RYR2",
        "source_table": SOURCE_TABLE_LABEL,
        "source_row_ordinal": ordinal,
        "source_row_id": ordinal,  # xlsx has no natural PK; use ordinal
        "source_type": source_type,
        "lit_flag": None,
        "cohort_flag": None,
        "pmid": pmid,
        "cohort_label": cohort_label,
        "year": _maybe_year_from_pmid(pmid),
        "variant_raw": variant_raw,
        "variant_normalized": variant_norm,
        "variant_normalize_method": method,
        "affected_lqt": lqts,
        "affected_lqt3": None,
        "affected_brs1": None,
        "affected_sqts": None,
        "affected_other": None,
        "scd": scd,
        "sids": None,
        "unaffected": unaffected,
        "ambiguous": None,
        "unknown_phenotype": None,
        "homozygous": None,
        "compound_hets": None,
        "het_count": None,
        "total_stated": None,
        "carriers_total_computed": carriers_total_computed,
        "affected_total": affected_total,
        "count_balance_delta": None,
        "other_diag": _coerce_str(raw.get("proband diagnosis")) or None,
        "quality_flags": "|".join(flags),
        "crosswalk_match_method": None,  # no RYR2 browser table to crosswalk against
        "browser_variant": None,
    }

    extras = {
        "affected_cpvt": cpvt,
        "affected_arvc": arvc,
        "cardiac_arrest": cardiac_arrest,
        "syncope": syncope,
        "seizure": seizure,
        "vf": vf,
        "vt": vt,
        "neural_phenotype": neural,
        "epilepsy": epilepsy,
        "adhd": adhd,
        "developmental_delay": devdelay,
        "num_other_diag": num_other,
        "exon": _coerce_str(raw.get("exon")),
        "variant_codon": _coerce_str(raw.get("Variant Codon")),
    }
    return long_row, extras


def _coerce_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value != value:
        return ""
    return str(value).strip()


def _git_commit() -> Optional[str]:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(GVF_REPO),
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return out.strip() or None
    except Exception:
        return None


def _git_branch() -> Optional[str]:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=str(GVF_REPO),
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return out.strip() or None
    except Exception:
        return None


def _portable_source_path(path: Path) -> str:
    """Keep generated manifests portable across machines."""
    try:
        return str(path.resolve().relative_to(GVF_REPO))
    except ValueError:
        return path.name


def _read_existing_long(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """Read the current clinical_counts_long.csv, drop any RYR2 rows so we
    can rewrite them with fresh data while preserving KCNH2/KCNQ1/SCN5A.
    """
    if not path.is_file():
        return list(NORMALIZED_COLUMNS), []
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or NORMALIZED_COLUMNS)
        rows = [r for r in reader if (r.get("gene") or "").upper() != "RYR2"]
    return fieldnames, rows


def _write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Path to RYR2 xlsx")
    parser.add_argument(
        "--out", required=True, type=Path, help="Gold-standard output root"
    )
    args = parser.parse_args()

    xlsx_path: Path = args.input.expanduser().resolve()
    out_root: Path = args.out.expanduser().resolve()

    if not xlsx_path.is_file():
        raise SystemExit(f"ERROR: input xlsx not found: {xlsx_path}")

    started_at = dt.datetime.now(dt.timezone.utc).isoformat()
    src_mtime = dt.datetime.fromtimestamp(
        xlsx_path.stat().st_mtime, tz=dt.timezone.utc
    ).isoformat()
    src_size = xlsx_path.stat().st_size

    print(f"Reading {xlsx_path}")
    df = pd.read_excel(xlsx_path, sheet_name=0, dtype=object)
    sheet_name = pd.ExcelFile(xlsx_path).sheet_names[0]
    print(f"  sheet={sheet_name!r}  rows={len(df)}  cols={len(df.columns)}")

    raw_dir = out_root / "raw"
    norm_dir = out_root / "normalized"
    qc_dir = out_root / "qc"
    raw_dir.mkdir(parents=True, exist_ok=True)
    norm_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Verbatim raw dump (audit trail) -----------------------------------
    raw_path = raw_dir / "RYR2_clinical_counts.csv"
    raw_cols = ["source_row_ordinal"] + list(df.columns)
    with raw_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(raw_cols)
        for ordinal, (_, row) in enumerate(df.iterrows(), start=1):
            out_row = [ordinal]
            for col in df.columns:
                val = row[col]
                if val is None or (isinstance(val, float) and val != val):
                    out_row.append("")
                else:
                    out_row.append(val)
            writer.writerow(out_row)
    print(f"  raw/{raw_path.name}: {len(df)} rows")

    # ---- 2. Normalize per-row -------------------------------------------------
    long_rows: List[Dict[str, Any]] = []
    per_gene_rows: List[Dict[str, Any]] = []
    gnomad_rows: List[Tuple[int, Dict[str, Any]]] = []
    qc_duplicates_counter: Counter = Counter()
    source_type_counter: Counter = Counter()

    for ordinal, (_, row) in enumerate(df.iterrows(), start=1):
        raw = row.to_dict()
        long_row, extras = _row_ryr2(raw, ordinal)
        source_type_counter[long_row["source_type"]] += 1

        if long_row["source_type"] == "blank":
            gnomad_rows.append((ordinal, raw))
            # also keep them out of normalized output — they're not gold standard
            continue

        long_rows.append(long_row)
        per_gene_row = dict(long_row)
        per_gene_row.update(extras)
        per_gene_rows.append(per_gene_row)

        key = (
            "RYR2",
            long_row["variant_normalized"] or long_row["variant_raw"],
            long_row["pmid"] or "",
        )
        qc_duplicates_counter[key] += 1

    # Annotate duplicates with quality_flag
    for r in long_rows:
        key = ("RYR2", r["variant_normalized"] or r["variant_raw"], r["pmid"] or "")
        if qc_duplicates_counter[key] > 1:
            existing = r.get("quality_flags") or ""
            r["quality_flags"] = "|".join(filter(None, [existing, "duplicate_pair"]))
    for r in per_gene_rows:
        key = ("RYR2", r["variant_normalized"] or r["variant_raw"], r["pmid"] or "")
        if qc_duplicates_counter[key] > 1:
            existing = r.get("quality_flags") or ""
            r["quality_flags"] = "|".join(filter(None, [existing, "duplicate_pair"]))

    # ---- 3. Per-gene normalized CSV ------------------------------------------
    per_gene_path = norm_dir / "RYR2_clinical_counts.csv"
    per_gene_cols = list(NORMALIZED_COLUMNS) + list(RYR2_EXTRA_COLUMNS)
    _write_csv(per_gene_path, per_gene_rows, per_gene_cols)
    print(
        f"  normalized/{per_gene_path.name}: {len(per_gene_rows)} rows ({len(per_gene_cols)} cols)"
    )

    # ---- 4. recall_input CSV --------------------------------------------------
    recall_path = norm_dir / "RYR2_recall_input.csv"
    recall_count = 0
    with recall_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["variant", "pmid", "carriers", "affected", "unaffected"])
        for r in long_rows:
            if r["source_type"] != "pubmed":
                continue
            if not r["variant_normalized"]:
                continue
            writer.writerow(
                [
                    r["variant_normalized"],
                    r["pmid"],
                    int(r["carriers_total_computed"] or 0),
                    int(r["affected_total"] or 0),
                    int(r["unaffected"] or 0),
                ]
            )
            recall_count += 1
    print(f"  normalized/{recall_path.name}: {recall_count} pubmed rows")

    # ---- 5. gnomAD-only raw stash --------------------------------------------
    gnomad_path = raw_dir / "RYR2_gnomad_rows.csv"
    with gnomad_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(raw_cols)
        for ordinal, raw in gnomad_rows:
            out_row = [ordinal]
            for col in df.columns:
                val = raw.get(col)
                if val is None or (isinstance(val, float) and val != val):
                    out_row.append("")
                else:
                    out_row.append(val)
            writer.writerow(out_row)
    print(f"  raw/{gnomad_path.name}: {len(gnomad_rows)} rows (no PMID)")

    # ---- 6. Append RYR2 to clinical_counts_long.csv --------------------------
    long_path = norm_dir / "clinical_counts_long.csv"
    existing_fieldnames, existing_rows = _read_existing_long(long_path)
    # Existing schema is the canonical 34-col contract; ensure ours matches.
    if set(existing_fieldnames) != set(NORMALIZED_COLUMNS):
        # Be defensive — print but don't fail; the writer will use existing
        # column order and DictWriter ignores extras.
        print(
            f"  WARN: clinical_counts_long.csv columns differ from expected. "
            f"existing={len(existing_fieldnames)} expected={len(NORMALIZED_COLUMNS)}",
            file=sys.stderr,
        )
    fieldnames = existing_fieldnames if existing_rows else NORMALIZED_COLUMNS
    combined_rows = existing_rows + long_rows
    _write_csv(long_path, combined_rows, fieldnames)
    print(
        f"  normalized/{long_path.name}: {len(combined_rows)} rows "
        f"({len(existing_rows)} existing + {len(long_rows)} RYR2 appended)"
    )

    # ---- 7. Append duplicates to qc/duplicate_variant_source_pairs.csv -------
    dup_path = qc_dir / "duplicate_variant_source_pairs.csv"
    dup_fieldnames = ["gene", "variant_normalized", "pmid", "row_count"]
    existing_dup_rows: List[Dict[str, str]] = []
    if dup_path.is_file():
        with dup_path.open("r", newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            existing_dup_rows = [
                r for r in reader if (r.get("gene") or "").upper() != "RYR2"
            ]
    new_dup_rows = [
        {"gene": "RYR2", "variant_normalized": k[1], "pmid": k[2], "row_count": n}
        for k, n in qc_duplicates_counter.items()
        if n > 1
    ]
    _write_csv(dup_path, existing_dup_rows + new_dup_rows, dup_fieldnames)
    print(f"  qc/{dup_path.name}: appended {len(new_dup_rows)} RYR2 duplicate pairs")

    # ---- 8. Update qc/summary.json --------------------------------------------
    summary_path = qc_dir / "summary.json"
    summary: Dict[str, Any] = {}
    if summary_path.is_file():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except Exception as exc:
            print(
                f"  WARN: could not parse existing summary.json ({exc!r}); rebuilding minimally",
                file=sys.stderr,
            )
            summary = {}
    summary.setdefault("row_counts", {})
    unique_pmids = sorted({r["pmid"] for r in long_rows if r.get("pmid")})
    summary["row_counts"]["RYR2"] = {
        "normalized_rows": len(long_rows),
        "recall_input_rows": recall_count,
        "source_type_breakdown": dict(source_type_counter),
        "unique_pmid_count": len(unique_pmids),
    }
    # Update overall anomaly counts (only the duplicate-pairs metric is affected
    # by RYR2; unmatched/balance/browser-recon stay as-is because RYR2 has no
    # crosswalk and no stated totals).
    anomaly = summary.setdefault("anomaly_counts", {})
    # Recompute duplicate count from the file we just wrote.
    with dup_path.open("r", newline="", encoding="utf-8") as fh:
        anomaly["duplicate_variant_source_pairs"] = sum(1 for _ in csv.DictReader(fh))
    summary["generated_at_utc_ryr2"] = started_at
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"  qc/summary.json: RYR2 entry added")

    # ---- 9. Update manifest.json ----------------------------------------------
    manifest_path = out_root / "manifest.json"
    manifest: Dict[str, Any] = {}
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception as exc:
            print(
                f"  WARN: could not parse existing manifest.json ({exc!r}); rebuilding minimally",
                file=sys.stderr,
            )
            manifest = {}
    manifest.setdefault("genes", {})
    manifest["genes"]["RYR2"] = {
        "source_kind": "xlsx",
        "source_path": _portable_source_path(xlsx_path),
        "source_filename": xlsx_path.name,
        "source_mtime_utc": src_mtime,
        "source_size_bytes": src_size,
        "sheet_name": sheet_name,
        "raw_row_count_total": int(len(df)),
        "raw_row_count_literature": len(long_rows),
        "raw_row_count_gnomad_only": len(gnomad_rows),
        "normalized_row_count": len(long_rows),
        "recall_input_row_count": recall_count,
        "source_type_breakdown": dict(source_type_counter),
        "unique_pmid_count": len(unique_pmids),
        "browser_table_fqn": None,
        "browser_filter": None,
        "note": (
            "Loaded from a hand-curated OneDrive xlsx (RYR2_20241129.xlsx). "
            "No Variant_Browser per-PMID table exists for RYR2. The xlsx ships "
            "rich per-disease columns (CPVT, ARVC, CA, syncope, seizure, VF, VT, "
            "neural phenotype, epilepsy, ADHD, developmental delay) preserved in "
            "the per-gene RYR2_clinical_counts.csv (suffix columns); only the "
            "subset that fits the unified long-format schema is propagated to "
            "clinical_counts_long.csv. Population-only rows (no PMID, ~5094 "
            "gnomAD entries) are stashed in raw/RYR2_gnomad_rows.csv."
        ),
        "exporter": "scripts/build_ryr2_gold_standard_from_xlsx.py",
        "exporter_run_at_utc": started_at,
        "exporter_git_commit": _git_commit(),
        "exporter_git_branch": _git_branch(),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"  manifest.json: RYR2 entry updated")

    print()
    print("=== summary ===")
    print(f"  raw rows:                {len(df)}")
    print(f"  literature (pubmed) rows: {source_type_counter.get('pubmed', 0)}")
    print(
        f"  personal_communication:   {source_type_counter.get('personal_communication', 0)}"
    )
    print(f"  gnomAD-only / blank:      {source_type_counter.get('blank', 0)}")
    print(f"  recall_input rows:        {recall_count}")
    print(f"  unique PMIDs:             {len(unique_pmids)}")
    print(
        f"  duplicate (variant,PMID): {sum(1 for n in qc_duplicates_counter.values() if n > 1)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
