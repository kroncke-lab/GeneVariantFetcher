#!/usr/bin/env python3
"""Export the GeneVariantFetcher gold-standard counts from Variant_Browser's Azure SQL.

Run from inside the Variant_Browser venv with its .env sourced so Django + ODBC
credentials are available:

    cd /path/to/Variant_Browser
    set -a; source .env; set +a
    venv/bin/python /path/to/GeneVariantFetcher/scripts/build_gold_standard_from_varbrowser.py \
        --out /path/to/GeneVariantFetcher/gene_variant_fetcher_gold_standard

The script:
  1. Bootstraps Django (DJANGO_SETTINGS_MODULE=variantbrowser.settings).
  2. Dumps each source table verbatim into <out>/raw/<table>.csv.
  3. Re-iterates rows with tolerant int parsing, classifies source_type, normalizes
     variants, and writes <out>/normalized/{<GENE>_clinical_counts.csv, clinical_counts_long.csv,
     variant_crosswalk.csv, <GENE>_recall_input.csv}.
  4. Emits QC artifacts under <out>/qc/ and the run manifest at <out>/manifest.json.

Idempotent: re-running overwrites every output in place.
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
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

GVF_REPO = Path(__file__).resolve().parents[1]
if str(GVF_REPO) not in sys.path:
    sys.path.insert(0, str(GVF_REPO))


# Locate the Variant_Browser project root so Django can import `variantbrowser.settings`.
# Strategy 1: env var (set by caller if needed).
# Strategy 2: derive from the active venv (venv/bin/python -> venv -> Variant_Browser root).
# Strategy 3: cwd if it contains `manage.py` and `variantbrowser/settings.py`.
def _locate_variantbrowser_root() -> Path:
    env_hint = os.environ.get("VARIANTBROWSER_ROOT")
    candidates: List[Path] = []
    if env_hint:
        candidates.append(Path(env_hint).resolve())
    try:
        candidates.append(Path(sys.executable).resolve().parents[2])
    except IndexError:
        pass
    candidates.append(Path.cwd().resolve())
    for cand in candidates:
        if (cand / "manage.py").is_file() and (
            cand / "variantbrowser" / "settings.py"
        ).is_file():
            return cand
    raise SystemExit(
        "ERROR: could not locate the Variant_Browser project root.\n"
        "Tried: " + ", ".join(str(c) for c in candidates) + "\n"
        "Hint: set VARIANTBROWSER_ROOT, or run from inside the Variant_Browser checkout."
    )


VARIANTBROWSER_ROOT = _locate_variantbrowser_root()
if str(VARIANTBROWSER_ROOT) not in sys.path:
    sys.path.insert(0, str(VARIANTBROWSER_ROOT))

import django  # noqa: E402

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "variantbrowser.settings")
django.setup()

from django.db import connection  # noqa: E402


# Load utils/variant_normalizer.py directly (bypass utils/__init__.py whose
# transitive imports pull in `requests`, which isn't installed in the
# Variant_Browser venv where this script runs).
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
    normalize_variant = _vn_mod.normalize_variant  # type: ignore[assignment]
except Exception as exc:  # pragma: no cover - import guard
    print(
        f"WARNING: could not import variant_normalizer ({exc!r}); "
        "variant_normalized will fall back to variant_raw.",
        file=sys.stderr,
    )

    def normalize_variant(variant: str, gene_symbol: str = "KCNH2") -> str:  # type: ignore[misc]
        return variant or ""


# ------------------------------------------------------------------- helpers --


def _ensure_dirs(out_root: Path) -> None:
    for sub in ("raw", "normalized", "qc", "schema"):
        (out_root / sub).mkdir(parents=True, exist_ok=True)


def _fetch_table_ordered(
    table_fqn: str, order_by_sql: str
) -> Tuple[List[str], List[List[Any]], List[Dict[str, Any]]]:
    """Single canonical fetch for a table.

    Returns (column_names, row_lists, row_dicts) where row_lists preserves
    SELECT order column-for-column (for the raw CSV dump) and row_dicts is the
    same rows reshaped for the normalization mappers. Same iteration order
    both ways guarantees that the ordinal index in normalized output matches
    the row position in the raw CSV.
    """
    sql = f"SELECT * FROM {table_fqn} ORDER BY {order_by_sql}"
    with connection.cursor() as cur:
        cur.execute(sql)
        cols = [c[0] for c in cur.description]
        raw_rows: List[List[Any]] = []
        dict_rows: List[Dict[str, Any]] = []
        while True:
            batch = cur.fetchmany(1000)
            if not batch:
                break
            for row in batch:
                raw_rows.append(list(row))
                dict_rows.append(dict(zip(cols, row)))
    return cols, raw_rows, dict_rows


def _write_raw_csv(raw_filename: Path, cols: List[str], rows: List[List[Any]]) -> int:
    with raw_filename.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(cols)
        for row in rows:
            writer.writerow(["" if v is None else v for v in row])
    return len(rows)


def _independent_count(table_fqn: str) -> int:
    with connection.cursor() as cur:
        cur.execute(f"SELECT COUNT(*) FROM {table_fqn}")
        return int(cur.fetchone()[0])


_TRUE_BLANKS = {"", "-", "—", "?", "n/a", "na", "none", "null"}


def _parse_int(value: Any, *, col_name: str, flags: List[str]) -> Optional[int]:
    """Coerce a possibly-text value to int. Records parse failures in flags."""
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
    if text == "" or text.lower() in _TRUE_BLANKS:
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


def _maybe_year(value: Any) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, int):
        return value if 1900 < value < 2100 else None
    text = str(value).strip()
    if not text:
        return None
    m = re.search(r"(19|20)\d{2}", text)
    if m:
        try:
            return int(m.group(0))
        except ValueError:
            return None
    return None


def _classify_pmid_text(raw_pmid: Any) -> Tuple[str, Optional[str], Optional[str]]:
    """Returns (source_type, pmid_digits_or_None, cohort_label_or_None)."""
    if raw_pmid is None:
        return "blank", None, None
    text = str(raw_pmid).strip()
    if not text:
        return "blank", None, None
    if text.isdigit():
        return "pubmed", text, None
    lowered = text.lower()
    if "personal communication" in lowered:
        return "personal_communication", None, None
    if "cohort" in lowered:
        return "cohort", None, text
    digits = re.findall(r"\d+", text)
    if digits and len(digits[0]) >= 6:
        # PubMed IDs are 6+ digits; pick the first plausible run, but flag the
        # caller so they know it wasn't a pristine numeric cell.
        return "pubmed", digits[0], None
    return "cohort", None, text


def _sum_non_null(*values: Optional[int]) -> int:
    return sum(v for v in values if v is not None)


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


def _portable_path(path: Path) -> str:
    """Return a manifest-safe path without machine-specific absolute prefixes."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(GVF_REPO))
    except ValueError:
        return resolved.name


# ---------------------------------------------------------- normalization ---

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

GENE_TO_NORMALIZER_SYMBOL = {"KCNH2": "KCNH2", "KCNQ1": "KCNQ1", "SCN5A": "SCN5A"}


def _normalize(variant: str, gene: str, flags: List[str]) -> Tuple[str, str]:
    raw = (variant or "").strip()
    if not raw:
        return "", "raw_fallback"
    try:
        norm = normalize_variant(raw, GENE_TO_NORMALIZER_SYMBOL.get(gene, gene))
        if norm:
            if norm != raw:
                pass  # silent; we record the method below
            return norm, "normalizer"
    except Exception:
        pass
    flags.append("variant_normalize_fallback")
    return raw, "raw_fallback"


# ----------------------------------------------------------- KCNH2 mapping ---


def _row_kcnh2(raw: Dict[str, Any]) -> Dict[str, Any]:
    flags: List[str] = []
    source_type, pmid, cohort_label = _classify_pmid_text(raw.get("PMID"))
    variant_raw = (raw.get("Variant") or "").strip()
    variant_norm, method = _normalize(variant_raw, "KCNH2", flags)
    if method == "raw_fallback":
        pass  # _normalize already flagged

    lqt = _parse_int(raw.get("LQT"), col_name="LQT", flags=flags)
    car = _parse_int(raw.get("CAR"), col_name="CAR", flags=flags)
    amb = _parse_int(raw.get("AMB"), col_name="AMB", flags=flags)
    una = _parse_int(raw.get("UNA"), col_name="UNA", flags=flags)
    amb_una = _parse_int(raw.get("AMB+UNA"), col_name="AMB+UNA", flags=flags)
    unk = _parse_int(raw.get("UNK"), col_name="UNK", flags=flags)
    hom = _parse_int(raw.get("HOM"), col_name="HOM", flags=flags)
    scd = _parse_int(raw.get("SCD"), col_name="SCD", flags=flags)

    # Prefer the computed AMB+UNA column as `unaffected` if present, else AMB+UNA components.
    if amb_una is not None:
        unaffected = amb_una
    else:
        unaffected = _sum_non_null(amb, una)

    # KCNH2 carrier total is the sum of all phenotype buckets we have.
    # CAR is the heterozygous "any carrier" count and is sometimes the cohort total.
    components = [lqt, unaffected, unk, hom]
    carriers_total_computed = _sum_non_null(*components)
    affected_total = _sum_non_null(lqt)

    return {
        "gene": "KCNH2",
        "source_table": "[dbo].[KCNH2_clinical]",
        "source_row_id": raw.get("clinical_papers_id"),
        "source_type": source_type,
        "lit_flag": None,
        "cohort_flag": None,
        "pmid": pmid,
        "cohort_label": cohort_label,
        "year": _maybe_year(raw.get("Year")),
        "variant_raw": variant_raw,
        "variant_normalized": variant_norm,
        "variant_normalize_method": method,
        "affected_lqt": lqt,
        "affected_lqt3": None,
        "affected_brs1": None,
        "affected_sqts": None,
        "affected_other": None,
        "scd": scd,
        "sids": None,
        "unaffected": unaffected,
        "ambiguous": amb,
        "unknown_phenotype": unk,
        "homozygous": hom,
        "compound_hets": None,
        "het_count": car,
        "total_stated": None,
        "carriers_total_computed": carriers_total_computed,
        "affected_total": affected_total,
        "count_balance_delta": None,
        "other_diag": raw.get("Other DIAG"),
        "quality_flags": "|".join(flags),
    }


# ----------------------------------------------------------- KCNQ1 mapping ---


def _bool_flag(raw_value: Any) -> Optional[bool]:
    if raw_value is None:
        return None
    if isinstance(raw_value, bool):
        return raw_value
    try:
        return bool(int(raw_value))
    except (ValueError, TypeError):
        return None


def _row_kcnq1(raw: Dict[str, Any]) -> Dict[str, Any]:
    flags: List[str] = []
    lit_flag = _bool_flag(raw.get("Lit_Y_N"))
    cohort_flag = _bool_flag(raw.get("Cohort_Y_N"))
    pmid_int = raw.get("PMID")
    if pmid_int is None:
        source_type = "cohort" if cohort_flag else "blank"
        pmid = None
        cohort_label = "KCNQ1 cohort (no PMID)" if source_type == "cohort" else None
    else:
        pmid = str(int(pmid_int))
        if lit_flag:
            source_type = "pubmed"
            cohort_label = None
        elif cohort_flag:
            source_type = "cohort"
            cohort_label = f"PMID {pmid} cohort"
        else:
            # No explicit flag, but PMID is set — treat as pubmed.
            source_type = "pubmed"
            cohort_label = None

    variant_raw = (raw.get("var") or "").strip()
    variant_norm, method = _normalize(variant_raw, "KCNQ1", flags)

    lqt1 = _parse_int(raw.get("LQT1"), col_name="LQT1", flags=flags)
    sqts = _parse_int(raw.get("SQTS"), col_name="SQTS", flags=flags)
    asymptomatic = _parse_int(
        raw.get("Asymptomatic"), col_name="Asymptomatic", flags=flags
    )
    het = _parse_int(raw.get("Het"), col_name="Het", flags=flags)
    homozygous = _parse_int(
        raw.get("Homozygous_Carriers"), col_name="Homozygous_Carriers", flags=flags
    )
    compound_hets = _parse_int(
        raw.get("Heterozygous_JLNS"), col_name="Heterozygous_JLNS", flags=flags
    )
    ambiguous = _parse_int(
        raw.get("Ambiguous_Phenotype_Hets"),
        col_name="Ambiguous_Phenotype_Hets",
        flags=flags,
    )
    unknown = _parse_int(
        raw.get("Unknown_Phenotype_Hets"),
        col_name="Unknown_Phenotype_Hets",
        flags=flags,
    )
    scd = _parse_int(
        raw.get("Sudden_Cardiac_Death"), col_name="Sudden_Cardiac_Death", flags=flags
    )
    sids = _parse_int(raw.get("SIDS"), col_name="SIDS", flags=flags)

    carriers_total_computed = _sum_non_null(
        lqt1, sqts, asymptomatic, ambiguous, unknown, homozygous, compound_hets
    )
    affected_total = _sum_non_null(lqt1, sqts)

    return {
        "gene": "KCNQ1",
        "source_table": "[dbo].[KCNQ1_clinical_v10]",
        "source_row_id": raw.get("Columna_1"),
        "source_type": source_type,
        "lit_flag": lit_flag,
        "cohort_flag": cohort_flag,
        "pmid": pmid,
        "cohort_label": cohort_label,
        "year": _maybe_year(raw.get("Year")),
        "variant_raw": variant_raw,
        "variant_normalized": variant_norm,
        "variant_normalize_method": method,
        "affected_lqt": lqt1,
        "affected_lqt3": None,
        "affected_brs1": None,
        "affected_sqts": sqts,
        "affected_other": None,
        "scd": scd,
        "sids": sids,
        "unaffected": asymptomatic,
        "ambiguous": ambiguous,
        "unknown_phenotype": unknown,
        "homozygous": homozygous,
        "compound_hets": compound_hets,
        "het_count": het,
        "total_stated": None,
        "carriers_total_computed": carriers_total_computed,
        "affected_total": affected_total,
        "count_balance_delta": None,
        "other_diag": raw.get("Other_Diagnoses_of_Interest"),
        "quality_flags": "|".join(flags),
    }


# ----------------------------------------------------------- SCN5A mapping ---


def _row_scn5a(raw: Dict[str, Any]) -> Dict[str, Any]:
    flags: List[str] = []
    pmid_value = raw.get("PMID")
    if pmid_value is None:
        source_type = "blank"
        pmid = None
    else:
        pmid = str(int(pmid_value))
        source_type = "pubmed"

    variant_raw = (raw.get("Variant") or "").strip()
    variant_norm, method = _normalize(variant_raw, "SCN5A", flags)

    lqt3 = _parse_int(raw.get("LQT3"), col_name="LQT3", flags=flags)
    brs1 = _parse_int(raw.get("BrS"), col_name="BrS", flags=flags)
    other = _parse_int(raw.get("Other"), col_name="Other", flags=flags)
    normal = _parse_int(raw.get("Normal"), col_name="Normal", flags=flags)
    total_stated = _parse_int(raw.get("Total"), col_name="Total", flags=flags)

    carriers_total_computed = _sum_non_null(lqt3, brs1, other, normal)
    affected_total = _sum_non_null(lqt3, brs1, other)
    count_balance_delta = None
    if total_stated is not None:
        count_balance_delta = total_stated - carriers_total_computed

    return {
        "gene": "SCN5A",
        "source_table": "[SCN5A].[SCN5A_papers]",
        "source_row_id": raw.get("id"),
        "source_type": source_type,
        "lit_flag": None,
        "cohort_flag": None,
        "pmid": pmid,
        "cohort_label": None,
        "year": _maybe_year(raw.get("Year")),
        "variant_raw": variant_raw,
        "variant_normalized": variant_norm,
        "variant_normalize_method": method,
        "affected_lqt": None,
        "affected_lqt3": lqt3,
        "affected_brs1": brs1,
        "affected_sqts": None,
        "affected_other": other,
        "scd": None,
        "sids": None,
        "unaffected": normal,
        "ambiguous": None,
        "unknown_phenotype": None,
        "homozygous": None,
        "compound_hets": None,
        "het_count": None,
        "total_stated": total_stated,
        "carriers_total_computed": carriers_total_computed,
        "affected_total": affected_total,
        "count_balance_delta": count_balance_delta,
        "other_diag": raw.get("Other_Disease"),
        "quality_flags": "|".join(flags),
    }


# ----------------------------------------------------------- gene configs ---

GENES = [
    {
        "gene": "KCNH2",
        "table_fqn": "[dbo].[KCNH2_clinical]",
        "raw_filename": "KCNH2_clinical.csv",
        "order_by_sql": "clinical_papers_id",
        "row_mapper": _row_kcnh2,
        "browser_table_fqn": "[dbo].[all_vars_annotated]",
        "browser_filter_sql": "isoform = 'A'",
        "browser_var_column": "var",
        "browser_hgvsc_column": "HGVSc",
        "browser_total_column": "total_carriers_adj",
    },
    {
        "gene": "KCNQ1",
        "table_fqn": "[dbo].[KCNQ1_clinical_v10]",
        "raw_filename": "KCNQ1_clinical_v10.csv",
        # Columna_1 is not unique and has NULLs; var disambiguates remaining ties.
        "order_by_sql": "ISNULL([Columna_1], 9999999), [var]",
        "row_mapper": _row_kcnq1,
        "browser_table_fqn": "[dbo].[KCNQ1_browser_main]",
        "browser_filter_sql": None,
        "browser_var_column": "var",
        "browser_hgvsc_column": "HGVSc",
        "browser_total_column": "total_carriers",
    },
    {
        "gene": "SCN5A",
        "table_fqn": "[SCN5A].[SCN5A_papers]",
        "raw_filename": "SCN5A_papers.csv",
        "order_by_sql": "[id]",
        "row_mapper": _row_scn5a,
        "browser_table_fqn": "[dbo].[scn5a_dataset]",
        "browser_filter_sql": None,
        "browser_var_column": "var",
        "browser_hgvsc_column": "HGVSc",
        "browser_total_column": "total_carriers",
    },
]


# ----------------------------------------------------------- crosswalk ---


def _load_browser_variants(
    gene_cfg: Dict[str, Any],
) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """Returns (browser_rows, normalized_lookup) for a gene's browser table."""
    fqn = gene_cfg["browser_table_fqn"]
    where = (
        f"WHERE {gene_cfg['browser_filter_sql']}"
        if gene_cfg["browser_filter_sql"]
        else ""
    )
    var_col = gene_cfg["browser_var_column"]
    hgvsc_col = gene_cfg["browser_hgvsc_column"]
    total_col = gene_cfg["browser_total_column"]
    sql = (
        f"SELECT [{var_col}] AS var, [{hgvsc_col}] AS hgvsc, "
        f"[{total_col}] AS total_carriers FROM {fqn} {where}"
    )
    with connection.cursor() as cur:
        cur.execute(sql)
        cols = [c[0] for c in cur.description]
        rows = [dict(zip(cols, r)) for r in cur.fetchall()]
    lookup: Dict[str, Dict[str, Any]] = {}
    for r in rows:
        for key in (r.get("var"), r.get("hgvsc")):
            if not key:
                continue
            try:
                norm = normalize_variant(str(key), gene_cfg["gene"])
            except Exception:
                norm = str(key)
            if norm and norm not in lookup:
                lookup[norm] = r
    return rows, lookup


def _try_crosswalk(
    normalized_variant: str,
    raw_variant: str,
    lookup: Dict[str, Dict[str, Any]],
) -> Tuple[Optional[str], str]:
    """Returns (browser_variant, match_method)."""
    if not normalized_variant:
        return None, "unmatched"
    match = lookup.get(normalized_variant)
    if match:
        return match.get("var"), "normalized"
    match = lookup.get(raw_variant)
    if match:
        return match.get("var"), "exact"
    upper = normalized_variant.upper()
    match = lookup.get(upper)
    if match:
        return match.get("var"), "normalized"
    return None, "unmatched"


# ----------------------------------------------------------- main run ---


def _write_csv(
    path: Path, rows: Iterable[Dict[str, Any]], fieldnames: List[str]
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            n += 1
    return n


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path, help="Output root directory")
    args = parser.parse_args()

    out_root: Path = args.out.resolve()
    _ensure_dirs(out_root)

    started_at = dt.datetime.now(dt.timezone.utc).isoformat()

    print(f"Bootstrapping export to {out_root}")
    print(
        f"  Source DB: {os.environ.get('DJANGO_DATABASE_SERVER', '?')}::"
        f"{os.environ.get('DJANGO_DATABASE_NAME', '?')}"
    )

    manifest_genes: Dict[str, Any] = {}

    all_normalized_rows: List[Dict[str, Any]] = []
    qc_unmatched: List[Dict[str, Any]] = []
    qc_count_flags: List[Dict[str, Any]] = []
    qc_browser_recon: List[Dict[str, Any]] = []
    qc_duplicates: List[Dict[str, Any]] = []

    for gene_cfg in GENES:
        gene = gene_cfg["gene"]
        table_fqn = gene_cfg["table_fqn"]
        print(f"\n=== {gene} :: {table_fqn} ===")

        raw_path = out_root / "raw" / gene_cfg["raw_filename"]
        _cols, raw_lists, raw_rows = _fetch_table_ordered(
            table_fqn, gene_cfg["order_by_sql"]
        )
        raw_count = _write_raw_csv(raw_path, _cols, raw_lists)
        independent_count = _independent_count(table_fqn)
        print(
            f"  raw/{raw_path.name}: {raw_count} rows (independent SELECT COUNT(*)={independent_count})"
        )
        if raw_count != independent_count:
            print(
                f"  WARN: row-count mismatch! raw={raw_count} count={independent_count}",
                file=sys.stderr,
            )
        print(
            f"  Loaded {len(raw_rows)} rows for normalization (ordered by {gene_cfg['order_by_sql']})"
        )

        browser_rows, browser_lookup = _load_browser_variants(gene_cfg)
        print(
            f"  Loaded {len(browser_rows)} browser rows; {len(browser_lookup)} normalized keys"
        )

        gene_normalized: List[Dict[str, Any]] = []
        per_gene_pair_counter: Counter = Counter()
        lit_carriers_per_variant: Dict[str, Dict[str, Any]] = defaultdict(
            lambda: {"sum_carriers": 0, "rows": 0}
        )

        for ordinal, raw_row in enumerate(raw_rows, start=1):
            row = gene_cfg["row_mapper"](raw_row)
            # 1-based ordinal in the SELECT * iteration order. Combined with `gene`
            # this is a guaranteed-unique trace key, even when the natural PK
            # column has duplicates or NULLs (KCNQ1.Columna_1 has both).
            row["source_row_ordinal"] = ordinal
            # Crosswalk
            browser_variant, match_method = _try_crosswalk(
                row["variant_normalized"], row["variant_raw"], browser_lookup
            )
            row["crosswalk_match_method"] = match_method
            row["browser_variant"] = browser_variant
            if match_method == "unmatched":
                existing_flags = row.get("quality_flags") or ""
                row["quality_flags"] = "|".join(
                    filter(None, [existing_flags, "crosswalk_unmatched"])
                )
                qc_unmatched.append(
                    {
                        "gene": gene,
                        "source_table": row["source_table"],
                        "source_row_id": row["source_row_id"],
                        "variant_raw": row["variant_raw"],
                        "variant_normalized": row["variant_normalized"],
                        "pmid": row["pmid"],
                    }
                )

            # Count balance check (SCN5A is the only one with a stated total)
            if (
                row["count_balance_delta"] is not None
                and row["count_balance_delta"] != 0
            ):
                qc_count_flags.append(
                    {
                        "gene": gene,
                        "source_table": row["source_table"],
                        "source_row_id": row["source_row_id"],
                        "variant_raw": row["variant_raw"],
                        "pmid": row["pmid"],
                        "total_stated": row["total_stated"],
                        "carriers_total_computed": row["carriers_total_computed"],
                        "count_balance_delta": row["count_balance_delta"],
                    }
                )

            # Duplicate tracking
            key = (
                gene,
                row["variant_normalized"] or row["variant_raw"],
                row["pmid"] or "",
            )
            per_gene_pair_counter[key] += 1

            # Accumulate browser reconciliation only for literature rows.
            if row["source_type"] == "pubmed" and browser_variant:
                bucket = lit_carriers_per_variant[browser_variant]
                bucket["sum_carriers"] += int(row["carriers_total_computed"] or 0)
                bucket["rows"] += 1

            gene_normalized.append(row)

        # Per-gene normalized CSV
        per_gene_path = out_root / "normalized" / f"{gene}_clinical_counts.csv"
        _write_csv(per_gene_path, gene_normalized, NORMALIZED_COLUMNS)
        print(f"  normalized/{per_gene_path.name}: {len(gene_normalized)} rows")

        # Per-gene recall_input CSV (pubmed only)
        recall_path = out_root / "normalized" / f"{gene}_recall_input.csv"
        recall_rows = 0
        with recall_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["variant", "pmid", "carriers", "affected", "unaffected"])
            for r in gene_normalized:
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
                recall_rows += 1
        print(f"  normalized/{recall_path.name}: {recall_rows} pubmed rows")

        # Crosswalk CSV (gene-scoped — gets merged later)
        crosswalk_path = out_root / "normalized" / f"{gene}_variant_crosswalk.csv"
        with crosswalk_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                [
                    "gene",
                    "clinical_variant_raw",
                    "clinical_variant_normalized",
                    "browser_variant",
                    "match_method",
                    "source_row_id",
                    "pmid",
                ]
            )
            for r in gene_normalized:
                writer.writerow(
                    [
                        gene,
                        r["variant_raw"],
                        r["variant_normalized"],
                        r["browser_variant"] or "",
                        r["crosswalk_match_method"],
                        r["source_row_id"],
                        r["pmid"] or "",
                    ]
                )

        # Duplicate-pairs surface
        for key, n in per_gene_pair_counter.items():
            if n > 1:
                qc_duplicates.append(
                    {
                        "gene": key[0],
                        "variant_normalized": key[1],
                        "pmid": key[2],
                        "row_count": n,
                    }
                )

        # Browser reconciliation
        browser_by_var = {r.get("var"): r for r in browser_rows if r.get("var")}
        for var_name, bucket in lit_carriers_per_variant.items():
            br = browser_by_var.get(var_name) or {}
            browser_total = br.get("total_carriers")
            qc_browser_recon.append(
                {
                    "gene": gene,
                    "browser_variant": var_name,
                    "literature_carrier_sum": bucket["sum_carriers"],
                    "literature_row_count": bucket["rows"],
                    "browser_total_carriers": browser_total,
                    "delta_browser_minus_lit": (
                        (int(browser_total) - bucket["sum_carriers"])
                        if browser_total is not None
                        else None
                    ),
                }
            )

        # Manifest entry
        gene_summary = Counter(r["source_type"] for r in gene_normalized)
        manifest_genes[gene] = {
            "source_table": table_fqn,
            "raw_row_count": raw_count,
            "independent_count_query": independent_count,
            "normalized_row_count": len(gene_normalized),
            "recall_input_row_count": recall_rows,
            "source_type_breakdown": dict(gene_summary),
            "browser_table_fqn": gene_cfg["browser_table_fqn"],
            "browser_filter": gene_cfg["browser_filter_sql"],
            "browser_row_count": len(browser_rows),
            "browser_normalized_keys": len(browser_lookup),
        }

        all_normalized_rows.extend(gene_normalized)

    # RYR2 manifest stub
    manifest_genes["RYR2"] = {
        "source_table": None,
        "note": (
            "No per-PMID clinical paper table exists in Variant_Browser. "
            "Only [RyR2].[sub_tmp_AM_REVEL_againstALL_2_adj] (aggregate variant table) is present. "
            "Skipped: no rows to emit."
        ),
    }

    # Unified long-form CSV
    long_path = out_root / "normalized" / "clinical_counts_long.csv"
    _write_csv(long_path, all_normalized_rows, NORMALIZED_COLUMNS)
    print(f"\nnormalized/{long_path.name}: {len(all_normalized_rows)} rows")

    # Unified crosswalk CSV (concat of per-gene)
    crosswalk_all = out_root / "normalized" / "variant_crosswalk.csv"
    with crosswalk_all.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "gene",
                "clinical_variant_raw",
                "clinical_variant_normalized",
                "browser_variant",
                "match_method",
                "source_row_id",
                "pmid",
            ]
        )
        for r in all_normalized_rows:
            writer.writerow(
                [
                    r["gene"],
                    r["variant_raw"],
                    r["variant_normalized"],
                    r["browser_variant"] or "",
                    r["crosswalk_match_method"],
                    r["source_row_id"],
                    r["pmid"] or "",
                ]
            )

    # QC outputs
    qc_dir = out_root / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(
        qc_dir / "unmatched_variants.csv",
        qc_unmatched,
        [
            "gene",
            "source_table",
            "source_row_id",
            "variant_raw",
            "variant_normalized",
            "pmid",
        ],
    )
    _write_csv(
        qc_dir / "count_balance_flags.csv",
        qc_count_flags,
        [
            "gene",
            "source_table",
            "source_row_id",
            "variant_raw",
            "pmid",
            "total_stated",
            "carriers_total_computed",
            "count_balance_delta",
        ],
    )
    _write_csv(
        qc_dir / "duplicate_variant_source_pairs.csv",
        qc_duplicates,
        ["gene", "variant_normalized", "pmid", "row_count"],
    )
    _write_csv(
        qc_dir / "browser_total_reconciliation.csv",
        qc_browser_recon,
        [
            "gene",
            "browser_variant",
            "literature_carrier_sum",
            "literature_row_count",
            "browser_total_carriers",
            "delta_browser_minus_lit",
        ],
    )

    # Headline summary.json
    summary = {
        "generated_at_utc": started_at,
        "row_counts": {
            gene: {
                "normalized_rows": manifest_genes[gene].get("normalized_row_count", 0),
                "recall_input_rows": manifest_genes[gene].get(
                    "recall_input_row_count", 0
                ),
                "source_type_breakdown": manifest_genes[gene].get(
                    "source_type_breakdown", {}
                ),
            }
            for gene in ("KCNH2", "KCNQ1", "SCN5A")
        },
        "anomaly_counts": {
            "unmatched_variants": len(qc_unmatched),
            "count_balance_flags": len(qc_count_flags),
            "duplicate_variant_source_pairs": len(qc_duplicates),
            "browser_total_reconciliation_rows": len(qc_browser_recon),
        },
        "unique_pmids_by_gene": {
            gene: sorted(
                {
                    r["pmid"]
                    for r in all_normalized_rows
                    if r["gene"] == gene and r["pmid"]
                }
            )
            for gene in ("KCNH2", "KCNQ1", "SCN5A")
        },
    }
    # Replace the lists with counts (the long lists bloat summary.json).
    for gene in ("KCNH2", "KCNQ1", "SCN5A"):
        summary["row_counts"][gene]["unique_pmid_count"] = len(
            summary["unique_pmids_by_gene"][gene]
        )
    del summary["unique_pmids_by_gene"]

    (qc_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"qc/summary.json: {json.dumps(summary['row_counts'], indent=2)}")

    # Manifest
    manifest = {
        "generated_at_utc": started_at,
        "generator_script": str(Path(__file__).resolve().relative_to(GVF_REPO)),
        "gvf_repo_path": ".",
        "gvf_git_commit": _git_commit(),
        "gvf_git_branch": _git_branch(),
        "source_db_server": os.environ.get("DJANGO_DATABASE_SERVER"),
        "source_db_name": os.environ.get("DJANGO_DATABASE_NAME"),
        "source_db_driver": os.environ.get("DJANGO_DATABASE_DRIVER"),
        "kcnh2_browser_isoform_filter": "isoform='A'",
        "genes": manifest_genes,
        "output_root": _portable_path(out_root),
    }
    (out_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"manifest.json written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
