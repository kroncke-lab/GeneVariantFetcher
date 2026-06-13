#!/usr/bin/env python3
"""
Variant Comparison Tool for GeneVariantFetcher

Compares a manually curated Excel sheet of variant carriers against
an automated literature-curation SQLite database produced by GeneVariantFetcher.

Usage:
    python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db

For fuzzy variant matching:
    python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db \
        --variant_match_mode fuzzy --fuzzy_threshold 0.85

Author: GeneVariantFetcher Team
"""

import argparse
import json
import logging
import re
import sqlite3
import unicodedata
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Set, Tuple

import pandas as pd
from utils.source_layers import (
    count_reasons,
    junk_notation_reason,
    normalize_source_layer,
    source_layer_sql_case,
)

# Try to import rapidfuzz for better fuzzy matching; fallback to difflib
try:
    from rapidfuzz import fuzz as rapidfuzz_fuzz

    RAPIDFUZZ_AVAILABLE = True
except ImportError:
    RAPIDFUZZ_AVAILABLE = False
    import difflib

# Optional YAML support for mapping files
try:
    import yaml

    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

# Configure logging
logger = logging.getLogger(__name__)


# =============================================================================
# COLUMN SYNONYM DETECTION
# =============================================================================

COLUMN_SYNONYMS: Dict[str, List[str]] = {
    "pmid": [
        "pmid",
        "pubmed_id",
        "pubmed",
        "pm_id",
        "pub_id",
        "pubmedid",
        "pmid_id",
        "paper_id",
        "article_id",
        "reference_id",
    ],
    "variant": [
        "hgvs",
        "hgvs_p",
        "hgvs_c",
        "protein_change",
        "cdna_change",
        "variant",
        "aa_change",
        "mutation",
        "amino_acid",
        "protein_variant",
        "cdna_variant",
        "nucleotide_change",
        "dna_change",
        "protein_mutation",
        "variant_name",
        "variant_id",
        "hgvs_protein",
        "hgvs_cdna",
        "p_notation",
        "c_notation",
        "protein_notation",
        "cdna_notation",
    ],
    "rsid": ["rsid", "rs_id", "dbsnp", "rs", "rs_number", "dbsnp_id", "snp_id"],
    "carriers_total": [
        "carriers",
        "n_carriers",
        "carrier_count",
        "total_carriers",
        "n_total",
        "subjects",
        "probands",
        "individuals",
        "n_individuals",
        "total",
        "n",
        "count",
        "sample_size",
        "n_subjects",
        "carrier_n",
        "total_n",
        "num_carriers",
        "number_carriers",
        "car",  # LQTS-specific: carrier count
    ],
    "affected_count": [
        "affected",
        "cases",
        "n_affected",
        "symptomatic",
        "patients",
        "n_cases",
        "case_count",
        "affected_count",
        "n_symptomatic",
        "diseased",
        "n_diseased",
        "num_affected",
        "number_affected",
        "lqt",  # LQTS-specific: Long QT syndrome affected
    ],
    "unaffected_count": [
        "unaffected",
        "controls",
        "n_unaffected",
        "asymptomatic",
        "n_controls",
        "control_count",
        "unaffected_count",
        "n_asymptomatic",
        "healthy",
        "n_healthy",
        "num_unaffected",
        "number_unaffected",
        "amb+una",
        "ambuna",  # LQTS-specific: ambiguous + unaffected combined
    ],
    "phenotype": [
        "phenotype",
        "condition",
        "disease",
        "diagnosis",
        "clinical",
        "clinical_phenotype",
        "phenotype_details",
        "disorder",
        "syndrome",
    ],
    "notes": [
        "notes",
        "comments",
        "remarks",
        "additional_info",
        "description",
        "free_text",
        "annotation",
        "additional_notes",
    ],
}

# LQTS-specific affected phenotype columns to sum for total affected count
AFFECTED_PHENOTYPE_COLUMNS: List[str] = [
    "lqt",
    "hom",
    "smpt",
    "asd",
    "scd",
    "drg",
    "stqs",
]


def normalize_column_name(col: Any) -> str:
    """Normalize column name for matching (lowercase, strip, remove underscores/spaces)."""
    if col is None:
        return ""
    col_str = str(col)
    return re.sub(r"[\s_\-]+", "", col_str.lower().strip())


def normalize_pmid(pmid: Any) -> str:
    """
    Normalize a PMID value to a consistent string format.

    Handles:
    - Float values like 16470702.0 -> "16470702"
    - Integer values -> "16470702"
    - String values -> strip and convert
    - NaN/None -> ""
    """
    if pmid is None or (isinstance(pmid, float) and pd.isna(pmid)):
        return ""

    # Try to convert to int first (handles floats like 16470702.0)
    try:
        pmid_int = int(float(pmid))
        return str(pmid_int)
    except (ValueError, TypeError):
        # Fall back to string conversion
        return str(pmid).strip()


def detect_columns(
    df: pd.DataFrame, mapping: Optional[Dict[str, str]] = None
) -> Dict[str, Any]:
    """
    Detect column mappings from DataFrame using synonym matching.

    Args:
        df: Input DataFrame
        mapping: Optional explicit column mapping from user

    Returns:
        Dict mapping semantic field names to actual column names.
        Special key 'affected_phenotype_columns' contains list of detected
        phenotype-specific affected columns to sum.
    """
    detected: Dict[str, Any] = {field: None for field in COLUMN_SYNONYMS.keys()}
    detected["affected_phenotype_columns"] = []  # List of columns to sum for affected

    # If explicit mapping provided, use it first
    if mapping:
        for field, col_name in mapping.items():
            if col_name in df.columns:
                detected[field] = col_name
                logger.info(f"Using explicit mapping: {field} -> {col_name}")
            else:
                logger.warning(f"Explicit mapping column not found: {col_name}")

    # Auto-detect remaining fields
    normalized_cols = {normalize_column_name(c): c for c in df.columns}

    for field, synonyms in COLUMN_SYNONYMS.items():
        if detected[field] is not None:
            continue  # Already mapped

        for synonym in synonyms:
            norm_syn = normalize_column_name(synonym)
            if norm_syn in normalized_cols:
                detected[field] = normalized_cols[norm_syn]
                logger.info(f"Auto-detected column: {field} -> {detected[field]}")
                break

        # Also try partial matching for compound names
        if detected[field] is None:
            for norm_col, orig_col in normalized_cols.items():
                for synonym in synonyms:
                    norm_syn = normalize_column_name(synonym)
                    # Only allow the synonym to be contained in a longer
                    # column name. The reverse direction mis-maps affected as
                    # unaffected, and short synonyms like "rs" match carriers.
                    if len(norm_syn) >= 4 and norm_syn in norm_col:
                        detected[field] = orig_col
                        logger.info(f"Partial match column: {field} -> {orig_col}")
                        break
                if detected[field] is not None:
                    break

    # Detect LQTS-specific affected phenotype columns to sum
    for pheno_col in AFFECTED_PHENOTYPE_COLUMNS:
        norm_pheno = normalize_column_name(pheno_col)
        if norm_pheno in normalized_cols:
            actual_col = normalized_cols[norm_pheno]
            detected["affected_phenotype_columns"].append(actual_col)
            logger.info(f"Detected affected phenotype column: {actual_col}")

    if detected["affected_phenotype_columns"]:
        logger.info(
            f"Will sum {len(detected['affected_phenotype_columns'])} affected phenotype columns: {detected['affected_phenotype_columns']}"
        )

    return detected


# =============================================================================
# SQLITE INTROSPECTION
# =============================================================================


@dataclass
class TableInfo:
    """Information about a database table."""

    name: str
    columns: List[str]
    has_pmid: bool = False
    has_variant: bool = False
    has_counts: bool = False
    pmid_column: Optional[str] = None
    variant_columns: List[str] = field(default_factory=list)
    count_columns: Dict[str, str] = field(default_factory=dict)


def introspect_sqlite(conn: sqlite3.Connection) -> Dict[str, TableInfo]:
    """
    Introspect SQLite database to discover schema.

    Args:
        conn: SQLite connection

    Returns:
        Dict of table name -> TableInfo
    """
    cursor = conn.cursor()

    # Get all tables
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [row[0] for row in cursor.fetchall()]

    logger.info(f"Discovered tables: {tables}")

    pmid_patterns = ["pmid", "pubmed", "paper_id", "article_id"]
    variant_patterns = [
        "protein_notation",
        "cdna_notation",
        "hgvs",
        "variant",
        "mutation",
        "genomic_position",
        "protein_change",
    ]
    count_patterns = {
        "carriers_total": [
            "total_carriers",
            "carriers",
            "total_carriers_observed",
            "n_total",
        ],
        "affected_count": ["affected_count", "affected", "n_affected", "cases"],
        "unaffected_count": [
            "unaffected_count",
            "unaffected",
            "n_unaffected",
            "controls",
        ],
        "uncertain_count": ["uncertain_count", "uncertain", "n_uncertain"],
    }

    table_info = {}

    for table in tables:
        cursor.execute(f"PRAGMA table_info({table})")
        columns = [row[1] for row in cursor.fetchall()]

        info = TableInfo(name=table, columns=columns)

        # Check for PMID column
        for col in columns:
            col_lower = col.lower()
            for pattern in pmid_patterns:
                if pattern in col_lower:
                    info.has_pmid = True
                    info.pmid_column = col
                    break

        # Check for variant columns
        for col in columns:
            col_lower = col.lower()
            for pattern in variant_patterns:
                if pattern in col_lower:
                    info.has_variant = True
                    info.variant_columns.append(col)
                    break

        # Check for count columns
        for col in columns:
            col_lower = col.lower()
            for count_type, patterns in count_patterns.items():
                for pattern in patterns:
                    if pattern in col_lower:
                        info.has_counts = True
                        info.count_columns[count_type] = col
                        break

        table_info[table] = info
        logger.debug(
            f"Table {table}: pmid={info.has_pmid}, variant={info.has_variant}, counts={info.has_counts}"
        )

    return table_info


def find_best_data_source(table_info: Dict[str, TableInfo]) -> Tuple[str, TableInfo]:
    """
    Find the best table for extracting variant-paper-count data.

    Priority:
    1. union_all — UNION of variant_papers + penetrance_data + individual_records
       joined to variants. Preserves every (pmid, variant) link from any source
       (variant_papers alone covers ~397 PMIDs / 2.9k links that penetrance_data
       misses). Counts come from penetrance_data when present, else from
       aggregated individual_records, else NULL.
    2. penetrance_data (legacy single-table)
    3. Combined pmid+variant+counts table
    4. individual_records aggregation
    5. variant_papers + variants join

    Args:
        table_info: Dict of TableInfo from introspection

    Returns:
        Tuple of (strategy_name, primary_table_info)
    """
    # Preferred: UNION across all three sources when variants table exists
    if "variants" in table_info and (
        "variant_papers" in table_info
        or "penetrance_data" in table_info
        or "individual_records" in table_info
    ):
        logger.info(
            "Using union_all strategy: variant_papers + penetrance_data + individual_records"
        )
        # Use whichever of the three exists as the "primary" reference
        for candidate in ("variant_papers", "penetrance_data", "individual_records"):
            if candidate in table_info:
                return ("union_all", table_info[candidate])

    # Legacy fallbacks (kept for non-GVF databases)
    if "penetrance_data" in table_info:
        info = table_info["penetrance_data"]
        if info.has_pmid and info.has_counts:
            logger.info("Using penetrance_data table as primary source")
            return ("penetrance_data", info)

    for name, info in table_info.items():
        if info.has_pmid and info.has_variant and info.has_counts:
            logger.info(f"Using table {name} with pmid+variant+counts")
            return ("combined", info)

    if "individual_records" in table_info:
        info = table_info["individual_records"]
        if info.has_pmid:
            logger.info("Using individual_records table (will aggregate counts)")
            return ("individual_records", info)

    if "variant_papers" in table_info and "variants" in table_info:
        logger.info("Using variant_papers + variants join strategy")
        return ("variant_papers_join", table_info["variant_papers"])

    raise ValueError(
        "Cannot find suitable tables for comparison.\n"
        f"Available tables: {list(table_info.keys())}\n"
        "Expected: penetrance_data, individual_records, or variant_papers with variants\n"
        "Consider using --mapping to specify column mappings."
    )


def extract_sqlite_data(
    conn: sqlite3.Connection, table_info: Dict[str, TableInfo]
) -> pd.DataFrame:
    """
    Extract variant data from SQLite database.

    Args:
        conn: SQLite connection
        table_info: Schema information from introspection

    Returns:
        DataFrame with columns: pmid, variant, carriers_total, affected_count, unaffected_count
    """
    strategy, primary_table = find_best_data_source(table_info)

    if strategy == "union_all":
        # Union every (pmid, variant_id) link from variant_papers, penetrance_data,
        # and individual_records. Counts: prefer penetrance_data > aggregated
        # individual_records > NULL. This recovers links that exist in
        # variant_papers but never made it into a penetrance_data row.
        union_parts = []
        if "variant_papers" in table_info:
            union_parts.append("SELECT pmid, variant_id FROM variant_papers")
        if "penetrance_data" in table_info:
            union_parts.append("SELECT pmid, variant_id FROM penetrance_data")
        if "individual_records" in table_info:
            union_parts.append("SELECT pmid, variant_id FROM individual_records")
        union_sql = "\n  UNION\n  ".join(union_parts)

        has_pd = "penetrance_data" in table_info
        has_ir = "individual_records" in table_info
        has_vp = "variant_papers" in table_info
        has_vp_layer = has_vp and "source_layer" in table_info["variant_papers"].columns
        has_vp_source_location = (
            has_vp and "source_location" in table_info["variant_papers"].columns
        )
        has_vp_additional_notes = (
            has_vp and "additional_notes" in table_info["variant_papers"].columns
        )

        pd_join = (
            "LEFT JOIN penetrance_data pd "
            "ON pd.pmid = al.pmid AND pd.variant_id = al.variant_id"
            if has_pd
            else ""
        )
        ir_join = (
            "LEFT JOIN ir_agg ir ON ir.pmid = al.pmid AND ir.variant_id = al.variant_id"
            if has_ir
            else ""
        )
        ir_cte = (
            """,
        ir_agg AS (
            SELECT
                pmid,
                variant_id,
                COUNT(*) AS carriers_total,
                SUM(CASE WHEN affected_status='affected' THEN 1 ELSE 0 END) AS affected_count,
                SUM(CASE WHEN affected_status='unaffected' THEN 1 ELSE 0 END) AS unaffected_count,
                SUM(CASE WHEN affected_status='uncertain' THEN 1 ELSE 0 END) AS uncertain_count
            FROM individual_records
            GROUP BY pmid, variant_id
        )"""
            if has_ir
            else ""
        )
        vp_layer_expr = source_layer_sql_case(
            "source_location" if has_vp_source_location else "NULL",
            "source_layer" if has_vp_layer else None,
            additional_notes_expr=(
                "additional_notes" if has_vp_additional_notes else None
            ),
        )
        vp_layer_cte = (
            f""",
        vp_layer AS (
            SELECT
                pmid,
                variant_id,
                GROUP_CONCAT(DISTINCT {vp_layer_expr}) AS source_layer
            FROM variant_papers
            GROUP BY pmid, variant_id
        )"""
            if has_vp
            else ""
        )
        vp_join = (
            "LEFT JOIN vp_layer vp ON vp.pmid = al.pmid AND vp.variant_id = al.variant_id"
            if has_vp
            else ""
        )
        layer_select_expr = (
            "COALESCE(vp.source_layer, 'llm_text')" if has_vp else "'llm_text'"
        )

        carriers_expr = []
        affected_expr = []
        unaffected_expr = []
        uncertain_expr = []
        if has_pd:
            carriers_expr.append("pd.total_carriers_observed")
            affected_expr.append("pd.affected_count")
            unaffected_expr.append("pd.unaffected_count")
            uncertain_expr.append("pd.uncertain_count")
        if has_ir:
            carriers_expr.append("ir.carriers_total")
            affected_expr.append("ir.affected_count")
            unaffected_expr.append("ir.unaffected_count")
            uncertain_expr.append("ir.uncertain_count")
        carriers_expr.append("NULL")
        affected_expr.append("NULL")
        unaffected_expr.append("NULL")
        uncertain_expr.append("NULL")

        def coalesce_sql(values: list[str]) -> str:
            return values[0] if len(values) == 1 else f"COALESCE({', '.join(values)})"

        query = f"""
        WITH all_links AS (
            {union_sql}
        ){ir_cte}{vp_layer_cte}
        SELECT
            al.pmid,
            COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position) AS variant,
            v.gene_symbol,
            v.protein_notation,
            v.cdna_notation,
            {layer_select_expr} AS source_layer,
            {coalesce_sql(carriers_expr)} AS carriers_total,
            {coalesce_sql(affected_expr)} AS affected_count,
            {coalesce_sql(unaffected_expr)} AS unaffected_count,
            {coalesce_sql(uncertain_expr)} AS uncertain_count
        FROM all_links al
        JOIN variants v ON al.variant_id = v.variant_id
        {pd_join}
        {ir_join}
        {vp_join}
        """
        logger.info(
            "Executing union_all query across variant_papers + penetrance_data + individual_records"
        )
        df = pd.read_sql_query(query, conn)
        logger.info(
            f"union_all: {len(df)} rows from "
            f"{df['pmid'].nunique()} PMIDs, "
            f"{df['variant'].nunique()} variants"
        )

    elif strategy == "penetrance_data":
        # Join penetrance_data with variants to get variant notation
        query = """
            SELECT
                pd.pmid,
                COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position) as variant,
                v.gene_symbol,
                v.protein_notation,
                v.cdna_notation,
                'llm_text' as source_layer,
                pd.total_carriers_observed as carriers_total,
                pd.affected_count,
                pd.unaffected_count,
                pd.uncertain_count
            FROM penetrance_data pd
            JOIN variants v ON pd.variant_id = v.variant_id
        """
        logger.info("Executing penetrance_data query")
        df = pd.read_sql_query(query, conn)

    elif strategy == "individual_records":
        # Aggregate individual records by variant+pmid
        query = """
            SELECT
                ir.pmid,
                COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position) as variant,
                v.gene_symbol,
                v.protein_notation,
                v.cdna_notation,
                'llm_text' as source_layer,
                COUNT(*) as carriers_total,
                SUM(CASE WHEN ir.affected_status = 'affected' THEN 1 ELSE 0 END) as affected_count,
                SUM(CASE WHEN ir.affected_status = 'unaffected' THEN 1 ELSE 0 END) as unaffected_count,
                SUM(CASE WHEN ir.affected_status = 'uncertain' THEN 1 ELSE 0 END) as uncertain_count
            FROM individual_records ir
            JOIN variants v ON ir.variant_id = v.variant_id
            GROUP BY ir.pmid, v.variant_id
        """
        logger.info("Executing individual_records aggregation query")
        df = pd.read_sql_query(query, conn)

    elif strategy == "variant_papers_join":
        # Just get variant-paper associations (no counts available)
        has_vp_layer = "source_layer" in table_info["variant_papers"].columns
        has_vp_source_location = (
            "source_location" in table_info["variant_papers"].columns
        )
        has_vp_additional_notes = (
            "additional_notes" in table_info["variant_papers"].columns
        )
        layer_expr = source_layer_sql_case(
            "vp.source_location" if has_vp_source_location else "NULL",
            "vp.source_layer" if has_vp_layer else None,
            additional_notes_expr=(
                "vp.additional_notes" if has_vp_additional_notes else None
            ),
        )
        query = f"""
            SELECT
                vp.pmid,
                COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position) as variant,
                v.gene_symbol,
                v.protein_notation,
                v.cdna_notation,
                {layer_expr} as source_layer,
                NULL as carriers_total,
                NULL as affected_count,
                NULL as unaffected_count,
                NULL as uncertain_count
            FROM variant_papers vp
            JOIN variants v ON vp.variant_id = v.variant_id
        """
        logger.info("Executing variant_papers query (no counts available)")
        df = pd.read_sql_query(query, conn)
        logger.warning("No count data available in SQLite - using variant_papers only")

    else:
        # Fallback: try to query the combined table directly
        cols = primary_table.columns
        pmid_col = primary_table.pmid_column
        variant_cols = primary_table.variant_columns
        count_cols = primary_table.count_columns

        variant_select = (
            f"COALESCE({', '.join(variant_cols)})" if variant_cols else "NULL"
        )

        query = f"""
            SELECT
                {pmid_col} as pmid,
                {variant_select} as variant,
                NULL as gene_symbol,
                'llm_text' as source_layer,
                {count_cols.get("carriers_total", "NULL")} as carriers_total,
                {count_cols.get("affected_count", "NULL")} as affected_count,
                {count_cols.get("unaffected_count", "NULL")} as unaffected_count,
                {count_cols.get("uncertain_count", "NULL")} as uncertain_count
            FROM {primary_table.name}
        """
        logger.info(f"Executing custom query on {primary_table.name}")
        df = pd.read_sql_query(query, conn)

    # Normalize PMIDs (handles float conversion like 16470702.0 -> "16470702")
    df["pmid"] = df["pmid"].apply(normalize_pmid)
    df = filter_junk_sqlite_rows(df)

    logger.info(f"Extracted {len(df)} records from SQLite")
    return df


def filter_junk_sqlite_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Drop known-invalid figure/regex-table notation artifacts before scoring."""

    if df.empty:
        return df

    reasons: list[str | None] = []
    for _, row in df.iterrows():
        reasons.append(
            junk_notation_reason(
                source_layer=row.get("source_layer"),
                protein_notation=row.get("protein_notation"),
                cdna_notation=row.get("cdna_notation"),
                variant=row.get("variant"),
                gene_symbol=row.get("gene_symbol"),
            )
        )
    mask = pd.Series([reason is not None for reason in reasons], index=df.index)
    if not mask.any():
        return df

    counts = count_reasons(reason for reason in reasons if reason)
    breakdown = ", ".join(
        f"{reason}={count}" for reason, count in sorted(counts.items())
    )
    logger.info(
        "Dropped %d junk figure/regex_table SQLite rows before scoring (%s)",
        int(mask.sum()),
        breakdown,
    )
    return df.loc[~mask].copy()


# =============================================================================
# VARIANT NORMALIZATION
# =============================================================================

# Amino acid 3-letter to 1-letter mapping
AA_3_TO_1 = {
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
    "Ter": "*",
    "X": "*",
}

AA_1_TO_3 = {v: k for k, v in AA_3_TO_1.items() if k != "X"}
AA_1_TO_3["*"] = "Ter"


def normalize_unicode(s: str) -> str:
    """Normalize unicode characters (dashes, quotes, etc.)."""
    if not s:
        return s
    # Normalize to NFKC form
    s = unicodedata.normalize("NFKC", s)
    # Replace various dash characters with standard hyphen
    s = re.sub(r"[\u2010\u2011\u2012\u2013\u2014\u2015\u2212]", "-", s)
    # Replace fancy quotes
    s = re.sub(r"[\u2018\u2019]", "'", s)
    s = re.sub(r"[\u201c\u201d]", '"', s)
    return s


def normalize_variant(variant: str) -> str:
    """
    Normalize a variant string for comparison.

    NOTE: This is a simplified local normalizer for comparison purposes.
    For full variant normalization with gene-specific validation, use:
        from utils.variant_normalizer import normalize_variant as normalize_variant_full

    - Strips whitespace
    - Normalizes unicode
    - Collapses multiple spaces
    - Standardizes case for common prefixes
    - Normalizes termination symbols (X -> *)
    """
    if not variant or pd.isna(variant):
        return ""

    variant = str(variant).strip()
    variant = normalize_unicode(variant)
    variant = re.sub(r"\s+", " ", variant)

    # Standardize common prefixes
    variant = re.sub(r"^p\.\s*", "p.", variant, flags=re.IGNORECASE)
    variant = re.sub(r"^c\.\s*", "c.", variant, flags=re.IGNORECASE)

    # Normalize termination symbols: X, Ter -> *
    # Handle fsX, fsX123, X at end of variant
    variant = re.sub(r"fsX(\d*)$", r"fs*\1", variant, flags=re.IGNORECASE)
    variant = re.sub(r"Ter(\d*)$", r"*\1", variant)
    # Handle X alone at end (but not in middle of word)
    variant = re.sub(r"(\d)X$", r"\1*", variant)

    return variant


def to_canonical_form(variant: str) -> Optional[str]:
    """
    Reduce a variant to a single canonical form for matching.

    The canonical form uses:
    - Single-letter amino acid codes
    - No p. prefix
    - 'X' for stop codons
    - 'fsX' for frameshifts (no terminal position)
    - 'Del' for deletions
    - 'Ins' for insertions

    This matches the gold-standard heterozygote database format.

    Examples:
        p.Ala561Val     -> A561V
        p.A561V         -> A561V
        Ala561Val       -> A561V
        A561V           -> A561V
        p.Ala121Leufs*12 -> A121fsX
        p.H302fsX339    -> H302fsX
        A282fs          -> A282fsX
        p.Glu807*       -> E807X
        p.C39*          -> C39X
        p.Cys108Tyr*    -> C108Y  (trailing * on missense is artifact)
        p.Gly628del     -> G628Del
        Y475del         -> Y475Del
    """
    if not variant or pd.isna(variant):
        return None

    v = str(variant).strip()

    # Strip p. prefix
    if v.lower().startswith("p."):
        v = v[2:]

    # --- Numeric in-frame deletion with deleted residues: 73-73 DEL AAP -> P73Del ---
    m = re.match(
        r"^(\d+)(?:[-_](\d+))?\s*del\s*([A-Z]+)$",
        v,
        re.IGNORECASE,
    )
    if m:
        position = m.group(2) or m.group(1)
        deleted = m.group(3).upper()
        return f"{deleted[-1]}{position}Del"

    # --- Numeric in-frame insertion with inserted residues: 392INSW -> W392Ins ---
    m = re.match(r"^(\d+)\s*ins\s*([A-Z]+)$", v, re.IGNORECASE)
    if m:
        inserted = m.group(2).upper()
        return f"{inserted[0]}{m.group(1)}Ins"

    # --- Three-letter range deletion: Gln1507_Pro1509del -> Q1507_P1509Del ---
    m = re.match(
        r"^([A-Z][a-z]{2})(\d+)[_-]([A-Z][a-z]{2})(\d+)del([A-Z]*)\d*$",
        v,
        re.IGNORECASE,
    )
    if m:
        ref1 = AA_3_TO_1.get(m.group(1).title())
        ref2 = AA_3_TO_1.get(m.group(3).title())
        if ref1 and ref2:
            return f"{ref1}{m.group(2)}_{ref2}{m.group(4)}Del{m.group(5).upper()}"

    # --- Single-letter range deletion: K1505_Q1507DEL -> K1505_Q1507Del ---
    m = re.match(
        r"^([A-Z])(\d+)[_-]([A-Z])(\d+)del([A-Z]*)\d*$",
        v,
        re.IGNORECASE,
    )
    if m:
        return (
            f"{m.group(1).upper()}{m.group(2)}_"
            f"{m.group(3).upper()}{m.group(4)}Del{m.group(5).upper()}"
        )

    # --- Numeric protein insertion: 4944_4945INSH -> 4944_4945InsH ---
    m = re.match(r"^(\d+)_(\d+)ins([A-Z][A-Za-z]{0,2})$", v, re.IGNORECASE)
    if m:
        inserted = m.group(3)
        if len(inserted) == 3 and inserted.title() in AA_3_TO_1:
            inserted = AA_3_TO_1[inserted.title()]
        else:
            inserted = inserted.upper()
        return f"{m.group(1)}_{m.group(2)}Ins{inserted}"

    # --- Three-letter missense with trailing *: Ala429Pro* -> A429P ---
    m = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\*$", v)
    if m:
        ref = AA_3_TO_1.get(m.group(1))
        alt = AA_3_TO_1.get(m.group(3))
        if ref and alt:
            return f"{ref}{m.group(2)}{alt}"

    # --- Three-letter missense: Ala561Val -> A561V ---
    m = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", v)
    if m:
        ref = AA_3_TO_1.get(m.group(1))
        alt = AA_3_TO_1.get(m.group(3))
        if ref and alt:
            return f"{ref}{m.group(2)}{alt}"

    # --- Three-letter frameshift: Ala121Leufs*12, Pro926AlafsTer14 -> A121fsX ---
    m = re.match(
        r"^([A-Z][a-z]{2})(\d+)(?:[A-Z][a-z]{2})?"
        r"fs(?:Ter\d*|[/+]\*?\d+|[\*X]?\d*)$",
        v,
    )
    if m:
        ref = AA_3_TO_1.get(m.group(1))
        if ref:
            return f"{ref}{m.group(2)}fsX"

    # --- Single-letter frameshift: A282fs, V131fs/185, K897fs+*49 -> A282fsX ---
    m = re.match(
        r"^([A-Z])(\d+)[A-Z]?fs(?:Ter\d*|[/+]\*?\d+|[\*X]?\d*)$",
        v,
        re.IGNORECASE,
    )
    if m:
        return f"{m.group(1).upper()}{m.group(2)}fsX"

    # --- Malformed frameshift: Pro926AlaX14 -> P926fsX ---
    m = re.match(r"^([A-Z][a-z]{2})(\d+)[A-Z][a-z]{2}X\d+$", v)
    if m:
        ref = AA_3_TO_1.get(m.group(1))
        if ref:
            return f"{ref}{m.group(2)}fsX"

    # --- Three-letter stop: Glu807Ter, Arg1014* -> E807X ---
    m = re.match(r"^([A-Z][a-z]{2})(\d+)(?:Ter|\*|Stop)$", v)
    if m:
        ref = AA_3_TO_1.get(m.group(1))
        if ref:
            return f"{ref}{m.group(2)}X"

    # --- Single-letter stop: E807*, C39* -> E807X ---
    m = re.match(r"^([A-Z])(\d+)[\*]$", v)
    if m:
        return f"{m.group(1)}{m.group(2)}X"

    # --- Three-letter del: Gly628del -> G628Del ---
    m = re.match(r"^([A-Z][a-z]{2})(\d+)(del)$", v, re.IGNORECASE)
    if m:
        ref = AA_3_TO_1.get(m.group(1))
        if ref:
            return f"{ref}{m.group(2)}Del"

    # --- Single-letter del: A671del -> A671Del ---
    m = re.match(r"^([A-Z])(\d+)(del)$", v, re.IGNORECASE)
    if m:
        return f"{m.group(1)}{m.group(2)}Del"

    # --- Three-letter dup: Ile82dup -> I82dup ---
    m = re.match(r"^([A-Z][a-z]{2})(\d+)(dup)$", v, re.IGNORECASE)
    if m:
        ref = AA_3_TO_1.get(m.group(1))
        if ref:
            return f"{ref}{m.group(2)}dup"

    # --- Single-letter dup: I82dup -> I82dup ---
    m = re.match(r"^([A-Z])(\d+)(dup)$", v, re.IGNORECASE)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}dup"

    # --- Three-letter ins: Gly189ins -> G189Ins ---
    m = re.match(r"^([A-Z][a-z]{2})(\d+)(ins)$", v, re.IGNORECASE)
    if m:
        ref = AA_3_TO_1.get(m.group(1))
        if ref:
            return f"{ref}{m.group(2)}Ins"

    # --- Single-letter ins: G189ins -> G189Ins ---
    m = re.match(r"^([A-Z])(\d+)(ins)$", v, re.IGNORECASE)
    if m:
        return f"{m.group(1)}{m.group(2)}Ins"

    # --- Single-letter stop with X (already canonical) ---
    m = re.match(r"^([A-Z])(\d+)X$", v)
    if m:
        return v

    # --- Single-letter missense: A561V (already canonical) ---
    m = re.match(r"^[A-Z]\d+[A-Z]$", v)
    if m:
        return v

    # --- Splice: A715sp (already canonical) ---
    m = re.match(r"^([A-Z])(\d+)sp$", v, re.IGNORECASE)
    if m:
        return f"{m.group(1)}{m.group(2)}sp"

    # --- Three-letter splice: Ala715sp -> A715sp ---
    m = re.match(r"^([A-Z][a-z]{2})(\d+)sp$", v, re.IGNORECASE)
    if m:
        ref = AA_3_TO_1.get(m.group(1))
        if ref:
            return f"{ref}{m.group(2)}sp"

    # --- IVS notation: pass through ---
    if v.startswith("IVS"):
        return v

    # Could not canonicalize
    return None


def convert_aa_3_to_1(variant: str) -> Optional[str]:
    """
    Convert 3-letter amino acid codes to 1-letter in a variant string.

    e.g., p.Arg123His -> p.R123H
    Also handles: p.Arg123del, p.Arg123fs, p.Arg123*, p.Arg123_Lys130del
    """
    if not variant or not variant.lower().startswith("p."):
        return None

    result = variant[:2]  # Keep "p."
    remaining = variant[2:]

    # Pattern 1: Handle range deletions like p.Arg123_Lys130del before the
    # generic single-position pattern can consume only p.Arg123.
    range_pattern = r"([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)(.*)"
    range_match = re.match(range_pattern, remaining)
    if range_match:
        aa1, pos1, aa2, pos2, rest = range_match.groups()
        aa1_short = AA_3_TO_1.get(aa1, aa1)
        aa2_short = AA_3_TO_1.get(aa2, aa2)
        return f"p.{aa1_short}{pos1}_{aa2_short}{pos2}{rest}"

    # Pattern 2: 3-letter AA, position, 3-letter AA (with optional fs, del, etc.)
    pattern = r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})?(.*)$"
    match = re.match(pattern, remaining)

    if match:
        aa1, pos, aa2, rest = match.groups()
        aa1_short = AA_3_TO_1.get(aa1, aa1)
        aa2_short = AA_3_TO_1.get(aa2, aa2) if aa2 else ""
        result = f"p.{aa1_short}{pos}{aa2_short}{rest or ''}"
        return result

    return None


def convert_aa_1_to_3(variant: str) -> Optional[str]:
    """
    Convert 1-letter amino acid codes to 3-letter in a variant string.

    e.g., p.R123H -> p.Arg123His
    Also handles: R123H (without p. prefix), R123del, R123fs, R123*
    """
    # Handle variants without p. prefix
    if (
        variant
        and not variant.lower().startswith("p.")
        and not variant.lower().startswith("c.")
    ):
        # Check if it looks like a protein variant (letter + number + letter/suffix)
        if re.match(r"^[A-Z]\d+[A-Z*]", variant, re.IGNORECASE):
            variant = f"p.{variant}"

    if not variant or not variant.lower().startswith("p."):
        return None

    # Pattern 1: Handle range deletions like p.R123_K130del before the generic
    # single-position pattern can consume only p.R123.
    range_pattern = r"^p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)(.*)$"
    range_match = re.match(range_pattern, variant)
    if range_match:
        aa1, pos1, aa2, pos2, rest = range_match.groups()
        aa1_long = AA_1_TO_3.get(aa1, aa1)
        aa2_long = AA_1_TO_3.get(aa2, aa2)
        return f"p.{aa1_long}{pos1}_{aa2_long}{pos2}{rest}"

    # Pattern 2: single letter AA, position, single letter AA (with optional suffix)
    pattern = r"^p\.([A-Z\*])(\d+)([A-Z\*])?(.*)$"
    match = re.match(pattern, variant)

    if match:
        aa1, pos, aa2, rest = match.groups()
        aa1_long = AA_1_TO_3.get(aa1, aa1)
        aa2_long = AA_1_TO_3.get(aa2, aa2) if aa2 else ""
        result = f"p.{aa1_long}{pos}{aa2_long}{rest or ''}"
        return result

    return None


def add_p_prefix(variant: str) -> Optional[str]:
    """
    Add p. prefix to a variant that looks like a protein change but lacks the prefix.

    e.g., A123V -> p.A123V, Ala123Val -> p.Ala123Val
    """
    if not variant:
        return None
    if variant.lower().startswith("p.") or variant.lower().startswith("c."):
        return None  # Already has prefix

    # Check if it looks like a protein variant:
    # - Single letter + number + single letter (A123V)
    # - Three letter + number + three letter (Ala123Val)
    # - Single/three letter + number + fs* (A123fs*)
    # - Single/three letter + number + * (A123*)
    aa = r"(?:[A-Z][a-z]{2}|[A-Z])"
    protein_pattern = (
        rf"^(?:{aa}\d+(?:_{aa}\d+)?"
        rf"(?:{aa}|fs\*?\d*|\*\d*|del\d*|dup|ins{aa}*)"
        rf"|\d+_\d+ins{aa}+)$"
    )
    if re.match(protein_pattern, variant, re.IGNORECASE):
        return f"p.{variant}"

    return None


def remove_p_prefix(variant: str) -> Optional[str]:
    """
    Remove p. prefix from a protein variant.

    e.g., p.A123V -> A123V
    """
    if variant and variant.lower().startswith("p."):
        return variant[2:]
    return None


def get_variant_forms(variant: str) -> Set[str]:
    """
    Get all equivalent forms of a variant for matching.

    Returns set of normalized variant strings including:
    - Original normalized
    - Canonical form (single-letter, no prefix, gold-standard format)
    - 3-letter to 1-letter conversion
    - 1-letter to 3-letter conversion
    - With and without p. prefix
    """
    forms = set()

    normalized = normalize_variant(variant)
    if normalized:
        forms.add(normalized)

    # Add canonical form (the most important for matching)
    canonical = to_canonical_form(variant)
    if canonical:
        forms.add(canonical)
        # Also add p.-prefixed canonical
        forms.add(f"p.{canonical}")

    # Try amino acid conversions on normalized form
    converted_1 = convert_aa_3_to_1(normalized)
    if converted_1:
        forms.add(normalize_variant(converted_1))

    converted_3 = convert_aa_1_to_3(normalized)
    if converted_3:
        forms.add(normalize_variant(converted_3))

    # Try adding p. prefix if missing
    with_prefix = add_p_prefix(normalized)
    if with_prefix:
        forms.add(normalize_variant(with_prefix))
        # Also convert this form
        converted_1_p = convert_aa_3_to_1(with_prefix)
        if converted_1_p:
            forms.add(normalize_variant(converted_1_p))
        converted_3_p = convert_aa_1_to_3(with_prefix)
        if converted_3_p:
            forms.add(normalize_variant(converted_3_p))

    # Try removing p. prefix
    without_prefix = remove_p_prefix(normalized)
    if without_prefix:
        forms.add(normalize_variant(without_prefix))
        # Also convert this form
        converted_1_np = convert_aa_3_to_1(f"p.{without_prefix}")
        if converted_1_np:
            # Add both with and without prefix
            forms.add(normalize_variant(converted_1_np))
            forms.add(
                normalize_variant(converted_1_np[2:])
                if converted_1_np.startswith("p.")
                else converted_1_np
            )
        converted_3_np = convert_aa_1_to_3(f"p.{without_prefix}")
        if converted_3_np:
            forms.add(normalize_variant(converted_3_np))
            forms.add(
                normalize_variant(converted_3_np[2:])
                if converted_3_np.startswith("p.")
                else converted_3_np
            )

    return forms


# =============================================================================
# FUZZY MATCHING
# =============================================================================


def compute_similarity(s1: str, s2: str) -> float:
    """
    Compute similarity ratio between two strings.

    Uses rapidfuzz if available, otherwise difflib.
    Returns value between 0 and 1.
    """
    if not s1 or not s2:
        return 0.0

    if RAPIDFUZZ_AVAILABLE:
        return rapidfuzz_fuzz.ratio(s1, s2) / 100.0
    else:
        return difflib.SequenceMatcher(None, s1, s2).ratio()


def _positional_digits(variant: str) -> Tuple[int, ...]:
    """
    Extract every run of digits from a variant string, as a tuple of ints.

    Two variants whose digit tuples differ refer to different positions
    (G572S = (572,) vs G628S = (628,)) and must never be fuzzy-matched.
    Returns an empty tuple for variants with no extractable digits.
    """
    return tuple(int(d) for d in re.findall(r"\d+", str(variant or "")))


def _positions_compatible(a: str, b: str) -> bool:
    """
    True if the position-digit tuples of ``a`` and ``b`` overlap as subsets.

    Used as a guard on fuzzy matches:

    - Same-position variants: G572S = (572,), G572R = (572,) → compatible.
    - Different-position variants: G572S = (572,), G628S = (628,) →
      incompatible. (Levenshtein happily calls these 81.8% similar
      because position digits carry no semantic weight to the distance
      metric.)
    - Frameshift with offset: p.Pro926AlafsTer14 = (926, 14) vs P926fsX
      = (926,). The smaller set is a subset of the larger, so the guard
      allows the match. This is the cross-notation case the bridge needs.
    - Either side with no digits (IVS splice notation): return True and
      defer to similarity scoring.
    """
    da, db = _positional_digits(a), _positional_digits(b)
    if not da or not db:
        return True
    # cDNA coordinates carry intronic offsets and ranges where subset matching
    # is unsafe: c.1890+5G>A and c.1890G>A share 1890 but are different
    # variants. Cross-notation indel recovery is handled by the explicit
    # cDNA/protein bridge instead of fuzzy matching.
    a_is_cdna = str(a or "").strip().lower().startswith("c.")
    b_is_cdna = str(b or "").strip().lower().startswith("c.")
    if a_is_cdna or b_is_cdna:
        return da == db
    sa, sb = set(da), set(db)
    if len(sa) <= len(sb):
        return sa.issubset(sb)
    return sb.issubset(sa)


def find_best_match(
    query_variant: str,
    candidates: List[str],
    threshold: float = 0.85,
    consumed: Optional[Set[str]] = None,
) -> Tuple[Optional[str], float, str]:
    """
    Find the best matching variant from candidates.

    Args:
        query_variant: Variant to match.
        candidates: List of candidate variants (raw strings).
        threshold: Minimum similarity threshold for fuzzy matches.
        consumed: Optional set of candidate raw strings already claimed by
            an earlier query in this PMID — those are skipped to enforce
            1-to-1 assignment (prevents one sqlite variant matching
            multiple gold rows).

    Returns:
        Tuple of (best_match, score, match_type).
    """
    consumed = consumed or set()
    query_forms = get_variant_forms(query_variant)

    # First try exact match on any form.
    for candidate in candidates:
        if candidate in consumed:
            continue
        candidate_forms = get_variant_forms(candidate)
        if query_forms & candidate_forms:
            return (candidate, 1.0, "exact")

    # Fuzzy matching with a positional-digit guard. Only consider candidates
    # whose position digits overlap with the query — this prevents
    # Levenshtein from declaring G572S ≈ G628S as 81.8% similar.
    best_match = None
    best_score = 0.0

    for candidate in candidates:
        if candidate in consumed:
            continue
        if not _positions_compatible(query_variant, candidate):
            continue

        cand_forms = get_variant_forms(candidate)
        for q_form in query_forms:
            for c_form in cand_forms:
                if not _positions_compatible(q_form, c_form):
                    continue
                score = compute_similarity(q_form, c_form)
                if score > best_score:
                    best_score = score
                    best_match = candidate

    if best_score >= threshold:
        return (best_match, best_score, "fuzzy")

    return (None, best_score, "none")


def _cdna_indel_protein_positions(cdna: str) -> Optional[Tuple[int, int]]:
    """
    From a cDNA indel notation, compute the (start, end) protein-aa positions
    in 1-based amino-acid coordinates via standard codon math (aa = (nt-1)//3 + 1).

    Examples:
        c.842dupG          → (281, 281)
        c.1631_1632delAG   → (544, 544)
        c.221_242del       → (74, 81)
        c.547_553delGGCGG  → (183, 185)

    Returns None for non-indel cDNA notations (e.g. ``c.1558-1G>C`` splice
    sites) where there's no clean protein-position equivalent.
    """
    if not cdna or not isinstance(cdna, str):
        return None
    m = re.match(
        r"^c\.(\d+)(?:[+\-]\d+)?(?:_(\d+)(?:[+\-]\d+)?)?\s*(del|dup|ins)",
        cdna.strip(),
        re.IGNORECASE,
    )
    if not m:
        return None
    start = int(m.group(1))
    end = int(m.group(2)) if m.group(2) else start
    aa_start = (start - 1) // 3 + 1
    aa_end = (end - 1) // 3 + 1
    return (aa_start, aa_end)


def _protein_indel_position(variant: str) -> Optional[Tuple[str, int, str]]:
    """
    From a protein-level indel/frameshift/stop canonical form, return
    ``(ref_aa, position, op)`` where ``op`` is one of: ``"fs"``, ``"del"``,
    ``"ins"``, ``"X"``.

    Examples:
        R281fsX  → ('R', 281, 'fs')
        T74fsX   → ('T', 74, 'fs')
        V215X    → ('V', 215, 'X')
        G628Del  → ('G', 628, 'del')

    Returns None when the variant isn't a recognisable protein-level
    indel/stop form.
    """
    if not variant:
        return None
    m = re.match(
        r"^([A-Z])(\d+)(?:_[A-Z]\d+)?(fsX\d*|fs|Del|Ins|X)$",
        str(variant).strip(),
    )
    if not m:
        return None
    aa, pos, op_raw = m.group(1), int(m.group(2)), m.group(3)
    op = op_raw.lower()
    if op.startswith("fs"):
        op = "fs"
    elif op == "del":
        op = "del"
    elif op == "ins":
        op = "ins"
    elif op == "x":
        op = "X"
    return (aa, pos, op)


# =============================================================================
# COMPARISON LOGIC
# =============================================================================


@dataclass
class ComparisonRow:
    """Single row of comparison results."""

    pmid: str
    excel_variant_raw: str
    excel_variant_norm: str
    sqlite_variant_raw: Optional[str]
    sqlite_variant_norm: Optional[str]
    match_type: str  # 'exact', 'fuzzy', 'none'
    match_score: Optional[float]

    excel_carriers_total: Optional[int]
    sqlite_carriers_total: Optional[int]
    carriers_diff: Optional[int]

    excel_affected: Optional[int]
    sqlite_affected: Optional[int]
    affected_diff: Optional[int]

    excel_unaffected: Optional[int]
    sqlite_unaffected: Optional[int]
    unaffected_diff: Optional[int]

    phenotype: str = "ALL"

    missing_in_sqlite: bool = False
    missing_in_excel: bool = False
    count_mismatch: bool = False
    sqlite_source_layer: Optional[str] = None


def safe_int(value: Any) -> Optional[int]:
    """Safely convert a value to int, returning None for invalid values."""
    if value is None or pd.isna(value):
        return None
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return None


def compute_diff(a: Optional[int], b: Optional[int]) -> Optional[int]:
    """Compute difference between two optional integers."""
    if a is None or b is None:
        return None
    return a - b


def load_excel_data(
    excel_path: Path,
    sheet_name: Optional[str],
    column_mapping: Optional[Dict[str, str]],
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Load and parse curated gold-standard data.

    Returns:
        Tuple of (DataFrame, detected_columns)
    """
    logger.info(f"Loading gold-standard file: {excel_path}")

    suffix = excel_path.suffix.lower()

    try:
        if suffix == ".csv":
            df = pd.read_csv(excel_path)
        else:
            # Select engine based on file suffix so legacy .xls files are supported
            engine = "xlrd" if suffix == ".xls" else "openpyxl"
            # If sheet_name is None, default to first sheet (index 0)
            read_sheet = sheet_name if sheet_name is not None else 0
            df = pd.read_excel(excel_path, sheet_name=read_sheet, engine=engine)
            # Handle case where read_excel returns a dict (when sheet_name is a list or None with multiple sheets)
            if isinstance(df, dict):
                # Get the first sheet
                first_sheet = list(df.keys())[0]
                logger.info(f"Multiple sheets found, using first sheet: {first_sheet}")
                df = df[first_sheet]
    except ImportError as exc:
        if suffix == ".xls":
            raise ImportError(
                "Reading .xls files requires the optional dependency xlrd. "
                "Install it (e.g., `pip install xlrd`) or convert the file to "
                ".xlsx/.csv before retrying."
            ) from exc
        raise
    logger.info(f"Loaded {len(df)} rows from gold-standard file")
    logger.info(f"Columns: {list(df.columns)}")

    detected = detect_columns(df, column_mapping)

    # Validate required columns
    if detected["pmid"] is None:
        raise ValueError(
            f"Cannot find PMID column. Available columns: {list(df.columns)}\n"
            "Use --mapping to specify column mappings."
        )

    if detected["variant"] is None:
        raise ValueError(
            f"Cannot find variant column. Available columns: {list(df.columns)}\n"
            "Use --mapping to specify column mappings."
        )

    return df, detected


def aggregate_excel_data(
    df: pd.DataFrame, detected: Dict[str, Any]
) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """
    Aggregate Excel data by (pmid, variant) key.

    Returns:
        Dict mapping (pmid, variant_norm) -> aggregated data
    """
    pmid_col = detected["pmid"]
    variant_col = detected["variant"]
    carriers_col = detected.get("carriers_total")
    affected_col = detected.get("affected_count")
    unaffected_col = detected.get("unaffected_count")
    phenotype_col = detected.get("phenotype")
    affected_phenotype_cols = detected.get("affected_phenotype_columns", [])

    aggregated = {}

    for _, row in df.iterrows():
        pmid = normalize_pmid(row[pmid_col])
        variant_raw = str(row[variant_col]) if pd.notna(row[variant_col]) else ""
        variant_norm = normalize_variant(variant_raw)
        # Use canonical form as primary key for better cross-notation matching
        # e.g., "p.Ala561Val" and "A561V" both canonicalize to "A561V"
        canonical = to_canonical_form(variant_raw)
        variant_key = canonical if canonical else variant_norm

        if not pmid or not variant_key:
            continue

        key = (pmid, variant_key)

        if key not in aggregated:
            aggregated[key] = {
                "pmid": pmid,
                "variant_raw": variant_raw,
                "variant_norm": variant_key,
                "carriers_total": 0,
                "affected_count": 0,
                "unaffected_count": 0,
                "phenotypes": set(),
            }

        # Aggregate counts
        if carriers_col and pd.notna(row.get(carriers_col)):
            aggregated[key]["carriers_total"] += safe_int(row[carriers_col]) or 0

        # Sum affected from multiple phenotype columns if available
        if affected_phenotype_cols:
            row_affected = 0
            for pheno_col in affected_phenotype_cols:
                if pd.notna(row.get(pheno_col)):
                    row_affected += safe_int(row[pheno_col]) or 0
            aggregated[key]["affected_count"] += row_affected
        elif affected_col and pd.notna(row.get(affected_col)):
            # Fallback to single affected column
            aggregated[key]["affected_count"] += safe_int(row[affected_col]) or 0

        if unaffected_col and pd.notna(row.get(unaffected_col)):
            aggregated[key]["unaffected_count"] += safe_int(row[unaffected_col]) or 0

        if phenotype_col and pd.notna(row.get(phenotype_col)):
            aggregated[key]["phenotypes"].add(str(row[phenotype_col]))

    logger.info(
        f"Aggregated Excel data into {len(aggregated)} unique (pmid, variant) pairs"
    )
    return aggregated


def aggregate_sqlite_data(df: pd.DataFrame) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """
    Aggregate SQLite data by (pmid, variant) key.

    Returns:
        Dict mapping (pmid, variant_norm) -> aggregated data
    """
    aggregated = {}

    for _, row in df.iterrows():
        pmid = normalize_pmid(row["pmid"])
        variant_raw = str(row["variant"]) if pd.notna(row["variant"]) else ""
        variant_norm = normalize_variant(variant_raw)
        # Use canonical form as primary key (matches Excel aggregation)
        canonical = to_canonical_form(variant_raw)
        variant_key = canonical if canonical else variant_norm

        if not pmid or not variant_key:
            continue

        key = (pmid, variant_key)

        raw_source_layer = str(row.get("source_layer") or "").strip()
        source_layer = normalize_source_layer(raw_source_layer) or "llm_text"
        if "," in raw_source_layer:
            source_layer = "mixed"

        if key not in aggregated:
            aggregated[key] = {
                "pmid": pmid,
                "variant_raw": variant_raw,
                "variant_norm": variant_key,
                "carriers_total": 0,
                "affected_count": 0,
                "unaffected_count": 0,
                "protein_notation": row.get("protein_notation"),
                "cdna_notation": row.get("cdna_notation"),
                "source_layer": source_layer,
            }
        elif source_layer != aggregated[key].get("source_layer"):
            aggregated[key]["source_layer"] = "mixed"

        # Aggregate counts (handle None values)
        aggregated[key]["carriers_total"] += safe_int(row.get("carriers_total")) or 0
        aggregated[key]["affected_count"] += safe_int(row.get("affected_count")) or 0
        aggregated[key]["unaffected_count"] += (
            safe_int(row.get("unaffected_count")) or 0
        )

    logger.info(
        f"Aggregated SQLite data into {len(aggregated)} unique (pmid, variant) pairs"
    )
    return aggregated


_INDEL_CANONICAL_RE = re.compile(
    r"^(?:[A-Z]\d+(?:_[A-Z]\d+)?(?:fsX|Del|Ins|dup)|\d+_\d+Ins[A-Z]+)$",
    re.IGNORECASE,
)
# Match indel tokens in cDNA notation. "dup", "del", "ins", and "fs" may be
# followed immediately by a base letter (e.g. c.842dupG, c.362delA), so word
# boundaries on the right side are wrong here — keep only the left boundary.
_CDNA_INDEL_TOKEN_RE = re.compile(r"(?<![A-Za-z])(?:dup|del|ins|fs)", re.IGNORECASE)


def _is_indel_canonical(canonical: Optional[str]) -> bool:
    """Return True if a canonical form represents a protein indel/frameshift."""
    if not canonical:
        return False
    return bool(_INDEL_CANONICAL_RE.match(canonical))


def _cdna_implied_codon(cdna: Optional[str]) -> Optional[int]:
    """
    Map a cDNA position to its implied protein codon (1-indexed, ceil division).

    Returns the codon of the first cDNA position cited. Returns None if the
    string has no parseable cDNA position or is not an indel-type cDNA.
    """
    if not cdna:
        return None
    if not _CDNA_INDEL_TOKEN_RE.search(cdna):
        return None
    m = re.search(r"c\.?\s*(\d+)", cdna)
    if not m:
        m = re.search(r"(\d+)", cdna)
    if not m:
        return None
    try:
        cdna_pos = int(m.group(1))
    except ValueError:
        return None
    if cdna_pos <= 0:
        return None
    return (cdna_pos + 2) // 3


def _protein_codon_from_canonical(canonical: Optional[str]) -> Optional[int]:
    """Extract the protein codon position from a canonical form like R281fsX."""
    if not canonical:
        return None
    m = re.match(r"^[A-Z](\d+)", canonical)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


# Coding-region cDNA substitution (requires literal c. prefix + base>base, so
# intronic forms like c.703+1G>A are not mistaken for a coding position).
_CDNA_SUBSTITUTION_RE = re.compile(r"c\.\s*(\d+)\s*[ACGTacgt]\s*>\s*[ACGTacgt]")
# Canonical protein substitution / nonsense like Y51X, Y111C, P117L (single
# residue, not an indel/frameshift, which are handled by the indel bridge).
_PROTEIN_SUBSTITUTION_CANON_RE = re.compile(r"^[A-Z](\d+)[A-Z*X]$")


def _cdna_substitution_codon(cdna: Optional[str]) -> Optional[int]:
    """Implied protein codon for a coding cDNA substitution (c.153C>A -> 51)."""
    if not cdna:
        return None
    m = _CDNA_SUBSTITUTION_RE.search(cdna)
    if not m:
        return None
    try:
        pos = int(m.group(1))
    except ValueError:
        return None
    return (pos + 2) // 3 if pos > 0 else None


def _protein_substitution_codon(canonical: Optional[str]) -> Optional[int]:
    """Codon of a protein substitution/nonsense canonical form (Y51X -> 51).

    Returns None for indel/frameshift forms (handled by the indel bridge).
    """
    if not canonical or _is_indel_canonical(canonical):
        return None
    m = _PROTEIN_SUBSTITUTION_CANON_RE.match(canonical)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def _find_cdna_protein_bridge(
    excel_entry: Dict[str, Any],
    available: List[Tuple[Tuple[str, str], Dict[str, Any]]],
) -> Optional[Tuple[Tuple[str, str], Dict[str, Any]]]:
    """
    Find a SQLite candidate whose cDNA indel notation implies the same protein
    codon as the gold protein indel.

    Gold protein notation like R281fsX maps to codon 281; an LLM-extracted
    cDNA dup/del/ins at base 841-843 maps to the same codon. The matcher
    otherwise sees only the raw cDNA string and misses the equivalence.

    Returns the matching (sqlite_key, sqlite_entry) or None.
    """
    excel_raw = excel_entry.get("variant_raw", "")
    excel_canonical = to_canonical_form(excel_raw)

    # Protein gold -> cDNA SQLite bridge.
    target_codon = None
    if _is_indel_canonical(excel_canonical):
        target_codon = _protein_codon_from_canonical(excel_canonical)

        if target_codon is not None:
            for sqlite_key, entry in available:
                cdna = entry.get("cdna_notation")
                if not cdna or not isinstance(cdna, str):
                    continue
                protein_range = _cdna_indel_protein_positions(cdna)
                if protein_range:
                    start, end = protein_range
                    if start <= target_codon <= end:
                        return (sqlite_key, entry)

                implied = _cdna_implied_codon(cdna)
                if implied is None:
                    continue
                # Exact codon match. A frameshift caused by a single-base indel at
                # codon N starts at codon N; tolerate ±0 by default to keep precision
                # high — the diagnosis only reports recoverable matches at the right
                # codon.
                if implied == target_codon:
                    return (sqlite_key, entry)

    # Protein-substitution / nonsense gold -> cDNA-substitution SQLite bridge.
    # Gold like Y51X or P117L maps to a codon; a stored cDNA SNV (e.g. c.153C>A)
    # implies the same codon and is the SAME variant in cDNA notation. Require a
    # UNIQUE cDNA candidate at that codon on this PMID to avoid mispairing the
    # two distinct substitutions that can share a codon (e.g. R190W vs R190Q).
    if target_codon is None:
        sub_codon = _protein_substitution_codon(excel_canonical)
        if sub_codon is not None:
            sub_cands = [
                (sqlite_key, entry)
                for sqlite_key, entry in available
                if isinstance(entry.get("cdna_notation"), str)
                and _cdna_substitution_codon(entry.get("cdna_notation")) == sub_codon
            ]
            if len(sub_cands) == 1:
                return sub_cands[0]

    # cDNA gold -> protein SQLite bridge. This is the reverse case: SQLite may
    # prefer protein_notation as its display key even though the stored cDNA
    # notation, or an equivalent protein frameshift, is what the gold row uses.
    excel_range = _cdna_indel_protein_positions(str(excel_raw))
    if not excel_range:
        return None
    start, end = excel_range
    for sqlite_key, entry in available:
        for protein in (entry.get("protein_notation"), entry.get("variant_raw")):
            protein_canonical = to_canonical_form(str(protein or ""))
            protein_info = _protein_indel_position(protein_canonical or "")
            if not protein_info:
                continue
            _, protein_pos, _ = protein_info
            # HGVS protein frameshift positions can be one codon before the
            # first changed cDNA codon when the nucleotide event falls at a
            # codon boundary. Keep this narrow to avoid cross-position matches.
            if start - 1 <= protein_pos <= end + 1:
                return (sqlite_key, entry)
    return None


def _entry_variant_forms(entry: Dict[str, Any]) -> Set[str]:
    """Return every comparable variant form stored on a SQLite aggregate row."""

    forms: Set[str] = set()
    for variant_field in ("variant_raw", "protein_notation", "cdna_notation"):
        value = entry.get(variant_field)
        if value is None or pd.isna(value):
            continue
        text = str(value).strip()
        if not text:
            continue
        forms.update(get_variant_forms(text))
    return forms


def compare_data(
    excel_data: Dict[Tuple[str, str], Dict[str, Any]],
    sqlite_data: Dict[Tuple[str, str], Dict[str, Any]],
    match_mode: str,
    fuzzy_threshold: float,
) -> List[ComparisonRow]:
    """
    Compare Excel and SQLite data with greedy 1-to-1 assignment.

    Two-pass matching ensures each SQLite variant is claimed by at most one
    Excel row: pass 1 takes direct canonical-key matches (highest confidence);
    pass 2 handles form-intersection, fuzzy, and cDNA->protein bridge matches
    against the *remaining* SQLite candidates.

    Args:
        excel_data: Aggregated Excel data, keyed by (pmid, canonical_variant)
        sqlite_data: Aggregated SQLite data, keyed by (pmid, canonical_variant)
        match_mode: 'exact' or 'fuzzy'
        fuzzy_threshold: Threshold for fuzzy matching

    Returns:
        List of ComparisonRow objects
    """
    results: List[ComparisonRow] = []
    consumed_sqlite_keys: Set[Tuple[str, str]] = set()
    matched_excel_keys: Set[Tuple[str, str]] = set()

    # Index SQLite entries by PMID for per-PMID candidate lookup
    sqlite_by_pmid: Dict[str, List[Tuple[Tuple[str, str], Dict[str, Any]]]] = {}
    for sqlite_key, entry in sqlite_data.items():
        sqlite_by_pmid.setdefault(sqlite_key[0], []).append((sqlite_key, entry))

    def _available_for(pmid: str) -> List[Tuple[Tuple[str, str], Dict[str, Any]]]:
        return [
            (sqlite_key, entry)
            for (sqlite_key, entry) in sqlite_by_pmid.get(pmid, [])
            if sqlite_key not in consumed_sqlite_keys
        ]

    # Pass 1: direct canonical-key exact matches.
    # Both sides aggregate by canonical form (see aggregate_*_data), so an
    # identical key on both sides is the strongest possible match.
    for excel_key, excel_entry in excel_data.items():
        if excel_key in sqlite_data and excel_key not in consumed_sqlite_keys:
            sqlite_entry = sqlite_data[excel_key]
            consumed_sqlite_keys.add(excel_key)
            matched_excel_keys.add(excel_key)
            results.append(
                create_comparison_row(excel_entry, sqlite_entry, "exact", 1.0)
            )

    # Pass 2: remaining excel items — form-intersection, fuzzy, or cDNA bridge
    for excel_key, excel_entry in excel_data.items():
        if excel_key in matched_excel_keys:
            continue
        excel_pmid, _ = excel_key
        available = _available_for(excel_pmid)

        row: Optional[ComparisonRow] = None
        excel_forms = get_variant_forms(excel_entry["variant_raw"])

        # Exact form-intersection against all notations carried by the SQLite
        # row. The aggregate display key prefers protein_notation, but a row can
        # also carry an exact cDNA match that would otherwise be hidden.
        for sqlite_key, entry in available:
            if excel_forms & _entry_variant_forms(entry):
                consumed_sqlite_keys.add(sqlite_key)
                matched_excel_keys.add(excel_key)
                row = create_comparison_row(excel_entry, entry, "exact", 1.0)
                break

        if row is None:
            bridge = _find_cdna_protein_bridge(excel_entry, available)
            if bridge:
                sqlite_key, entry = bridge
                consumed_sqlite_keys.add(sqlite_key)
                matched_excel_keys.add(excel_key)
                row = create_comparison_row(
                    excel_entry,
                    entry,
                    (
                        "fuzzy_cdna_bridge"
                        if match_mode == "fuzzy"
                        else "exact_cdna_bridge"
                    ),
                    1.0,
                )

        if row is None and match_mode == "fuzzy":
            # Fuzzy mode: compare every form stored on each available SQLite
            # row, not just the preferred display variant. This preserves the
            # positional guard while allowing cDNA/protein alternate notations
            # from the same row to participate in matching.
            best_key: Optional[Tuple[str, str]] = None
            best_entry: Optional[Dict[str, Any]] = None
            best_score = 0.0
            for sqlite_key, entry in available:
                for q_form in excel_forms:
                    for c_form in _entry_variant_forms(entry):
                        if not _positions_compatible(q_form, c_form):
                            continue
                        score = compute_similarity(q_form, c_form)
                        if score > best_score:
                            best_score = score
                            best_key = sqlite_key
                            best_entry = entry
            if (
                best_entry is not None
                and best_key is not None
                and best_score >= fuzzy_threshold
            ):
                consumed_sqlite_keys.add(best_key)
                matched_excel_keys.add(excel_key)
                row = create_comparison_row(
                    excel_entry, best_entry, "fuzzy", best_score
                )

            # Retain the legacy raw-string matcher as a fallback in case a
            # caller relies on its exact return semantics for unusual inputs.
            if row is None:
                candidate_raws = [entry["variant_raw"] for (_, entry) in available]
                best_match, score, match_type = find_best_match(
                    excel_entry["variant_raw"], candidate_raws, fuzzy_threshold
                )
                if best_match:
                    for sqlite_key, entry in available:
                        if entry["variant_raw"] == best_match:
                            consumed_sqlite_keys.add(sqlite_key)
                            matched_excel_keys.add(excel_key)
                            row = create_comparison_row(
                                excel_entry, entry, match_type, score
                            )
                            break

        if row is None and match_mode == "exact":
            # No further exact-only work needed; form and cDNA bridge attempts
            # already ran above.
            pass
        elif row is None and match_mode != "fuzzy":
            # Preserve future non-fuzzy modes by falling back to the legacy
            # exact raw-string helper.
            candidate_raws = [entry["variant_raw"] for (_, entry) in available]
            best_match, score, match_type = find_best_match(
                excel_entry["variant_raw"], candidate_raws, fuzzy_threshold
            )
            if best_match:
                for sqlite_key, entry in available:
                    if entry["variant_raw"] == best_match:
                        consumed_sqlite_keys.add(sqlite_key)
                        matched_excel_keys.add(excel_key)
                        row = create_comparison_row(
                            excel_entry, entry, match_type, score
                        )
                        break

        if row is None:
            row = create_comparison_row(excel_entry, None, "none", None)
            row.missing_in_sqlite = True

        results.append(row)

    # Add SQLite entries not consumed by any Excel row
    for key, sqlite_entry in sqlite_data.items():
        if key in consumed_sqlite_keys:
            continue
        results.append(
            ComparisonRow(
                pmid=sqlite_entry["pmid"],
                excel_variant_raw="",
                excel_variant_norm="",
                sqlite_variant_raw=sqlite_entry["variant_raw"],
                sqlite_variant_norm=sqlite_entry["variant_norm"],
                match_type="none",
                match_score=None,
                excel_carriers_total=None,
                sqlite_carriers_total=sqlite_entry["carriers_total"] or None,
                carriers_diff=None,
                excel_affected=None,
                sqlite_affected=sqlite_entry["affected_count"] or None,
                affected_diff=None,
                excel_unaffected=None,
                sqlite_unaffected=sqlite_entry["unaffected_count"] or None,
                unaffected_diff=None,
                missing_in_excel=True,
                sqlite_source_layer=sqlite_entry.get("source_layer"),
            )
        )

    logger.info(f"Comparison complete: {len(results)} total rows")
    return results


def create_comparison_row(
    excel_entry: Dict[str, Any],
    sqlite_entry: Optional[Dict[str, Any]],
    match_type: str,
    match_score: Optional[float],
) -> ComparisonRow:
    """Create a ComparisonRow from Excel and SQLite entries."""

    excel_carriers = excel_entry["carriers_total"] or None
    excel_affected = excel_entry["affected_count"] or None
    excel_unaffected = excel_entry["unaffected_count"] or None

    if sqlite_entry:
        sqlite_carriers = sqlite_entry["carriers_total"] or None
        sqlite_affected = sqlite_entry["affected_count"] or None
        sqlite_unaffected = sqlite_entry["unaffected_count"] or None

        carriers_diff = compute_diff(excel_carriers, sqlite_carriers)
        affected_diff = compute_diff(excel_affected, sqlite_affected)
        unaffected_diff = compute_diff(excel_unaffected, sqlite_unaffected)

        # Determine if there's a count mismatch
        count_mismatch = any(
            [
                carriers_diff is not None and carriers_diff != 0,
                affected_diff is not None and affected_diff != 0,
                unaffected_diff is not None and unaffected_diff != 0,
            ]
        )
    else:
        sqlite_carriers = None
        sqlite_affected = None
        sqlite_unaffected = None
        carriers_diff = None
        affected_diff = None
        unaffected_diff = None
        count_mismatch = False

    phenotypes = excel_entry.get("phenotypes", set())
    phenotype = ", ".join(sorted(phenotypes)) if phenotypes else "ALL"

    return ComparisonRow(
        pmid=excel_entry["pmid"],
        excel_variant_raw=excel_entry["variant_raw"],
        excel_variant_norm=excel_entry["variant_norm"],
        sqlite_variant_raw=sqlite_entry["variant_raw"] if sqlite_entry else None,
        sqlite_variant_norm=sqlite_entry["variant_norm"] if sqlite_entry else None,
        match_type=match_type,
        match_score=match_score,
        excel_carriers_total=excel_carriers,
        sqlite_carriers_total=sqlite_carriers,
        carriers_diff=carriers_diff,
        excel_affected=excel_affected,
        sqlite_affected=sqlite_affected,
        affected_diff=affected_diff,
        excel_unaffected=excel_unaffected,
        sqlite_unaffected=sqlite_unaffected,
        unaffected_diff=unaffected_diff,
        phenotype=phenotype,
        count_mismatch=count_mismatch,
        sqlite_source_layer=sqlite_entry.get("source_layer") if sqlite_entry else None,
    )


def _ratio(numerator: int, denominator: int) -> Optional[float]:
    if denominator == 0:
        return None
    return numerator / denominator


def _counted_value(*values: Optional[int]) -> int:
    return sum(v for v in values if v is not None)


def _row_carrier_count(row: ComparisonRow) -> Optional[int]:
    if row.excel_carriers_total is not None:
        return row.excel_carriers_total
    if row.excel_affected is not None or row.excel_unaffected is not None:
        return _counted_value(row.excel_affected, row.excel_unaffected)
    return None


def compute_recall_summary(results: List[ComparisonRow]) -> Dict[str, Any]:
    """Compute gold-side recall metrics from comparison rows.

    The denominator is the curated/Excel side. Extra SQLite-only rows are
    precision findings and do not contribute to recall denominators.
    """
    gold_rows = [r for r in results if not r.missing_in_excel]
    matched_rows = [
        r for r in gold_rows if not r.missing_in_sqlite and r.match_type != "none"
    ]

    gold_pmids = {r.pmid for r in gold_rows if r.pmid}
    matched_pmids = {r.pmid for r in matched_rows if r.pmid}

    gold_variants = {r.excel_variant_norm for r in gold_rows if r.excel_variant_norm}
    matched_variants = {
        r.excel_variant_norm for r in matched_rows if r.excel_variant_norm
    }

    gold_carriers = _counted_value(*(_row_carrier_count(r) for r in gold_rows))
    matched_carriers = _counted_value(*(_row_carrier_count(r) for r in matched_rows))

    gold_affected = _counted_value(*(r.excel_affected for r in gold_rows))
    matched_affected = _counted_value(*(r.excel_affected for r in matched_rows))

    gold_unaffected = _counted_value(*(r.excel_unaffected for r in gold_rows))
    matched_unaffected = _counted_value(*(r.excel_unaffected for r in matched_rows))

    return {
        "pmids": {
            "matched": len(matched_pmids),
            "gold": len(gold_pmids),
            "recall": _ratio(len(matched_pmids), len(gold_pmids)),
        },
        "variant_rows": {
            "matched": len(matched_rows),
            "gold": len(gold_rows),
            "recall": _ratio(len(matched_rows), len(gold_rows)),
        },
        "unique_variants": {
            "matched": len(matched_variants),
            "gold": len(gold_variants),
            "recall": _ratio(len(matched_variants), len(gold_variants)),
        },
        "patients": {
            "matched": matched_carriers,
            "gold": gold_carriers,
            "recall": _ratio(matched_carriers, gold_carriers),
        },
        "affected": {
            "matched": matched_affected,
            "gold": gold_affected,
            "recall": _ratio(matched_affected, gold_affected),
        },
        "unaffected": {
            "matched": matched_unaffected,
            "gold": gold_unaffected,
            "recall": _ratio(matched_unaffected, gold_unaffected),
        },
    }


def compute_rows_mae(results: List[ComparisonRow]) -> Dict[str, Any]:
    """Compute per-row mean absolute error for count fields on matched rows.

    Only matched (PMID x variant) rows with both gold and extracted counts
    present contribute. Lower is better; target is `MAE -> 0`. Returned
    aggregate is gold-dependent and complements the recall metrics.

    Per-field output: ``{"sum": <int>, "n": <int>, "mae": <float | None>}``.
    """
    matched_rows = [
        r
        for r in results
        if not r.missing_in_excel and not r.missing_in_sqlite and r.match_type != "none"
    ]
    field_pairs = [
        ("carriers", "excel_carriers_total", "sqlite_carriers_total"),
        ("affected", "excel_affected", "sqlite_affected"),
        ("unaffected", "excel_unaffected", "sqlite_unaffected"),
    ]
    out: Dict[str, Any] = {}
    for label, gold_attr, ext_attr in field_pairs:
        total_abs = 0
        n = 0
        for r in matched_rows:
            gv = getattr(r, gold_attr, None)
            sv = getattr(r, ext_attr, None)
            if gv is None or sv is None:
                continue
            try:
                total_abs += abs(int(gv) - int(sv))
                n += 1
            except (TypeError, ValueError):
                continue
        out[label] = {
            "sum_abs_error": total_abs,
            "n_matched": n,
            "mae": (total_abs / n) if n else None,
        }
    return out


def compute_precision_summary(results: List[ComparisonRow]) -> Dict[str, Any]:
    """Compute a gold-PMID-restricted, extra-rows-relative-to-gold rate.

    This is NOT clean precision. The gold standard is a curator-selected
    subset, not a paper-exhaustive enumeration of every variant in a paper, so
    a DB row that has no gold match might still be a true positive the curator
    simply did not record. To keep the denominator judgeable, we restrict to
    PMIDs the gold standard actually curated: an unmatched DB row on a PMID
    that gold never touched cannot be adjudicated and is excluded (in one
    sample, ~81% of unmatched DB rows were on non-gold PMIDs). The same
    gold-PMID-restriction is applied in
    ``scripts/recall_audit/paper_disagreement_report.py``.

    Interpret the result as a false-positive UPPER BOUND on gold-curated
    papers, not as precision: even within gold PMIDs, an "extra" DB row may be
    a real variant the curator omitted. Hence the key name
    ``precision_vs_gold_pmids`` and the ``note`` caveat.

    Computation:
        matched_db          = gold rows that matched a DB row
                              (not missing_in_sqlite, match_type != "none")
        extra_on_gold_pmids = missing_in_excel (DB-only) rows whose PMID is in
                              the gold PMID set
        counted_extra_on_gold_pmids = extra rows on gold PMIDs carrying at least
                                      one extracted count field
        precision_vs_gold_pmids = matched_db / (matched_db + extra_on_gold_pmids)
                                  (None when the denominator is 0)
        precision_vs_counted_gold_pmids = matched_db /
                                          (matched_db + counted_extra_on_gold_pmids)
    """
    gold_rows = [r for r in results if not r.missing_in_excel]
    matched_rows = [
        r for r in gold_rows if not r.missing_in_sqlite and r.match_type != "none"
    ]
    gold_pmids = {r.pmid for r in gold_rows if r.pmid}

    matched_db = len(matched_rows)
    extras_on_gold = [r for r in results if r.missing_in_excel and r.pmid in gold_pmids]
    # Extra DB-only rows count only on PMIDs gold curated; rows on non-gold
    # PMIDs are not judgeable by the gold standard and are excluded.
    extra_on_gold_pmids = len(extras_on_gold)
    counted_extra_on_gold_pmids = sum(
        1
        for r in extras_on_gold
        if any(
            value is not None
            for value in (
                r.sqlite_carriers_total,
                r.sqlite_affected,
                r.sqlite_unaffected,
            )
        )
    )
    by_layer: Dict[str, Dict[str, Any]] = {}
    layers = {
        (normalize_source_layer(r.sqlite_source_layer) or "llm_text")
        for r in [*matched_rows, *extras_on_gold]
    }
    for layer in sorted(layers):
        layer_matched = sum(
            1
            for r in matched_rows
            if (normalize_source_layer(r.sqlite_source_layer) or "llm_text") == layer
        )
        layer_extras = [
            r
            for r in extras_on_gold
            if (normalize_source_layer(r.sqlite_source_layer) or "llm_text") == layer
        ]
        layer_counted = sum(
            1
            for r in layer_extras
            if any(
                value is not None
                for value in (
                    r.sqlite_carriers_total,
                    r.sqlite_affected,
                    r.sqlite_unaffected,
                )
            )
        )
        by_layer[layer] = {
            "matched_db": layer_matched,
            "extra_on_gold_pmids": len(layer_extras),
            "counted_extra_on_gold_pmids": layer_counted,
            "precision_vs_gold_pmids": _ratio(
                layer_matched, layer_matched + len(layer_extras)
            ),
            "precision_vs_counted_gold_pmids": _ratio(
                layer_matched, layer_matched + layer_counted
            ),
        }

    return {
        "matched_db": matched_db,
        "extra_on_gold_pmids": extra_on_gold_pmids,
        "counted_extra_on_gold_pmids": counted_extra_on_gold_pmids,
        "precision_vs_gold_pmids": _ratio(matched_db, matched_db + extra_on_gold_pmids),
        "precision_vs_counted_gold_pmids": _ratio(
            matched_db, matched_db + counted_extra_on_gold_pmids
        ),
        "by_source_layer": by_layer,
        "note": (
            "Upper bound on false-positive rate, restricted to gold-curated "
            "PMIDs. Gold is a curator-selected subset, not paper-exhaustive, "
            "so 'extra' DB rows on gold PMIDs may still be true positives the "
            "curator omitted. counted_extra_on_gold_pmids restricts that "
            "denominator to extra rows carrying extracted counts. These are "
            "extra-rows-relative-to-gold rates, NOT clean precision."
        ),
    }


# =============================================================================
# OUTPUT GENERATION
# =============================================================================


def generate_outputs(
    results: List[ComparisonRow], outdir: Path, excel_path: Path, sqlite_path: Path
) -> Dict[str, Any]:
    """
    Generate all output files.

    Args:
        results: List of ComparisonRow objects
        outdir: Output directory
        excel_path: Path to input Excel file
        sqlite_path: Path to input SQLite file

    Returns:
        Summary statistics
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Convert results to DataFrame
    df = pd.DataFrame([asdict(r) for r in results])

    # Calculate summary statistics
    summary = {
        "input_files": {"excel": str(excel_path), "sqlite": str(sqlite_path)},
        "total_rows": len(results),
        "matched_exact": sum(1 for r in results if r.match_type == "exact"),
        "matched_fuzzy": sum(1 for r in results if r.match_type == "fuzzy"),
        "matched_cdna_bridge": sum(
            1 for r in results if r.match_type.endswith("_cdna_bridge")
        ),
        "unmatched": sum(1 for r in results if r.match_type == "none"),
        "missing_in_sqlite": sum(1 for r in results if r.missing_in_sqlite),
        "missing_in_excel": sum(1 for r in results if r.missing_in_excel),
        "count_mismatches": sum(1 for r in results if r.count_mismatch),
        "unique_pmids": len(set(r.pmid for r in results)),
        "recall": compute_recall_summary(results),
        "mae": compute_rows_mae(results),
        "precision": compute_precision_summary(results),
        "top_mismatches": [],
    }

    # Get top mismatches (by absolute diff)
    mismatched = [r for r in results if r.count_mismatch]
    mismatched.sort(
        key=lambda r: abs(r.carriers_diff or 0) + abs(r.affected_diff or 0),
        reverse=True,
    )
    summary["top_mismatches"] = [asdict(r) for r in mismatched[:20]]

    # Write full comparison CSV
    df.to_csv(outdir / "discrepancies.csv", index=False)
    logger.info(f"Wrote {outdir / 'discrepancies.csv'}")

    # Write missing_in_sqlite.csv
    missing_sqlite = df[df["missing_in_sqlite"]]
    if len(missing_sqlite) > 0:
        missing_sqlite.to_csv(outdir / "missing_in_sqlite.csv", index=False)
        logger.info(
            f"Wrote {outdir / 'missing_in_sqlite.csv'} ({len(missing_sqlite)} rows)"
        )

    # Write missing_in_excel.csv
    missing_excel = df[df["missing_in_excel"]]
    if len(missing_excel) > 0:
        missing_excel.to_csv(outdir / "missing_in_excel.csv", index=False)
        logger.info(
            f"Wrote {outdir / 'missing_in_excel.csv'} ({len(missing_excel)} rows)"
        )

        # Write the gold-PMID-restricted subset of unmatched DB rows for
        # manual adjudication. These are the rows that feed the
        # precision_vs_gold_pmids denominator (extra rows on gold PMIDs).
        gold_pmids = {r.pmid for r in results if not r.missing_in_excel and r.pmid}
        extra_on_gold = missing_excel[missing_excel["pmid"].isin(gold_pmids)]
        if len(extra_on_gold) > 0:
            extra_on_gold.to_csv(
                outdir / "unmatched_db_rows_on_gold_pmids.csv", index=False
            )
            logger.info(
                f"Wrote {outdir / 'unmatched_db_rows_on_gold_pmids.csv'} "
                f"({len(extra_on_gold)} rows)"
            )

    # Write summary.json
    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Wrote {outdir / 'summary.json'}")

    # Write markdown report
    write_markdown_report(results, summary, outdir / "report.md")
    logger.info(f"Wrote {outdir / 'report.md'}")

    return summary


def write_markdown_report(
    results: List[ComparisonRow], summary: Dict[str, Any], output_path: Path
) -> None:
    """Write a human-readable markdown report."""

    def fmt_recall(metric: Dict[str, Any]) -> str:
        recall = metric.get("recall")
        recall_text = "n/a" if recall is None else f"{recall:.1%}"
        return f"{metric.get('matched', 0)}/{metric.get('gold', 0)} ({recall_text})"

    lines = [
        "# Variant Comparison Report",
        "",
        "## Overview",
        "",
        f"- **Excel file**: `{summary['input_files']['excel']}`",
        f"- **SQLite file**: `{summary['input_files']['sqlite']}`",
        f"- **Total comparisons**: {summary['total_rows']}",
        f"- **Unique PMIDs**: {summary['unique_pmids']}",
        "",
        "## Match Summary",
        "",
        "| Metric | Count |",
        "|--------|-------|",
        f"| Exact matches | {summary['matched_exact']} |",
        f"| Fuzzy matches | {summary['matched_fuzzy']} |",
        f"| cDNA/protein bridge matches | {summary.get('matched_cdna_bridge', 0)} |",
        f"| Unmatched | {summary['unmatched']} |",
        f"| Missing in SQLite | {summary['missing_in_sqlite']} |",
        f"| Missing in Excel | {summary['missing_in_excel']} |",
        f"| Count mismatches | {summary['count_mismatches']} |",
        "",
    ]

    recall = summary.get("recall", {})
    if recall:
        lines.extend(
            [
                "## Recall",
                "",
                "| Dimension | Matched / Gold |",
                "|-----------|----------------|",
                f"| PMIDs | {fmt_recall(recall.get('pmids', {}))} |",
                f"| Variant rows | {fmt_recall(recall.get('variant_rows', {}))} |",
                f"| Unique variants | {fmt_recall(recall.get('unique_variants', {}))} |",
                f"| Patients/carriers | {fmt_recall(recall.get('patients', {}))} |",
                f"| Affected | {fmt_recall(recall.get('affected', {}))} |",
                f"| Unaffected | {fmt_recall(recall.get('unaffected', {}))} |",
                "",
            ]
        )

    precision = summary.get("precision", {})
    if precision:
        pvg = precision.get("precision_vs_gold_pmids")
        pvg_text = "n/a" if pvg is None else f"{pvg:.1%}"
        pcg = precision.get("precision_vs_counted_gold_pmids")
        pcg_text = "n/a" if pcg is None else f"{pcg:.1%}"
        lines.extend(
            [
                "## Precision (counted extras vs gold PMIDs)",
                "",
                "Headline precision uses only count-bearing extra rows on "
                "gold-curated PMIDs. The raw gold-PMID rate is a loose "
                "false-positive upper bound dominated by zero-count variant "
                "mentions.",
                "",
                f"- Matched DB rows: {precision.get('matched_db', 0)}",
                f"- Counted extra DB rows on gold PMIDs: "
                f"{precision.get('counted_extra_on_gold_pmids', 0)}",
                f"- precision_vs_counted_gold_pmids: {pcg_text}",
                f"- Loose extra DB rows on gold PMIDs: "
                f"{precision.get('extra_on_gold_pmids', 0)}",
                f"- loose precision_vs_gold_pmids: {pvg_text}",
                "",
            ]
        )
        by_layer = precision.get("by_source_layer") or {}
        if by_layer:
            lines.extend(
                [
                    "| Source layer | Matched DB rows | Extra rows | Counted extra rows | precision_vs_gold_pmids | precision_vs_counted_gold_pmids |",
                    "|--------------|-----------------|------------|--------------------|-------------------------|--------------------------------|",
                ]
            )
            for layer, block in sorted(by_layer.items()):
                layer_p = block.get("precision_vs_gold_pmids")
                layer_pc = block.get("precision_vs_counted_gold_pmids")
                lines.append(
                    f"| {layer} | {block.get('matched_db', 0)} | "
                    f"{block.get('extra_on_gold_pmids', 0)} | "
                    f"{block.get('counted_extra_on_gold_pmids', 0)} | "
                    f"{'n/a' if layer_p is None else f'{layer_p:.1%}'} | "
                    f"{'n/a' if layer_pc is None else f'{layer_pc:.1%}'} |"
                )
            lines.append("")

    # Top mismatches table
    if summary["count_mismatches"] > 0:
        lines.extend(
            [
                "## Top Count Mismatches",
                "",
                "| PMID | Excel Variant | SQLite Variant | Carriers (E/S/Δ) | Affected (E/S/Δ) |",
                "|------|---------------|----------------|------------------|------------------|",
            ]
        )

        mismatched = [r for r in results if r.count_mismatch]
        mismatched.sort(
            key=lambda r: abs(r.carriers_diff or 0) + abs(r.affected_diff or 0),
            reverse=True,
        )

        for r in mismatched[:50]:
            carriers = f"{r.excel_carriers_total or '-'}/{r.sqlite_carriers_total or '-'}/{r.carriers_diff or '-'}"
            affected = f"{r.excel_affected or '-'}/{r.sqlite_affected or '-'}/{r.affected_diff or '-'}"
            lines.append(
                f"| {r.pmid} | {r.excel_variant_raw[:30]} | {r.sqlite_variant_raw or '-'} | {carriers} | {affected} |"
            )

        lines.append("")

    # Missing in SQLite
    missing_sqlite = [r for r in results if r.missing_in_sqlite]
    if missing_sqlite:
        lines.extend(
            [
                "## Variants Missing in SQLite (Top 50)",
                "",
                "| PMID | Excel Variant |",
                "|------|---------------|",
            ]
        )
        for r in missing_sqlite[:50]:
            lines.append(f"| {r.pmid} | {r.excel_variant_raw} |")
        lines.append("")

    # Missing in Excel
    missing_excel = [r for r in results if r.missing_in_excel]
    if missing_excel:
        lines.extend(
            [
                "## Variants Missing in Excel (Top 50)",
                "",
                "| PMID | SQLite Variant |",
                "|------|----------------|",
            ]
        )
        for r in missing_excel[:50]:
            lines.append(f"| {r.pmid} | {r.sqlite_variant_raw} |")
        lines.append("")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))


# =============================================================================
# CLI ENTRYPOINT
# =============================================================================


def load_mapping_file(mapping_path: Path) -> Dict[str, str]:
    """Load column mapping from JSON or YAML file."""
    with open(mapping_path, "r") as f:
        content = f.read()

    # Try JSON first
    try:
        return json.loads(content)
    except json.JSONDecodeError:
        pass

    # Try YAML
    if YAML_AVAILABLE:
        try:
            return yaml.safe_load(content)
        except yaml.YAMLError:
            pass

    raise ValueError(f"Cannot parse mapping file: {mapping_path}")


def main():
    """Main CLI entrypoint."""
    parser = argparse.ArgumentParser(
        description="Compare curated gold-standard variants against a GVF SQLite database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic comparison
  python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db

  # Normalized recall CSVs are also supported
  python compare_variants.py --excel gene_variant_fetcher_gold_standard/normalized/KCNQ1_recall_input.csv \\
      --sqlite KCNQ1.db

  # Fuzzy matching with lower threshold
  python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db \\
      --variant_match_mode fuzzy --fuzzy_threshold 0.80

  # Specify sheet and output directory
  python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db \\
      --sheet "Variants" --outdir ./results

  # Use explicit column mapping
  python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db \\
      --mapping column_map.yaml
""",
    )

    parser.add_argument(
        "--excel",
        "-e",
        type=Path,
        required=True,
        help="Path to Excel/CSV file with curated variant data",
    )

    parser.add_argument(
        "--sqlite",
        "-s",
        type=Path,
        required=True,
        help="Path to SQLite database (e.g., KCNH2.db)",
    )

    parser.add_argument(
        "--sheet",
        type=str,
        default=None,
        help="Excel sheet name (default: first sheet)",
    )

    parser.add_argument(
        "--outdir",
        "-o",
        type=Path,
        default=Path("./compare_out"),
        help="Output directory (default: ./compare_out)",
    )

    parser.add_argument(
        "--mapping",
        "-m",
        type=Path,
        default=None,
        help="JSON or YAML file with explicit column mappings",
    )

    parser.add_argument(
        "--variant_match_mode",
        choices=["exact", "fuzzy"],
        default="fuzzy",
        help="Variant matching mode (default: fuzzy for better cross-notation matching)",
    )

    parser.add_argument(
        "--fuzzy_threshold",
        type=float,
        default=0.80,
        help="Minimum similarity for fuzzy matching (default: 0.80 for better cross-notation matching)",
    )

    parser.add_argument(
        "--log_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)",
    )

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Validate inputs
    if not args.excel.exists():
        logger.error(f"Excel file not found: {args.excel}")
        return 1

    if not args.sqlite.exists():
        logger.error(f"SQLite file not found: {args.sqlite}")
        return 1

    # Load optional mapping
    column_mapping = None
    if args.mapping:
        if not args.mapping.exists():
            logger.error(f"Mapping file not found: {args.mapping}")
            return 1
        column_mapping = load_mapping_file(args.mapping)
        logger.info(f"Loaded column mapping: {column_mapping}")

    logger.info("=" * 80)
    logger.info("VARIANT COMPARISON TOOL")
    logger.info("=" * 80)
    logger.info(f"Excel: {args.excel}")
    logger.info(f"SQLite: {args.sqlite}")
    logger.info(f"Match mode: {args.variant_match_mode}")
    if args.variant_match_mode == "fuzzy":
        logger.info(f"Fuzzy threshold: {args.fuzzy_threshold}")
    logger.info("=" * 80)

    try:
        # Load Excel data
        excel_df, detected_columns = load_excel_data(
            args.excel, args.sheet, column_mapping
        )

        # Connect to SQLite
        conn = sqlite3.connect(f"file:{args.sqlite}?mode=ro", uri=True)
        logger.info(f"Connected to SQLite database: {args.sqlite}")

        # Introspect SQLite schema
        table_info = introspect_sqlite(conn)

        # Extract SQLite data
        sqlite_df = extract_sqlite_data(conn, table_info)
        conn.close()

        # Aggregate data
        excel_aggregated = aggregate_excel_data(excel_df, detected_columns)
        sqlite_aggregated = aggregate_sqlite_data(sqlite_df)

        # Compare
        results = compare_data(
            excel_aggregated,
            sqlite_aggregated,
            args.variant_match_mode,
            args.fuzzy_threshold,
        )

        # Generate outputs
        summary = generate_outputs(results, args.outdir, args.excel, args.sqlite)

        # Print summary
        logger.info("=" * 80)
        logger.info("COMPARISON COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Total rows: {summary['total_rows']}")
        logger.info(f"Exact matches: {summary['matched_exact']}")
        logger.info(f"Fuzzy matches: {summary['matched_fuzzy']}")
        logger.info(f"Missing in SQLite: {summary['missing_in_sqlite']}")
        logger.info(f"Missing in Excel: {summary['missing_in_excel']}")
        logger.info(f"Count mismatches: {summary['count_mismatches']}")
        logger.info(f"Output directory: {args.outdir}")
        logger.info("=" * 80)

        return 0

    except Exception as e:
        logger.error(f"Comparison failed: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    import sys

    sys.exit(main())
