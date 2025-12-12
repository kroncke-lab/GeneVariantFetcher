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
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Set, Tuple

import pandas as pd

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
    'pmid': [
        'pmid', 'pubmed_id', 'pubmed', 'pm_id', 'pub_id', 'pubmedid',
        'pmid_id', 'paper_id', 'article_id', 'reference_id'
    ],
    'variant': [
        'hgvs', 'hgvs_p', 'hgvs_c', 'protein_change', 'cdna_change',
        'variant', 'aa_change', 'mutation', 'amino_acid', 'protein_variant',
        'cdna_variant', 'nucleotide_change', 'dna_change', 'protein_mutation',
        'variant_name', 'variant_id', 'hgvs_protein', 'hgvs_cdna', 'p_notation',
        'c_notation', 'protein_notation', 'cdna_notation'
    ],
    'rsid': [
        'rsid', 'rs_id', 'dbsnp', 'rs', 'rs_number', 'dbsnp_id', 'snp_id'
    ],
    'carriers_total': [
        'carriers', 'n_carriers', 'carrier_count', 'total_carriers',
        'n_total', 'subjects', 'probands', 'individuals', 'n_individuals',
        'total', 'n', 'count', 'sample_size', 'n_subjects', 'carrier_n',
        'total_n', 'num_carriers', 'number_carriers'
    ],
    'affected_count': [
        'affected', 'cases', 'n_affected', 'symptomatic', 'patients',
        'n_cases', 'case_count', 'affected_count', 'n_symptomatic',
        'diseased', 'n_diseased', 'num_affected', 'number_affected'
    ],
    'unaffected_count': [
        'unaffected', 'controls', 'n_unaffected', 'asymptomatic',
        'n_controls', 'control_count', 'unaffected_count', 'n_asymptomatic',
        'healthy', 'n_healthy', 'num_unaffected', 'number_unaffected'
    ],
    'phenotype': [
        'phenotype', 'condition', 'disease', 'diagnosis', 'clinical',
        'clinical_phenotype', 'phenotype_details', 'disorder', 'syndrome'
    ],
    'notes': [
        'notes', 'comments', 'remarks', 'additional_info', 'description',
        'free_text', 'annotation', 'additional_notes'
    ]
}


def normalize_column_name(col: str) -> str:
    """Normalize column name for matching (lowercase, strip, remove underscores/spaces)."""
    return re.sub(r'[\s_\-]+', '', col.lower().strip())


def detect_columns(df: pd.DataFrame, mapping: Optional[Dict[str, str]] = None) -> Dict[str, Optional[str]]:
    """
    Detect column mappings from DataFrame using synonym matching.

    Args:
        df: Input DataFrame
        mapping: Optional explicit column mapping from user

    Returns:
        Dict mapping semantic field names to actual column names
    """
    detected = {field: None for field in COLUMN_SYNONYMS.keys()}

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
                    if norm_syn in norm_col or norm_col in norm_syn:
                        detected[field] = orig_col
                        logger.info(f"Partial match column: {field} -> {orig_col}")
                        break
                if detected[field] is not None:
                    break

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

    pmid_patterns = ['pmid', 'pubmed', 'paper_id', 'article_id']
    variant_patterns = [
        'protein_notation', 'cdna_notation', 'hgvs', 'variant',
        'mutation', 'genomic_position', 'protein_change'
    ]
    count_patterns = {
        'carriers_total': ['total_carriers', 'carriers', 'total_carriers_observed', 'n_total'],
        'affected_count': ['affected_count', 'affected', 'n_affected', 'cases'],
        'unaffected_count': ['unaffected_count', 'unaffected', 'n_unaffected', 'controls'],
        'uncertain_count': ['uncertain_count', 'uncertain', 'n_uncertain']
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
        logger.debug(f"Table {table}: pmid={info.has_pmid}, variant={info.has_variant}, counts={info.has_counts}")

    return table_info


def find_best_data_source(table_info: Dict[str, TableInfo]) -> Tuple[str, TableInfo]:
    """
    Find the best table for extracting variant-paper-count data.

    Priority:
    1. Tables with pmid + variant + counts (penetrance_data joined with variants)
    2. Tables with pmid + variant (variant_papers)
    3. Individual records (for aggregation)

    Args:
        table_info: Dict of TableInfo from introspection

    Returns:
        Tuple of (strategy_name, primary_table_info)
    """
    # Check for penetrance_data (ideal case)
    if 'penetrance_data' in table_info:
        info = table_info['penetrance_data']
        if info.has_pmid and info.has_counts:
            logger.info("Using penetrance_data table as primary source")
            return ('penetrance_data', info)

    # Check for tables with all three: pmid, variant, counts
    for name, info in table_info.items():
        if info.has_pmid and info.has_variant and info.has_counts:
            logger.info(f"Using table {name} with pmid+variant+counts")
            return ('combined', info)

    # Check for individual_records (aggregate affected_status)
    if 'individual_records' in table_info:
        info = table_info['individual_records']
        if info.has_pmid:
            logger.info("Using individual_records table (will aggregate counts)")
            return ('individual_records', info)

    # Fallback: variant_papers + variants
    if 'variant_papers' in table_info and 'variants' in table_info:
        logger.info("Using variant_papers + variants join strategy")
        return ('variant_papers_join', table_info['variant_papers'])

    # No suitable source found
    raise ValueError(
        "Cannot find suitable tables for comparison.\n"
        f"Available tables: {list(table_info.keys())}\n"
        "Expected: penetrance_data, individual_records, or variant_papers with variants\n"
        "Consider using --mapping to specify column mappings."
    )


def extract_sqlite_data(conn: sqlite3.Connection, table_info: Dict[str, TableInfo]) -> pd.DataFrame:
    """
    Extract variant data from SQLite database.

    Args:
        conn: SQLite connection
        table_info: Schema information from introspection

    Returns:
        DataFrame with columns: pmid, variant, carriers_total, affected_count, unaffected_count
    """
    strategy, primary_table = find_best_data_source(table_info)

    if strategy == 'penetrance_data':
        # Join penetrance_data with variants to get variant notation
        query = """
            SELECT
                pd.pmid,
                COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position) as variant,
                v.protein_notation,
                v.cdna_notation,
                pd.total_carriers_observed as carriers_total,
                pd.affected_count,
                pd.unaffected_count,
                pd.uncertain_count
            FROM penetrance_data pd
            JOIN variants v ON pd.variant_id = v.variant_id
        """
        logger.info("Executing penetrance_data query")
        df = pd.read_sql_query(query, conn)

    elif strategy == 'individual_records':
        # Aggregate individual records by variant+pmid
        query = """
            SELECT
                ir.pmid,
                COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position) as variant,
                v.protein_notation,
                v.cdna_notation,
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

    elif strategy == 'variant_papers_join':
        # Just get variant-paper associations (no counts available)
        query = """
            SELECT
                vp.pmid,
                COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position) as variant,
                v.protein_notation,
                v.cdna_notation,
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

        variant_select = f"COALESCE({', '.join(variant_cols)})" if variant_cols else "NULL"

        query = f"""
            SELECT
                {pmid_col} as pmid,
                {variant_select} as variant,
                {count_cols.get('carriers_total', 'NULL')} as carriers_total,
                {count_cols.get('affected_count', 'NULL')} as affected_count,
                {count_cols.get('unaffected_count', 'NULL')} as unaffected_count,
                {count_cols.get('uncertain_count', 'NULL')} as uncertain_count
            FROM {primary_table.name}
        """
        logger.info(f"Executing custom query on {primary_table.name}")
        df = pd.read_sql_query(query, conn)

    # Ensure pmid is string
    df['pmid'] = df['pmid'].astype(str)

    logger.info(f"Extracted {len(df)} records from SQLite")
    return df


# =============================================================================
# VARIANT NORMALIZATION
# =============================================================================

# Amino acid 3-letter to 1-letter mapping
AA_3_TO_1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*', 'X': '*'
}

AA_1_TO_3 = {v: k for k, v in AA_3_TO_1.items() if k != 'X'}
AA_1_TO_3['*'] = 'Ter'


def normalize_unicode(s: str) -> str:
    """Normalize unicode characters (dashes, quotes, etc.)."""
    if not s:
        return s
    # Normalize to NFKC form
    s = unicodedata.normalize('NFKC', s)
    # Replace various dash characters with standard hyphen
    s = re.sub(r'[\u2010\u2011\u2012\u2013\u2014\u2015\u2212]', '-', s)
    # Replace fancy quotes
    s = re.sub(r'[\u2018\u2019]', "'", s)
    s = re.sub(r'[\u201c\u201d]', '"', s)
    return s


def normalize_variant(variant: str) -> str:
    """
    Normalize a variant string for comparison.

    - Strips whitespace
    - Normalizes unicode
    - Collapses multiple spaces
    - Standardizes case for common prefixes
    """
    if not variant or pd.isna(variant):
        return ''

    variant = str(variant).strip()
    variant = normalize_unicode(variant)
    variant = re.sub(r'\s+', ' ', variant)

    # Standardize common prefixes
    variant = re.sub(r'^p\.\s*', 'p.', variant, flags=re.IGNORECASE)
    variant = re.sub(r'^c\.\s*', 'c.', variant, flags=re.IGNORECASE)

    return variant


def convert_aa_3_to_1(variant: str) -> Optional[str]:
    """
    Convert 3-letter amino acid codes to 1-letter in a variant string.

    e.g., p.Arg123His -> p.R123H
    """
    if not variant or not variant.lower().startswith('p.'):
        return None

    result = variant[:2]  # Keep "p."
    remaining = variant[2:]

    # Pattern: 3-letter AA, position, 3-letter AA (with optional fs, del, etc.)
    pattern = r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})?(.*)$'
    match = re.match(pattern, remaining)

    if match:
        aa1, pos, aa2, rest = match.groups()
        aa1_short = AA_3_TO_1.get(aa1, aa1)
        aa2_short = AA_3_TO_1.get(aa2, aa2) if aa2 else ''
        result = f"p.{aa1_short}{pos}{aa2_short}{rest or ''}"
        return result

    return None


def convert_aa_1_to_3(variant: str) -> Optional[str]:
    """
    Convert 1-letter amino acid codes to 3-letter in a variant string.

    e.g., p.R123H -> p.Arg123His
    """
    if not variant or not variant.lower().startswith('p.'):
        return None

    # Pattern: single letter AA, position, single letter AA (with optional suffix)
    pattern = r'^p\.([A-Z\*])(\d+)([A-Z\*])?(.*)$'
    match = re.match(pattern, variant)

    if match:
        aa1, pos, aa2, rest = match.groups()
        aa1_long = AA_1_TO_3.get(aa1, aa1)
        aa2_long = AA_1_TO_3.get(aa2, aa2) if aa2 else ''
        result = f"p.{aa1_long}{pos}{aa2_long}{rest or ''}"
        return result

    return None


def get_variant_forms(variant: str) -> Set[str]:
    """
    Get all equivalent forms of a variant for matching.

    Returns set of normalized variant strings including:
    - Original normalized
    - 3-letter to 1-letter conversion
    - 1-letter to 3-letter conversion
    """
    forms = set()

    normalized = normalize_variant(variant)
    if normalized:
        forms.add(normalized)

    # Try amino acid conversions
    converted_1 = convert_aa_3_to_1(normalized)
    if converted_1:
        forms.add(normalize_variant(converted_1))

    converted_3 = convert_aa_1_to_3(normalized)
    if converted_3:
        forms.add(normalize_variant(converted_3))

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


def find_best_match(
    query_variant: str,
    candidates: List[str],
    threshold: float = 0.85
) -> Tuple[Optional[str], float, str]:
    """
    Find the best matching variant from candidates.

    Args:
        query_variant: Variant to match
        candidates: List of candidate variants
        threshold: Minimum similarity threshold

    Returns:
        Tuple of (best_match, score, match_type)
    """
    query_norm = normalize_variant(query_variant)
    query_forms = get_variant_forms(query_variant)

    # First try exact match on any form
    for candidate in candidates:
        candidate_norm = normalize_variant(candidate)
        candidate_forms = get_variant_forms(candidate)

        if query_forms & candidate_forms:  # Intersection
            return (candidate, 1.0, 'exact')

    # Fuzzy matching
    best_match = None
    best_score = 0.0

    for candidate in candidates:
        candidate_norm = normalize_variant(candidate)

        # Compare all forms
        for q_form in query_forms:
            for c_form in get_variant_forms(candidate):
                score = compute_similarity(q_form, c_form)
                if score > best_score:
                    best_score = score
                    best_match = candidate

    if best_score >= threshold:
        return (best_match, best_score, 'fuzzy')

    return (None, best_score, 'none')


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

    phenotype: str = 'ALL'

    missing_in_sqlite: bool = False
    missing_in_excel: bool = False
    count_mismatch: bool = False


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
    column_mapping: Optional[Dict[str, str]]
) -> Tuple[pd.DataFrame, Dict[str, Optional[str]]]:
    """
    Load and parse Excel data.

    Returns:
        Tuple of (DataFrame, detected_columns)
    """
    logger.info(f"Loading Excel file: {excel_path}")

    df = pd.read_excel(excel_path, sheet_name=sheet_name, engine='openpyxl')
    logger.info(f"Loaded {len(df)} rows from Excel")
    logger.info(f"Columns: {list(df.columns)}")

    detected = detect_columns(df, column_mapping)

    # Validate required columns
    if detected['pmid'] is None:
        raise ValueError(
            f"Cannot find PMID column. Available columns: {list(df.columns)}\n"
            "Use --mapping to specify column mappings."
        )

    if detected['variant'] is None:
        raise ValueError(
            f"Cannot find variant column. Available columns: {list(df.columns)}\n"
            "Use --mapping to specify column mappings."
        )

    return df, detected


def aggregate_excel_data(
    df: pd.DataFrame,
    detected: Dict[str, Optional[str]]
) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """
    Aggregate Excel data by (pmid, variant) key.

    Returns:
        Dict mapping (pmid, variant_norm) -> aggregated data
    """
    pmid_col = detected['pmid']
    variant_col = detected['variant']
    carriers_col = detected.get('carriers_total')
    affected_col = detected.get('affected_count')
    unaffected_col = detected.get('unaffected_count')
    phenotype_col = detected.get('phenotype')

    aggregated = {}

    for _, row in df.iterrows():
        pmid = str(row[pmid_col]).strip()
        variant_raw = str(row[variant_col]) if pd.notna(row[variant_col]) else ''
        variant_norm = normalize_variant(variant_raw)

        if not pmid or not variant_norm:
            continue

        key = (pmid, variant_norm)

        if key not in aggregated:
            aggregated[key] = {
                'pmid': pmid,
                'variant_raw': variant_raw,
                'variant_norm': variant_norm,
                'carriers_total': 0,
                'affected_count': 0,
                'unaffected_count': 0,
                'phenotypes': set()
            }

        # Aggregate counts
        if carriers_col and pd.notna(row.get(carriers_col)):
            aggregated[key]['carriers_total'] += safe_int(row[carriers_col]) or 0

        if affected_col and pd.notna(row.get(affected_col)):
            aggregated[key]['affected_count'] += safe_int(row[affected_col]) or 0

        if unaffected_col and pd.notna(row.get(unaffected_col)):
            aggregated[key]['unaffected_count'] += safe_int(row[unaffected_col]) or 0

        if phenotype_col and pd.notna(row.get(phenotype_col)):
            aggregated[key]['phenotypes'].add(str(row[phenotype_col]))

    logger.info(f"Aggregated Excel data into {len(aggregated)} unique (pmid, variant) pairs")
    return aggregated


def aggregate_sqlite_data(df: pd.DataFrame) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """
    Aggregate SQLite data by (pmid, variant) key.

    Returns:
        Dict mapping (pmid, variant_norm) -> aggregated data
    """
    aggregated = {}

    for _, row in df.iterrows():
        pmid = str(row['pmid']).strip()
        variant_raw = str(row['variant']) if pd.notna(row['variant']) else ''
        variant_norm = normalize_variant(variant_raw)

        if not pmid or not variant_norm:
            continue

        key = (pmid, variant_norm)

        if key not in aggregated:
            aggregated[key] = {
                'pmid': pmid,
                'variant_raw': variant_raw,
                'variant_norm': variant_norm,
                'carriers_total': 0,
                'affected_count': 0,
                'unaffected_count': 0,
                'protein_notation': row.get('protein_notation'),
                'cdna_notation': row.get('cdna_notation')
            }

        # Aggregate counts (handle None values)
        aggregated[key]['carriers_total'] += safe_int(row.get('carriers_total')) or 0
        aggregated[key]['affected_count'] += safe_int(row.get('affected_count')) or 0
        aggregated[key]['unaffected_count'] += safe_int(row.get('unaffected_count')) or 0

    logger.info(f"Aggregated SQLite data into {len(aggregated)} unique (pmid, variant) pairs")
    return aggregated


def compare_data(
    excel_data: Dict[Tuple[str, str], Dict[str, Any]],
    sqlite_data: Dict[Tuple[str, str], Dict[str, Any]],
    match_mode: str,
    fuzzy_threshold: float
) -> List[ComparisonRow]:
    """
    Compare Excel and SQLite data.

    Args:
        excel_data: Aggregated Excel data
        sqlite_data: Aggregated SQLite data
        match_mode: 'exact' or 'fuzzy'
        fuzzy_threshold: Threshold for fuzzy matching

    Returns:
        List of ComparisonRow objects
    """
    results = []
    matched_sqlite_keys = set()

    # Group SQLite variants by PMID for matching
    sqlite_by_pmid: Dict[str, List[Tuple[str, str]]] = {}
    for (pmid, variant_norm), data in sqlite_data.items():
        if pmid not in sqlite_by_pmid:
            sqlite_by_pmid[pmid] = []
        sqlite_by_pmid[pmid].append((variant_norm, data['variant_raw']))

    # Process Excel entries
    for (excel_pmid, excel_variant_norm), excel_entry in excel_data.items():
        # Get SQLite candidates for this PMID
        sqlite_candidates = sqlite_by_pmid.get(excel_pmid, [])
        candidate_variants = [raw for (norm, raw) in sqlite_candidates]

        # Try to match
        if match_mode == 'exact':
            # Exact match only
            match_key = (excel_pmid, excel_variant_norm)
            if match_key in sqlite_data:
                sqlite_entry = sqlite_data[match_key]
                matched_sqlite_keys.add(match_key)

                row = create_comparison_row(
                    excel_entry, sqlite_entry,
                    'exact', 1.0
                )
            else:
                row = create_comparison_row(
                    excel_entry, None,
                    'none', None
                )
                row.missing_in_sqlite = True
        else:
            # Fuzzy matching
            best_match, score, match_type = find_best_match(
                excel_entry['variant_raw'],
                candidate_variants,
                fuzzy_threshold
            )

            if best_match:
                # Find the sqlite entry for this match
                best_norm = normalize_variant(best_match)
                match_key = (excel_pmid, best_norm)

                # Try exact key first, then search
                if match_key in sqlite_data:
                    sqlite_entry = sqlite_data[match_key]
                else:
                    # Search for matching variant in this PMID
                    sqlite_entry = None
                    for (pmid, norm), entry in sqlite_data.items():
                        if pmid == excel_pmid and entry['variant_raw'] == best_match:
                            sqlite_entry = entry
                            match_key = (pmid, norm)
                            break

                if sqlite_entry:
                    matched_sqlite_keys.add(match_key)
                    row = create_comparison_row(
                        excel_entry, sqlite_entry,
                        match_type, score
                    )
                else:
                    row = create_comparison_row(
                        excel_entry, None,
                        'none', None
                    )
                    row.missing_in_sqlite = True
            else:
                row = create_comparison_row(
                    excel_entry, None,
                    'none', None
                )
                row.missing_in_sqlite = True

        results.append(row)

    # Add SQLite entries not matched
    for key, sqlite_entry in sqlite_data.items():
        if key not in matched_sqlite_keys:
            row = ComparisonRow(
                pmid=sqlite_entry['pmid'],
                excel_variant_raw='',
                excel_variant_norm='',
                sqlite_variant_raw=sqlite_entry['variant_raw'],
                sqlite_variant_norm=sqlite_entry['variant_norm'],
                match_type='none',
                match_score=None,
                excel_carriers_total=None,
                sqlite_carriers_total=sqlite_entry['carriers_total'] or None,
                carriers_diff=None,
                excel_affected=None,
                sqlite_affected=sqlite_entry['affected_count'] or None,
                affected_diff=None,
                excel_unaffected=None,
                sqlite_unaffected=sqlite_entry['unaffected_count'] or None,
                unaffected_diff=None,
                missing_in_excel=True
            )
            results.append(row)

    logger.info(f"Comparison complete: {len(results)} total rows")
    return results


def create_comparison_row(
    excel_entry: Dict[str, Any],
    sqlite_entry: Optional[Dict[str, Any]],
    match_type: str,
    match_score: Optional[float]
) -> ComparisonRow:
    """Create a ComparisonRow from Excel and SQLite entries."""

    excel_carriers = excel_entry['carriers_total'] or None
    excel_affected = excel_entry['affected_count'] or None
    excel_unaffected = excel_entry['unaffected_count'] or None

    if sqlite_entry:
        sqlite_carriers = sqlite_entry['carriers_total'] or None
        sqlite_affected = sqlite_entry['affected_count'] or None
        sqlite_unaffected = sqlite_entry['unaffected_count'] or None

        carriers_diff = compute_diff(excel_carriers, sqlite_carriers)
        affected_diff = compute_diff(excel_affected, sqlite_affected)
        unaffected_diff = compute_diff(excel_unaffected, sqlite_unaffected)

        # Determine if there's a count mismatch
        count_mismatch = any([
            carriers_diff is not None and carriers_diff != 0,
            affected_diff is not None and affected_diff != 0,
            unaffected_diff is not None and unaffected_diff != 0
        ])
    else:
        sqlite_carriers = None
        sqlite_affected = None
        sqlite_unaffected = None
        carriers_diff = None
        affected_diff = None
        unaffected_diff = None
        count_mismatch = False

    phenotypes = excel_entry.get('phenotypes', set())
    phenotype = ', '.join(sorted(phenotypes)) if phenotypes else 'ALL'

    return ComparisonRow(
        pmid=excel_entry['pmid'],
        excel_variant_raw=excel_entry['variant_raw'],
        excel_variant_norm=excel_entry['variant_norm'],
        sqlite_variant_raw=sqlite_entry['variant_raw'] if sqlite_entry else None,
        sqlite_variant_norm=sqlite_entry['variant_norm'] if sqlite_entry else None,
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
        count_mismatch=count_mismatch
    )


# =============================================================================
# OUTPUT GENERATION
# =============================================================================

def generate_outputs(
    results: List[ComparisonRow],
    outdir: Path,
    excel_path: Path,
    sqlite_path: Path
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
        'input_files': {
            'excel': str(excel_path),
            'sqlite': str(sqlite_path)
        },
        'total_rows': len(results),
        'matched_exact': sum(1 for r in results if r.match_type == 'exact'),
        'matched_fuzzy': sum(1 for r in results if r.match_type == 'fuzzy'),
        'unmatched': sum(1 for r in results if r.match_type == 'none'),
        'missing_in_sqlite': sum(1 for r in results if r.missing_in_sqlite),
        'missing_in_excel': sum(1 for r in results if r.missing_in_excel),
        'count_mismatches': sum(1 for r in results if r.count_mismatch),
        'unique_pmids': len(set(r.pmid for r in results)),
        'top_mismatches': []
    }

    # Get top mismatches (by absolute diff)
    mismatched = [r for r in results if r.count_mismatch]
    mismatched.sort(key=lambda r: abs(r.carriers_diff or 0) + abs(r.affected_diff or 0), reverse=True)
    summary['top_mismatches'] = [asdict(r) for r in mismatched[:20]]

    # Write full comparison CSV
    df.to_csv(outdir / 'discrepancies.csv', index=False)
    logger.info(f"Wrote {outdir / 'discrepancies.csv'}")

    # Write missing_in_sqlite.csv
    missing_sqlite = df[df['missing_in_sqlite'] == True]
    if len(missing_sqlite) > 0:
        missing_sqlite.to_csv(outdir / 'missing_in_sqlite.csv', index=False)
        logger.info(f"Wrote {outdir / 'missing_in_sqlite.csv'} ({len(missing_sqlite)} rows)")

    # Write missing_in_excel.csv
    missing_excel = df[df['missing_in_excel'] == True]
    if len(missing_excel) > 0:
        missing_excel.to_csv(outdir / 'missing_in_excel.csv', index=False)
        logger.info(f"Wrote {outdir / 'missing_in_excel.csv'} ({len(missing_excel)} rows)")

    # Write summary.json
    with open(outdir / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Wrote {outdir / 'summary.json'}")

    # Write markdown report
    write_markdown_report(results, summary, outdir / 'report.md')
    logger.info(f"Wrote {outdir / 'report.md'}")

    return summary


def write_markdown_report(
    results: List[ComparisonRow],
    summary: Dict[str, Any],
    output_path: Path
) -> None:
    """Write a human-readable markdown report."""

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
        f"| Metric | Count |",
        f"|--------|-------|",
        f"| Exact matches | {summary['matched_exact']} |",
        f"| Fuzzy matches | {summary['matched_fuzzy']} |",
        f"| Unmatched | {summary['unmatched']} |",
        f"| Missing in SQLite | {summary['missing_in_sqlite']} |",
        f"| Missing in Excel | {summary['missing_in_excel']} |",
        f"| Count mismatches | {summary['count_mismatches']} |",
        "",
    ]

    # Top mismatches table
    if summary['count_mismatches'] > 0:
        lines.extend([
            "## Top Count Mismatches",
            "",
            "| PMID | Excel Variant | SQLite Variant | Carriers (E/S/Δ) | Affected (E/S/Δ) |",
            "|------|---------------|----------------|------------------|------------------|",
        ])

        mismatched = [r for r in results if r.count_mismatch]
        mismatched.sort(key=lambda r: abs(r.carriers_diff or 0) + abs(r.affected_diff or 0), reverse=True)

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
        lines.extend([
            "## Variants Missing in SQLite (Top 50)",
            "",
            "| PMID | Excel Variant |",
            "|------|---------------|",
        ])
        for r in missing_sqlite[:50]:
            lines.append(f"| {r.pmid} | {r.excel_variant_raw} |")
        lines.append("")

    # Missing in Excel
    missing_excel = [r for r in results if r.missing_in_excel]
    if missing_excel:
        lines.extend([
            "## Variants Missing in Excel (Top 50)",
            "",
            "| PMID | SQLite Variant |",
            "|------|----------------|",
        ])
        for r in missing_excel[:50]:
            lines.append(f"| {r.pmid} | {r.sqlite_variant_raw} |")
        lines.append("")

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))


# =============================================================================
# CLI ENTRYPOINT
# =============================================================================

def load_mapping_file(mapping_path: Path) -> Dict[str, str]:
    """Load column mapping from JSON or YAML file."""
    with open(mapping_path, 'r') as f:
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
        description="Compare variant data between Excel curation and SQLite database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic comparison
  python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db

  # Fuzzy matching with lower threshold
  python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db \\
      --variant_match_mode fuzzy --fuzzy_threshold 0.80

  # Specify sheet and output directory
  python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db \\
      --sheet "Variants" --outdir ./results

  # Use explicit column mapping
  python compare_variants.py --excel curated.xlsx --sqlite KCNH2.db \\
      --mapping column_map.yaml
"""
    )

    parser.add_argument(
        '--excel', '-e',
        type=Path,
        required=True,
        help="Path to Excel file with curated variant data"
    )

    parser.add_argument(
        '--sqlite', '-s',
        type=Path,
        required=True,
        help="Path to SQLite database (e.g., KCNH2.db)"
    )

    parser.add_argument(
        '--sheet',
        type=str,
        default=None,
        help="Excel sheet name (default: first sheet)"
    )

    parser.add_argument(
        '--outdir', '-o',
        type=Path,
        default=Path('./compare_out'),
        help="Output directory (default: ./compare_out)"
    )

    parser.add_argument(
        '--mapping', '-m',
        type=Path,
        default=None,
        help="JSON or YAML file with explicit column mappings"
    )

    parser.add_argument(
        '--variant_match_mode',
        choices=['exact', 'fuzzy'],
        default='exact',
        help="Variant matching mode (default: exact)"
    )

    parser.add_argument(
        '--fuzzy_threshold',
        type=float,
        default=0.85,
        help="Minimum similarity for fuzzy matching (default: 0.85)"
    )

    parser.add_argument(
        '--log_level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help="Logging level (default: INFO)"
    )

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
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
    if args.variant_match_mode == 'fuzzy':
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
            args.fuzzy_threshold
        )

        # Generate outputs
        summary = generate_outputs(
            results, args.outdir, args.excel, args.sqlite
        )

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
