#!/usr/bin/env python3
"""
SQLite Migration Script for Gene Variant Fetcher

Migrates file-based extraction data to a normalized SQLite database.
Includes cleanup and archival functions for the file system.

The script automatically detects the gene symbol from the directory path
(e.g., /output/TTR/...) or extraction JSON files and creates a
gene-specific database (e.g., TTR.db).

Usage:
    # Point to parent directory containing extractions/ subdirectory
    # This will create TTR.db automatically
    python migrate_to_sqlite.py --data-dir /path/to/output/TTR/20251125_114028

    # Point directly to directory containing JSON files
    # Database name will be auto-detected (e.g., TTR.db)
    python migrate_to_sqlite.py --data-dir /path/to/output/TTR/20251125_114028/extractions
    python migrate_to_sqlite.py --data-dir /path/to/output/TTR/20251125_114028/extractions_rerun/20251125_151454

    # Specify custom database path (overrides auto-detection)
    python migrate_to_sqlite.py --data-dir /path/to/your/data --db custom_name.db

The script will automatically find JSON extraction files in:
1. The specified directory itself (if it contains *_PMID_*.json files)
2. An 'extractions' subdirectory
3. Alternative subdirectories (extractions_rerun, extraction, etc.)
4. Nested timestamped subdirectories
"""

import argparse
import hashlib
import json
import logging
import os
import re
import shutil
import sqlite3
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

# Configure logging using centralized utility
from utils.logging_utils import get_logger, setup_logging
from utils.pmid_utils import extract_gene_from_filename, extract_pmid_from_filename
from utils.geo_ancestry import enrich_ethnicity_origin
from utils.source_layers import (
    infer_source_layer_from_text,
    junk_notation_reason,
)

setup_logging(level=logging.INFO)
logger = get_logger(__name__)

PROTEIN_NOTATION_RE = re.compile(
    r"^(?:p\.)?(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])"
    r"\d{1,4}"
    # Optional second residue for HGVS range notations like p.Asp2_Arg135del
    # or p.Lys100_Glu105delinsX. Without this, valid multi-residue dels/ins
    # are silently dropped at migration. delins must precede del.
    r"(?:[_-](?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])\d{1,4})?"
    r"(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY*X?]|fs(?:X|\*)?\d*"
    r"|delins[ACDEFGHIKLMNPQRSTVWY]+"
    r"|del[ACDEFGHIKLMNPQRSTVWY]*|dup|ins[ACDEFGHIKLMNPQRSTVWY]*)",
    re.IGNORECASE,
)
CDNA_NOTATION_RE = re.compile(
    r"^c\.\d+(?:[+-]\d+)?[ACGT]>[ACGT]$"
    r"|^c\.\d+(?:[+-]\d+)?(?:del|dup|ins)[ACGT]*$"
    r"|^c\.\d+(?:_\d+)?(?:del|dup|ins)[ACGT]*$"
    r"|^c\.\d+(?:_\d+)?del[ACGT]*ins[ACGT]+$"
    r"|^c\.\d+(?:[+-]\d+)?del[ACGT]*ins[ACGT]+$"
    r"|^IVS\d+[+-]\d+[ACGT]>[ACGT]$"
    r"|^IVS\d+[+-]\d+(?:del|dup|ins)[ACGT]*$",
    re.IGNORECASE,
)

# Closed vocabulary for per-variant structural classification (Stage 5 B1).
VARIANT_CLASS_VALUES = frozenset(
    {
        "missense",
        "nonsense",
        "frameshift",
        "inframe_indel",
        "splice",
        "deep_intronic",
        "large_deletion",
        "large_duplication",
        "cnv",
        "exon_deletion",
        "exon_duplication",
        "complex",
        "other",
    }
)

FACT_PROVENANCE_FIELDS: Tuple[str, ...] = (
    "variant_id",
    "pmid",
    "fact_type",
    "fact_value",
    "individual_id",
    "source_location",
    "source_section",
    "source_paragraph",
    "source_table",
    "source_row",
    "source_column",
    "evidence_quote",
    "count_type",
    "source_layer",
    "provenance_kind",
)

OBSERVATION_PROVENANCE_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("source_container", "TEXT"),
    ("source_kind", "TEXT"),
    ("source_ref", "TEXT"),
    ("page_label", "TEXT"),
    ("pdf_page", "INTEGER"),
    ("row_label", "TEXT"),
    ("row_ordinal", "INTEGER"),
    ("column_ref", "TEXT"),
    ("figure_panel", "TEXT"),
    ("source_record_id", "TEXT"),
    ("locator_extra", "TEXT"),
)

OBSERVATION_PROVENANCE_KEYS: Tuple[str, ...] = tuple(
    col for col, _decl in OBSERVATION_PROVENANCE_COLUMNS
)

# =============================================================================
# INPUT VALIDATION
# =============================================================================


class ValidationError(Exception):
    """Raised when input validation fails."""

    pass


def sanitize_variant_notation(variant_data: Dict[str, Any]) -> bool:
    """
    Drop malformed protein/cDNA notation before SQLite insertion.

    Returns True if the variant still has at least one usable notation. This is
    a migration-time backstop for older extraction JSON files that may contain
    table artifacts such as single alleles (A/C/G/T), p-values, or cohort sizes.

    Structural events with a valid ``variant_class`` + ``structural_description``
    are kept even without point-form protein/cDNA notation.
    """
    protein = (variant_data.get("protein_notation") or "").strip().replace(" ", "")
    cdna = (variant_data.get("cdna_notation") or "").strip().replace(" ", "")
    genomic = (variant_data.get("genomic_position") or "").strip()
    vclass = (variant_data.get("variant_class") or "").strip().lower()
    structural = (variant_data.get("structural_description") or "").strip()

    if protein and not PROTEIN_NOTATION_RE.match(protein):
        variant_data["protein_notation"] = None
        protein = ""
    if cdna and not CDNA_NOTATION_RE.match(cdna):
        variant_data["cdna_notation"] = None
        cdna = ""

    if vclass and vclass not in VARIANT_CLASS_VALUES:
        variant_data["variant_class"] = None
        vclass = ""

    has_structural = bool(vclass and structural)
    return bool(protein or cdna or genomic or has_structural)


def infer_source_layer(variant_data: Dict[str, Any]) -> str:
    """Return a stable source-layer label for a variant-paper link."""

    return infer_source_layer_from_text(
        source_location=variant_data.get("source_location"),
        additional_notes=variant_data.get("additional_notes"),
        extraction_source=variant_data.get("extraction_source"),
        source_layer=variant_data.get("source_layer"),
    )


def normalize_affected_status(status: Any) -> Optional[str]:
    """Normalize extracted individual status values to the SQLite enum."""

    if status is None:
        return None
    text = str(status).strip().lower()
    if not text:
        return None
    if text in {"affected", "case", "symptomatic", "diseased", "patient"}:
        return "affected"
    if text in {
        "unaffected",
        "control",
        "asymptomatic",
        "healthy",
        "normal",
        "nonaffected",
        "non-affected",
    }:
        return "unaffected"
    if text in {
        "uncertain",
        "unknown",
        "ambiguous",
        "carrier",
        "carrier only",
        "not reported",
        "na",
        "n/a",
    }:
        return "uncertain"
    return "uncertain"


def validate_input_directory(data_dir: Path) -> None:
    """
    Validate that input directory exists and is readable.

    Args:
        data_dir: Path to data directory

    Raises:
        ValidationError: If directory is invalid
    """
    if not data_dir.exists():
        raise ValidationError(f"Input directory does not exist: {data_dir}")

    if not data_dir.is_dir():
        raise ValidationError(f"Input path is not a directory: {data_dir}")

    if not os.access(data_dir, os.R_OK):
        raise ValidationError(f"Input directory is not readable: {data_dir}")


def validate_has_extraction_files(extraction_dir: Path) -> int:
    """
    Validate that directory contains extraction JSON files.

    Args:
        extraction_dir: Path to directory with extraction files

    Returns:
        Number of JSON files found

    Raises:
        ValidationError: If no extraction files found
    """
    json_files = find_extraction_json_files(extraction_dir)

    if not json_files:
        raise ValidationError(
            f"No extraction JSON files found in {extraction_dir}. "
            f"Expected files matching pattern: *_PMID_*.json or *_extraction.json"
        )

    return len(json_files)


_EXTRACTION_JSON_NAME_RE = re.compile(
    r"^(?:[A-Za-z0-9_-]+_PMID_\d+|[A-Za-z0-9_-]+_PMID_FULL|\d+_extraction)\.json$"
)


def is_extraction_json_file(path: Path) -> bool:
    """Return True for canonical extraction JSON filenames only."""
    return bool(_EXTRACTION_JSON_NAME_RE.match(path.name))


def find_extraction_json_files(extraction_dir: Path) -> List[Path]:
    """Return extraction JSON files in all supported GVF naming schemes.

    Keep this intentionally stricter than ``*_PMID_*.json`` so timestamped
    backups such as ``KCNH2_PMID_123.20260514_pre.json`` are not migrated as
    duplicate papers.
    """
    return sorted(
        (p for p in extraction_dir.glob("*.json") if is_extraction_json_file(p)),
        key=lambda p: p.name,
    )


def validate_db_path_writable(db_path: str) -> None:
    """
    Validate that database path is writable.

    Args:
        db_path: Path to SQLite database file

    Raises:
        ValidationError: If database path is not writable
    """
    db_path = Path(db_path)

    if db_path.exists():
        # Check if existing file is writable
        if not os.access(db_path, os.W_OK):
            raise ValidationError(
                f"Database file exists but is not writable: {db_path}"
            )
    else:
        # Check if parent directory is writable
        parent = db_path.parent if db_path.parent.exists() else Path.cwd()
        if not os.access(parent, os.W_OK):
            raise ValidationError(
                f"Cannot create database file (directory not writable): {db_path}"
            )


def validate_migrate_inputs(
    data_dir: Path,
    extraction_dir: Optional[Path],
    db_path: str,
) -> Path:
    """
    Validate all inputs for the migrate operation.

    Args:
        data_dir: Path to data directory
        extraction_dir: Path to directory with extraction files (may be None if not yet found)
        db_path: Path to SQLite database file

    Returns:
        Validated extraction directory path

    Raises:
        ValidationError: If validation fails
    """
    # Validate input directory
    validate_input_directory(data_dir)

    # Validate database path is writable
    validate_db_path_writable(db_path)

    # If extraction_dir is provided, validate it has files
    if extraction_dir:
        validate_input_directory(extraction_dir)
        validate_has_extraction_files(extraction_dir)
        return extraction_dir

    return data_dir


def validate_extraction_data(
    data: Dict[str, Any], filename: str = ""
) -> Tuple[bool, List[str], List[str]]:
    """
    Validate extraction data before migration.

    Args:
        data: Extraction JSON data
        filename: Source filename (for error messages)

    Returns:
        Tuple of (is_valid, errors, warnings)
    """
    errors = []
    warnings = []

    # Check paper_metadata
    paper_meta = data.get("paper_metadata", {})
    if not paper_meta:
        errors.append("Missing paper_metadata section")
    else:
        pmid = paper_meta.get("pmid")
        if not pmid or pmid == "UNKNOWN":
            errors.append(f"Missing or invalid paper_metadata.pmid: {pmid}")
        if not paper_meta.get("title"):
            warnings.append("Missing paper_metadata.title")

    # Check variants
    variants = data.get("variants", [])
    if not variants:
        warnings.append("No variants in extraction")
    else:
        for i, variant in enumerate(variants):
            gene = variant.get("gene_symbol")
            if not gene:
                errors.append(f"Variant {i}: missing required field 'gene_symbol'")

            # Check for at least one notation
            has_notation = any(
                [
                    variant.get("cdna_notation"),
                    variant.get("protein_notation"),
                    variant.get("genomic_position"),
                ]
            )
            if not has_notation:
                warnings.append(f"Variant {i}: no notation (cdna/protein/genomic)")

            # Check individuals for penetrance data
            individuals = variant.get("individuals", [])
            for j, ind in enumerate(individuals):
                if not ind.get("affected_status"):
                    warnings.append(
                        f"Variant {i}, Individual {j}: missing affected_status"
                    )

    is_valid = len(errors) == 0
    return is_valid, errors, warnings


def repair_extraction_data(
    data: Dict[str, Any], filename: str, gene_symbol_override: Optional[str] = None
) -> Tuple[Dict[str, Any], List[str]]:
    """
    Attempt to repair common issues in extraction data.

    Args:
        data: Extraction JSON data
        filename: Source filename (for PMID/gene extraction)
        gene_symbol_override: Force this gene symbol if provided

    Returns:
        Tuple of (repaired_data, list of repairs made)
    """
    repairs = []
    data = json.loads(json.dumps(data))  # Deep copy

    # Extract info from filename
    filename_pmid = extract_pmid_from_filename(filename)
    filename_gene = extract_gene_from_filename(filename)
    default_gene = gene_symbol_override or filename_gene

    # Ensure paper_metadata exists
    if "paper_metadata" not in data:
        data["paper_metadata"] = {
            k: data[k]
            for k in ("pmid", "title", "extraction_summary")
            if data.get(k) is not None
        }
        repairs.append("Created missing paper_metadata section")

    paper_meta = data["paper_metadata"]

    # Fix missing PMID
    if not paper_meta.get("pmid") or paper_meta.get("pmid") == "UNKNOWN":
        if filename_pmid:
            old_val = paper_meta.get("pmid", "MISSING")
            paper_meta["pmid"] = filename_pmid
            repairs.append(f"Set pmid from filename: {old_val} -> {filename_pmid}")

    # Fix missing title
    if not paper_meta.get("title"):
        paper_meta["title"] = f"Paper {paper_meta.get('pmid', 'Unknown')}"
        repairs.append("Set default title")

    # Store gene in paper_metadata for fallback
    if default_gene and not paper_meta.get("gene_symbol"):
        paper_meta["gene_symbol"] = default_gene
        repairs.append(f"Set paper_metadata.gene_symbol: {default_gene}")

    # Repair variants
    variants = data.get("variants", [])
    for i, variant in enumerate(variants):
        if not variant.get("gene_symbol"):
            gene = paper_meta.get("gene_symbol") or default_gene or "UNKNOWN_GENE"
            variant["gene_symbol"] = gene
            repairs.append(f"Variant {i}: set gene_symbol to '{gene}'")

        # Ensure individuals have affected_status
        repaired_individuals = []
        for j, ind in enumerate(variant.get("individuals", [])):
            if not isinstance(ind, dict):
                repairs.append(f"Variant {i}, Individual {j}: dropped non-object")
                continue
            if not ind.get("affected_status"):
                ind["affected_status"] = "uncertain"
                repairs.append(f"Variant {i}, Individual {j}: set affected_status")
            else:
                normalized = normalize_affected_status(ind.get("affected_status"))
                if normalized != ind.get("affected_status"):
                    repairs.append(
                        f"Variant {i}, Individual {j}: normalized affected_status"
                    )
                ind["affected_status"] = normalized
            repaired_individuals.append(ind)
        if len(repaired_individuals) != len(variant.get("individuals", [])):
            variant["individuals"] = repaired_individuals

        repaired_records = []
        for j, record in enumerate(variant.get("individual_records", [])):
            if not isinstance(record, dict):
                repairs.append(f"Variant {i}, Record {j}: dropped non-object")
                continue
            normalized = normalize_affected_status(record.get("affected_status"))
            if normalized != record.get("affected_status"):
                repairs.append(f"Variant {i}, Record {j}: normalized affected_status")
            record["affected_status"] = normalized or "uncertain"
            repaired_records.append(record)
        if len(repaired_records) != len(variant.get("individual_records", [])):
            variant["individual_records"] = repaired_records

    # Mark as repaired
    if repairs:
        if "extraction_metadata" not in data:
            data["extraction_metadata"] = {}
        data["extraction_metadata"]["was_repaired"] = True
        data["extraction_metadata"]["repair_count"] = len(repairs)

    return data, repairs


# ============================================================================
# DATABASE SCHEMA INITIALIZATION
# ============================================================================


def create_database_schema(db_path: str) -> sqlite3.Connection:
    """
    Create the SQLite database with a normalized schema.

    Args:
        db_path: Path to the SQLite database file

    Returns:
        Database connection
    """
    logger.info(f"Creating database schema at: {db_path}")

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Enable foreign keys
    cursor.execute("PRAGMA foreign_keys = ON")

    # ========================================================================
    # Papers table: Stores paper metadata
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS papers (
            pmid TEXT PRIMARY KEY,
            title TEXT,
            first_author TEXT,
            journal TEXT,
            publication_date TEXT,
            doi TEXT,
            pmc_id TEXT,
            gene_symbol TEXT,
            extraction_summary TEXT,
            extraction_timestamp TEXT,
            created_at TEXT DEFAULT CURRENT_TIMESTAMP
        )
    """)

    # ========================================================================
    # Variants table: Stores unique variants with normalized identifiers
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT NOT NULL,
            cdna_notation TEXT,
            protein_notation TEXT,
            genomic_position TEXT,
            clinical_significance TEXT,
            evidence_level TEXT,
            variant_class TEXT,  -- missense, splice, exon_deletion, cnv, ...
            structural_description TEXT,  -- free text e.g. "deletion of exons 3-5"

            -- Composite unique constraint (flexible nulls)
            UNIQUE(gene_symbol, cdna_notation, protein_notation, genomic_position)
        )
    """)

    # Create index for faster lookups
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_variants_gene
        ON variants(gene_symbol)
    """)

    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_variants_protein
        ON variants(protein_notation)
    """)

    # ========================================================================
    # Variant-Paper association: Many-to-many relationship
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS variant_papers (
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            source_location TEXT,
            additional_notes TEXT,
            key_quotes TEXT,  -- JSON array of quotes
            count_provenance TEXT,  -- JSON: which column/count-type each carrier/affected count came from (the "why")
            source_layer TEXT,  -- stable provenance layer: llm_table, clinvar, pubtator, figure, etc.

            PRIMARY KEY (variant_id, pmid),
            FOREIGN KEY (variant_id) REFERENCES variants(variant_id) ON DELETE CASCADE,
            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    # ========================================================================
    # Penetrance data: Cohort-level statistics per variant per paper
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS penetrance_data (
            penetrance_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            total_carriers_observed INTEGER,
            affected_count INTEGER,
            unaffected_count INTEGER,
            uncertain_count INTEGER,
            penetrance_percentage REAL,
            trust_tier TEXT DEFAULT 'trusted',
            trust_reasons TEXT,
            trust_rule_version TEXT,

            FOREIGN KEY (variant_id) REFERENCES variants(variant_id) ON DELETE CASCADE,
            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    # ========================================================================
    # Age-dependent penetrance: Age-stratified penetrance data
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS age_dependent_penetrance (
            age_penetrance_id INTEGER PRIMARY KEY AUTOINCREMENT,
            penetrance_id INTEGER NOT NULL,
            age_range TEXT NOT NULL,
            penetrance_percentage REAL,
            carriers_in_range INTEGER,
            affected_in_range INTEGER,

            FOREIGN KEY (penetrance_id) REFERENCES penetrance_data(penetrance_id) ON DELETE CASCADE
        )
    """)

    # ========================================================================
    # Individual records: Person-level carrier and affected status
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS individual_records (
            record_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            individual_id TEXT,
            age_at_evaluation INTEGER,
            age_at_onset INTEGER,
            age_at_diagnosis INTEGER,
            sex TEXT,
            affected_status TEXT CHECK(affected_status IN ('affected', 'unaffected', 'uncertain')),
            phenotype_details TEXT,
            evidence_sentence TEXT,
            ethnicity TEXT,
            geographic_origin TEXT,
            source_container TEXT,
            source_kind TEXT,
            source_ref TEXT,
            page_label TEXT,
            pdf_page INTEGER,
            row_label TEXT,
            row_ordinal INTEGER,
            column_ref TEXT,
            figure_panel TEXT,
            source_record_id TEXT,
            locator_extra TEXT,

            FOREIGN KEY (variant_id) REFERENCES variants(variant_id) ON DELETE CASCADE,
            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    # ========================================================================
    # Fact provenance: exact source pointers for variant/count/person facts
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS fact_provenance (
            provenance_id INTEGER PRIMARY KEY AUTOINCREMENT,
            fact_hash TEXT NOT NULL UNIQUE,
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            fact_type TEXT NOT NULL,
            fact_value TEXT,
            individual_id TEXT,
            source_location TEXT,
            source_section TEXT,
            source_paragraph TEXT,
            source_table TEXT,
            source_row TEXT,
            source_column TEXT,
            evidence_quote TEXT,
            count_type TEXT,
            source_layer TEXT,
            provenance_kind TEXT,
            created_at TEXT DEFAULT CURRENT_TIMESTAMP,

            FOREIGN KEY (variant_id) REFERENCES variants(variant_id) ON DELETE CASCADE,
            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_fact_provenance_variant_paper
        ON fact_provenance(variant_id, pmid, fact_type)
    """)

    # Backfill new columns onto pre-existing DBs (CREATE IF NOT EXISTS won't add them).
    for table, col, decl in (
        ("papers", "first_author", "TEXT"),
        ("variant_papers", "source_location", "TEXT"),
        ("variant_papers", "additional_notes", "TEXT"),
        ("variant_papers", "key_quotes", "TEXT"),
        ("variant_papers", "count_provenance", "TEXT"),
        ("variant_papers", "source_layer", "TEXT"),
        ("penetrance_data", "trust_tier", "TEXT DEFAULT 'trusted'"),
        ("penetrance_data", "trust_reasons", "TEXT"),
        ("penetrance_data", "trust_rule_version", "TEXT"),
        ("individual_records", "ethnicity", "TEXT"),
        ("individual_records", "geographic_origin", "TEXT"),
        ("extraction_metadata", "study_type", "TEXT"),
        ("extraction_metadata", "study_design", "TEXT"),
        ("extraction_metadata", "ascertainment", "TEXT"),
        ("extraction_metadata", "cohort_source", "TEXT"),
        ("extraction_metadata", "population", "TEXT"),
        ("extraction_metadata", "study_summary", "TEXT"),
        ("variants", "variant_class", "TEXT"),
        ("variants", "structural_description", "TEXT"),
        *(
            ("individual_records", col, decl)
            for col, decl in OBSERVATION_PROVENANCE_COLUMNS
        ),
        *(("phenotypes", col, decl) for col, decl in OBSERVATION_PROVENANCE_COLUMNS),
    ):
        table_exists = cursor.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
            (table,),
        ).fetchone()
        if not table_exists:
            continue
        existing = {r[1] for r in cursor.execute(f"PRAGMA table_info({table})")}
        if col not in existing:
            cursor.execute(f"ALTER TABLE {table} ADD COLUMN {col} {decl}")

    # Index for querying individual records by affected status
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_individual_records_affected
        ON individual_records(variant_id, affected_status)
    """)

    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_individual_records_source_record
        ON individual_records(source_record_id)
    """)

    # ========================================================================
    # Functional data: In vitro and functional assay data
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS functional_data (
            functional_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            summary TEXT,
            assays TEXT,  -- JSON array of assays

            FOREIGN KEY (variant_id) REFERENCES variants(variant_id) ON DELETE CASCADE,
            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    # ========================================================================
    # Phenotypes: Patient phenotype descriptions
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS phenotypes (
            phenotype_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            patient_count INTEGER,
            demographics TEXT,
            phenotype_description TEXT,
            source_container TEXT,
            source_kind TEXT,
            source_ref TEXT,
            page_label TEXT,
            pdf_page INTEGER,
            row_label TEXT,
            row_ordinal INTEGER,
            column_ref TEXT,
            figure_panel TEXT,
            source_record_id TEXT,
            locator_extra TEXT,

            FOREIGN KEY (variant_id) REFERENCES variants(variant_id) ON DELETE CASCADE,
            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_phenotypes_source_record
        ON phenotypes(source_record_id)
    """)

    # ========================================================================
    # Additional variant metadata
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS variant_metadata (
            metadata_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            pmid TEXT NOT NULL,
            segregation_data TEXT,
            population_frequency TEXT,

            FOREIGN KEY (variant_id) REFERENCES variants(variant_id) ON DELETE CASCADE,
            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    # ========================================================================
    # Extraction metadata: Track extraction quality and processing
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS extraction_metadata (
            extraction_id INTEGER PRIMARY KEY AUTOINCREMENT,
            pmid TEXT NOT NULL,
            total_variants_found INTEGER,
            extraction_confidence TEXT,
            study_type TEXT,  -- clinical/functional/mixed (legacy coarse label)
            study_design TEXT,  -- case_report, cohort_population, gwas, ...
            ascertainment TEXT,  -- proband_referral, biobank, family_cascade, ...
            cohort_source TEXT,  -- free-text cohort origin
            population TEXT,  -- ancestry/geography/founder population
            study_summary TEXT,  -- 1-3 sentence study narrative
            challenges TEXT,  -- JSON array
            notes TEXT,
            model_used TEXT,
            extraction_timestamp TEXT,
            source_type TEXT,  -- 'fulltext', 'abstract_only', or NULL
            abstract_only INTEGER DEFAULT 0,  -- 1 if extracted from abstract only
            source_file TEXT,  -- Path to source file (DATA_ZONES.md, FULL_CONTEXT.md, etc.)

            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    # ========================================================================
    # Tables processed: Track which tables were processed in each paper
    # ========================================================================
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS tables_processed (
            table_id INTEGER PRIMARY KEY AUTOINCREMENT,
            pmid TEXT NOT NULL,
            table_name TEXT,
            table_caption TEXT,
            variants_extracted INTEGER,

            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    conn.commit()
    logger.info("✓ Database schema created successfully")

    # Upgrade existing schema if needed
    upgrade_database_schema(conn)

    return conn


def upgrade_database_schema(conn: sqlite3.Connection) -> None:
    """
    Upgrade existing database schema to add new columns.

    Handles backward compatibility for databases created before
    new columns were added.

    Args:
        conn: Database connection
    """
    cursor = conn.cursor()

    def ensure_column(table: str, column: str, declaration: str) -> None:
        table_exists = cursor.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
            (table,),
        ).fetchone()
        if not table_exists:
            return
        cursor.execute(f"PRAGMA table_info({table})")
        columns = {row[1] for row in cursor.fetchall()}
        if column in columns:
            return
        logger.info("Upgrading schema: adding %s column to %s", column, table)
        cursor.execute(f"ALTER TABLE {table} ADD COLUMN {column} {declaration}")

    # Check if source_file column exists in extraction_metadata
    cursor.execute("PRAGMA table_info(extraction_metadata)")
    columns = {row[1] for row in cursor.fetchall()}

    if "source_file" not in columns:
        logger.info(
            "Upgrading schema: adding source_file column to extraction_metadata"
        )
        cursor.execute("""
            ALTER TABLE extraction_metadata
            ADD COLUMN source_file TEXT
        """)
        conn.commit()
        logger.info("✓ Schema upgrade complete: source_file column added")

    for table in ("individual_records", "phenotypes"):
        for column, declaration in OBSERVATION_PROVENANCE_COLUMNS:
            ensure_column(table, column, declaration)
    conn.commit()


# ============================================================================
# ETL FUNCTIONS: Extract, Transform, Load
# ============================================================================


def get_or_create_variant(cursor: sqlite3.Cursor, variant_data: Dict[str, Any]) -> int:
    """
    Get existing variant ID or create new variant entry.

    Args:
        cursor: Database cursor
        variant_data: Variant data dictionary

    Returns:
        variant_id
    """
    gene_symbol = variant_data.get("gene_symbol")
    cdna = variant_data.get("cdna_notation")
    protein = variant_data.get("protein_notation")
    genomic = variant_data.get("genomic_position")
    vclass = (variant_data.get("variant_class") or "").strip().lower() or None
    if vclass and vclass not in VARIANT_CLASS_VALUES:
        vclass = None
    structural = (variant_data.get("structural_description") or "").strip() or None

    # Point-form identity first; structural-only events key on description.
    if cdna or protein or genomic:
        cursor.execute(
            """
            SELECT variant_id FROM variants
            WHERE gene_symbol = ?
            AND (cdna_notation = ? OR (cdna_notation IS NULL AND ? IS NULL))
            AND (protein_notation = ? OR (protein_notation IS NULL AND ? IS NULL))
            AND (genomic_position = ? OR (genomic_position IS NULL AND ? IS NULL))
        """,
            (gene_symbol, cdna, cdna, protein, protein, genomic, genomic),
        )
    elif structural:
        cursor.execute(
            """
            SELECT variant_id FROM variants
            WHERE gene_symbol = ?
            AND structural_description = ?
            AND cdna_notation IS NULL
            AND protein_notation IS NULL
            AND genomic_position IS NULL
        """,
            (gene_symbol, structural),
        )
    else:
        cursor.execute("SELECT variant_id FROM variants WHERE 0")

    result = cursor.fetchone()
    if result:
        vid = result[0]
        if vclass or structural:
            cursor.execute(
                """
                UPDATE variants
                SET variant_class = COALESCE(variant_class, ?),
                    structural_description = COALESCE(structural_description, ?)
                WHERE variant_id = ?
                """,
                (vclass, structural, vid),
            )
        return vid

    cursor.execute(
        """
        INSERT INTO variants (
            gene_symbol, cdna_notation, protein_notation,
            genomic_position, clinical_significance, evidence_level,
            variant_class, structural_description
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """,
        (
            gene_symbol,
            cdna,
            protein,
            genomic,
            variant_data.get("clinical_significance"),
            variant_data.get("evidence_level"),
            vclass,
            structural,
        ),
    )

    return cursor.lastrowid


def insert_paper_metadata(
    cursor: sqlite3.Cursor,
    extraction_data: Dict[str, Any],
    replace_existing: bool = True,
) -> None:
    """
    Insert or update paper metadata.

    Args:
        cursor: Database cursor
        extraction_data: Extraction JSON data
        replace_existing: If True, use the legacy replace behavior. If False,
            upsert metadata without deleting existing PMID-linked evidence.
    """
    paper_meta = extraction_data.get("paper_metadata", {})
    pmid = paper_meta.get("pmid")

    if not pmid:
        logger.warning("No PMID found in extraction data")
        return

    gene_symbol = (
        extraction_data.get("variants", [{}])[0].get("gene_symbol")
        if extraction_data.get("variants")
        else None
    )
    values = (
        pmid,
        paper_meta.get("title"),
        paper_meta.get("first_author"),
        paper_meta.get("journal"),
        paper_meta.get("publication_date"),
        paper_meta.get("doi"),
        paper_meta.get("pmc_id"),
        gene_symbol,
        paper_meta.get("extraction_summary"),
        datetime.now().isoformat(),
    )
    cols = (
        "pmid, title, first_author, journal, publication_date, doi, pmc_id, "
        "gene_symbol, extraction_summary, extraction_timestamp"
    )
    placeholders = ", ".join(["?"] * len(values))

    if replace_existing:
        cursor.execute(
            f"INSERT OR REPLACE INTO papers ({cols}) VALUES ({placeholders})", values
        )
        return

    cursor.execute(
        f"""
        INSERT INTO papers ({cols}) VALUES ({placeholders})
        ON CONFLICT(pmid) DO UPDATE SET
            title = COALESCE(excluded.title, papers.title),
            first_author = COALESCE(excluded.first_author, papers.first_author),
            journal = COALESCE(excluded.journal, papers.journal),
            publication_date = COALESCE(excluded.publication_date, papers.publication_date),
            doi = COALESCE(excluded.doi, papers.doi),
            pmc_id = COALESCE(excluded.pmc_id, papers.pmc_id),
            gene_symbol = COALESCE(papers.gene_symbol, excluded.gene_symbol),
            extraction_summary = COALESCE(
                excluded.extraction_summary, papers.extraction_summary
            ),
            extraction_timestamp = excluded.extraction_timestamp
    """,
        values,
    )


def _matching_row_id(
    cursor: sqlite3.Cursor,
    table: str,
    id_column: str,
    match: Dict[str, Any],
) -> Optional[int]:
    """Return the first matching synthetic id, or None.

    NULL-aware (uses ``col IS NULL`` rather than ``col = NULL``). This keeps the
    child-row inserts idempotent: the same fact extracted from several table
    representations of one variant (e.g. multi-cohort / multi-classification
    supplement tables) — or re-migrated on a refresh — must not create a second
    identical row, which would make carrier/affected counts sum N-fold at scoring
    time. Table and column names are code-controlled (never user input).
    """
    clauses: List[str] = []
    params: List[Any] = []
    for col, val in match.items():
        if val is None:
            clauses.append(f"{col} IS NULL")
        else:
            clauses.append(f"{col} = ?")
            params.append(val)
    sql = (  # noqa: S608
        f"SELECT {id_column} FROM {table} "
        f"WHERE {' AND '.join(clauses)} ORDER BY {id_column} LIMIT 1"
    )
    row = cursor.execute(sql, params).fetchone()
    return int(row[0]) if row else None


def _row_exists(cursor: sqlite3.Cursor, table: str, match: Dict[str, Any]) -> bool:
    """Return True if *table* already holds a row matching every column in *match*."""
    return _matching_row_id(cursor, table, "rowid", match) is not None


def _insert_age_dependent_penetrance(
    cursor: sqlite3.Cursor,
    penetrance_id: int,
    age_dependent_rows: Iterable[Dict[str, Any]],
) -> None:
    """Insert unique age-bin rows under an existing or newly-created parent."""
    for age_dep in age_dependent_rows:
        if not isinstance(age_dep, dict):
            logger.warning(
                "Skipping non-object age_dependent_penetrance for penetrance_id %s",
                penetrance_id,
            )
            continue
        row = {
            "penetrance_id": penetrance_id,
            "age_range": age_dep.get("age_range"),
            "penetrance_percentage": age_dep.get("penetrance_percentage"),
            "carriers_in_range": age_dep.get("carriers_in_range"),
            "affected_in_range": age_dep.get("affected_in_range"),
        }
        if _row_exists(cursor, "age_dependent_penetrance", row):
            continue
        cursor.execute(
            """
            INSERT INTO age_dependent_penetrance (
                penetrance_id, age_range, penetrance_percentage,
                carriers_in_range, affected_in_range
            ) VALUES (?, ?, ?, ?, ?)
        """,
            (
                penetrance_id,
                row["age_range"],
                row["penetrance_percentage"],
                row["carriers_in_range"],
                row["affected_in_range"],
            ),
        )


def _clean_optional_text(value: Any) -> Optional[str]:
    """Return a stripped text value, preserving None for missing/blank fields."""
    if value is None:
        return None
    if isinstance(value, (dict, list)):
        text = json.dumps(value, sort_keys=True, ensure_ascii=False)
    else:
        text = str(value)
    text = text.strip()
    return text or None


def _first_quote(variant_data: Dict[str, Any]) -> Optional[str]:
    quotes = variant_data.get("key_quotes")
    if isinstance(quotes, list):
        for quote in quotes:
            text = _clean_optional_text(quote)
            if text:
                return text
        # An empty / quote-less list must yield None, not the JSON-dumped "[]"
        # (that literal used to land in evidence_quote for every quote-less fact).
        return None
    return _clean_optional_text(quotes)


def _normalize_count_provenance(raw: Any) -> Dict[str, Any]:
    if isinstance(raw, dict):
        return raw
    if isinstance(raw, str) and raw.strip():
        try:
            parsed = json.loads(raw)
        except json.JSONDecodeError:
            return {}
        return parsed if isinstance(parsed, dict) else {}
    return {}


def _infer_source_components(
    source_location: Optional[str],
) -> Dict[str, Optional[str]]:
    """Best-effort split of a free-text source pointer into queryable pieces."""
    text = source_location or ""
    table = None
    row = None
    paragraph = None
    section = None

    table_match = re.search(
        r"\b(?:(?:Supplementary|Supplemental)\s+)?e?Table\s+[A-Za-z0-9._-]+",
        text,
        re.IGNORECASE,
    )
    if table_match:
        table = table_match.group(0)

    row_match = re.search(
        r"\brows?\s+[A-Za-z0-9._-]+(?:\s*(?:-|–|,|and)\s*[A-Za-z0-9._-]+)?",
        text,
        re.IGNORECASE,
    )
    if row_match:
        row = row_match.group(0)

    para_match = re.search(
        r"\b(?:paragraph|para\.?)\s+[A-Za-z0-9._-]+",
        text,
        re.IGNORECASE,
    )
    if para_match:
        paragraph = para_match.group(0)

    section_match = re.search(
        r"\b(Abstract|Introduction|Methods?|Results?|Discussion|Conclusion|"
        r"Case(?: presentation| report)?|Supplementary Material)\b",
        text,
        re.IGNORECASE,
    )
    if section_match:
        section = section_match.group(0)

    return {
        "source_table": _clean_optional_text(table),
        "source_row": _clean_optional_text(row),
        "source_paragraph": _clean_optional_text(paragraph),
        "source_section": _clean_optional_text(section),
    }


def _nested_provenance(raw: Dict[str, Any]) -> Dict[str, Any]:
    """Return flat per-observation provenance keys from a possibly nested row."""
    merged: Dict[str, Any] = {}
    for key in ("provenance", "source_provenance", "source_locator", "locator"):
        nested = raw.get(key)
        if isinstance(nested, dict):
            merged.update(nested)
    merged.update(raw)
    return merged


def _normalize_source_container(value: Any) -> Optional[str]:
    text = _clean_optional_text(value)
    if not text:
        return None
    lowered = text.lower()
    if lowered in {"supplement", "supplemental", "supplementary", "supp"}:
        return "supplement"
    if lowered == "main":
        return "main"
    return None


def _normalize_source_kind(value: Any) -> Optional[str]:
    text = _clean_optional_text(value)
    if not text:
        return None
    lowered = text.lower()
    if lowered in {"table", "figure", "text", "abstract"}:
        return lowered
    if lowered in {"section", "paragraph", "body"}:
        return "text"
    return None


def _coerce_optional_int(value: Any) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    text = _clean_optional_text(value)
    if not text:
        return None
    match = re.search(r"\d+", text)
    return int(match.group(0)) if match else None


def _infer_observation_source_from_location(
    source_location: Any,
) -> Dict[str, Optional[Any]]:
    """Best-effort structured locator from legacy free-text source_location."""
    text = _clean_optional_text(source_location) or ""
    if not text:
        return {}
    lowered = text.lower()

    source_container = (
        "supplement"
        if re.search(
            r"\bsupp(?:lement(?:ary|al)?)?\b|etable|e-table|table\s+s\d", lowered
        )
        else "main"
    )
    source_kind = None
    if re.search(r"\b(?:e?table|tab\.?)\b", lowered):
        source_kind = "table"
    elif re.search(r"\b(?:fig(?:ure)?\.?)\b", lowered):
        source_kind = "figure"
    elif "abstract" in lowered:
        source_kind = "abstract"
    elif re.search(
        r"\b(results?|methods?|discussion|conclusions?|case(?: presentation| report)?|paragraph|text scan)\b",
        lowered,
    ):
        source_kind = "text"

    ref = None
    figure_panel = None
    ref_match = re.search(
        r"\b(?:(?:Supplementary|Supplemental|Supp\.?)\s+)?e?"
        r"(?:Table|Fig(?:ure)?\.?)\s*[A-Za-z0-9][A-Za-z0-9._-]*[A-Za-z]?",
        text,
        re.IGNORECASE,
    )
    if ref_match:
        ref = ref_match.group(0).strip(" .;,")
        if source_kind == "figure":
            panel_match = re.search(r"\d+\s*([A-Za-z])$", ref)
            if panel_match:
                figure_panel = panel_match.group(1).upper()
    elif source_kind in {"table", "figure"}:
        ref = text
    else:
        section_match = re.search(
            r"\b(Abstract|Introduction|Methods?|Results?|Discussion|Conclusions?|"
            r"Case(?: presentation| report)?|Supplementary Material)\b",
            text,
            re.IGNORECASE,
        )
        if section_match:
            ref = section_match.group(0)

    row_ordinal = None
    row_match = re.search(r"\brows?\s+(\d+)\b", text, re.IGNORECASE)
    if row_match:
        row_ordinal = int(row_match.group(1))

    return {
        "source_container": source_container,
        "source_kind": source_kind,
        "source_ref": _clean_optional_text(ref),
        "row_ordinal": row_ordinal,
        "figure_panel": figure_panel,
    }


def _serialize_locator_extra(raw: Any, provenance: Dict[str, Any]) -> Optional[str]:
    value = raw
    if value is None:
        extra: Dict[str, Any] = {}
        for out_key, in_keys in {
            "bbox": ("bbox", "bounding_box"),
            "pmc_xpath": ("pmc_xpath", "xpath"),
            "cell_coords": ("cell_coords", "cell_coordinates"),
            "parser_confidence": ("parser_confidence", "confidence"),
        }.items():
            for in_key in in_keys:
                if provenance.get(in_key) is not None:
                    extra[out_key] = provenance[in_key]
                    break
        value = extra or None
    if value is None:
        return None
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        try:
            value = json.loads(text)
        except json.JSONDecodeError:
            value = text
    return json.dumps(value, sort_keys=True, ensure_ascii=False)


def _stable_source_record_id(
    *,
    pmid: Any,
    source_notation: Any,
    individual_id: Any = None,
    row_label: Any = None,
    source_ref: Any = None,
) -> str:
    person_key = _clean_optional_text(individual_id) or _clean_optional_text(row_label)
    payload = {
        "pmid": _clean_optional_text(pmid) or "",
        "protein_notation": _clean_optional_text(source_notation) or "",
        "person": person_key or "",
        "source_ref": _clean_optional_text(source_ref) or "",
    }
    encoded = json.dumps(
        {k: v.casefold() for k, v in payload.items()},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _observation_provenance_values(
    *,
    observation: Dict[str, Any],
    variant_data: Dict[str, Any],
    pmid: str,
    individual_id: Any = None,
    default_column_ref: Any = None,
) -> Dict[str, Any]:
    """Normalize structured per-observation provenance for DB insertion."""
    obs_prov = _nested_provenance(observation)
    variant_prov = _nested_provenance(variant_data)
    source_location = (
        obs_prov.get("source_location")
        or variant_prov.get("source_location")
        or variant_data.get("source_location")
    )
    inferred = _infer_observation_source_from_location(source_location)

    def pick(*keys: str, default: Any = None) -> Any:
        for key in keys:
            if obs_prov.get(key) is not None:
                return obs_prov.get(key)
        for key in keys:
            if variant_prov.get(key) is not None:
                return variant_prov.get(key)
        return inferred.get(keys[0], default)

    source_container = _normalize_source_container(
        pick("source_container", "container")
    )
    source_kind = _normalize_source_kind(pick("source_kind", "kind"))
    source_ref = _clean_optional_text(
        pick(
            "source_ref",
            "source_table",
            "table",
            "source_figure",
            "figure",
            "source_section",
            "section",
        )
    )
    row_label = _clean_optional_text(
        pick("row_label", "source_row", "row", "case_label")
    )
    if not row_label:
        row_label = _clean_optional_text(individual_id)
    column_ref = _clean_optional_text(
        pick("column_ref", "source_column", "column", default=default_column_ref)
    )
    figure_panel = _clean_optional_text(pick("figure_panel", "panel"))
    if not figure_panel and source_kind == "figure" and source_ref:
        panel_match = re.search(r"\d+\s*([A-Za-z])$", source_ref)
        if panel_match:
            figure_panel = panel_match.group(1).upper()

    source_notation = (
        variant_data.get("protein_notation")
        or variant_data.get("cdna_notation")
        or variant_data.get("genomic_position")
    )
    return {
        "source_container": source_container,
        "source_kind": source_kind,
        "source_ref": source_ref,
        "page_label": _clean_optional_text(
            pick("page_label", "printed_page", "source_page")
        ),
        "pdf_page": _coerce_optional_int(pick("pdf_page", "pdf_page_number")),
        "row_label": row_label,
        "row_ordinal": _coerce_optional_int(pick("row_ordinal", "row_index")),
        "column_ref": column_ref,
        "figure_panel": figure_panel,
        "source_record_id": _stable_source_record_id(
            pmid=pmid,
            source_notation=source_notation,
            individual_id=individual_id,
            row_label=row_label,
            source_ref=source_ref,
        ),
        "locator_extra": _serialize_locator_extra(
            obs_prov.get("locator_extra") or variant_prov.get("locator_extra"),
            obs_prov,
        ),
    }


def _fact_hash(row: Dict[str, Any]) -> str:
    payload = {
        key: _clean_optional_text(row.get(key)) for key in FACT_PROVENANCE_FIELDS
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def insert_fact_provenance(
    cursor: sqlite3.Cursor,
    *,
    variant_id: int,
    pmid: str,
    fact_type: str,
    fact_value: Any = None,
    individual_id: Any = None,
    source_location: Any = None,
    source_section: Any = None,
    source_paragraph: Any = None,
    source_table: Any = None,
    source_row: Any = None,
    source_column: Any = None,
    evidence_quote: Any = None,
    count_type: Any = None,
    source_layer: Any = None,
    provenance_kind: str = "synthesized",
) -> None:
    """Insert one fact-level provenance row, idempotently.

    ``fact_hash`` is a stable natural key over the semantic fields. It makes
    re-migration, refresh replay, and layered recovery additive without creating
    duplicate evidence rows.
    """
    source_location_text = _clean_optional_text(source_location)
    inferred = _infer_source_components(source_location_text)
    row = {
        "variant_id": variant_id,
        "pmid": _clean_optional_text(pmid),
        "fact_type": _clean_optional_text(fact_type),
        "fact_value": _clean_optional_text(fact_value),
        "individual_id": _clean_optional_text(individual_id),
        "source_location": source_location_text,
        "source_section": _clean_optional_text(source_section)
        or inferred["source_section"],
        "source_paragraph": _clean_optional_text(source_paragraph)
        or inferred["source_paragraph"],
        "source_table": _clean_optional_text(source_table) or inferred["source_table"],
        "source_row": _clean_optional_text(source_row) or inferred["source_row"],
        "source_column": _clean_optional_text(source_column),
        "evidence_quote": _clean_optional_text(evidence_quote),
        "count_type": _clean_optional_text(count_type),
        "source_layer": _clean_optional_text(source_layer),
        "provenance_kind": _clean_optional_text(provenance_kind),
    }
    if not row["fact_type"]:
        return
    row["fact_hash"] = _fact_hash(row)
    cursor.execute(
        """
        INSERT OR IGNORE INTO fact_provenance (
            fact_hash, variant_id, pmid, fact_type, fact_value, individual_id,
            source_location, source_section, source_paragraph, source_table,
            source_row, source_column, evidence_quote, count_type, source_layer,
            provenance_kind
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            row["fact_hash"],
            variant_id,
            row["pmid"],
            row["fact_type"],
            row["fact_value"],
            row["individual_id"],
            row["source_location"],
            row["source_section"],
            row["source_paragraph"],
            row["source_table"],
            row["source_row"],
            row["source_column"],
            row["evidence_quote"],
            row["count_type"],
            row["source_layer"],
            row["provenance_kind"],
        ),
    )


def _count_provenance_for(
    count_provenance: Dict[str, Any], fact_type: str
) -> Tuple[Optional[str], Optional[str]]:
    prefix = {
        "patient_count": "carriers",
        "total_carriers_observed": "carriers",
        "affected_count": "affected",
        "unaffected_count": "unaffected",
    }.get(fact_type)
    if not prefix:
        return None, None
    return (
        _clean_optional_text(count_provenance.get(f"{prefix}_column_label")),
        _clean_optional_text(count_provenance.get(f"{prefix}_count_type")),
    )


def _insert_explicit_fact_provenance(
    cursor: sqlite3.Cursor,
    *,
    pmid: str,
    variant_id: int,
    variant_data: Dict[str, Any],
    source_layer: str,
) -> None:
    raw_rows = variant_data.get("fact_provenance") or variant_data.get(
        "fact_provenance_records"
    )
    if not isinstance(raw_rows, list):
        return
    for raw in raw_rows:
        if not isinstance(raw, dict):
            logger.warning(
                "Skipping non-object fact_provenance row for PMID %s variant_id %s",
                pmid,
                variant_id,
            )
            continue
        insert_fact_provenance(
            cursor,
            variant_id=variant_id,
            pmid=pmid,
            fact_type=raw.get("fact_type") or raw.get("field"),
            fact_value=raw.get("fact_value")
            if "fact_value" in raw
            else raw.get("value"),
            individual_id=raw.get("individual_id"),
            source_location=raw.get("source_location")
            or variant_data.get("source_location"),
            source_section=raw.get("source_section") or raw.get("section"),
            source_paragraph=raw.get("source_paragraph") or raw.get("paragraph"),
            source_table=raw.get("source_table") or raw.get("table"),
            source_row=raw.get("source_row") or raw.get("row"),
            source_column=raw.get("source_column") or raw.get("column"),
            evidence_quote=raw.get("evidence_quote")
            or raw.get("quote")
            or raw.get("evidence_sentence"),
            count_type=raw.get("count_type"),
            source_layer=raw.get("source_layer") or source_layer,
            provenance_kind="explicit",
        )


def _insert_standard_fact_provenance(
    cursor: sqlite3.Cursor,
    *,
    pmid: str,
    variant_id: int,
    variant_data: Dict[str, Any],
    source_layer: str,
) -> None:
    """Write the repo-standard provenance rows for every migrated variant."""
    source_location = variant_data.get("source_location")
    quote = _first_quote(variant_data)

    _insert_explicit_fact_provenance(
        cursor,
        pmid=pmid,
        variant_id=variant_id,
        variant_data=variant_data,
        source_layer=source_layer,
    )

    # Structured locators the extractor may attach at the variant top level
    # (e.g. the table-regex pass records source_table/source_row/source_column +
    # the verbatim row). Thread them into the synthesized identity fact so the
    # within-paper coordinate survives to fact_provenance.
    top_source_table = variant_data.get("source_table")
    top_source_row = variant_data.get("source_row")
    top_source_column = variant_data.get("source_column")
    top_evidence_quote = variant_data.get("evidence_quote") or quote

    for notation_key in ("protein_notation", "cdna_notation", "genomic_position"):
        notation = variant_data.get(notation_key)
        if notation:
            insert_fact_provenance(
                cursor,
                variant_id=variant_id,
                pmid=pmid,
                fact_type="variant_identity",
                fact_value=notation,
                source_location=source_location,
                source_table=top_source_table,
                source_row=top_source_row,
                source_column=top_source_column,
                evidence_quote=top_evidence_quote,
                source_layer=source_layer,
                provenance_kind="synthesized",
            )

    if variant_data.get("clinical_significance") is not None:
        insert_fact_provenance(
            cursor,
            variant_id=variant_id,
            pmid=pmid,
            fact_type="clinical_significance",
            fact_value=variant_data.get("clinical_significance"),
            source_location=source_location,
            evidence_quote=quote,
            source_layer=source_layer,
            provenance_kind="synthesized",
        )

    count_provenance = _normalize_count_provenance(variant_data.get("count_provenance"))
    patients = variant_data.get("patients") or {}
    if patients.get("count") is not None:
        column, count_type = _count_provenance_for(count_provenance, "patient_count")
        insert_fact_provenance(
            cursor,
            variant_id=variant_id,
            pmid=pmid,
            fact_type="patient_count",
            fact_value=patients.get("count"),
            source_location=source_location,
            source_column=column,
            evidence_quote=quote,
            count_type=count_type,
            source_layer=source_layer,
            provenance_kind="synthesized",
        )

    penetrance = variant_data.get("penetrance_data") or {}
    for fact_type in (
        "total_carriers_observed",
        "affected_count",
        "unaffected_count",
        "uncertain_count",
    ):
        if penetrance.get(fact_type) is None:
            continue
        column, count_type = _count_provenance_for(count_provenance, fact_type)
        insert_fact_provenance(
            cursor,
            variant_id=variant_id,
            pmid=pmid,
            fact_type=fact_type,
            fact_value=penetrance.get(fact_type),
            source_location=source_location,
            source_column=column,
            evidence_quote=quote,
            count_type=count_type,
            source_layer=source_layer,
            provenance_kind="synthesized",
        )

    for record in variant_data.get("individual_records", []):
        if not isinstance(record, dict):
            continue
        individual_id = record.get("individual_id")
        evidence = record.get("evidence_sentence") or quote
        insert_fact_provenance(
            cursor,
            variant_id=variant_id,
            pmid=pmid,
            fact_type="individual_carrier_status",
            fact_value="carrier",
            individual_id=individual_id,
            source_location=source_location,
            evidence_quote=evidence,
            source_layer=source_layer,
            provenance_kind="synthesized",
        )
        affected_status = normalize_affected_status(record.get("affected_status"))
        if affected_status:
            insert_fact_provenance(
                cursor,
                variant_id=variant_id,
                pmid=pmid,
                fact_type="individual_affected_status",
                fact_value=affected_status,
                individual_id=individual_id,
                source_location=source_location,
                evidence_quote=evidence,
                source_layer=source_layer,
                provenance_kind="synthesized",
            )


def insert_variant_data(
    cursor: sqlite3.Cursor,
    pmid: str,
    variant_data: Dict[str, Any],
    preserve_existing_evidence: bool = False,
) -> Optional[int]:
    """
    Insert variant and all associated data.

    Args:
        cursor: Database cursor
        pmid: Paper PMID
        variant_data: Variant data dictionary
        preserve_existing_evidence: If True, do not duplicate existing
            penetrance or individual rows for the same PMID/variant.

    Returns:
        variant_id, or None if the variant had no usable notation
    """
    source_layer = infer_source_layer(variant_data)
    junk_reason = junk_notation_reason(
        source_layer=source_layer,
        protein_notation=variant_data.get("protein_notation"),
        cdna_notation=variant_data.get("cdna_notation"),
        variant=(
            variant_data.get("protein_notation")
            or variant_data.get("cdna_notation")
            or variant_data.get("genomic_position")
        ),
        gene_symbol=variant_data.get("gene_symbol"),
    )
    if junk_reason:
        logger.info(
            "Dropping junk %s variant in PMID %s "
            "(reason=%s, gene=%s, protein=%r, cdna=%r)",
            source_layer,
            pmid,
            junk_reason,
            variant_data.get("gene_symbol"),
            variant_data.get("protein_notation"),
            variant_data.get("cdna_notation"),
        )
        return None

    original_notation = {
        "protein_notation": variant_data.get("protein_notation"),
        "cdna_notation": variant_data.get("cdna_notation"),
        "genomic_position": variant_data.get("genomic_position"),
    }
    if not sanitize_variant_notation(variant_data):
        # Promote to WARNING so cohort-summary hallucinations are visible in
        # run summaries — Tier 3 sometimes emits group-level entries with
        # all-null notation (e.g., "34 patients with pore-region mutations").
        # The prompt rejects these but model drift can still produce them.
        patients = variant_data.get("patients") or {}
        logger.warning(
            "Dropping variant with no usable notation in PMID %s "
            "(gene=%s, patients.count=%s, clinical_significance=%s): %s",
            pmid,
            variant_data.get("gene_symbol"),
            patients.get("count"),
            variant_data.get("clinical_significance"),
            original_notation,
        )
        return None

    # Get or create variant
    variant_id = get_or_create_variant(cursor, variant_data)

    # Insert variant-paper association
    cursor.execute(
        """
        INSERT OR IGNORE INTO variant_papers (
            variant_id, pmid, source_location, additional_notes, key_quotes,
            count_provenance, source_layer
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
    """,
        (
            variant_id,
            pmid,
            variant_data.get("source_location"),
            variant_data.get("additional_notes"),
            json.dumps(variant_data.get("key_quotes", [])),
            json.dumps(variant_data["count_provenance"])
            if variant_data.get("count_provenance")
            else None,
            source_layer,
        ),
    )
    cursor.execute(
        """
        UPDATE variant_papers
        SET source_layer = ?
        WHERE variant_id = ?
          AND pmid = ?
          AND (source_layer IS NULL OR TRIM(source_layer) = '')
        """,
        (source_layer, variant_id, pmid),
    )

    _insert_standard_fact_provenance(
        cursor,
        pmid=pmid,
        variant_id=variant_id,
        variant_data=variant_data,
        source_layer=source_layer,
    )

    # Insert penetrance data
    penetrance = variant_data.get("penetrance_data", {})
    if penetrance and any(
        penetrance.get(k) is not None
        for k in ["total_carriers_observed", "affected_count", "unaffected_count"]
    ):
        penetrance_values = (
            variant_id,
            pmid,
            penetrance.get("total_carriers_observed"),
            penetrance.get("affected_count"),
            penetrance.get("unaffected_count"),
            penetrance.get("uncertain_count"),
            penetrance.get("penetrance_percentage"),
        )
        # Idempotency guard: an identical penetrance row for this (variant_id,
        # pmid) must never be written twice. Exact duplicates come from one
        # variant represented across several supplement-table cells (multi-cohort
        # / multi-classification) or from re-migrating a paper on a refresh;
        # without this guard they are summed N-fold into the carrier/affected
        # counts at scoring time (the regression that took KCNQ1 PMID 32893267
        # to 4x gold).
        exact_penetrance_match = {
            "variant_id": variant_id,
            "pmid": pmid,
            "total_carriers_observed": penetrance.get("total_carriers_observed"),
            "affected_count": penetrance.get("affected_count"),
            "unaffected_count": penetrance.get("unaffected_count"),
            "uncertain_count": penetrance.get("uncertain_count"),
            "penetrance_percentage": penetrance.get("penetrance_percentage"),
        }
        exact_penetrance_id = _matching_row_id(
            cursor,
            "penetrance_data",
            "penetrance_id",
            exact_penetrance_match,
        )
        is_exact_penetrance_dup = exact_penetrance_id is not None

        existing_penetrance = None
        if preserve_existing_evidence and not is_exact_penetrance_dup:
            existing_penetrance = cursor.execute(
                """
                SELECT penetrance_id, total_carriers_observed, affected_count,
                       unaffected_count, uncertain_count, penetrance_percentage
                FROM penetrance_data
                WHERE variant_id = ? AND pmid = ?
                ORDER BY penetrance_id
                LIMIT 1
                """,
                (variant_id, pmid),
            ).fetchone()

        if is_exact_penetrance_dup:
            # Identical parent already stored. Still merge any unique
            # age-stratified children carried by this duplicate representation.
            _insert_age_dependent_penetrance(
                cursor,
                exact_penetrance_id,
                penetrance.get("age_dependent_penetrance", []),
            )
        elif existing_penetrance:
            penetrance_id = existing_penetrance[0]
            current_values = existing_penetrance[1:]
            incoming_values = penetrance_values[2:]
            merged_values = tuple(
                incoming if current is None and incoming is not None else current
                for current, incoming in zip(current_values, incoming_values)
            )
            if merged_values != current_values:
                cursor.execute(
                    """
                    UPDATE penetrance_data
                    SET total_carriers_observed = ?,
                        affected_count = ?,
                        unaffected_count = ?,
                        uncertain_count = ?,
                        penetrance_percentage = ?
                    WHERE penetrance_id = ?
                    """,
                    (*merged_values, penetrance_id),
                )
            _insert_age_dependent_penetrance(
                cursor,
                penetrance_id,
                penetrance.get("age_dependent_penetrance", []),
            )
        else:
            cursor.execute(
                """
                INSERT INTO penetrance_data (
                    variant_id, pmid, total_carriers_observed, affected_count,
                    unaffected_count, uncertain_count, penetrance_percentage
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
                penetrance_values,
            )

            penetrance_id = cursor.lastrowid

            # Insert age-dependent penetrance only for a newly inserted parent
            # row; duplicates are folded into the stored parent above.
            _insert_age_dependent_penetrance(
                cursor,
                penetrance_id,
                penetrance.get("age_dependent_penetrance", []),
            )

    # Insert individual records
    for record in variant_data.get("individual_records", []):
        if not isinstance(record, dict):
            logger.warning(
                "Skipping non-object individual_record for PMID %s variant_id %s",
                pmid,
                variant_id,
            )
            continue
        affected_status = normalize_affected_status(record.get("affected_status"))
        individual_id = record.get("individual_id")
        individual_provenance = _observation_provenance_values(
            observation=record,
            variant_data=variant_data,
            pmid=pmid,
            individual_id=individual_id,
        )
        # Deterministic ethnicity / cohort-origin backfill from the surrounding
        # patient text when the extractor didn't set them. Author-affiliation
        # fallback is applied later (network) via the metadata backfill, so migrate
        # stays offline.
        patients_block = variant_data.get("patients") or {}
        record_ethnicity, record_origin = enrich_ethnicity_origin(
            stated_ethnicity=record.get("ethnicity"),
            stated_origin=record.get("geographic_origin"),
            text_fields=[
                record.get("phenotype_details"),
                record.get("evidence_sentence"),
                patients_block.get("demographics"),
                variant_data.get("additional_notes"),
            ],
        )
        if preserve_existing_evidence:
            if individual_id:
                existing_record = cursor.execute(
                    """
                    SELECT 1 FROM individual_records
                    WHERE variant_id = ? AND pmid = ? AND individual_id = ?
                    LIMIT 1
                    """,
                    (variant_id, pmid, individual_id),
                ).fetchone()
            else:
                existing_record = cursor.execute(
                    """
                    SELECT 1 FROM individual_records
                    WHERE variant_id = ? AND pmid = ?
                    LIMIT 1
                    """,
                    (variant_id, pmid),
                ).fetchone()
            if existing_record:
                continue
        # Idempotency guard: skip an individual record identical to one already
        # stored for this (variant_id, pmid) (applies in replace mode too).
        if _row_exists(
            cursor,
            "individual_records",
            {
                "variant_id": variant_id,
                "pmid": pmid,
                "individual_id": individual_id,
                "age_at_evaluation": record.get("age_at_evaluation"),
                "age_at_onset": record.get("age_at_onset"),
                "age_at_diagnosis": record.get("age_at_diagnosis"),
                "sex": record.get("sex"),
                "affected_status": affected_status,
                "phenotype_details": record.get("phenotype_details"),
                "evidence_sentence": record.get("evidence_sentence"),
                "ethnicity": record_ethnicity,
                "geographic_origin": record_origin,
                **individual_provenance,
            },
        ):
            continue
        cursor.execute(
            """
            INSERT INTO individual_records (
                variant_id, pmid, individual_id, age_at_evaluation,
                age_at_onset, age_at_diagnosis, sex, affected_status,
                phenotype_details, evidence_sentence, ethnicity, geographic_origin,
                source_container, source_kind, source_ref, page_label, pdf_page,
                row_label, row_ordinal, column_ref, figure_panel, source_record_id,
                locator_extra
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
            (
                variant_id,
                pmid,
                individual_id,
                record.get("age_at_evaluation"),
                record.get("age_at_onset"),
                record.get("age_at_diagnosis"),
                record.get("sex"),
                affected_status,
                record.get("phenotype_details"),
                record.get("evidence_sentence"),
                record_ethnicity,
                record_origin,
                individual_provenance["source_container"],
                individual_provenance["source_kind"],
                individual_provenance["source_ref"],
                individual_provenance["page_label"],
                individual_provenance["pdf_page"],
                individual_provenance["row_label"],
                individual_provenance["row_ordinal"],
                individual_provenance["column_ref"],
                individual_provenance["figure_panel"],
                individual_provenance["source_record_id"],
                individual_provenance["locator_extra"],
            ),
        )

    # Insert functional data (exact-dup guard keeps re-migration idempotent —
    # _land does not pre-delete these child rows the way it does penetrance).
    functional = variant_data.get("functional_data", {})
    if functional and (functional.get("summary") or functional.get("assays")):
        assays_json = json.dumps(functional.get("assays", []))
        if not _row_exists(
            cursor,
            "functional_data",
            {
                "variant_id": variant_id,
                "pmid": pmid,
                "summary": functional.get("summary"),
                "assays": assays_json,
            },
        ):
            cursor.execute(
                """
                INSERT INTO functional_data (
                    variant_id, pmid, summary, assays
                ) VALUES (?, ?, ?, ?)
            """,
                (variant_id, pmid, functional.get("summary"), assays_json),
            )

    # Insert phenotype data
    patients = variant_data.get("patients", {})
    if patients and (patients.get("count") or patients.get("phenotype")):
        count_provenance = _normalize_count_provenance(
            variant_data.get("count_provenance")
        )
        phenotype_provenance = _observation_provenance_values(
            observation=patients,
            variant_data=variant_data,
            pmid=pmid,
            individual_id=patients.get("individual_id"),
            default_column_ref=count_provenance.get("carriers_column_label"),
        )
        if not _row_exists(
            cursor,
            "phenotypes",
            {
                "variant_id": variant_id,
                "pmid": pmid,
                "patient_count": patients.get("count"),
                "demographics": patients.get("demographics"),
                "phenotype_description": patients.get("phenotype"),
                **phenotype_provenance,
            },
        ):
            cursor.execute(
                """
                INSERT INTO phenotypes (
                    variant_id, pmid, patient_count, demographics,
                    phenotype_description, source_container, source_kind, source_ref,
                    page_label, pdf_page, row_label, row_ordinal, column_ref,
                    figure_panel, source_record_id, locator_extra
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
                (
                    variant_id,
                    pmid,
                    patients.get("count"),
                    patients.get("demographics"),
                    patients.get("phenotype"),
                    phenotype_provenance["source_container"],
                    phenotype_provenance["source_kind"],
                    phenotype_provenance["source_ref"],
                    phenotype_provenance["page_label"],
                    phenotype_provenance["pdf_page"],
                    phenotype_provenance["row_label"],
                    phenotype_provenance["row_ordinal"],
                    phenotype_provenance["column_ref"],
                    phenotype_provenance["figure_panel"],
                    phenotype_provenance["source_record_id"],
                    phenotype_provenance["locator_extra"],
                ),
            )

    # Insert variant metadata
    if variant_data.get("segregation_data") or variant_data.get("population_frequency"):
        if not _row_exists(
            cursor,
            "variant_metadata",
            {
                "variant_id": variant_id,
                "pmid": pmid,
                "segregation_data": variant_data.get("segregation_data"),
                "population_frequency": variant_data.get("population_frequency"),
            },
        ):
            cursor.execute(
                """
                INSERT INTO variant_metadata (
                    variant_id, pmid, segregation_data, population_frequency
                ) VALUES (?, ?, ?, ?)
            """,
                (
                    variant_id,
                    pmid,
                    variant_data.get("segregation_data"),
                    variant_data.get("population_frequency"),
                ),
            )

    return variant_id


def _resolve_reference_validation_policy() -> str:
    """Return the reference-validation policy from Settings (default 'off').

    Falls back to 'off' (strict no-op) if Settings can't load, mirroring the
    count-hygiene policy resolution in pipeline/steps.py.
    """
    try:
        from config.settings import get_settings

        return getattr(get_settings(), "reference_validation_policy", "off") or "off"
    except Exception:
        return "off"


def _apply_reference_validation(
    variants: List[Dict[str, Any]], *, policy: str
) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    """Apply the B1 reference-transcript gate to a variant list before insert.

    Validates each variant's protein notation against its gene's cached reference
    protein (``pipeline.reference_validation.validate_protein_variant``). Behavior
    by policy:
      * ``off``  — strict no-op; returns the list unchanged (imports nothing).
      * ``flag`` — annotate every rejected variant with a ``reference_validation``
        block (status/reason/expected/observed/high_confidence) and keep it.
      * ``drop`` — annotate every reject, but only *omit* HIGH-CONFIDENCE rejects
        (clean substitutions/nonsense beyond the N-terminal isoform window).
        Frameshift/indel and N-terminal rejects are annotated but kept, because
        their reference residue is notation-convention / isoform ambiguous and a
        hard drop there loses real variants (measured: 7 of 8 KCNH2 recall losses
        were exactly this class).

    A variant is only ever ``reject`` when its gene has a cached reference
    sequence AND the notation names a concrete reference residue+position that
    mismatches (or is out of range). cDNA-only / uncached-gene variants resolve to
    ``unknown`` and always pass, so the gate is turnkey-safe on new genes.

    Returns ``(kept_variants, stats)`` where stats counts ok/unknown/reject plus
    flagged/dropped.
    """
    stats = {"ok": 0, "unknown": 0, "reject": 0, "flagged": 0, "dropped": 0}
    if policy == "off" or not variants:
        return variants, stats

    from pipeline.reference_validation import (
        is_high_confidence_reject,
        validate_protein_variant,
    )

    kept: List[Dict[str, Any]] = []
    for variant in variants:
        gene = variant.get("gene_symbol") or ""
        protein = variant.get("protein_notation") or ""
        verdict = validate_protein_variant(gene, protein)
        stats[verdict.status] = stats.get(verdict.status, 0) + 1
        if not verdict.is_reject:
            kept.append(variant)
            continue
        # Drop only high-confidence rejects; flag mode never drops.
        high_conf = is_high_confidence_reject(protein, verdict)
        will_drop = policy == "drop" and high_conf
        variant["reference_validation"] = {
            "status": verdict.status,
            "reason": verdict.reason,
            "expected": verdict.expected,
            "observed": verdict.observed,
            "policy": policy,
            "high_confidence": high_conf,
            "dropped": will_drop,
        }
        stats["flagged"] += 1
        if will_drop:
            stats["dropped"] += 1
            continue
        kept.append(variant)
    return kept, stats


def migrate_extraction_file(
    cursor: sqlite3.Cursor,
    json_file: Path,
    auto_repair: bool = True,
    replace_existing_paper: bool = True,
) -> Tuple[bool, str]:
    """
    Migrate a single extraction JSON file to the database.

    Args:
        cursor: Database cursor
        json_file: Path to JSON file
        auto_repair: If True, attempt to repair common issues before migration
        replace_existing_paper: If True, use the legacy paper replace behavior.
            If False, preserve existing PMID-linked evidence and add new data.

    Returns:
        Tuple of (success, message)
    """
    try:
        with open(json_file, "r", encoding="utf-8") as f:
            extraction_data = json.load(f)

        # Validate before migration
        is_valid, errors, warnings = validate_extraction_data(
            extraction_data, json_file.name
        )

        # Log warnings
        for w in warnings:
            logger.warning(f"{json_file.name}: {w}")

        # Attempt repair if invalid and auto_repair is enabled
        if not is_valid and auto_repair:
            logger.info(f"Attempting auto-repair for {json_file.name}: {errors}")
            extraction_data, repairs = repair_extraction_data(
                extraction_data, json_file.name
            )
            if repairs:
                logger.info(f"Applied {len(repairs)} repairs: {repairs}")

            # Re-validate after repair
            is_valid, errors, _ = validate_extraction_data(
                extraction_data, json_file.name
            )

        # If still invalid after repair, fail
        if not is_valid:
            return False, f"Validation failed for {json_file.name}: {errors}"

        # Insert paper metadata
        insert_paper_metadata(
            cursor, extraction_data, replace_existing=replace_existing_paper
        )

        pmid = extraction_data.get("paper_metadata", {}).get("pmid", "UNKNOWN")

        # Insert variants
        variants = extraction_data.get("variants", [])
        # B1 reference-transcript validation gate. Strict no-op unless
        # REFERENCE_VALIDATION_POLICY is flag/drop; gold-free and turnkey-safe
        # (cDNA-only notations and uncached genes always pass).
        refval_policy = _resolve_reference_validation_policy()
        if refval_policy != "off":
            variants, refval_stats = _apply_reference_validation(
                variants, policy=refval_policy
            )
            if refval_stats["flagged"]:
                logger.info(
                    "%s: reference-validation (%s) flagged %d variant(s), "
                    "dropped %d (ok=%d unknown=%d reject=%d)",
                    json_file.name,
                    refval_policy,
                    refval_stats["flagged"],
                    refval_stats["dropped"],
                    refval_stats["ok"],
                    refval_stats["unknown"],
                    refval_stats["reject"],
                )
        for variant in variants:
            insert_variant_data(
                cursor,
                pmid,
                variant,
                preserve_existing_evidence=not replace_existing_paper,
            )

        # Insert extraction metadata
        extraction_meta = extraction_data.get("extraction_metadata", {})
        # Determine source_file: prefer from metadata, fallback to extraction JSON filename
        source_file = extraction_meta.get("source_file")
        if not source_file:
            # Fallback: use the JSON file path itself as the source reference
            source_file = str(json_file)

        # Always insert extraction metadata (at minimum we have pmid and source_file)
        if extraction_meta or source_file:
            cursor.execute(
                """
                INSERT INTO extraction_metadata (
                    pmid, total_variants_found, extraction_confidence,
                    study_type, study_design, ascertainment, cohort_source,
                    population, study_summary,
                    challenges, notes, extraction_timestamp, source_type, abstract_only,
                    source_file, model_used
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
                (
                    pmid,
                    extraction_meta.get("total_variants_found"),
                    extraction_meta.get("extraction_confidence"),
                    extraction_meta.get("study_type"),
                    extraction_meta.get("study_design"),
                    extraction_meta.get("ascertainment"),
                    extraction_meta.get("cohort_source"),
                    extraction_meta.get("population"),
                    extraction_meta.get("study_summary"),
                    json.dumps(extraction_meta.get("challenges", [])),
                    extraction_meta.get("notes"),
                    datetime.now().isoformat(),
                    extraction_meta.get("source_type"),
                    1 if extraction_meta.get("abstract_only") else 0,
                    source_file,
                    extraction_meta.get("model_used"),
                ),
            )

        # Insert tables processed
        for table in extraction_data.get("tables_processed", []):
            cursor.execute(
                """
                INSERT INTO tables_processed (
                    pmid, table_name, table_caption, variants_extracted
                ) VALUES (?, ?, ?, ?)
            """,
                (
                    pmid,
                    table.get("table_name"),
                    table.get("table_caption"),
                    table.get("variants_extracted"),
                ),
            )

        return True, f"Successfully migrated {json_file.name}"

    except Exception as e:
        error_msg = f"Failed to migrate {json_file.name}: {str(e)}"
        logger.error(error_msg)
        return False, error_msg


def migrate_extraction_directory(
    conn: sqlite3.Connection,
    extraction_dir: Path,
    replace_existing_paper: bool = True,
) -> Dict[str, Any]:
    """
    Migrate all extraction JSON files from a directory.

    Args:
        conn: Database connection
        extraction_dir: Directory containing extraction JSON files
        replace_existing_paper: If True, use the legacy paper replace behavior.
            If False, preserve existing PMID-linked evidence and add new data.

    Returns:
        Migration statistics
    """
    logger.info(f"Migrating extraction files from: {extraction_dir}")

    # Find all JSON files matching extraction pattern
    json_files = find_extraction_json_files(extraction_dir)

    if not json_files:
        logger.warning(f"No extraction JSON files found in {extraction_dir}")
        return {"total_files": 0, "successful": 0, "failed": 0, "errors": []}

    logger.info(f"Found {len(json_files)} extraction files to migrate")

    cursor = conn.cursor()
    successful = 0
    failed = 0
    errors = []

    # Per-file atomicity. Each file migrates inside its own SAVEPOINT so a
    # mid-file failure (e.g. a constraint error raised after paper metadata was
    # already inserted) rolls back only that file's partial rows instead of
    # riding along on the directory-wide commit below. Without this, a failed
    # file leaves half-written paper/variant rows -- and under
    # replace_existing_paper=True (delete-then-insert) it can delete a paper's
    # good rows and leave only the partial replacement. Transactions are driven
    # explicitly so SAVEPOINT semantics do not depend on sqlite3's
    # implicit-transaction mode.
    prev_isolation = conn.isolation_level
    conn.isolation_level = None
    try:
        cursor.execute("BEGIN")
        for idx, json_file in enumerate(json_files):
            savepoint = f"gvf_migrate_{idx}"
            cursor.execute(f"SAVEPOINT {savepoint}")
            success, message = migrate_extraction_file(
                cursor,
                json_file,
                replace_existing_paper=replace_existing_paper,
            )
            if success:
                cursor.execute(f"RELEASE SAVEPOINT {savepoint}")
                successful += 1
                logger.info(f"✓ [{successful}/{len(json_files)}] {message}")
            else:
                # Discard any partial writes this file made before it failed,
                # then drop the savepoint marker.
                cursor.execute(f"ROLLBACK TO SAVEPOINT {savepoint}")
                cursor.execute(f"RELEASE SAVEPOINT {savepoint}")
                failed += 1
                errors.append(message)
                logger.error(f"✗ [{failed} failures] {message}")
        # conn.commit() persists here despite isolation_level=None: the explicit
        # BEGIN above opens a transaction, so SQLite is not in autocommit at the C
        # level and commit() issues a real COMMIT. rollback() stays a method call
        # -- a safe no-op if a sqlite error already aborted the transaction, where
        # a bare cursor "ROLLBACK" would instead raise and mask the original error.
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.isolation_level = prev_isolation

    stats = {
        "total_files": len(json_files),
        "successful": successful,
        "failed": failed,
        "errors": errors,
    }

    logger.info(f"✓ Migration complete: {successful}/{len(json_files)} successful")

    return stats


# ============================================================================
# CLEANUP AND ARCHIVAL FUNCTIONS
# ============================================================================


_DEDUP_CHILD_TABLES: Tuple[str, ...] = (
    "penetrance_data",
    "individual_records",
    "fact_provenance",
    "functional_data",
    "phenotypes",
    "variant_metadata",
)


def _reparent_age_bins_for_duplicate_penetrance(cursor: sqlite3.Cursor) -> None:
    """Move unique age-bin children onto the kept duplicate parent row.

    ``dedup_existing_rows`` removes exact-duplicate ``penetrance_data`` parents.
    Before that delete runs, preserve any nonduplicate
    ``age_dependent_penetrance`` rows hanging off the parent rows that will be
    cascaded away.
    """
    parent_info = cursor.execute("PRAGMA table_info(penetrance_data)").fetchall()
    if not parent_info:
        return
    if not cursor.execute("PRAGMA table_info(age_dependent_penetrance)").fetchall():
        return
    parent_pk = next((row[1] for row in parent_info if row[5]), "penetrance_id")
    parent_value_cols = [row[1] for row in parent_info if row[1] != parent_pk]
    parent_rows = cursor.execute(
        f"SELECT {parent_pk}, {', '.join(parent_value_cols)} "  # noqa: S608
        f"FROM penetrance_data ORDER BY {parent_pk}"  # noqa: S608
    ).fetchall()

    grouped: Dict[Tuple[Any, ...], List[int]] = {}
    for row in parent_rows:
        grouped.setdefault(tuple(row[1:]), []).append(int(row[0]))

    for ids in grouped.values():
        if len(ids) < 2:
            continue
        keep_id = ids[0]
        for duplicate_id in ids[1:]:
            age_rows = cursor.execute(
                """
                SELECT age_range, penetrance_percentage, carriers_in_range,
                       affected_in_range
                FROM age_dependent_penetrance
                WHERE penetrance_id = ?
                ORDER BY age_penetrance_id
                """,
                (duplicate_id,),
            ).fetchall()
            _insert_age_dependent_penetrance(
                cursor,
                keep_id,
                [
                    {
                        "age_range": row[0],
                        "penetrance_percentage": row[1],
                        "carriers_in_range": row[2],
                        "affected_in_range": row[3],
                    }
                    for row in age_rows
                ],
            )


def _delete_orphan_age_dependent_penetrance(cursor: sqlite3.Cursor) -> int:
    """Remove age-bin rows whose penetrance parent is gone.

    This is normally redundant with the ON DELETE CASCADE constraint. It is
    intentionally kept as a belt-and-suspenders cleanup because SQLite ignores
    ``PRAGMA foreign_keys = ON`` inside an already-open transaction, and older
    ad hoc DB scripts may also have disabled FK enforcement.
    """
    if not cursor.execute("PRAGMA table_info(age_dependent_penetrance)").fetchall():
        return 0
    if not cursor.execute("PRAGMA table_info(penetrance_data)").fetchall():
        return 0
    before = cursor.execute("SELECT COUNT(*) FROM age_dependent_penetrance").fetchone()[
        0
    ]
    cursor.execute(
        """
        DELETE FROM age_dependent_penetrance
        WHERE NOT EXISTS (
            SELECT 1
            FROM penetrance_data
            WHERE penetrance_data.penetrance_id =
                  age_dependent_penetrance.penetrance_id
        )
        """
    )
    after = cursor.execute("SELECT COUNT(*) FROM age_dependent_penetrance").fetchone()[
        0
    ]
    return before - after


def dedup_existing_rows(
    conn: sqlite3.Connection,
    tables: Tuple[str, ...] = _DEDUP_CHILD_TABLES,
) -> Dict[str, int]:
    """Collapse exact-duplicate child rows already stored in a DB (idempotent).

    For each table, rows that are identical on every column except the synthetic
    primary key are reduced to one (the lowest ``rowid`` is kept).
    ``age_dependent_penetrance`` rows under removed penetrance parents are
    either dropped by FK cascade or cleaned up explicitly when FK enforcement is
    off. Genuinely different rows (e.g. real sub-cohort splits) are untouched.
    Re-running on an already-clean DB removes nothing.

    Back-fill companion to the exact-dup insert guards in ``insert_variant_data``:
    the guards stop new duplicates being written; this removes the ones that
    accumulated before the guards existed (carrier/affected counts summed N-fold
    from variants represented across multiple cohort/supplement tables).

    Returns ``{table: rows_removed}``.
    """
    removed: Dict[str, int] = {}
    conn.execute("PRAGMA foreign_keys = ON")
    cur = conn.cursor()
    for table in tables:
        info = cur.execute(f"PRAGMA table_info({table})").fetchall()  # noqa: S608
        if not info:
            continue
        pk_cols = {row[1] for row in info if row[5]}
        value_cols = [row[1] for row in info if row[1] not in pk_cols]
        if not value_cols:
            continue
        if table == "penetrance_data":
            _reparent_age_bins_for_duplicate_penetrance(cur)
        # IFNULL sentinel so NULLs group together (NULL == NULL for dedup); the
        # char(0) byte never appears in real extracted values.
        group_expr = ", ".join(f"IFNULL({c}, char(0))" for c in value_cols)
        before = cur.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]  # noqa: S608
        cur.execute(  # noqa: S608 - table/column names are code-controlled
            f"DELETE FROM {table} WHERE rowid NOT IN "
            f"(SELECT MIN(rowid) FROM {table} GROUP BY {group_expr})"
        )
        if table == "penetrance_data":
            orphan_removed = _delete_orphan_age_dependent_penetrance(cur)
            if orphan_removed:
                removed["age_dependent_penetrance_orphans"] = (
                    removed.get("age_dependent_penetrance_orphans", 0) + orphan_removed
                )
        after = cur.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]  # noqa: S608
        removed[table] = before - after
    conn.commit()
    return removed


def find_and_delete_empty_directories(
    root_dir: Path, dry_run: bool = False
) -> List[Path]:
    """
    Recursively find and delete all empty directories.

    Args:
        root_dir: Root directory to search
        dry_run: If True, only report without deleting

    Returns:
        List of deleted directories
    """
    logger.info(f"Searching for empty directories in: {root_dir}")

    deleted_dirs = []

    # Walk bottom-up to handle nested empty directories
    for dirpath, dirnames, filenames in os.walk(root_dir, topdown=False):
        current_dir = Path(dirpath)

        # Check if directory is empty (no files and no subdirectories)
        try:
            if not any(current_dir.iterdir()):
                if dry_run:
                    logger.info(
                        f"[DRY RUN] Would delete empty directory: {current_dir}"
                    )
                else:
                    current_dir.rmdir()
                    logger.info(f"✓ Deleted empty directory: {current_dir}")
                deleted_dirs.append(current_dir)
        except Exception as e:
            logger.warning(f"Could not process directory {current_dir}: {e}")

    logger.info(
        f"{'[DRY RUN] Found' if dry_run else 'Deleted'} {len(deleted_dirs)} empty directories"
    )

    return deleted_dirs


def archive_pmc_fulltext(
    pmc_dir: Path, archive_path: Optional[Path] = None, delete_after_zip: bool = False
) -> Tuple[bool, str]:
    """
    Compress pmc_fulltext directory into a ZIP archive.

    Args:
        pmc_dir: Path to pmc_fulltext directory
        archive_path: Output ZIP path (defaults to pmc_dir.zip)
        delete_after_zip: If True, delete original after successful compression

    Returns:
        Tuple of (success, message)
    """
    if not pmc_dir.exists():
        return False, f"Directory not found: {pmc_dir}"

    if not pmc_dir.is_dir():
        return False, f"Not a directory: {pmc_dir}"

    if archive_path is None:
        archive_path = pmc_dir.parent / f"{pmc_dir.name}.zip"

    logger.info(f"Archiving {pmc_dir} to {archive_path}")

    try:
        # Create ZIP archive
        with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as zipf:
            file_count = 0
            for file_path in pmc_dir.rglob("*"):
                if file_path.is_file():
                    arcname = file_path.relative_to(pmc_dir.parent)
                    zipf.write(file_path, arcname)
                    file_count += 1
                    if file_count % 100 == 0:
                        logger.info(f"  Archived {file_count} files...")

        # Verify archive was created
        if not archive_path.exists():
            return False, "Archive file was not created"

        archive_size_mb = archive_path.stat().st_size / (1024 * 1024)
        logger.info(f"✓ Archive created: {archive_path} ({archive_size_mb:.2f} MB)")

        # Optionally delete original directory
        if delete_after_zip:
            logger.info(f"Deleting original directory: {pmc_dir}")
            shutil.rmtree(pmc_dir)
            logger.info(f"✓ Deleted original directory: {pmc_dir}")

        return True, f"Successfully archived to {archive_path}"

    except Exception as e:
        error_msg = f"Failed to archive {pmc_dir}: {str(e)}"
        logger.error(error_msg)
        return False, error_msg


def cleanup_data_directory(
    data_dir: Path,
    delete_empty_dirs: bool = True,
    archive_pmc: bool = True,
    delete_pmc_after_archive: bool = False,
    dry_run: bool = False,
) -> Dict[str, Any]:
    """
    Comprehensive cleanup of data directory.

    Args:
        data_dir: Root data directory
        delete_empty_dirs: Delete empty directories (especially *_supplements)
        archive_pmc: Archive pmc_fulltext directory
        delete_pmc_after_archive: Delete pmc_fulltext after successful archival
        dry_run: If True, report actions without executing

    Returns:
        Cleanup statistics
    """
    logger.info(
        f"{'[DRY RUN] ' if dry_run else ''}Cleaning up data directory: {data_dir}"
    )

    results = {"empty_dirs_deleted": [], "archives_created": [], "errors": []}

    # Find and delete empty directories
    if delete_empty_dirs:
        try:
            deleted = find_and_delete_empty_directories(data_dir, dry_run=dry_run)
            results["empty_dirs_deleted"] = [str(d) for d in deleted]
        except Exception as e:
            error_msg = f"Error deleting empty directories: {e}"
            logger.error(error_msg)
            results["errors"].append(error_msg)

    # Archive pmc_fulltext directory
    if archive_pmc:
        pmc_dir = data_dir / "pmc_fulltext"
        if pmc_dir.exists():
            if not dry_run:
                success, message = archive_pmc_fulltext(
                    pmc_dir, delete_after_zip=delete_pmc_after_archive
                )
                if success:
                    results["archives_created"].append(
                        str(pmc_dir.parent / f"{pmc_dir.name}.zip")
                    )
                else:
                    results["errors"].append(message)
            else:
                logger.info(f"[DRY RUN] Would archive: {pmc_dir}")
                results["archives_created"].append(f"[DRY RUN] {pmc_dir}.zip")
        else:
            logger.info(f"pmc_fulltext directory not found in {data_dir}")

    return results


# ============================================================================
# MAIN CLI
# ============================================================================


def extract_gene_from_path(data_dir: Path) -> Optional[str]:
    """
    Extract gene symbol from directory path.

    Expected path structure: {OUTPUT_DIR}/{GENE}/timestamp/...
    where GENE is an uppercase symbol like BRCA1, SCN5A, etc.

    Args:
        data_dir: Path to data directory

    Returns:
        Gene symbol or None if not found
    """
    parts = data_dir.parts

    # Look for uppercase gene symbol pattern followed by timestamp
    try:
        for i, part in enumerate(parts):
            # Check if this looks like a gene symbol (uppercase, reasonable length)
            if part.isupper() and 2 <= len(part) <= 10:
                # Check if next part looks like a timestamp (YYYYMMDD_HHMMSS)
                if i + 1 < len(parts):
                    next_part = parts[i + 1]
                    if (
                        len(next_part) == 15
                        and "_" in next_part
                        and next_part.replace("_", "").isdigit()
                    ):
                        return part
                # If no timestamp but path looks reasonable, still accept it
                return part
    except Exception as e:
        logger.debug(f"Could not extract gene from path: {e}")

    return None


def extract_gene_from_json(json_file: Path) -> Optional[str]:
    """
    Extract gene symbol from a JSON extraction file.

    Args:
        json_file: Path to JSON file

    Returns:
        Gene symbol or None if not found
    """
    try:
        with open(json_file, "r", encoding="utf-8") as f:
            data = json.load(f)

        # Try to get gene from variants
        variants = data.get("variants", [])
        if variants and len(variants) > 0:
            gene = variants[0].get("gene_symbol")
            if gene:
                return gene

    except Exception as e:
        logger.debug(f"Could not extract gene from JSON: {e}")

    return None


def determine_database_name(
    data_dir: Path, extraction_dir: Optional[Path] = None
) -> str:
    """
    Determine the appropriate database name based on gene symbol.

    Args:
        data_dir: Data directory path
        extraction_dir: Directory containing JSON files (if found)

    Returns:
        Database filename (e.g., "TTR.db" or "variants.db" as fallback)
    """
    # Strategy 1: Try to extract from path
    gene = extract_gene_from_path(data_dir)

    if gene:
        logger.info(f"Detected gene '{gene}' from directory path")
        return f"{gene}.db"

    # Strategy 2: Try to extract from JSON file
    if extraction_dir:
        json_files = find_extraction_json_files(extraction_dir)
        if json_files:
            gene = extract_gene_from_json(json_files[0])
            if gene:
                logger.info(f"Detected gene '{gene}' from extraction data")
                return f"{gene}.db"

    # Fallback to generic name
    logger.info("Could not detect gene symbol, using default database name")
    return "variants.db"


def main():
    """Main CLI entrypoint."""
    parser = argparse.ArgumentParser(
        description="Migrate Gene Variant Fetcher data to SQLite and clean up file system"
    )

    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help=(
            "Path to data directory. Can point to parent dir with extractions/ "
            "subdir, or directly to dir containing *_PMID_*.json or "
            "*_extraction.json files"
        ),
    )

    parser.add_argument(
        "--db",
        type=str,
        default=None,
        help="SQLite database path (default: auto-detect based on gene symbol, e.g., TTR.db)",
    )

    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Run cleanup and archival after migration",
    )

    parser.add_argument(
        "--delete-pmc-after-archive",
        action="store_true",
        help="Delete pmc_fulltext directory after successful archival",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Dry run mode: report actions without executing",
    )

    parser.add_argument(
        "--preserve-existing-paper",
        action="store_true",
        help=(
            "Do not replace an existing paper row during migration. Use this "
            "when importing targeted recovery extractions into an existing DB."
        ),
    )

    parser.add_argument(
        "--extractions-subdir",
        type=str,
        default="extractions",
        help="Name of extractions subdirectory (default: extractions, could be extractions_rerun)",
    )

    args = parser.parse_args()

    data_dir = args.data_dir

    # ========================================================================
    # INPUT VALIDATION
    # ========================================================================
    try:
        validate_input_directory(data_dir)
    except ValidationError as e:
        logger.error(f"Input validation failed: {e}")
        return 2

    # ========================================================================
    # STEP 1: Find extraction directory
    # ========================================================================
    # Find extraction directory - try multiple strategies
    extraction_dir = None

    # Strategy 1: Check if data_dir itself contains JSON files
    json_files_in_data_dir = find_extraction_json_files(data_dir)
    if json_files_in_data_dir:
        extraction_dir = data_dir
        logger.info(
            f"✓ Found {len(json_files_in_data_dir)} JSON files directly in: {extraction_dir}"
        )
    else:
        # Strategy 2: Check specified extractions subdirectory
        extraction_subdir = data_dir / args.extractions_subdir
        if extraction_subdir.exists() and extraction_subdir.is_dir():
            json_files_in_subdir = find_extraction_json_files(extraction_subdir)
            if json_files_in_subdir:
                extraction_dir = extraction_subdir
                logger.info(
                    f"✓ Found {len(json_files_in_subdir)} JSON files in: {extraction_dir}"
                )
            else:
                logger.warning(
                    f"Subdirectory exists but no JSON files found: {extraction_subdir}"
                )
        else:
            logger.warning(f"Extraction subdirectory not found: {extraction_subdir}")

        # Strategy 3: Search for common alternative subdirectories
        if not extraction_dir:
            logger.info("Searching for alternative extraction directories...")
            alternatives = ["extractions_rerun", "extractions", "extraction"]

            # Look for timestamped subdirectories
            for subdir in data_dir.iterdir():
                if subdir.is_dir():
                    # Check timestamped subdirectories (e.g., 20251125_151454)
                    for potential_dir in [subdir] + list(subdir.iterdir()):
                        if potential_dir.is_dir():
                            json_files = find_extraction_json_files(potential_dir)
                            if json_files:
                                extraction_dir = potential_dir
                                logger.info(
                                    f"✓ Found {len(json_files)} JSON files in: {extraction_dir}"
                                )
                                break
                    if extraction_dir:
                        break

            # Try predefined alternatives
            if not extraction_dir:
                for alt in alternatives:
                    alt_dir = data_dir / alt
                    if alt_dir.exists() and alt_dir.is_dir():
                        json_files_in_alt = find_extraction_json_files(alt_dir)
                        if json_files_in_alt:
                            extraction_dir = alt_dir
                            logger.info(
                                f"✓ Found {len(json_files_in_alt)} JSON files in: {extraction_dir}"
                            )
                            break

    # Final check
    if not extraction_dir:
        logger.error(
            f"No extraction JSON files found in {data_dir} or its subdirectories"
        )
        logger.error(
            "Expected files matching pattern: *_PMID_*.json or *_extraction.json"
        )
        return 1

    # ========================================================================
    # STEP 2: Determine database name
    # ========================================================================
    if args.db:
        db_path = args.db
        logger.info(f"Using user-specified database: {db_path}")
    else:
        db_path = determine_database_name(data_dir, extraction_dir)

    # ========================================================================
    # STEP 2.5: Validate all inputs
    # ========================================================================
    try:
        validate_migrate_inputs(data_dir, extraction_dir, db_path)
    except ValidationError as e:
        logger.error(f"Input validation failed: {e}")
        return 2

    logger.info("=" * 80)
    logger.info("GENE VARIANT FETCHER: SQLite MIGRATION")
    logger.info("=" * 80)
    logger.info(f"Data directory: {data_dir}")
    logger.info(f"Database: {db_path}")
    logger.info(f"Cleanup: {args.cleanup}")
    logger.info(f"Dry run: {args.dry_run}")
    logger.info(f"Preserve existing paper evidence: {args.preserve_existing_paper}")
    logger.info("=" * 80)

    # ========================================================================
    # STEP 3: Create database schema
    # ========================================================================
    if not args.dry_run:
        conn = create_database_schema(db_path)
    else:
        logger.info("[DRY RUN] Would create database schema")
        conn = None

    # ========================================================================
    # STEP 4: Migrate extraction data
    # ========================================================================
    if not args.dry_run:
        migration_stats = migrate_extraction_directory(
            conn,
            extraction_dir,
            replace_existing_paper=not args.preserve_existing_paper,
        )

        logger.info("\n" + "=" * 80)
        logger.info("MIGRATION STATISTICS")
        logger.info("=" * 80)
        logger.info(f"Total files: {migration_stats['total_files']}")
        logger.info(f"Successful: {migration_stats['successful']}")
        logger.info(f"Failed: {migration_stats['failed']}")

        if migration_stats["errors"]:
            logger.info("\nErrors:")
            for error in migration_stats["errors"]:
                logger.info(f"  - {error}")
    else:
        logger.info(
            f"[DRY RUN] Would migrate {len(list(extraction_dir.glob('*.json')))} JSON files"
        )

    # ========================================================================
    # STEP 5: Database statistics
    # ========================================================================
    if not args.dry_run and conn:
        cursor = conn.cursor()

        cursor.execute("SELECT COUNT(*) FROM papers")
        paper_count = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM variants")
        variant_count = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM individual_records")
        individual_count = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM penetrance_data")
        penetrance_count = cursor.fetchone()[0]

        logger.info("\n" + "=" * 80)
        logger.info("DATABASE STATISTICS")
        logger.info("=" * 80)
        logger.info(f"Papers: {paper_count}")
        logger.info(f"Variants: {variant_count}")
        logger.info(f"Individual records: {individual_count}")
        logger.info(f"Penetrance data points: {penetrance_count}")
        logger.info("=" * 80)

        conn.close()

    # ========================================================================
    # STEP 6: Cleanup and archival
    # ========================================================================
    if args.cleanup:
        logger.info("\n" + "=" * 80)
        logger.info("CLEANUP AND ARCHIVAL")
        logger.info("=" * 80)

        cleanup_results = cleanup_data_directory(
            data_dir,
            delete_empty_dirs=True,
            archive_pmc=True,
            delete_pmc_after_archive=args.delete_pmc_after_archive,
            dry_run=args.dry_run,
        )

        logger.info(
            f"Empty directories {'found' if args.dry_run else 'deleted'}: {len(cleanup_results['empty_dirs_deleted'])}"
        )
        logger.info(f"Archives created: {len(cleanup_results['archives_created'])}")

        if cleanup_results["errors"]:
            logger.info(f"Errors: {len(cleanup_results['errors'])}")
            for error in cleanup_results["errors"]:
                logger.info(f"  - {error}")

    logger.info("\n" + "=" * 80)
    logger.info("✓ MIGRATION COMPLETE")
    logger.info("=" * 80)

    return 0


if __name__ == "__main__":
    import os
    import sys

    sys.exit(main())
