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
import json
import logging
import os
import re
import shutil
import sqlite3
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Configure logging using centralized utility
from utils.logging_utils import get_logger, setup_logging
from utils.pmid_utils import extract_gene_from_filename, extract_pmid_from_filename

setup_logging(level=logging.INFO)
logger = get_logger(__name__)

PROTEIN_NOTATION_RE = re.compile(
    r"^(?:p\.)?(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])"
    r"\d{1,4}"
    # Optional second residue for HGVS range notations like p.Asp2_Arg135del
    # or p.Lys100_Glu105delinsX. Without this, valid multi-residue dels/ins
    # are silently dropped at migration.
    r"(?:[_-](?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])\d{1,4})?"
    r"(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY*X?]|fs(?:X|\*)?\d*|del[ACDEFGHIKLMNPQRSTVWY]*|dup|ins[ACDEFGHIKLMNPQRSTVWY]*)",
    re.IGNORECASE,
)
CDNA_NOTATION_RE = re.compile(
    r"^c\.\d+(?:[+-]\d+)?[ACGT]>[ACGT]$"
    r"|^c\.\d+(?:[+-]\d+)?(?:del|dup|ins)[ACGT]*$"
    r"|^c\.\d+(?:_\d+)?(?:del|dup|ins)[ACGT]*$",
    re.IGNORECASE,
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
    """
    protein = (variant_data.get("protein_notation") or "").strip().replace(" ", "")
    cdna = (variant_data.get("cdna_notation") or "").strip().replace(" ", "")
    genomic = (variant_data.get("genomic_position") or "").strip()

    if protein and not PROTEIN_NOTATION_RE.match(protein):
        variant_data["protein_notation"] = None
        protein = ""
    if cdna and not CDNA_NOTATION_RE.match(cdna):
        variant_data["cdna_notation"] = None
        cdna = ""

    return bool(protein or cdna or genomic)


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

            FOREIGN KEY (variant_id) REFERENCES variants(variant_id) ON DELETE CASCADE,
            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
    """)

    # Index for querying individual records by affected status
    cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_individual_records_affected
        ON individual_records(variant_id, affected_status)
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

            FOREIGN KEY (variant_id) REFERENCES variants(variant_id) ON DELETE CASCADE,
            FOREIGN KEY (pmid) REFERENCES papers(pmid) ON DELETE CASCADE
        )
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

    # Try to find existing variant
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

    result = cursor.fetchone()
    if result:
        return result[0]

    # Create new variant
    cursor.execute(
        """
        INSERT INTO variants (
            gene_symbol, cdna_notation, protein_notation,
            genomic_position, clinical_significance, evidence_level
        ) VALUES (?, ?, ?, ?, ?, ?)
    """,
        (
            gene_symbol,
            cdna,
            protein,
            genomic,
            variant_data.get("clinical_significance"),
            variant_data.get("evidence_level"),
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
        gene_symbol,
        paper_meta.get("extraction_summary"),
        datetime.now().isoformat(),
    )

    if replace_existing:
        cursor.execute(
            """
            INSERT OR REPLACE INTO papers (
                pmid, title, gene_symbol, extraction_summary, extraction_timestamp
            ) VALUES (?, ?, ?, ?, ?)
        """,
            values,
        )
        return

    cursor.execute(
        """
        INSERT INTO papers (
            pmid, title, gene_symbol, extraction_summary, extraction_timestamp
        ) VALUES (?, ?, ?, ?, ?)
        ON CONFLICT(pmid) DO UPDATE SET
            title = COALESCE(excluded.title, papers.title),
            gene_symbol = COALESCE(papers.gene_symbol, excluded.gene_symbol),
            extraction_summary = COALESCE(
                excluded.extraction_summary,
                papers.extraction_summary
            ),
            extraction_timestamp = excluded.extraction_timestamp
    """,
        values,
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
            variant_id, pmid, source_location, additional_notes, key_quotes
        ) VALUES (?, ?, ?, ?, ?)
    """,
        (
            variant_id,
            pmid,
            variant_data.get("source_location"),
            variant_data.get("additional_notes"),
            json.dumps(variant_data.get("key_quotes", [])),
        ),
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
        existing_penetrance = None
        if preserve_existing_evidence:
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

        if existing_penetrance:
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
            # row. Existing parent rows are kept intact in preserve mode.
            for age_dep in penetrance.get("age_dependent_penetrance", []):
                cursor.execute(
                    """
                    INSERT INTO age_dependent_penetrance (
                        penetrance_id, age_range, penetrance_percentage,
                        carriers_in_range, affected_in_range
                    ) VALUES (?, ?, ?, ?, ?)
                """,
                    (
                        penetrance_id,
                        age_dep.get("age_range"),
                        age_dep.get("penetrance_percentage"),
                        age_dep.get("carriers_in_range"),
                        age_dep.get("affected_in_range"),
                    ),
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
        cursor.execute(
            """
            INSERT INTO individual_records (
                variant_id, pmid, individual_id, age_at_evaluation,
                age_at_onset, age_at_diagnosis, sex, affected_status,
                phenotype_details, evidence_sentence
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
            ),
        )

    # Insert functional data
    functional = variant_data.get("functional_data", {})
    if functional and (functional.get("summary") or functional.get("assays")):
        cursor.execute(
            """
            INSERT INTO functional_data (
                variant_id, pmid, summary, assays
            ) VALUES (?, ?, ?, ?)
        """,
            (
                variant_id,
                pmid,
                functional.get("summary"),
                json.dumps(functional.get("assays", [])),
            ),
        )

    # Insert phenotype data
    patients = variant_data.get("patients", {})
    if patients and (patients.get("count") or patients.get("phenotype")):
        cursor.execute(
            """
            INSERT INTO phenotypes (
                variant_id, pmid, patient_count, demographics, phenotype_description
            ) VALUES (?, ?, ?, ?, ?)
        """,
            (
                variant_id,
                pmid,
                patients.get("count"),
                patients.get("demographics"),
                patients.get("phenotype"),
            ),
        )

    # Insert variant metadata
    if variant_data.get("segregation_data") or variant_data.get("population_frequency"):
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
                    challenges, notes, extraction_timestamp, source_type, abstract_only,
                    source_file, model_used
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
                (
                    pmid,
                    extraction_meta.get("total_variants_found"),
                    extraction_meta.get("extraction_confidence"),
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

    for json_file in json_files:
        success, message = migrate_extraction_file(
            cursor,
            json_file,
            replace_existing_paper=replace_existing_paper,
        )
        if success:
            successful += 1
            logger.info(f"✓ [{successful}/{len(json_files)}] {message}")
        else:
            failed += 1
            errors.append(message)
            logger.error(f"✗ [{failed} failures] {message}")

    conn.commit()

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
