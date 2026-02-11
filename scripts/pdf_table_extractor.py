#!/usr/bin/env python3
"""
PDF Table Extraction for GVF Pipeline

Extracts tables from PDFs and Excel files, identifies variant-containing tables,
normalizes variant notation, and outputs in GVF-compatible format.

Author: Boswell (AI Assistant)
Date: 2026-02-10
"""

import json
import logging
import re
import sys
from pathlib import Path
from typing import Any, Optional

import openpyxl
import pandas as pd
import pdfplumber
import xlrd

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ============================================================================
# Variant Pattern Definitions
# ============================================================================

# cDNA patterns: c.123A>G, c.123_125del, c.123dup, etc.
CDNA_PATTERN = re.compile(
    r"c\.[\-\*]?\d+(?:[_+\-]\d+)?(?:[ACGT]+>[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+|[ACGT]+)?",
    re.IGNORECASE,
)

# Protein patterns: p.Ala123Val, p.A123V, p.Gly628Ser, etc.
PROTEIN_PATTERN = re.compile(
    r"p\.(?:[A-Z][a-z]{2}\d+[A-Z][a-z]{2,}|[A-Z]\d+[A-Z]|[A-Z][a-z]{2}\d+(?:fs|del|dup|\*|Ter))",
    re.IGNORECASE,
)

# Genomic patterns: g.12345A>G, NC_000007.14:g.150974740T>C
GENOMIC_PATTERN = re.compile(
    r"(?:NC_\d+\.\d+[:\s]*)?(g\.\d+[ACGT]+>[ACGT]+|g\.\d+(?:del|dup|ins)[ACGT]*)",
    re.IGNORECASE,
)

# Generic variant pattern for table detection
VARIANT_KEYWORDS = [
    "variant",
    "mutation",
    "cdna",
    "c.",
    "p.",
    "amino acid",
    "protein change",
    "nucleotide",
    "substitution",
    "deletion",
    "missense",
    "nonsense",
    "frameshift",
    "splice",
    "patient",
    "proband",
    "carrier",
    "genotype",
    "allele",
]

# Column name patterns that indicate variant tables
VARIANT_COLUMN_PATTERNS = [
    r"variant",
    r"mutation",
    r"c\.",
    r"cdna",
    r"protein",
    r"p\.",
    r"amino\s*acid",
    r"nucleotide",
    r"change",
    r"patient",
    r"proband",
    r"family",
    r"phenotype",
    r"diagnosis",
    r"genotype",
    r"allele",
]


# ============================================================================
# Variant Normalization
# ============================================================================

# Three-letter to one-letter amino acid mapping
AA_THREE_TO_ONE = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
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
    "Stop": "*",
}

AA_ONE_TO_THREE = {v: k for k, v in AA_THREE_TO_ONE.items() if len(k) == 3}


def normalize_protein_notation(protein: str) -> dict:
    """
    Normalize protein notation to standard format.
    Returns dict with both one-letter and three-letter formats.
    """
    if not protein:
        return {"one_letter": None, "three_letter": None, "original": protein}

    protein = protein.strip()

    # Match three-letter format: p.Gly628Ser
    match_3 = re.match(
        r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Ter|del|dup|fs)",
        protein,
        re.IGNORECASE,
    )
    if match_3:
        aa1, pos, aa2 = match_3.groups()
        aa1 = aa1.capitalize()
        aa2 = aa2.capitalize() if aa2 not in ["*", "Ter", "del", "dup", "fs"] else aa2

        one_aa1 = AA_THREE_TO_ONE.get(aa1, aa1)
        if aa2 in ["*", "Ter"]:
            one_aa2 = "*"
            three_aa2 = "Ter"
        elif aa2 in ["del", "dup", "fs"]:
            one_aa2 = aa2
            three_aa2 = aa2
        else:
            one_aa2 = AA_THREE_TO_ONE.get(aa2, aa2)
            three_aa2 = aa2

        return {
            "one_letter": f"p.{one_aa1}{pos}{one_aa2}",
            "three_letter": f"p.{aa1}{pos}{three_aa2}",
            "original": protein,
        }

    # Match one-letter format: p.G628S
    match_1 = re.match(r"p\.([A-Z])(\d+)([A-Z]|\*|del|dup|fs)", protein, re.IGNORECASE)
    if match_1:
        aa1, pos, aa2 = match_1.groups()
        aa1 = aa1.upper()
        aa2 = aa2.upper() if aa2 not in ["*", "del", "dup", "fs"] else aa2

        three_aa1 = AA_ONE_TO_THREE.get(aa1, aa1)
        if aa2 == "*":
            three_aa2 = "Ter"
        elif aa2 in ["del", "dup", "fs"]:
            three_aa2 = aa2
        else:
            three_aa2 = AA_ONE_TO_THREE.get(aa2, aa2)

        return {
            "one_letter": f"p.{aa1}{pos}{aa2}",
            "three_letter": f"p.{three_aa1}{pos}{three_aa2}",
            "original": protein,
        }

    return {"one_letter": None, "three_letter": protein, "original": protein}


def normalize_cdna_notation(cdna: str) -> str:
    """Normalize cDNA notation to standard format."""
    if not cdna:
        return None

    cdna = cdna.strip()

    # Ensure it starts with c.
    if not cdna.lower().startswith("c."):
        cdna = "c." + cdna

    # Standardize case for nucleotides
    return cdna


def extract_variants_from_text(text: str) -> list:
    """Extract variant notations from free text."""
    variants = []

    # Find cDNA variants
    cdna_matches = CDNA_PATTERN.findall(text)
    for match in cdna_matches:
        variants.append(
            {
                "type": "cdna",
                "notation": normalize_cdna_notation(match),
                "original": match,
            }
        )

    # Find protein variants
    protein_matches = PROTEIN_PATTERN.findall(text)
    for match in protein_matches:
        normalized = normalize_protein_notation(match)
        variants.append(
            {
                "type": "protein",
                "notation": normalized["three_letter"],
                "one_letter": normalized["one_letter"],
                "original": match,
            }
        )

    return variants


# ============================================================================
# Table Detection and Classification
# ============================================================================


def is_variant_table(df: pd.DataFrame) -> tuple[bool, float]:
    """
    Determine if a DataFrame likely contains variant data.
    Returns (is_variant_table, confidence_score).
    """
    if df is None or df.empty:
        return False, 0.0

    score = 0.0

    # Check column names
    columns_str = " ".join(str(c).lower() for c in df.columns)
    for pattern in VARIANT_COLUMN_PATTERNS:
        if re.search(pattern, columns_str, re.IGNORECASE):
            score += 0.15

    # Check cell contents for variant patterns
    sample_text = df.head(10).to_string()

    cdna_count = len(CDNA_PATTERN.findall(sample_text))
    protein_count = len(PROTEIN_PATTERN.findall(sample_text))
    genomic_count = len(GENOMIC_PATTERN.findall(sample_text))

    if cdna_count > 0:
        score += min(0.3, cdna_count * 0.1)
    if protein_count > 0:
        score += min(0.3, protein_count * 0.1)
    if genomic_count > 0:
        score += min(0.3, genomic_count * 0.1)

    # Check for variant keywords
    for keyword in VARIANT_KEYWORDS:
        if keyword.lower() in sample_text.lower():
            score += 0.05

    score = min(1.0, score)
    return score >= 0.3, score


def classify_columns(df: pd.DataFrame) -> dict:
    """Classify columns in a DataFrame to identify variant-related fields."""
    classification = {
        "cdna_columns": [],
        "protein_columns": [],
        "genomic_columns": [],
        "patient_columns": [],
        "phenotype_columns": [],
        "gene_columns": [],
        "other_columns": [],
    }

    for col in df.columns:
        col_str = str(col).lower()
        sample_values = df[col].head(10).astype(str).str.cat(sep=" ")

        # Gene column
        if any(x in col_str for x in ["gene", "symbol"]):
            classification["gene_columns"].append(col)
        # cDNA column
        elif "cdna" in col_str or "c." in col_str or "nucleotide" in col_str:
            classification["cdna_columns"].append(col)
        elif len(CDNA_PATTERN.findall(sample_values)) > 2:
            classification["cdna_columns"].append(col)
        # Protein column
        elif (
            "protein" in col_str
            or "p." in col_str
            or "amino" in col_str
            or "change" in col_str
        ):
            classification["protein_columns"].append(col)
        elif len(PROTEIN_PATTERN.findall(sample_values)) > 2:
            classification["protein_columns"].append(col)
        # Genomic column
        elif "genomic" in col_str or "grch" in col_str or "g." in col_str:
            classification["genomic_columns"].append(col)
        elif len(GENOMIC_PATTERN.findall(sample_values)) > 2:
            classification["genomic_columns"].append(col)
        # Patient column
        elif any(
            x in col_str for x in ["patient", "proband", "family", "carrier", "subject"]
        ):
            classification["patient_columns"].append(col)
        # Phenotype column
        elif any(
            x in col_str for x in ["phenotype", "diagnosis", "clinical", "symptom"]
        ):
            classification["phenotype_columns"].append(col)
        else:
            classification["other_columns"].append(col)

    return classification


# ============================================================================
# PDF Table Extraction
# ============================================================================


def extract_tables_from_pdf(pdf_path: Path) -> list[dict]:
    """Extract all tables from a PDF file using pdfplumber."""
    tables = []

    try:
        with pdfplumber.open(pdf_path) as pdf:
            for page_num, page in enumerate(pdf.pages, 1):
                page_tables = page.extract_tables()

                for table_idx, table in enumerate(page_tables):
                    if table and len(table) > 1:
                        # Convert to DataFrame
                        try:
                            # First row as header
                            df = pd.DataFrame(table[1:], columns=table[0])
                            df = df.dropna(how="all")

                            if not df.empty:
                                is_variant, confidence = is_variant_table(df)
                                tables.append(
                                    {
                                        "page": page_num,
                                        "table_index": table_idx,
                                        "dataframe": df,
                                        "is_variant_table": is_variant,
                                        "confidence": confidence,
                                        "source": str(pdf_path),
                                    }
                                )
                        except Exception as e:
                            logger.warning(
                                f"Error converting table on page {page_num}: {e}"
                            )
                            continue
    except Exception as e:
        logger.error(f"Error processing PDF {pdf_path}: {e}")

    return tables


# ============================================================================
# Excel Extraction
# ============================================================================


def extract_tables_from_excel(excel_path: Path) -> list[dict]:
    """Extract all sheets from an Excel file as tables."""
    tables = []

    try:
        # Try openpyxl first (for .xlsx)
        if excel_path.suffix.lower() == ".xlsx":
            xlsx = pd.ExcelFile(excel_path, engine="openpyxl")
        else:
            xlsx = pd.ExcelFile(excel_path, engine="xlrd")

        for sheet_name in xlsx.sheet_names:
            try:
                # First read with no header to inspect structure
                df_raw = pd.read_excel(xlsx, sheet_name=sheet_name, header=None)
                df_raw = df_raw.dropna(how="all")

                if df_raw.empty:
                    continue

                # Find the best header row (look for variant keywords in first 5 rows)
                header_row = 0
                for i in range(min(5, len(df_raw))):
                    row_text = " ".join(
                        str(v).lower() for v in df_raw.iloc[i] if pd.notna(v)
                    )
                    if any(
                        kw in row_text
                        for kw in [
                            "variant",
                            "protein",
                            "cdna",
                            "change",
                            "mutation",
                            "gene",
                        ]
                    ):
                        header_row = i
                        break

                # Re-read with correct header
                df = pd.read_excel(xlsx, sheet_name=sheet_name, header=header_row)
                df = df.dropna(how="all")

                # Clean column names
                df.columns = [
                    str(c).strip() if pd.notna(c) else f"col_{i}"
                    for i, c in enumerate(df.columns)
                ]

                if not df.empty:
                    is_variant, confidence = is_variant_table(df)
                    tables.append(
                        {
                            "sheet": sheet_name,
                            "dataframe": df,
                            "is_variant_table": is_variant,
                            "confidence": confidence,
                            "source": str(excel_path),
                            "header_row": header_row,
                        }
                    )
            except Exception as e:
                logger.warning(f"Error reading sheet {sheet_name}: {e}")
                continue
    except Exception as e:
        logger.error(f"Error processing Excel {excel_path}: {e}")

    return tables


# ============================================================================
# Variant Extraction from Tables
# ============================================================================


def extract_variants_from_table(table_info: dict) -> list[dict]:
    """Extract structured variant data from a classified table."""
    df = table_info["dataframe"]
    variants = []

    col_class = classify_columns(df)

    # Get column names for extraction
    cdna_col = col_class["cdna_columns"][0] if col_class["cdna_columns"] else None
    protein_col = (
        col_class["protein_columns"][0] if col_class["protein_columns"] else None
    )
    genomic_col = (
        col_class["genomic_columns"][0] if col_class["genomic_columns"] else None
    )
    gene_col = col_class["gene_columns"][0] if col_class["gene_columns"] else None
    patient_col = (
        col_class["patient_columns"][0] if col_class["patient_columns"] else None
    )
    phenotype_col = (
        col_class["phenotype_columns"][0] if col_class["phenotype_columns"] else None
    )

    for idx, row in df.iterrows():
        variant = {
            "gene_symbol": None,
            "cdna_notation": None,
            "protein_notation": None,
            "genomic_notation": None,
            "patient_info": None,
            "phenotype": None,
            "source_table": table_info.get("page") or table_info.get("sheet"),
            "source_file": table_info["source"],
            "row_data": {k: str(v) for k, v in row.to_dict().items() if pd.notna(v)},
        }

        # Extract gene symbol
        if gene_col and pd.notna(row.get(gene_col)):
            variant["gene_symbol"] = str(row[gene_col]).strip()

        # Extract cDNA
        if cdna_col and pd.notna(row.get(cdna_col)):
            cdna_text = str(row[cdna_col])
            cdna_matches = CDNA_PATTERN.findall(cdna_text)
            if cdna_matches:
                variant["cdna_notation"] = normalize_cdna_notation(cdna_matches[0])

        # Extract protein
        if protein_col and pd.notna(row.get(protein_col)):
            protein_text = str(row[protein_col])
            protein_matches = PROTEIN_PATTERN.findall(protein_text)
            if protein_matches:
                normalized = normalize_protein_notation(protein_matches[0])
                variant["protein_notation"] = normalized["three_letter"]
                variant["protein_one_letter"] = normalized["one_letter"]

        # Extract genomic
        if genomic_col and pd.notna(row.get(genomic_col)):
            genomic_text = str(row[genomic_col])
            genomic_matches = GENOMIC_PATTERN.findall(genomic_text)
            if genomic_matches:
                variant["genomic_notation"] = genomic_matches[0]

        # Extract patient info
        if patient_col and pd.notna(row.get(patient_col)):
            variant["patient_info"] = str(row[patient_col])

        # Extract phenotype
        if phenotype_col and pd.notna(row.get(phenotype_col)):
            variant["phenotype"] = str(row[phenotype_col])

        # Also try to extract variants from other columns if not found
        if not variant["cdna_notation"] and not variant["protein_notation"]:
            row_text = " ".join(str(v) for v in row.values if pd.notna(v))
            text_variants = extract_variants_from_text(row_text)

            for tv in text_variants:
                if tv["type"] == "cdna" and not variant["cdna_notation"]:
                    variant["cdna_notation"] = tv["notation"]
                elif tv["type"] == "protein" and not variant["protein_notation"]:
                    variant["protein_notation"] = tv["notation"]
                    variant["protein_one_letter"] = tv.get("one_letter")

        # Only add if we found at least one variant notation
        if (
            variant["cdna_notation"]
            or variant["protein_notation"]
            or variant["genomic_notation"]
        ):
            variants.append(variant)

    return variants


# ============================================================================
# GVF-Compatible Output
# ============================================================================


def convert_to_gvf_format(
    variants: list[dict], source_file: str, gene_symbol: str = None
) -> dict:
    """Convert extracted variants to GVF pipeline-compatible JSON format."""
    gvf_variants = []

    for v in variants:
        # Use gene from variant if available, else use provided gene_symbol
        var_gene = v.get("gene_symbol") or gene_symbol

        gvf_variant = {
            "gene_symbol": var_gene,
            "cdna_notation": v.get("cdna_notation"),
            "protein_notation": v.get("protein_notation"),
            "genomic_notation": v.get("genomic_notation"),
            "genomic_position": None,
            "clinical_significance": None,
            "patients": {
                "count": 1 if v.get("patient_info") else None,
                "demographics": v.get("patient_info"),
                "phenotype": v.get("phenotype"),
            },
            "penetrance_data": {
                "total_carriers_observed": None,
                "affected_count": None,
                "unaffected_count": None,
                "uncertain_count": None,
                "penetrance_percentage": None,
                "age_dependent_penetrance": [],
            },
            "individual_records": [],
            "functional_data": None,
            "segregation_data": None,
            "population_frequency": None,
            "evidence_level": "table_extraction",
            "source_location": f"Table: {v.get('source_table')} from {Path(source_file).name}",
            "additional_notes": "Extracted via PDF/Excel table extraction",
            "key_quotes": [],
            "raw_row_data": v.get("row_data", {}),
        }
        gvf_variants.append(gvf_variant)

    return {
        "paper_metadata": {
            "pmid": None,
            "title": Path(source_file).stem,
            "extraction_summary": f"Table extraction found {len(gvf_variants)} variants",
        },
        "variants": gvf_variants,
        "tables_processed": [source_file],
        "extraction_metadata": {
            "total_variants_found": len(gvf_variants),
            "extraction_confidence": "medium",
            "study_type": "table_extraction",
            "challenges": [],
            "notes": "Automated extraction from PDF/Excel tables",
        },
    }


# ============================================================================
# Main Processing Functions
# ============================================================================


def process_pdf(pdf_path: Path, gene_symbol: str = None) -> dict:
    """Process a single PDF and extract variant tables."""
    logger.info(f"Processing PDF: {pdf_path}")

    tables = extract_tables_from_pdf(pdf_path)
    variant_tables = [t for t in tables if t["is_variant_table"]]

    logger.info(
        f"  Found {len(tables)} tables, {len(variant_tables)} appear to contain variants"
    )

    all_variants = []
    for table in variant_tables:
        variants = extract_variants_from_table(table)
        all_variants.extend(variants)
        logger.info(
            f"    Page {table['page']}: {len(variants)} variants (confidence: {table['confidence']:.2f})"
        )

    return {
        "source_file": str(pdf_path),
        "total_tables": len(tables),
        "variant_tables": len(variant_tables),
        "variants": all_variants,
        "gvf_output": convert_to_gvf_format(all_variants, str(pdf_path), gene_symbol),
    }


def process_excel(excel_path: Path, gene_symbol: str = None) -> dict:
    """Process a single Excel file and extract variant tables."""
    logger.info(f"Processing Excel: {excel_path}")

    tables = extract_tables_from_excel(excel_path)
    variant_tables = [t for t in tables if t["is_variant_table"]]

    logger.info(
        f"  Found {len(tables)} sheets, {len(variant_tables)} appear to contain variants"
    )

    all_variants = []
    for table in variant_tables:
        variants = extract_variants_from_table(table)
        all_variants.extend(variants)
        logger.info(
            f"    Sheet '{table['sheet']}': {len(variants)} variants (confidence: {table['confidence']:.2f})"
        )

    return {
        "source_file": str(excel_path),
        "total_sheets": len(tables),
        "variant_sheets": len(variant_tables),
        "variants": all_variants,
        "gvf_output": convert_to_gvf_format(all_variants, str(excel_path), gene_symbol),
    }


def process_directory(
    directory: Path,
    output_dir: Path = None,
    gene_symbol: str = None,
    recursive: bool = True,
) -> dict:
    """Process all PDFs and Excel files in a directory."""
    results = {
        "pdfs_processed": 0,
        "excels_processed": 0,
        "total_variants": 0,
        "files": [],
    }

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Find all files
    if recursive:
        pdfs = list(directory.rglob("*.pdf"))
        excels = list(directory.rglob("*.xlsx")) + list(directory.rglob("*.xls"))
    else:
        pdfs = list(directory.glob("*.pdf"))
        excels = list(directory.glob("*.xlsx")) + list(directory.glob("*.xls"))

    logger.info(f"Found {len(pdfs)} PDFs and {len(excels)} Excel files")

    # Process PDFs
    for pdf_path in pdfs:
        try:
            result = process_pdf(pdf_path, gene_symbol)
            results["pdfs_processed"] += 1
            results["total_variants"] += len(result["variants"])
            results["files"].append(result)

            if output_dir and result["variants"]:
                output_file = output_dir / f"{pdf_path.stem}_table_extraction.json"
                with open(output_file, "w") as f:
                    json.dump(result["gvf_output"], f, indent=2)
        except Exception as e:
            logger.error(f"Failed to process {pdf_path}: {e}")

    # Process Excel files
    for excel_path in excels:
        try:
            result = process_excel(excel_path, gene_symbol)
            results["excels_processed"] += 1
            results["total_variants"] += len(result["variants"])
            results["files"].append(result)

            if output_dir and result["variants"]:
                output_file = output_dir / f"{excel_path.stem}_table_extraction.json"
                with open(output_file, "w") as f:
                    json.dump(result["gvf_output"], f, indent=2)
        except Exception as e:
            logger.error(f"Failed to process {excel_path}: {e}")

    return results


# ============================================================================
# CLI Entry Point
# ============================================================================


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract variant tables from PDFs and Excel files"
    )
    parser.add_argument("input", help="Input PDF/Excel file or directory")
    parser.add_argument("-o", "--output", help="Output directory for JSON files")
    parser.add_argument("-g", "--gene", help="Gene symbol (e.g., KCNH2)")
    parser.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        help="Recursively process subdirectories",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    input_path = Path(args.input)
    output_dir = Path(args.output) if args.output else None

    if input_path.is_file():
        if input_path.suffix.lower() == ".pdf":
            result = process_pdf(input_path, args.gene)
        elif input_path.suffix.lower() in [".xlsx", ".xls"]:
            result = process_excel(input_path, args.gene)
        else:
            logger.error(f"Unsupported file type: {input_path.suffix}")
            sys.exit(1)

        print(json.dumps(result["gvf_output"], indent=2))

    elif input_path.is_dir():
        results = process_directory(input_path, output_dir, args.gene, args.recursive)

        print(f"\n{'=' * 60}")
        print("SUMMARY")
        print(f"{'=' * 60}")
        print(f"PDFs processed: {results['pdfs_processed']}")
        print(f"Excel files processed: {results['excels_processed']}")
        print(f"Total variants extracted: {results['total_variants']}")

        if output_dir:
            print(f"Output written to: {output_dir}")

    else:
        logger.error(f"Input not found: {input_path}")
        sys.exit(1)


if __name__ == "__main__":
    main()
