#!/usr/bin/env python3
"""
Extract variant data from TTR.db SQLite database to CSV format.

Extracts:
- Variant identity (cDNA and protein notation)
- Number of affected
- Number of unaffected
- Notes
- Individual affected data with onset information
- Onset age (when available)
"""

import argparse
import csv
import sqlite3
import sys
from pathlib import Path


def extract_variants_to_csv(
    db_path: str, output_csv: str, include_individual_records: bool = True
) -> None:
    """
    Extract variant data from TTR.db to CSV.

    Args:
        db_path: Path to TTR.db SQLite database
        output_csv: Path to output CSV file
        include_individual_records: If True, include individual-level records in separate CSV
    """
    # Check if database exists
    if not Path(db_path).exists():
        print(f"Error: Database file not found: {db_path}")
        sys.exit(1)

    print(f"Connecting to database: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row  # Enable column access by name
    cursor = conn.cursor()

    # Query to extract variant-level data with aggregated penetrance statistics
    # This aggregates data across all papers for each variant
    query = """
    SELECT
        v.variant_id,
        v.gene_symbol,
        v.cdna_notation,
        v.protein_notation,
        v.genomic_position,
        v.clinical_significance,
        v.evidence_level,

        -- Aggregated penetrance data (sum across all papers, using subqueries to avoid JOIN multiplication)
        COALESCE((
            SELECT SUM(total_carriers_observed)
            FROM penetrance_data pd2
            WHERE pd2.variant_id = v.variant_id
        ), 0) as total_carriers,
        COALESCE((
            SELECT SUM(affected_count)
            FROM penetrance_data pd3
            WHERE pd3.variant_id = v.variant_id
        ), 0) as total_affected,
        COALESCE((
            SELECT SUM(unaffected_count)
            FROM penetrance_data pd4
            WHERE pd4.variant_id = v.variant_id
        ), 0) as total_unaffected,
        COALESCE((
            SELECT SUM(uncertain_count)
            FROM penetrance_data pd5
            WHERE pd5.variant_id = v.variant_id
        ), 0) as total_uncertain,

        -- Individual records counts (from subqueries to avoid JOIN multiplication)
        COALESCE((
            SELECT COUNT(*)
            FROM individual_records ir2
            WHERE ir2.variant_id = v.variant_id
            AND ir2.affected_status = 'affected'
        ), 0) as individual_affected_count,
        COALESCE((
            SELECT COUNT(*)
            FROM individual_records ir3
            WHERE ir3.variant_id = v.variant_id
            AND ir3.affected_status = 'unaffected'
        ), 0) as individual_unaffected_count,
        COALESCE((
            SELECT COUNT(*)
            FROM individual_records ir4
            WHERE ir4.variant_id = v.variant_id
            AND ir4.affected_status = 'uncertain'
        ), 0) as individual_uncertain_count,

        -- Aggregate notes from all papers (concatenated, using subqueries for distinct values)
        -- Note: SQLite doesn't support GROUP_CONCAT(DISTINCT expr, separator), so we use REPLACE
        REPLACE((
            SELECT GROUP_CONCAT(DISTINCT additional_notes)
            FROM variant_papers vp2
            WHERE vp2.variant_id = v.variant_id
            AND additional_notes IS NOT NULL
            AND additional_notes != ''
        ), ',', ' | ') as notes,

        -- Source information
        (
            SELECT GROUP_CONCAT(DISTINCT pmid)
            FROM variant_papers vp3
            WHERE vp3.variant_id = v.variant_id
        ) as source_pmids,
        REPLACE((
            SELECT GROUP_CONCAT(DISTINCT source_location)
            FROM variant_papers vp4
            WHERE vp4.variant_id = v.variant_id
            AND source_location IS NOT NULL
            AND source_location != ''
        ), ',', ' | ') as source_locations,

        -- Phenotype information
        REPLACE((
            SELECT GROUP_CONCAT(DISTINCT phenotype_description)
            FROM phenotypes ph2
            WHERE ph2.variant_id = v.variant_id
            AND phenotype_description IS NOT NULL
            AND phenotype_description != ''
        ), ',', ' | ') as phenotype_descriptions,

        -- Count of papers reporting this variant
        (
            SELECT COUNT(DISTINCT pmid)
            FROM variant_papers vp5
            WHERE vp5.variant_id = v.variant_id
        ) as paper_count

    FROM variants v
    WHERE v.gene_symbol = 'TTR' OR v.gene_symbol IS NULL
    ORDER BY v.cdna_notation, v.protein_notation
    """

    print("Executing query...")
    cursor.execute(query)
    rows = cursor.fetchall()

    if not rows:
        print("Warning: No variants found in database.")
        conn.close()
        return

    print(f"Found {len(rows)} variants. Writing to CSV...")

    # Write main variant CSV
    with open(output_csv, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "variant_id",
            "gene_symbol",
            "cdna_notation",
            "protein_notation",
            "genomic_position",
            "clinical_significance",
            "evidence_level",
            "total_carriers",
            "total_affected",
            "total_unaffected",
            "total_uncertain",
            "individual_affected_count",
            "individual_unaffected_count",
            "individual_uncertain_count",
            "notes",
            "source_pmids",
            "source_locations",
            "phenotype_descriptions",
            "paper_count",
        ]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    "variant_id": row["variant_id"],
                    "gene_symbol": row["gene_symbol"] or "",
                    "cdna_notation": row["cdna_notation"] or "",
                    "protein_notation": row["protein_notation"] or "",
                    "genomic_position": row["genomic_position"] or "",
                    "clinical_significance": row["clinical_significance"] or "",
                    "evidence_level": row["evidence_level"] or "",
                    "total_carriers": row["total_carriers"],
                    "total_affected": row["total_affected"],
                    "total_unaffected": row["total_unaffected"],
                    "total_uncertain": row["total_uncertain"],
                    "individual_affected_count": row["individual_affected_count"],
                    "individual_unaffected_count": row["individual_unaffected_count"],
                    "individual_uncertain_count": row["individual_uncertain_count"],
                    "notes": row["notes"] or "",
                    "source_pmids": row["source_pmids"] or "",
                    "source_locations": row["source_locations"] or "",
                    "phenotype_descriptions": row["phenotype_descriptions"] or "",
                    "paper_count": row["paper_count"],
                }
            )

    print(f"✓ Main variant data written to: {output_csv}")

    # Extract individual records with onset data if requested
    if include_individual_records:
        individual_csv = output_csv.replace(".csv", "_individual_records.csv")

        individual_query = """
        SELECT
            ir.record_id,
            v.variant_id,
            v.cdna_notation,
            v.protein_notation,
            ir.individual_id,
            ir.age_at_evaluation,
            ir.age_at_onset,
            ir.age_at_diagnosis,
            ir.sex,
            ir.affected_status,
            ir.phenotype_details,
            ir.evidence_sentence,
            ir.pmid,
            p.title as paper_title
        FROM individual_records ir
        JOIN variants v ON ir.variant_id = v.variant_id
        LEFT JOIN papers p ON ir.pmid = p.pmid
        WHERE v.gene_symbol = 'TTR' OR v.gene_symbol IS NULL
        ORDER BY v.cdna_notation, v.protein_notation, ir.affected_status, ir.age_at_onset
        """

        cursor.execute(individual_query)
        individual_rows = cursor.fetchall()

        if individual_rows:
            print(f"Found {len(individual_rows)} individual records. Writing to CSV...")

            with open(individual_csv, "w", newline="", encoding="utf-8") as csvfile:
                fieldnames = [
                    "record_id",
                    "variant_id",
                    "cdna_notation",
                    "protein_notation",
                    "individual_id",
                    "age_at_evaluation",
                    "age_at_onset",
                    "age_at_diagnosis",
                    "sex",
                    "affected_status",
                    "phenotype_details",
                    "evidence_sentence",
                    "pmid",
                    "paper_title",
                ]

                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()

                for row in individual_rows:
                    writer.writerow(
                        {
                            "record_id": row["record_id"],
                            "variant_id": row["variant_id"],
                            "cdna_notation": row["cdna_notation"] or "",
                            "protein_notation": row["protein_notation"] or "",
                            "individual_id": row["individual_id"] or "",
                            "age_at_evaluation": row["age_at_evaluation"] or "",
                            "age_at_onset": row["age_at_onset"] or "",
                            "age_at_diagnosis": row["age_at_diagnosis"] or "",
                            "sex": row["sex"] or "",
                            "affected_status": row["affected_status"] or "",
                            "phenotype_details": row["phenotype_details"] or "",
                            "evidence_sentence": row["evidence_sentence"] or "",
                            "pmid": row["pmid"] or "",
                            "paper_title": row["paper_title"] or "",
                        }
                    )

            print(f"✓ Individual records written to: {individual_csv}")
        else:
            print("No individual records found in database.")

    conn.close()
    print("\n✓ Extraction complete!")


def main():
    parser = argparse.ArgumentParser(
        description="Extract variant data from TTR.db to CSV format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract to default output file
  python extract_ttr_to_csv.py --db "/path/to/TTR.db"

  # Specify custom output file
  python extract_ttr_to_csv.py --db "/path/to/TTR.db" --output variants.csv

  # Skip individual records
  python extract_ttr_to_csv.py --db "/path/to/TTR.db" --no-individual-records
        """,
    )

    parser.add_argument(
        "--db", type=str, required=True, help="Path to TTR.db SQLite database file"
    )

    parser.add_argument(
        "--output",
        type=str,
        default="ttr_variants.csv",
        help="Output CSV file path (default: ttr_variants.csv)",
    )

    parser.add_argument(
        "--no-individual-records",
        action="store_true",
        help="Skip extraction of individual records CSV",
    )

    args = parser.parse_args()

    extract_variants_to_csv(
        db_path=args.db,
        output_csv=args.output,
        include_individual_records=not args.no_individual_records,
    )


if __name__ == "__main__":
    main()
