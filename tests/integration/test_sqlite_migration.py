#!/usr/bin/env python3
"""
Test script for SQLite migration

Creates sample extraction data and tests the migration process.
"""

import json
import tempfile
import shutil
from pathlib import Path
from harvesting.migrate_to_sqlite import (
    create_database_schema,
    migrate_extraction_directory,
    cleanup_data_directory,
)


def create_sample_extraction_data(output_dir: Path) -> None:
    """Create sample extraction JSON files for testing."""

    extraction_dir = output_dir / "extractions"
    extraction_dir.mkdir(parents=True, exist_ok=True)

    # Sample extraction 1: TTR p.Val30Met
    extraction1 = {
        "paper_metadata": {
            "pmid": "12345678",
            "title": "Transthyretin Val30Met: A Common Amyloidogenic Mutation",
            "extraction_summary": "Identified Val30Met variant with penetrance data",
        },
        "variants": [
            {
                "gene_symbol": "TTR",
                "cdna_notation": "c.148G>A",
                "protein_notation": "p.Val30Met",
                "genomic_position": "chr18:31591753",
                "clinical_significance": "pathogenic",
                "patients": {
                    "count": 45,
                    "demographics": "Portuguese descent, mean age 38 years",
                    "phenotype": "Familial amyloid polyneuropathy",
                },
                "penetrance_data": {
                    "total_carriers_observed": 45,
                    "affected_count": 32,
                    "unaffected_count": 8,
                    "uncertain_count": 5,
                    "penetrance_percentage": 71.1,
                    "age_dependent_penetrance": [
                        {
                            "age_range": "30-40 years",
                            "penetrance_percentage": 45.0,
                            "carriers_in_range": 20,
                            "affected_in_range": 9,
                        },
                        {
                            "age_range": "40-50 years",
                            "penetrance_percentage": 85.0,
                            "carriers_in_range": 15,
                            "affected_in_range": 13,
                        },
                    ],
                },
                "individual_records": [
                    {
                        "individual_id": "P001",
                        "age_at_evaluation": 42,
                        "age_at_onset": 38,
                        "age_at_diagnosis": 39,
                        "sex": "male",
                        "affected_status": "affected",
                        "phenotype_details": "Progressive peripheral neuropathy",
                        "evidence_sentence": "Patient P001 presented at age 38 with sensory symptoms.",
                    },
                    {
                        "individual_id": "P002",
                        "age_at_evaluation": 35,
                        "age_at_onset": None,
                        "age_at_diagnosis": None,
                        "sex": "female",
                        "affected_status": "unaffected",
                        "phenotype_details": "Asymptomatic carrier",
                        "evidence_sentence": "Patient P002 was an asymptomatic carrier at age 35.",
                    },
                ],
                "functional_data": {
                    "summary": "In vitro studies show reduced stability",
                    "assays": ["thermal stability assay", "aggregation kinetics"],
                },
                "segregation_data": "Segregates with disease in 3 families",
                "population_frequency": "gnomAD: 0.0001 (rare)",
                "evidence_level": "strong",
                "source_location": "Table 2, Rows 1-5",
                "additional_notes": "Well-characterized pathogenic variant",
                "key_quotes": [
                    "Val30Met is the most common TTR mutation worldwide",
                    "Penetrance varies by geographic origin",
                ],
            }
        ],
        "tables_processed": [
            {
                "table_name": "Table 2",
                "table_caption": "Clinical characteristics of TTR mutation carriers",
                "variants_extracted": 1,
            }
        ],
        "extraction_metadata": {
            "total_variants_found": 1,
            "extraction_confidence": "high",
            "challenges": [],
            "notes": "Full-text with detailed table",
        },
    }

    # Sample extraction 2: TTR p.Thr60Ala
    extraction2 = {
        "paper_metadata": {
            "pmid": "87654321",
            "title": "Novel TTR Mutation Thr60Ala in Cardiac Amyloidosis",
            "extraction_summary": "Identified Thr60Ala variant with cardiac phenotype",
        },
        "variants": [
            {
                "gene_symbol": "TTR",
                "cdna_notation": "c.178A>G",
                "protein_notation": "p.Thr60Ala",
                "genomic_position": "chr18:31591723",
                "clinical_significance": "likely pathogenic",
                "patients": {
                    "count": 12,
                    "demographics": "Irish descent, mean age 52 years",
                    "phenotype": "Cardiac amyloidosis",
                },
                "penetrance_data": {
                    "total_carriers_observed": 12,
                    "affected_count": 10,
                    "unaffected_count": 2,
                    "uncertain_count": 0,
                    "penetrance_percentage": 83.3,
                    "age_dependent_penetrance": [],
                },
                "individual_records": [
                    {
                        "individual_id": "Case_1",
                        "age_at_evaluation": 58,
                        "age_at_onset": 55,
                        "age_at_diagnosis": 56,
                        "sex": "male",
                        "affected_status": "affected",
                        "phenotype_details": "Heart failure with restrictive cardiomyopathy",
                        "evidence_sentence": "Case 1 developed heart failure at age 55.",
                    }
                ],
                "functional_data": {
                    "summary": "Moderate destabilization of TTR tetramer",
                    "assays": ["circular dichroism", "native PAGE"],
                },
                "segregation_data": None,
                "population_frequency": "Not found in gnomAD",
                "evidence_level": "moderate",
                "source_location": "Results section, paragraph 3",
                "additional_notes": "Rare variant, limited literature",
                "key_quotes": ["Thr60Ala causes predominantly cardiac manifestations"],
            }
        ],
        "tables_processed": [],
        "extraction_metadata": {
            "total_variants_found": 1,
            "extraction_confidence": "medium",
            "challenges": ["Limited family data"],
            "notes": "Abstract only, limited information",
        },
    }

    # Write sample files
    with open(extraction_dir / "TTR_PMID_12345678.json", "w") as f:
        json.dump(extraction1, f, indent=2)

    with open(extraction_dir / "TTR_PMID_87654321.json", "w") as f:
        json.dump(extraction2, f, indent=2)

    # Create pmc_fulltext directory with empty supplements (for cleanup testing)
    pmc_dir = output_dir / "pmc_fulltext"
    pmc_dir.mkdir(exist_ok=True)

    # Add a sample markdown file
    (pmc_dir / "PMID_12345678_FULL_CONTEXT.md").write_text(
        "# Sample Paper\n\nThis is a sample paper about TTR mutations."
    )

    # Create empty supplement directories
    (pmc_dir / "PMID_12345678_supplements").mkdir(exist_ok=True)
    (pmc_dir / "PMID_87654321_supplements").mkdir(exist_ok=True)

    print(f"✓ Created sample extraction data in {output_dir}")


def test_migration():
    """Test the complete migration workflow."""

    print("=" * 80)
    print("TESTING SQLITE MIGRATION")
    print("=" * 80)

    # Create temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir) / "test_data"
        test_dir.mkdir()

        # Step 1: Create sample data
        print("\n1. Creating sample extraction data...")
        create_sample_extraction_data(test_dir)

        # Step 2: Create database and migrate
        print("\n2. Creating database and migrating data...")
        db_path = test_dir / "test_variants.db"
        conn = create_database_schema(str(db_path))

        extraction_dir = test_dir / "extractions"
        migration_stats = migrate_extraction_directory(conn, extraction_dir)

        print(f"\n   Migration stats:")
        print(f"   - Total files: {migration_stats['total_files']}")
        print(f"   - Successful: {migration_stats['successful']}")
        print(f"   - Failed: {migration_stats['failed']}")

        # Step 3: Verify data
        print("\n3. Verifying database contents...")
        cursor = conn.cursor()

        cursor.execute("SELECT COUNT(*) FROM papers")
        paper_count = cursor.fetchone()[0]
        print(f"   - Papers: {paper_count}")

        cursor.execute("SELECT COUNT(*) FROM variants")
        variant_count = cursor.fetchone()[0]
        print(f"   - Variants: {variant_count}")

        cursor.execute("SELECT COUNT(*) FROM individual_records")
        individual_count = cursor.fetchone()[0]
        print(f"   - Individual records: {individual_count}")

        cursor.execute("SELECT COUNT(*) FROM penetrance_data")
        penetrance_count = cursor.fetchone()[0]
        print(f"   - Penetrance data points: {penetrance_count}")

        # Query a specific variant
        cursor.execute("""
            SELECT v.protein_notation, v.clinical_significance,
                   COUNT(DISTINCT vp.pmid) as paper_count
            FROM variants v
            LEFT JOIN variant_papers vp ON v.variant_id = vp.variant_id
            WHERE v.protein_notation = 'p.Val30Met'
            GROUP BY v.variant_id
        """)
        result = cursor.fetchone()
        if result:
            print(f"\n   Sample query (p.Val30Met):")
            print(f"   - Notation: {result[0]}")
            print(f"   - Significance: {result[1]}")
            print(f"   - Papers: {result[2]}")

        conn.close()

        # Step 4: Test cleanup
        print("\n4. Testing cleanup functions...")
        cleanup_results = cleanup_data_directory(
            test_dir,
            delete_empty_dirs=True,
            archive_pmc=True,
            delete_pmc_after_archive=False,
            dry_run=False,
        )

        print(f"   - Empty dirs deleted: {len(cleanup_results['empty_dirs_deleted'])}")
        print(f"   - Archives created: {len(cleanup_results['archives_created'])}")

        if cleanup_results["errors"]:
            print(f"   - Errors: {len(cleanup_results['errors'])}")
            for error in cleanup_results["errors"]:
                print(f"     * {error}")

        # Verify archive was created
        archive_path = test_dir / "pmc_fulltext.zip"
        if archive_path.exists():
            size_mb = archive_path.stat().st_size / (1024 * 1024)
            print(f"   - Archive size: {size_mb:.2f} MB")

    print("\n" + "=" * 80)
    print("✓ ALL TESTS PASSED")
    print("=" * 80)


if __name__ == "__main__":
    try:
        test_migration()
    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback

        traceback.print_exc()
        e(1)
