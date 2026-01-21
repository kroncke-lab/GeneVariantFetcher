#!/usr/bin/env python3
"""
Tests for the variant comparison tool.

Run with: pytest tests/test_compare_variants.py -v
"""

import json
import sqlite3
import tempfile
from pathlib import Path

import pandas as pd
import pytest

# Import the module under test
from compare_variants import (
    AA_3_TO_1,
    AA_1_TO_3,
    normalize_column_name,
    detect_columns,
    normalize_unicode,
    normalize_variant,
    convert_aa_3_to_1,
    convert_aa_1_to_3,
    get_variant_forms,
    compute_similarity,
    find_best_match,
    safe_int,
    compute_diff,
    aggregate_excel_data,
    aggregate_sqlite_data,
    compare_data,
    introspect_sqlite,
    extract_sqlite_data,
    ComparisonRow,
)


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def sample_excel_df():
    """Create a sample DataFrame mimicking Excel input."""
    data = {
        "PMID": ["12345678", "12345678", "87654321", "11111111"],
        "protein_change": ["p.Arg123His", "p.Val456Met", "p.Gly789Ser", "p.Leu100Pro"],
        "carriers": [10, 5, 8, None],
        "affected": [7, 3, 6, 2],
        "unaffected": [3, 2, 2, 1],
        "phenotype": ["LQT2", "LQT2", "BrS", "LQT2"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_sqlite_db():
    """Create an in-memory SQLite database with test data."""
    conn = sqlite3.connect(":memory:")
    cursor = conn.cursor()

    # Create schema matching GeneVariantFetcher
    cursor.execute("""
        CREATE TABLE papers (
            pmid TEXT PRIMARY KEY,
            title TEXT,
            gene_symbol TEXT
        )
    """)

    cursor.execute("""
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT,
            cdna_notation TEXT,
            protein_notation TEXT,
            genomic_position TEXT
        )
    """)

    cursor.execute("""
        CREATE TABLE penetrance_data (
            penetrance_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER,
            pmid TEXT,
            total_carriers_observed INTEGER,
            affected_count INTEGER,
            unaffected_count INTEGER,
            uncertain_count INTEGER
        )
    """)

    cursor.execute("""
        CREATE TABLE individual_records (
            record_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER,
            pmid TEXT,
            affected_status TEXT
        )
    """)

    # Insert test data
    cursor.execute("INSERT INTO papers VALUES ('12345678', 'Test Paper 1', 'KCNH2')")
    cursor.execute("INSERT INTO papers VALUES ('87654321', 'Test Paper 2', 'KCNH2')")

    cursor.execute("""
        INSERT INTO variants (gene_symbol, protein_notation, cdna_notation)
        VALUES ('KCNH2', 'p.Arg123His', 'c.368G>A')
    """)
    cursor.execute("""
        INSERT INTO variants (gene_symbol, protein_notation, cdna_notation)
        VALUES ('KCNH2', 'p.Val456Met', 'c.1366G>A')
    """)
    cursor.execute("""
        INSERT INTO variants (gene_symbol, protein_notation, cdna_notation)
        VALUES ('KCNH2', 'p.Gly789Ser', 'c.2365G>A')
    """)
    cursor.execute("""
        INSERT INTO variants (gene_symbol, protein_notation, cdna_notation)
        VALUES ('KCNH2', 'p.Trp999*', 'c.2997G>A')
    """)

    # Insert penetrance data with slightly different counts for testing mismatches
    cursor.execute("""
        INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed, affected_count, unaffected_count)
        VALUES (1, '12345678', 10, 7, 3)
    """)
    cursor.execute("""
        INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed, affected_count, unaffected_count)
        VALUES (2, '12345678', 6, 4, 2)
    """)  # Different from Excel!
    cursor.execute("""
        INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed, affected_count, unaffected_count)
        VALUES (3, '87654321', 8, 6, 2)
    """)
    cursor.execute("""
        INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed, affected_count, unaffected_count)
        VALUES (4, '99999999', 3, 2, 1)
    """)  # Not in Excel

    conn.commit()
    return conn


# =============================================================================
# COLUMN DETECTION TESTS
# =============================================================================


class TestColumnDetection:
    def test_normalize_column_name(self):
        assert normalize_column_name("PMID") == "pmid"
        assert normalize_column_name("PubMed_ID") == "pubmedid"
        assert normalize_column_name("protein change") == "proteinchange"
        assert normalize_column_name("  N_Affected  ") == "naffected"

    def test_detect_columns_basic(self, sample_excel_df):
        detected = detect_columns(sample_excel_df)

        assert detected["pmid"] == "PMID"
        assert detected["variant"] == "protein_change"
        assert detected["carriers_total"] == "carriers"
        assert detected["affected_count"] == "affected"
        assert detected["unaffected_count"] == "unaffected"
        assert detected["phenotype"] == "phenotype"

    def test_detect_columns_with_mapping(self, sample_excel_df):
        mapping = {"pmid": "PMID", "variant": "protein_change"}
        detected = detect_columns(sample_excel_df, mapping)

        assert detected["pmid"] == "PMID"
        assert detected["variant"] == "protein_change"

    def test_detect_columns_alternative_names(self):
        df = pd.DataFrame(
            {
                "PubMed ID": ["123"],
                "hgvs_p": ["p.A1B"],
                "n_carriers": [5],
                "cases": [3],
                "controls": [2],
            }
        )
        detected = detect_columns(df)

        assert detected["pmid"] == "PubMed ID"
        assert detected["variant"] == "hgvs_p"
        assert detected["carriers_total"] == "n_carriers"
        assert detected["affected_count"] == "cases"
        assert detected["unaffected_count"] == "controls"


# =============================================================================
# VARIANT NORMALIZATION TESTS
# =============================================================================


class TestVariantNormalization:
    def test_normalize_unicode(self):
        # Test various unicode dash characters
        assert normalize_unicode("p.Arg123\u2013His") == "p.Arg123-His"  # en-dash
        assert normalize_unicode("p.Arg123\u2014His") == "p.Arg123-His"  # em-dash
        assert normalize_unicode("p.Arg123\u2212His") == "p.Arg123-His"  # minus sign

    def test_normalize_variant_basic(self):
        assert normalize_variant("  p.Arg123His  ") == "p.Arg123His"
        assert normalize_variant("P. Arg123His") == "p.Arg123His"
        assert normalize_variant("c. 368G>A") == "c.368G>A"
        assert (
            normalize_variant("p.  Arg123His") == "p. Arg123His"
        )  # collapses multiple spaces

    def test_normalize_variant_none_handling(self):
        assert normalize_variant(None) == ""
        assert normalize_variant("") == ""

    def test_convert_aa_3_to_1(self):
        assert convert_aa_3_to_1("p.Arg123His") == "p.R123H"
        assert convert_aa_3_to_1("p.Val456Met") == "p.V456M"
        assert convert_aa_3_to_1("p.Gly789Ser") == "p.G789S"
        assert convert_aa_3_to_1("p.Trp999Ter") == "p.W999*"

    def test_convert_aa_3_to_1_with_suffix(self):
        # Frameshift
        assert convert_aa_3_to_1("p.Arg123fs") == "p.R123fs"

    def test_convert_aa_1_to_3(self):
        assert convert_aa_1_to_3("p.R123H") == "p.Arg123His"
        assert convert_aa_1_to_3("p.V456M") == "p.Val456Met"
        assert convert_aa_1_to_3("p.G789S") == "p.Gly789Ser"
        assert convert_aa_1_to_3("p.W999*") == "p.Trp999Ter"

    def test_convert_aa_non_protein(self):
        # Should return None for non-protein variants
        assert convert_aa_3_to_1("c.368G>A") is None
        assert convert_aa_1_to_3("c.368G>A") is None

    def test_get_variant_forms(self):
        forms = get_variant_forms("p.Arg123His")
        assert "p.Arg123His" in forms
        assert "p.R123H" in forms

        forms = get_variant_forms("p.R123H")
        assert "p.R123H" in forms
        assert "p.Arg123His" in forms


# =============================================================================
# MATCHING TESTS
# =============================================================================


class TestMatching:
    def test_compute_similarity_exact(self):
        score = compute_similarity("p.Arg123His", "p.Arg123His")
        assert score == 1.0

    def test_compute_similarity_different(self):
        score = compute_similarity("p.Arg123His", "p.Val456Met")
        assert score < 0.5

    def test_compute_similarity_similar(self):
        score = compute_similarity("p.Arg123His", "p.Arg123Gln")
        assert 0.7 < score < 1.0

    def test_find_best_match_exact(self):
        candidates = ["p.Arg123His", "p.Val456Met", "p.Gly789Ser"]

        match, score, match_type = find_best_match("p.Arg123His", candidates)
        assert match == "p.Arg123His"
        assert score == 1.0
        assert match_type == "exact"

    def test_find_best_match_aa_conversion(self):
        candidates = ["p.Arg123His", "p.Val456Met"]

        # 1-letter should match 3-letter
        match, score, match_type = find_best_match("p.R123H", candidates)
        assert match == "p.Arg123His"
        assert score == 1.0
        assert match_type == "exact"

    def test_find_best_match_fuzzy(self):
        candidates = ["p.Arg123His", "p.Val456Met"]

        # Similar but not exact
        match, score, match_type = find_best_match(
            "p.Arg123Gln", candidates, threshold=0.7
        )
        assert match == "p.Arg123His"
        assert match_type == "fuzzy"
        assert 0.7 <= score < 1.0

    def test_find_best_match_no_match(self):
        candidates = ["p.Arg123His", "p.Val456Met"]

        match, score, match_type = find_best_match(
            "p.Trp999Ter", candidates, threshold=0.9
        )
        assert match is None
        assert match_type == "none"


# =============================================================================
# UTILITY FUNCTION TESTS
# =============================================================================


class TestUtilities:
    def test_safe_int(self):
        assert safe_int(5) == 5
        assert safe_int(5.7) == 5
        assert safe_int("10") == 10
        assert safe_int(None) is None
        assert safe_int("invalid") is None
        assert safe_int(float("nan")) is None

    def test_compute_diff(self):
        assert compute_diff(10, 7) == 3
        assert compute_diff(5, 5) == 0
        assert compute_diff(None, 5) is None
        assert compute_diff(10, None) is None
        assert compute_diff(None, None) is None


# =============================================================================
# DATA AGGREGATION TESTS
# =============================================================================


class TestAggregation:
    def test_aggregate_excel_data(self, sample_excel_df):
        detected = {
            "pmid": "PMID",
            "variant": "protein_change",
            "carriers_total": "carriers",
            "affected_count": "affected",
            "unaffected_count": "unaffected",
            "phenotype": "phenotype",
        }

        aggregated = aggregate_excel_data(sample_excel_df, detected)

        assert len(aggregated) == 4

        key = ("12345678", "p.Arg123His")
        assert key in aggregated
        assert aggregated[key]["carriers_total"] == 10
        assert aggregated[key]["affected_count"] == 7
        assert aggregated[key]["unaffected_count"] == 3

    def test_aggregate_sqlite_data(self, sample_sqlite_db):
        # First introspect and extract
        table_info = introspect_sqlite(sample_sqlite_db)
        df = extract_sqlite_data(sample_sqlite_db, table_info)

        aggregated = aggregate_sqlite_data(df)

        assert len(aggregated) == 4

        # Check specific entry
        key = ("12345678", "p.Arg123His")
        assert key in aggregated
        assert aggregated[key]["carriers_total"] == 10
        assert aggregated[key]["affected_count"] == 7


# =============================================================================
# SQLITE INTROSPECTION TESTS
# =============================================================================


class TestSQLiteIntrospection:
    def test_introspect_sqlite(self, sample_sqlite_db):
        table_info = introspect_sqlite(sample_sqlite_db)

        assert "papers" in table_info
        assert "variants" in table_info
        assert "penetrance_data" in table_info

        assert table_info["penetrance_data"].has_pmid
        assert table_info["penetrance_data"].has_counts
        assert "affected_count" in table_info["penetrance_data"].count_columns

    def test_extract_sqlite_data(self, sample_sqlite_db):
        table_info = introspect_sqlite(sample_sqlite_db)
        df = extract_sqlite_data(sample_sqlite_db, table_info)

        assert len(df) == 4
        assert "pmid" in df.columns
        assert "variant" in df.columns
        assert "carriers_total" in df.columns
        assert "affected_count" in df.columns


# =============================================================================
# COMPARISON TESTS
# =============================================================================


class TestComparison:
    def test_compare_data_exact_match(self, sample_excel_df, sample_sqlite_db):
        # Prepare Excel data
        detected = {
            "pmid": "PMID",
            "variant": "protein_change",
            "carriers_total": "carriers",
            "affected_count": "affected",
            "unaffected_count": "unaffected",
            "phenotype": "phenotype",
        }
        excel_data = aggregate_excel_data(sample_excel_df, detected)

        # Prepare SQLite data
        table_info = introspect_sqlite(sample_sqlite_db)
        sqlite_df = extract_sqlite_data(sample_sqlite_db, table_info)
        sqlite_data = aggregate_sqlite_data(sqlite_df)

        # Compare
        results = compare_data(excel_data, sqlite_data, "exact", 0.85)

        # Should have results for all entries
        assert len(results) > 0

        # Check for expected matches and mismatches
        exact_matches = [r for r in results if r.match_type == "exact"]
        assert len(exact_matches) >= 2  # At least Arg123His and Gly789Ser

        # Check for count mismatch (Val456Met has different counts)
        mismatches = [r for r in results if r.count_mismatch]
        assert len(mismatches) >= 1

        # Check for missing in sqlite (Leu100Pro not in SQLite)
        missing_sqlite = [r for r in results if r.missing_in_sqlite]
        assert len(missing_sqlite) >= 1

        # Check for missing in excel (Trp999* in SQLite but not Excel)
        missing_excel = [r for r in results if r.missing_in_excel]
        assert len(missing_excel) >= 1

    def test_compare_data_fuzzy_match(self, sample_excel_df, sample_sqlite_db):
        # Modify Excel to have 1-letter codes
        df = sample_excel_df.copy()
        df.loc[0, "protein_change"] = "p.R123H"  # 1-letter version

        detected = {
            "pmid": "PMID",
            "variant": "protein_change",
            "carriers_total": "carriers",
            "affected_count": "affected",
            "unaffected_count": "unaffected",
            "phenotype": "phenotype",
        }
        excel_data = aggregate_excel_data(df, detected)

        # Prepare SQLite data
        table_info = introspect_sqlite(sample_sqlite_db)
        sqlite_df = extract_sqlite_data(sample_sqlite_db, table_info)
        sqlite_data = aggregate_sqlite_data(sqlite_df)

        # Compare with fuzzy matching
        results = compare_data(excel_data, sqlite_data, "fuzzy", 0.85)

        # p.R123H should match p.Arg123His
        r123h_result = [r for r in results if r.excel_variant_raw == "p.R123H"]
        assert len(r123h_result) == 1
        assert r123h_result[0].match_type == "exact"  # AA conversion counts as exact
        assert r123h_result[0].sqlite_variant_raw == "p.Arg123His"


# =============================================================================
# INTEGRATION TESTS
# =============================================================================


class TestIntegration:
    def test_full_workflow(self, tmp_path):
        """Test the complete workflow with temp files."""
        # Create Excel file
        excel_data = {
            "PMID": ["11111111", "22222222"],
            "Variant": ["p.Arg100His", "p.Val200Met"],
            "Carriers": [5, 3],
            "Affected": [4, 2],
            "Unaffected": [1, 1],
        }
        excel_df = pd.DataFrame(excel_data)
        excel_path = tmp_path / "test.xlsx"
        excel_df.to_excel(excel_path, index=False)

        # Create SQLite database
        db_path = tmp_path / "test.db"
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()

        cursor.execute("""
            CREATE TABLE variants (
                variant_id INTEGER PRIMARY KEY,
                protein_notation TEXT
            )
        """)
        cursor.execute("""
            CREATE TABLE penetrance_data (
                penetrance_id INTEGER PRIMARY KEY,
                variant_id INTEGER,
                pmid TEXT,
                total_carriers_observed INTEGER,
                affected_count INTEGER,
                unaffected_count INTEGER
            )
        """)

        cursor.execute("INSERT INTO variants VALUES (1, 'p.Arg100His')")
        cursor.execute("INSERT INTO variants VALUES (2, 'p.Val200Met')")
        cursor.execute("INSERT INTO penetrance_data VALUES (1, 1, '11111111', 5, 4, 1)")
        cursor.execute("INSERT INTO penetrance_data VALUES (2, 2, '22222222', 3, 2, 1)")
        conn.commit()
        conn.close()

        # Now run comparison (simplified - just test loading works)
        from compare_variants import load_excel_data

        df, detected = load_excel_data(excel_path, None, None)
        assert len(df) == 2
        assert detected["pmid"] == "PMID"
        assert detected["variant"] == "Variant"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
