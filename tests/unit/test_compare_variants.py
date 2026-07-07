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
from cli.compare_variants import (
    AA_3_TO_1,
    AA_1_TO_3,
    normalize_column_name,
    detect_columns,
    normalize_unicode,
    normalize_variant,
    convert_aa_3_to_1,
    convert_aa_1_to_3,
    get_variant_forms,
    to_canonical_form,
    compute_similarity,
    find_best_match,
    safe_int,
    compute_diff,
    load_excel_data,
    aggregate_excel_data,
    aggregate_sqlite_data,
    compare_data,
    introspect_sqlite,
    extract_sqlite_data,
    compute_recall_summary,
    compute_precision_summary,
    generate_outputs,
    ComparisonRow,
)


def _make_row(
    *,
    pmid: str,
    matched: bool,
    missing_in_excel: bool = False,
    sqlite_source_layer: str | None = None,
    sqlite_carriers_total: int | None = None,
) -> ComparisonRow:
    """Build a minimal ComparisonRow for precision testing.

    matched=True  -> a gold row that matched a DB row (exact).
    matched=False & missing_in_excel=False -> a gold row with no DB match.
    missing_in_excel=True -> a DB-only (extra) row.
    """
    if missing_in_excel:
        return ComparisonRow(
            pmid=pmid,
            excel_variant_raw="",
            excel_variant_norm="",
            sqlite_variant_raw="p.Xxx999Yyy",
            sqlite_variant_norm="p.Xxx999Yyy",
            match_type="none",
            match_score=None,
            excel_carriers_total=None,
            sqlite_carriers_total=sqlite_carriers_total,
            carriers_diff=None,
            excel_affected=None,
            sqlite_affected=None,
            affected_diff=None,
            excel_unaffected=None,
            sqlite_unaffected=None,
            unaffected_diff=None,
            missing_in_excel=True,
            sqlite_source_layer=sqlite_source_layer,
        )
    return ComparisonRow(
        pmid=pmid,
        excel_variant_raw="p.Arg1His",
        excel_variant_norm="p.Arg1His",
        sqlite_variant_raw="p.Arg1His" if matched else None,
        sqlite_variant_norm="p.Arg1His" if matched else None,
        match_type="exact" if matched else "none",
        match_score=1.0 if matched else None,
        excel_carriers_total=None,
        sqlite_carriers_total=None,
        carriers_diff=None,
        excel_affected=None,
        sqlite_affected=None,
        affected_diff=None,
        excel_unaffected=None,
        sqlite_unaffected=None,
        unaffected_diff=None,
        missing_in_sqlite=not matched,
        sqlite_source_layer=sqlite_source_layer,
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

    def test_detect_columns_does_not_invent_unaffected_from_affected(self):
        df = pd.DataFrame(
            {
                "variant": ["G148R"],
                "pmid": ["34135346"],
                "carriers": [1],
                "affected": [0],
            }
        )
        detected = detect_columns(df)

        assert detected["affected_count"] == "affected"
        assert detected["unaffected_count"] is None
        assert detected["rsid"] is None

    def test_load_csv_recall_input(self, tmp_path):
        csv_path = tmp_path / "KCNQ1_recall_input.csv"
        csv_path.write_text(
            "variant,pmid,carriers,affected\nG148R,34135346,1,0\nD202N,34135346,2,2\n",
            encoding="utf-8",
        )

        df, detected = load_excel_data(csv_path, None, None)
        aggregated = aggregate_excel_data(df, detected)

        assert detected["pmid"] == "pmid"
        assert detected["variant"] == "variant"
        assert detected["carriers_total"] == "carriers"
        assert detected["affected_count"] == "affected"
        assert ("34135346", "G148R") in aggregated
        assert aggregated[("34135346", "D202N")]["carriers_total"] == 2


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
        assert normalize_variant("p.  Arg123His") == "p.Arg123His"  # collapses spaces

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

    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ("p.Pro926AlafsTer14", "P926fsX"),
            ("p.Pro926AlaX14", "P926fsX"),
            ("p.V131fs/185", "V131fsX"),
            ("R513VfsTer8", "R513fsX"),
            ("p.K897fs+*49", "K897fsX"),
            ("I82dup", "I82dup"),
            ("p.K1505_Q1507DEL", "K1505_Q1507Del"),
            ("K1505_Q1507del", "K1505_Q1507Del"),
            ("p.A178-G189del", "A178_G189Del"),
            ("Q521-Y522DELT", "Q521_Y522DelT"),
            ("p.Gln1507_Pro1509del", "Q1507_P1509Del"),
            ("4944_4945INSH", "4944_4945InsH"),
            ("73-73 DEL AAP", "P73Del"),
            ("339DELF", "F339Del"),
            ("392INSW", "W392Ins"),
            ("614DELH", "H614Del"),
        ],
    )
    def test_to_canonical_form_frameshift_variants_seen_in_extractions(
        self, raw, expected
    ):
        assert to_canonical_form(raw) == expected

    def test_range_indel_forms_match_case_and_prefix_variants(self):
        forms = get_variant_forms("p.K1505_Q1507DEL")

        assert "K1505_Q1507Del" in forms
        assert "p.K1505_Q1507Del" in forms
        assert "p.Lys1505_Gln1507DEL" in forms

        match, score, match_type = find_best_match(
            "p.K1505_Q1507DEL", ["K1505_Q1507del"]
        )
        assert match == "K1505_Q1507del"
        assert score == 1.0
        assert match_type == "exact"

        match, score, match_type = find_best_match("p.A178-G189del", ["A178-G189DEL"])
        assert match == "A178-G189DEL"
        assert score == 1.0
        assert match_type == "exact"


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

        # Keys now use canonical form: p.Arg123His -> R123H
        key = ("12345678", "R123H")
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

        # Check specific entry - keys now use canonical form
        key = ("12345678", "R123H")
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

    def test_extract_sqlite_data_prefers_explicit_source_layer(self):
        conn = sqlite3.connect(":memory:")
        cur = conn.cursor()
        cur.executescript(
            """
            CREATE TABLE variants (
                variant_id INTEGER PRIMARY KEY,
                gene_symbol TEXT,
                cdna_notation TEXT,
                protein_notation TEXT,
                genomic_position TEXT
            );
            CREATE TABLE variant_papers (
                variant_id INTEGER,
                pmid TEXT,
                source_location TEXT,
                source_layer TEXT
            );
            INSERT INTO variants VALUES (1, 'KCNH2', NULL, 'p.Arg123His', NULL);
            INSERT INTO variant_papers
            VALUES (1, '12345678', 'ClinVar (PMID citation)', 'figure');
            """
        )

        df = extract_sqlite_data(conn, introspect_sqlite(conn))

        assert len(df) == 1
        assert df.iloc[0]["source_layer"] == "figure"
        conn.close()

    def test_extract_sqlite_data_rejects_junk_in_figure_and_regex_table_layers(self):
        conn = sqlite3.connect(":memory:")
        cur = conn.cursor()
        cur.executescript(
            """
            CREATE TABLE variants (
                variant_id INTEGER PRIMARY KEY,
                gene_symbol TEXT,
                cdna_notation TEXT,
                protein_notation TEXT,
                genomic_position TEXT
            );
            CREATE TABLE variant_papers (
                variant_id INTEGER,
                pmid TEXT,
                source_location TEXT,
                source_layer TEXT
            );
            INSERT INTO variants VALUES (1, 'SCN5A', NULL, 'SCN5A', NULL);
            INSERT INTO variant_papers VALUES (1, '111', 'table regex', 'regex_table');
            INSERT INTO variants VALUES (2, 'KCNH2', NULL, 'RQ', NULL);
            INSERT INTO variant_papers VALUES (2, '222', 'figure-reader', 'figure');
            INSERT INTO variants VALUES (3, 'KCNQ1', NULL, 'D4S6 residue 56', NULL);
            INSERT INTO variant_papers VALUES (3, '333', 'figure-reader', 'figure');
            INSERT INTO variants VALUES (4, 'KCNH2', NULL, 'RQ', NULL);
            INSERT INTO variant_papers VALUES (4, '444', 'scanner text scan', 'regex_text');
            INSERT INTO variants VALUES (5, 'KCNH2', NULL, 'p.Arg123His', NULL);
            INSERT INTO variant_papers VALUES (5, '555', 'Supplement Table S4', 'regex_table');
            """
        )

        df = extract_sqlite_data(conn, introspect_sqlite(conn))

        assert set(df["pmid"]) == {"444", "555"}
        assert "regex_text" in set(df["source_layer"])
        assert "regex_table" in set(df["source_layer"])
        conn.close()


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

    def test_recall_summary_counts_gold_dimensions(
        self, sample_excel_df, sample_sqlite_db
    ):
        detected = {
            "pmid": "PMID",
            "variant": "protein_change",
            "carriers_total": "carriers",
            "affected_count": "affected",
            "unaffected_count": "unaffected",
            "phenotype": "phenotype",
        }
        excel_data = aggregate_excel_data(sample_excel_df, detected)
        table_info = introspect_sqlite(sample_sqlite_db)
        sqlite_df = extract_sqlite_data(sample_sqlite_db, table_info)
        sqlite_data = aggregate_sqlite_data(sqlite_df)

        results = compare_data(excel_data, sqlite_data, "exact", 0.85)
        recall = compute_recall_summary(results)

        assert recall["pmids"] == {"matched": 2, "gold": 3, "recall": 2 / 3}
        assert recall["variant_rows"] == {"matched": 3, "gold": 4, "recall": 0.75}
        assert recall["unique_variants"] == {"matched": 3, "gold": 4, "recall": 0.75}
        assert recall["patients"] == {"matched": 23, "gold": 26, "recall": 23 / 26}
        assert recall["affected"] == {"matched": 16, "gold": 18, "recall": 16 / 18}
        assert recall["unaffected"] == {"matched": 7, "gold": 8, "recall": 7 / 8}

    def test_generate_outputs_writes_recall_metrics(
        self, sample_excel_df, sample_sqlite_db, tmp_path
    ):
        detected = {
            "pmid": "PMID",
            "variant": "protein_change",
            "carriers_total": "carriers",
            "affected_count": "affected",
            "unaffected_count": "unaffected",
            "phenotype": "phenotype",
        }
        excel_data = aggregate_excel_data(sample_excel_df, detected)
        table_info = introspect_sqlite(sample_sqlite_db)
        sqlite_df = extract_sqlite_data(sample_sqlite_db, table_info)
        sqlite_data = aggregate_sqlite_data(sqlite_df)
        results = compare_data(excel_data, sqlite_data, "exact", 0.85)

        summary = generate_outputs(
            results,
            tmp_path,
            Path("gold.xlsx"),
            Path("gvf.db"),
        )

        assert summary["recall"]["unique_variants"]["matched"] == 3
        assert summary["recall"]["patients"]["gold"] == 26
        report = (tmp_path / "report.md").read_text(encoding="utf-8")
        assert "## Recall" in report
        assert "| Affected | 16/18 (88.9%) |" in report

    def test_precision_summary_restricts_extra_rows_to_gold_pmids(self):
        """Extra DB rows count only when their PMID was curated by gold.

        3 matched gold rows + 1 extra DB row on a gold PMID gives
        3 / (3 + 1) = 0.75. Two further extra DB rows on a NON-gold PMID must
        be excluded from the denominator and must not change the ratio.
        """
        results = [
            _make_row(pmid="111", matched=True),
            _make_row(pmid="111", matched=True),
            _make_row(pmid="222", matched=True),
            # Extra DB-only row on a gold PMID -> counts in denominator.
            _make_row(pmid="222", matched=False, missing_in_excel=True),
            # Extra DB-only rows on a PMID gold never curated -> excluded.
            _make_row(pmid="999", matched=False, missing_in_excel=True),
            _make_row(pmid="999", matched=False, missing_in_excel=True),
        ]

        precision = compute_precision_summary(results)

        assert precision["matched_db"] == 3
        assert precision["extra_on_gold_pmids"] == 1
        assert precision["counted_extra_on_gold_pmids"] == 0
        assert precision["precision_vs_gold_pmids"] == 0.75
        assert precision["precision_vs_counted_gold_pmids"] == 1.0
        assert "note" in precision and "NOT clean precision" in precision["note"]

    def test_precision_summary_decomposes_counted_extras_by_source_layer(self):
        """Zero-count linker rows and counted extraction rows stay separable."""
        results = [
            _make_row(
                pmid="111",
                matched=True,
                sqlite_source_layer="llm_table",
            ),
            _make_row(
                pmid="111",
                matched=True,
                sqlite_source_layer="clinvar",
            ),
            # Zero-count ClinVar linker row: included in the upper-bound proxy,
            # excluded from the counted-extra denominator.
            _make_row(
                pmid="111",
                matched=False,
                missing_in_excel=True,
                sqlite_source_layer="clinvar",
            ),
            # Count-bearing table extraction row: included in both denominators.
            _make_row(
                pmid="111",
                matched=False,
                missing_in_excel=True,
                sqlite_source_layer="llm_table",
                sqlite_carriers_total=4,
            ),
        ]

        precision = compute_precision_summary(results)

        assert precision["matched_db"] == 2
        assert precision["extra_on_gold_pmids"] == 2
        assert precision["counted_extra_on_gold_pmids"] == 1
        assert precision["precision_vs_gold_pmids"] == 0.5
        assert precision["precision_vs_counted_gold_pmids"] == pytest.approx(2 / 3)

        by_layer = precision["by_source_layer"]
        assert by_layer["clinvar"]["extra_on_gold_pmids"] == 1
        assert by_layer["clinvar"]["counted_extra_on_gold_pmids"] == 0
        assert by_layer["clinvar"]["precision_vs_counted_gold_pmids"] == 1.0
        assert by_layer["llm_table"]["extra_on_gold_pmids"] == 1
        assert by_layer["llm_table"]["counted_extra_on_gold_pmids"] == 1
        assert by_layer["llm_table"]["precision_vs_counted_gold_pmids"] == 0.5

    def test_precision_summary_non_gold_extras_do_not_change_ratio(self):
        """Adding only non-gold-PMID extra rows leaves precision at 1.0."""
        base = [
            _make_row(pmid="111", matched=True),
            _make_row(pmid="222", matched=True),
        ]
        with_non_gold_extras = base + [
            _make_row(pmid="888", matched=False, missing_in_excel=True),
            _make_row(pmid="999", matched=False, missing_in_excel=True),
        ]

        baseline = compute_precision_summary(base)
        widened = compute_precision_summary(with_non_gold_extras)

        assert baseline["precision_vs_gold_pmids"] == 1.0
        assert widened["extra_on_gold_pmids"] == 0
        assert widened["precision_vs_gold_pmids"] == 1.0

    def test_precision_summary_none_when_no_judgeable_rows(self):
        """No matched rows and no extras on gold PMIDs -> ratio is None."""
        results = [
            # Gold row that was missed entirely (no DB match) -> not in denom.
            _make_row(pmid="111", matched=False),
            # Extra DB row on a non-gold PMID -> excluded.
            _make_row(pmid="999", matched=False, missing_in_excel=True),
        ]

        precision = compute_precision_summary(results)

        assert precision["matched_db"] == 0
        assert precision["extra_on_gold_pmids"] == 0
        assert precision["counted_extra_on_gold_pmids"] == 0
        assert precision["precision_vs_gold_pmids"] is None
        assert precision["precision_vs_counted_gold_pmids"] is None

    def test_generate_outputs_emits_precision_and_gold_pmid_csv(self, tmp_path):
        """generate_outputs surfaces precision + the gold-PMID extras CSV."""
        results = [
            _make_row(pmid="111", matched=True),
            _make_row(pmid="111", matched=True),
            _make_row(pmid="222", matched=True),
            _make_row(pmid="222", matched=False, missing_in_excel=True),
            _make_row(pmid="999", matched=False, missing_in_excel=True),
        ]

        summary = generate_outputs(results, tmp_path, Path("gold.xlsx"), Path("gvf.db"))

        assert summary["precision"]["matched_db"] == 3
        assert summary["precision"]["extra_on_gold_pmids"] == 1
        assert summary["precision"]["counted_extra_on_gold_pmids"] == 0
        assert summary["precision"]["precision_vs_gold_pmids"] == 0.75
        assert summary["precision"]["precision_vs_counted_gold_pmids"] == 1.0

        # Recall/MAE blocks unchanged in shape by the precision addition.
        assert "recall" in summary and "mae" in summary

        gold_pmid_csv = tmp_path / "unmatched_db_rows_on_gold_pmids.csv"
        assert gold_pmid_csv.exists()
        rows = pd.read_csv(gold_pmid_csv, dtype=str)
        # Only the gold-PMID extra (222), not the non-gold extra (999).
        assert set(rows["pmid"]) == {"222"}

        report = (tmp_path / "report.md").read_text(encoding="utf-8")
        assert "## Precision (counted extras vs gold PMIDs)" in report


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
        from cli.compare_variants import load_excel_data

        df, detected = load_excel_data(excel_path, None, None)
        assert len(df) == 2
        assert detected["pmid"] == "PMID"
        assert detected["variant"] == "Variant"


# =============================================================================
# REGRESSION TESTS FOR v4 RECALL FIXES
# (positional-digit guard, greedy 1-to-1, cDNA->protein bridge)
# =============================================================================


class TestPositionalDigitGuard:
    """Fuzzy matcher must reject candidates at a different protein position."""

    def test_rejects_different_position_missense(self):
        # G572S vs G628S: Levenshtein on canonical strings reports ~0.82
        # similarity. The guard must reject — different codon = different
        # variant, regardless of string similarity.
        match, score, match_type = find_best_match(
            "p.Gly572Ser", ["p.Gly628Ser"], threshold=0.7
        )
        assert match is None
        assert match_type == "none"

    def test_rejects_different_position_frameshift(self):
        # T152fsX vs T152H: same position but different variant type.
        # Position guard does not reject this; threshold + similarity decide.
        # The bigger concern is e.g. R582C vs R620C — different positions.
        match, _, match_type = find_best_match(
            "p.Arg582Cys", ["p.Arg620Cys"], threshold=0.7
        )
        assert match is None
        assert match_type == "none"

    def test_accepts_same_position_different_alt(self):
        # D501N vs D501R: same position, different alt residue.
        # Genuinely uncertain — the guard should not reject; threshold decides.
        # Here similarity is high enough at 0.7 threshold to count as fuzzy.
        match, score, match_type = find_best_match(
            "p.Asp501Asn", ["p.Asp501Arg"], threshold=0.7
        )
        assert match == "p.Asp501Arg"
        assert match_type == "fuzzy"
        assert score >= 0.7

    def test_accepts_same_position_aa_conversion(self):
        # 1-letter <-> 3-letter notation at same position is exact, not fuzzy.
        match, score, match_type = find_best_match(
            "p.R123H", ["p.Arg123His"], threshold=0.7
        )
        assert match == "p.Arg123His"
        assert match_type == "exact"
        assert score == 1.0

    def test_does_not_block_close_positions(self):
        # Even adjacent codons must be rejected when positions differ. A 561V
        # vs A562V are different variants; the matcher should not collapse them.
        match, _, match_type = find_best_match(
            "p.Ala561Val", ["p.Ala562Val"], threshold=0.7
        )
        assert match is None
        assert match_type == "none"

    def test_handles_unparseable_positions(self):
        # IVS notation has no extractable protein position; the guard should
        # fall through to similarity rather than reject blindly.
        match, _, _ = find_best_match("IVS3-2A>G", ["IVS3-2A>T"], threshold=0.7)
        # Both have position digits in the canonical/normalized form ("3", "2")
        # so the guard sees them as matching — similarity then decides.
        assert match == "IVS3-2A>T"

    def test_rejects_cdna_intronic_offset_subset_match(self):
        # c.1890+5G>A and c.1890G>A share the genomic base prefix but are not
        # the same cDNA coordinate. Subset positional matching is only safe for
        # protein frameshift suffixes, not cDNA offsets.
        match, _, match_type = find_best_match(
            "c.1890+5G>A", ["c.1890G>A"], threshold=0.7
        )
        assert match is None
        assert match_type == "none"


class TestGreedyOneToOne:
    """Each SQLite variant must be claimed by at most one Excel row."""

    def _build_inputs(self, excel_rows, sqlite_rows):
        excel_df = pd.DataFrame(excel_rows)
        detected = {
            "pmid": "PMID",
            "variant": "protein_change",
            "carriers_total": "carriers",
            "affected_count": "affected",
            "unaffected_count": "unaffected",
            "phenotype": "phenotype",
        }
        excel_data = aggregate_excel_data(excel_df, detected)

        sqlite_df = pd.DataFrame(sqlite_rows)
        sqlite_data = aggregate_sqlite_data(sqlite_df)
        return excel_data, sqlite_data

    def test_one_sqlite_not_double_counted(self):
        # Same PMID has gold G572S and G628S; SQLite has only G628S. Before the
        # fix, G628S matched G628S exactly AND was fuzzy-claimed by G572S, so
        # the same SQLite row counted as two matches. After: exactly one.
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["29650123", "29650123"],
                "protein_change": ["p.Gly572Ser", "p.Gly628Ser"],
                "carriers": [1, 1],
                "affected": [1, 1],
                "unaffected": [0, 0],
                "phenotype": ["LQT2", "LQT2"],
            },
            sqlite_rows={
                "pmid": ["29650123"],
                "variant": ["p.Gly628Ser"],
                "protein_notation": ["p.Gly628Ser"],
                "cdna_notation": [None],
                "carriers_total": [1],
                "affected_count": [1],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "fuzzy", 0.70)
        matched = [
            r
            for r in results
            if r.match_type in ("exact", "fuzzy") and not r.missing_in_excel
        ]
        # Exactly one Excel row should match (G628S), the other should be
        # missing-in-sqlite. No double-count.
        assert len(matched) == 1
        assert matched[0].excel_variant_norm == "G628S"

        missing = [r for r in results if r.missing_in_sqlite]
        assert len(missing) == 1
        assert missing[0].excel_variant_norm == "G572S"

    def test_exact_wins_over_fuzzy_for_same_candidate(self):
        # Gold has both R123H (exact match in sqlite) and R123Q (fuzzy candidate
        # against the same sqlite row). Pass-1 must claim the exact match first
        # so the fuzzy attempt sees an empty pool.
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["111", "111"],
                "protein_change": ["p.Arg123Gln", "p.Arg123His"],
                "carriers": [1, 1],
                "affected": [1, 1],
                "unaffected": [0, 0],
                "phenotype": ["LQT2", "LQT2"],
            },
            sqlite_rows={
                "pmid": ["111"],
                "variant": ["p.Arg123His"],
                "protein_notation": ["p.Arg123His"],
                "cdna_notation": [None],
                "carriers_total": [1],
                "affected_count": [1],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "fuzzy", 0.70)
        # R123H matched exactly; R123Q has no remaining candidate → missing.
        exact = [r for r in results if r.match_type == "exact"]
        assert len(exact) == 1
        assert exact[0].excel_variant_norm == "R123H"

        missing = [r for r in results if r.missing_in_sqlite]
        assert len(missing) == 1
        assert missing[0].excel_variant_norm == "R123Q"


class TestCDNAProteinBridge:
    """Gold protein indels can bridge to SQLite-only cDNA notation."""

    def _build_inputs(self, excel_rows, sqlite_rows):
        excel_df = pd.DataFrame(excel_rows)
        detected = {
            "pmid": "PMID",
            "variant": "protein_change",
            "carriers_total": "carriers",
            "affected_count": "affected",
            "unaffected_count": "unaffected",
            "phenotype": "phenotype",
        }
        excel_data = aggregate_excel_data(excel_df, detected)

        sqlite_df = pd.DataFrame(sqlite_rows)
        sqlite_data = aggregate_sqlite_data(sqlite_df)
        return excel_data, sqlite_data

    def test_matches_cdna_notation_hidden_by_sqlite_protein_display_key(self):
        # SQLite aggregation uses protein_notation as the display variant when
        # both protein and cDNA are present. The comparison still needs to see
        # the stored cDNA notation so gold cDNA rows do not become false misses.
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["19279983"],
                "protein_change": ["c.3480delT"],
                "carriers": [3],
                "affected": [2],
                "unaffected": [1],
                "phenotype": ["BrS"],
            },
            sqlite_rows={
                "pmid": ["19279983"],
                "variant": ["p.Phe1160Leufs*"],
                "protein_notation": ["p.Phe1160Leufs*"],
                "cdna_notation": ["c.3480delT"],
                "carriers_total": [3],
                "affected_count": [2],
                "unaffected_count": [1],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "fuzzy", 0.80)
        matched = [r for r in results if r.match_type == "exact"]
        assert len(matched) == 1
        assert matched[0].excel_variant_raw == "c.3480delT"
        assert matched[0].sqlite_variant_raw == "p.Phe1160Leufs*"

    def test_bridges_gold_cdna_to_sqlite_protein_frameshift(self):
        # Gold cDNA c.5464_5467delTCTG maps to the Leu1821/1822 frameshift
        # boundary. SQLite may only expose the protein notation as its primary
        # variant display key.
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["17897635"],
                "protein_change": ["c.5464_5467delTCTG"],
                "carriers": [1],
                "affected": [1],
                "unaffected": [0],
                "phenotype": ["BrS"],
            },
            sqlite_rows={
                "pmid": ["17897635"],
                "variant": ["p.Leu1821fs*10"],
                "protein_notation": ["p.Leu1821fs*10"],
                "cdna_notation": [None],
                "carriers_total": [1],
                "affected_count": [1],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "exact", 0.85)
        bridged = [r for r in results if r.match_type == "exact_cdna_bridge"]
        assert len(bridged) == 1
        assert bridged[0].sqlite_variant_raw == "p.Leu1821fs*10"

    def test_bridges_frameshift_via_cdna_position(self):
        # Gold: R281fsX (protein codon 281).
        # SQLite: only c.842dupG stored — implies codon ceil(842/3) = 281.
        # Bridge should match these even though the strings share nothing.
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["29622001"],
                "protein_change": ["R281fsX"],
                "carriers": [1],
                "affected": [1],
                "unaffected": [0],
                "phenotype": ["LQT2"],
            },
            sqlite_rows={
                "pmid": ["29622001"],
                # variant column falls back to cdna_notation when protein is NULL
                "variant": ["c.842dupG"],
                "protein_notation": [None],
                "cdna_notation": ["c.842dupG"],
                "carriers_total": [1],
                "affected_count": [1],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "fuzzy", 0.80)
        matched = [r for r in results if r.match_type.endswith("cdna_bridge")]
        assert len(matched) == 1
        assert matched[0].excel_variant_norm == "R281fsX"
        assert matched[0].sqlite_variant_raw == "c.842dupG"

    def test_bridge_works_in_exact_mode(self):
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["111"],
                "protein_change": ["A121fsX"],
                "carriers": [1],
                "affected": [1],
                "unaffected": [0],
                "phenotype": ["LQT2"],
            },
            sqlite_rows={
                "pmid": ["111"],
                # codon 121 -> cdna 361-363
                "variant": ["c.362delA"],
                "protein_notation": [None],
                "cdna_notation": ["c.362delA"],
                "carriers_total": [1],
                "affected_count": [1],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "exact", 0.85)
        matched = [r for r in results if r.match_type == "exact_cdna_bridge"]
        assert len(matched) == 1

    def test_bridge_matches_protein_position_inside_cdna_indel_range(self):
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["29622001"],
                "protein_change": ["P1034fsX"],
                "carriers": [1],
                "affected": [1],
                "unaffected": [0],
                "phenotype": ["LQT2"],
            },
            sqlite_rows={
                "pmid": ["29622001"],
                # c.3093_3106 spans protein positions 1031-1036.
                "variant": ["c.3093_3106del"],
                "protein_notation": [None],
                "cdna_notation": ["c.3093_3106del"],
                "carriers_total": [1],
                "affected_count": [1],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "exact", 0.85)
        matched = [r for r in results if r.match_type == "exact_cdna_bridge"]
        assert len(matched) == 1
        assert matched[0].sqlite_variant_raw == "c.3093_3106del"

    def test_bridge_rejects_wrong_codon(self):
        # Gold codon 281; cDNA position 100 → codon 34. Must not bridge.
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["111"],
                "protein_change": ["R281fsX"],
                "carriers": [1],
                "affected": [1],
                "unaffected": [0],
                "phenotype": ["LQT2"],
            },
            sqlite_rows={
                "pmid": ["111"],
                "variant": ["c.100dupG"],
                "protein_notation": [None],
                "cdna_notation": ["c.100dupG"],
                "carriers_total": [1],
                "affected_count": [1],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "fuzzy", 0.80)
        bridged = [r for r in results if r.match_type.endswith("cdna_bridge")]
        assert len(bridged) == 0
        # R281fsX should be marked missing-in-sqlite
        missing = [r for r in results if r.missing_in_sqlite]
        assert any(r.excel_variant_norm == "R281fsX" for r in missing)

    def test_bridge_skips_non_indel_cdna(self):
        # Gold is an indel; SQLite cDNA is a substitution (c.XXXG>A). Bridge
        # is for indel<->indel only — substitutions cause missense/nonsense,
        # not frameshifts.
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["111"],
                "protein_change": ["R281fsX"],
                "carriers": [1],
                "affected": [1],
                "unaffected": [0],
                "phenotype": ["LQT2"],
            },
            sqlite_rows={
                "pmid": ["111"],
                "variant": ["c.842G>A"],
                "protein_notation": [None],
                "cdna_notation": ["c.842G>A"],
                "carriers_total": [1],
                "affected_count": [1],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "fuzzy", 0.80)
        bridged = [r for r in results if r.match_type.endswith("cdna_bridge")]
        assert len(bridged) == 0

    def test_bridge_does_not_fire_for_missense_gold(self):
        # Gold A561V (missense) must not be bridged to any cDNA — only indel
        # canonical forms (fsX, Del, Ins, dup) can use the bridge.
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["111"],
                "protein_change": ["A561V"],
                "carriers": [1],
                "affected": [1],
                "unaffected": [0],
                "phenotype": ["LQT2"],
            },
            sqlite_rows={
                "pmid": ["111"],
                "variant": ["c.1682delC"],
                "protein_notation": [None],
                "cdna_notation": ["c.1682delC"],
                "carriers_total": [1],
                "affected_count": [1],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "fuzzy", 0.80)
        bridged = [r for r in results if r.match_type.endswith("cdna_bridge")]
        assert len(bridged) == 0

    def test_bridge_rejects_intronic_offset_deletion(self):
        # Gold c.5256_5278-2757del is a multi-kb intron-spanning structural
        # deletion; the -2757 offset makes codon math meaningless. It must NOT
        # bridge to an unrelated coding frameshift (c.5266dupC = p.Gln1756fs)
        # that merely happens to sit inside the naively-computed codon window.
        # (Real BRCA1 33468216 case: the false match inflated carriers MAE by
        # 210 by pairing a 1-carrier deletion with a 211-carrier duplication.)
        excel_data, sqlite_data = self._build_inputs(
            excel_rows={
                "PMID": ["111"],
                "protein_change": ["c.5256_5278-2757del"],
                "carriers": [1],
                "affected": [1],
                "unaffected": [0],
                "phenotype": ["HBOC"],
            },
            sqlite_rows={
                "pmid": ["111"],
                "variant": ["p.Gln1756fs"],
                "protein_notation": ["p.Gln1756fs"],
                "cdna_notation": ["c.5266dupC"],
                "carriers_total": [211],
                "affected_count": [211],
                "unaffected_count": [0],
                "uncertain_count": [0],
            },
        )

        results = compare_data(excel_data, sqlite_data, "fuzzy", 0.80)
        bridged = [r for r in results if r.match_type.endswith("cdna_bridge")]
        assert len(bridged) == 0
        gold_row = [r for r in results if r.excel_variant_raw == "c.5256_5278-2757del"]
        assert gold_row and gold_row[0].missing_in_sqlite is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
