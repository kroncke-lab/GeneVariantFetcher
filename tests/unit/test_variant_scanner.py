"""
Tests for the variant scanner module.

Run with: pytest tests/unit/test_variant_scanner.py -v
"""

import pytest

from utils.variant_scanner import (
    ScannedVariant,
    ScanResult,
    VariantScanner,
    merge_scanner_results,
    scan_document_for_variants,
)


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def scanner():
    """Default KCNH2 scanner."""
    return VariantScanner("KCNH2")


@pytest.fixture
def scn5a_scanner():
    """SCN5A scanner for cross-gene tests."""
    return VariantScanner("SCN5A")


# Base test text from variant_scanner.py __main__ block
TEST_TEXT = """
The patient carried the p.Arg534Cys mutation, also known as R534C.
This variant c.1600C>T was previously reported in LQT2 families.
We also identified A561V, Gly628Ser, and the frameshift mutation L987fsX.
The IVS9+1G>A splice variant was found in 3 families.
The intronic variant c.2398+1G>A causes splicing defects.
The W1001X nonsense mutation leads to truncation.
In Table 2, we list the mutations: T613M, N470D, and the deletion p.Leu552del.
Also found: c.526C>T (R176W), p.Thr613Met, and c.1234del.
"""


# =============================================================================
# TestScannedVariant
# =============================================================================


class TestScannedVariant:
    """Tests for the ScannedVariant dataclass."""

    def test_to_dict_contains_expected_keys(self):
        sv = ScannedVariant(
            raw_text="p.Arg534Cys",
            normalized="R534C",
            variant_type="missense",
            notation_type="protein",
            position=534,
            context="the p.Arg534Cys mutation",
            confidence=0.95,
            source="protein_hgvs_full",
        )
        d = sv.to_dict()
        assert d["raw_text"] == "p.Arg534Cys"
        assert d["normalized"] == "R534C"
        assert d["variant_type"] == "missense"
        assert d["notation_type"] == "protein"
        assert d["position"] == 534
        assert d["confidence"] == 0.95
        assert d["source"] == "protein_hgvs_full"

    def test_to_dict_excludes_context(self):
        sv = ScannedVariant(
            raw_text="R534C",
            normalized="R534C",
            variant_type="missense",
            notation_type="protein",
            position=534,
            context="some context text here",
            confidence=0.70,
            source="protein_single_letter",
        )
        d = sv.to_dict()
        assert "context" not in d


# =============================================================================
# TestScanResult
# =============================================================================


class TestScanResult:
    """Tests for the ScanResult dataclass."""

    def test_empty_hints(self):
        result = ScanResult()
        assert result.get_hints_for_prompt() == ""

    def test_hints_format(self):
        result = ScanResult()
        result.variants = [
            ScannedVariant("R534C", "R534C", "missense", "protein", 534, "", 0.95, "test"),
            ScannedVariant("A561V", "A561V", "missense", "protein", 561, "", 0.70, "test"),
        ]
        result.stats["gene"] = "KCNH2"
        hints = result.get_hints_for_prompt()
        assert "PRE-SCANNED VARIANT HINTS" in hints
        assert "R534C" in hints
        assert "A561V" in hints
        assert "HIGH" in hints  # R534C has 0.95 confidence
        assert "KCNH2" in hints

    def test_hints_max_limit(self):
        result = ScanResult()
        result.stats["gene"] = "KCNH2"
        # Create 60 variants
        for i in range(60):
            result.variants.append(
                ScannedVariant(
                    f"A{100+i}V", f"A{100+i}V", "missense", "protein",
                    100 + i, "", 0.70, "test",
                )
            )
        hints = result.get_hints_for_prompt(max_hints=10)
        # Should only contain 10 numbered entries
        assert "10." in hints
        assert "11." not in hints

    def test_hints_dedup(self):
        result = ScanResult()
        result.stats["gene"] = "KCNH2"
        result.variants = [
            ScannedVariant("R534C", "R534C", "missense", "protein", 534, "", 0.95, "a"),
            ScannedVariant("p.R534C", "R534C", "missense", "protein", 534, "", 0.90, "b"),
        ]
        hints = result.get_hints_for_prompt()
        # After dedup, only 1 unique R534C should appear in the list
        assert "  1. R534C" in hints
        assert "  2. R534C" not in hints

    def test_to_variant_dicts(self):
        result = ScanResult()
        result.variants = [
            ScannedVariant("R534C", "R534C", "missense", "protein", 534, "", 0.95, "test"),
            ScannedVariant("c.1600C>T", "c.1600C>T", "substitution", "cdna", 1600, "", 0.95, "test"),
        ]
        dicts = result.to_variant_dicts("KCNH2")
        assert len(dicts) == 2
        assert dicts[0]["gene_symbol"] == "KCNH2"
        assert dicts[0]["protein_notation"] == "R534C"
        assert dicts[0]["cdna_notation"] is None
        assert dicts[1]["cdna_notation"] == "c.1600C>T"
        assert dicts[1]["protein_notation"] is None
        assert dicts[0]["evidence_level"] == "scanner"

    def test_to_variant_dicts_dedup(self):
        result = ScanResult()
        result.variants = [
            ScannedVariant("R534C", "R534C", "missense", "protein", 534, "", 0.95, "a"),
            ScannedVariant("p.R534C", "R534C", "missense", "protein", 534, "", 0.90, "b"),
        ]
        dicts = result.to_variant_dicts("KCNH2")
        assert len(dicts) == 1


# =============================================================================
# TestProteinVariantScanning
# =============================================================================


class TestProteinVariantScanning:
    """Parametrized tests for each protein pattern."""

    @pytest.mark.parametrize(
        "text, expected_normalized, expected_confidence",
        [
            ("p.Arg534Cys", "R534C", 0.95),
            ("p.Ala561Val", "A561V", 0.95),
            ("p.Gly628Ser", "G628S", 0.95),
            ("p.Thr613Met", "T613M", 0.95),
        ],
        ids=["Arg534Cys", "Ala561Val", "Gly628Ser", "Thr613Met"],
    )
    def test_full_hgvs(self, scanner, text, expected_normalized, expected_confidence):
        result = scanner.scan(f"found the {text} mutation")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms
        matching = [v for v in result.variants if v.normalized == expected_normalized]
        assert any(v.confidence >= expected_confidence for v in matching)

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("p.R534C", "R534C"),
            ("p.A561V", "A561V"),
            ("p.L987fs", "L987fsX"),
        ],
        ids=["pR534C", "pA561V", "pL987fs"],
    )
    def test_short_hgvs(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the {text} variant was detected")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("Arg534Cys", "R534C"),
            ("Gly628Ser", "G628S"),
        ],
        ids=["Arg534Cys_nop", "Gly628Ser_nop"],
    )
    def test_three_letter(self, scanner, text, expected_normalized):
        result = scanner.scan(f"we identified {text} in families")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("R534C", "R534C"),
            ("A561V", "A561V"),
            ("T613M", "T613M"),
            ("N470D", "N470D"),
        ],
        ids=["R534C", "A561V", "T613M", "N470D"],
    )
    def test_single_letter(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the {text} mutation was found in this study")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("HERGG604S", "G604S"),
            ("KCNH2A561V", "A561V"),
            ("hERGT613M", "T613M"),
            ("KCNH2-R534C", "R534C"),
        ],
        ids=["HERGG604S", "KCNH2A561V", "hERGT613M", "KCNH2-R534C"],
    )
    def test_concatenated_gene_variant(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the {text} mutant channel")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("L987fsX", "L987fsX"),
            ("p.Leu987fs", "L987fsX"),
            ("L987fs*10", "L987fsX"),
        ],
        ids=["L987fsX", "pLeu987fs", "L987fs_star10"],
    )
    def test_frameshifts(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the frameshift {text} was identified")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("W1001X", "W1001X"),
            ("p.Trp1001Ter", "W1001X"),
        ],
        ids=["W1001X", "pTrp1001Ter"],
    )
    def test_nonsense(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the nonsense mutation {text} was found")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("L552del", "L552del"),
            ("p.Leu552del", "L552del"),
        ],
        ids=["L552del", "pLeu552del"],
    )
    def test_deletions(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the deletion {text} was identified")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms


# =============================================================================
# TestCdnaVariantScanning
# =============================================================================


class TestCdnaVariantScanning:
    """Parametrized tests for cDNA patterns."""

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("c.1600C>T", "c.1600C>T"),
            ("c.526C>T", "c.526C>T"),
        ],
        ids=["c1600CT", "c526CT"],
    )
    def test_standard_cdna(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the variant {text} was detected")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("c.1234+1G>A", "c.1234+1G>A"),
            ("c.2398+1G>A", "c.2398+1G>A"),
        ],
        ids=["c1234_plus1GA", "c2398_plus1GA"],
    )
    def test_intronic_cdna(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the splice variant {text} causes skipping")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("c.1234del", "c.1234del"),
            ("c.1234_1235delAG", "c.1234_1235delAG"),
            ("c.1234dup", "c.1234dup"),
        ],
        ids=["c1234del", "c1234_1235delAG", "c1234dup"],
    )
    def test_cdna_indels(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the {text} variant was identified")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms


# =============================================================================
# TestSpliceVariantScanning
# =============================================================================


class TestSpliceVariantScanning:
    """Tests for IVS splice variant patterns."""

    @pytest.mark.parametrize(
        "text, expected_normalized",
        [
            ("IVS9+1G>A", "IVS9+1G>A"),
            ("IVS5-2A>G", "IVS5-2A>G"),
        ],
        ids=["IVS9_plus1GA", "IVS5_minus2AG"],
    )
    def test_ivs_patterns(self, scanner, text, expected_normalized):
        result = scanner.scan(f"the splice variant {text} was found")
        norms = {v.normalized for v in result.variants}
        assert expected_normalized in norms


# =============================================================================
# TestNarrativeScanning
# =============================================================================


class TestNarrativeScanning:
    """Tests for variants in narrative context."""

    def test_the_x_mutation(self, scanner):
        result = scanner.scan("the R534C mutation was found in three families")
        norms = {v.normalized for v in result.variants}
        assert "R534C" in norms

    def test_carrying_variant(self, scanner):
        result = scanner.scan("patients carrying the p.Arg534Cys variant were studied")
        norms = {v.normalized for v in result.variants}
        assert "R534C" in norms

    def test_identified_variant(self, scanner):
        result = scanner.scan("we identified the A561V in probands")
        norms = {v.normalized for v in result.variants}
        assert "A561V" in norms

    def test_cdna_mutation_context(self, scanner):
        result = scanner.scan("the c.1600C>T mutation was identified in the proband")
        norms = {v.normalized for v in result.variants}
        assert "c.1600C>T" in norms


# =============================================================================
# TestFalsePositiveFiltering
# =============================================================================


class TestFalsePositiveFiltering:
    """Tests that false positives are filtered out."""

    @pytest.mark.parametrize(
        "text",
        [
            "Table 1 shows results",
            "see Figure 2 for pedigree",
            "HEK293 cells were transfected",
            "37C incubation temperature",
        ],
        ids=["table_ref", "figure_ref", "cell_line", "temperature"],
    )
    def test_false_positives_excluded(self, scanner, text):
        result = scanner.scan(text)
        # These shouldn't produce any protein variants from the reference text
        protein_variants = [v for v in result.variants if v.notation_type == "protein"]
        for v in protein_variants:
            # None of these should be the false positive patterns themselves
            assert v.normalized not in {"S1", "T1", "F2", "T2"}

    def test_short_codes_filtered(self, scanner):
        s = VariantScanner("KCNH2")
        assert s._is_false_positive("S1")
        assert s._is_false_positive("T2")
        assert s._is_false_positive("A1")

    def test_valid_variants_not_filtered(self, scanner):
        assert not scanner._is_false_positive("R534C")
        assert not scanner._is_false_positive("p.Arg534Cys")


# =============================================================================
# TestNonTargetHotspotFiltering
# =============================================================================


class TestNonTargetHotspotFiltering:
    """Tests that known non-target gene hotspots are filtered."""

    @pytest.mark.parametrize(
        "variant",
        ["R175H", "R248W", "G12D", "G12V", "V600E", "E545K"],
        ids=["TP53_R175H", "TP53_R248W", "KRAS_G12D", "KRAS_G12V", "BRAF_V600E", "PIK3CA_E545K"],
    )
    def test_hotspots_filtered(self, scanner, variant):
        result = scanner.scan(f"the {variant} mutation was identified")
        norms = {v.normalized for v in result.variants}
        assert variant not in norms


# =============================================================================
# TestPositionValidation
# =============================================================================


class TestPositionValidation:
    """Tests that variants beyond protein length are filtered."""

    def test_valid_position_kept(self, scanner):
        # KCNH2 protein length is 1159
        result = scanner.scan("the R534C mutation was found")
        norms = {v.normalized for v in result.variants}
        assert "R534C" in norms

    def test_position_at_boundary_kept(self, scanner):
        result = scanner.scan("the L1159X nonsense mutation truncates the protein")
        norms = {v.normalized for v in result.variants}
        assert "L1159X" in norms

    def test_position_beyond_length_filtered(self, scanner):
        # Position 2000 > KCNH2 length of 1159
        result = scanner.scan("the R2000C mutation was found")
        norms = {v.normalized for v in result.variants}
        assert "R2000C" not in norms

    def test_cdna_position_not_filtered_by_protein_length(self, scanner):
        # cDNA positions can be larger than protein length
        result = scanner.scan("c.3000C>T was detected")
        norms = {v.normalized for v in result.variants}
        assert "c.3000C>T" in norms


# =============================================================================
# TestUnicodeNormalization
# =============================================================================


class TestUnicodeNormalization:
    """Tests that Unicode arrow characters are normalized."""

    def test_right_arrow_normalized(self, scanner):
        # \u2192 = →
        result = scanner.scan("c.1600C\u2192T was found in the proband")
        norms = {v.normalized for v in result.variants}
        assert "c.1600C>T" in norms

    def test_heavy_arrow_normalized(self, scanner):
        # \u21d2 = ⇒
        result = scanner.scan("c.526C\u21d2T change was identified")
        norms = {v.normalized for v in result.variants}
        assert "c.526C>T" in norms


# =============================================================================
# TestMergeScannerResults
# =============================================================================


class TestMergeScannerResults:
    """Tests for merge_scanner_results()."""

    def test_adds_new_variants(self):
        extracted = {
            "variants": [
                {"gene_symbol": "KCNH2", "protein_notation": "R534C", "cdna_notation": None},
            ],
            "extraction_metadata": {"total_variants_found": 1},
        }
        scan_result = ScanResult()
        scan_result.variants = [
            ScannedVariant("A561V", "A561V", "missense", "protein", 561, "", 0.80, "test"),
        ]

        merged = merge_scanner_results(extracted, scan_result, "KCNH2", min_confidence=0.5)
        protein_notations = [v.get("protein_notation") for v in merged["variants"]]
        assert "A561V" in protein_notations
        assert len(merged["variants"]) == 2

    def test_skips_existing_variants(self):
        extracted = {
            "variants": [
                {"gene_symbol": "KCNH2", "protein_notation": "R534C", "cdna_notation": None},
            ],
            "extraction_metadata": {"total_variants_found": 1},
        }
        scan_result = ScanResult()
        scan_result.variants = [
            ScannedVariant("R534C", "R534C", "missense", "protein", 534, "", 0.95, "test"),
        ]

        merged = merge_scanner_results(extracted, scan_result, "KCNH2", min_confidence=0.5)
        assert len(merged["variants"]) == 1

    def test_respects_min_confidence(self):
        extracted = {
            "variants": [],
            "extraction_metadata": {"total_variants_found": 0},
        }
        scan_result = ScanResult()
        scan_result.variants = [
            ScannedVariant("A561V", "A561V", "missense", "protein", 561, "", 0.30, "test"),
        ]

        merged = merge_scanner_results(extracted, scan_result, "KCNH2", min_confidence=0.5)
        assert len(merged["variants"]) == 0

    def test_updates_metadata(self):
        extracted = {
            "variants": [],
            "extraction_metadata": {"total_variants_found": 0},
        }
        scan_result = ScanResult()
        scan_result.variants = [
            ScannedVariant("A561V", "A561V", "missense", "protein", 561, "", 0.80, "test"),
        ]

        merged = merge_scanner_results(extracted, scan_result, "KCNH2", min_confidence=0.5)
        assert merged["extraction_metadata"]["scanner_added"] == 1
        assert merged["extraction_metadata"]["total_variants_found"] == 1


# =============================================================================
# TestScanDocumentConvenience
# =============================================================================


class TestScanDocumentConvenience:
    """Tests for scan_document_for_variants() wrapper."""

    def test_basic_scan(self):
        result = scan_document_for_variants(TEST_TEXT, "KCNH2")
        assert len(result.variants) > 0
        assert result.stats["gene"] == "KCNH2"

    def test_finds_expected_variants(self):
        result = scan_document_for_variants(TEST_TEXT, "KCNH2")
        norms = {v.normalized for v in result.variants}
        # These should all be found in the base test text
        assert "R534C" in norms
        assert "A561V" in norms
        assert "T613M" in norms
        assert "c.1600C>T" in norms

    def test_stats_populated(self):
        result = scan_document_for_variants(TEST_TEXT, "KCNH2")
        assert "text_length" in result.stats
        assert "unique_variants" in result.stats
        assert result.stats["unique_variants"] > 0

    def test_custom_source_label(self):
        result = scan_document_for_variants("A561V mutation", "KCNH2", source="supplement")
        assert result.stats.get("source") == "supplement"

    def test_different_gene(self):
        text = "the p.Arg1193Gln variant in SCN5A"
        result = scan_document_for_variants(text, "SCN5A")
        norms = {v.normalized for v in result.variants}
        assert "R1193Q" in norms
