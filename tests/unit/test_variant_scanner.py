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
            ScannedVariant(
                "R534C", "R534C", "missense", "protein", 534, "", 0.95, "test"
            ),
            ScannedVariant(
                "A561V", "A561V", "missense", "protein", 561, "", 0.70, "test"
            ),
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
                    f"A{100 + i}V",
                    f"A{100 + i}V",
                    "missense",
                    "protein",
                    100 + i,
                    "",
                    0.70,
                    "test",
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
            ScannedVariant(
                "p.R534C", "R534C", "missense", "protein", 534, "", 0.90, "b"
            ),
        ]
        hints = result.get_hints_for_prompt()
        # After dedup, only 1 unique R534C should appear in the list
        assert "  1. R534C" in hints
        assert "  2. R534C" not in hints

    def test_to_variant_dicts(self):
        result = ScanResult()
        result.variants = [
            ScannedVariant(
                "R534C", "R534C", "missense", "protein", 534, "", 0.95, "test"
            ),
            ScannedVariant(
                "c.1600C>T", "c.1600C>T", "substitution", "cdna", 1600, "", 0.95, "test"
            ),
        ]
        dicts = result.to_variant_dicts("KCNH2")
        assert len(dicts) == 2
        assert dicts[0]["gene_symbol"] == "KCNH2"
        assert dicts[0]["protein_notation"] == "R534C"
        assert dicts[0]["cdna_notation"] is None
        assert dicts[1]["cdna_notation"] == "c.1600C>T"
        assert dicts[1]["protein_notation"] is None
        assert dicts[0]["evidence_level"] == "scanner"
        assert dicts[0]["key_quotes"] == ["R534C"]
        assert dicts[1]["key_quotes"] == ["c.1600C>T"]

    def test_to_variant_dicts_dedup(self):
        result = ScanResult()
        result.variants = [
            ScannedVariant("R534C", "R534C", "missense", "protein", 534, "", 0.95, "a"),
            ScannedVariant(
                "p.R534C", "R534C", "missense", "protein", 534, "", 0.90, "b"
            ),
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

    @pytest.mark.parametrize(
        "text",
        ["p.Lys1505_Gln1507del", "K1505_Q1507del"],
        ids=["pLys1505_Gln1507del", "K1505_Q1507del"],
    )
    def test_scn5a_range_deletions(self, scn5a_scanner, text):
        result = scn5a_scanner.scan(f"SCN5A deletion {text} was identified")
        norms = {v.normalized for v in result.variants}
        assert "K1505_Q1507del" in norms

    def test_range_deletion_does_not_emit_suffix_only_deletion(self, scn5a_scanner):
        result = scn5a_scanner.scan(
            "SCN5A p.Lys1505_Gln1507del was identified in Table 2."
        )
        norms = {v.normalized for v in result.variants}
        assert "K1505_Q1507del" in norms
        assert "Q1507del" not in norms


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
# TestGeneContextFiltering
# =============================================================================


class TestGeneContextFiltering:
    """Tests that multi-gene table rows do not bleed into the target gene."""

    def test_scanner_filters_rows_explicitly_labeled_with_other_gene(
        self, scn5a_scanner
    ):
        text = """
        | Patient | Gene | Aa change (p.) |
        | 5 | KCNH2 | Ala561Val |
        | 19 | SCN5A | Glu1784Lys |
        """
        result = scn5a_scanner.scan(text)
        norms = {v.normalized for v in result.variants}
        assert "E1784K" in norms
        assert "A561V" not in norms

    def test_target_gene_context_keeps_variant_when_other_gene_also_named(
        self, scn5a_scanner
    ):
        text = "SCN5A and KCNH2 were screened; the SCN5A Glu1784Lys row is listed."
        result = scn5a_scanner.scan(text)
        norms = {v.normalized for v in result.variants}
        assert "E1784K" in norms

    def test_nearest_gene_assignment_filters_other_gene_variant_in_sentence(
        self, scn5a_scanner
    ):
        text = (
            "KCNQ1 p.(Pro73Thr) was previously described as a VUS, "
            "and SCN5A p.(Pro1008Ser) had been associated with arrhythmia."
        )
        result = scn5a_scanner.scan(text)
        norms = {v.normalized for v in result.variants}
        assert "P1008S" in norms
        assert "P73T" not in norms

    def test_nearest_gene_assignment_handles_wrapped_multigene_table_rows(
        self, scn5a_scanner
    ):
        text = """
        | 33 | KCNQ1 | 1 | Cytoplasmic |  |  |  |
        | N-term | Pro73Thr | 217C>A | Missense | het |
        | 46 | SCN5A | 12 | Cytoplasmic |  |  |  |
        | C-term | Arg620Cys | 1858 C>T | Missense | het |
        """
        result = scn5a_scanner.scan(text)
        norms = {v.normalized for v in result.variants}
        assert "R620C" in norms
        assert "P73T" not in norms


# =============================================================================
# TestNonTargetHotspotFiltering
# =============================================================================


class TestNonTargetHotspotFiltering:
    """Tests that known non-target gene hotspots are filtered."""

    @pytest.mark.parametrize(
        "variant",
        ["R175H", "R248W", "G12D", "G12V", "V600E", "E545K"],
        ids=[
            "TP53_R175H",
            "TP53_R248W",
            "KRAS_G12D",
            "KRAS_G12V",
            "BRAF_V600E",
            "PIK3CA_E545K",
        ],
    )
    def test_hotspots_filtered(self, scanner, variant):
        result = scanner.scan(f"the {variant} mutation was identified")
        norms = {v.normalized for v in result.variants}
        assert variant not in norms

    @pytest.mark.parametrize(
        ("gene", "variant"),
        [("TP53", "R248Q"), ("KRAS", "G12D"), ("BRAF", "V600E"), ("PIK3CA", "H1047R")],
    )
    def test_hotspots_kept_when_their_gene_is_target(self, gene, variant):
        scanner = VariantScanner(gene)
        result = scanner.scan(f"the {variant} mutation was identified")
        norms = {v.normalized for v in result.variants}
        assert variant in norms


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

    def test_scanner_gene_context_uses_default_aliases(self):
        scanner = VariantScanner("MYBPC3")
        context = "Table 1 | cMyBP-C | p.Arg502Trp | 4 carriers |"

        assert scanner._context_mentions_gene(context, "MYBPC3") is True
        assert scanner._gene_assigned_to_variant(context, "p.Arg502Trp") == "MYBPC3"


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
                {
                    "gene_symbol": "KCNH2",
                    "protein_notation": "R534C",
                    "cdna_notation": None,
                },
            ],
            "extraction_metadata": {"total_variants_found": 1},
        }
        scan_result = ScanResult()
        scan_result.variants = [
            ScannedVariant(
                "A561V", "A561V", "missense", "protein", 561, "", 0.80, "test"
            ),
        ]

        merged = merge_scanner_results(
            extracted,
            scan_result,
            "KCNH2",
            min_confidence=0.5,
            document_text="In our KCNH2 cohort, A561V was identified in a proband.",
        )
        protein_notations = [v.get("protein_notation") for v in merged["variants"]]
        assert "A561V" in protein_notations
        assert len(merged["variants"]) == 2

    def test_skips_existing_variants(self):
        extracted = {
            "variants": [
                {
                    "gene_symbol": "KCNH2",
                    "protein_notation": "R534C",
                    "cdna_notation": None,
                },
            ],
            "extraction_metadata": {"total_variants_found": 1},
        }
        scan_result = ScanResult()
        scan_result.variants = [
            ScannedVariant(
                "R534C", "R534C", "missense", "protein", 534, "", 0.95, "test"
            ),
        ]

        merged = merge_scanner_results(
            extracted, scan_result, "KCNH2", min_confidence=0.5
        )
        assert len(merged["variants"]) == 1

    def test_respects_min_confidence(self):
        extracted = {
            "variants": [],
            "extraction_metadata": {"total_variants_found": 0},
        }
        scan_result = ScanResult()
        scan_result.variants = [
            ScannedVariant(
                "A561V", "A561V", "missense", "protein", 561, "", 0.30, "test"
            ),
        ]

        merged = merge_scanner_results(
            extracted, scan_result, "KCNH2", min_confidence=0.5
        )
        assert len(merged["variants"]) == 0

    def test_updates_metadata(self):
        extracted = {
            "variants": [],
            "extraction_metadata": {"total_variants_found": 0},
        }
        scan_result = ScanResult()
        scan_result.variants = [
            ScannedVariant(
                "A561V", "A561V", "missense", "protein", 561, "", 0.80, "test"
            ),
        ]

        merged = merge_scanner_results(
            extracted,
            scan_result,
            "KCNH2",
            min_confidence=0.5,
            document_text="In our KCNH2 cohort, A561V was identified in a proband.",
        )
        assert merged["extraction_metadata"]["scanner_added"] == 1
        assert merged["extraction_metadata"]["total_variants_found"] == 1

    def test_fails_closed_without_gene_or_current_study_support(self):
        extracted = {"variants": [], "extraction_metadata": {}}
        scan_result = ScanResult(
            variants=[
                ScannedVariant(
                    "A561V", "A561V", "missense", "protein", 561, "", 0.90, "test"
                )
            ]
        )

        merged = merge_scanner_results(
            extracted,
            scan_result,
            "KCNH2",
            document_text="A561V is a well-known comparator mutation.",
        )

        assert merged["variants"] == []
        assert merged["extraction_metadata"]["scanner_merge"]["skipped"][0][
            "reason"
        ] in {"gene_unattributed", "study_unattributed"}

    def test_rejects_wrong_gene_even_when_study_observed_it(self):
        text = (
            "We examined BRCA1 and BRCA2. "
            "The BRCA1 variant c.181T>G was identified in two patients."
        )
        scan_result = ScanResult(
            variants=[
                ScannedVariant(
                    "c.181T>G",
                    "c.181T>G",
                    "substitution",
                    "cdna",
                    181,
                    "",
                    0.95,
                    "test",
                )
            ]
        )

        merged = merge_scanner_results(
            {"variants": [], "extraction_metadata": {}},
            scan_result,
            "BRCA2",
            document_text=text,
        )

        assert merged["variants"] == []
        assert (
            merged["extraction_metadata"]["scanner_merge"]["skipped"][0]["reason"]
            == "wrong_gene"
        )

    def test_rejects_embedded_other_gene_prefix(self):
        text = (
            "Four variants found in three subjects included KCNH2_R176W and "
            "SCN5A_S1103Y."
        )
        scan_result = scan_document_for_variants(text, "SCN5A")

        merged = merge_scanner_results(
            {"variants": [], "extraction_metadata": {}},
            scan_result,
            "SCN5A",
            document_text=text,
        )

        assert not any(v.get("protein_notation") == "R176W" for v in merged["variants"])
        rejected = {
            item["variant"]: item["reason"]
            for item in merged["extraction_metadata"]["scanner_merge"]["skipped"]
        }
        assert rejected["R176W"] == "wrong_gene"

    def test_rejects_bibliography_and_unconfirmed_artifact(self):
        text = "\n".join(
            [
                "Christiansen M, Hedley P, et al. A founder family with p.F29L in KCNH2.",
                "Five RYR2 mutations (H469Y, L2299F) were identified in cases; repeated sequencing could not confirm these DNA artifacts.",
            ]
        )
        for gene in ("KCNH2", "RYR2"):
            scan_result = scan_document_for_variants(text, gene)
            merged = merge_scanner_results(
                {"variants": [], "extraction_metadata": {}},
                scan_result,
                gene,
                document_text=text,
            )
            assert merged["variants"] == []

    def test_rejects_background_and_compilation_mentions(self):
        text = "\n".join(
            [
                "BRCA2 recurrent mutations have been reported (c.658_659delGT).",
                "Table 3 Overview of BRCA2 variants | Study | Variant |",
                "| Smith et al. 2010 | c.3847_3848delGT | 1 |",
            ]
        )
        scan_result = scan_document_for_variants(text, "BRCA2")

        merged = merge_scanner_results(
            {"variants": [], "extraction_metadata": {}},
            scan_result,
            "BRCA2",
            document_text=text,
        )

        assert merged["variants"] == []
        reasons = {
            item["reason"]
            for item in merged["extraction_metadata"]["scanner_merge"]["skipped"]
        }
        assert reasons & {"background_mention", "table_like", "reference_list"}

    def test_any_independently_supported_mention_can_pass(self):
        text = "\n".join(
            [
                "The KCNH2 A561V variant was previously reported by Smith et al.",
                "In our KCNH2 cohort, A561V was identified in two probands.",
            ]
        )
        scan_result = scan_document_for_variants(text, "KCNH2")

        merged = merge_scanner_results(
            {"variants": [], "extraction_metadata": {}},
            scan_result,
            "KCNH2",
            document_text=text,
        )

        assert {v["protein_notation"] for v in merged["variants"]} == {"A561V"}

    def test_current_study_two_variant_list_can_pass(self):
        text = "In our RYR2 cohort, H469Y and L2299F were identified in two probands."
        scan_result = scan_document_for_variants(text, "RYR2")

        merged = merge_scanner_results(
            {"variants": [], "extraction_metadata": {}},
            scan_result,
            "RYR2",
            document_text=text,
        )

        assert {v["protein_notation"] for v in merged["variants"]} == {
            "H469Y",
            "L2299F",
        }

    def test_normalized_existing_identity_is_not_readded(self):
        extracted = {
            "variants": [
                {
                    "gene_symbol": "BRCA2",
                    "protein_notation": "p.F2638*",
                    "cdna_notation": "c.7913_7917delTTCCT",
                }
            ],
            "extraction_metadata": {},
        }
        scan_result = ScanResult(
            variants=[
                ScannedVariant(
                    "F2638*",
                    "F2638*",
                    "nonsense",
                    "protein",
                    2638,
                    "",
                    0.95,
                    "test",
                )
            ]
        )

        merged = merge_scanner_results(
            extracted,
            scan_result,
            "BRCA2",
            document_text="BRCA2 F2638* was identified in one patient.",
        )

        assert len(merged["variants"]) == 1
        assert (
            merged["extraction_metadata"]["scanner_merge"]["skipped"][0]["reason"]
            == "already_extracted"
        )

    def test_braca2_compilation_regression_keeps_only_authoritative_rows(self):
        text = "\n".join(
            [
                "BRCA1, BRCA2 and PALB2 were examined. The PALB2 mutation c.509_510delGA was identified in two patients.",
                "Some recurrent mutations of BRCA2 have been reported (c.658_659delGT, c.3847_3848delGT, c.5946delT).",
                "Table 1 Molecular characteristics BRCA1 c.181T>G p.C61G BRCA2 c.1310_1313delAAGA c.9371A>T c.9403delC PALB2 c.509_510delGA",
                "In ten BRCA2 patients, c.1310_1313delAAGA, c.6267_6269delinsC, c.7913_7917delTTCCT, c.9027delT, c.9371A>T, c.9403delC, and c.10095delCinsGAATTATATCT were identified.",
                "Table 3 Overview of BRCA2 mutations Smith et al. 2000 Jones et al. 2004 c.658_659delGT c.3847_3848delGT c.5239_5240insT c.5946delT Current study",
            ]
        )
        authoritative = [
            "c.1310_1313delAAGA",
            "c.6267_6269delinsC",
            "c.7913_7917delTTCCT",
            "c.9027delT",
            "c.9371A>T",
            "c.9403delC",
            "c.10095delCinsGAATTATATCT",
        ]
        extracted = {
            "variants": [
                {"gene_symbol": "BRCA2", "cdna_notation": notation}
                for notation in authoritative
            ],
            "extraction_metadata": {},
        }
        scan_result = scan_document_for_variants(text, "BRCA2")

        merged = merge_scanner_results(
            extracted,
            scan_result,
            "BRCA2",
            document_text=text,
        )

        assert len(merged["variants"]) == 7
        assert merged["extraction_metadata"]["scanner_added"] == 0
        merged_notations = {v.get("cdna_notation") for v in merged["variants"]}
        assert "c.509_510delGA" not in merged_notations
        assert "c.658_659delGT" not in merged_notations


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
        result = scan_document_for_variants(
            "A561V mutation", "KCNH2", source="supplement"
        )
        assert result.stats.get("source") == "supplement"

    def test_different_gene(self):
        text = "the p.Arg1193Gln variant in SCN5A"
        result = scan_document_for_variants(text, "SCN5A")
        norms = {v.normalized for v in result.variants}
        assert "R1193Q" in norms


# =============================================================================
# DOCUMENT-LEVEL GENE ATTRIBUTION
# =============================================================================


# Real papers put whole paragraphs between mentions. These fixtures pad to more
# than the ±240-character attribution window so each mention is judged locally,
# the way it is in a real full text.
_FILLER = (
    "Clinical evaluation included resting electrocardiography, exercise "
    "testing and Holter monitoring for every family member who consented. "
    "Corrected QT intervals were measured with Bazett's formula by two "
    "independent readers who were blinded to genotype, and disagreements "
    "were resolved by consensus review of the tracings. "
)


class TestDocumentGeneAttribution:
    """A ±50-char window often misses the gene label; the whole document does not.

    Papers that sequence two genes name the gene once per sentence rather than
    beside every repetition of the variant, so the per-hit conflict filter sees
    a window with no gene in it at all and keeps the token for whichever gene
    happens to be running.
    """

    TWO_GENE_TEXT = (
        "All exons of KCNH2 and KCNQ1 were sequenced in the probands. "
        + _FILLER
        + "DNA sequence analysis of the proband revealed a heterozygous "
        "transition at nucleotide 560 of KCNQ1, which results in a "
        "substitution of amino acid residue leucine by proline (L187P). "
        + _FILLER
        + "The L187P mutation is located between domains S2 and S3 of KCNQ1."
    )

    def test_other_genes_variant_is_not_scanned_for_the_target(self):
        norms = {
            v.normalized
            for v in scan_document_for_variants(self.TWO_GENE_TEXT, "KCNH2").variants
        }
        assert "L187P" not in norms

    def test_the_owning_gene_still_scans_it(self):
        norms = {
            v.normalized
            for v in scan_document_for_variants(self.TWO_GENE_TEXT, "KCNQ1").variants
        }
        assert "L187P" in norms

    def test_attribution_is_collected_across_every_mention(self):
        scanner = VariantScanner("KCNH2")
        assert scanner._document_gene_attribution(self.TWO_GENE_TEXT, "L187P") == {
            "KCNQ1"
        }
        assert scanner._document_assigns_variant_elsewhere(self.TWO_GENE_TEXT, "L187P")

    def test_ambiguous_window_abstains_instead_of_voting(self):
        """ "P297S KCNH2 and P1177L SCN5A" must not assign P1177L to KCNH2.

        The per-hit assigner prefers the nearest mention *before* the token, so
        this window reads as KCNH2 even though the label sits immediately after
        the variant. A window that names the target gene too abstains, which
        leaves the candidate in place.
        """
        text = (
            "In total, 5 of 44 cases carried a mutation in 1 of the 3 genes "
            "(R190W KCNQ1, F29L KCNH2, P297S KCNH2 and P1177L SCN5A)."
        )
        scanner = VariantScanner("SCN5A")
        assert scanner._document_gene_attribution(text, "P1177L") == set()
        assert not scanner._document_assigns_variant_elsewhere(text, "P1177L")

    def test_document_with_no_gene_label_keeps_everything(self):
        text = "The proband carried L187P and a second sequence change."
        scanner = VariantScanner("KCNH2")
        assert scanner._document_gene_attribution(text, "L187P") == set()
        assert not scanner._document_assigns_variant_elsewhere(text, "L187P")

    def test_target_gene_mention_keeps_its_own_variant(self):
        text = (
            "Sequencing of KCNH2 identified the A561V substitution in the "
            "proband. " + _FILLER + "The A561V carrier had a prolonged QTc."
        )
        scanner = VariantScanner("KCNH2")
        assert "KCNH2" in scanner._document_gene_attribution(text, "A561V")
        assert not scanner._document_assigns_variant_elsewhere(text, "A561V")

    def test_scan_caches_document_attribution_for_repeated_literal(self, monkeypatch):
        scanner = VariantScanner("KCNH2")
        calls = 0

        def attributed_elsewhere(_text: str, raw_text: str) -> set[str]:
            nonlocal calls
            calls += 1
            assert raw_text == "L187P"
            return {"KCNQ1"}

        monkeypatch.setattr(scanner, "_document_gene_attribution", attributed_elsewhere)
        result = scanner.scan("L187P was found. Later L187P was confirmed.")

        assert calls == 1
        assert "L187P" not in {variant.normalized for variant in result.variants}

    def test_scan_skips_document_attribution_for_accepted_normalized_duplicate(
        self, monkeypatch
    ):
        scanner = VariantScanner("KCNH2")
        calls = 0

        def target_or_unassigned(_text: str, _raw_text: str) -> set[str]:
            nonlocal calls
            calls += 1
            return set()

        monkeypatch.setattr(scanner, "_document_gene_attribution", target_or_unassigned)
        result = scanner.scan("p.Arg534Cys was confirmed as R534C.")

        assert calls == 1
        assert [
            variant.normalized
            for variant in result.variants
            if variant.normalized == "R534C"
        ] == ["R534C"]
