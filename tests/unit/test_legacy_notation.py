"""Source-only legacy indels stay useful without becoming invented HGVS."""

from __future__ import annotations

import sqlite3

from harvesting.migrate_to_sqlite import (
    create_database_schema,
    insert_variant_data,
    sanitize_variant_notation,
)
from pipeline.extraction import ExpertExtractor
from utils.legacy_notation import (
    normalize_legacy_notation,
    preserve_source_only_legacy_identity,
)
from utils.variant_normalizer import VariantNormalizer
from utils.variant_scanner import VariantScanner, merge_scanner_results


def test_strict_legacy_grammar_is_narrow_and_normalized():
    assert normalize_legacy_notation(" 4321DeLAC ") == "4321delAC"
    assert normalize_legacy_notation("7456insT") == "7456insT"
    assert normalize_legacy_notation("120delta") is None
    assert normalize_legacy_notation("4321delac") is None
    assert normalize_legacy_notation("42delAC") is None
    assert normalize_legacy_notation("4321del") is None


def test_cdna_normalizer_never_invents_prefix_for_legacy_identity():
    normalizer = VariantNormalizer("BRCA2")
    assert normalizer.normalize_cdna("4321delAC") is None
    assert normalizer.normalize_cdna("c.4321delAC") == "c.4321delAC"


def test_extraction_filter_reverses_fabricated_cdna_using_verbatim_source():
    extractor = ExpertExtractor(models=["test-model"])
    payload = {
        "variants": [
            {
                "gene_symbol": "BRCA2",
                "cdna_notation": "c.4321delAC",
                "protein_notation": None,
                "source_notation": "4321delAC",
            }
        ]
    }

    filtered = extractor._filter_extraction_artifacts(payload, "BRCA2")

    assert len(filtered["variants"]) == 1
    assert filtered["variants"][0]["cdna_notation"] is None
    assert filtered["variants"][0]["legacy_notation"] == "4321delAC"


def test_fabricated_brca_cdna_is_removed_when_real_protein_is_present():
    extractor = ExpertExtractor(models=["test-model"])
    payload = {
        "variants": [
            {
                "gene_symbol": "BRCA2",
                "cdna_notation": "c.4321delAC",
                "protein_notation": "p.Gln1441fs",
                "source_notation": "4321delAC",
            }
        ]
    }

    (variant,) = extractor._filter_extraction_artifacts(payload, "BRCA2")["variants"]

    assert variant["cdna_notation"] is None
    assert variant["legacy_notation"] == "4321delAC"
    assert variant["protein_notation"] == "p.Gln1441fs"


def test_non_brca_prefixless_indel_remains_cdna():
    extractor = ExpertExtractor(models=["test-model"])
    payload = {
        "variants": [
            {
                "gene_symbol": "BMPR2",
                "cdna_notation": "1234delA",
                "source_notation": "1234delA",
            }
        ]
    }

    (variant,) = extractor._filter_extraction_artifacts(payload, "BMPR2")["variants"]

    assert variant["cdna_notation"] == "c.1234delA"
    assert variant.get("legacy_notation") is None


def test_legacy_field_wins_before_scanner_merge_without_duplicate():
    source = "In BRCA2, 4321delAC was identified in a participant."
    payload = {
        "variants": [
            {
                "gene_symbol": "BRCA2",
                "cdna_notation": "c.4321delAC",
                "legacy_notation": "4321delAC",
                "source_notation": "c.4321delAC",
            }
        ]
    }

    merged = merge_scanner_results(
        payload,
        VariantScanner("BRCA2").scan(source),
        "BRCA2",
        document_text=source,
    )

    assert len(merged["variants"]) == 1
    assert merged["variants"][0]["cdna_notation"] is None
    assert merged["variants"][0]["legacy_notation"] == "4321delAC"


def test_migration_preserves_and_deduplicates_source_only_legacy(tmp_path):
    db = tmp_path / "legacy.db"
    conn = create_database_schema(str(db))
    conn.row_factory = sqlite3.Row
    conn.execute("INSERT INTO papers (pmid, title) VALUES ('1', 'Synthetic')")
    row = {
        "gene_symbol": "BRCA2",
        "cdna_notation": "c.4321delAC",
        "source_notation": "4321delAC",
        "patients": {},
        "penetrance_data": {},
    }

    first = insert_variant_data(conn.cursor(), "1", dict(row))
    second = insert_variant_data(conn.cursor(), "1", dict(row))
    conn.commit()

    assert first == second
    stored = conn.execute(
        "SELECT cdna_notation, legacy_notation FROM variants"
    ).fetchone()
    assert stored["cdna_notation"] is None
    assert stored["legacy_notation"] == "4321delAC"
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_migration_reuses_proven_old_fabricated_cdna_row(tmp_path):
    db = tmp_path / "legacy-backcompat.db"
    conn = create_database_schema(str(db))
    conn.row_factory = sqlite3.Row
    conn.execute("INSERT INTO papers (pmid, title) VALUES ('1', 'Synthetic')")
    cursor = conn.execute(
        """
        INSERT INTO variants (gene_symbol, cdna_notation)
        VALUES ('BRCA2', 'c.4321delAC')
        """
    )
    old_id = cursor.lastrowid
    conn.execute(
        """
        INSERT INTO variant_papers (variant_id, pmid, source_notation)
        VALUES (?, '1', ?)
        """,
        (old_id, "4321\t delAC\r\n"),
    )

    reused = insert_variant_data(
        conn.cursor(),
        "1",
        {
            "gene_symbol": "BRCA2",
            "legacy_notation": "4321delAC",
            "source_notation": "4321delAC",
            "patients": {},
            "penetrance_data": {},
        },
    )
    conn.commit()

    assert reused == old_id
    stored = conn.execute(
        "SELECT cdna_notation, legacy_notation FROM variants WHERE variant_id = ?",
        (old_id,),
    ).fetchone()
    assert tuple(stored) == (None, "4321delAC")
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_sanitizer_rejects_loose_source_text_but_keeps_strict_legacy():
    valid = {"source_notation": "7456insT"}
    invalid = {"source_notation": "approximately 7456 inserted bases"}
    assert sanitize_variant_notation(valid) is True
    assert valid["legacy_notation"] == "7456insT"
    assert sanitize_variant_notation(invalid) is False


def test_scanner_emits_legacy_identity_and_verbatim_source():
    result = VariantScanner(gene_symbol="BRCA2").scan(
        "In BRCA2, 7456insT was identified in a participant."
    )
    legacy = next(v for v in result.variants if v.source == "bic_prefixless")
    assert legacy.normalized == "7456insT"
    assert legacy.notation_type == "legacy"
    (row,) = [
        item for item in result.to_variant_dicts("BRCA2") if item.get("legacy_notation")
    ]
    assert row["cdna_notation"] is None
    assert row["legacy_notation"] == "7456insT"
    assert row["source_notation"] == "7456insT"


def test_legacy_grammar_accepts_bic_deleted_base_counts():
    """BIC writes the deleted bases or their count; both are the same class.

    ``4184del4`` is standard BRCA2 BIC notation. The extraction prompt now
    routes bare labels into ``legacy_notation`` instead of ``cdna_notation``,
    so a grammar that accepted only the base form silently destroyed these
    identities rather than merely declining to normalize them.
    """

    assert normalize_legacy_notation("4184del4") == "4184del4"
    assert normalize_legacy_notation("1294del40") == "1294del40"
    assert normalize_legacy_notation("300dup12") == "300dup12"
    assert normalize_legacy_notation("5213_5216del4") == "5213_5216del4"
    assert normalize_legacy_notation("1773-1776del4") == "1773-1776del4"
    assert normalize_legacy_notation("7734_7740del6ins9") == "7734_7740del6ins9"
    # Still narrow: no payload at all, and over-long counts, stay out.
    assert normalize_legacy_notation("4321del") is None
    assert normalize_legacy_notation("4321del1234") is None


def test_prefixed_bic_count_form_is_demoted_instead_of_dropped():
    extractor = ExpertExtractor(models=["test-model"])
    payload = {
        "variants": [
            {
                "gene_symbol": "BRCA2",
                "cdna_notation": "c.4184del4",
                "source_notation": "c.4184del4",
            }
        ]
    }

    (variant,) = extractor._filter_extraction_artifacts(payload, "BRCA2")["variants"]

    assert variant["cdna_notation"] is None
    assert variant["legacy_notation"] == "4184del4"


def test_prefixed_bic_range_count_forms_are_demoted_instead_of_dropped():
    extractor = ExpertExtractor(models=["test-model"])
    payload = {
        "variants": [
            {
                "gene_symbol": "BRCA2",
                "cdna_notation": notation,
                "source_notation": notation,
            }
            for notation in (
                "c.5213_5216del4",
                "c.1773-1776del4",
                "c.7734_7740del6ins9",
                "c.8374_8384del11ins3",
                "c.9063_9078del16",
            )
        ]
    }

    variants = extractor._filter_extraction_artifacts(payload, "BRCA2")["variants"]

    assert [variant["legacy_notation"] for variant in variants] == [
        "5213_5216del4",
        "1773-1776del4",
        "7734_7740del6ins9",
        "8374_8384del11ins3",
        "9063_9078del16",
    ]
    assert all(variant["cdna_notation"] is None for variant in variants)


def test_intronic_range_with_source_bases_remains_cdna():
    extractor = ExpertExtractor(models=["test-model"])
    payload = {
        "variants": [
            {
                "gene_symbol": "BRCA2",
                "cdna_notation": "c.7436_7437-2delAGAT",
                "source_notation": "c.7436_7437-2delAGAT",
            }
        ]
    }

    (variant,) = extractor._filter_extraction_artifacts(payload, "BRCA2")["variants"]

    assert variant["cdna_notation"] == "c.7436_7437-2delAGAT"
    assert variant.get("legacy_notation") is None


def test_scanner_routes_count_form_and_non_brca_prefixless_indel_by_gene():
    brca = VariantScanner("BRCA2").scan(
        "In BRCA2, 4184del4 was identified in a participant."
    )
    bmpr2 = VariantScanner("BMPR2").scan(
        "In BMPR2, 1234delA was identified in a participant."
    )

    brca_row = next(v for v in brca.variants if v.raw_text == "4184del4")
    bmpr2_row = next(v for v in bmpr2.variants if v.raw_text == "1234delA")
    assert (brca_row.notation_type, brca_row.normalized) == ("legacy", "4184del4")
    assert (bmpr2_row.notation_type, bmpr2_row.normalized) == ("cdna", "c.1234delA")


def test_source_proven_legacy_survives_a_junk_cdna_value():
    """A transcript accession in the cDNA slot must not delete the variant.

    ``_filter_extraction_artifacts`` nulls malformed cDNA *after* legacy
    preservation runs, so discarding a source-proven legacy label merely
    because some cDNA string is present left the row with no identity at all
    and the nameless-row guard dropped it.
    """

    variant = {
        "cdna_notation": "NM_000059.3",
        "legacy_notation": "4321delAC",
        "source_notation": "4321delAC",
    }
    assert preserve_source_only_legacy_identity(variant) == "4321delAC"
    assert variant["legacy_notation"] == "4321delAC"

    # An alias that the source does NOT independently print is still dropped
    # when a real HGVS identity is present.
    genuine = {
        "cdna_notation": "c.68_69delAG",
        "legacy_notation": "185delAG",
        "source_notation": "c.68_69delAG",
    }
    assert preserve_source_only_legacy_identity(genuine) is None
    assert genuine["cdna_notation"] == "c.68_69delAG"
    assert genuine["legacy_notation"] is None


def test_scanner_legacy_does_not_duplicate_a_fabricated_cdna_row():
    """The two lanes must not emit the same mutation twice.

    Deterministic table lanes still prefix a bare BIC label with ``c.`` when
    they cannot see the verbatim cell, and that row carries the counts. A
    scanner legacy identity for the same label would otherwise be appended as
    a second, count-less row for one physical mutation.
    """

    text = "In BRCA2, 4321delAC was identified in seven carriers of our cohort."
    payload = {
        "variants": [
            {
                "gene_symbol": "BRCA2",
                "cdna_notation": "c.4321delAC",
                "patients": {"total_carriers": 7},
            }
        ],
        "extraction_metadata": {},
    }
    merged = merge_scanner_results(
        payload, VariantScanner("BRCA2").scan(text), "BRCA2", document_text=text
    )
    assert len(merged["variants"]) == 1
    assert merged["variants"][0]["cdna_notation"] == "c.4321delAC"

    # An unrelated existing row must not suppress a genuine legacy identity.
    other = {
        "variants": [{"gene_symbol": "BRCA2", "cdna_notation": "c.9999A>G"}],
        "extraction_metadata": {},
    }
    assert (
        len(
            merge_scanner_results(
                other,
                VariantScanner("BRCA2").scan(text),
                "BRCA2",
                document_text=text,
            )["variants"]
        )
        == 2
    )
