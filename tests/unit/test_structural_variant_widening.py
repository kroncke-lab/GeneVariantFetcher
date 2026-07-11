"""Stage 5 B: structural / splice / delins widening unit coverage."""

from __future__ import annotations

import csv
from pathlib import Path

from cli.compare_variants import to_canonical_form
from harvesting.migrate_to_sqlite import (
    create_database_schema,
    get_or_create_variant,
    sanitize_variant_notation,
)
from pipeline.extraction import ExpertExtractor
from utils.variant_normalizer import (
    VariantNormalizer,
    structural_variant_identity,
)
from utils.variant_scanner import VariantScanner, merge_scanner_results


STRUCTURAL_FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "benchmarks"
    / "curated_extraction_eval"
    / "fixtures"
    / "structural_widening_cases.csv"
)


def test_filter_keeps_structural_only_variant():
    extractor = ExpertExtractor(models=["gpt-4"])
    data = {
        "extraction_metadata": {},
        "variants": [
            {
                "gene_symbol": "BRCA1",
                "variant_class": "exon_deletion",
                "structural_description": "deletion of exons 3-5",
            },
            {
                "gene_symbol": "BRCA1",
                "protein_notation": "0.99",  # still dropped
            },
        ],
    }
    filtered = extractor._filter_extraction_artifacts(data, "BRCA1")
    assert len(filtered["variants"]) == 1
    kept = filtered["variants"][0]
    assert kept["variant_class"] == "exon_deletion"
    assert "exons 3-5" in kept["structural_description"]


def test_sanitize_keeps_structural_identity():
    v = {
        "gene_symbol": "LDLR",
        "variant_class": "large_deletion",
        "structural_description": "whole-gene deletion",
    }
    assert sanitize_variant_notation(v) is True


def test_notation_free_identity_requires_genuine_structural_class():
    extractor = ExpertExtractor(models=["gpt-4"])
    rows = [
        {
            "gene_symbol": "BRCA1",
            "variant_class": "missense",
            "structural_description": "a described missense change",
        },
        {
            "gene_symbol": "BRCA1",
            "variant_class": "other",
            "structural_description": "N-terminal mutation group",
        },
        {
            "gene_symbol": "BRCA1",
            "variant_class": "exon_duplication",
            "structural_description": "duplication of exons 8-10",
        },
    ]

    filtered = extractor._filter_extraction_artifacts(
        {"extraction_metadata": {}, "variants": rows}, "BRCA1"
    )

    assert filtered["variants"] == [rows[2]]
    assert sanitize_variant_notation(dict(rows[0])) is False
    assert sanitize_variant_notation(dict(rows[1])) is False
    assert sanitize_variant_notation(dict(rows[2])) is True


def test_protein_validation_fullmatches_and_supports_three_letter_delins():
    extractor = ExpertExtractor(models=["gpt-4"])
    valid = {
        "gene_symbol": "KCNH2",
        "protein_notation": "p.Arg176delinsLysGly",
    }
    trailing_junk = {
        "gene_symbol": "KCNH2",
        "protein_notation": "p.Arg176Trp-not-a-variant",
    }

    filtered = extractor._filter_extraction_artifacts(
        {
            "extraction_metadata": {},
            "variants": [valid, trailing_junk],
        },
        "KCNH2",
    )

    assert filtered["variants"] == [valid]
    assert sanitize_variant_notation(dict(valid)) is True
    invalid = dict(trailing_junk)
    assert sanitize_variant_notation(invalid) is False
    assert invalid["protein_notation"] is None


def test_to_canonical_form_ivs_delins_structural():
    assert to_canonical_form("IVS9+1G>A") == "IVS9+1G>A"
    assert to_canonical_form("ivs5-2a>g") == "IVS5-2A>G"
    assert to_canonical_form("c.4513_4521delAAGCAGinsT") == "c.4513_4521delinsT"
    assert to_canonical_form("del:exon3-5") == "del:exon3-5"
    assert to_canonical_form("deletion of exons 3-5") == "del:exon3-5"
    assert to_canonical_form("whole-gene deletion") == "del:wholegene"


def test_scanner_finds_delins_and_exon_events():
    text = (
        "In BRCA1 we identified c.123_125delAGinsT and a deletion of exons 3-5. "
        "A BRCA1 whole-gene deletion was also reported."
    )
    scanner = VariantScanner(gene_symbol="BRCA1")
    result = scanner.scan(text)
    norms = {v.normalized for v in result.variants}
    sources = {v.source for v in result.variants}
    assert any("del" in n and "ins" in n for n in norms) or any(
        "delins" in (v.variant_class or "") or v.source == "cdna_delins"
        for v in result.variants
    )
    assert "exon_event" in sources or any(
        v.variant_class == "exon_deletion" for v in result.variants
    )
    assert any(v.normalized == "del:wholegene" for v in result.variants)


def test_scanner_bic_prefixless_requires_gene_context():
    scanner = VariantScanner(gene_symbol="BRCA1")
    # Without gene mention nearby, bare lab number should not match
    no_gene = scanner.scan("The assay used 185delAG reagent lot 12.")
    with_gene = scanner.scan("BRCA1 185delAG was found in the kindred.")
    assert not any(v.source == "bic_prefixless" for v in no_gene.variants)
    assert any(v.source == "bic_prefixless" for v in with_gene.variants)


def test_scanner_bic_ignores_lowercase_delta_word():
    """'120delta' (del + lowercase 'ta') must not match as c.120delTA; uppercase
    BIC bases like 185delAG still do."""
    scanner = VariantScanner(gene_symbol="SCN5A")
    delta = scanner.scan("The SCN5A ECG showed a 120delta wave morphology.")
    assert not any(v.source == "bic_prefixless" for v in delta.variants)
    real = scanner.scan("SCN5A 185delAG segregated in the family.")
    assert any(v.source == "bic_prefixless" for v in real.variants)


def test_scanner_skips_negated_structural_events():
    """Negation before or after an event must suppress the structural hint."""
    scanner = VariantScanner(gene_symbol="BRCA1")
    negated_before = scanner.scan(
        "In BRCA1 there was no deletion of exons 3-5 on MLPA."
    )
    negated_after = scanner.scan(
        "In BRCA1 deletion of exons 3-5 was not detected by MLPA."
    )
    negated_whole_gene = scanner.scan(
        "A BRCA1 whole-gene deletion was not identified in the proband."
    )
    positive = scanner.scan("In BRCA1 a deletion of exons 3-5 was identified.")
    assert not any(v.source == "exon_event" for v in negated_before.variants)
    assert not any(v.source == "exon_event" for v in negated_after.variants)
    assert not any(v.source == "whole_gene_del" for v in negated_whole_gene.variants)
    assert any(v.source == "exon_event" for v in positive.variants)


def test_scanner_requires_target_gene_and_positive_structural_finding():
    scanner = VariantScanner(gene_symbol="BRCA1")
    no_gene = scanner.scan("A deletion of exons 3-5 was identified by MLPA.")
    wrong_gene = scanner.scan("A BRCA2 deletion of exons 3-5 was identified.")
    methods_only = scanner.scan(
        "BRCA1 deletion of exons 3-5 was included in the assay design."
    )
    assay_capability = scanner.scan(
        "The MLPA assay was designed to detect a BRCA1 deletion of exons 3-5."
    )
    positive = scanner.scan(
        "BRCA1 deletion of exons 3-5 was confirmed in the affected proband."
    )

    assert not any(v.source == "exon_event" for v in no_gene.variants)
    assert not any(v.source == "exon_event" for v in wrong_gene.variants)
    assert not any(v.source == "exon_event" for v in methods_only.variants)
    assert not any(v.source == "exon_event" for v in assay_capability.variants)
    assert any(v.source == "exon_event" for v in positive.variants)


def test_structural_identity_drives_scanner_table_and_continuation_dedup():
    extractor = ExpertExtractor(models=["gpt-4"])
    exon = {
        "gene_symbol": "BRCA1",
        "variant_class": "exon_deletion",
        "structural_description": "Deletion of exons 3 – 5",
    }
    exon_alias = {
        "gene_symbol": "BRCA1",
        "variant_class": "exon_deletion",
        "structural_description": "del:exon3-5",
    }
    whole_gene = {
        "gene_symbol": "BRCA1",
        "variant_class": "large_deletion",
        "structural_description": "whole-gene deletion",
    }

    scanned = VariantScanner("BRCA1").scan(
        "BRCA1 deletion of exons 3-5 was identified in the proband."
    )
    merged_scanner = merge_scanner_results(
        {"variants": [exon.copy()], "extraction_metadata": {}},
        scanned,
        "BRCA1",
    )
    assert len(merged_scanner["variants"]) == 1

    deduped = extractor._dedupe_table_variants([exon, exon_alias, whole_gene])
    assert len(deduped) == 2
    assert {
        structural_variant_identity(v["structural_description"]) for v in deduped
    } == {"del:exon3-5", "del:wholegene"}

    merged_table = extractor._merge_table_variants(
        {"variants": [exon.copy()], "extraction_metadata": {}},
        [exon_alias, whole_gene],
    )
    assert len(merged_table["variants"]) == 2

    merged_continuation = extractor._merge_continuation_results(
        {"variants": [exon.copy()], "extraction_metadata": {}},
        {"variants": [exon_alias, whole_gene], "extraction_metadata": {}},
    )
    assert len(merged_continuation["variants"]) == 2


def test_scanner_variant_dict_preserves_structural_and_splice_identity():
    result = VariantScanner("BRCA1").scan(
        "BRCA1 deletion of exons 3-5 was identified together with IVS9+1G>A."
    )
    rows = result.to_variant_dicts("BRCA1")

    structural = next(row for row in rows if row["variant_class"] == "exon_deletion")
    assert structural["protein_notation"] is None
    assert structural["cdna_notation"] is None
    assert structural["structural_description"] == "deletion of exons 3-5"

    splice = next(row for row in rows if row["variant_class"] == "splice")
    assert splice["cdna_notation"] == "IVS9+1G>A"


def test_structural_identity_is_stable_in_normalizer_and_sqlite(tmp_path):
    normalizer = VariantNormalizer("BRCA1")
    canonical = normalizer.get_canonical_form(
        {"structural_description": "Deletion of exons 3 – 5"}
    )
    assert canonical["canonical_key"] == "BRCA1:del:exon3-5"

    conn = create_database_schema(str(tmp_path / "structural.db"))
    cursor = conn.cursor()
    first = get_or_create_variant(
        cursor,
        {
            "gene_symbol": "BRCA1",
            "variant_class": "exon_deletion",
            "structural_description": "Deletion of exons 3 – 5",
        },
    )
    replay = get_or_create_variant(
        cursor,
        {
            "gene_symbol": "BRCA1",
            "variant_class": "exon_deletion",
            "structural_description": "del:exon3-5",
        },
    )
    assert replay == first
    assert cursor.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_synthetic_structural_fixture_is_exercised_deterministically():
    with STRUCTURAL_FIXTURE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 4
    for row in rows:
        assert int(row["affected"]) + int(row["unaffected"]) == int(row["carriers"])
        result = VariantScanner(row["gene"]).scan(
            f"{row['gene']} analysis identified {row['mention']} "
            "in an affected carrier."
        )
        matches = [
            variant
            for variant in result.variants
            if (to_canonical_form(variant.normalized) or variant.normalized)
            == row["expected_canonical"]
        ]
        assert matches, row["case_id"]
        assert any(v.variant_class == row["variant_class"] for v in matches)


def test_extract_sqlite_data_surfaces_structural_description(tmp_path):
    """A structural-only variant (point-form NULL) must reach the scorer's
    variant key via structural_description, so it can canonicalize and match."""
    import sqlite3

    from cli.compare_variants import extract_sqlite_data, introspect_sqlite

    db = tmp_path / "s.db"
    conn = sqlite3.connect(db)
    conn.executescript(
        """
        CREATE TABLE variants (variant_id INTEGER PRIMARY KEY, gene_symbol TEXT,
            protein_notation TEXT, cdna_notation TEXT, genomic_position TEXT,
            structural_description TEXT);
        CREATE TABLE variant_papers (variant_id INTEGER, pmid TEXT);
        CREATE TABLE penetrance_data (penetrance_id INTEGER PRIMARY KEY,
            variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
            affected_count INTEGER, unaffected_count INTEGER, uncertain_count INTEGER);
        """
    )
    conn.execute(
        "INSERT INTO variants VALUES (1,'BRCA1',NULL,NULL,NULL,'deletion of exons 3-5')"
    )
    conn.execute("INSERT INTO variant_papers VALUES (1, '111')")
    conn.execute("INSERT INTO penetrance_data VALUES (1,1,'111',4,4,0,0)")
    conn.commit()
    table_info = introspect_sqlite(conn)
    df = extract_sqlite_data(conn, table_info)
    conn.close()

    assert "deletion of exons 3-5" in set(df["variant"])
    assert to_canonical_form("deletion of exons 3-5") == "del:exon3-5"
