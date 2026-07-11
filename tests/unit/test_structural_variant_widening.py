"""Stage 5 B: structural / splice / delins widening unit coverage."""

from __future__ import annotations

from cli.compare_variants import to_canonical_form
from harvesting.migrate_to_sqlite import sanitize_variant_notation
from pipeline.extraction import ExpertExtractor
from utils.variant_scanner import VariantScanner


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
        "A whole-gene deletion was also reported."
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
    """A ruled-out exon deletion must not be emitted as a positive finding."""
    scanner = VariantScanner(gene_symbol="BRCA1")
    negated = scanner.scan("In BRCA1 there was no deletion of exons 3-5 on MLPA.")
    positive = scanner.scan("In BRCA1 a deletion of exons 3-5 was identified.")
    assert not any(v.source == "exon_event" for v in negated.variants)
    assert any(v.source == "exon_event" for v in positive.variants)


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
