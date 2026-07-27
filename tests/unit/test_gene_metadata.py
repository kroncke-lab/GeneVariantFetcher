import sqlite3

import utils.gene_metadata as gene_metadata
from pipeline.steps import discover_synonyms
from utils.gene_metadata import (
    clear_gene_metadata_cache,
    get_gene_aliases,
    get_gene_metadata,
    lookup_variantfeatures_residue,
)


def _write_case_fixture_db(db_path):
    """Build a VariantFeatures-shaped slice covering both gene-symbol casings.

    ``XORF1`` exists only in mixed case, mirroring the HGNC lowercase ``orf``
    convention seen in real slices (``C19orf25``). ``COLL1`` is stored under both
    casings with deliberately different payloads so the exact-match-first
    ordering is observable.
    """

    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        CREATE TABLE genes (
            symbol TEXT,
            canonical_transcript TEXT,
            ncbi_id TEXT,
            ensembl_id TEXT
        );
        CREATE TABLE variant_consequences (
            gene_symbol TEXT,
            transcript_id TEXT,
            hgvs_p TEXT,
            hgvs_c TEXT,
            aa_pos INTEGER,
            aa_ref TEXT,
            aa_alt TEXT
        );
        CREATE TABLE transcripts (
            transcript_id TEXT,
            gene_symbol TEXT,
            refseq_match TEXT,
            protein_id TEXT,
            cds_length INTEGER,
            is_canonical INTEGER,
            is_mane_select INTEGER
        );
        INSERT INTO genes VALUES
            ('Xorf1', 'ENST000LOWER', '999', 'ENSG000LOWER'),
            ('COLL1', 'ENST000UPPER', '111', 'ENSG000UPPER');
        INSERT INTO transcripts VALUES
            ('ENST000LOWER', 'Xorf1', 'NM_000999.3', 'ENSP000LOWER', 303, 1, 1),
            ('ENST000UPPER', 'COLL1', 'NM_000111.2', 'ENSP000UPPER', 63, 1, 1);
        INSERT INTO variant_consequences VALUES
            ('Xorf1', 'ENST000LOWER', 'ENSP000LOWER:p.Cys10Arg', 'ENST000LOWER:c.28T>C', 10, 'C', 'R'),
            ('Xorf1', 'ENST000LOWER', 'ENSP000LOWER:p.Ter101Ter', NULL, 101, '*', '*'),
            ('COLL1', 'ENST000UPPER', 'ENSP000UPPER:p.Ala20Gly', NULL, 20, 'A', 'G'),
            ('Coll1', 'ENST000COLLLOWER', 'ENSP000COLLLOWER:p.Trp900Tyr', NULL, 900, 'W', 'Y');
        """
    )
    conn.commit()
    conn.close()


def test_builtin_gene_aliases_include_common_protein_names_and_typos():
    aliases = get_gene_aliases("MYBPC3", include_query_aliases=True)
    metadata = get_gene_metadata("MYBPC3")

    assert "MYPBC3" in aliases
    assert "cMyBP-C" in aliases
    assert "cardiac myosin-binding protein C" in aliases
    assert metadata.protein_length == 1274


def test_variantfeatures_metadata_and_residue_lookup(tmp_path, monkeypatch):
    db_path = tmp_path / "variants.db"
    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        CREATE TABLE genes (
            symbol TEXT,
            canonical_transcript TEXT,
            ncbi_id TEXT,
            ensembl_id TEXT
        );
        CREATE TABLE variant_consequences (
            gene_symbol TEXT,
            transcript_id TEXT,
            hgvs_p TEXT,
            hgvs_c TEXT,
            aa_pos INTEGER,
            aa_ref TEXT,
            aa_alt TEXT
        );
        INSERT INTO genes VALUES ('TEST1', 'ENST000TEST', '1234', 'ENSG000TEST');
        INSERT INTO variant_consequences VALUES
            ('TEST1', 'ENST000TEST', 'ENSP000TEST:p.Cys10Arg', 'ENST000TEST:c.28T>C', 10, 'C', 'R'),
            ('TEST1', 'ENST000TEST', 'ENSP000TEST:p.Ter101Ter', NULL, 101, '*', '*');
        """
    )
    conn.close()
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    clear_gene_metadata_cache()

    metadata = get_gene_metadata("TEST1")
    residue = lookup_variantfeatures_residue(
        "TEST1",
        position=10,
        protein_notation="p.Cys10Arg",
        cdna_notation="c.28T>C",
    )

    assert metadata.protein_length == 100
    assert metadata.canonical_transcript == "ENST000TEST"
    assert metadata.ensembl_id == "ENSG000TEST"
    assert residue is not None
    assert residue.reference_residues == ("C",)
    assert residue.alternate_residues == ("R",)
    assert residue.matched_hgvs_p is True
    assert residue.matched_hgvs_c is True

    clear_gene_metadata_cache()


def test_mixed_case_gene_symbol_still_resolves_via_upper_fallback(
    tmp_path, monkeypatch
):
    """A gene stored only in mixed case must stay reachable.

    The exact-match probe misses, and because ``SELECT MAX(aa_pos)`` returns a
    single all-NULL row rather than no rows, the fallback has to trigger on the
    NULL payload instead of on an empty result set.
    """

    db_path = tmp_path / "variants.db"
    _write_case_fixture_db(db_path)
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    clear_gene_metadata_cache()

    metadata = get_gene_metadata("XORF1")
    residue = lookup_variantfeatures_residue(
        "XORF1",
        position=10,
        protein_notation="p.Cys10Arg",
        cdna_notation="c.28T>C",
    )

    assert metadata.protein_length == 100
    assert metadata.canonical_transcript == "ENST000LOWER"
    assert metadata.ensembl_id == "ENSG000LOWER"
    assert metadata.ncbi_id == "999"
    # Proves VariantFeatures actually resolved; the alias assertion below cannot,
    # because base metadata already carries ('XORF1',) when the lookup returns None.
    assert "variantfeatures" in metadata.sources
    # _dedupe folds aliases case-insensitively, so the stored 'Xorf1' collapses
    # into the normalized symbol rather than appearing alongside it.
    assert metadata.aliases == ("XORF1",)
    # The transcripts query is gene-scoped too, so the fallback has to carry it.
    assert metadata.refseq_transcripts == ("NM_000999.3",)
    assert metadata.protein_ids == ("ENSP000LOWER",)
    assert residue is not None
    assert residue.reference_residues == ("C",)
    assert residue.alternate_residues == ("R",)
    assert residue.matched_hgvs_p is True
    assert residue.matched_hgvs_c is True

    clear_gene_metadata_cache()


def test_exact_case_match_wins_over_upper_fallback(tmp_path, monkeypatch):
    """When the exact casing matches, the mixed-case twin is not folded in.

    This is the deliberate narrowing that buys the index: the old
    ``UPPER(gene_symbol) = ?`` predicate unioned both casings, so it would have
    reported the ``Coll1`` row's aa_pos 900 here.
    """

    db_path = tmp_path / "variants.db"
    _write_case_fixture_db(db_path)
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    clear_gene_metadata_cache()

    metadata = get_gene_metadata("COLL1")

    assert metadata.protein_length == 19
    assert metadata.canonical_transcript == "ENST000UPPER"
    assert metadata.refseq_transcripts == ("NM_000111.2",)
    assert metadata.protein_ids == ("ENSP000UPPER",)

    clear_gene_metadata_cache()


def test_exact_case_lookup_does_not_emit_an_upper_scan(tmp_path, monkeypatch):
    """Pin the performance property, not just the result.

    Wrapping the column in ``UPPER()`` makes ``idx_consequences_gene``
    unusable and turns each lookup into a full index scan, so the hot path must
    not reach for that predicate when the exact form already matched.
    """

    db_path = tmp_path / "variants.db"
    _write_case_fixture_db(db_path)
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))

    statements: list[str] = []
    real_connect = gene_metadata._connect_readonly

    def traced_connect(path):
        conn = real_connect(path)
        conn.set_trace_callback(statements.append)
        return conn

    monkeypatch.setattr(gene_metadata, "_connect_readonly", traced_connect)

    # sqlite3 hands the trace callback the *expanded* SQL, so the bound symbol
    # appears as a literal. That lets each gene-scoped table be pinned positively
    # rather than relying on "some statement ran" — the four sqlite_master probes
    # from _has_table would satisfy a bare non-empty check on their own.
    def emitted(table: str, predicate: str) -> bool:
        return any(f"FROM {table}" in sql and predicate in sql for sql in statements)

    clear_gene_metadata_cache()
    statements.clear()
    get_gene_metadata("COLL1")
    assert emitted("genes", "symbol = 'COLL1'")
    assert emitted("transcripts", "gene_symbol = 'COLL1'")
    assert emitted("variant_consequences", "gene_symbol = 'COLL1'")
    assert not any("UPPER(" in sql for sql in statements)

    # The per-variant residue lookup is the hottest converted query — it has no
    # lru_cache and runs once per parsed protein position, so a regression here
    # is paid thousands of times per paper rather than once per gene.
    statements.clear()
    residue = lookup_variantfeatures_residue("COLL1", position=20)
    assert residue is not None
    assert residue.reference_residues == ("A",)
    assert emitted("variant_consequences", "gene_symbol = 'COLL1' AND aa_pos = 20")
    assert not any("UPPER(" in sql for sql in statements)

    # The mixed-case-only gene is the one case that still has to pay the scan.
    clear_gene_metadata_cache()
    statements.clear()
    get_gene_metadata("XORF1")
    assert any("UPPER(" in sql for sql in statements)

    clear_gene_metadata_cache()


def test_legacy_variants_schema_covers_both_casings(tmp_path, monkeypatch):
    """The legacy ``variants`` path was converted too, and had no test at all.

    It is reached only when ``variant_consequences`` is absent, so the modern
    fixture never exercises it. Both the metadata branch and
    ``_lookup_legacy_variantfeatures_residue`` bind the gene as the first
    parameter, which is what makes a WHERE-clause reorder fail silently here.
    """

    db_path = tmp_path / "variants.db"
    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        CREATE TABLE variants (
            gene TEXT,
            resnum INTEGER,
            uniprot_id TEXT,
            var TEXT,
            var_hgvs_p TEXT,
            var_hgvs_c TEXT,
            wt_aa TEXT,
            mut_aa TEXT
        );
        INSERT INTO variants VALUES
            ('LEG1', 10, 'P12345', 'C10R', 'p.Cys10Arg', 'c.28T>C', 'C', 'R'),
            ('LEG1', 250, 'P12345', 'A250T', 'p.Ala250Thr', 'c.748G>A', 'A', 'T'),
            ('Leg2', 20, 'Q67890', 'W20Y', 'p.Trp20Tyr', 'c.59G>A', 'W', 'Y'),
            ('Leg2', 400, 'Q67890', 'G400S', 'p.Gly400Ser', 'c.1198G>A', 'G', 'S');
        """
    )
    conn.commit()
    conn.close()
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    clear_gene_metadata_cache()

    # Exact casing: indexed predicate is enough.
    upper = get_gene_metadata("LEG1")
    assert upper.protein_length == 250
    assert upper.protein_ids == ("P12345",)
    residue = lookup_variantfeatures_residue(
        "LEG1", position=10, protein_notation="p.Cys10Arg", cdna_notation="c.28T>C"
    )
    assert residue is not None
    assert residue.reference_residues == ("C",)
    assert residue.alternate_residues == ("R",)
    assert residue.transcripts == ("P12345",)
    assert residue.matched_hgvs_p is True
    assert residue.matched_hgvs_c is True
    # Position is the second bound parameter; a transposed binding returns None.
    assert lookup_variantfeatures_residue("LEG1", position=999) is None

    # Mixed case only: the UPPER() fallback has to carry the legacy path too.
    clear_gene_metadata_cache()
    lower = get_gene_metadata("LEG2")
    assert lower.protein_length == 400
    assert lower.protein_ids == ("Q67890",)
    legacy_residue = lookup_variantfeatures_residue(
        "LEG2", position=20, protein_notation="p.Trp20Tyr"
    )
    assert legacy_residue is not None
    assert legacy_residue.reference_residues == ("W",)
    assert legacy_residue.matched_hgvs_p is True

    clear_gene_metadata_cache()


def test_discover_synonyms_keeps_builtin_aliases_when_ncbi_fails(monkeypatch):
    import gene_literature.synonym_finder as synonym_finder

    class FailingSynonymFinder:
        def __init__(self, **kwargs):
            pass

        def find_gene_synonyms(self, *args, **kwargs):
            raise RuntimeError("offline")

    monkeypatch.setattr(synonym_finder, "SynonymFinder", FailingSynonymFinder)

    result = discover_synonyms("MYBPC3", email="nobody@example.org")

    assert result.success is False
    assert "cMyBP-C" in result.data["synonyms"]
    assert "cardiac myosin-binding protein C" in result.data["synonyms"]
