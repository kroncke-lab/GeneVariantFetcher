import sqlite3

from pipeline.steps import discover_synonyms
from utils.gene_metadata import (
    clear_gene_metadata_cache,
    get_gene_aliases,
    get_gene_metadata,
    lookup_variantfeatures_residue,
)


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
