import sqlite3

import pytest

import utils.gene_metadata as gene_metadata
from pipeline.steps import discover_synonyms
from utils.gene_metadata import (
    clear_gene_metadata_cache,
    get_gene_aliases,
    get_gene_metadata,
    lookup_variantfeatures_residue,
    resolve_variantfeatures_gene_symbols,
)


def _write_case_fixture_db(db_path):
    """Build a VariantFeatures-shaped slice covering both gene-symbol casings.

    ``XORF1`` exists only in mixed case, mirroring the HGNC lowercase ``orf``
    convention seen in real slices (``C19orf25``). ``COLL1`` is stored under both
    casings with deliberately different payloads so the exact-match-first
    ordering is observable. ``YORF2`` is stored under *two* spellings, neither of
    them the uppercase form, which is the only shape that forces the stored
    spellings to be unioned.

    The indexes include the composite residue key required by production
    lookups, plus ``idx_transcripts_gene`` and ``genes.symbol UNIQUE``. Without them the
    stored-spelling census takes its unindexed ``SELECT DISTINCT`` branch, so the
    loose index scan that production actually uses would run in zero tests.
    """

    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        CREATE TABLE genes (
            symbol TEXT UNIQUE,
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
        CREATE INDEX idx_consequences_gene_pos
            ON variant_consequences(gene_symbol, aa_pos);
        CREATE INDEX idx_transcripts_gene ON transcripts(gene_symbol);
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
            ('Coll1', 'ENST000COLLLOWER', 'ENSP000COLLLOWER:p.Trp900Tyr', NULL, 900, 'W', 'Y'),
            ('Yorf2', 'ENST000YLOW', 'ENSP000YLOW:p.Ser50Pro', NULL, 50, 'S', 'P'),
            ('yorf2', 'ENST000YMIX', 'ENSP000YMIX:p.Gly300Asp', NULL, 300, 'G', 'D');
        """
    )
    conn.commit()
    conn.close()


def _trace_statements(monkeypatch) -> list[str]:
    """Capture every SQL statement the module executes.

    sqlite3 hands the trace callback the *expanded* SQL, so bound symbols appear
    as literals and each predicate can be pinned positively.
    """

    statements: list[str] = []
    real_connect = gene_metadata._connect_readonly

    def traced_connect(path):
        conn = real_connect(path)
        conn.set_trace_callback(statements.append)
        return conn

    monkeypatch.setattr(gene_metadata, "_connect_readonly", traced_connect)
    return statements


def _census_statements(statements, table: str, column: str) -> list[str]:
    """The census passes over one column -- the walk or the DISTINCT fallback.

    Matched on both table and column: the module's own ``SELECT DISTINCT
    transcript_id`` query is not a census pass, and ``transcripts`` carries a
    ``gene_symbol`` column of its own that would otherwise be double-counted.
    """

    return [
        sql
        for sql in statements
        if (
            f"MIN({column}) FROM {table}" in sql
            and "RECURSIVE" in sql
            and f"MIN({column}) FROM {table} WHERE" in sql
        )
        or f"SELECT DISTINCT {column} FROM {table}" in sql
    ]


def test_builtin_gene_aliases_include_common_protein_names_and_typos():
    aliases = get_gene_aliases("MYBPC3", include_query_aliases=True)
    metadata = get_gene_metadata("MYBPC3")

    assert "MYPBC3" in aliases
    assert "cMyBP-C" in aliases
    assert "cardiac myosin-binding protein C" in aliases
    assert metadata.protein_length == 1274


def test_bmpr2_registration_turns_position_validation_back_on():
    """An unregistered gene validates nothing, and says so silently.

    ``validate_position`` returns True whenever ``protein_length`` is None and
    ``variant_scanner`` falls back to a 9999 ceiling, so before BMPR2 was
    registered every parsed position was accepted -- including positions carried
    in from a co-tabulated gene. 1038 is NP_001195.2, the protein of the MANE
    Select transcript NM_001204.7.
    """

    from utils.variant_normalizer import PROTEIN_LENGTHS, VariantNormalizer

    metadata = get_gene_metadata("BMPR2")
    assert metadata.protein_length == 1038
    assert PROTEIN_LENGTHS["BMPR2"] == 1038
    assert metadata.canonical_transcript == "NM_001204.7"

    normalizer = VariantNormalizer("BMPR2")
    assert normalizer.validate_position(1038) is True
    assert normalizer.validate_position(1039) is False
    # A RYR2-scale position is the shape a multi-gene PAH panel table leaks in.
    assert normalizer.validate_position(2870) is False


def test_bmpr2_aliases_do_not_match_the_other_pah_panel_genes():
    """BMPR2 is tabulated next to the genes it must not absorb.

    Heritable-PAH papers list BMPR2 alongside ACVRL1/ALK1, ENG, SMAD9 and CAV1,
    so an over-broad alias here reads a neighbouring row as a BMPR2 mention. The
    historical locus names PPH1/POVD1 still have to retrieve founder-era papers,
    which is why they are query-only rather than gene-mention aliases.
    """

    mention = gene_metadata.gene_alias_regex("BMPR2", include_query_aliases=False)
    query = gene_metadata.gene_alias_regex("BMPR2", include_query_aliases=True)

    for text in ("BMPR2 mutation carriers", "BMPR-II protein", "BMPR3"):
        assert mention.search(text), text

    for text in (
        "ACVRL1 (ALK1) variant",
        "ENG mutation",
        "SMAD9 carriers",
        "CAV1 c.474delA",
        "anaplastic lymphoma kinase (ALK)",
    ):
        assert not mention.search(text), text
        assert not query.search(text), text

    # Query-only: retrievable as literature terms, never a gene mention.
    assert query.search("linkage to the PPH1 locus")
    assert not mention.search("linkage to the PPH1 locus")


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
        CREATE INDEX idx_consequences_gene_pos
            ON variant_consequences(gene_symbol, aa_pos);
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


def test_mixed_case_gene_symbol_still_resolves_via_stored_spelling(
    tmp_path, monkeypatch
):
    """A gene stored only in mixed case must stay reachable.

    The exact-match probe misses, and because ``SELECT MAX(aa_pos)`` returns a
    single all-NULL row rather than no rows, the recovery has to trigger on the
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


def test_builtin_length_and_transcript_avoid_consequence_table_scan(
    tmp_path, monkeypatch
):
    db_path = tmp_path / "variants.db"
    _write_case_fixture_db(db_path)
    conn = sqlite3.connect(db_path)
    conn.execute(
        "INSERT INTO genes VALUES ('BRCA2', 'NM_000059.4', '675', 'ENSG00000139618')"
    )
    conn.execute(
        "INSERT INTO transcripts VALUES "
        "('NM_000059.4', 'BRCA2', 'NM_000059.4', 'NP_000050.3', 10257, 1, 1)"
    )
    conn.execute(
        "INSERT INTO variant_consequences VALUES "
        "('BRCA2', 'NM_000059.4', 'NP_000050.3:p.Ter3419Ter', NULL, 3419, '*', '*')"
    )
    conn.commit()
    conn.close()
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    statements = _trace_statements(monkeypatch)

    clear_gene_metadata_cache()
    metadata = get_gene_metadata("BRCA2")

    assert metadata.protein_length == 3418
    assert metadata.canonical_transcript == "NM_000059.4"
    assert not any(
        "SELECT MAX(aa_pos)" in sql or "SELECT DISTINCT transcript_id" in sql
        for sql in statements
    )
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
    statements = _trace_statements(monkeypatch)

    # Each gene-scoped table is pinned positively rather than by "some statement
    # ran" — the four sqlite_master probes from _has_table would satisfy a bare
    # non-empty check on their own.
    def emitted(table: str, predicate: str) -> bool:
        return any(f"FROM {table}" in sql and predicate in sql for sql in statements)

    clear_gene_metadata_cache()
    statements.clear()
    get_gene_metadata("COLL1")
    assert emitted("genes", "symbol = 'COLL1'")
    assert emitted("transcripts", "gene_symbol = 'COLL1'")
    assert emitted("variant_consequences", "gene_symbol = 'COLL1'")
    assert not any("UPPER(" in sql for sql in statements)
    # The census is lazy: every exact probe matched, so nothing should have
    # enumerated the stored spellings.
    assert _census_statements(statements, "variant_consequences", "gene_symbol") == []
    assert _census_statements(statements, "genes", "symbol") == []

    # The per-variant residue lookup is the hottest converted query — it runs
    # once per parsed protein position, so a regression here is paid thousands of
    # times per paper rather than once per gene.
    statements.clear()
    residue = lookup_variantfeatures_residue("COLL1", position=20)
    assert residue is not None
    assert residue.reference_residues == ("A",)
    assert emitted(
        "variant_consequences", "gene_symbol = 'COLL1' AND aa_pos IS NOT NULL"
    )
    assert not any("UPPER(" in sql for sql in statements)
    assert _census_statements(statements, "variant_consequences", "gene_symbol") == []

    # The mixed-case-only gene used to be the one case that still paid a full
    # scan. It now resolves through the census by querying the stored spelling
    # exactly, so no UPPER() predicate is emitted for it either.
    clear_gene_metadata_cache()
    statements.clear()
    get_gene_metadata("XORF1")
    assert not any("UPPER(" in sql for sql in statements)
    assert emitted("variant_consequences", "gene_symbol IN ('Xorf1')")
    assert emitted("genes", "symbol IN ('Xorf1')")
    assert emitted("transcripts", "gene_symbol IN ('Xorf1')")
    # Indexed columns take the loose index scan, not a whole-column walk.
    assert all(
        "RECURSIVE" in sql
        for sql in _census_statements(statements, "variant_consequences", "gene_symbol")
    )

    clear_gene_metadata_cache()


def test_absent_gene_resolves_from_the_census_without_any_scan(tmp_path, monkeypatch):
    """The half of the index fix that 88d486c deliberately left open.

    A gene that is absent misses the exact probe, so the old code fell back to
    ``UPPER(<column>) = ?`` and paid a full index scan — 8.5-24s per call on the
    real 38 GB slice, twice per ``get_gene_metadata`` and once per protein
    position in ``lookup_variantfeatures_residue``. A novel gene is by definition
    absent, which makes this the case a new-gene run pays on every lookup.
    """

    db_path = tmp_path / "variants.db"
    _write_case_fixture_db(db_path)
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    statements = _trace_statements(monkeypatch)
    clear_gene_metadata_cache()
    statements.clear()

    metadata = get_gene_metadata("ZZZZ9")

    # Absence still reports absence: base metadata only, no VariantFeatures.
    assert metadata.sources == ()
    assert metadata.protein_length is None
    assert metadata.canonical_transcript is None

    assert not any("UPPER(" in sql for sql in statements)
    # Positive proof the verdict came from the census rather than from a lucky
    # exact miss: the census ran, and it ran once per column despite two
    # gene-scoped queries missing on variant_consequences.
    assert (
        len(_census_statements(statements, "variant_consequences", "gene_symbol")) == 1
    )
    # Only the two exact probes touched the data table; no second query per miss.
    data_reads = [
        sql
        for sql in statements
        if "FROM variant_consequences" in sql
        and sql
        not in _census_statements(statements, "variant_consequences", "gene_symbol")
    ]
    assert len(data_reads) == 2
    assert all("gene_symbol = 'ZZZZ9'" in sql for sql in data_reads)

    # Same for the per-position lookup, and the census is already warm.
    statements.clear()
    assert lookup_variantfeatures_residue("ZZZZ9", position=20) is None
    assert not any("UPPER(" in sql for sql in statements)
    assert _census_statements(statements, "variant_consequences", "gene_symbol") == []

    clear_gene_metadata_cache()


def test_every_stored_casing_is_unioned_when_the_exact_spelling_has_no_payload(
    tmp_path, monkeypatch
):
    """``IN (<spellings>)`` has to select exactly what ``UPPER() = ?`` did.

    ``YORF2`` is stored as both ``Yorf2`` and ``yorf2``, so the exact probe finds
    nothing and the census owns the whole answer. Querying only the first stored
    spelling would report aa_pos 50 instead of 300 — a silent metadata loss that
    the two casing tests above cannot catch, because there each gene has a single
    non-uppercase spelling.
    """

    db_path = tmp_path / "variants.db"
    _write_case_fixture_db(db_path)
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    statements = _trace_statements(monkeypatch)
    clear_gene_metadata_cache()
    statements.clear()

    metadata = get_gene_metadata("YORF2")

    assert metadata.protein_length == 299
    assert "variantfeatures" in metadata.sources
    assert not any("UPPER(" in sql for sql in statements)
    unions = [
        sql
        for sql in statements
        if "FROM variant_consequences" in sql and "gene_symbol IN (" in sql
    ]
    assert unions
    assert all("'Yorf2'" in sql and "'yorf2'" in sql for sql in unions)

    clear_gene_metadata_cache()


def test_residue_lookup_is_memoized_per_gene_range(tmp_path, monkeypatch):
    """``target_gene_specificity`` calls this once per parsed protein position.

    Distinct positions recur across a paper's variants and tables. The database
    has no ``(gene_symbol, aa_pos)`` index, so the first lookup loads the gene's
    indexed range and later positions must resolve without another query.
    """

    db_path = tmp_path / "variants.db"
    _write_case_fixture_db(db_path)
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    statements = _trace_statements(monkeypatch)
    clear_gene_metadata_cache()

    first = lookup_variantfeatures_residue("COLL1", position=20)
    assert first is not None
    assert statements, "the first call must reach the database"

    statements.clear()
    assert lookup_variantfeatures_residue("COLL1", position=20) is first
    assert statements == []

    # A different position is a distinct result key but shares the gene range.
    statements.clear()
    lookup_variantfeatures_residue("COLL1", position=900)
    assert statements == []

    # Misses are cached too — an absent gene is the hot path this protects.
    statements.clear()
    assert lookup_variantfeatures_residue("ZZZZ9", position=20) is None
    assert statements
    statements.clear()
    assert lookup_variantfeatures_residue("ZZZZ9", position=20) is None
    assert statements == []

    clear_gene_metadata_cache()
    statements.clear()
    assert lookup_variantfeatures_residue("COLL1", position=20) is not None
    assert statements, "clear_gene_metadata_cache must drop the residue cache"

    clear_gene_metadata_cache()


def test_census_cache_is_keyed_by_path_and_dropped_on_clear(tmp_path, monkeypatch):
    """The census is a process-lifetime cache, so a stale one must be droppable.

    Clearing only the metadata cache leaves the census in place and the gene still
    reads as absent; ``clear_gene_metadata_cache`` has to drop both, or a fixture
    rewritten in place answers from the previous database.
    """

    db_path = tmp_path / "variants.db"

    def write(rows):
        db_path.unlink(missing_ok=True)
        conn = sqlite3.connect(db_path)
        conn.executescript(
            """
            CREATE TABLE variant_consequences (
                gene_symbol TEXT,
                transcript_id TEXT,
                hgvs_p TEXT,
                hgvs_c TEXT,
                aa_pos INTEGER,
                aa_ref TEXT,
                aa_alt TEXT
            );
            CREATE INDEX idx_consequences_gene ON variant_consequences(gene_symbol);
            """
        )
        conn.executemany(
            "INSERT INTO variant_consequences VALUES (?, 'ENST000X', NULL, NULL, ?, 'A', 'G')",
            rows,
        )
        conn.commit()
        conn.close()

    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    write([("AAA1", 100)])
    clear_gene_metadata_cache()

    assert get_gene_metadata("BBB2").protein_length is None

    write([("AAA1", 100), ("Bbb2", 400)])
    gene_metadata.get_gene_metadata.cache_clear()
    assert get_gene_metadata("BBB2").protein_length is None, (
        "the census is expected to be cached, not re-read per lookup"
    )

    clear_gene_metadata_cache()
    assert get_gene_metadata("BBB2").protein_length == 399

    clear_gene_metadata_cache()


def test_census_avoids_the_loose_walk_without_a_usable_index(tmp_path, monkeypatch):
    """Pick the census strategy from the indexes, not unconditionally.

    The recursive ``MIN()`` walk is one index seek per distinct value on an
    indexed column, but one *full scan* per distinct value without an index — so
    on an unindexed column a single ``SELECT DISTINCT`` pass is the cheaper
    census. A partial index does not count: SQLite only uses one when the query
    implies its WHERE clause, which the unconstrained walk never does.
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
        CREATE INDEX idx_variants_partial ON variants(gene) WHERE resnum > 1000;
        INSERT INTO variants VALUES
            ('Leg3', 20, 'Q1', 'W20Y', 'p.Trp20Tyr', 'c.59G>A', 'W', 'Y');
        """
    )
    conn.commit()
    conn.close()
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    statements = _trace_statements(monkeypatch)
    clear_gene_metadata_cache()
    statements.clear()

    # Resolves through the census even though only a partial index exists.
    assert get_gene_metadata("LEG3").protein_length == 20

    census = _census_statements(statements, "variants", "gene")
    assert census, "the census must still run"
    assert not any("RECURSIVE" in sql for sql in census)
    assert not any("UPPER(" in sql for sql in statements)

    clear_gene_metadata_cache()


def test_gene_rows_keeps_the_scanning_predicate_when_no_census_is_possible():
    """A connection with no backing file cannot be censused.

    ``_resolve_variantfeatures_db`` only ever hands back real files, so this guard
    exists for a caller that supplies its own connection. Reading "the census
    returned nothing" as "the gene is absent" is exactly the silent
    metadata-zeroing failure the census path is built to avoid, so the scanning
    predicate has to stay reachable rather than being assumed unnecessary.
    """

    conn = sqlite3.connect(":memory:")
    conn.execute("CREATE TABLE variant_consequences (gene_symbol TEXT, aa_pos INTEGER)")
    conn.execute("INSERT INTO variant_consequences VALUES ('C19orf25', 142)")
    statements: list[str] = []
    conn.set_trace_callback(statements.append)

    rows = gene_metadata._gene_rows(
        conn,
        "SELECT MAX(aa_pos) FROM variant_consequences WHERE {gene_predicate}",
        "variant_consequences",
        "gene_symbol",
        "C19ORF25",
    )

    assert rows == [(142,)]
    assert any("UPPER(gene_symbol)" in sql for sql in statements)
    conn.close()


def test_census_folds_case_like_sqlite_upper_not_like_python(tmp_path, monkeypatch):
    """The census stands in for ``UPPER(column) = ?``, so it must fold identically.

    SQLite's ``upper()`` only folds ASCII; Python's ``str.upper()`` is
    Unicode-aware. ``UPPER('Xé1')`` is ``'Xé1'``, which never equalled the
    ``'XÉ1'`` that ``normalize_gene_symbol`` produces — so this symbol did not
    resolve before and must not start resolving now. Real HGNC symbols are ASCII,
    which is precisely why a Unicode fold here would go unnoticed.
    """

    db_path = tmp_path / "variants.db"
    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        CREATE TABLE variant_consequences (
            gene_symbol TEXT,
            transcript_id TEXT,
            hgvs_p TEXT,
            hgvs_c TEXT,
            aa_pos INTEGER,
            aa_ref TEXT,
            aa_alt TEXT
        );
        INSERT INTO variant_consequences VALUES
            ('Xé1', 'ENST000ACC', NULL, NULL, 200, 'A', 'G');
        """
    )
    conn.commit()
    conn.close()
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    clear_gene_metadata_cache()

    assert get_gene_metadata("xé1").protein_length is None
    assert gene_metadata._sqlite_upper("Xé1") == "Xé1"
    assert "Xé1".upper() == "XÉ1"

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


def test_resolve_gene_symbols_probes_exact_before_scanning(tmp_path):
    """The multi-query helper: one casing decision, reused as an equality test.

    ``_gene_rows`` fixes one query at a time. A caller running several gene-scoped
    queries on one connection should not let each fall back on its own -- the
    ``UPPER()`` form is exactly what costs the full scan.
    """

    db_path = tmp_path / "variants.db"
    _write_case_fixture_db(db_path)
    conn = sqlite3.connect(db_path)
    statements: list[str] = []
    conn.set_trace_callback(statements.append)
    try:
        # Stored uppercase: resolved from the indexed probe, no scan at all.
        assert resolve_variantfeatures_gene_symbols(conn, "COLL1") == ("COLL1",)
        assert not any("UPPER(" in sql for sql in statements)

        # Stored mixed-case only: the real casing comes back, so the caller can
        # still bind an equality test. A bare exact match would return ().
        statements.clear()
        assert resolve_variantfeatures_gene_symbols(conn, "XORF1") == ("Xorf1",)
        assert sum("UPPER(" in sql for sql in statements) == 1

        # Absent: empty tuple, which is the caller's signal to skip its queries.
        assert resolve_variantfeatures_gene_symbols(conn, "NOPE") == ()

        # Lowercase input is normalized like every other gene entry point.
        assert resolve_variantfeatures_gene_symbols(conn, "coll1") == ("COLL1",)
        assert resolve_variantfeatures_gene_symbols(conn, "") == ()

        # Same deliberate narrowing as _gene_rows: when the exact casing matches,
        # the mixed-case twin ('Coll1', aa_pos 900) is not folded in.
        assert "Coll1" not in resolve_variantfeatures_gene_symbols(conn, "COLL1")

        # Other gene-scoped tables are reachable through the same helper.
        assert resolve_variantfeatures_gene_symbols(
            conn, "XORF1", table="transcripts"
        ) == ("Xorf1",)
        assert resolve_variantfeatures_gene_symbols(
            conn, "XORF1", table="genes", column="symbol"
        ) == ("Xorf1",)
    finally:
        conn.close()


def test_resolve_gene_symbols_rejects_unsafe_identifiers(tmp_path):
    """Table/column are interpolated, not bound, so they are validated."""

    db_path = tmp_path / "variants.db"
    _write_case_fixture_db(db_path)
    conn = sqlite3.connect(db_path)
    try:
        with pytest.raises(ValueError):
            resolve_variantfeatures_gene_symbols(
                conn, "COLL1", table="variant_consequences; DROP TABLE genes"
            )
        with pytest.raises(ValueError):
            resolve_variantfeatures_gene_symbols(conn, "COLL1", column="gene_symbol)--")
        assert _has_rows(conn, "genes")
    finally:
        conn.close()


def _has_rows(conn, table):
    return conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0] > 0


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


def test_residue_support_survives_a_missing_composite_index(tmp_path, monkeypatch):
    """A gene-only index must still yield residue evidence, not silent absence.

    The production VariantFeatures slice indexes ``gene_symbol`` alone. Gating
    residue support on a ``(gene_symbol, aa_pos)`` index disabled the lookup in
    exactly that configuration, and an absent residue reads downstream as
    "nothing contradicts this variant" — turning a wrong-residue quarantine
    signal into silent support.
    """

    db_path = tmp_path / "gene_only_index.db"
    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        CREATE TABLE variant_consequences (
            gene_symbol TEXT,
            transcript_id TEXT,
            hgvs_p TEXT,
            hgvs_c TEXT,
            aa_pos INTEGER,
            aa_ref TEXT,
            aa_alt TEXT
        );
        CREATE INDEX idx_vc_gene ON variant_consequences(gene_symbol);
        INSERT INTO variant_consequences VALUES
            ('TESTGENE', 'ENST1', 'ENSP1:p.Gly10Val', 'ENST1:c.29G>T', 10, 'G', 'V');
        """
    )
    conn.close()
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    clear_gene_metadata_cache()

    residue = lookup_variantfeatures_residue("TESTGENE", position=10)
    assert residue is not None, "gene-only index must not disable residue support"
    assert residue.reference_residues == ("G",)

    clear_gene_metadata_cache()


def test_unreadable_variantfeatures_is_not_memoized_as_absent(tmp_path, monkeypatch):
    """A read failure must not be cached as "this gene has no residues"."""

    db_path = tmp_path / "residues.db"
    _write_case_fixture_db(db_path)
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    clear_gene_metadata_cache()

    real_connect = gene_metadata._connect_readonly
    calls = {"n": 0}

    def failing_connect(path):
        calls["n"] += 1
        if calls["n"] == 1:
            raise sqlite3.OperationalError("disk I/O error")
        return real_connect(path)

    monkeypatch.setattr(gene_metadata, "_connect_readonly", failing_connect)
    assert lookup_variantfeatures_residue("COLL1", position=20) is None

    # The failure must not have poisoned the per-gene map: a later healthy
    # read has to reach the database again and produce the real residue.
    monkeypatch.setattr(gene_metadata, "_connect_readonly", real_connect)
    clear_gene_metadata_cache()
    assert lookup_variantfeatures_residue("COLL1", position=20) is not None

    clear_gene_metadata_cache()


def test_oversized_gene_falls_back_to_point_lookups_not_silence(tmp_path, monkeypatch):
    """Too large to prefetch must mean slower, never "no residue evidence"."""

    db_path = tmp_path / "big.db"
    conn = sqlite3.connect(db_path)
    conn.executescript(
        """
        CREATE TABLE variant_consequences (
            gene_symbol TEXT,
            transcript_id TEXT,
            hgvs_p TEXT,
            hgvs_c TEXT,
            aa_pos INTEGER,
            aa_ref TEXT,
            aa_alt TEXT
        );
        CREATE INDEX idx_vc_gene ON variant_consequences(gene_symbol);
        INSERT INTO variant_consequences VALUES
            ('BIGGENE', 'ENST1', 'ENSP1:p.Gly10Val', 'ENST1:c.29G>T', 10, 'G', 'V');
        """
    )
    conn.close()
    monkeypatch.setenv("VARIANTFEATURES_DB", str(db_path))
    # Force the oversized branch without materializing a huge fixture.
    monkeypatch.setattr(gene_metadata, "_RESIDUE_INDEX_MAX_ROWS", 0)
    clear_gene_metadata_cache()

    residue = lookup_variantfeatures_residue("BIGGENE", position=10)
    assert residue is not None, "oversized genes must still resolve residues"
    assert residue.reference_residues == ("G",)
    assert lookup_variantfeatures_residue("BIGGENE", position=9999) is None

    clear_gene_metadata_cache()
