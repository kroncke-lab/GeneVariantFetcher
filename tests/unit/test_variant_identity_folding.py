"""Cross-route identity folding: spelling only, never completion.

The relation must fold two spellings of one fully specified allele and refuse
everything else — conflicts, sparse/rich pairs, cross-gene coincidences, and
any chain that would only connect through a field one side does not supply.

Every fixture is synthetic. Gene names are placeholders so the behaviour cannot
be satisfied by a gene-specific special case.
"""

from __future__ import annotations

import sqlite3

import pytest

from pipeline.variant_identity import (
    FoldVerdict,
    VariantIdentity,
    fold_decision,
    find_identity_conflicts,
)


def ident(**kwargs) -> VariantIdentity:
    return VariantIdentity(**kwargs)


# --- Folds ------------------------------------------------------------------


@pytest.mark.parametrize(
    "left_protein,right_protein",
    [
        ("p.Q403*", "p.Gln403Ter"),
        ("p.Y589C", "p.Tyr589Cys"),
        ("p.Arg147Ter", "p.Arg147*"),
        ("p.Gln42Arg", "p.Q42R"),
    ],
)
def test_pure_spelling_difference_folds(left_protein, right_protein):
    decision = fold_decision(
        ident(
            gene_symbol="GENEA", cdna_notation="c.100A>G", protein_notation=left_protein
        ),
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.100A>G",
            protein_notation=right_protein,
        ),
    )
    assert decision.verdict is FoldVerdict.FOLD, decision.reason


def test_identical_rows_fold():
    row = ident(gene_symbol="GENEA", cdna_notation="c.1A>G", protein_notation="p.M1V")
    assert fold_decision(row, row).folds


# --- Refusals ---------------------------------------------------------------


def test_same_cdna_different_protein_consequence_refuses():
    """A real source contradiction must stay two identities.

    One paper reported the same nucleotide change with two different protein
    consequences. Merging them would launder a contradiction into a single
    confident row.
    """
    decision = fold_decision(
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.1156G>A",
            protein_notation="p.Glu386Gln",
        ),
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.1156G>A",
            protein_notation="p.Glu386Lys",
        ),
    )
    assert decision.verdict is FoldVerdict.REFUSE_CONFLICT


def test_sparse_and_rich_do_not_fold():
    """Completing a missing field asserts evidence the source never gave."""
    decision = fold_decision(
        ident(
            gene_symbol="GENEA", cdna_notation="c.461C>G", protein_notation="p.A154G"
        ),
        ident(gene_symbol="GENEA", protein_notation="p.A154G"),
    )
    assert decision.verdict is FoldVerdict.REFUSE_SPARSE


def test_missing_protein_does_not_bridge_two_conflicting_rows():
    """The NULL-bridge case.

    A row with no protein is 'compatible' with two rows that conflict with each
    other. If the relation were transitive through a missing field, those two
    would end up sharing an identity through the third.
    """
    a = ident(gene_symbol="GENEA", cdna_notation="c.1G>A", protein_notation="p.Arg1Cys")
    b = ident(gene_symbol="GENEA", cdna_notation="c.1G>A", protein_notation="p.Arg1His")
    bridge = ident(gene_symbol="GENEA", cdna_notation="c.1G>A")

    assert fold_decision(a, b).verdict is FoldVerdict.REFUSE_CONFLICT
    assert fold_decision(a, bridge).verdict is FoldVerdict.REFUSE_SPARSE
    assert fold_decision(b, bridge).verdict is FoldVerdict.REFUSE_SPARSE


def test_identical_notation_under_different_genes_never_folds():
    """Cross-gene adversarial case: same string, different biology."""
    decision = fold_decision(
        ident(gene_symbol="GENEA", cdna_notation="c.100A>G", protein_notation="p.M1V"),
        ident(gene_symbol="GENEB", cdna_notation="c.100A>G", protein_notation="p.M1V"),
    )
    assert decision.verdict is FoldVerdict.REFUSE_GENE


def test_missing_gene_never_folds():
    assert not fold_decision(
        ident(cdna_notation="c.100A>G"), ident(cdna_notation="c.100A>G")
    ).folds


def test_variant_class_disagreement_refuses():
    decision = fold_decision(
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.100A>G",
            protein_notation="p.M1V",
            variant_class="missense",
        ),
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.100A>G",
            protein_notation="p.M1V",
            variant_class="nonsense",
        ),
    )
    assert decision.verdict is FoldVerdict.REFUSE_CONFLICT


def test_genomic_position_conflict_refuses():
    decision = fold_decision(
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.100A>G",
            protein_notation="p.M1V",
            genomic_position="chr1:100",
        ),
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.100A>G",
            protein_notation="p.M1V",
            genomic_position="chr1:999",
        ),
    )
    assert decision.verdict is FoldVerdict.REFUSE_CONFLICT


def test_legacy_label_never_aliases_onto_a_coordinate_identity():
    decision = fold_decision(
        ident(gene_symbol="GENEA", legacy_notation="4321delAC"),
        ident(gene_symbol="GENEA", cdna_notation="c.4321delAC"),
    )
    assert decision.verdict is FoldVerdict.REFUSE_SPARSE


def test_relation_is_symmetric():
    a = ident(gene_symbol="GENEA", cdna_notation="c.1A>G", protein_notation="p.M1V")
    b = ident(gene_symbol="GENEA", cdna_notation="c.1A>G")
    assert fold_decision(a, b).verdict is fold_decision(b, a).verdict


# --- Insert-time behaviour --------------------------------------------------


def _db(tmp_path):
    from harvesting.migrate_to_sqlite import create_database_schema

    path = tmp_path / "variants.db"
    conn = create_database_schema(str(path))
    return conn, path


def test_migration_does_not_create_two_rows_for_one_spelling(tmp_path):
    from harvesting.migrate_to_sqlite import get_or_create_variant

    conn, _ = _db(tmp_path)
    cur = conn.cursor()
    first = get_or_create_variant(
        cur,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.1207C>T",
            "protein_notation": "p.Q403*",
        },
    )
    second = get_or_create_variant(
        cur,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.1207C>T",
            "protein_notation": "p.Gln403Ter",
        },
    )
    assert first == second
    assert cur.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_migration_keeps_conflicting_consequences_separate(tmp_path):
    from harvesting.migrate_to_sqlite import get_or_create_variant

    conn, _ = _db(tmp_path)
    cur = conn.cursor()
    first = get_or_create_variant(
        cur,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.1156G>A",
            "protein_notation": "p.Glu386Gln",
        },
    )
    second = get_or_create_variant(
        cur,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.1156G>A",
            "protein_notation": "p.Glu386Lys",
        },
    )
    assert first != second
    assert cur.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 2
    conn.close()


def test_migration_keeps_same_notation_under_two_genes_separate(tmp_path):
    from harvesting.migrate_to_sqlite import get_or_create_variant

    conn, _ = _db(tmp_path)
    cur = conn.cursor()
    a = get_or_create_variant(
        cur,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.100A>G",
            "protein_notation": "p.M1V",
        },
    )
    b = get_or_create_variant(
        cur,
        {
            "gene_symbol": "GENEB",
            "cdna_notation": "c.100A>G",
            "protein_notation": "p.M1V",
        },
    )
    assert a != b
    conn.close()


def test_fold_order_does_not_change_the_outcome(tmp_path):
    """Insertion order must not decide how many identities exist."""
    from harvesting.migrate_to_sqlite import get_or_create_variant

    spellings = [
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.1A>G",
            "protein_notation": "p.Gln403Ter",
        },
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.1A>G",
            "protein_notation": "p.Q403*",
        },
    ]
    totals = []
    for index, order in enumerate((spellings, list(reversed(spellings)))):
        case_dir = tmp_path / f"order{index}"
        case_dir.mkdir()
        conn, _ = _db(case_dir)
        cur = conn.cursor()
        for payload in order:
            get_or_create_variant(cur, dict(payload))
        totals.append(cur.execute("SELECT COUNT(*) FROM variants").fetchone()[0])
        conn.close()
    assert totals[0] == totals[1] == 1


# --- Detector ---------------------------------------------------------------


def test_detector_reports_rather_than_merges(tmp_path):
    """The detector is a worklist. It must never mutate the database."""
    path = tmp_path / "legacy.db"
    conn = sqlite3.connect(path)
    conn.executescript(
        """
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY,
            gene_symbol TEXT, cdna_notation TEXT, protein_notation TEXT,
            legacy_notation TEXT, genomic_position TEXT,
            variant_class TEXT, structural_description TEXT
        );
        CREATE TABLE penetrance_data (
            variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
            affected_count INTEGER, unaffected_count INTEGER
        );
        """
    )
    conn.executemany(
        "INSERT INTO variants (variant_id, gene_symbol, cdna_notation, "
        "protein_notation) VALUES (?,?,?,?)",
        [
            (1, "GENEA", "c.1207C>T", "p.Q403*"),
            (2, "GENEA", "c.1207C>T", "p.Gln403Ter"),
            (3, "GENEA", "c.1156G>A", "p.Glu386Gln"),
            (4, "GENEA", "c.1156G>A", "p.Glu386Lys"),
            (5, "GENEA", "c.461C>G", "p.A154G"),
            (6, "GENEA", None, "p.A154G"),
        ],
    )
    conn.commit()
    before = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    conn.close()

    report = find_identity_conflicts(str(path))

    assert report["foldable_spelling_duplicates"] == 1
    assert report["identity_conflicts"] >= 1
    assert report["sparse_partial_observations"] >= 1

    conn = sqlite3.connect(path)
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == before
    conn.close()


def test_detector_refuses_a_fold_when_counts_disagree(tmp_path):
    """Same spelling, contradictory counts for one paper: adjudicate, not merge."""
    path = tmp_path / "counts.db"
    conn = sqlite3.connect(path)
    conn.executescript(
        """
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY,
            gene_symbol TEXT, cdna_notation TEXT, protein_notation TEXT,
            legacy_notation TEXT, genomic_position TEXT,
            variant_class TEXT, structural_description TEXT
        );
        CREATE TABLE penetrance_data (
            variant_id INTEGER, pmid TEXT, total_carriers_observed INTEGER,
            affected_count INTEGER, unaffected_count INTEGER
        );
        """
    )
    conn.executemany(
        "INSERT INTO variants (variant_id, gene_symbol, cdna_notation, "
        "protein_notation) VALUES (?,?,?,?)",
        [(1, "GENEA", "c.1A>G", "p.Q403*"), (2, "GENEA", "c.1A>G", "p.Gln403Ter")],
    )
    conn.executemany(
        "INSERT INTO penetrance_data VALUES (?,?,?,?,?)",
        [(1, "P1", 3, None, None), (2, "P1", 5, None, None)],
    )
    conn.commit()
    conn.close()

    report = find_identity_conflicts(str(path))

    assert report["foldable_spelling_duplicates"] == 0
    assert any(
        entry["verdict"] == "refuse_count_conflict" for entry in report["conflicts"]
    )


# --- Shared cross-route resolver (database-linkage ingests) -----------------


def _linkage_db(tmp_path):
    path = tmp_path / "linkage.db"
    conn = sqlite3.connect(path)
    conn.executescript(
        """
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT NOT NULL,
            cdna_notation TEXT, protein_notation TEXT, legacy_notation TEXT,
            genomic_position TEXT, variant_class TEXT,
            structural_description TEXT
        );
        """
    )
    return conn


def test_linkage_row_folds_onto_a_paper_row_spelled_differently(tmp_path):
    """The largest measured duplicate source: linkage bypassed identity.

    Each linkage ingest carried its own ``ensure_variant`` matching on the raw
    string, so a database row spelling a call ``p.Gln403Ter`` created a second
    identity beside the paper row's ``p.Q403*``.
    """
    from pipeline.variant_identity import resolve_variant_identity

    conn = _linkage_db(tmp_path)
    paper = resolve_variant_identity(conn, "GENEA", "c.1207C>T", "p.Q403*")
    linkage = resolve_variant_identity(conn, "GENEA", "c.1207C>T", "p.Gln403Ter")

    assert paper == linkage
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_resolver_refuses_to_complete_a_protein_only_linkage_row(tmp_path):
    """A protein-only row must not acquire a neighbouring row's cDNA."""
    from pipeline.variant_identity import resolve_variant_identity

    conn = _linkage_db(tmp_path)
    rich = resolve_variant_identity(conn, "GENEA", "c.461C>G", "p.A154G")
    sparse = resolve_variant_identity(conn, "GENEA", None, "p.A154G")

    assert rich != sparse
    assert (
        conn.execute(
            "SELECT cdna_notation FROM variants WHERE variant_id = ?", (sparse,)
        ).fetchone()[0]
        is None
    )
    conn.close()


def test_resolver_keeps_conflicting_consequences_separate(tmp_path):
    from pipeline.variant_identity import resolve_variant_identity

    conn = _linkage_db(tmp_path)
    a = resolve_variant_identity(conn, "GENEA", "c.1156G>A", "p.Glu386Gln")
    b = resolve_variant_identity(conn, "GENEA", "c.1156G>A", "p.Glu386Lys")
    assert a != b
    conn.close()


def test_resolver_never_folds_across_genes(tmp_path):
    from pipeline.variant_identity import resolve_variant_identity

    conn = _linkage_db(tmp_path)
    a = resolve_variant_identity(conn, "GENEA", "c.100A>G", "p.M1V")
    b = resolve_variant_identity(conn, "GENEB", "c.100A>G", "p.M1V")
    assert a != b
    conn.close()


def test_resolver_is_idempotent(tmp_path):
    from pipeline.variant_identity import resolve_variant_identity

    conn = _linkage_db(tmp_path)
    ids = {
        resolve_variant_identity(conn, "GENEA", "c.100A>G", "p.M1V") for _ in range(4)
    }
    assert len(ids) == 1
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_resolver_refuses_an_ambiguous_match(tmp_path):
    """Two stored rows both claiming a spelling: insert rather than guess."""
    from pipeline.variant_identity import resolve_variant_identity

    conn = _linkage_db(tmp_path)
    conn.executemany(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) "
        "VALUES (?,?,?)",
        [
            ("GENEA", "c.1A>G", "p.Q403*"),
            ("GENEA", "c.1A>G", "p.Gln403Ter"),
        ],
    )
    conn.commit()

    # A third spelling of the same call matches neither stored row exactly but
    # folds with both, so neither may be chosen. An exact match would be
    # unambiguous and is handled by the fast path; this is the case where a
    # winner could only be picked arbitrarily.
    resolved = resolve_variant_identity(conn, "GENEA", "c.1A>G", "p.Gln403*")

    assert resolved not in (1, 2), "resolver picked one of two ambiguous rows"
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 3
    conn.close()


def test_resolver_exact_match_wins_over_ambiguity(tmp_path):
    """An exact spelling match is unambiguous even when a fold rival exists."""
    from pipeline.variant_identity import resolve_variant_identity

    conn = _linkage_db(tmp_path)
    conn.executemany(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) "
        "VALUES (?,?,?)",
        [
            ("GENEA", "c.1A>G", "p.Q403*"),
            ("GENEA", "c.1A>G", "p.Gln403Ter"),
        ],
    )
    conn.commit()

    assert resolve_variant_identity(conn, "GENEA", "c.1A>G", "p.Q403*") == 1
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 2
    conn.close()


# --- Regressions found by adversarial review -------------------------------


def test_two_unparseable_notations_do_not_compare_equal():
    """``None == None`` must not read as "matching cDNA".

    The canonical helper returned ``None`` both for an absent value and for one
    the normalizer could not parse, so two different unparseable strings on the
    same gene satisfied the equality check and folded.
    """
    decision = fold_decision(
        ident(gene_symbol="GENEA", cdna_notation="c.foo", protein_notation="p.M1V"),
        ident(gene_symbol="GENEA", cdna_notation="c.bar", protein_notation="p.M1V"),
    )
    assert decision.verdict is FoldVerdict.REFUSE_CONFLICT


def test_protein_only_spelling_duplicates_fold_across_routes(tmp_path):
    """The duplicate class the raw-string candidate lookup could never see.

    A linkage row carrying only a protein consequence, spelled differently from
    the paper row's, has no byte-identical field to join on. Selecting
    candidates by raw equality therefore offered the predicate nothing and a
    second identity was minted -- the commonest real duplicate.
    """
    from pipeline.variant_identity import resolve_variant_identity

    conn = _linkage_db(tmp_path)
    paper = resolve_variant_identity(conn, "GENEA", None, "p.Q403*")
    linkage = resolve_variant_identity(conn, "GENEA", None, "p.Gln403Ter")

    assert paper == linkage
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_case_only_cdna_difference_folds(tmp_path):
    from pipeline.variant_identity import resolve_variant_identity

    conn = _linkage_db(tmp_path)
    a = resolve_variant_identity(conn, "GENEA", "c.1207C>T", "p.Q403*")
    b = resolve_variant_identity(conn, "GENEA", "c.1207C>T", "p.Gln403Ter")
    assert a == b
    conn.close()


def test_writer_and_detector_implement_the_same_relation(tmp_path):
    """A worklist that reports folds the writer would not make is misleading."""
    from pipeline.variant_identity import (
        find_identity_conflicts,
        resolve_variant_identity,
    )

    conn = _linkage_db(tmp_path)
    conn.execute(
        "CREATE TABLE penetrance_data (variant_id INTEGER, pmid TEXT, "
        "total_carriers_observed INTEGER, affected_count INTEGER, "
        "unaffected_count INTEGER)"
    )
    for protein in ("p.Q403*", "p.Gln403Ter", "p.Gln403*"):
        resolve_variant_identity(conn, "GENEA", None, protein)
    conn.commit()
    path = tmp_path / "linkage.db"
    conn.close()

    report = find_identity_conflicts(str(path))
    assert report["foldable_spelling_duplicates"] == 0, (
        "detector still reports spelling folds the writer already performed"
    )


def test_positions_are_spelling_independent():
    from pipeline.variant_identity import notation_positions

    assert notation_positions("p.Gln403Ter") == ["403"]
    assert notation_positions("p.Q403*") == ["403"]
    assert notation_positions("c.1207C>T", "p.Q403*") == ["1207", "403"]
    assert notation_positions(None, None) == []


def test_a_one_sided_genomic_coordinate_does_not_block_a_fold():
    """A derived coordinate is an annotation, not an identity axis.

    Two rows with byte-identical cDNA *and* protein were being kept apart
    because one carried a genomic coordinate and the other did not. That is
    the opposite of the duplicate this relation exists to prevent: the allele
    is already pinned by gene plus cDNA plus consequence.
    """
    decision = fold_decision(
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.1471C>T",
            protein_notation="p.Arg491Trp",
            genomic_position="chr2:203417496",
        ),
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.1471C>T",
            protein_notation="p.Arg491Trp",
        ),
    )
    assert decision.verdict is FoldVerdict.FOLD, decision.reason


def test_two_different_genomic_coordinates_still_refuse():
    decision = fold_decision(
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.1471C>T",
            protein_notation="p.Arg491Trp",
            genomic_position="chr2:203417496",
        ),
        ident(
            gene_symbol="GENEA",
            cdna_notation="c.1471C>T",
            protein_notation="p.Arg491Trp",
            genomic_position="chr2:999999999",
        ),
    )
    assert decision.verdict is FoldVerdict.REFUSE_CONFLICT


def test_a_missing_cdna_or_protein_still_refuses():
    """The relaxation must not leak into the identity-defining fields."""
    assert (
        fold_decision(
            ident(
                gene_symbol="GENEA", cdna_notation="c.1A>G", protein_notation="p.M1V"
            ),
            ident(gene_symbol="GENEA", protein_notation="p.M1V"),
        ).verdict
        is FoldVerdict.REFUSE_SPARSE
    )
    assert (
        fold_decision(
            ident(
                gene_symbol="GENEA", cdna_notation="c.1A>G", protein_notation="p.M1V"
            ),
            ident(gene_symbol="GENEA", cdna_notation="c.1A>G"),
        ).verdict
        is FoldVerdict.REFUSE_SPARSE
    )


def test_folding_keeps_the_genomic_coordinate(tmp_path):
    """The coordinate must survive the fold, in either insertion order."""
    from harvesting.migrate_to_sqlite import get_or_create_variant

    for order in (0, 1):
        case = tmp_path / f"coord{order}"
        case.mkdir()
        conn, _ = _db(case)
        cur = conn.cursor()
        rows = [
            {
                "gene_symbol": "GENEA",
                "cdna_notation": "c.1471C>T",
                "protein_notation": "p.Arg491Trp",
                "genomic_position": "chr2:203417496",
            },
            {
                "gene_symbol": "GENEA",
                "cdna_notation": "c.1471C>T",
                "protein_notation": "p.Arg491Trp",
            },
        ]
        if order:
            rows.reverse()
        ids = {get_or_create_variant(cur, dict(r)) for r in rows}
        assert len(ids) == 1, "identical cDNA+protein produced two identities"
        stored = cur.execute("SELECT genomic_position FROM variants").fetchone()[0]
        assert stored == "chr2:203417496", "the coordinate was discarded"
        conn.close()
