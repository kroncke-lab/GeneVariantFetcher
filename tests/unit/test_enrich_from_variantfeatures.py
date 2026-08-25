"""variantFeatures enrichment: gene-scoped queries must stay index-friendly.

``upper(gene_symbol) = ?`` cannot use ``idx_consequences_gene``, so it turned
every gene lookup into a full scan. On the real warehouse that cost minutes per
run for byte-identical output, the worst of it on the ``annotations_pathogenicity``
join (43.5M rows). These tests pin the emitted SQL, not just the results, because
the results are identical either way -- the predicate shape *is* the fix.

The mixed-case cases matter as much as the fast path: a bare swap to
``gene_symbol = ?`` silently returns nothing for HGNC symbols that keep a
lowercase ``orf`` (``C19orf25``), which would look like a clean speedup while
dropping the gene's entire enrichment.
"""

import sqlite3
import sys

import pytest

import scripts.enrich_from_variantfeatures as enrich

VF_SCHEMA = """
CREATE TABLE variant_consequences (
    variant_id INTEGER NOT NULL,
    transcript_id TEXT NOT NULL,
    gene_symbol TEXT,
    consequence TEXT,
    hgvs_c TEXT,
    hgvs_p TEXT,
    aa_pos INTEGER,
    aa_ref TEXT,
    aa_alt TEXT,
    is_canonical INTEGER DEFAULT 0,
    is_mane_select INTEGER DEFAULT 0
);
CREATE INDEX idx_consequences_gene ON variant_consequences(gene_symbol, consequence);
CREATE TABLE annotations_pathogenicity (
    variant_id INTEGER NOT NULL,
    predictor TEXT NOT NULL,
    predictor_version TEXT NOT NULL DEFAULT '',
    score REAL,
    PRIMARY KEY (variant_id, predictor, predictor_version)
);
"""

# 'KCNH2' is stored uppercase; 'C19orf25' only in mixed case, as in the real
# warehouse. The second KCNH2 row for variant 1 is the non-MANE twin, so the
# representative-consequence ranking stays observable.
VF_ROWS = """
INSERT INTO variant_consequences
    (variant_id, transcript_id, gene_symbol, consequence, hgvs_c, hgvs_p,
     aa_pos, aa_ref, aa_alt, is_canonical, is_mane_select)
VALUES
    (1, 'ENST_MANE', 'KCNH2', 'missense_variant',
     'ENST_MANE:c.334T>C', 'ENSP_MANE:p.Cys112Arg', 112, 'C', 'R', 1, 1),
    (1, 'ENST_ALT', 'KCNH2', 'missense_variant',
     'ENST_ALT:c.400T>C', 'ENSP_ALT:p.Cys200Arg', 200, 'C', 'R', 0, 0),
    (2, 'ENST_MANE', 'KCNH2', 'stop_gained',
     'ENST_MANE:c.3580A>T', 'ENSP_MANE:p.Lys1194Ter', 1194, 'K', '*', 1, 1),
    (10, 'ENST_ORF', 'C19orf25', 'missense_variant',
     'ENST_ORF:c.89C>G', 'ENSP_ORF:p.Ala30Gly', 30, 'A', 'G', 1, 1);
INSERT INTO annotations_pathogenicity (variant_id, predictor, score) VALUES
    (1, 'alphamissense', 0.91),
    (1, 'revel', 0.77),
    (1, 'cadd_phred', 28.5),
    (1, 'not_a_headline_predictor', 1.0),
    (2, 'revel', 0.55),
    (10, 'alphamissense', 0.42),
    (10, 'eve', 0.31);
"""


@pytest.fixture
def vf_db(tmp_path):
    path = tmp_path / "variants.db"
    conn = sqlite3.connect(path)
    conn.executescript(VF_SCHEMA + VF_ROWS)
    conn.commit()
    conn.close()
    return path


@pytest.fixture
def traced(monkeypatch):
    """Record SQL from the variantFeatures connection only.

    The run-DB query in ``main`` keeps its ``upper()`` fold on purpose, so a
    global trace would conflate the two databases.
    """

    statements: list[str] = []
    real_connect = sqlite3.connect

    def connect(target, *args, **kwargs):
        conn = real_connect(target, *args, **kwargs)
        if "variants.db" in str(target):
            conn.set_trace_callback(statements.append)
        return conn

    monkeypatch.setattr(enrich.sqlite3, "connect", connect)
    return statements


def _hit(statements, table, needle):
    return any(table in sql and needle in sql for sql in statements)


def test_uppercase_gene_uses_the_indexed_predicate_everywhere(vf_db, traced):
    prot_map, cdna_map, meta, max_pos, positions, _residues = enrich.load_vf(
        vf_db, "KCNH2"
    )

    # Results first: the predicate change must not move any output.
    assert prot_map == {"C112R": 1, "K1194*": 2}
    assert cdna_map == {"c.334T>C": 1, "c.3580A>T": 2}
    assert max_pos == 1194
    assert positions == {112, 1194}
    # MANE-select row won over its non-MANE twin for variant 1.
    assert meta[1]["hgvs_p"] == "ENSP_MANE:p.Cys112Arg"
    assert meta[1]["aa_pos"] == 112
    # Only the requested predictors attach; the unlisted one is filtered out.
    assert meta[1]["scores"] == {
        "alphamissense": 0.91,
        "revel": 0.77,
        "cadd_phred": 28.5,
    }
    assert meta[2]["scores"] == {"revel": 0.55}

    # ...and the shape that buys the index.
    assert _hit(traced, "FROM variant_consequences", "gene_symbol IN ('KCNH2')")
    assert _hit(traced, "annotations_pathogenicity", "gene_symbol IN ('KCNH2')")
    assert not any("UPPER(" in sql.upper() for sql in traced)


def test_uppercase_gene_join_plan_uses_the_gene_index(vf_db):
    """The predicate is only worth pinning if SQLite really picks the index."""

    con = sqlite3.connect(f"file:{vf_db}?mode=ro", uri=True)
    try:
        plan = " ".join(
            str(row)
            for row in con.execute(
                "EXPLAIN QUERY PLAN SELECT DISTINCT variant_id "
                "FROM variant_consequences WHERE gene_symbol IN ('KCNH2')"
            )
        )
        assert "idx_consequences_gene" in plan
        upper_plan = " ".join(
            str(row)
            for row in con.execute(
                "EXPLAIN QUERY PLAN SELECT DISTINCT variant_id "
                "FROM variant_consequences WHERE upper(gene_symbol) = 'KCNH2'"
            )
        )
        assert "idx_consequences_gene" not in upper_plan
    finally:
        con.close()


def test_mixed_case_only_gene_still_resolves_and_keeps_the_join_indexed(vf_db, traced):
    """The C19orf25 case: reachable, and the scan is paid once on the small table.

    A bare swap to ``gene_symbol = ?`` returns nothing here. Resolving the casing
    once is what keeps the 43.5M-row pathogenicity join on the index even on this
    fallback path -- a per-query fallback would scan twice.
    """

    prot_map, cdna_map, meta, max_pos, positions, _residues = enrich.load_vf(
        vf_db, "C19ORF25"
    )

    assert prot_map == {"A30G": 10}
    assert cdna_map == {"c.89C>G": 10}
    assert max_pos == 30
    assert positions == {30}
    assert meta[10]["scores"] == {"alphamissense": 0.42, "eve": 0.31}

    # The one scan is the casing probe, and it is on variant_consequences.
    upper_stmts = [sql for sql in traced if "UPPER(" in sql.upper()]
    assert len(upper_stmts) == 1
    assert "DISTINCT gene_symbol" in upper_stmts[0]
    # Both payload queries bind the real stored casing as an equality test.
    assert _hit(traced, "FROM variant_consequences", "gene_symbol IN ('C19orf25')")
    assert _hit(traced, "annotations_pathogenicity", "gene_symbol IN ('C19orf25')")
    assert not any(
        "annotations_pathogenicity" in sql and "UPPER(" in sql.upper() for sql in traced
    )


def test_absent_gene_never_touches_the_pathogenicity_join(vf_db, traced):
    """Every gene without a slice used to pay both full scans for zero rows."""

    prot_map, cdna_map, meta, max_pos, positions, _residues = enrich.load_vf(
        vf_db, "NOSUCHGENE"
    )

    # Same empty shape the callers already handle: nothing matches, nothing is
    # quarantined because max_pos of 0 disables the out-of-range test.
    assert (prot_map, cdna_map, meta, max_pos, positions, _residues) == (
        {},
        {},
        {},
        0,
        set(),
        {},
    )
    assert not any("annotations_pathogenicity" in sql for sql in traced)


def test_enrichment_end_to_end_matches_and_flags(tmp_path, vf_db, capsys):
    """Guards the whole script, so an import or arity slip cannot pass unnoticed.

    ``load_vf`` returns six values against a three-value annotation; only a real
    invocation catches an unpacking regression.
    """

    run_db = tmp_path / "run.db"
    con = sqlite3.connect(run_db)
    con.executescript(
        """
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY,
            gene_symbol TEXT NOT NULL,
            cdna_notation TEXT,
            protein_notation TEXT
        );
        CREATE TABLE variant_papers (variant_id INTEGER, pmid TEXT);
        INSERT INTO variants VALUES
            (1, 'KCNH2', 'c.334T>C', 'p.Cys112Arg'),
            (2, 'KCNH2', 'c.3580A>T', NULL),
            (3, 'kcnh2', NULL, 'p.Ala99999Val'),
            (4, 'KCNH2', NULL, 'p.Arg500Gln');
        INSERT INTO variant_papers VALUES (3, '111'), (3, '222'), (4, '333');
        """
    )
    con.commit()
    con.close()

    fp_out = tmp_path / "fp.csv"
    argv = [
        "enrich_from_variantfeatures.py",
        "--gene",
        "KCNH2",
        "--db",
        str(run_db),
        "--vf",
        str(vf_db),
        "--fp-out",
        str(fp_out),
    ]
    old_argv = sys.argv
    sys.argv = argv
    try:
        assert enrich.main() == 0
    finally:
        sys.argv = old_argv

    con = sqlite3.connect(run_db)
    rows = {
        r[0]: r
        for r in con.execute(
            "SELECT variant_id, matched, match_method, canonical_aa_key, "
            "alphamissense, revel, fp_class FROM vf_enrichment"
        )
    }
    con.close()

    # Protein match carries the in-silico scores through.
    assert rows[1][1:] == (1, "protein", "C112R", 0.91, 0.77, None)
    # cDNA fallback matched the stop variant, which has only one predictor.
    assert rows[2][1:5] == (1, "cdna", "K1194*", None)
    assert rows[2][5] == 0.55
    # Row 3 is the mixed-case run-DB symbol -- it must still be picked up, and its
    # residue is past max_pos, so it is the wrong-gene signal.
    assert rows[3][1] == 0
    assert rows[3][6] == "misparse_out_of_range"
    # In range but not in the warehouse slice: not a false positive.
    assert rows[4][6] == "novel_in_range"

    report = fp_out.read_text().splitlines()
    assert report[0].startswith("gvf_variant_id")
    # Only the high-confidence class is reported, ordered by paper count.
    assert len(report) == 2
    assert report[1].startswith("3,")
    assert report[1].endswith(",misparse_out_of_range,2")


def test_residue_mismatch_splits_wrong_gene_from_numbering_and_unknown():
    """Out-of-range was the only wrong-gene signal, and it is the weaker half.

    BRCA1's P871/E1038/K1183 haplotype sits comfortably inside BRCA2's 3,418
    residues, so under a BRCA2 run every one scored "novel_in_range" and none
    was ever flagged. On the 150-paper re-extraction the residue check flags
    504 BRCA1 / 230 BRCA2 / 64 BMPR2 rows; BMPR2 had ZERO out-of-range, so the
    old classifier reported nothing gene-related for it at all.
    """
    residues = {871: {"L"}, 872: {"P"}, 1038: {"K"}}
    others = {"BRCA1": {871: {"P"}}}

    # disagrees here AND positively matches another gene -> confident wrong gene
    assert (
        enrich.classify_unmatched("p.P871L", "", 3418, {871: {"L"}}, others)
        == "wrong_gene_residue_mismatch"
    )
    # disagrees here, fits this gene one residue over, and NO other gene claims
    # the exact position -> legacy BIC numbering, not a wrong-gene row.
    # Quarantining these deletes real data: 125 BRCA1, 29 BRCA2 and 24 BMPR2
    # rows in the 150-paper re-extraction.
    assert (
        enrich.classify_unmatched(
            "p.P871L", "", 3418, residues, {"BRCA1": {999: {"P"}}}
        )
        == "residue_offset_suspect"
    )

    # PRECEDENCE: an exact match in another gene outranks a fuzzy offset match
    # in this one. On a 3,418-residue protein a +/-3 window finds some amino
    # acid by coincidence almost every time -- BRCA2 has E at 1035/1036, so
    # BRCA1's E1038G was rescued as "numbering" and published under BRCA2
    # across five papers, while BRCA1 has E at exactly 1038.
    assert (
        enrich.classify_unmatched(
            "p.E1038G",
            "",
            3418,
            {1035: {"E"}, 1036: {"E"}, 1038: {"Y"}},
            {"BRCA1": {1038: {"E"}}},
        )
        == "wrong_gene_residue_mismatch"
    )
    # disagrees here and matches no gene we know -> suspect, but not shown to be
    # misattributed. Reported, never quarantined.
    assert (
        enrich.classify_unmatched("p.W500A", "", 3418, {500: {"K"}}, others)
        == "residue_unverified"
    )
    # in range, reference residue agrees -> merely unseen by the warehouse
    assert enrich.classify_unmatched("p.L871P", "", 3418, residues) == "novel_in_range"
    # position the warehouse knows nothing about stays novel, never a false flag
    assert enrich.classify_unmatched("p.M999T", "", 3418, residues) == "novel_in_range"
    # out of range still wins
    assert (
        enrich.classify_unmatched("p.P9999L", "", 3418, residues)
        == "misparse_out_of_range"
    )
    # no residue map (gene absent from the warehouse) must never flag
    assert enrich.classify_unmatched("p.P871L", "", 3418, {}) == "novel_in_range"


def test_classify_unmatched_preserves_strict_source_only_legacy_identity():
    assert (
        enrich.classify_unmatched("", "", 3418, {}, {}, legacy="4321delAC")
        == "legacy_source_notation"
    )


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("P1730H", "P1730H"),
        ("Pro871Leu", "P871L"),
        ("p.P1730H", "P1730H"),
        ("p.Pro871Leu", "P871L"),
        ("pPro871Leu", "P871L"),
    ],
)
def test_proline_is_not_consumed_as_an_hgvs_prefix(source, expected):
    assert enrich.parse_protein_key(source) == expected


@pytest.mark.parametrize(
    ("protein", "position", "reference"),
    [
        ("S1140T", 1140, "S"),
        ("G1319V", 1319, "G"),
        ("D1243N", 1243, "D"),
        ("L1346P", 1346, "L"),
        ("S1382I", 1382, "S"),
        ("E1574K", 1574, "E"),
        ("R1583C", 1583, "R"),
        ("N1722D", 1722, "N"),
    ],
)
def test_scn5a_q1077_long_isoform_uses_lookup_offset_only(protein, position, reference):
    source_identity = protein
    residues = {
        1077: {"E"},
        2015: {"V"},
        position: {"A" if reference != "A" else "C"},
        position - 1: {reference},
    }

    shifted = enrich.scn5a_q1077_isoform_key("SCN5A", protein, 2015, residues)

    # This is a warehouse lookup key. The source/literature identity remains
    # byte-for-byte unchanged in the run DB and scored projection.
    assert shifted == f"{reference}{position - 1}{protein[-1]}"
    assert protein == source_identity


def test_scn5a_q1077_isoform_bridge_is_narrow_and_fail_closed():
    residues = {
        1077: {"E"},
        1078: {"S"},
        1139: {"S"},
        1140: {"T"},
        1318: {"G"},
        1319: {"M"},
        2015: {"W"},
    }

    assert enrich.scn5a_q1077_isoform_key("SCN5A", "G1319V", 2015, residues) == "G1318V"
    assert (
        enrich.scn5a_q1077_isoform_key("SCN5A", "P.G1319V", 2015, residues) == "G1318V"
    )
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "G1319X", 2015, residues) == "G1318*"
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "E1078V", 2015, residues) == "E1077V"
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "W2016X", 2015, residues) == "W2015*"
    # Exact same-coordinate reference wins; do not double-shift a MANE call.
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "M1319V", 2015, residues) is None
    # The insertion residue, out-of-range calls, synonymous and structural
    # alleles are not eligible for the simple missense/stop bridge.
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "E1077V", 2015, residues) is None
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "W2017X", 2015, residues) is None
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "G1319G", 2015, residues) is None
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "G1319del", 2015, residues) is None
    # No DNA-coordinate shifting, no other-gene generalization, and no use
    # after a warehouse reference update to the 2,016-aa isoform.
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "c.3956G>T", 2015, residues) is None
    assert enrich.scn5a_q1077_isoform_key("KCNH2", "G1319V", 2015, residues) is None
    assert enrich.scn5a_q1077_isoform_key("SCN5A", "G1319V", 2016, residues) is None
    long_length_reference = {**residues, 2016: {"A"}}
    assert (
        enrich.scn5a_q1077_isoform_key("SCN5A", "G1319V", 2016, long_length_reference)
        is None
    )
    long_reference = {**residues, 1077: {"Q"}}
    assert (
        enrich.scn5a_q1077_isoform_key("SCN5A", "G1319V", 2015, long_reference) is None
    )


def test_scn5a_q1077_enrichment_preserves_identity_and_other_gene_precedence(
    tmp_path, capsys
):
    vf = tmp_path / "scn5a_variants.db"
    con = sqlite3.connect(vf)
    con.executescript(
        VF_SCHEMA
        + """
        INSERT INTO variant_consequences
          (variant_id, transcript_id, gene_symbol, consequence, hgvs_c, hgvs_p,
           aa_pos, aa_ref, aa_alt, is_canonical, is_mane_select)
        VALUES
          (1, 'SCN5A_MANE', 'SCN5A', 'missense_variant', 'c.3230A>C',
           'p.Glu1077Ala', 1077, 'E', 'A', 1, 1),
          (2, 'SCN5A_MANE', 'SCN5A', 'missense_variant', 'c.3953G>T',
           'p.Gly1318Val', 1318, 'G', 'V', 1, 1),
          (3, 'SCN5A_MANE', 'SCN5A', 'missense_variant', 'c.3956A>G',
           'p.Met1319Val', 1319, 'M', 'V', 1, 1),
          (4, 'SCN5A_MANE', 'SCN5A', 'missense_variant', 'c.6044A>G',
           'p.Ile2015Val', 2015, 'I', 'V', 1, 1),
          (5, 'OTHER_MANE', 'KCNH2', 'missense_variant', 'c.3956G>T',
           'p.Gly1319Val', 1319, 'G', 'V', 1, 1);
        """
    )
    con.commit()
    con.close()

    run_db = tmp_path / "SCN5A.db"
    con = sqlite3.connect(run_db)
    con.executescript(
        """
        CREATE TABLE variants (
          variant_id INTEGER PRIMARY KEY, gene_symbol TEXT NOT NULL,
          cdna_notation TEXT, protein_notation TEXT
        );
        CREATE TABLE variant_papers (variant_id INTEGER, pmid TEXT);
        INSERT INTO variants VALUES (1, 'SCN5A', NULL, 'G1319V');
        INSERT INTO variants VALUES (2, 'SCN5A', NULL, 'G1319W');
        INSERT INTO variants VALUES (3, 'SCN5A', NULL, 'G1319del');
        INSERT INTO variants VALUES (4, 'SCN5A', NULL, 'M1319A');
        INSERT INTO variant_papers VALUES (1, '1'), (2, '1'), (3, '1'), (4, '1');
        """
    )
    con.commit()
    con.close()

    old_argv = sys.argv
    sys.argv = [
        "enrich_from_variantfeatures.py",
        "--gene",
        "SCN5A",
        "--db",
        str(run_db),
        "--vf",
        str(vf),
    ]
    try:
        assert enrich.main() == 0
    finally:
        sys.argv = old_argv

    con = sqlite3.connect(run_db)
    identities = dict(con.execute("SELECT variant_id, protein_notation FROM variants"))
    rows = {
        row[0]: row[1:]
        for row in con.execute(
            "SELECT variant_id, matched, match_method, canonical_aa_key, fp_class "
            "FROM vf_enrichment ORDER BY variant_id"
        )
    }
    con.close()

    assert identities == {1: "G1319V", 2: "G1319W", 3: "G1319del", 4: "M1319A"}
    # Exact shifted variant: canonical VF identity is attached, but the source
    # identity is not rewritten. The other-gene G1319 match cannot preempt the
    # biologically proven same-gene coordinate bridge.
    assert rows[1] == (1, "protein_known_isoform_offset", "G1318V", None)
    # A novel alternate is trusted by the proven reference coordinate without
    # pretending it existed in the warehouse alt catalog.
    assert rows[2] == (
        0,
        "protein_known_isoform_offset",
        None,
        "known_isoform_offset",
    )
    # Structural calls remain in the generic held class; an exact current-
    # coordinate reference remains an ordinary novel in-range allele.
    assert rows[3][3] == "residue_offset_suspect"
    assert rows[4][3] == "novel_in_range"


# --- cDNA coding-range bound ------------------------------------------------
#
# The protein branch of ``classify_unmatched`` has always bounded residues
# against the gene's length, but the cDNA-only branch returned
# ``cdna_only_unmatched`` with no bound -- and that class is on the Variant
# Browser trusted importer's admit-list. A BRCA2 allele written only in c.
# notation therefore walked straight into a BRCA1 collaborator queue: 152 of
# them were live on 2026-08-25, across 13 of the 50 BRCA1 papers, every one a
# combined BRCA1/BRCA2 study.

BRCA1_MAX_AA = 1863  # -> c. ceiling 5592
BRCA2_MAX_AA = 3418  # -> c. ceiling 10257


def test_cdna_coding_ceiling_is_three_bases_plus_the_stop_codon():
    assert enrich.cdna_coding_ceiling(BRCA1_MAX_AA) == 5592
    assert enrich.cdna_coding_ceiling(BRCA2_MAX_AA) == 10257
    # An unknown gene has no ceiling and must not be bounded into quarantine.
    assert enrich.cdna_coding_ceiling(0) == 0


@pytest.mark.parametrize(
    "cdna",
    [
        "c.5946delT",  # BRCA2 Ashkenazi founder allele
        "c.6275_6276delTT",  # BRCA2 founder allele
        "c.9485-1G>A",  # BRCA2 splice acceptor
        "c.5722_5723delCT",
    ],
)
def test_brca2_cdna_in_a_brca1_run_is_quarantined_not_admitted(cdna):
    assert enrich.cdna_out_of_range(cdna, BRCA1_MAX_AA) is True
    assert (
        enrich.classify_unmatched("", cdna, BRCA1_MAX_AA, {}, {}) == "cdna_out_of_range"
    )
    # ...and the same allele is unremarkable in its own gene.
    assert enrich.cdna_out_of_range(cdna, BRCA2_MAX_AA) is False
    assert (
        enrich.classify_unmatched("", cdna, BRCA2_MAX_AA, {}, {})
        == "cdna_only_unmatched"
    )


@pytest.mark.parametrize(
    "cdna",
    [
        "c.5592A>G",  # the last coding base
        "c.5592+34G>T",  # intronic offset from an in-range base
        "c.181T>G",
        "c.-19G>A",  # 5' UTR: legitimately outside the CDS
        "c.*103A>C",  # 3' UTR: likewise
    ],
)
def test_in_range_and_utr_cdna_stays_admissible_for_brca1(cdna):
    assert enrich.cdna_out_of_range(cdna, BRCA1_MAX_AA) is False
    assert (
        enrich.classify_unmatched("", cdna, BRCA1_MAX_AA, {}, {})
        == "cdna_only_unmatched"
    )


def test_cdna_bound_never_fires_without_a_known_protein_length():
    """A gene variantFeatures does not cover must keep its prior behaviour."""
    assert enrich.cdna_out_of_range("c.99999A>G", 0) is False
    assert (
        enrich.classify_unmatched("", "c.99999A>G", 0, {}, {}) == "cdna_only_unmatched"
    )


def test_source_only_legacy_identity_still_outranks_the_cdna_bound():
    """Legacy notation is resolved before the range check, as it always was."""
    assert (
        enrich.classify_unmatched(
            "", "c.9485-1G>A", BRCA1_MAX_AA, {}, {}, legacy="4321delAC"
        )
        == "legacy_source_notation"
    )


def test_a_row_with_no_notation_is_unaffected_by_the_cdna_bound():
    assert (
        enrich.classify_unmatched("", "", BRCA1_MAX_AA, {}, {}) == "no_notation_suspect"
    )
