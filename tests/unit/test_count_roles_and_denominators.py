"""Typed count roles, zero provenance, and the finalized-list denominator.

Two independent guarantees are pinned here.

*Counts.* A number's role (individual, family, proband, cohort total) becomes a
persisted column so a downstream consumer can no longer read 41 families as 41
patients. An unsupported zero is never invented, and a historical zero whose
provenance was never recorded is never rewritten in either direction.

*Denominators.* Run artifacts are rendered from the papers extraction actually
produced records for, not from what an acquisition step intended to fetch.

Fixtures are synthetic. Identifiers are placeholders.
"""

from __future__ import annotations

import sqlite3

import pytest

from pipeline.source_ledger import (
    CLASS_ABSTRACT_ONLY,
    CLASS_FULLTEXT,
    CLASS_UNVERIFIED,
    build_source_ledger,
)

# --- Count roles ------------------------------------------------------------


def _schema(tmp_path):
    from harvesting.migrate_to_sqlite import create_database_schema

    return create_database_schema(str(tmp_path / "counts.db"))


def _insert(conn, variant_data, pmid="P1"):
    from harvesting.migrate_to_sqlite import insert_paper_metadata, insert_variant_data

    cur = conn.cursor()
    insert_paper_metadata(
        cur,
        {
            "paper_metadata": {"pmid": pmid, "title": "Synthetic fixture"},
            "variants": [variant_data],
        },
    )
    insert_variant_data(cur, pmid, variant_data)
    conn.commit()
    return cur


def _penetrance(conn, pmid="P1"):
    conn.row_factory = sqlite3.Row
    return conn.execute(
        "SELECT * FROM penetrance_data WHERE pmid = ?", (pmid,)
    ).fetchone()


def test_declared_family_count_is_persisted_as_its_own_role(tmp_path):
    """A family count must not be readable as an individual count."""
    conn = _schema(tmp_path)
    _insert(
        conn,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.100A>G",
            "penetrance_data": {"total_carriers_observed": 41},
            "count_provenance": {
                "carriers_column_label": "Families",
                "carriers_count_type": "family_count",
            },
        },
    )
    row = _penetrance(conn)
    assert row["carriers_role"] == "family_count"
    assert row["total_carriers_observed"] == 41
    conn.close()


def test_undeclared_role_is_unknown_not_per_variant_carrier(tmp_path):
    """Absence of a declaration must never default to an individual count."""
    conn = _schema(tmp_path)
    _insert(
        conn,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.100A>G",
            "penetrance_data": {"total_carriers_observed": 7},
        },
    )
    assert _penetrance(conn)["carriers_role"] == "unknown"
    conn.close()


def test_sourced_zero_is_marked_sourced_and_kept(tmp_path):
    conn = _schema(tmp_path)
    _insert(
        conn,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.100A>G",
            "penetrance_data": {
                "total_carriers_observed": 5,
                "affected_count": 5,
                "unaffected_count": 0,
            },
            "count_provenance": {
                "unaffected_column_label": "Unaffected carriers",
                "unaffected_count_type": "per_variant_carrier",
            },
        },
    )
    row = _penetrance(conn)
    assert row["unaffected_count"] == 0
    assert row["unaffected_zero_provenance"] == "sourced"
    conn.close()


def test_zero_without_a_named_column_is_unknown_not_nulled(tmp_path):
    """Unknown-provenance zero is a third state, not NULL and not a claim.

    Nulling it would delete a possibly real complete-penetrance observation;
    trusting it would report a possibly defaulted value as a finding. It is
    preserved exactly and labelled.
    """
    conn = _schema(tmp_path)
    _insert(
        conn,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.100A>G",
            "penetrance_data": {
                "total_carriers_observed": 5,
                "unaffected_count": 0,
            },
        },
    )
    row = _penetrance(conn)
    assert row["unaffected_count"] == 0, "a stored zero must not be rewritten"
    assert row["unaffected_zero_provenance"] == "unknown"
    conn.close()


def test_unobserved_partition_stays_null(tmp_path):
    conn = _schema(tmp_path)
    _insert(
        conn,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.100A>G",
            "penetrance_data": {"total_carriers_observed": 5},
        },
    )
    row = _penetrance(conn)
    assert row["affected_count"] is None
    assert row["unaffected_count"] is None
    assert row["affected_zero_provenance"] is None
    conn.close()


def test_impossible_partition_is_recorded_not_repaired(tmp_path):
    """Reconciliation detects. It must not make the arithmetic close."""
    conn = _schema(tmp_path)
    _insert(
        conn,
        {
            "gene_symbol": "GENEA",
            "cdna_notation": "c.100A>G",
            "penetrance_data": {
                "total_carriers_observed": 4,
                "affected_count": 3,
                "unaffected_count": 3,
            },
            "count_provenance": {
                "affected_column_label": "Affected",
                "affected_count_type": "per_variant_carrier",
                "unaffected_column_label": "Unaffected",
                "unaffected_count_type": "per_variant_carrier",
            },
        },
    )
    row = _penetrance(conn)
    if row is not None:
        # If the row is stored, the integers must be untouched and the
        # contradiction recorded rather than silently balanced.
        assert (row["affected_count"], row["unaffected_count"]) == (3, 3)
        assert row["total_carriers_observed"] == 4
        assert "partitions_exceed_total" in (row["count_reconciliation"] or "")
    conn.close()


def test_upgrade_adds_role_columns_without_touching_stored_zeros(tmp_path):
    """A pre-existing database gains the axis; its data is not reinterpreted."""
    from harvesting.migrate_to_sqlite import upgrade_database_schema

    path = tmp_path / "legacy.db"
    conn = sqlite3.connect(path)
    conn.executescript(
        """
        CREATE TABLE penetrance_data (
            penetrance_id INTEGER PRIMARY KEY,
            variant_id INTEGER, pmid TEXT,
            total_carriers_observed INTEGER, affected_count INTEGER,
            unaffected_count INTEGER, uncertain_count INTEGER,
            penetrance_percentage REAL
        );
        CREATE TABLE extraction_metadata (
            extraction_id INTEGER PRIMARY KEY, pmid TEXT
        );
        """
    )
    conn.execute(
        "INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed,"
        " affected_count, unaffected_count) VALUES (1, 'P1', 9, 9, 0)"
    )
    conn.commit()

    upgrade_database_schema(conn)

    conn.row_factory = sqlite3.Row
    row = conn.execute("SELECT * FROM penetrance_data").fetchone()
    assert row["unaffected_count"] == 0, "historical zero must survive the upgrade"
    assert row["carriers_role"] == "unknown"
    conn.close()


# --- Finalized-list denominators -------------------------------------------


def _record(pmid, *, declared, source_file=None, variants=None, size=100_000):
    meta = {"source_type": declared, "source_size_bytes": size}
    if declared == CLASS_ABSTRACT_ONLY:
        meta["abstract_only"] = 1
    if source_file:
        meta["source_file"] = str(source_file)
    return {
        "paper_metadata": {"pmid": pmid},
        "extraction_metadata": meta,
        "variants": variants or [],
    }


def _fulltext_file(tmp_path, pmid):
    path = tmp_path / f"{pmid}_FULL_CONTEXT.md"
    path.write_text("# Article\n\n" + ("body text " * 400), encoding="utf-8")
    return path


def test_denominator_follows_extraction_not_the_download_step(tmp_path):
    """The measured defect: a resume run reported 0/50 full text.

    ``download_fulltext`` labels every paper it did not itself fetch as
    abstract-only, so a run served from the corpus cache reported zero full
    text while its own per-paper records said every paper was full text.
    """
    records = [
        _record(
            f"P{i}",
            declared=CLASS_FULLTEXT,
            source_file=_fulltext_file(tmp_path, f"P{i}"),
        )
        for i in range(5)
    ]
    ledger = build_source_ledger(
        records,
        requested_pmids=[f"P{i}" for i in range(5)],
        search_dirs=[tmp_path],
    )
    summary = ledger.as_completeness_dict()
    assert summary["papers_with_fulltext"] == 5
    assert summary["papers_abstract_only"] == 0
    assert summary["denominator_source"] == "finalized_extraction_records"


def test_every_reported_count_equals_the_length_of_its_list(tmp_path):
    """Counts were reported for the population while lists were truncated."""
    records = [
        _record(
            f"P{i}",
            declared=CLASS_FULLTEXT,
            source_file=_fulltext_file(tmp_path, f"P{i}"),
            variants=[{"penetrance_data": {"total_carriers_observed": 1}}],
        )
        for i in range(60)
    ]
    summary = build_source_ledger(
        records, search_dirs=[tmp_path]
    ).as_completeness_dict()
    assert len(summary["fulltext_pmids"]) == summary["papers_with_fulltext"]
    assert len(summary["abstract_only_pmids"]) == summary["papers_abstract_only"]
    assert len(summary["single_carrier_pmids"]) == summary["single_carrier_papers"]
    assert len(summary["zero_variant_pmids"]) == summary["zero_variant_papers"]


def test_single_carrier_predicate_reads_a_real_count_key(tmp_path):
    """The predicate read two keys that do not exist in the schema.

    ``variant["carriers_total"]`` and ``variant["total_carriers"]`` are absent,
    so ``... or 1`` made the expression constant and every paper in every run
    was flagged as single-carrier.
    """
    records = [
        _record(
            "P_MANY",
            declared=CLASS_FULLTEXT,
            source_file=_fulltext_file(tmp_path, "P_MANY"),
            variants=[{"penetrance_data": {"total_carriers_observed": 7}}],
        ),
        _record(
            "P_ONE",
            declared=CLASS_FULLTEXT,
            source_file=_fulltext_file(tmp_path, "P_ONE"),
            variants=[{"penetrance_data": {"total_carriers_observed": 1}}],
        ),
        _record(
            "P_NOCOUNT",
            declared=CLASS_FULLTEXT,
            source_file=_fulltext_file(tmp_path, "P_NOCOUNT"),
            variants=[{"penetrance_data": {}}],
        ),
    ]
    ledger = build_source_ledger(records, search_dirs=[tmp_path])
    assert ledger.single_carrier_pmids == ["P_ONE"]


def test_declared_fulltext_with_unreadable_source_is_unverified(tmp_path):
    """The inverse bug: a declaration is not proof of article bytes."""
    record = _record("P1", declared=CLASS_FULLTEXT, source_file=tmp_path / "gone.md")
    ledger = build_source_ledger([record], search_dirs=[tmp_path])
    assert ledger.rows[0].effective_class == CLASS_UNVERIFIED
    assert ledger.rows[0].discrepancy


def test_declared_fulltext_over_stub_bytes_is_reported_as_a_discrepancy(tmp_path):
    stub = tmp_path / "P1_FULL_CONTEXT.md"
    stub.write_text(
        "# ABSTRACT-ONLY FALLBACK\n\nfull text could not be retrieved; "
        "contains only the pubmed abstract\n" + "x" * 2000,
        encoding="utf-8",
    )
    record = _record("P1", declared=CLASS_FULLTEXT, source_file=stub)
    ledger = build_source_ledger([record], search_dirs=[tmp_path])
    assert ledger.rows[0].effective_class == CLASS_ABSTRACT_ONLY
    assert ledger.rows[0].discrepancy


def test_requested_paper_with_no_extraction_record_is_reported_missing(tmp_path):
    """Real coverage loss must not be hidden by a finalized-list denominator."""
    records = [
        _record(
            "P1", declared=CLASS_FULLTEXT, source_file=_fulltext_file(tmp_path, "P1")
        )
    ]
    ledger = build_source_ledger(
        records, requested_pmids=["P1", "P2"], search_dirs=[tmp_path]
    )
    assert ledger.missing_pmids == ["P2"]


@pytest.mark.parametrize("gene", ["GENEA", "GENEB"])
def test_ledger_is_gene_agnostic(tmp_path, gene):
    record = _record(
        f"{gene}_P1",
        declared=CLASS_FULLTEXT,
        source_file=_fulltext_file(tmp_path, f"{gene}_P1"),
    )
    ledger = build_source_ledger([record], search_dirs=[tmp_path])
    assert ledger.rows[0].effective_class == CLASS_FULLTEXT


def test_trust_gate_reads_the_persisted_role_when_the_json_is_gone(tmp_path):
    """The role columns must be a protection, not only a record.

    The trust gate reads count roles from a JSON blob on ``variant_papers``,
    reached through a LEFT JOIN. A missing or malformed blob silently dropped
    the role, so a family count read as an individual count. The persisted
    column sits on the fact row itself and recovers that.
    """
    from harvesting.migrate_to_sqlite import create_database_schema
    from pipeline.trust_gate import apply_trust_gate

    db = tmp_path / "roles.db"
    conn = create_database_schema(str(db))
    conn.execute("INSERT INTO papers (pmid, gene_symbol) VALUES ('P1', 'GENEA')")
    conn.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation) VALUES ('GENEA', 'c.1A>G')"
    )
    # variant_papers row deliberately absent: the JSON provenance is gone.
    conn.execute(
        "INSERT INTO penetrance_data (variant_id, pmid, total_carriers_observed, "
        "carriers_role) VALUES (1, 'P1', 41, 'family_count')"
    )
    conn.commit()
    conn.close()

    apply_trust_gate(str(db))

    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    row = conn.execute("SELECT * FROM penetrance_data").fetchone()
    assert "family_count_not_carrier" in (row["trust_reasons"] or ""), (
        "a declared family count was not recognised without the JSON blob"
    )
    assert row["total_carriers_observed"] == 41, "the raw count must be preserved"
    conn.close()
