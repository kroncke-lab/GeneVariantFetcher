"""Unit tests for the paper-metadata backfill (offline; no network)."""

import sqlite3

from scripts.backfill_paper_metadata import (
    _layer,
    backfill_db,
    esummary_to_columns,
    pmids_missing_after_local,
)


def test_esummary_to_columns_maps_and_normalizes_pmc():
    md = {
        "Title": "A paper",
        "FullJournalName": "American journal of obstetrics and gynecology",
        "Source": "Am J Obstet Gynecol",
        "PubDate": "1999 Apr",
        "DOI": "10.1016/s0002-9378(99)70663-0",
        "AuthorList": ["Karlan BY", "Platt LD"],
        "ArticleIds": {"pmc": "123456", "doi": "10.x/ignored"},
    }
    cols = esummary_to_columns(md)
    assert cols["journal"] == "American journal of obstetrics and gynecology"
    assert cols["publication_date"] == "1999 Apr"
    assert cols["doi"] == "10.1016/s0002-9378(99)70663-0"
    assert cols["first_author"] == "Karlan BY"
    assert cols["pmc_id"] == "PMC123456"  # bare accession normalized


def test_esummary_doi_falls_back_to_articleids():
    md = {"ArticleIds": {"doi": "10.1/xyz", "pmc": "PMC9"}, "Source": "J"}
    cols = esummary_to_columns(md)
    assert cols["doi"] == "10.1/xyz"
    assert cols["pmc_id"] == "PMC9"  # already prefixed → unchanged
    assert cols["journal"] == "J"


def test_layer_local_wins_and_none_never_clobbers():
    fetched = {"journal": "Fetched J", "doi": "10.1/a", "pmc_id": "PMC1"}
    abstracts = {"journal": "Local J", "first_author": None}  # None must not clobber
    merged = _layer(fetched, abstracts)
    assert merged["journal"] == "Local J"  # local (later) wins
    assert merged["doi"] == "10.1/a"  # only fetched had it
    assert "first_author" not in merged  # None dropped, not written as empty


def _make_db(path, rows):
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE papers (pmid TEXT PRIMARY KEY, title TEXT)")
    con.executemany("INSERT INTO papers (pmid, title) VALUES (?, ?)", rows)
    con.commit()
    con.close()


def test_backfill_fills_missing_and_is_idempotent(tmp_path):
    db = tmp_path / "G.db"
    _make_db(db, [("111", "P1"), ("222", "P2")])
    abstracts = {
        "111": {"journal": "J1", "first_author": "A1", "publication_date": "2001"}
    }
    fetched = {"111": {"doi": "10.1/a", "pmc_id": "PMC1"}, "222": {"journal": "J2"}}

    filled = backfill_db(db, abstracts, {}, fetched)
    assert filled["journal"] == 2  # J1 (local) + J2 (fetched)
    assert filled["doi"] == 1 and filled["pmc_id"] == 1

    con = sqlite3.connect(db)
    got = dict(con.execute("SELECT pmid, journal FROM papers").fetchall())
    assert got == {"111": "J1", "222": "J2"}
    con.close()

    # Second run is a no-op (nothing empty left to fill).
    again = backfill_db(db, abstracts, {}, fetched)
    assert sum(again.values()) == 0


def test_backfill_never_overwrites_existing(tmp_path):
    db = tmp_path / "G.db"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE papers (pmid TEXT, title TEXT, journal TEXT)")
    con.execute("INSERT INTO papers VALUES ('1', 'P', 'Real Journal')")
    con.commit()
    con.close()

    filled = backfill_db(db, {"1": {"journal": "Wrong"}}, {}, {})
    assert filled["journal"] == 0  # existing value protected
    con = sqlite3.connect(db)
    assert con.execute("SELECT journal FROM papers").fetchone()[0] == "Real Journal"
    con.close()


def test_pmids_missing_after_local():
    # "1" has all ESummary columns locally, but author_affiliation/author_country
    # are never cached locally, so BOTH PMIDs are fetch candidates (affiliations
    # are pulled for every paper as the cohort-origin last resort).
    abstracts = {"1": {"journal": "J", "first_author": "A", "publication_date": "2000"}}
    artifacts = {"1": {"doi": "10.1/a", "pmc_id": "PMC1"}}
    missing = pmids_missing_after_local(["1", "2"], abstracts, artifacts)
    assert missing == ["1", "2"]
