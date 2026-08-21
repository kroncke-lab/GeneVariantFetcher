import json
import sqlite3
from pathlib import Path

from cli.automated_workflow import _apply_paper_scope_gate
from harvesting.migrate_to_sqlite import (
    create_database_schema,
    migrate_extraction_directory,
)
from scripts.extract_figure_variants import ingest_cached_variants
from scripts.recall_recovery import ingest_clinvar
from scripts.recall_recovery import ingest_pubtator
from utils.paper_scope import (
    NONHUMAN_ORTHOLOG_MODEL,
    NONHUMAN_ORTHOLOG_REASON,
    db_paper_scope_exclusions,
    metadata_paper_scope_exclusion_reason,
    purge_excluded_paper_evidence,
)


def _scope_db(path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(path)
    con.executescript(
        """
        CREATE TABLE extraction_metadata (pmid TEXT, model_used TEXT);
        CREATE TABLE variant_papers (variant_id INTEGER, pmid TEXT);
        CREATE TABLE penetrance_data (variant_id INTEGER, pmid TEXT);
        CREATE TABLE papers (pmid TEXT PRIMARY KEY, gene_symbol TEXT);
        """
    )
    con.execute(
        "INSERT INTO extraction_metadata VALUES (?, ?)",
        ("19944633", NONHUMAN_ORTHOLOG_MODEL),
    )
    con.execute("INSERT INTO papers VALUES ('19944633', 'BRCA2')")
    con.execute("INSERT INTO variant_papers VALUES (1, '19944633')")
    con.execute("INSERT INTO penetrance_data VALUES (1, '19944633')")
    con.commit()
    return con


def test_db_scope_exclusion_purges_evidence_but_keeps_audit_record(tmp_path: Path):
    con = _scope_db(tmp_path / "scope.db")

    assert db_paper_scope_exclusions(con) == {"19944633": NONHUMAN_ORTHOLOG_REASON}
    assert purge_excluded_paper_evidence(con) == 2
    assert con.execute("SELECT COUNT(*) FROM variant_papers").fetchone()[0] == 0
    assert con.execute("SELECT COUNT(*) FROM penetrance_data").fetchone()[0] == 0
    assert con.execute("SELECT COUNT(*) FROM extraction_metadata").fetchone()[0] == 1
    assert con.execute("SELECT COUNT(*) FROM papers").fetchone()[0] == 1


def test_late_nonhuman_drop_is_persisted_as_durable_db_exclusion(tmp_path: Path):
    extraction_dir = tmp_path / "extractions"
    extraction_dir.mkdir()
    (extraction_dir / "BRCA2_PMID_19944633.json").write_text(
        json.dumps(
            {
                "paper_metadata": {
                    "pmid": "19944633",
                    "title": "Canine BRCA2",
                    "journal": "Veterinary Journal",
                },
                "variants": [],
                "extraction_metadata": {
                    "model_used": "azure_ai/grok-4.3",
                    "dropped_nonhuman_ortholog": 19,
                },
            }
        ),
        encoding="utf-8",
    )
    con = create_database_schema(str(tmp_path / "late-drop.db"))
    stats = migrate_extraction_directory(con, extraction_dir)

    assert stats["successful"] == 1
    assert db_paper_scope_exclusions(con) == {"19944633": NONHUMAN_ORTHOLOG_REASON}
    assert (
        con.execute(
            "SELECT paper_scope_exclusion_reason FROM extraction_metadata"
        ).fetchone()[0]
        == NONHUMAN_ORTHOLOG_REASON
    )
    assert (
        metadata_paper_scope_exclusion_reason({"dropped_nonhuman_ortholog": 19})
        == NONHUMAN_ORTHOLOG_REASON
    )


def test_explicit_manifest_scope_gate_rejects_canine_title_before_models(
    tmp_path: Path,
):
    abstract = tmp_path / "19944633.json"
    abstract.write_text(
        json.dumps(
            {
                "metadata": {
                    "title": "Single nucleotide variation in exon 11 of canine BRCA2"
                },
                "abstract": "Canine mammary tissue was sequenced.",
            }
        ),
        encoding="utf-8",
    )

    kept, dropped = _apply_paper_scope_gate(
        pmids=["19944633"],
        abstract_records={"19944633": str(abstract)},
        gene_symbol="BRCA2",
        output_path=tmp_path,
    )

    assert kept == []
    assert dropped == [
        ("19944633", "paper studies a non-human ortholog of the target gene")
    ]
    assert (
        "19944633"
        in (tmp_path / "pmid_status" / "paper_scope_exclusions.csv").read_text()
    )


def test_figure_ingest_honors_persisted_paper_scope_exclusion(tmp_path: Path):
    db = tmp_path / "scope.db"
    con = _scope_db(db)
    con.close()

    assert (
        ingest_cached_variants(
            pmid="19944633",
            gene="BRCA2",
            distinct=[{"protein": "p.Lys805Arg"}],
            db_path=db,
        )
        == 0
    )
    with sqlite3.connect(db) as check:
        assert check.execute("SELECT COUNT(*) FROM variant_papers").fetchone()[0] == 0


def test_pubtator_rejects_nonhuman_document_even_without_db_metadata(monkeypatch):
    class Response:
        status_code = 200
        text = json.dumps(
            {
                "PubTator3": [
                    {
                        "id": "19944633",
                        "passages": [
                            {
                                "text": "Single nucleotide variation in canine BRCA2",
                                "annotations": [
                                    {
                                        "text": "K805R",
                                        "infons": {
                                            "type": "ProteinMutation",
                                            "identifier": "CorrespondingGene:675;HGVS:p.Lys805Arg",
                                        },
                                    }
                                ],
                            }
                        ],
                    }
                ]
            }
        )

    monkeypatch.setattr(
        ingest_pubtator.requests, "get", lambda *args, **kwargs: Response()
    )
    monkeypatch.setattr(ingest_pubtator.time, "sleep", lambda *_: None)

    assert (
        ingest_pubtator.fetch_pubtator_mutations(["19944633"], "675", gene="BRCA2")
        == {}
    )


def test_clinvar_preflight_rejects_nonhuman_pubmed_title(monkeypatch):
    class Response:
        status_code = 200

        @staticmethod
        def json():
            return {
                "result": {
                    "19944633": {
                        "title": "Single nucleotide variation in exon 11 of canine BRCA2"
                    }
                }
            }

    monkeypatch.setattr(
        ingest_clinvar.requests, "get", lambda *args, **kwargs: Response()
    )
    monkeypatch.setattr(ingest_clinvar, "_throttle", lambda *_: None)

    assert ingest_clinvar.find_nonhuman_pubmed_pmids(
        {"19944633"}, "BRCA2", "test@example.org", ""
    ) == {"19944633"}
