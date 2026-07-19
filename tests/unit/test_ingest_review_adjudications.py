"""Tests for scripts/ingest_review_adjudications.py.

Exercises the Variant_Browser → GVF adjudication round-trip: a lead-approved CSV
(the schema ``manage.py export_gold_standard`` emits) is matched against a real run DB
built through the actual ``migrate_to_sqlite`` path, so the (pmid, canonical
notation) match contract is checked against the same schema the recall scorer
reads. Covers all six verdicts, the correction overlay (extracted + corrected
values kept side by side), the follow-up queue, the summary deltas, idempotency,
and the no-DB path.
"""

import csv
import hashlib
import io
import json
import sqlite3
import sys

import pytest

import harvesting.migrate_to_sqlite as migrate
import scripts.ingest_review_adjudications as ingest

PMID = "10086971"

# protein_notation -> (carriers, affected, unaffected) seeded into the run DB.
DB_VARIANTS = {
    "p.Ser818Leu": (2, 2, 0),
    "p.Ala561Val": (1, 1, 0),
    "p.Gly628Ser": (4, 3, 1),
    "p.Asn588Lys": (2, 1, 1),
}

# Export rows: one per verdict. (variant_label, source_notation, verdict,
# corrected_affected, corrected_unaffected, corrected_total, classification, comment)
EXPORT_ROWS = [
    ("S818L", "p.Ser818Leu", "confirm", "", "", "", "", "looks right"),
    ("A561V", "p.Ala561Val", "correct_counts", "5", "2", "7", "Pathogenic", "miscount"),
    ("G628S", "p.Gly628Ser", "wrong_variant", "", "", "", "", "actually G628A"),
    ("N588K", "p.Asn588Lys", "wrong_paper", "", "", "", "", "wrong PMID linked"),
    ("R100Q", "p.Arg100Gln", "missing", "3", "1", "4", "", "GVF missed this one"),
    ("T200M", "p.Thr200Met", "other", "", "", "", "", "needs a domain expert"),
]


def _build_db(tmp_path):
    """Build a real KCNH2 run DB through the actual migration path."""
    ext_dir = tmp_path / "extractions"
    ext_dir.mkdir()
    doc = {
        "paper_metadata": {"pmid": PMID, "title": "Test paper", "gene_symbol": "KCNH2"},
        "variants": [
            {
                "gene_symbol": "KCNH2",
                "protein_notation": protein,
                "penetrance_data": {
                    "total_carriers_observed": carriers,
                    "affected_count": affected,
                    "unaffected_count": unaffected,
                },
            }
            for protein, (carriers, affected, unaffected) in DB_VARIANTS.items()
        ],
    }
    (ext_dir / f"{PMID}.json").write_text(json.dumps(doc))
    db_path = tmp_path / "KCNH2.db"
    conn = migrate.create_database_schema(str(db_path))
    try:
        cur = conn.cursor()
        ok, msg = migrate.migrate_extraction_file(cur, ext_dir / f"{PMID}.json")
        assert ok, msg
        conn.commit()
    finally:
        conn.close()
    return db_path


def _write_export(tmp_path):
    path = tmp_path / "adjudications.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "gene",
                "record_key",
                "status",
                "revision",
                "source_reviewer_user_id",
                "source_reviewer",
                "decided_by_user_id",
                "decided_by",
                "variant_label",
                "source_notation",
                "pmid",
                "verdict",
                "corrected_affected",
                "corrected_unaffected",
                "corrected_total",
                "corrected_classification",
                "comment",
                "adjudicator",
                "updated_at",
            ]
        )
        for i, (label, src, verdict, c_aff, c_unaff, c_tot, cls, comment) in enumerate(
            EXPORT_ROWS
        ):
            writer.writerow(
                [
                    "KCNH2",
                    f"gold-record-{i}",
                    "gold_standard",
                    "1",
                    f"account-{i}",
                    f"reviewer{i}",
                    "lead-account",
                    "lead",
                    label,
                    src,
                    PMID,
                    verdict,
                    c_aff,
                    c_unaff,
                    c_tot,
                    cls,
                    comment,
                    f"reviewer{i}",
                    "2026-06-13T10:00:00",
                ]
            )
    return path


def _run(tmp_path, *extra_argv):
    out_dir = tmp_path / "adj_out"
    argv = [
        "ingest_review_adjudications.py",
        "--export-csv",
        str(_write_export(tmp_path)),
        "--out-dir",
        str(out_dir),
        *extra_argv,
    ]
    old = sys.argv
    sys.argv = argv
    try:
        assert ingest.main() == 0
    finally:
        sys.argv = old
    return out_dir


def _read_overlay(out_dir):
    with (out_dir / "KCNH2_review_adjudications.csv").open(encoding="utf-8") as handle:
        return {r["source_notation"]: r for r in csv.DictReader(handle)}


def test_canonicalization_matches_scorer():
    # The overlay key must be identical for 3-letter and 1-letter spellings,
    # which is how it round-trips to the DB's aggregated (pmid, variant) key.
    assert ingest._variant_key("p.Ser818Leu") == ingest._variant_key("S818L")


def test_verdicts_map_and_match_against_real_db(tmp_path):
    db = _build_db(tmp_path)
    out_dir = _run(tmp_path, "--db", f"KCNH2={db}")
    rows = _read_overlay(out_dir)
    assert len(rows) == len(EXPORT_ROWS)

    confirm = rows["p.Ser818Leu"]
    assert confirm["record_key"] == "gold-record-0"
    assert confirm["status"] == "gold_standard"
    assert confirm["revision"] == "1"
    assert confirm["source_reviewer_user_id"] == "account-0"
    assert confirm["source_reviewer"] == "reviewer0"
    assert confirm["decided_by_user_id"] == "lead-account"
    assert confirm["decided_by"] == "lead"
    assert confirm["action"] == "gold_confirmed"
    assert confirm["match_status"] == "matched"
    # Both extracted and corrected survive; confirm keeps the extracted value.
    assert confirm["extracted_affected"] == "2"

    override = rows["p.Ala561Val"]
    assert override["action"] == "count_override"
    assert override["match_status"] == "matched"
    assert override["extracted_affected"] == "1"  # original kept
    assert override["corrected_affected"] == "5"  # adjudication kept

    assert rows["p.Gly628Ser"]["action"] == "false_positive"
    assert rows["p.Gly628Ser"]["match_status"] == "matched"
    assert rows["p.Asn588Lys"]["action"] == "excluded"
    assert rows["p.Asn588Lys"]["match_status"] == "matched"

    # missing/other reference variants GVF never extracted -> unmatched is expected.
    assert rows["p.Arg100Gln"]["action"] == "followup_missing"
    assert rows["p.Arg100Gln"]["match_status"] == "unmatched"
    assert rows["p.Thr200Met"]["action"] == "followup_other"


def test_followup_queue_and_summary(tmp_path):
    db = _build_db(tmp_path)
    out_dir = _run(tmp_path, "--db", f"KCNH2={db}")

    with (out_dir / "review_followup_queue.csv").open(encoding="utf-8") as handle:
        queue = list(csv.DictReader(handle))
    # Only missing + other queue (the four actionable verdicts all matched).
    queued_verdicts = sorted(r["verdict"] for r in queue)
    assert queued_verdicts == ["missing", "other"]

    summary = json.loads((out_dir / "review_adjudications_summary.json").read_text())
    kcnh2 = summary["KCNH2"]
    assert kcnh2["total"] == 6
    assert kcnh2["actions"]["count_override"] == 1
    assert kcnh2["count_override_matched"] == 1
    # A561V: affected 1->5 (+4), unaffected 0->2 (+2). Only matched overrides count.
    assert kcnh2["net_affected_delta"] == 4
    assert kcnh2["net_unaffected_delta"] == 2


def test_idempotent(tmp_path):
    db = _build_db(tmp_path)
    out_dir = _run(tmp_path, "--db", f"KCNH2={db}")
    first = (out_dir / "KCNH2_review_adjudications.csv").read_bytes()
    out_dir2 = _run(tmp_path, "--db", f"KCNH2={db}")
    assert (out_dir2 / "KCNH2_review_adjudications.csv").read_bytes() == first


def test_no_db_mode_queues_everything(tmp_path):
    out_dir = _run(tmp_path, "--no-db")
    rows = _read_overlay(out_dir)
    assert all(r["match_status"] == "no_db" for r in rows.values())
    assert all(r["extracted_affected"] == "" for r in rows.values())
    with (out_dir / "review_followup_queue.csv").open(encoding="utf-8") as handle:
        queue = list(csv.DictReader(handle))
    # With no DB, nothing can be verified, so every row needs human follow-up.
    assert len(queue) == len(EXPORT_ROWS)


def test_missing_export_errors(tmp_path):
    sys_argv = sys.argv
    sys.argv = [
        "ingest_review_adjudications.py",
        "--export-csv",
        str(tmp_path / "does_not_exist.csv"),
    ]
    try:
        with pytest.raises(SystemExit):
            ingest.main()
    finally:
        sys.argv = sys_argv


def test_multi_reviewer_audit_export_is_rejected(tmp_path):
    path = tmp_path / "unsafe_adjudications.csv"
    path.write_text("gene,pmid,verdict,adjudicator\nKCNH2,1,confirm,nate\n")

    with pytest.raises(SystemExit, match="Do not ingest the multi-reviewer"):
        ingest._read_export(path)


def test_disputed_gold_audit_export_is_rejected(tmp_path):
    path = _write_export(tmp_path)
    rows = list(csv.DictReader(path.open()))
    rows[0]["status"] = "disputed"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    with pytest.raises(SystemExit, match="non-accepted status"):
        ingest._read_export(path)


class _FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class _FakeOpener:
    def __init__(self, payload):
        self.payload = payload
        self.request = None

    def open(self, request, timeout):
        self.request = request
        return _FakeResponse(json.dumps(self.payload).encode())


def _live_payload(rows):
    canonical = json.dumps(
        rows,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode()
    return {
        "ok": True,
        "schema_version": 1,
        "source": "variant_browser_azure_review",
        "generated_at": "2026-07-19T20:00:00+00:00",
        "dataset_revision": hashlib.sha256(canonical).hexdigest(),
        "columns": list(rows[0]) if rows else sorted(ingest.GOLD_EXPORT_MARKERS),
        "record_count": len(rows),
        "records": rows,
    }


def test_live_sync_populates_atomic_sqlite_cache_with_identity_audit(tmp_path):
    source_rows = list(csv.DictReader(_write_export(tmp_path).open()))
    payload = _live_payload(source_rows)
    opener = _FakeOpener(payload)
    cache_db = tmp_path / "review_gold.sqlite3"
    out_dir = tmp_path / "live_out"

    state = ingest.sync_live_gold(
        source_url="https://variantbrowser.org/review/api/gold-standard/",
        token="secret-token-" * 4,
        cache_db=cache_db,
        out_dir=out_dir,
        no_db=True,
        opener=opener,
        announce=False,
    )

    assert state["record_count"] == len(EXPORT_ROWS)
    assert state["dataset_revision"] == payload["dataset_revision"]
    assert opener.request.get_header("Authorization").startswith("Bearer ")
    conn = sqlite3.connect(cache_db)
    try:
        gold = conn.execute(
            "SELECT record_key, source_reviewer_user_id, decided_by_user_id "
            "FROM review_gold_records ORDER BY record_key"
        ).fetchall()
        overlays = conn.execute("SELECT COUNT(*) FROM review_gold_overlays").fetchone()[
            0
        ]
        sync_state = conn.execute(
            "SELECT dataset_revision, record_count FROM review_gold_sync_state"
        ).fetchone()
    finally:
        conn.close()
    assert gold[0] == ("gold-record-0", "account-0", "lead-account")
    assert overlays == len(EXPORT_ROWS)
    assert sync_state == (payload["dataset_revision"], len(EXPORT_ROWS))

    # A second complete snapshot replaces rather than duplicates the cache.
    ingest.sync_live_gold(
        source_url="https://variantbrowser.org/review/api/gold-standard/",
        token="secret-token-" * 4,
        cache_db=cache_db,
        out_dir=out_dir,
        no_db=True,
        opener=_FakeOpener(payload),
        announce=False,
    )
    conn = sqlite3.connect(cache_db)
    try:
        assert conn.execute("SELECT COUNT(*) FROM review_gold_records").fetchone()[
            0
        ] == len(EXPORT_ROWS)
    finally:
        conn.close()


def test_live_sync_rejects_disputed_or_checksum_tampered_payload(tmp_path):
    source_rows = list(csv.DictReader(_write_export(tmp_path).open()))
    disputed = [dict(row) for row in source_rows]
    disputed[0]["status"] = "disputed"
    with pytest.raises(ingest.GoldSyncError, match="non-accepted status"):
        ingest.fetch_live_gold(
            "https://variantbrowser.org/review/api/gold-standard/",
            "secret-token-" * 4,
            opener=_FakeOpener(_live_payload(disputed)),
        )

    tampered = _live_payload(source_rows)
    tampered["dataset_revision"] = "0" * 64
    with pytest.raises(ingest.GoldSyncError, match="checksum"):
        ingest.fetch_live_gold(
            "https://variantbrowser.org/review/api/gold-standard/",
            "secret-token-" * 4,
            opener=_FakeOpener(tampered),
        )

    missing_identity = [dict(row) for row in source_rows]
    missing_identity[0]["decided_by_user_id"] = ""
    with pytest.raises(ingest.GoldSyncError, match="approving lead"):
        ingest.fetch_live_gold(
            "https://variantbrowser.org/review/api/gold-standard/",
            "secret-token-" * 4,
            opener=_FakeOpener(_live_payload(missing_identity)),
        )


def test_live_sync_requires_https_and_never_accepts_short_token():
    with pytest.raises(ingest.GoldSyncError, match="HTTPS"):
        ingest.fetch_live_gold(
            "http://variantbrowser.org/review/api/gold-standard/", "x" * 64
        )
    with pytest.raises(ingest.GoldSyncError, match="too short"):
        ingest.fetch_live_gold(
            "https://variantbrowser.org/review/api/gold-standard/", "short"
        )


def test_live_sync_redirect_handler_fails_explicitly():
    request = ingest.urllib.request.Request(
        "https://variantbrowser.org/review/api/gold-standard/"
    )

    with pytest.raises(ingest.urllib.error.HTTPError, match="blocked for security"):
        ingest._NoRedirect().redirect_request(
            request,
            None,
            302,
            "Found",
            {},
            "https://other.invalid/gold",
        )


def test_read_live_sync_state_handles_paths_with_spaces(tmp_path):
    path = tmp_path / "directory with spaces" / "review_gold.sqlite3"
    path.parent.mkdir()
    conn = sqlite3.connect(path)
    try:
        conn.execute(
            "CREATE TABLE review_gold_sync_state ("
            "singleton INTEGER PRIMARY KEY, dataset_revision TEXT)"
        )
        conn.execute("INSERT INTO review_gold_sync_state VALUES (1, 'abc123')")
        conn.commit()
    finally:
        conn.close()

    assert ingest.read_live_sync_state(path)["dataset_revision"] == "abc123"
