"""Tests for harvesting persistence helpers."""

import json

from harvesting.persistence import (
    append_paywalled_entry,
    append_success_entry,
    initialize_harvest_logs,
    write_pmid_status_file,
)


def test_initialize_harvest_logs_writes_headers(tmp_path):
    paywalled_log = tmp_path / "paywalled_missing.csv"
    success_log = tmp_path / "successful_downloads.csv"

    initialize_harvest_logs(paywalled_log, success_log)

    paywalled_header = paywalled_log.read_text(encoding="utf-8").splitlines()[0]
    success_header = success_log.read_text(encoding="utf-8").splitlines()[0]
    assert paywalled_header.startswith("PMID,Reason,URL")
    assert success_header == "PMID,PMCID,Supplements_Downloaded"


def test_append_entries_and_status_file(tmp_path):
    paywalled_log = tmp_path / "paywalled_missing.csv"
    success_log = tmp_path / "successful_downloads.csv"
    initialize_harvest_logs(paywalled_log, success_log)

    append_paywalled_entry(paywalled_log, "123", "reason", "http://example.org")
    append_success_entry(success_log, "456", "PMC456", 2)

    paywalled_rows = paywalled_log.read_text(encoding="utf-8").splitlines()
    success_rows = success_log.read_text(encoding="utf-8").splitlines()
    assert len(paywalled_rows) == 2
    assert "123,reason,http://example.org" in paywalled_rows[1]
    assert success_rows[1] == "456,PMC456,2"

    status_file = write_pmid_status_file(
        output_dir=tmp_path,
        pmid="456",
        status="extracted",
        details={"variant_count": 4, "source": "pmc"},
    )
    payload = json.loads(status_file.read_text(encoding="utf-8"))
    assert payload["pmid"] == "456"
    assert payload["status"] == "extracted"
    assert payload["variant_count"] == 4
    assert payload["source"] == "pmc"
