"""Tests for shared free-text output writer helpers."""

import csv

from harvesting.free_text_output_service import (
    publisher_api_fallback_source,
    source_from_free_text_flags,
    write_free_text_output,
)


def test_source_from_free_text_flags():
    elsevier_source = source_from_free_text_flags(
        used_elsevier_api=True, used_wiley_api=False
    )
    assert elsevier_source.success_marker == "ELSEVIER_API"
    assert elsevier_source.status_source == "elsevier-api"
    assert elsevier_source.source_tag == "[via Elsevier API]"

    wiley_source = source_from_free_text_flags(
        used_elsevier_api=False, used_wiley_api=True
    )
    assert wiley_source.success_marker == "WILEY_API"
    assert wiley_source.status_source == "wiley-api"
    assert wiley_source.source_tag == "[via Wiley API]"

    publisher_source = source_from_free_text_flags(
        used_elsevier_api=False, used_wiley_api=False
    )
    assert publisher_source.success_marker == "PUBLISHER_FREE"
    assert publisher_source.status_source == "publisher-free"
    assert publisher_source.source_tag == "[from publisher]"


def test_write_free_text_output_writes_file_log_and_status(tmp_path):
    success_log = tmp_path / "success.csv"
    success_log.write_text("PMID,PMCID,Supplements_Downloaded\n", encoding="utf-8")

    status_calls = []
    source = publisher_api_fallback_source()
    output_file, unified_content = write_free_text_output(
        output_dir=tmp_path,
        success_log=success_log,
        pmid="123456",
        main_markdown="main",
        supplement_markdown="\n\nsupp",
        downloaded_count=2,
        source=source,
        write_pmid_status=lambda pmid, status, details: status_calls.append(
            (pmid, status, details)
        ),
    )

    assert output_file.exists()
    assert output_file.read_text(encoding="utf-8") == "main\n\nsupp"
    assert unified_content == "main\n\nsupp"

    with open(success_log, "r", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[-1] == ["123456", "publisher-api", "2"]

    assert len(status_calls) == 1
    pmid, status, details = status_calls[0]
    assert pmid == "123456"
    assert status == "extracted"
    assert details["variant_count"] == 0
    assert details["source"] == "publisher-api-fallback"
