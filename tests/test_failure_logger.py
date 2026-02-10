"""Tests for the failure taxonomy logging system."""

import json
from pathlib import Path

import pytest

from pipeline.failure_logger import ExtractionFailure, FailureCode, FailureLogger


@pytest.fixture
def failure_log_path(tmp_path):
    """Provide a temp path for the failure log JSON."""
    return tmp_path / "extraction_failures.json"


@pytest.fixture
def logger(failure_log_path):
    """Provide a fresh FailureLogger writing to a temp file."""
    return FailureLogger(output_path=failure_log_path)


class TestFailureCode:
    def test_all_codes_exist(self):
        expected = [
            "PAYWALL", "NO_FULL_TEXT", "XML_PARSE_ERROR", "SUPPLEMENT_MISSING",
            "LLM_TIMEOUT", "NO_VARIANTS_FOUND", "LLM_INPUT_ERROR", "INCORRECT_JSON",
            "TEXT_TOO_SHORT", "METADATA_ERROR", "HTML_GARBAGE",
        ]
        actual = [c.value for c in FailureCode]
        assert set(expected) == set(actual)

    def test_code_values_match_names(self):
        for code in FailureCode:
            assert code.name == code.value


class TestExtractionFailure:
    def test_round_trip_dict(self):
        failure = ExtractionFailure(
            pmid="12345678",
            failure_code=FailureCode.PAYWALL,
            description="Could not access full text",
            model_used="gpt-4o",
            additional_info={"publisher": "Elsevier"},
        )
        d = failure.to_dict()
        assert d["failure_code"] == "PAYWALL"
        assert d["pmid"] == "12345678"

        restored = ExtractionFailure.from_dict(d)
        assert restored.failure_code == FailureCode.PAYWALL
        assert restored.pmid == failure.pmid
        assert restored.description == failure.description
        assert restored.model_used == failure.model_used
        assert restored.additional_info == failure.additional_info

    def test_timestamp_auto_populated(self):
        failure = ExtractionFailure(
            pmid="1", failure_code=FailureCode.LLM_TIMEOUT, description="timed out"
        )
        assert failure.timestamp  # non-empty


class TestFailureLogger:
    def test_log_single_failure(self, logger, failure_log_path):
        logger.log_failure("111", FailureCode.PAYWALL, "Behind paywall")
        assert len(logger.failures) == 1
        assert failure_log_path.exists()

        data = json.loads(failure_log_path.read_text())
        assert len(data) == 1
        assert data[0]["failure_code"] == "PAYWALL"

    def test_log_multiple_failures(self, logger):
        logger.log_failure("111", FailureCode.PAYWALL, "paywall")
        logger.log_failure("222", FailureCode.LLM_TIMEOUT, "timeout")
        logger.log_failure("333", FailureCode.PAYWALL, "paywall again")
        assert len(logger.failures) == 3

    def test_log_failure_with_kwargs(self, logger):
        f = logger.log_failure(
            "111",
            FailureCode.INCORRECT_JSON,
            "LLM returned invalid JSON",
            model_used="gpt-4o",
            raw_response="not json",
        )
        assert f.model_used == "gpt-4o"
        assert f.additional_info["raw_response"] == "not json"

    def test_get_failures_by_code(self, logger):
        logger.log_failure("111", FailureCode.PAYWALL, "paywall")
        logger.log_failure("222", FailureCode.LLM_TIMEOUT, "timeout")
        logger.log_failure("333", FailureCode.PAYWALL, "paywall 2")

        paywalls = logger.get_failures_by_code(FailureCode.PAYWALL)
        assert len(paywalls) == 2
        assert all(f.failure_code == FailureCode.PAYWALL for f in paywalls)

        timeouts = logger.get_failures_by_code(FailureCode.LLM_TIMEOUT)
        assert len(timeouts) == 1

    def test_get_failures_by_code_empty(self, logger):
        result = logger.get_failures_by_code(FailureCode.HTML_GARBAGE)
        assert result == []

    def test_get_failure_stats(self, logger):
        logger.log_failure("1", FailureCode.PAYWALL, "p")
        logger.log_failure("2", FailureCode.PAYWALL, "p")
        logger.log_failure("3", FailureCode.LLM_TIMEOUT, "t")
        logger.log_failure("4", FailureCode.TEXT_TOO_SHORT, "s")

        stats = logger.get_failure_stats()
        assert stats == {"PAYWALL": 2, "LLM_TIMEOUT": 1, "TEXT_TOO_SHORT": 1}

    def test_get_failure_stats_empty(self, logger):
        assert logger.get_failure_stats() == {}

    def test_get_failures_for_pmid(self, logger):
        logger.log_failure("111", FailureCode.PAYWALL, "paywall")
        logger.log_failure("111", FailureCode.LLM_TIMEOUT, "then timeout")
        logger.log_failure("222", FailureCode.PAYWALL, "different paper")

        pmid_failures = logger.get_failures_for_pmid("111")
        assert len(pmid_failures) == 2

    def test_persistence_across_instances(self, failure_log_path):
        logger1 = FailureLogger(output_path=failure_log_path)
        logger1.log_failure("111", FailureCode.PAYWALL, "paywall")
        logger1.log_failure("222", FailureCode.LLM_TIMEOUT, "timeout")

        logger2 = FailureLogger(output_path=failure_log_path)
        assert len(logger2.failures) == 2
        assert logger2.failures[0].pmid == "111"

    def test_handles_missing_file(self, tmp_path):
        path = tmp_path / "nonexistent" / "failures.json"
        fl = FailureLogger(output_path=path)
        assert fl.failures == []
        fl.log_failure("1", FailureCode.PAYWALL, "test")
        assert path.exists()

    def test_handles_corrupt_file(self, failure_log_path):
        failure_log_path.write_text("not valid json{{{")
        fl = FailureLogger(output_path=failure_log_path)
        assert fl.failures == []
