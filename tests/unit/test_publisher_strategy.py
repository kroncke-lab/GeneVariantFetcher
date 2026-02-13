"""Tests for publisher fallback strategy planning."""

from harvesting.publisher_strategy import PublisherAttempt, build_publisher_attempt_plan


def _noop_fetch(doi: str, pmid: str):
    return None, None


def test_build_plan_prefers_doi_matches():
    attempts = [
        PublisherAttempt("Elsevier", should_try=False, try_fetch=_noop_fetch),
        PublisherAttempt("Springer", should_try=True, try_fetch=_noop_fetch),
        PublisherAttempt("Wiley", should_try=True, try_fetch=_noop_fetch),
    ]

    plan = build_publisher_attempt_plan("10.1007/test", attempts)

    assert [attempt.name for attempt in plan] == ["Springer", "Wiley", "Elsevier"]


def test_build_plan_keeps_unmatched_when_needed():
    attempts = [
        PublisherAttempt("Elsevier", should_try=False, try_fetch=_noop_fetch),
        PublisherAttempt("Wiley", should_try=False, try_fetch=_noop_fetch),
    ]

    plan = build_publisher_attempt_plan("10.0000/unknown", attempts)

    assert [attempt.name for attempt in plan] == ["Elsevier", "Wiley"]
