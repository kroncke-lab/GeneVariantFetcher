"""Unit coverage for utils.resilience.CircuitBreaker half-open transitions.

The pre-existing suites never exercised the half-open state: tests/test_circuit_breaker.py
covers an unrelated ExpertExtractor input-quality guard, and the integration test only
drives closed -> open. These tests pin the half-open trial behavior, including the
regression fixed here: a failed half-open trial must re-open the breaker instead of
leaving it half-open (which would keep letting traffic through to a broken dependency).
"""

from utils.resilience import CircuitBreaker


def _open_then_half_open(max_failures: int = 3) -> CircuitBreaker:
    """Return a breaker that has tripped open and been promoted to half-open."""
    cb = CircuitBreaker("test", max_failures=max_failures, reset_timeout=0)
    for _ in range(max_failures):
        cb.record_failure()
    assert cb.state == "open"
    # reset_timeout=0 => the next is_open() check promotes open -> half-open.
    assert cb.is_open() is False
    assert cb.state == "half-open"
    return cb


def test_failed_half_open_trial_reopens():
    """A failing trial while half-open re-opens the breaker (regression)."""
    cb = _open_then_half_open()

    cb.record_failure()

    # Must be back to "open", not lingering "half-open" (the old bug).
    assert cb.state == "open"


def test_half_open_closes_after_three_successes():
    """Three consecutive successful trials close the breaker and reset counts."""
    cb = _open_then_half_open()

    cb.record_success()
    cb.record_success()
    assert cb.state == "half-open"  # not yet three

    cb.record_success()
    assert cb.state == "closed"
    assert cb.failure_count == 0
    assert cb.consecutive_successes == 0


def test_closed_breaker_opens_after_max_failures():
    """Baseline closed -> open path is unchanged by the half-open fix."""
    cb = CircuitBreaker("test", max_failures=2, reset_timeout=60)
    assert cb.state == "closed"
    cb.record_failure()
    assert cb.state == "closed"
    cb.record_failure()
    assert cb.state == "open"
    assert cb.is_open() is True
