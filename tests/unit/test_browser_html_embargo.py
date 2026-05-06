"""Tests for Tier 3.5 embargo eligibility logic.

Per-strategy embargo policies + a shared comparator that respects the
approximate-month math used in production (30.44 days/month average).
"""

from __future__ import annotations

import datetime

from harvesting.browser_html.embargo import EmbargoChecker, _parse_pubdate_dict


def _today():
    return datetime.date(2026, 5, 6)


def test_no_embargo_policy_always_eligible():
    checker = EmbargoChecker(today=_today())
    eligible, reason = checker.is_eligible(pub_date=None, embargo_months=None)
    assert eligible is True
    assert "no embargo policy" in reason


def test_zero_embargo_always_eligible_even_without_pub_date():
    checker = EmbargoChecker(today=_today())
    eligible, reason = checker.is_eligible(pub_date=None, embargo_months=0)
    assert eligible is True
    assert "no embargo" in reason


def test_embargo_present_but_no_pub_date_blocks():
    checker = EmbargoChecker(today=_today())
    eligible, reason = checker.is_eligible(pub_date=None, embargo_months=12)
    assert eligible is False
    assert "pub_date unknown" in reason


def test_recent_paper_blocked_by_12mo_embargo():
    today = _today()
    pub_date = today - datetime.timedelta(days=30 * 6)  # 6 months ago
    checker = EmbargoChecker(today=today)
    eligible, _ = checker.is_eligible(pub_date=pub_date, embargo_months=12)
    assert eligible is False


def test_old_paper_passes_12mo_embargo():
    today = _today()
    pub_date = today - datetime.timedelta(days=30 * 14)  # 14 months ago
    checker = EmbargoChecker(today=today)
    eligible, _ = checker.is_eligible(pub_date=pub_date, embargo_months=12)
    assert eligible is True


def test_exactly_at_embargo_threshold_eligible():
    today = _today()
    # 12 months at 30.44 days = 365.28 days
    pub_date = today - datetime.timedelta(days=366)
    checker = EmbargoChecker(today=today)
    eligible, _ = checker.is_eligible(pub_date=pub_date, embargo_months=12)
    assert eligible is True


def test_parse_pubdate_year_only():
    pd = {"Year": "2014"}
    assert _parse_pubdate_dict(pd) == datetime.date(2014, 1, 1)


def test_parse_pubdate_year_and_month_short_string():
    pd = {"Year": "2014", "Month": "Jan"}
    assert _parse_pubdate_dict(pd) == datetime.date(2014, 1, 1)


def test_parse_pubdate_full():
    pd = {"Year": "2014", "Month": "06", "Day": "15"}
    assert _parse_pubdate_dict(pd) == datetime.date(2014, 6, 15)


def test_parse_pubdate_medline_fallback():
    pd = {"MedlineDate": "2014 Jan-Feb"}
    assert _parse_pubdate_dict(pd) == datetime.date(2014, 1, 1)


def test_parse_pubdate_missing_year_returns_none():
    assert _parse_pubdate_dict({}) is None
