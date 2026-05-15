"""Tests for browser-cookie domain coverage used by paywall recovery."""

from __future__ import annotations

from harvesting.browser_html.cookie_loader import DEFAULT_PUBLISHER_DOMAINS


def test_neurology_org_cookie_domain_is_loaded_for_wolters_kluwer_supplements():
    """PMID 19038855 supplementary files are served from neurology.org."""
    assert "neurology.org" in DEFAULT_PUBLISHER_DOMAINS


def test_mayo_clinic_proceedings_cookie_domain_is_loaded_for_elsevier_imprint():
    """PMID 14661677 lands on mayoclinicproceedings.org."""
    assert "mayoclinicproceedings.org" in DEFAULT_PUBLISHER_DOMAINS
