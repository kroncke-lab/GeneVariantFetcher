"""Tests for utils.env_utils.get_env_int."""

from __future__ import annotations

from utils.env_utils import get_env_int


def test_get_env_int_parses_valid_value(monkeypatch):
    monkeypatch.setenv("GVF_TEST_INT", "42")
    assert get_env_int("GVF_TEST_INT", 7) == 42


def test_get_env_int_uses_default_when_unset(monkeypatch):
    monkeypatch.delenv("GVF_TEST_INT", raising=False)
    assert get_env_int("GVF_TEST_INT", 7) == 7


def test_get_env_int_falls_back_on_blank_or_garbage(monkeypatch):
    monkeypatch.setenv("GVF_TEST_INT", "   ")
    assert get_env_int("GVF_TEST_INT", 7) == 7
    monkeypatch.setenv("GVF_TEST_INT", "not-a-number")
    assert get_env_int("GVF_TEST_INT", 7) == 7
