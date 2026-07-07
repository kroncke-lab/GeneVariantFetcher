"""Small helpers for reading configuration from environment variables."""

from __future__ import annotations

import os


def get_env_int(name: str, default: int) -> int:
    """Return the int value of env var ``name``, or ``default``.

    Falls back to ``default`` when the variable is unset, blank, or not a
    valid integer, so a typo like ``FOO=abc`` cannot crash module import.
    """
    raw = os.environ.get(name)
    if raw is None or not raw.strip():
        return default
    try:
        return int(raw)
    except ValueError:
        return default
