"""Small helpers for reading configuration from environment variables."""

from __future__ import annotations

import os

DISABLE_LOCAL_DATA_ENV = "GVF_DISABLE_LOCAL_DATA"

_TRUTHY = frozenset({"1", "true", "yes", "on"})
_FALSEY = frozenset({"0", "false", "no", "off"})


def get_env_bool(name: str, default: bool = False) -> bool:
    """Return the boolean value of env var ``name``, or ``default``.

    Accepts ``1/true/yes/on`` and ``0/false/no/off`` (case-insensitive). An
    unset, blank, or unrecognized value falls back to ``default`` so a typo
    cannot crash module import.
    """
    raw = os.environ.get(name)
    if raw is None:
        return default
    value = raw.strip().lower()
    if value in _TRUTHY:
        return True
    if value in _FALSEY:
        return False
    return default


def local_data_discovery_disabled() -> bool:
    """True when *implicit* discovery of on-disk local data is switched off.

    Several helpers fall back to guessing a path when nothing is configured —
    a sibling ``../variantFeatures`` checkout, ``<repo>/corpus``, and so on.
    Those guesses resolve in a developer's main checkout but not in a side
    worktree or in CI, so the same code silently takes different paths
    depending on where it runs. The offline unit suite sets
    ``GVF_DISABLE_LOCAL_DATA=1`` to make that fallback layer inert.

    This only suppresses *guessing*. An explicitly configured path — a
    function argument or an env var such as ``VARIANTFEATURES_DB`` /
    ``GVF_CORPUS_DIR`` — is still honoured, so a test that builds its own
    fixture database keeps working.
    """
    return get_env_bool(DISABLE_LOCAL_DATA_ENV, False)


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
