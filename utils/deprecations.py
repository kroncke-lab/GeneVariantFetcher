"""Small, consistent deprecation notices for legacy importable helpers."""

from __future__ import annotations

import warnings


def warn_deprecated(name: str, replacement: str) -> None:
    """Warn callers while preserving one compatibility window."""

    warnings.warn(
        f"{name} is deprecated; use {replacement}.",
        DeprecationWarning,
        stacklevel=3,
    )
