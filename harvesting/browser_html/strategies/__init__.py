"""Per-publisher strategies for Tier 3.5 browser HTML harvesting.

To add a new publisher:

1. Create ``my_publisher.py`` in this directory
2. Define a class inheriting from ``PublisherStrategy``
3. Set ``NAME``, ``DOI_PREFIXES``, ``DOMAINS``, ``EMBARGO_MONTHS``
4. Implement ``fetch(page, ctx) -> FetchResult``

The registry below auto-imports every module in this package and registers
any concrete subclass of ``PublisherStrategy``. No manual registration is
needed — drop the file in and it's available.
"""

from __future__ import annotations

import importlib
import logging
import pkgutil
from typing import Dict, List, Optional, Type

from ..base import PublisherStrategy

logger = logging.getLogger(__name__)


_REGISTRY: Dict[str, Type[PublisherStrategy]] = {}
_LOADED = False


def _autodiscover() -> None:
    """Import every sibling module so subclasses register themselves."""
    global _LOADED
    if _LOADED:
        return
    package = __name__
    package_path = __path__  # type: ignore[name-defined]
    for _finder, modname, _ispkg in pkgutil.iter_modules(package_path):
        if modname.startswith("_"):
            continue
        try:
            importlib.import_module(f"{package}.{modname}")
        except Exception as e:
            logger.warning(f"Strategy module {modname} failed to import: {e}")
    _LOADED = True


def register(cls: Type[PublisherStrategy]) -> Type[PublisherStrategy]:
    """Decorator to register a strategy class under its NAME."""
    name = getattr(cls, "NAME", "").strip().lower()
    if not name or name == "base":
        raise ValueError(
            f"Strategy {cls.__name__} must set a non-empty NAME (got {name!r})"
        )
    if name in _REGISTRY:
        # Last writer wins; warn so collisions are visible.
        logger.warning(
            "Strategy NAME collision for %s: replacing %s with %s",
            name,
            _REGISTRY[name].__name__,
            cls.__name__,
        )
    _REGISTRY[name] = cls
    return cls


def all_strategies() -> List[PublisherStrategy]:
    """Instantiate every registered strategy in registration order."""
    _autodiscover()
    return [cls() for cls in _REGISTRY.values()]


def find_strategy(
    doi: str,
    url: str = "",
    allowlist: Optional[List[str]] = None,
) -> Optional[PublisherStrategy]:
    """Resolve the first matching strategy for a DOI/URL.

    The generic strategy is excluded from the prefix-match round and only
    selected when no other strategy claims the DOI.
    """
    _autodiscover()

    allowed = {n.lower() for n in allowlist} if allowlist is not None else None

    generic: Optional[PublisherStrategy] = None
    for name, cls in _REGISTRY.items():
        if allowed is not None and name not in allowed:
            continue
        if name == "generic":
            generic = cls()
            continue
        instance = cls()
        if instance.matches(doi, url):
            return instance

    return generic


def registered_names() -> List[str]:
    """List registered strategy NAMEs (after autodiscovery)."""
    _autodiscover()
    return list(_REGISTRY.keys())
