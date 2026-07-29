"""Fail-closed guard for the repo paths that are backed by external storage.

On Brett's workstation the repo path ``corpus/`` is an absolute symlink to an
external APFS volume (see CLAUDE.md > "Operating Shape"); the multi-GB fetched
source lives there, not on the internal disk. Any write that would otherwise
``mkdir(parents=True)`` its way into ``corpus/`` has to answer one question
first: is the link actually present?

A *dangling* link already fails on its own — ``mkdir`` refuses to follow one. A
*missing* link does not: ``mkdir`` cheerfully creates a real local ``corpus/``,
and the run then builds a second corpus on the internal disk that no later run
reads, re-fetching paywalled source it already has. That split is silent, slow
to notice, and expensive to undo, so the guard refuses the write instead.

This is the write-side companion to ``local_data_discovery_disabled`` in
``utils.env_utils``: that one stops the offline suite from *reading* whatever is
on disk, this one stops a real run from *creating* storage in the wrong place.
"""

from __future__ import annotations

import os
from pathlib import Path

from utils.env_utils import get_env_bool

ALLOW_LOCAL_ENV = "GVF_ALLOW_LOCAL_CORPUS"

# Repo-relative paths that belong on external storage and must never be created.
EXTERNAL_PATHS = ("corpus",)

# Both anchors for the same tree: the unresolved one matches how callers build
# their paths, the resolved one matches the same location reached through a
# symlinked checkout.
_ANCHORS = tuple(
    dict.fromkeys(
        (
            Path(__file__).parent.parent,
            Path(__file__).resolve().parent.parent,
        )
    )
)


class LocalStorageError(Exception):
    """Raised when a write would create local storage that belongs on an external volume."""


def local_storage_allowed() -> bool:
    """True when the operator has opted into plain local storage for ``corpus/``.

    The escape hatch exists for a machine that legitimately has no external
    volume (CI, a collaborator's laptop, a throwaway checkout).
    """
    return get_env_bool(ALLOW_LOCAL_ENV, False)


def external_path_state(name: str) -> str:
    """Classify repo-relative path ``name`` as ``linked``, ``dangling``, or ``absent``.

    ``linked`` covers both a symlink whose target is reachable and a real local
    directory — either way the storage exists and a caller may write into it.
    """
    for anchor in _ANCHORS:
        path = anchor / name
        if path.is_dir():
            return "linked"
        if path.is_symlink():
            return "dangling"
    return "absent"


def _guarded_root(target: Path) -> str | None:
    """Return the guarded repo-relative name containing ``target``, if any.

    ``abspath`` normalizes without resolving symlinks, so a path under the repo's
    ``corpus/`` link keeps its repo-relative shape and still matches. Callers
    must therefore guard *before* calling ``Path.resolve()``, which would rewrite
    the path to its on-volume target and slip past this check.
    """
    try:
        absolute = Path(os.path.abspath(target))
    except (OSError, ValueError):
        return None
    for name in EXTERNAL_PATHS:
        for anchor in _ANCHORS:
            if absolute.is_relative_to(anchor / name):
                return name
    return None


def require_external_storage(target: Path) -> None:
    """Refuse to create ``target`` when its external-storage link is not in place.

    Paths outside the guarded trees — a pytest ``tmp_path``, an explicit
    ``GVF_CORPUS_DIR`` somewhere else, an already-resolved path on the volume
    itself — are not this guard's business and pass straight through.
    """
    name = _guarded_root(target)
    if name is None:
        return

    state = external_path_state(name)
    if state == "linked" or local_storage_allowed():
        return

    found = {
        "absent": f"no {name!r} entry in the repo root at all",
        "dangling": f"{name!r} is a symlink whose target is unreachable",
    }[state]
    raise LocalStorageError(
        f"refusing to create {target} — {found}.\n"
        f"\n"
        f"{name}/ is expected to be an absolute symlink to external storage, so "
        f"creating it here would put a second copy on the internal disk that no "
        f"later run reads — re-fetching source you already have.\n"
        f"\n"
        f"Mount the APFS volume 'Ezekers' at /Volumes/Ezekers (and recreate the "
        f"{name}/ symlink if this is a fresh checkout), then retry. See "
        f"CLAUDE.md > 'Operating Shape'.\n"
        f"To use plain local storage on purpose, set {ALLOW_LOCAL_ENV}=1."
    )
