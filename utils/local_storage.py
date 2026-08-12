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
    """Classify a guarded repo path as linked, local, dangling, or absent."""
    for anchor in _ANCHORS:
        path = anchor / name
        if path.is_symlink():
            return "linked" if path.is_dir() else "dangling"
        if path.is_dir():
            return "local"
    return "absent"


def external_link_target(name: str) -> Path | None:
    """The symlink target of repo-relative ``name``, whether or not it resolves.

    Read with ``os.readlink`` rather than ``resolve()`` so an unmounted volume
    still reports where it was pointed — that string is what tells an operator
    which drive to attach.
    """
    for anchor in _ANCHORS:
        path = anchor / name
        if path.is_symlink():
            try:
                return Path(os.readlink(path))
            except OSError:
                return None
    return None


def external_volume_name(name: str) -> str | None:
    """The macOS volume a repo-relative link points into, e.g. ``Ezekers``.

    Derived from the link itself rather than hardcoded, so this says the right
    thing on a machine whose drive is named something else — and says nothing at
    all on a collaborator's checkout that has no link to read.
    """
    target = external_link_target(name)
    if target is None:
        return None
    parts = target.parts
    try:
        return parts[parts.index("Volumes") + 1]
    except (ValueError, IndexError):
        return None


def _remedy_lines(name: str) -> list:
    """Guidance for restoring ``name``, matched to what the checkout looks like.

    A dangling link means a known drive is detached — say which one. Nothing to
    read means a fresh clone, where telling someone to mount a volume they have
    never heard of is worse than saying nothing: they need the setup path.
    """
    volume = external_volume_name(name)
    if volume:
        return [
            f"Attach the volume {volume!r} and mount it at /Volumes/{volume}, then retry.",
            f"The {name}/ symlink points at {external_link_target(name)}.",
        ]
    return [
        f"This checkout has no {name}/ at all, which is the normal state of a "
        f"fresh clone. The corpus is a cache of fetched source, so it is safe to "
        f"start empty — it just has to be somewhere deliberate.",
        "Set one up as either an absolute symlink to your own storage:",
        f"    ln -s /path/to/your/storage {name}",
        f"or a plain local directory, with {ALLOW_LOCAL_ENV}=1 set:",
        f"    mkdir {name}",
    ]


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
        "absent": f"there is no {name}/ in the repo root",
        "dangling": f"{name}/ is a symlink whose target is unreachable",
        "local": f"{name}/ is a real local directory rather than an external link",
    }[state]
    lines = [
        f"cannot use {name}/ — {found}.",
        "",
        f"Refusing to create it, because {target} on the internal disk would be a "
        f"second copy that no later run reads — re-fetching source you already have.",
        "",
        *_remedy_lines(name),
    ]
    raise LocalStorageError("\n".join(lines))
