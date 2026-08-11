"""Hermeticity guards for the offline unit suite.

`pytest tests/unit` must behave identically in the main checkout, in a
`.claude/worktrees/<name>` side worktree, and in CI. Several helpers fall back
to *guessing* a local-data path when nothing is configured — most notably
`utils.gene_metadata.default_variantfeatures_db_path`, whose sibling-checkout
guess resolves to a real multi-gigabyte VariantFeatures database on a
developer's machine and to nothing anywhere else. That one path difference is
enough to route the suite through a different code path depending on where it
runs, and it is what previously made a unit test appear to hang for over ten
minutes locally while passing in seconds in CI.

So the offline suite switches implicit local-data discovery off and strips the
env vars that point at on-disk local data, which a developer's shell or `.env`
(loaded in `tests/conftest.py`) may have set. Tests that genuinely want a
database build their own fixture and set `VARIANTFEATURES_DB` via monkeypatch —
an explicitly configured path still wins over the guard. A test that wants the
guessing layer back asks for the `allow_local_data_discovery` fixture.
"""

from __future__ import annotations

import os

import pytest

from utils.env_utils import DISABLE_LOCAL_DATA_ENV
from utils.gene_metadata import clear_gene_metadata_cache

# Env vars that point production code at on-disk local data or a sibling repo.
LOCAL_DATA_ENV_VARS = (
    "VARIANTFEATURES_DB",
    "VARIANT_FEATURES_DB",
    "GVF_CORPUS_DIR",
    "GVF_REVIEW_REPO",
    "VARIANT_BROWSER_DIR",
    "GVF_EZPROXY_PROFILE_DIR",
)


def pytest_configure(config):
    """Disable implicit local-data discovery for the whole offline suite.

    Set here rather than only in a fixture so module-level constants evaluated
    during collection see the guard too. `setdefault` leaves an explicit
    `GVF_DISABLE_LOCAL_DATA=0` alone, which is the escape hatch for reproducing
    the old discover-whatever-is-on-disk behaviour.
    """
    os.environ.setdefault(DISABLE_LOCAL_DATA_ENV, "1")


@pytest.fixture(autouse=True)
def hermetic_local_data(monkeypatch):
    """Keep each offline test independent of a developer's on-disk local data."""
    for name in LOCAL_DATA_ENV_VARS:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv(
        DISABLE_LOCAL_DATA_ENV, os.environ.get(DISABLE_LOCAL_DATA_ENV) or "1"
    )
    # get_gene_metadata is lru_cached on (gene, db path), so a test that points
    # it at a fixture database would otherwise leak those entries to its
    # successors — the same cross-test coupling the env guard removes.
    clear_gene_metadata_cache()
    yield
    clear_gene_metadata_cache()


@pytest.fixture
def allow_local_data_discovery(monkeypatch):
    """Opt back into implicit local-data discovery for one test.

    Only for tests that assert on the discovery logic itself; point the
    discovery at fixture paths rather than at whatever the developer happens to
    have checked out.
    """
    monkeypatch.setenv(DISABLE_LOCAL_DATA_ENV, "0")
