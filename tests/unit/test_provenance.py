"""Run provenance is captured and is content-addressed + fail-safe."""

from utils import provenance


def test_hash_files_is_stable_and_content_addressed(tmp_path):
    (tmp_path / "a.txt").write_text("alpha")
    (tmp_path / "b.txt").write_text("beta")

    h1 = provenance.hash_files(["a.txt", "b.txt"], root=tmp_path)
    h2 = provenance.hash_files(["a.txt", "b.txt"], root=tmp_path)
    assert h1 == h2 and h1 is not None

    (tmp_path / "b.txt").write_text("beta-changed")
    h3 = provenance.hash_files(["a.txt", "b.txt"], root=tmp_path)
    assert h3 != h1


def test_hash_files_folds_missing_as_sentinel(tmp_path):
    """A missing listed file is folded in (not skipped), so its absence changes
    the hash and can't silently vanish. Only an empty list returns None."""
    (tmp_path / "a.txt").write_text("alpha")
    present = provenance.hash_files(["a.txt"], root=tmp_path)
    with_missing = provenance.hash_files(["a.txt", "nope.txt"], root=tmp_path)
    assert with_missing is not None
    assert with_missing != present
    assert provenance.hash_files([], root=tmp_path) is None


def test_collect_provenance_has_expected_keys_and_never_raises():
    prov = provenance.collect_provenance()
    for key in (
        "git_sha",
        "git_dirty",
        "prompt_extractor_sha256",
        "dependency_lock_sha256",
        "model_routing",
    ):
        assert key in prov
    # model_routing is a dict (possibly empty if settings unavailable).
    assert isinstance(prov["model_routing"], dict)


def test_prompt_extractor_hash_tracks_real_files():
    # These files exist in the repo, so the hash must be populated.
    assert provenance.prompt_extractor_hash() is not None
    assert provenance.dependency_lock_hash() is not None


def test_collect_provenance_records_missing_hash_inputs():
    prov = provenance.collect_provenance()
    assert "prompt_extractor_files_missing" in prov
    # Every expected extractor file exists in the repo, so none are missing.
    assert prov["prompt_extractor_files_missing"] == []
