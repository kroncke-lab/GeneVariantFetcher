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


def test_hash_files_none_when_absent(tmp_path):
    assert provenance.hash_files(["nope.txt"], root=tmp_path) is None


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
