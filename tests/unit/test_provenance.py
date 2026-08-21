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
        "extractor_hash_algorithm",
        "extractor_code_sha256",
        "extractor_config_sha256",
        "prompt_extractor_sha256",
        "legacy_prompt_extractor_sha256",
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


def test_python_token_hash_ignores_comments_and_cosmetic_whitespace(tmp_path):
    path = tmp_path / "extractor.py"
    path.write_text("value = 1  # original comment\n", encoding="utf-8")
    before = provenance.hash_python_tokens(["extractor.py"], root=tmp_path)

    path.write_text("value=1 # rewritten comment\n", encoding="utf-8")
    cosmetic = provenance.hash_python_tokens(["extractor.py"], root=tmp_path)

    path.write_text("value = 2\n", encoding="utf-8")
    changed = provenance.hash_python_tokens(["extractor.py"], root=tmp_path)

    assert before == cosmetic
    assert changed != before


def test_resolved_config_changes_only_config_and_combined_hash():
    first = {"tier3_model": "model-a", "tier3_reasoning_effort": "high"}
    second = {"tier3_model": "model-b", "tier3_reasoning_effort": "high"}

    assert provenance.extractor_config_hash(first) != provenance.extractor_config_hash(
        second
    )
    assert provenance.prompt_extractor_hash(first) != provenance.prompt_extractor_hash(
        second
    )


def test_resolved_routing_is_a_closed_scientific_allowlist():
    routing = provenance.resolved_model_routing()
    allowed = {
        "model_provider",
        *(key for key, _resolver in provenance.EXTRACTOR_MODEL_RESOLVERS),
        *provenance.EXTRACTOR_CONFIG_FIELDS,
        *provenance.EXTRACTOR_RUNTIME_CONSTANTS,
    }

    assert set(routing) <= allowed
    assert {"output_dir", "email", "pmids", "timestamp"}.isdisjoint(routing)


def test_provenance_records_extraction_runtime_constants_not_false_knobs():
    routing = provenance.resolved_model_routing()

    assert {
        key: routing[key] for key in provenance.EXTRACTOR_RUNTIME_CONSTANTS
    } == provenance.EXTRACTOR_RUNTIME_CONSTANTS
    assert {
        "extraction_max_chars",
        "scanner_merge_confidence",
        "scanner_max_hints",
    }.isdisjoint(routing)


def test_legacy_digest_retains_the_pre_v2_raw_file_set():
    assert provenance.legacy_prompt_extractor_hash() == provenance.hash_files(
        provenance.LEGACY_PROMPT_EXTRACTOR_FILES
    )


def test_provenance_tracks_post_extraction_scientific_mutators():
    required = {
        "pipeline/claim_verifier.py",
        "pipeline/count_recovery.py",
        "pipeline/count_repair.py",
        "pipeline/trust_gate.py",
        "harvesting/migrate_to_sqlite.py",
        "utils/llm_utils.py",
        "utils/source_layers.py",
    }
    assert required <= set(provenance.EXTRACTOR_CODE_FILES)
    assert all(
        path.split("/", 1)[0]
        in {"cli", "config", "gene_literature", "harvesting", "pipeline", "utils"}
        for path in provenance.EXTRACTOR_CODE_FILES
    )
    assert "config/settings.py" not in provenance.EXTRACTOR_CODE_FILES
    assert "utils/provenance.py" not in provenance.EXTRACTOR_CODE_FILES


def test_collect_provenance_records_missing_hash_inputs():
    prov = provenance.collect_provenance()
    assert "extractor_code_files_missing" in prov
    assert "prompt_extractor_files_missing" in prov
    # Every expected extractor file exists in the repo, so none are missing.
    assert prov["prompt_extractor_files_missing"] == []
    assert prov["extractor_code_files_missing"] == []
