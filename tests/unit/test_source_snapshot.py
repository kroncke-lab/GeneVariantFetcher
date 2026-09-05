"""Frozen comparison sources must never silently turn into live acquisition."""

from unittest.mock import Mock

import pytest

from pipeline.source_snapshot import freeze_sources, stage_frozen_sources


def snapshot(tmp_path):
    source = tmp_path / "original" / "123_FULL_CONTEXT.md"
    source.parent.mkdir()
    source.write_text(
        "<!-- GVF_SUPPLEMENT_SURFACE_STATUS: unavailable -->\nAbstract-only source."
    )
    (source.parent / "123_extraction.json").write_text('{"variants": ["forbidden"]}')
    return freeze_sources(
        [{"gene": "SCN5A", "pmid": "123", "source": str(source)}], tmp_path / "frozen"
    )


def test_preserves_incomplete_source_and_excludes_predictions(tmp_path):
    manifest = snapshot(tmp_path)
    target = tmp_path / "run"
    assert stage_frozen_sources(manifest, "SCN5A", ["123"], target) == {"123"}
    assert "unavailable" in (target / "123_FULL_CONTEXT.md").read_text()
    assert not (target / "123_extraction.json").exists()


def test_modified_or_missing_source_cannot_fall_back(tmp_path):
    manifest = snapshot(tmp_path)
    with pytest.raises(ValueError, match="missing frozen source"):
        stage_frozen_sources(manifest, "SCN5A", ["999"], tmp_path / "run")
    (manifest.parent / "SCN5A/123/123_FULL_CONTEXT.md").write_text("different")
    with pytest.raises(ValueError, match="frozen source changed"):
        stage_frozen_sources(manifest, "SCN5A", ["123"], tmp_path / "run")


def test_frozen_download_never_calls_harvester_or_prior_cache(tmp_path, monkeypatch):
    import harvesting
    import pipeline.steps as steps

    manifest = snapshot(tmp_path)
    harvester = Mock()
    monkeypatch.setattr(harvesting, "PMCHarvester", Mock(return_value=harvester))
    monkeypatch.setattr(
        steps,
        "_consolidate_from_corpus",
        Mock(side_effect=AssertionError("live corpus")),
    )
    monkeypatch.setattr(
        steps,
        "_consolidate_prior_downloads",
        Mock(side_effect=AssertionError("cached predictions")),
    )
    monkeypatch.setenv("GVF_FROZEN_SOURCE_MANIFEST", str(manifest))
    monkeypatch.setenv("GVF_DISABLE_GOLD_DERIVED_ALIASES", "1")
    result = steps.download_fulltext(["123"], tmp_path / "run", "SCN5A")
    harvester.harvest.assert_not_called()
    assert result.stats["attempted"] == 1
