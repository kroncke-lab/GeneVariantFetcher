"""Nested-archive expansion: read the inner bundle, exactly once, safely.

A publisher bundle routinely contains the journal's own supplement archive. The
inner archive used to be written to disk and then skipped by every converter,
so its tables were never read by any extraction route.

All fixtures are built in-test. Nothing here depends on a real paper, a
publisher's directory layout, or a gene.
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path

import pytest

from harvesting.format_converters import (
    ARCHIVE_MAX_DEPTH,
    FormatConverter,
)


def _zip_bytes(members: dict[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, payload in members.items():
            zf.writestr(name, payload)
    return buffer.getvalue()


def _write_zip(path: Path, members: dict[str, bytes]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_zip_bytes(members))
    return path


def test_inner_archive_contents_are_converted(tmp_path):
    """The defect: a zip inside a zip reached no converter at all."""
    inner = _zip_bytes({"inner_table.txt": b"Variant\tCarriers\nc.100A>G\t4\n"})
    outer = _write_zip(
        tmp_path / "bundle.zip",
        {"figure.jpg": b"\xff\xd8\xff\xe0binary", "supplement-s001.zip": inner},
    )

    paths, markdown = FormatConverter().extract_zip_supplement(
        outer, dest_dir=tmp_path / "out"
    )

    assert "c.100A>G" in markdown, "inner archive contents never reached markdown"
    assert any(p.name == "inner_table.txt" for p in paths)


def test_inner_member_is_converted_exactly_once(tmp_path):
    """Double conversion would double every table the inner archive holds."""
    inner = _zip_bytes({"inner_table.txt": b"UNIQUE_ROW_MARKER\n"})
    outer = _write_zip(tmp_path / "bundle.zip", {"supplement-s001.zip": inner})

    _paths, markdown = FormatConverter().extract_zip_supplement(
        outer, dest_dir=tmp_path / "out"
    )

    assert markdown.count("UNIQUE_ROW_MARKER") == 1


def test_depth_limit_stops_descent(tmp_path):
    """A pathologically nested archive must terminate, not recurse forever."""
    payload = _zip_bytes({"deep.txt": b"DEEPEST_MARKER\n"})
    for _ in range(ARCHIVE_MAX_DEPTH + 2):
        payload = _zip_bytes({"nested.zip": payload})
    outer = _write_zip(tmp_path / "bomb.zip", {"nested.zip": payload})

    _paths, markdown = FormatConverter().extract_zip_supplement(
        outer, dest_dir=tmp_path / "out"
    )

    assert "DEEPEST_MARKER" not in markdown


def test_self_referential_archive_terminates(tmp_path):
    """Identical members are keyed by digest, so a cycle converts once."""
    inner = _zip_bytes({"leaf.txt": b"LEAF\n"})
    outer = _write_zip(
        tmp_path / "bundle.zip", {"a/copy.zip": inner, "b/copy.zip": inner}
    )

    _paths, markdown = FormatConverter().extract_zip_supplement(
        outer, dest_dir=tmp_path / "out"
    )

    assert markdown.count("LEAF") == 1, "duplicate archive members re-converted"


def test_zip_slip_member_is_not_written_outside_the_paper_directory(tmp_path):
    dest = tmp_path / "out"
    outer = _write_zip(
        tmp_path / "evil.zip",
        {"../escaped.txt": b"ESCAPED", "ok.txt": b"SAFE_CONTENT"},
    )

    paths, markdown = FormatConverter().extract_zip_supplement(outer, dest_dir=dest)

    assert not (tmp_path / "escaped.txt").exists()
    assert all(dest.resolve() in p.resolve().parents for p in paths)
    assert "SAFE_CONTENT" in markdown


def test_macos_metadata_members_are_not_treated_as_source(tmp_path):
    outer = _write_zip(
        tmp_path / "bundle.zip",
        {
            "__MACOSX/._table.txt": b"APPLEDOUBLE_JUNK",
            ".DS_Store": b"DSSTORE_JUNK",
            "table.txt": b"REAL_CONTENT",
        },
    )

    _paths, markdown = FormatConverter().extract_zip_supplement(
        outer, dest_dir=tmp_path / "out"
    )

    assert "APPLEDOUBLE_JUNK" not in markdown
    assert "DSSTORE_JUNK" not in markdown
    assert "REAL_CONTENT" in markdown


def test_oversized_declared_member_is_skipped(tmp_path, monkeypatch):
    """A member claiming an implausible expansion is refused, not expanded."""
    import harvesting.format_converters as fc

    monkeypatch.setattr(fc, "ARCHIVE_MAX_MEMBER_UNCOMPRESSED_BYTES", 16)
    outer = _write_zip(
        tmp_path / "bundle.zip",
        {"huge.txt": b"X" * 4096, "small.txt": b"TINY"},
    )

    _paths, markdown = FormatConverter().extract_zip_supplement(
        outer, dest_dir=tmp_path / "out"
    )

    assert "X" * 100 not in markdown
    assert "TINY" in markdown


def test_fold_expands_an_archive_harvest_never_unpacked(tmp_path):
    """The fold-path half of the defect.

    ``.zip`` is excluded from the fold's convertible suffixes so a harvest-time
    expansion is not double-counted. That exclusion also skipped archives
    harvest never expanded at all. The discriminator is whether the extraction
    directory already exists.
    """
    from harvesting.supplement_fold import fold_supplements_into_full_context

    pmid = "00000001"  # synthetic identifier, not a real record
    harvest = tmp_path / "harvest"
    supp = harvest / f"{pmid}_supplements"
    supp.mkdir(parents=True)
    (harvest / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# Body\n\nMain article text.\n", encoding="utf-8"
    )
    _write_zip(supp / "never_unpacked.zip", {"table_s1.txt": b"FOLDED_INNER_ROW\n"})

    result = fold_supplements_into_full_context(pmid, harvest)

    assert result is not None
    assert "FOLDED_INNER_ROW" in result.read_text(encoding="utf-8")


def test_fold_does_not_re_expand_an_already_extracted_archive(tmp_path):
    """An archive whose members are already on disk is left alone."""
    from harvesting.supplement_fold import _expand_pending_archives

    supp = tmp_path / "supplements"
    supp.mkdir()
    _write_zip(supp / "bundle.zip", {"t.txt": b"CONTENT\n"})
    (supp / "bundle").mkdir()
    (supp / "bundle" / "t.txt").write_text("CONTENT\n", encoding="utf-8")

    assert _expand_pending_archives(supp) == []


@pytest.mark.parametrize("gene", ["GENEA", "GENEB"])
def test_expansion_is_gene_agnostic(tmp_path, gene):
    """Nothing in the archive path may key on a gene or publisher layout."""
    inner = _zip_bytes({"t.txt": f"Gene\tVariant\n{gene}\tc.200G>A\n".encode()})
    outer = _write_zip(tmp_path / f"{gene}.zip", {"suppl.zip": inner})

    _paths, markdown = FormatConverter().extract_zip_supplement(
        outer, dest_dir=tmp_path / f"out_{gene}"
    )

    assert gene in markdown and "c.200G>A" in markdown


def test_refold_never_trades_real_text_for_a_converter_placeholder(tmp_path):
    """A converter that regresses must not delete source it can no longer read.

    A supplement that converted cleanly in an earlier environment can fail to
    convert in a later one. Placeholder labels are treated as droppable so a
    stale placeholder cannot block a newly recovered file, which means a
    re-fold would otherwise silently replace good text with an error marker.
    """
    from harvesting.supplement_fold import (
        FOLD_BEGIN,
        FOLD_END,
        fold_supplements_into_full_context,
    )

    pmid = "00000002"  # synthetic identifier
    harvest = tmp_path / "harvest"
    supp = harvest / f"{pmid}_supplements"
    supp.mkdir(parents=True)
    # A PDF this build cannot extract text from.
    (supp / "table_s1.pdf").write_bytes(b"%PDF-1.4\n" + b"\x00" * 4096)

    real_text = "IMPORTANT_TABLE_ROW c.100A>G 7 carriers\n" * 40
    (harvest / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# Body\n\nMain article text.\n\n"
        f"{FOLD_BEGIN}\n\n# FOLDED SUPPLEMENTS (re-extraction aid)\n"
        f"\n\n# SUPPLEMENTAL FILE 1: table_s1.pdf\n\n{real_text}\n\n{FOLD_END}\n",
        encoding="utf-8",
    )

    fold_supplements_into_full_context(pmid, harvest)

    surviving = (harvest / f"{pmid}_FULL_CONTEXT.md").read_text(encoding="utf-8")
    assert "IMPORTANT_TABLE_ROW" in surviving, (
        "existing converted text was replaced by a converter placeholder"
    )
    assert "text extraction failed" not in surviving


def test_placeholder_veto_is_per_label_and_does_not_block_new_content(tmp_path):
    """A veto for one file must not block every other change to the paper.

    An all-or-nothing refusal would let a single supplement this build can no
    longer convert suppress newly expanded nested-archive members, which is the
    one measured recall gain in this area.
    """
    from harvesting.supplement_fold import (
        FOLD_BEGIN,
        FOLD_END,
        fold_supplements_into_full_context,
    )

    pmid = "00000003"  # synthetic identifier
    harvest = tmp_path / "harvest"
    supp = harvest / f"{pmid}_supplements"
    supp.mkdir(parents=True)
    (supp / "old_table.pdf").write_bytes(b"%PDF-1.4\n" + b"\x00" * 4096)
    _write_zip(supp / "never_unpacked.zip", {"new_table.txt": b"NEW_INNER_ROW\n"})

    real_text = "PRESERVED_TABLE_ROW c.100A>G 7 carriers\n" * 40
    (harvest / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# Body\n\nMain article text.\n\n"
        f"{FOLD_BEGIN}\n\n# FOLDED SUPPLEMENTS (re-extraction aid)\n"
        f"\n\n# SUPPLEMENTAL FILE 1: old_table.pdf\n\n{real_text}\n\n{FOLD_END}\n",
        encoding="utf-8",
    )

    result = fold_supplements_into_full_context(pmid, harvest)

    assert result is not None
    text = result.read_text(encoding="utf-8")
    assert "PRESERVED_TABLE_ROW" in text, "existing converted text was dropped"
    assert "NEW_INNER_ROW" in text, "nested-archive content was blocked by the veto"


def test_empty_extraction_directory_does_not_mark_an_archive_expanded(tmp_path):
    """A failed extract leaves an empty directory; the archive is still unread."""
    from harvesting.supplement_fold import _expand_pending_archives

    supp = tmp_path / "supplements"
    supp.mkdir()
    _write_zip(supp / "bundle.zip", {"t.txt": b"CONTENT\n"})
    (supp / "bundle").mkdir()  # empty: a failed prior extraction

    assert _expand_pending_archives(supp), "empty dir wrongly counted as expanded"
