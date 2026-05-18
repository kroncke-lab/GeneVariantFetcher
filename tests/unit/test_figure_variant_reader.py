"""Unit tests for harvesting.figure_variant_reader (parser + helpers)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from harvesting.figure_variant_reader import (
    FigureReadResult,
    PMIDFigureReport,
    _parse_response,
    find_pmid_figures,
    is_image_path,
)


def test_is_image_path():
    assert is_image_path("fig1.jpg")
    assert is_image_path(Path("a/b/T1.PNG"))
    assert not is_image_path("paper.pdf")
    assert not is_image_path("note.md")


def test_find_pmid_figures(tmp_path: Path):
    pmc_dir = tmp_path / "pmc"
    pmc_dir.mkdir()
    figs = pmc_dir / "12345_figures"
    figs.mkdir()
    (figs / "fig_1.jpg").write_bytes(b"\x00")
    (figs / "fig_2.png").write_bytes(b"\x00")
    (figs / "captions.json").write_text("{}")  # not an image

    paths = find_pmid_figures(pmc_dir, "12345")
    names = {p.name for p in paths}
    assert names == {"fig_1.jpg", "fig_2.png"}

    # Missing PMID -> empty list, not error
    assert find_pmid_figures(pmc_dir, "99999") == []


def test_parse_response_object():
    body = json.dumps({"variants": [{"protein": "R176W"}, {"protein": "L552S"}]})
    out = _parse_response(body)
    assert [v["protein"] for v in out] == ["R176W", "L552S"]


def test_parse_response_list():
    body = json.dumps([{"protein": "G601S"}])
    out = _parse_response(body)
    assert out == [{"protein": "G601S"}]


def test_parse_response_fenced_json():
    body = '```json\n{"variants": [{"protein": "K897T"}]}\n```'
    out = _parse_response(body)
    assert out == [{"protein": "K897T"}]


def test_parse_response_garbage_returns_empty():
    assert _parse_response("no json here at all") == []
    assert _parse_response("") == []


def test_pmid_report_dedupes():
    report = PMIDFigureReport(
        pmid="1",
        gene="KCNH2",
        per_figure=[
            FigureReadResult(
                image_path="a.jpg",
                variants=[
                    {"protein": "R176W"},
                    {"protein": "L552S"},
                ],
            ),
            FigureReadResult(
                image_path="b.jpg",
                variants=[
                    {"protein": "R176W"},  # duplicate
                    {"protein": "G628S"},
                ],
            ),
        ],
    )
    proteins = sorted(v["protein"] for v in report.distinct_variants)
    assert proteins == ["G628S", "L552S", "R176W"]


def test_pmid_report_dedupe_skips_empty():
    report = PMIDFigureReport(
        pmid="1",
        gene="KCNH2",
        per_figure=[
            FigureReadResult(
                image_path="a.jpg",
                variants=[
                    {"protein": "", "cdna": ""},
                    {"protein": "K897T"},
                ],
            )
        ],
    )
    assert [v["protein"] for v in report.distinct_variants] == ["K897T"]
