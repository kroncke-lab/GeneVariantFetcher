"""Precision / guard regression cases for the curated extraction benchmark.

Complements the recall papers: those check we *find* gold variants; these check
we do NOT mint carrier evidence from annotation/predictor/population data or from
wrong-gene rows. Source data is real (PMID 33013630 Table 1 + variantFeatures
reference residues). Self-contained: no network, no LLM, no sibling DB needed.

Run:
    .venv/bin/python -m pytest benchmarks/curated_extraction_eval/negative_cases -q
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pipeline.table_router import MarkdownTable, parse_routed_table  # noqa: E402

_VARIANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z*]|[a-z]{2,3})$")


def _load_tsv(name: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    for line in (HERE / name).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        cells = line.split("\t")
        if header is None:
            header = cells
            continue
        rows.append(dict(zip(header, cells)))
    return rows


def _load_markdown_table(name: str) -> MarkdownTable:
    pipe_lines = [
        ln
        for ln in (HERE / name).read_text().splitlines()
        if ln.lstrip().startswith("|")
    ]
    header_line = pipe_lines[0]
    # drop the |---|---| separator row
    data_lines = [ln for ln in pipe_lines[2:]]
    header_cells = [c.strip() for c in header_line.strip().strip("|").split("|")]
    return MarkdownTable(
        table_id="T1",
        caption="Table 1. Nonsynonymous variants identified in SUDEP",
        header_line=header_line,
        header_cells=header_cells,
        data_lines=data_lines,
        char_start=0,
        char_end=0,
    )


def test_annotation_table_emits_no_carriers():
    """The PMID 33013630 gnomAD/SIFT/PolyPhen table -> 0 carriers under a KCNH2 job.

    Exercises the real extraction guard end-to-end on the deterministic table
    path: blank-cell rows inherit their off-target gene and are dropped, and the
    population/prediction-only columns suppress the infer-one-carrier-per-row
    fallback. If the protocol regresses, this returns variants instead of [].
    """
    table = _load_markdown_table("annotation_table_pmid33013630.md")
    variants = parse_routed_table(
        table, {"gene": 0, "protein": 1, "patient_count": 2}, "KCNH2"
    )
    assert variants == [], "annotation table minted carriers: " + ", ".join(
        v.get("protein_notation") or "?" for v in variants
    )


def test_variantfeatures_residues_flag_exactly_the_wrong_gene_rows():
    """The committed variantFeatures aa_ref slice catches every wrong-gene row
    and clears every genuine KCNH2 row — independent of any PMID narrative."""
    ref = {
        (r["gene"], int(r["aa_pos"])): r["aa_ref"]
        for r in _load_tsv("variantfeatures_residue_reference.tsv")
    }
    cases = _load_tsv("cases.tsv")
    assert cases, "no cases loaded"
    for c in cases:
        m = _VARIANT_RE.match(c["variant"])
        assert m, f"unparseable variant {c['variant']}"
        claimed_ref, pos = m.group(1), int(m.group(2))
        aa_ref = ref.get(("KCNH2", pos))
        assert aa_ref, f"no KCNH2 reference residue for position {pos}"
        residue_matches = claimed_ref == aa_ref
        if c["evidence_class"] == "wrong_gene":
            assert not residue_matches, (
                f"{c['variant']} should be wrong-gene: claims {claimed_ref}@{pos}, "
                f"KCNH2 ref is {aa_ref}"
            )
        else:
            assert residue_matches, (
                f"{c['variant']} should be a genuine KCNH2 residue: claims "
                f"{claimed_ref}@{pos}, KCNH2 ref is {aa_ref}"
            )


def test_only_text_supported_rows_expect_a_carrier():
    """Guard against silent re-introduction of the over-attribution bug: of all
    33013630 KCNH2 rows, only the two with explicit SUDEP-patient text support
    may carry a nonzero expected count."""
    cases = _load_tsv("cases.tsv")
    nonzero = {c["variant"] for c in cases if int(c["expected_carriers"]) > 0}
    assert nonzero == {"R744*", "G924A"}, nonzero
