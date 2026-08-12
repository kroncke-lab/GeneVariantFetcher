"""Index paths must survive a corpus that lives outside the repo.

`corpus/` is an absolute symlink to an external volume on Brett's workstation,
so `out.resolve()` leaves the repo and `Path.relative_to(REPO)` raises
ValueError. Before the fix, *every* invocation of build_source_corpus.py died
on the first indexed entry (`gvf-run` Step 4.5 logged
"corpus sync (build_source_corpus.py) exited 1"), and no new source was ever
folded into the external corpus.

These tests run the script against an out dir outside the repo — the same
shape as the symlinked volume — and pin both index fields that used to crash:
`fulltext` (falls back to the stable `corpus/<gene>/<pmid>/<file>` interface
path) and `chosen_from` (falls back to the absolute source dir).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "build_source_corpus.py"

# >= MIN_EXTRACTION_INPUT_SIZE (500) and not the abstract-only fallback shape.
USABLE_MD = (
    "# Novel BMPR2 variants in pulmonary arterial hypertension\n\n"
    "## Methods\n\nWe sequenced probands and relatives.\n\n"
    "| Variant | Carriers | Affected |\n|---|---|---|\n"
    "| p.Arg491Trp | 4 | 3 |\n| c.994C>T | 2 | 1 |\n\n"
) + ("Segregation and clinical detail. " * 20)


def run_builder(*argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        timeout=120,
    )


def read_index(out: Path) -> list[dict]:
    data = json.loads((out / "INDEX.json").read_text())
    return [entry for pmids in data["genes"].values() for entry in pmids.values()]


def test_out_dir_outside_repo_indexes_with_corpus_interface_paths(tmp_path):
    src = tmp_path / "runs" / "BMPR2" / "12345678"
    src.mkdir(parents=True)
    (src / "12345678_FULL_CONTEXT.md").write_text(USABLE_MD)
    out = tmp_path / "external_corpus"

    result = run_builder(
        "--apply", "--roots", str(src.parent.parent), "--out", str(out)
    )

    assert result.returncode == 0, result.stderr
    (entry,) = read_index(out)
    assert entry["fulltext"] == "corpus/BMPR2/12345678/12345678_FULL_CONTEXT.md"
    assert entry["chosen_from"] == str(src)
    assert (out / "BMPR2" / "12345678" / "12345678_FULL_CONTEXT.md").exists()


def test_assume_gene_files_run_dir_candidates_for_a_new_gene(tmp_path):
    """A scoped run-dir root has no gene component below it, so a NEW gene
    (not in KNOWN_GENES, empty corpus) is only reachable via --assume-gene —
    gvf-run's Step 4.5 passes the run's gene."""
    run_dir = tmp_path / "20260807_163246"
    (run_dir / "pmc_fulltext").mkdir(parents=True)
    (run_dir / "pmc_fulltext" / "87654321_FULL_CONTEXT.md").write_text(USABLE_MD)
    out = tmp_path / "external_corpus"

    skipped = run_builder("--apply", "--roots", str(run_dir), "--out", str(out))
    assert skipped.returncode == 0, skipped.stderr
    assert "skipped (gene not inferable)" in skipped.stdout
    assert not (out / "INDEX.json").exists()

    filed = run_builder(
        "--apply",
        "--roots",
        str(run_dir),
        "--out",
        str(out),
        "--assume-gene",
        "bmpr2",
    )
    assert filed.returncode == 0, filed.stderr
    data = json.loads((out / "INDEX.json").read_text())
    assert list(data["genes"]) == ["BMPR2"]
    assert list(data["genes"]["BMPR2"]) == ["87654321"]
    assert (out / "BMPR2" / "87654321" / "87654321_FULL_CONTEXT.md").exists()


def test_assume_gene_beats_existing_shared_pmid_membership(tmp_path):
    """A scoped run must not be filed under an arbitrary prior PMID gene."""
    run_dir = tmp_path / "run"
    (run_dir / "pmc_fulltext").mkdir(parents=True)
    (run_dir / "pmc_fulltext" / "24433082_FULL_CONTEXT.md").write_text(USABLE_MD)
    out = tmp_path / "external_corpus"
    for gene in ("BMPR2", "BRCA1"):
        seeded = out / gene / "24433082"
        seeded.mkdir(parents=True)
        (seeded / "24433082_FULL_CONTEXT.md").write_text(USABLE_MD)

    result = run_builder(
        "--apply",
        "--roots",
        str(run_dir),
        "--out",
        str(out),
        "--assume-gene",
        "KCNH2",
    )

    assert result.returncode == 0, result.stderr
    assert (out / "KCNH2" / "24433082" / "24433082_FULL_CONTEXT.md").exists()


def test_existing_corpus_only_entry_indexes_without_crashing(tmp_path):
    """The production crash: a corpus-resident paper with no source candidate."""
    out = tmp_path / "external_corpus"
    seeded = out / "APOE" / "10025785"
    seeded.mkdir(parents=True)
    (seeded / "10025785_FULL_CONTEXT.md").write_text(USABLE_MD)
    empty_root = tmp_path / "no_runs"
    empty_root.mkdir()

    result = run_builder("--apply", "--roots", str(empty_root), "--out", str(out))

    assert result.returncode == 0, result.stderr
    (entry,) = read_index(out)
    assert entry["fulltext"] == "corpus/APOE/10025785/10025785_FULL_CONTEXT.md"
    assert entry["chosen_from"] == "corpus"
