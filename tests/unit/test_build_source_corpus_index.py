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
