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

BODY_ONLY_MARKER = "<!-- GVF_SUPPLEMENT_SURFACE_STATUS: unavailable -->\n\n"


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


def test_body_only_source_is_indexed_incomplete_and_complete_copy_wins(tmp_path):
    roots = tmp_path / "runs" / "BMPR2"
    body_only = roots / "old" / "pmc_fulltext"
    complete = roots / "new" / "pmc_fulltext"
    body_only.mkdir(parents=True)
    complete.mkdir(parents=True)
    (body_only / "12345678_FULL_CONTEXT.md").write_text(
        BODY_ONLY_MARKER + USABLE_MD + ("larger body. " * 300),
        encoding="utf-8",
    )
    (complete / "12345678_FULL_CONTEXT.md").write_text(USABLE_MD, encoding="utf-8")
    out = tmp_path / "external_corpus"

    result = run_builder(
        "--apply", "--roots", str(roots), "--out", str(out), "--assume-gene", "BMPR2"
    )

    assert result.returncode == 0, result.stderr
    (entry,) = read_index(out)
    assert entry["full_text_status"] == "ok"
    stored = out / "BMPR2" / "12345678" / "12345678_FULL_CONTEXT.md"
    assert stored.read_text(encoding="utf-8") == USABLE_MD


def test_body_only_source_status_is_not_ok(tmp_path):
    source = tmp_path / "runs" / "BMPR2" / "12345679"
    source.mkdir(parents=True)
    (source / "12345679_FULL_CONTEXT.md").write_text(
        BODY_ONLY_MARKER + USABLE_MD,
        encoding="utf-8",
    )
    out = tmp_path / "external_corpus"

    result = run_builder(
        "--apply", "--roots", str(source.parent.parent), "--out", str(out)
    )

    assert result.returncode == 0, result.stderr
    (entry,) = read_index(out)
    assert entry["full_text_status"] == "body_only"


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


def test_scoped_noop_preserves_unrelated_index_manifest_and_does_not_bless_mutation(
    tmp_path,
):
    roots = tmp_path / "runs"
    for gene, pmid in (("BRCA1", "11111111"), ("BMPR2", "22222222")):
        source = roots / gene / pmid
        source.mkdir(parents=True)
        (source / f"{pmid}_FULL_CONTEXT.md").write_text(
            USABLE_MD.replace("BMPR2", gene),
            encoding="utf-8",
        )
    out = tmp_path / "external_corpus"
    seeded = run_builder("--apply", "--roots", str(roots), "--out", str(out))
    assert seeded.returncode == 0, seeded.stderr

    index_path = out / "INDEX.json"
    manifest_path = out / "MANIFEST.sha256"
    before_mtimes = (index_path.stat().st_mtime_ns, manifest_path.stat().st_mtime_ns)

    # Corrupt an unrelated corpus payload, then perform a BMPR2-only no-op sync.
    # The scoped operation must neither inspect/re-hash nor bless BRCA1.
    unrelated = out / "BRCA1" / "11111111" / "11111111_FULL_CONTEXT.md"
    unrelated.write_text(unrelated.read_text() + "\nunreviewed mutation\n")
    scoped = run_builder(
        "--apply",
        "--roots",
        str(roots / "BMPR2"),
        "--out",
        str(out),
        "--assume-gene",
        "BMPR2",
        "--gene",
        "BMPR2",
    )
    assert scoped.returncode == 0, scoped.stderr
    assert "upgrade=0 noop=1" in scoped.stdout
    assert before_mtimes == (
        index_path.stat().st_mtime_ns,
        manifest_path.stat().st_mtime_ns,
    )

    data = json.loads(index_path.read_text())
    assert set(data["genes"]) == {"BMPR2", "BRCA1"}
    verified = run_builder("--verify", "--out", str(out))
    assert verified.returncode == 1
    assert "VERIFY FAILED" in verified.stdout


def test_scoped_upgrade_merges_index_and_incrementally_updates_manifest(tmp_path):
    roots = tmp_path / "runs"
    for gene, pmid in (("BRCA1", "33333333"), ("BMPR2", "44444444")):
        source = roots / gene / pmid
        source.mkdir(parents=True)
        (source / f"{pmid}_FULL_CONTEXT.md").write_text(
            USABLE_MD.replace("BMPR2", gene),
            encoding="utf-8",
        )
    out = tmp_path / "external_corpus"
    seeded = run_builder("--apply", "--roots", str(roots), "--out", str(out))
    assert seeded.returncode == 0, seeded.stderr

    bmpr2_source = roots / "BMPR2" / "44444444" / "44444444_FULL_CONTEXT.md"
    bmpr2_source.write_text(
        bmpr2_source.read_text() + ("\nAdditional cohort table detail. " * 100)
    )
    upgraded = run_builder(
        "--apply",
        "--roots",
        str(roots / "BMPR2"),
        "--out",
        str(out),
        "--assume-gene",
        "BMPR2",
        "--gene",
        "BMPR2",
    )

    assert upgraded.returncode == 0, upgraded.stderr
    assert "upgrade=1" in upgraded.stdout
    data = json.loads((out / "INDEX.json").read_text())
    assert set(data["genes"]) == {"BMPR2", "BRCA1"}
    verified = run_builder("--verify", "--out", str(out))
    assert verified.returncode == 0, verified.stdout


def test_scoped_apply_refuses_corrupt_index_before_touching_payload(tmp_path):
    roots = tmp_path / "runs"
    source = roots / "BMPR2" / "55555555"
    source.mkdir(parents=True)
    source_file = source / "55555555_FULL_CONTEXT.md"
    source_file.write_text(USABLE_MD, encoding="utf-8")
    out = tmp_path / "external_corpus"
    seeded = run_builder("--apply", "--roots", str(roots), "--out", str(out))
    assert seeded.returncode == 0, seeded.stderr
    payload = out / "BMPR2" / "55555555" / "55555555_FULL_CONTEXT.md"
    before = payload.read_bytes()
    (out / "INDEX.json").write_text("{broken", encoding="utf-8")
    source_file.write_text(USABLE_MD + "\nnew source content\n", encoding="utf-8")

    scoped = run_builder(
        "--apply",
        "--gene",
        "BMPR2",
        "--assume-gene",
        "BMPR2",
        "--roots",
        str(roots / "BMPR2"),
        "--out",
        str(out),
    )

    assert scoped.returncode == 2
    assert "cannot safely load existing corpus index" in scoped.stderr
    assert payload.read_bytes() == before


def test_scoped_apply_refuses_incomplete_index_row_before_touching_payload(tmp_path):
    roots = tmp_path / "runs"
    source = roots / "BMPR2" / "55555556"
    source.mkdir(parents=True)
    source_file = source / "55555556_FULL_CONTEXT.md"
    source_file.write_text(USABLE_MD, encoding="utf-8")
    out = tmp_path / "external_corpus"
    seeded = run_builder("--apply", "--roots", str(roots), "--out", str(out))
    assert seeded.returncode == 0, seeded.stderr
    payload = out / "BMPR2" / "55555556" / "55555556_FULL_CONTEXT.md"
    before = payload.read_bytes()

    index_path = out / "INDEX.json"
    index = json.loads(index_path.read_text())
    del index["genes"]["BMPR2"]["55555556"]["chosen_from"]
    index_path.write_text(json.dumps(index), encoding="utf-8")
    source_file.write_text(USABLE_MD + "\nnew source content\n", encoding="utf-8")

    scoped = run_builder(
        "--apply",
        "--gene",
        "BMPR2",
        "--roots",
        str(roots / "BMPR2"),
        "--out",
        str(out),
    )

    assert scoped.returncode == 2
    assert "cannot safely load existing corpus index row" in scoped.stderr
    assert payload.read_bytes() == before


def test_scoped_apply_refuses_malformed_manifest_before_touching_payload(tmp_path):
    roots = tmp_path / "runs"
    source = roots / "BMPR2" / "55555557"
    source.mkdir(parents=True)
    source_file = source / "55555557_FULL_CONTEXT.md"
    source_file.write_text(USABLE_MD, encoding="utf-8")
    out = tmp_path / "external_corpus"
    seeded = run_builder("--apply", "--roots", str(roots), "--out", str(out))
    assert seeded.returncode == 0, seeded.stderr
    payload = out / "BMPR2" / "55555557" / "55555557_FULL_CONTEXT.md"
    before = payload.read_bytes()

    (out / "MANIFEST.sha256").write_text(
        "# total_files not-an-integer\n# total_bytes 10\n", encoding="utf-8"
    )
    source_file.write_text(USABLE_MD + "\nnew source content\n", encoding="utf-8")

    scoped = run_builder(
        "--apply",
        "--gene",
        "BMPR2",
        "--roots",
        str(roots / "BMPR2"),
        "--out",
        str(out),
    )

    assert scoped.returncode == 2
    assert "cannot safely load existing corpus manifest" in scoped.stderr
    assert payload.read_bytes() == before


def test_scoped_rebuild_is_rejected_without_deleting_other_gene(tmp_path):
    out = tmp_path / "external_corpus"
    unrelated = out / "BRCA1" / "66666666"
    unrelated.mkdir(parents=True)
    payload = unrelated / "66666666_FULL_CONTEXT.md"
    payload.write_text(USABLE_MD.replace("BMPR2", "BRCA1"), encoding="utf-8")

    result = run_builder(
        "--apply",
        "--rebuild",
        "--gene",
        "BMPR2",
        "--out",
        str(out),
    )

    assert result.returncode == 2
    assert "--rebuild cannot be combined with --gene" in result.stderr
    assert payload.exists()
