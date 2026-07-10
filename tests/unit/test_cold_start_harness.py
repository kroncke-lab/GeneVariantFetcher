"""Cold-start harness: gene gating and hermetic command construction."""

from pathlib import Path

from benchmarks.cold_start_eval.run_cold_start import (
    build_cold_start_command,
    cold_start_env,
    hermetic_env,
    is_cold_start_gene,
)


def test_known_genes_are_rejected():
    for gene in ("KCNH2", "BRCA1", "MYBPC3", "kcnq1"):
        ok, _ = is_cold_start_gene(gene)
        assert ok is False


def test_unseen_gene_is_accepted():
    ok, _ = is_cold_start_gene("LDLR")
    assert ok is True


def test_command_omits_pmid_file_and_isolates_corpus(tmp_path):
    cmd = build_cold_start_command(
        gene="LDLR",
        email="x@y.z",
        output_dir=tmp_path / "out",
        disease="familial hypercholesterolemia",
    )
    # No gold leak: discovery must run, so --pmid-file must be absent.
    assert "--pmid-file" not in cmd
    # No corpus write-back.
    assert "--no-corpus-sync" in cmd
    assert "gvf-run" in cmd and "LDLR" in cmd
    assert "--disease" in cmd
    # Corpus isolation goes through the environment, not a CLI flag.
    env = cold_start_env(tmp_path / "corpus")
    assert env["GVF_CORPUS_DIR"] == str(tmp_path / "corpus")


def test_source_recovery_toggle():
    cmd = build_cold_start_command(
        gene="LDLR",
        email="x@y.z",
        output_dir=Path("out"),
        source_recovery=False,
    )
    assert "--no-source-recovery" in cmd


def test_hermetic_env_strips_warm_start(tmp_path):
    """The live-run env drops warm-start GVF_* vars and pins the fresh corpus,
    while keeping credentials so the run can still authenticate."""
    base = {
        "PATH": "/usr/bin",
        "ANTHROPIC_API_KEY": "secret",
        "GVF_RESUME_DIR": "/old/run",
        "GVF_EXTRACTION_TOP_N": "50",
        "GVF_SOMETHING_DIR": "/cached",
        "GVF_CORPUS_DIR": "/stale/corpus",
    }
    env = hermetic_env(base, tmp_path / "fresh")
    assert env["GVF_CORPUS_DIR"] == str(tmp_path / "fresh")
    assert env["ANTHROPIC_API_KEY"] == "secret"
    assert env["PATH"] == "/usr/bin"
    assert "GVF_RESUME_DIR" not in env
    assert "GVF_EXTRACTION_TOP_N" not in env
    assert "GVF_SOMETHING_DIR" not in env  # GVF_*_DIR heuristic
