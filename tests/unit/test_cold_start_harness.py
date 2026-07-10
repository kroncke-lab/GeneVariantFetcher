"""Cold-start harness: gene gating and hermetic command construction."""

from pathlib import Path

from benchmarks.cold_start_eval.run_cold_start import (
    build_cold_start_command,
    cold_start_env,
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
        corpus_dir=tmp_path / "corpus",
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
        corpus_dir=Path("corpus"),
        source_recovery=False,
    )
    assert "--no-source-recovery" in cmd
