#!/usr/bin/env python3
"""Clean-room cold-start benchmark: any gene, empty corpus, no gold leak.

The curated extraction benchmark measures EXTRACTION quality on a frozen,
pre-fetched, PMID-fixed corpus -- it deliberately skips discovery, filtering,
and acquisition. This harness measures the opposite: whether ``gvf-run`` can
cold-start an *unseen* gene end to end (discover -> acquire -> extract) with:

  * an ISOLATED, empty corpus (``GVF_CORPUS_DIR`` -> a fresh dir) so no cached
    source warm-starts the run,
  * ``--no-corpus-sync`` so the cold gene is not folded back into the repo corpus,
  * a fresh ``--output`` dir so no prior run is scavenged,
  * NO ``--pmid-file`` so PMIDs come from live discovery, never a gold list.

Because an unseen gene has no gold CSV, recovery runs DB-observed and per-layer
recall scoring is skipped -- discovery/acquisition metrics come from the
gold-free source-QC bundle. To *score* extraction, curate a gold answer
separately (e.g. a ``<GENE>_recall_input.csv``) and run ``scripts/run_recall_suite.py``
against the produced DB AFTER the run -- never before (that would leak gold).

Default is ``--dry-run``: it prints the exact hermetic command and creates the
isolated corpus dir, but does NOT launch the (multi-hour/day) live run. Pass
``--run`` to execute.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Optional

BASE_DIR = Path(__file__).resolve().parents[2]

# Genes that already have a gold standard, are in the curated registry, are
# cached in corpus/, or are cardiac -- i.e. NOT a genuine cold start. Sourced
# from gene_variant_fetcher_gold_standard/, benchmarks/curated_extraction_eval/
# registry.tsv, corpus/, and config/cardiac_gene_synonyms.json.
NOT_COLD_START = {
    "KCNH2",
    "KCNQ1",
    "RYR2",
    "SCN5A",
    "APOE",
    "BRCA1",
    "BRCA2",
    "MYBPC3",
    "KCNE1",
    "KCNE2",
    "CACNA1C",
    "CALM1",
    "CALM2",
    "KCNJ2",
    "GJA5",
    "KCNA5",
    "NPPA",
    "PITX2",
}

# Suggested genuine cold-start candidates (unseen, non-cardiac, uncached).
SUGGESTED = {
    "LDLR": "familial hypercholesterolemia",
    "MLH1": "Lynch syndrome",
    "MSH2": "Lynch syndrome",
    "TP53": "Li-Fraumeni syndrome",
    "CFTR": "cystic fibrosis",
    "PAH": "phenylketonuria",
    "PTEN": "PTEN hamartoma tumor syndrome",
    "VHL": "von Hippel-Lindau disease",
}


def _covered_genes() -> set[str]:
    """Genes that are NOT a genuine cold start, derived at runtime so the guard
    can't drift from reality. Unions the hardcoded floor with the gold recall
    CSVs, curated gold overrides + registry, cardiac synonyms, and top-level
    corpus/ gene dirs. Best-effort per source (a missing/unreadable source is
    skipped, never fatal)."""
    covered = set(NOT_COLD_START)

    for d in (
        BASE_DIR / "gene_variant_fetcher_gold_standard" / "normalized",
        BASE_DIR / "benchmarks" / "curated_extraction_eval" / "gold_overrides",
    ):
        try:
            for csv_path in d.glob("*_recall_input.csv"):
                covered.add(csv_path.name[: -len("_recall_input.csv")].upper())
        except OSError:
            pass

    registry = BASE_DIR / "benchmarks" / "curated_extraction_eval" / "registry.tsv"
    try:
        lines = registry.read_text().splitlines()
        header = lines[0].split("\t") if lines else []
        gcol = header.index("gene") if "gene" in header else 0
        for line in lines[1:]:
            cells = line.split("\t")
            if len(cells) > gcol and cells[gcol].strip():
                covered.add(cells[gcol].strip().upper())
    except (OSError, ValueError, IndexError):
        pass

    try:
        data = json.loads(
            (BASE_DIR / "config" / "cardiac_gene_synonyms.json").read_text()
        )
        covered.update(str(k).upper() for k in data)
    except (OSError, json.JSONDecodeError):
        pass

    try:
        for child in (BASE_DIR / "corpus").iterdir():
            if child.is_dir() and not child.name.startswith("."):
                covered.add(child.name.upper())
    except OSError:
        pass

    return covered


def is_cold_start_gene(gene: str) -> tuple[bool, str]:
    """Return (is_genuine_cold_start, reason)."""
    g = gene.strip().upper()
    if g in _covered_genes():
        return (
            False,
            f"{g} already has gold/registry/corpus/cardiac coverage — not a "
            f"clean room. Pick an unseen gene, e.g. {', '.join(sorted(SUGGESTED))}.",
        )
    return (True, f"{g} looks like a genuine cold start.")


def build_cold_start_command(
    *,
    gene: str,
    email: str,
    output_dir: Path,
    disease: Optional[str] = None,
    source_recovery: bool = True,
) -> list[str]:
    """The hermetic ``gvf-run`` invocation for a cold-start measurement."""
    cmd = [
        sys.executable,
        "-m",
        "cli",
        "gvf-run",
        gene,
        "--email",
        email,
        "--output",
        str(output_dir),
        "--no-corpus-sync",
    ]
    if not source_recovery:
        cmd.append("--no-source-recovery")
    if disease:
        cmd.extend(["--disease", disease])
    return cmd


# GVF_* env vars that would warm-start a cold run from prior state / cached work
# and so defeat the clean room. Stripped from the subprocess env before launch;
# credentials (non-GVF) and the harness's own GVF_CORPUS_DIR are preserved.
WARM_START_ENV = (
    "GVF_RESUME_DIR",
    "GVF_REUSE_FULL_CONTEXT_BYTES",
    "GVF_EXTRACTION_TOP_N",
    "GVF_EXTRACTION_PRIORITY_OFFSET",
    "GVF_QA_MODEL",
    "GVF_MAX_WORKERS",
)


def cold_start_env(corpus_dir: Path) -> dict[str, str]:
    """The single explicit override: pin the corpus dir to the fresh isolated dir."""
    return {"GVF_CORPUS_DIR": str(corpus_dir)}


def hermetic_env(base: dict[str, str], corpus_dir: Path) -> dict[str, str]:
    """Build the cold-start subprocess env: strip warm-start/tuning GVF_* vars
    (and any GVF_*_DIR pointing at a cached tree) that would reuse prior state,
    then pin GVF_CORPUS_DIR to the fresh dir. Non-GVF env — PATH, API
    credentials — is preserved so the run can still authenticate and discover.
    """
    env = dict(base)
    for key in WARM_START_ENV:
        env.pop(key, None)
    for key in list(env):
        if key.startswith("GVF_") and key.endswith("_DIR") and key != "GVF_CORPUS_DIR":
            env.pop(key, None)
    env.update(cold_start_env(corpus_dir))
    return env


def summarize_acquisition(run_output_root: Path, gene: str) -> Optional[dict]:
    """Read the gold-free source-QC summary produced by the run, if present."""
    candidates = sorted(
        run_output_root.glob(f"{gene}/*/source_qc/source_acquisition_summary.json")
    )
    if not candidates:
        return None
    try:
        return json.loads(candidates[-1].read_text())
    except (OSError, json.JSONDecodeError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene", help="An UNSEEN, non-cardiac gene (e.g. LDLR).")
    parser.add_argument("--email", default="brett.kroncke@gmail.com")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Run output root (default: benchmarks/cold_start_eval/runs/<GENE>).",
    )
    parser.add_argument(
        "--corpus-dir",
        type=Path,
        default=None,
        help="Isolated corpus dir (default: <output>/.coldstart_corpus).",
    )
    parser.add_argument(
        "--disease",
        default=None,
        help="Optional phenotype clause (defaults to a suggested one if known).",
    )
    parser.add_argument("--no-source-recovery", action="store_true")
    parser.add_argument(
        "--run",
        action="store_true",
        help="Actually execute the (long) live run. Default is a dry run.",
    )
    parser.add_argument(
        "--allow-non-cold-gene",
        action="store_true",
        help="Proceed even if the gene is not a genuine cold start.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Wipe and recreate a non-empty isolated corpus dir instead of refusing.",
    )
    args = parser.parse_args()

    gene = args.gene.strip().upper()
    cold, reason = is_cold_start_gene(gene)
    print(reason)
    if not cold and not args.allow_non_cold_gene:
        print("Refusing non-cold-start gene; pass --allow-non-cold-gene to override.")
        return 2

    # Resolve to absolute so the harness (cwd) and the subprocess (cwd=BASE_DIR)
    # agree on which dirs are created/isolated vs. read.
    output_dir = (
        args.output or (BASE_DIR / "benchmarks" / "cold_start_eval" / "runs")
    ).resolve()
    corpus_dir = (
        args.corpus_dir or (output_dir / f".coldstart_corpus_{gene}")
    ).resolve()
    disease = args.disease or SUGGESTED.get(gene)

    # A non-empty isolated corpus is NOT a clean room — a prior cold run would
    # warm-start this one. Refuse (or --force wipe) instead of silently reusing.
    if corpus_dir.exists() and any(corpus_dir.iterdir()):
        if not args.force:
            print(
                f"Isolated corpus dir is not empty: {corpus_dir}\n"
                "It would warm-start the run. Remove it, pick a fresh --corpus-dir, "
                "or pass --force to wipe it."
            )
            return 2
        import shutil

        shutil.rmtree(corpus_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    corpus_dir.mkdir(parents=True, exist_ok=True)

    cmd = build_cold_start_command(
        gene=gene,
        email=args.email,
        output_dir=output_dir,
        disease=disease,
        source_recovery=not args.no_source_recovery,
    )
    env_overrides = cold_start_env(corpus_dir)

    env_prefix = " ".join(f"{k}={shlex.quote(v)}" for k, v in env_overrides.items())
    printable = f"{env_prefix} {shlex.join(cmd)}"
    print("\nHermetic cold-start command:")
    print(f"  {printable}\n")
    print("  (the live run also strips warm-start GVF_* env: resume/reuse/tuning)")
    print(f"Isolated corpus dir (empty): {corpus_dir}")
    print(f"Output root: {output_dir}")

    if not args.run:
        print("\n[dry-run] Not executing. Re-run with --run to launch the live run.")
        print("After a real run, acquisition metrics are read from:")
        print(
            f"  {output_dir}/{gene}/<timestamp>/source_qc/source_acquisition_summary.json"
        )
        return 0

    import os

    run_env = hermetic_env(dict(os.environ), corpus_dir)
    print("\n[run] Launching live cold-start run (this can take a long time)…")
    result = subprocess.run(cmd, env=run_env, cwd=str(BASE_DIR))
    acquisition = summarize_acquisition(output_dir, gene)
    if acquisition is not None:
        print("\nAcquisition summary (gold-free):")
        print(json.dumps(acquisition.get("pmid_coverage", acquisition), indent=2))
    else:
        print("\nNo source-QC summary found; inspect the run dir manually.")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
