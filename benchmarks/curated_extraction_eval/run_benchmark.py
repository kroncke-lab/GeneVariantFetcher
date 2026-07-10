#!/usr/bin/env python3
"""Run the curated-extraction-eval benchmark: recall + MAE on fixed papers.

WHAT THIS IS
  A cheap, FIXED regression harness. It scores extraction against the frozen
  gold subset in ./gold (hand-picked, strategy-diverse gold papers across the
  registry genes) and reports recall + rows-mode MAE +
  precision, plus a per-paper breakdown. Use it to tell whether a change to the
  prompts / harness / guardrails / matcher helped or regressed, WITHOUT paying to
  re-run the whole gold standard. See README.md for the full story.

TWO MODES
  score   (default, NO LLM tokens): score existing DB(s) on the subset. By
          default it scores the four canonical DBs from docs/RECALL_STATUS.md, so
          `python run_benchmark.py` works with zero args. Point it at your own DB
          with `--db GENE=/path/to.db` (repeatable) to check a matcher/scorer/
          normalizer change, or a DB you produced any other way.

  extract (opt-in, COSTS LLM TOKENS): process the set through the FULL, REGULAR,
          DEFAULT GVF pipeline and score the result. For each gene in the set it
          runs `python -m cli gvf-run <GENE> --pmid-file pmids/<GENE>.txt` with
          default settings — i.e. the exact scripts a normal run uses after selection:
          it downloads (or reuses cached) source, preprocesses/scouts, extracts,
          migrates, runs source recovery and the ClinVar/PubTator/figure recovery
          layers, applies the standard guardrails, builds per-gene DBs, then
          scores them. This is the representative path for prompt/harness/
          guardrail changes. Needs --email and a working .env (one LLM key +
          NCBI_EMAIL). Cost is bounded because the set is small. `--fast` adds
          --no-source-recovery for a quick, NON-representative sanity pass that
          skips source-acquisition refresh; do not use it for the benchmark
          number.

Mechanics: this is a thin wrapper around the production scorer
(scripts/run_recall_suite.py) pointed at ./gold, so scoring is byte-identical to
production — just restricted to these PMIDs because the gold subset only contains
them. Nothing about recall/MAE is reimplemented here.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
GOLD_DIR = HERE / "gold"
MANIFEST = HERE / "manifest.csv"
BASELINE = HERE / "expected_baseline.json"
SCORER = REPO / "scripts" / "run_recall_suite.py"
LOCAL_SOURCE_CORPUS = HERE / "sources"

# Canonical DBs — kept in sync with docs/RECALL_STATUS.md "Current Canonical
# Baseline" and build_fixture.py.
CANONICAL_DBS = {
    "KCNH2": REPO / "results/KCNH2/e2e_working_20260529_full/02_strict/KCNH2.db",
    "KCNQ1": REPO
    / "validation_runs/20260517_203904/results/KCNQ1/20260517_204424/KCNQ1.db",
    "SCN5A": REPO
    / "validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938/SCN5A.db",
    "RYR2": REPO
    / "validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938/RYR2.db",
    # Non-cardiac diversity genes (gold via gold_overrides/). Existence-gated in
    # resolve_dbs, so checkouts without these run dirs just skip them.
    "BRCA1": REPO / "results/BRCA1/20260616_132646/BRCA1.db",
    "BRCA2": REPO
    / "results/BRCA2/20260606_134517_hereditary_breast_cancer_500/BRCA2/20260606_134519/BRCA2.refresh_20260606_205358.db",
    "MYBPC3": REPO / "results/MYBPC3/20260616_132646/MYBPC3.refresh_20260617_091043.db",
    "APOE": REPO / "results/APOE/20260621_072155_full_redo/APOE.db",
}


def python_exe() -> str:
    venv = REPO / ".venv" / "bin" / "python"
    return str(venv) if venv.exists() else sys.executable


def load_manifest() -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with MANIFEST.open(newline="") as f:
        for r in csv.DictReader(f):
            rows[f"{r['gene']}:{r['pmid']}"] = r
    return rows


def genes_from_manifest(manifest: dict[str, dict[str, str]]) -> list[str]:
    """Return every registry gene represented in the generated manifest."""
    return sorted({r["gene"].strip().upper() for r in manifest.values()})


def complete_local_source_corpus(manifest: dict[str, dict[str, str]]) -> bool:
    """Return true when sources/ can serve as GVF_CORPUS_DIR for this fixture."""
    if not LOCAL_SOURCE_CORPUS.is_dir():
        return False
    for row in manifest.values():
        gene = row["gene"].strip().upper()
        pmid = row["pmid"].strip()
        if not (
            LOCAL_SOURCE_CORPUS / gene / pmid / f"{pmid}_FULL_CONTEXT.md"
        ).is_file():
            return False
    return True


def run_scorer(db_overrides: dict[str, Path], outdir: Path) -> dict:
    cmd = [
        python_exe(),
        str(SCORER),
        "--score",
        "--gold-dir",
        str(GOLD_DIR),
        "--genes",
        ",".join(sorted(db_overrides)),
        "--results-dir",
        str(REPO / "results"),
        "--outdir",
        str(outdir),
    ]
    for gene, path in db_overrides.items():
        cmd += ["--db", f"{gene}={path}"]
    res = subprocess.run(cmd, cwd=str(REPO), capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stdout + "\n" + res.stderr + "\n")
        raise SystemExit(f"scorer failed (exit {res.returncode})")
    summary = json.loads((outdir / "summary.json").read_text())
    return summary


def _is_true(v: str | None) -> bool:
    return (v or "").strip().lower() == "true"


def _fabs(v: str | None) -> float:
    try:
        return abs(float(v or 0))
    except (ValueError, TypeError):
        return 0.0


def per_paper_table(outdir: Path, manifest: dict[str, dict[str, str]]) -> list[dict]:
    """Roll up per-paper recall + count error from each gene's discrepancies.csv.

    discrepancies.csv is always written by the scorer (one row per gold/DB
    variant on gold PMIDs). A gold variant is a row with missing_in_excel=False;
    it is matched when missing_in_sqlite is also False.
    """
    want: dict[str, set[str]] = {}
    strat: dict[tuple[str, str], str] = {}
    for r in manifest.values():
        g, p = r["gene"].upper(), r["pmid"]
        want.setdefault(g, set()).add(p)
        strat[(g, p)] = r["strategy"]

    agg: dict[tuple[str, str], dict] = {}
    for gene, pmids in want.items():
        disc = outdir / gene / "discrepancies.csv"
        if not disc.exists():
            continue
        with disc.open(newline="") as f:
            for d in csv.DictReader(f):
                pmid = d.get("pmid", "")
                if pmid not in pmids:
                    continue
                rec = agg.setdefault(
                    (gene, pmid),
                    {
                        "gold": 0,
                        "matched": 0,
                        "cmm": 0,
                        "car": 0.0,
                        "aff": 0.0,
                        "unaf": 0.0,
                    },
                )
                if _is_true(d.get("missing_in_excel")):  # DB-only extra row
                    continue
                rec["gold"] += 1
                if not _is_true(d.get("missing_in_sqlite")):  # matched
                    rec["matched"] += 1
                    if _is_true(d.get("count_mismatch")):
                        rec["cmm"] += 1
                    rec["car"] += _fabs(d.get("carriers_diff"))
                    rec["aff"] += _fabs(d.get("affected_diff"))
                    rec["unaf"] += _fabs(d.get("unaffected_diff"))

    out: list[dict] = []
    for (gene, pmid), rec in agg.items():
        rr = rec["matched"] / rec["gold"] if rec["gold"] else None
        out.append(
            {
                "gene": gene,
                "pmid": pmid,
                "strategy": strat.get((gene, pmid), "?"),
                "gold": rec["gold"],
                "matched": rec["matched"],
                "row_recall": rr,
                "cmm": rec["cmm"],
                "cdiff": f"{rec['car']:.0f}/{rec['aff']:.0f}/{rec['unaf']:.0f}",
            }
        )
    out.sort(key=lambda d: (d["strategy"], d["gene"], d["pmid"]))
    return out


def fmt_pct(v) -> str:
    return "n/a" if v is None else f"{v:.1%}"


def display_path(path: Path) -> str:
    """Return a readable path whether callers supplied absolute or relative args."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO))
    except ValueError:
        return str(path)


def headline(summary: dict) -> dict:
    rec = summary.get("aggregate_recall", {})
    mae = summary.get("aggregate_mae", {})

    def r(name):
        b = rec.get(name, {})
        return {
            "matched": b.get("matched"),
            "gold": b.get("gold"),
            "recall": b.get("recall"),
        }

    def m(name):
        b = mae.get(name, {})
        return b.get("mae")

    e2e = summary.get("aggregate_count_error_end_to_end", {})

    def e(name):
        return e2e.get(name, {}).get("mae")

    return {
        "recall": {
            k: r(k)
            for k in (
                "pmids",
                "variant_rows",
                "unique_variants",
                "patients",
                "affected",
                "unaffected",
            )
        },
        "mae": {k: m(k) for k in ("carriers", "affected", "unaffected")},
        "count_error_end_to_end": {
            k: e(k) for k in ("carriers", "affected", "unaffected")
        },
    }


def print_scorecard(summary: dict, papers: list[dict], baseline: dict | None) -> None:
    h = headline(summary)
    print("\n" + "=" * 72)
    print("CURATED EXTRACTION EVAL — fixed curated gold subset")
    print("=" * 72)
    print("\nAggregate recall (matched / gold on the curated subset):")
    base_rec = (baseline or {}).get("recall", {})
    for name, b in h["recall"].items():
        delta = ""
        if (
            name in base_rec
            and b["recall"] is not None
            and base_rec[name].get("recall") is not None
        ):
            d = b["recall"] - base_rec[name]["recall"]
            if abs(d) >= 0.0005:
                delta = f"  ({'+' if d >= 0 else ''}{d * 100:.1f}pp vs baseline)"
        print(
            f"  {name:16s} {str(b['matched']):>4}/{str(b['gold']):<4} {fmt_pct(b['recall']):>7}{delta}"
        )

    print("\nRows-mode MAE (matched rows; lower is better):")
    base_mae = (baseline or {}).get("mae", {})
    for name, v in h["mae"].items():
        delta = ""
        if name in base_mae and v is not None and base_mae[name] is not None:
            d = v - base_mae[name]
            if abs(d) >= 0.005:
                delta = f"  ({'+' if d >= 0 else ''}{d:.3f} vs baseline {'WORSE' if d > 0 else 'better'})"
        print(f"  {name:16s} MAE={'n/a' if v is None else f'{v:.3f}'}{delta}")

    if papers:
        print("\nPer-paper (grouped by strategy):")
        print(
            f"  {'strat':7} {'gene':6} {'pmid':10} {'gold':>4} {'match':>5} {'recall':>7} "
            f"{'cmm':>3} {'|Δcar/aff/unaf|':>15}"
        )
        last = None
        for p in papers:
            if p["strategy"] != last:
                print(f"  --- {p['strategy']} ---")
                last = p["strategy"]
            rr = p["row_recall"]
            try:
                rr = f"{float(rr):.0%}"
            except (ValueError, TypeError):
                rr = str(rr)
            print(
                f"  {'':7} {p['gene']:6} {p['pmid']:10} {str(p['gold']):>4} "
                f"{str(p['matched']):>5} {rr:>7} {str(p['cmm']):>3} {p['cdiff']:>15}"
            )
    print()


def check_regression(
    summary: dict,
    baseline: dict | None,
    recall_tol: float,
    mae_tol: float,
) -> list[str]:
    """Return regression messages (empty if none) vs the frozen baseline.

    A recall dimension regresses if it drops more than ``recall_tol`` below
    baseline; a count MAE regresses if it rises more than ``mae_tol`` above
    baseline (higher MAE is worse).
    """
    if not baseline:
        return []
    h = headline(summary)
    problems: list[str] = []

    base_rec = baseline.get("recall", {})
    for name, block in h["recall"].items():
        base = (base_rec.get(name) or {}).get("recall")
        if base is None:
            continue  # dimension not tracked by the baseline
        cur = block["recall"]
        if cur is None:
            problems.append(
                f"recall.{name}: missing in current run (baseline {base:.3f})"
            )
        elif cur < base - recall_tol:
            problems.append(
                f"recall.{name}: {cur:.3f} < baseline {base:.3f} - {recall_tol}"
            )

    # Both matched-only MAE and the end-to-end count error (misses as zero) gate
    # on a rise. Including e2e is what stops the gate from re-flattering the
    # metric this branch just added; it activates once the baseline is
    # regenerated (--write-baseline) with the e2e block.
    for key, label in (("mae", "mae"), ("count_error_end_to_end", "e2e")):
        base_block = baseline.get(key, {})
        for name, cur in h.get(key, {}).items():
            base = base_block.get(name)
            if base is None:
                continue
            if cur is None:
                problems.append(
                    f"{label}.{name}: missing in current run (baseline {base:.3f})"
                )
            elif cur > base + mae_tol:
                problems.append(
                    f"{label}.{name}: {cur:.3f} > baseline {base:.3f} + {mae_tol}"
                )
    return problems


def gate_result(
    summary: dict,
    baseline: dict | None,
    recall_tol: float,
    mae_tol: float,
) -> list[str]:
    """The --fail-on-regression decision, fail-closed. A missing baseline is a
    hard failure: a gate that compared against nothing must not report success.
    """
    if baseline is None:
        return [
            "no baseline at expected_baseline.json — refusing to pass a gate that "
            "checked nothing (write one with --write-baseline)"
        ]
    return check_regression(summary, baseline, recall_tol, mae_tol)


def resolve_dbs(overrides: list[str] | None) -> dict[str, Path]:
    dbs: dict[str, Path] = {}
    for gene, path in CANONICAL_DBS.items():
        if path.exists():
            dbs[gene] = path
        else:
            print(f"  (skip {gene}: canonical DB not found at {path})")
    for item in overrides or []:
        if "=" not in item:
            raise SystemExit(f"--db must be GENE=/path/to.db, got {item!r}")
        gene, path = item.split("=", 1)
        dbs[gene.strip().upper()] = Path(path).expanduser()
    if not dbs:
        raise SystemExit(
            "No DBs to score. Pass --db GENE=/path/to.db, or run from a checkout "
            "that has the canonical DBs (see docs/RECALL_STATUS.md)."
        )
    return dbs


def do_extract(
    email: str,
    outroot: Path,
    genes: list[str],
    fast: bool,
    source_recovery_timeout_s: int,
    use_local_source_corpus: bool,
) -> dict[str, Path]:
    """Process the set through the standard `gvf gvf-run` pipeline, return DBs.

    Default = the full, regular, default pipeline. `fast` adds
    --no-source-recovery for a quick, non-representative pass that skips
    source-acquisition refresh.
    """
    outroot.mkdir(parents=True, exist_ok=True)
    produced: dict[str, Path] = {}
    if fast:
        print(
            "\n!! EXTRACT MODE (--fast): adds --no-source-recovery — a QUICK, "
            "NON-REPRESENTATIVE pass that skips source-acquisition refresh. "
            "Drop --fast for a real measurement.\n"
        )
    else:
        print(
            "\n!! EXTRACT MODE — running the FULL DEFAULT GVF pipeline "
            "(python -m cli gvf-run, default settings) on the set: "
            "downloads/reuses source, "
            "preprocesses/scouts, extracts, source recovery, and the figure/ClinVar/"
            "PubTator recovery layers. Representative; COSTS LLM TOKENS (bounded — "
            "the set is small).\n"
        )
    env = os.environ.copy()
    if use_local_source_corpus and "GVF_CORPUS_DIR" not in env:
        env["GVF_CORPUS_DIR"] = str(LOCAL_SOURCE_CORPUS)
        print(f"Using benchmark source cache: {LOCAL_SOURCE_CORPUS.relative_to(REPO)}")
    elif "GVF_CORPUS_DIR" in env:
        print(f"Using GVF_CORPUS_DIR from environment: {env['GVF_CORPUS_DIR']}")
    else:
        print("Using repo/global corpus cache (benchmark sources/ is incomplete).")

    for gene in genes:
        pmid_file = HERE / "pmids" / f"{gene}.txt"
        if not pmid_file.exists():
            print(f"  (skip {gene}: no pmids/{gene}.txt)")
            continue
        out = outroot / gene
        cmd = [
            python_exe(),
            "-m",
            "cli",
            "gvf-run",
            gene,
            "--email",
            email,
            "--pmid-file",
            str(pmid_file),
            "--source-recovery-timeout-s",
            str(source_recovery_timeout_s),
            "--output",
            str(out),
        ]
        if fast:
            cmd.append("--no-source-recovery")
        print("•", " ".join(cmd))
        res = subprocess.run(cmd, cwd=str(REPO), env=env)
        if res.returncode != 0:
            print(f"  (gvf-run failed for {gene}, exit {res.returncode}; skipping)")
            continue
        cands = sorted(out.glob(f"**/{gene}*.db"), key=lambda p: p.stat().st_mtime)
        if cands:
            produced[gene] = cands[-1]
            print(f"  -> {cands[-1]}")
        else:
            print(f"  (no {gene}*.db produced under {out})")
    return produced


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--mode", choices=["score", "extract"], default="score")
    ap.add_argument(
        "--db",
        action="append",
        help="Score this DB: GENE=/path/to.db (repeatable). "
        "Default: the four canonical DBs.",
    )
    ap.add_argument("--email", help="NCBI email (required for --mode extract).")
    ap.add_argument(
        "--fast",
        action="store_true",
        help="extract mode: use the quick --no-source-recovery path "
        "(NON-representative; skips source-acquisition refresh).",
    )
    ap.add_argument(
        "--source-recovery-timeout-s",
        type=int,
        default=120,
        help="extract mode: per-paper timeout forwarded to gvf-run source recovery.",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=HERE / ".last_run",
        help="Where scorer artifacts go (gitignored).",
    )
    ap.add_argument(
        "--extract-outdir",
        type=Path,
        default=HERE / ".extract_runs",
        help="Where --mode extract writes gvf-run output (gitignored).",
    )
    ap.add_argument(
        "--write-baseline",
        action="store_true",
        help="Save this run's headline numbers as expected_baseline.json.",
    )
    ap.add_argument(
        "--fail-on-regression",
        action="store_true",
        help="Exit nonzero if recall drops or count MAE rises beyond tolerance "
        "vs expected_baseline.json. Makes this runner a CI/nightly gate.",
    )
    ap.add_argument(
        "--regression-tol",
        type=float,
        default=0.005,
        help="Recall regression tolerance (fraction; default 0.005 = 0.5pp).",
    )
    ap.add_argument(
        "--mae-tol",
        type=float,
        default=0.05,
        help="Count-MAE regression tolerance (absolute; default 0.05, matching "
        "the refresh_recall land gate).",
    )
    args = ap.parse_args()

    manifest = load_manifest()
    genes_in_set = genes_from_manifest(manifest)

    if args.mode == "extract":
        if not args.email:
            raise SystemExit("--mode extract requires --email")
        dbs = do_extract(
            args.email,
            args.extract_outdir,
            genes_in_set,
            args.fast,
            args.source_recovery_timeout_s,
            complete_local_source_corpus(manifest),
        )
        if not dbs:
            raise SystemExit("extract produced no DBs to score.")
    else:
        dbs = resolve_dbs(args.db)

    print("Scoring DBs:")
    for g, p in sorted(dbs.items()):
        print(f"  {g}: {p}")

    summary = run_scorer(dbs, args.outdir)
    papers = per_paper_table(args.outdir, manifest)
    baseline = json.loads(BASELINE.read_text()) if BASELINE.exists() else None
    print_scorecard(summary, papers, baseline)

    if args.write_baseline:
        BASELINE.write_text(json.dumps(headline(summary), indent=2) + "\n")
        print(f"Wrote baseline -> {display_path(BASELINE)}")

    print(
        f"Full scorer artifacts: {display_path(args.outdir)}/ "
        "(summary.json, report.md, per-gene discrepancies.csv, "
        "paper_disagreement_report.csv)"
    )

    if args.fail_on_regression:
        problems = gate_result(summary, baseline, args.regression_tol, args.mae_tol)
        if problems:
            print("\n❌ REGRESSION gate failed:")
            for problem in problems:
                print(f"  - {problem}")
            return 1
        print("\n✓ No regression vs baseline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
