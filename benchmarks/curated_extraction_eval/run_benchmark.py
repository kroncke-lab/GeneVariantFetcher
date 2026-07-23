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
          default it scores every canonical DB registered below, so
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
import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
GOLD_DIR = HERE / "gold"
MANIFEST = HERE / "manifest.csv"
BASELINE = HERE / "expected_baseline.json"
SCORER = REPO / "scripts" / "run_recall_suite.py"
LOCAL_SOURCE_CORPUS = HERE / "sources"
BASELINE_SCHEMA_VERSION = 2
EXIT_STAGE_WARNINGS = 3

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
        # The benchmark scores its own frozen ./gold subset, which is independent
        # of the live Azure review-gold cache. Pin the scorer to the full gene
        # scope and disable live sync so the non-cardiac arm (BRCA1/BRCA2/MYBPC3/
        # APOE) is scored too; the scorer's default --review-gold-tier=cardiac
        # (added in PR #163) would otherwise silently drop every non-cardiac gene.
        "--review-gold-sync",
        "off",
        "--review-gold-tier",
        "all",
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


def gene_scope(genes: Iterable[str]) -> tuple[str, ...]:
    """Return the canonical identity for an exact benchmark gene set."""
    return tuple(sorted({gene.strip().upper() for gene in genes if gene.strip()}))


def baseline_scope_key(genes: Iterable[str]) -> str:
    """Stable, human-readable key used by expected_baseline.json schema v2."""
    return ",".join(gene_scope(genes))


def require_pmid_files_match_manifest(
    genes: Iterable[str],
    manifest: dict[str, dict[str, str]],
    pmid_dir: Path | None = None,
) -> None:
    """Refuse extraction when a selected PMID file has drifted from the manifest."""
    wanted = gene_scope(genes)
    expected = {gene: set() for gene in wanted}
    for row in manifest.values():
        gene = row.get("gene", "").strip().upper()
        if gene in expected:
            expected[gene].add(row.get("pmid", "").strip())

    root = pmid_dir or HERE / "pmids"
    problems: list[str] = []
    for gene in wanted:
        path = root / f"{gene}.txt"
        if not path.is_file():
            problems.append(f"{gene}: missing {path}")
            continue
        observed_list = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        observed: set[str] = set()
        duplicate_set: set[str] = set()
        for pmid in observed_list:
            if pmid in observed:
                duplicate_set.add(pmid)
            else:
                observed.add(pmid)
        duplicates = sorted(duplicate_set)
        missing = sorted(expected[gene] - observed)
        unexpected = sorted(observed - expected[gene])
        details = []
        if missing:
            details.append("missing=" + ",".join(missing))
        if unexpected:
            details.append("unexpected=" + ",".join(unexpected))
        if duplicates:
            details.append("duplicates=" + ",".join(duplicates))
        if details:
            problems.append(f"{gene}: " + " ".join(details))

    if problems:
        raise SystemExit(
            "benchmark PMID files do not match the selected manifest: "
            + "; ".join(problems)
        )


def fixture_sha256(genes: Iterable[str], manifest: dict[str, dict[str, str]]) -> str:
    """Fingerprint the selected manifest rows and frozen gold files.

    Gene scope alone is insufficient: adding a PMID or changing a gold row under
    the same genes changes the comparison population. Profiles record this hash
    so such fixture changes fail closed until their baseline is intentionally
    regenerated.
    """
    wanted = set(gene_scope(genes))
    hasher = hashlib.sha256()
    selected_rows = sorted(
        (
            key,
            row,
        )
        for key, row in manifest.items()
        if row.get("gene", "").strip().upper() in wanted
    )
    for key, row in selected_rows:
        hasher.update(key.encode("utf-8"))
        hasher.update(b"\0")
        hasher.update(
            json.dumps(row, sort_keys=True, separators=(",", ":")).encode("utf-8")
        )
        hasher.update(b"\0")
    for gene in sorted(wanted):
        gold_path = GOLD_DIR / "normalized" / f"{gene}_recall_input.csv"
        hasher.update(gene.encode("utf-8"))
        hasher.update(b"\0")
        hasher.update(gold_path.read_bytes() if gold_path.is_file() else b"<MISSING>")
        hasher.update(b"\0")
    return hasher.hexdigest()


def select_baseline_profile(
    document: dict | None,
    genes: Iterable[str],
    all_genes: Iterable[str],
    fixture_hash: str | None = None,
) -> tuple[dict | None, str | None]:
    """Select only a baseline whose gene set exactly matches ``genes``.

    Schema-v1 files had no scope metadata. They can still represent the full
    manifest, but must never be applied to a filtered run: that was the source
    of the misleading all-eight-vs-cardiac comparison this harness used to make.
    """
    if not isinstance(document, dict):
        return None, "expected_baseline.json is missing or invalid"

    wanted = gene_scope(genes)
    if document.get("schema_version") == BASELINE_SCHEMA_VERSION:
        profiles = document.get("profiles")
        if not isinstance(profiles, dict):
            return None, "expected_baseline.json schema v2 has no profiles object"

        profile = profiles.get(baseline_scope_key(wanted))
        if (
            not isinstance(profile, dict)
            or gene_scope(profile.get("genes", [])) != wanted
        ):
            return (
                None,
                "no exact baseline profile for genes " + ",".join(wanted),
            )
        stored_hash = profile.get("fixture_sha256")
        if fixture_hash is not None and stored_hash != fixture_hash:
            return (
                None,
                "baseline fixture fingerprint differs for genes " + ",".join(wanted),
            )
        return profile, None

    # Legacy baseline: its committed meaning was the entire manifest. Preserve
    # that compatibility only for that exact scope, never for a subset.
    if "recall" in document and wanted == gene_scope(all_genes):
        return (
            document,
            "using legacy all-gene baseline; rewrite it with --write-baseline",
        )
    if "recall" in document:
        return (
            None,
            "legacy baseline has no gene scope and cannot be compared with subset "
            + ",".join(wanted),
        )
    return None, "expected_baseline.json has an unsupported schema"


def update_baseline_document(
    document: dict | None,
    genes: Iterable[str],
    metrics: dict,
    all_genes: Iterable[str],
    fixture_hash: str | None = None,
) -> dict:
    """Return schema v2 with only the selected exact-scope profile replaced."""
    profiles: dict[str, dict] = {}
    if isinstance(document, dict) and document.get("schema_version") == 2:
        existing = document.get("profiles")
        if isinstance(existing, dict):
            profiles = {
                str(key): value
                for key, value in existing.items()
                if isinstance(value, dict)
            }
    elif isinstance(document, dict) and "recall" in document:
        # A v1 file was the full-manifest baseline. Retain it when a subset is
        # written so --write-baseline never erases a different scope.
        full_scope = gene_scope(all_genes)
        profiles[baseline_scope_key(full_scope)] = {
            "genes": list(full_scope),
            **document,
        }

    wanted = gene_scope(genes)
    profile = {"genes": list(wanted), **metrics}
    if fixture_hash is not None:
        profile["fixture_sha256"] = fixture_hash
    profiles[baseline_scope_key(wanted)] = profile
    return {
        "schema_version": BASELINE_SCHEMA_VERSION,
        "profiles": {key: profiles[key] for key in sorted(profiles)},
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
        base_block = base_rec.get(name) or {}
        denominator_matches = (
            base_block.get("gold") is None
            or b["gold"] is None
            or base_block["gold"] == b["gold"]
        )
        if not denominator_matches:
            delta = (
                f"  (NOT COMPARABLE: baseline gold={base_block['gold']}, "
                f"current gold={b['gold']})"
            )
        if (
            name in base_rec
            and denominator_matches
            and b["recall"] is not None
            and base_block.get("recall") is not None
        ):
            d = b["recall"] - base_block["recall"]
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
        base_block = base_rec.get(name) or {}
        base = base_block.get("recall")
        if base is None:
            continue  # dimension not tracked by the baseline
        if (
            base_block.get("gold") is not None
            and block["gold"] is not None
            and base_block["gold"] != block["gold"]
        ):
            problems.append(
                f"recall.{name}: current gold denominator {block['gold']} != "
                f"baseline {base_block['gold']}"
            )
            continue
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
    degraded_runs: dict[str, str] | None = None,
) -> list[str]:
    """The --fail-on-regression decision, fail-closed. A missing baseline is a
    hard failure: a gate that compared against nothing must not report success.
    """
    problems: list[str] = []
    if baseline is None:
        problems.append(
            "no baseline at expected_baseline.json — refusing to pass a gate that "
            "checked nothing (write one with --write-baseline)"
        )
    else:
        problems.extend(check_regression(summary, baseline, recall_tol, mae_tol))
    for gene, reason in sorted((degraded_runs or {}).items()):
        problems.append(f"extract.{gene}: degraded run ({reason})")
    return problems


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


@dataclass
class ExtractionBatch:
    """DBs that may be scored, plus explicit quality/failure state per gene."""

    dbs: dict[str, Path] = field(default_factory=dict)
    degraded: dict[str, str] = field(default_factory=dict)
    failures: dict[str, str] = field(default_factory=dict)


def _file_signature(path: Path, *, include_content: bool = False) -> tuple:
    stat = path.stat()
    digest = hashlib.sha256(path.read_bytes()).hexdigest() if include_content else None
    return (stat.st_mtime_ns, stat.st_size, stat.st_ino, digest)


def _snapshot_files(
    root: Path, pattern: str, *, include_content: bool = False
) -> dict[Path, tuple]:
    if not root.exists():
        return {}
    snapshot: dict[Path, tuple] = {}
    for path in root.glob(pattern):
        if not path.is_file():
            continue
        try:
            snapshot[path.resolve()] = _file_signature(
                path, include_content=include_content
            )
        except OSError:
            continue
    return snapshot


def _fresh_files(
    root: Path,
    pattern: str,
    before: dict[Path, tuple],
    *,
    include_content: bool = False,
) -> list[Path]:
    """Files newly created or modified after one subprocess invocation."""
    fresh: list[tuple[Path, int]] = []
    for path in root.glob(pattern) if root.exists() else []:
        if not path.is_file():
            continue
        resolved = path.resolve()
        try:
            signature = _file_signature(path, include_content=include_content)
        except OSError:
            continue
        if before.get(resolved) != signature:
            fresh.append((resolved, signature[0]))
    fresh.sort(key=lambda item: item[1], reverse=True)
    return [path for path, _ in fresh]


def _status_exit_code(payload: dict) -> int | None:
    raw = payload.get("exit_code")
    if raw is None:
        return None
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def _load_fresh_status(
    gene: str, status_paths: list[Path]
) -> tuple[dict | None, Path | None, str | None]:
    """Load the newest fresh, valid status for ``gene`` without using stale state."""
    malformed: list[str] = []
    for path in status_paths:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            malformed.append(f"{path}: {exc}")
            continue
        if not isinstance(payload, dict):
            malformed.append(f"{path}: root is not an object")
            continue
        status_gene = str(payload.get("gene") or "").strip().upper()
        if not status_gene:
            malformed.append(f"{path}: gene field is absent")
            continue
        if status_gene != gene:
            malformed.append(f"{path}: reports gene {status_gene}, expected {gene}")
            continue
        return payload, path, None
    if malformed:
        return (
            None,
            None,
            "no valid fresh RUN_STATUS.json (" + "; ".join(malformed) + ")",
        )
    return None, None, None


def _classify_status(
    payload: dict | None,
    actual_exit: int,
) -> tuple[bool, str | None, str | None]:
    """Return (scoreable, degraded reason, failure reason)."""
    if payload is None:
        if actual_exit != 0:
            return (
                False,
                None,
                f"gvf-run exited {actual_exit} without a fresh RUN_STATUS.json",
            )
        return (
            True,
            "exit 0 but no fresh RUN_STATUS.json; using fresh DB discovery",
            None,
        )

    status = str(payload.get("status") or "").strip().lower()
    severity = str(payload.get("severity") or "").strip().lower()
    declared_exit = _status_exit_code(payload)
    stage_failures = payload.get("stage_failures")

    if payload.get("exit_code") is not None and declared_exit is None:
        return False, None, f"RUN_STATUS exit_code is invalid: {payload['exit_code']!r}"
    if severity and severity not in {"ok", "warning"}:
        return False, None, f"RUN_STATUS severity is unscoreable: {severity!r}"
    if status and status not in {"completed", "completed_with_warnings"}:
        return False, None, f"RUN_STATUS status is unscoreable: {status!r}"
    if declared_exit is not None and declared_exit not in {0, EXIT_STAGE_WARNINGS}:
        return False, None, f"RUN_STATUS exit_code is fatal: {declared_exit}"
    if actual_exit not in {0, EXIT_STAGE_WARNINGS}:
        return False, None, f"gvf-run exited fatally with {actual_exit}"
    if declared_exit is not None and declared_exit != actual_exit:
        return (
            False,
            None,
            f"process exit {actual_exit} disagrees with RUN_STATUS exit {declared_exit}",
        )

    warning_signals = []
    if not status:
        warning_signals.append("status field absent")
    if not severity:
        warning_signals.append("severity field absent")
    if declared_exit is None:
        warning_signals.append("status exit_code absent")
    if status == "completed_with_warnings":
        warning_signals.append("status=completed_with_warnings")
    if severity == "warning":
        warning_signals.append("severity=warning")
    if declared_exit == EXIT_STAGE_WARNINGS:
        warning_signals.append(f"status exit={EXIT_STAGE_WARNINGS}")
    if actual_exit == EXIT_STAGE_WARNINGS:
        warning_signals.append(f"process exit={EXIT_STAGE_WARNINGS}")

    if actual_exit == EXIT_STAGE_WARNINGS and not (
        status == "completed_with_warnings" or severity == "warning"
    ):
        return (
            False,
            None,
            f"gvf-run exited {EXIT_STAGE_WARNINGS} but RUN_STATUS does not mark warnings",
        )
    if status == "completed" and severity == "warning":
        return False, None, "RUN_STATUS has contradictory completed/warning fields"
    if status == "completed_with_warnings" and severity == "ok":
        return False, None, "RUN_STATUS has contradictory warning/ok fields"
    warning_status = status == "completed_with_warnings" or severity == "warning"
    if warning_status and actual_exit != EXIT_STAGE_WARNINGS:
        return (
            False,
            None,
            "RUN_STATUS marks a warning completion but process exit is not "
            f"{EXIT_STAGE_WARNINGS}",
        )
    if isinstance(stage_failures, list) and stage_failures and not warning_status:
        return (
            False,
            None,
            "RUN_STATUS has stage_failures without warning completion fields",
        )

    if warning_signals:
        if isinstance(stage_failures, list) and stage_failures:
            detail = "; ".join(str(item) for item in stage_failures)
        else:
            detail = ", ".join(warning_signals)
        return True, detail, None
    return True, None, None


def _resolve_status_db(
    status_path: Path, payload: dict
) -> tuple[Path | None, str | None]:
    raw = payload.get("active_db")
    if raw is None or not str(raw).strip():
        return None, "fresh RUN_STATUS.json has no active_db"
    candidate = Path(str(raw)).expanduser()
    if not candidate.is_absolute():
        candidate = status_path.parent / candidate
    candidate = candidate.resolve()
    if not candidate.is_file():
        return None, f"RUN_STATUS active_db does not exist: {candidate}"
    return candidate, None


def require_complete_gene_set(dbs: dict[str, Path], requested: Iterable[str]) -> None:
    """Refuse partial aggregate scores, which change the denominator silently."""
    wanted = set(gene_scope(requested))
    missing = sorted(
        gene for gene in wanted if gene not in dbs or not dbs[gene].is_file()
    )
    if missing:
        raise SystemExit(
            "refusing to score a partial benchmark; missing requested gene DB(s): "
            + ", ".join(missing)
        )


def require_clean_baseline_write(degraded_runs: dict[str, str]) -> None:
    """Do not let an incomplete extraction run move a committed baseline."""
    if degraded_runs:
        genes = ", ".join(sorted(degraded_runs))
        raise SystemExit(
            "refusing to write baseline from degraded extraction run(s): " + genes
        )


def do_extract(
    email: str,
    outroot: Path,
    genes: list[str],
    fast: bool,
    source_recovery_timeout_s: int,
    use_local_source_corpus: bool,
) -> ExtractionBatch:
    """Process the set through the standard `gvf gvf-run` pipeline, return DBs.

    Default = the full, regular, default pipeline. `fast` adds
    --no-source-recovery for a quick, non-representative pass that skips
    source-acquisition refresh.
    """
    outroot.mkdir(parents=True, exist_ok=True)
    batch = ExtractionBatch()
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
            reason = f"no pmids/{gene}.txt"
            batch.failures[gene] = reason
            print(f"  ({gene} failed: {reason})")
            continue
        out = outroot / gene
        before_status = _snapshot_files(out, "**/RUN_STATUS.json", include_content=True)
        before_dbs = _snapshot_files(out, f"**/{gene}*.db")
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

        fresh_statuses = _fresh_files(
            out,
            "**/RUN_STATUS.json",
            before_status,
            include_content=True,
        )
        payload, status_path, status_error = _load_fresh_status(gene, fresh_statuses)
        if status_error:
            batch.failures[gene] = status_error
            print(f"  ({gene} failed: {status_error})")
            continue

        scoreable, degraded_reason, failure_reason = _classify_status(
            payload, res.returncode
        )
        if not scoreable:
            assert failure_reason is not None
            batch.failures[gene] = failure_reason
            print(f"  ({gene} failed: {failure_reason})")
            continue

        db: Path | None = None
        if payload is not None and status_path is not None:
            db, db_error = _resolve_status_db(status_path, payload)
            if db_error:
                batch.failures[gene] = db_error
                print(f"  ({gene} failed: {db_error})")
                continue
        if db is None:
            fresh_dbs = _fresh_files(out, f"**/{gene}*.db", before_dbs)
            if fresh_dbs:
                db = fresh_dbs[0]

        if db is not None:
            batch.dbs[gene] = db
            degradation = degraded_reason
            if fast:
                fast_reason = "--fast disables source recovery (non-representative)"
                degradation = (
                    f"{degradation}; {fast_reason}" if degradation else fast_reason
                )
            if degradation:
                batch.degraded[gene] = degradation
                print(f"  -> {db} [DEGRADED: {degradation}]")
            else:
                print(f"  -> {db}")
        else:
            reason = f"no fresh {gene}*.db produced under {out}"
            batch.failures[gene] = reason
            print(f"  ({gene} failed: {reason})")
    return batch


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--mode", choices=["score", "extract"], default="score")
    ap.add_argument(
        "--db",
        action="append",
        help="Score this DB: GENE=/path/to.db (repeatable). "
        "Default: every canonical DB registered for the selected gene scope.",
    )
    ap.add_argument("--email", help="NCBI email (required for --mode extract).")
    ap.add_argument(
        "--genes",
        help="Comma-separated gene subset to extract and/or score (e.g. "
        "KCNH2,KCNQ1,SCN5A,RYR2). Default: every gene in the set. Recall / MAE / "
        "precision are only trustworthy for genes with a fully human-curated gold "
        "standard (the four cardiac genes); the non-cardiac genes use "
        "curator/derived gold_overrides, so restrict scoring to the cardiac four "
        "for headline metrics.",
    )
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
    all_genes_in_set = genes_from_manifest(manifest)
    genes_in_set = list(all_genes_in_set)

    gene_filter: set[str] | None = None
    if args.genes:
        gene_filter = {g.strip().upper() for g in args.genes.split(",") if g.strip()}
        unknown = sorted(gene_filter - set(all_genes_in_set))
        if unknown:
            raise SystemExit(
                f"--genes {args.genes!r} includes gene(s) outside the set: "
                + ", ".join(unknown)
            )
        genes_in_set = [g for g in genes_in_set if g in gene_filter]
        if not genes_in_set:
            raise SystemExit(f"--genes {args.genes!r} matched no gene in the set")
    selected_genes = set(genes_in_set)
    manifest_in_scope = {
        key: row
        for key, row in manifest.items()
        if row["gene"].strip().upper() in selected_genes
    }
    current_fixture_hash = fixture_sha256(genes_in_set, manifest_in_scope)

    degraded_runs: dict[str, str] = {}
    if args.mode == "extract":
        if not args.email:
            raise SystemExit("--mode extract requires --email")
        require_pmid_files_match_manifest(genes_in_set, manifest_in_scope)
        extracted = do_extract(
            args.email,
            args.extract_outdir,
            genes_in_set,
            args.fast,
            args.source_recovery_timeout_s,
            complete_local_source_corpus(manifest_in_scope),
        )
        dbs = extracted.dbs
        degraded_runs = extracted.degraded
        if extracted.failures:
            details = "; ".join(
                f"{gene}: {reason}"
                for gene, reason in sorted(extracted.failures.items())
            )
            print(f"Extract failures: {details}")
    else:
        dbs = resolve_dbs(args.db)
        dbs = {g: p for g, p in dbs.items() if g in selected_genes}

    require_complete_gene_set(dbs, genes_in_set)

    print("Scoring DBs:")
    for g, p in sorted(dbs.items()):
        print(f"  {g}: {p}")

    summary = run_scorer(dbs, args.outdir)
    papers = per_paper_table(args.outdir, manifest_in_scope)
    baseline_document = (
        json.loads(BASELINE.read_text(encoding="utf-8")) if BASELINE.exists() else None
    )
    baseline, baseline_note = select_baseline_profile(
        baseline_document,
        genes_in_set,
        all_genes_in_set,
        current_fixture_hash,
    )
    if baseline_note:
        print(f"Baseline: {baseline_note}")
    print_scorecard(summary, papers, baseline)

    if degraded_runs:
        print("⚠️  Score includes degraded extraction run(s):")
        for gene, reason in sorted(degraded_runs.items()):
            print(f"  - {gene}: {reason}")

    if args.write_baseline:
        require_clean_baseline_write(degraded_runs)
        updated = update_baseline_document(
            baseline_document,
            genes_in_set,
            headline(summary),
            all_genes_in_set,
            current_fixture_hash,
        )
        BASELINE.write_text(json.dumps(updated, indent=2) + "\n", encoding="utf-8")
        print(f"Wrote baseline -> {display_path(BASELINE)}")

    print(
        f"Full scorer artifacts: {display_path(args.outdir)}/ "
        "(summary.json, report.md, per-gene discrepancies.csv, "
        "paper_disagreement_report.csv)"
    )

    if args.fail_on_regression:
        problems = gate_result(
            summary,
            baseline,
            args.regression_tol,
            args.mae_tol,
            degraded_runs,
        )
        if problems:
            print("\n❌ REGRESSION gate failed:")
            for problem in problems:
                print(f"  - {problem}")
            return 1
        print("\n✓ No regression vs baseline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
