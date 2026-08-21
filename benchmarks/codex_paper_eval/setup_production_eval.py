#!/usr/bin/env python3
"""Create and validate a blinded production ``gold_120`` evaluation scaffold.

The live cohort is 119 gene-paper attempts after the KCNH2 PMID erratum.  This
tool verifies the registry-pinned manifest, runs the existing gold-value-blind
``prepare`` step, writes exact per-gene PMID inputs, and emits separate scripts
for extraction and for the post-extraction lock/score boundary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import stat
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
TIERS = REPO / "benchmarks" / "evaluation_tiers"
DEFAULT_MANIFEST = TIERS / "tier2_gold_120.tsv"
DEFAULT_REGISTRY = TIERS / "registry.json"
DEFAULT_RUNS_DIR = HERE / "runs"
DEFAULT_CORPUS = REPO / "corpus"
DEFAULT_GOLD = REPO / "gene_variant_fetcher_gold_standard" / "normalized"
DEFAULT_SEED = 2026081301
PROVIDER_KEYS = ("OPENAI_API_KEY", "AZURE_AI_API_KEY", "ANTHROPIC_API_KEY")

RUNTIME_TREES = (
    "cli",
    "config",
    "gene_literature",
    "harvesting",
    "pipeline",
    "scripts",
    "utils",
)
RUNTIME_FILES = ("pyproject.toml", "requirements.lock", "uv.lock")


class SetupError(RuntimeError):
    """The evaluation scaffold is incomplete or no longer reproducible."""


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def read_manifest(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"[\t, ]+", line)
        if len(parts) != 2 or not parts[1].isdigit():
            raise SetupError(f"{path}:{line_number}: expected '<GENE> <PMID>'")
        rows.append((parts[0], parts[1]))
    if len(rows) != len(set(rows)):
        raise SetupError(f"{path}: duplicate gene/PMID attempts")
    return rows


def cohort_contract(manifest: Path, registry_path: Path) -> dict:
    registry = json.loads(registry_path.read_text())
    matches = [
        tier
        for tier in registry.get("tiers", [])
        if tier.get("manifest") == manifest.name
    ]
    if len(matches) != 1:
        raise SetupError(
            f"{manifest.name}: expected exactly one entry in {registry_path}"
        )
    tier = matches[0]
    rows = read_manifest(manifest)
    actual_digest = digest(manifest)
    actual_counts = Counter(gene for gene, _ in rows)
    errors: list[str] = []
    if actual_digest != tier.get("sha256"):
        errors.append(f"SHA-256 {actual_digest} != registry {tier.get('sha256')}")
    if len(rows) != tier.get("attempt_count"):
        errors.append(f"attempts {len(rows)} != registry {tier.get('attempt_count')}")
    unique_pmids = len({pmid for _, pmid in rows})
    if unique_pmids != tier.get("unique_pmid_count"):
        errors.append(
            f"unique PMIDs {unique_pmids} != registry {tier.get('unique_pmid_count')}"
        )
    if actual_counts != Counter(tier.get("gene_attempt_counts", {})):
        errors.append(
            f"gene counts {dict(actual_counts)} != registry "
            f"{tier.get('gene_attempt_counts')}"
        )
    if errors:
        raise SetupError(f"cohort registry mismatch: {'; '.join(errors)}")
    return {**tier, "rows": rows, "registry": str(registry_path.resolve())}


def runtime_source_files() -> list[Path]:
    files: set[Path] = set()
    for relative in RUNTIME_TREES:
        root = REPO / relative
        if not root.is_dir():
            raise SetupError(f"runtime source tree missing: {root}")
        files.update(
            path
            for path in root.rglob("*")
            if path.is_file()
            and "__pycache__" not in path.parts
            and path.suffix not in {".pyc", ".pyo"}
        )
    files.update(
        path for path in (REPO / name for name in RUNTIME_FILES) if path.is_file()
    )
    files.update(path for path in HERE.glob("*.py") if path.is_file())
    return sorted(files)


def runtime_fingerprint() -> dict:
    value = hashlib.sha256()
    files = runtime_source_files()
    for path in files:
        relative = path.relative_to(REPO).as_posix()
        value.update(relative.encode())
        value.update(b"\0")
        value.update(digest(path).encode())
        value.update(b"\n")
    return {"sha256": value.hexdigest(), "file_count": len(files)}


def git_state() -> dict:
    def run(*arguments: str) -> str | None:
        result = subprocess.run(
            ["git", *arguments],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        return result.stdout.strip() if result.returncode == 0 else None

    status = run("status", "--porcelain=v1", "--untracked-files=normal")
    return {
        "head": run("rev-parse", "HEAD"),
        "dirty": bool(status),
        "status_sha256": hashlib.sha256((status or "").encode()).hexdigest(),
        "runtime_source": runtime_fingerprint(),
    }


def configured_environment_names() -> set[str]:
    names = set(os.environ)
    dotenv = REPO / ".env"
    if dotenv.is_file():
        for raw in dotenv.read_text(errors="replace").splitlines():
            match = re.match(r"\s*(?:export\s+)?([A-Z][A-Z0-9_]*)\s*=", raw)
            if match:
                names.add(match.group(1))
    return names


def shell_quote(value: Path | str) -> str:
    return shlex.quote(str(value))


def runtime_python() -> Path:
    """Return the active environment entrypoint without dereferencing it."""
    return Path(sys.executable).absolute()


def make_extraction_script(
    *, run_dir: Path, genes: list[str], email: str, python: Path
) -> str:
    setup = HERE / "setup_production_eval.py"
    lines = [
        "#!/bin/zsh",
        "set -euo pipefail",
        f"cd {shell_quote(REPO)}",
        f"PYTHON={shell_quote(python)}",
        f"RUN_DIR={shell_quote(run_dir)}",
        f"EMAIL={shell_quote(email)}",
        f'"$PYTHON" {shell_quote(setup)} check --run-dir "$RUN_DIR"',
        'mkdir -p "$RUN_DIR/operator_logs"',
        "typeset -a pids",
        "",
    ]
    for gene in genes:
        lines.extend(
            [
                f'echo "Starting calibrated {gene} extraction"',
                "(",
                '  "$PYTHON" -m cli gvf-run '
                f'{gene} --email "$EMAIL" --output "$RUN_DIR/production_runs" '
                f'--pmid-file "$RUN_DIR/pmids/{gene}.txt" '
                "--no-source-recovery --no-corpus-sync --no-publish-review "
                f'2>&1 | tee "$RUN_DIR/operator_logs/{gene}.log"',
                ") &",
                "pids+=($!)",
                "",
            ]
        )
    lines.extend(
        [
            "failed=0",
            'for pid in "${pids[@]}"; do',
            '  if ! wait "$pid"; then',
            "    failed=1",
            "  fi",
            "done",
            "if (( failed )); then",
            '  echo "At least one gene extraction failed; do not lock or score." >&2',
            "  exit 1",
            "fi",
        ]
    )
    return "\n".join(lines).rstrip() + "\n"


def make_finalize_script(*, run_dir: Path, python: Path) -> str:
    setup = HERE / "setup_production_eval.py"
    rebind = HERE / "rebind_production_sources.py"
    converter = HERE / "db_to_predictions.py"
    evaluator = HERE / "run_eval.py"
    return f"""#!/bin/zsh
set -euo pipefail
cd {shell_quote(REPO)}
PYTHON={shell_quote(python)}
RUN_DIR={shell_quote(run_dir)}
"$PYTHON" {shell_quote(setup)} check --run-dir "$RUN_DIR"
"$PYTHON" {shell_quote(rebind)} \\
  --run-dir "$RUN_DIR" --production-root "$RUN_DIR/production_runs"
"$PYTHON" {shell_quote(converter)} \\
  --run-dir "$RUN_DIR" --production-root "$RUN_DIR/production_runs" \\
  --trust-mode trusted --identity-mode trusted \\
  --out "$RUN_DIR/predictions.json"
"$PYTHON" {shell_quote(evaluator)} lock --run-dir "$RUN_DIR"
"$PYTHON" {shell_quote(evaluator)} score --run-dir "$RUN_DIR"
"""


def write_executable(path: Path, content: str) -> None:
    path.write_text(content)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_runbook(run_dir: Path, contract: dict, genes: list[str]) -> None:
    counts = contract["gene_attempt_counts"]
    count_text = ", ".join(f"{gene} {counts[gene]}" for gene in genes)
    (run_dir / "RUNBOOK.md").write_text(
        f"""# Current-code gold_120 production evaluation

This scaffold is pinned to `{Path(contract["manifest"]).name}` at
`{contract["sha256"]}`: **{contract["attempt_count"]} gene-paper attempts** /
**{contract["unique_pmid_count"]} unique PMIDs** ({count_text}). The name remains
`gold_120` for continuity; KCNH2 PMID 10086972 is intentionally absent after the
erratum identified it as the wrong paper.

## 1. Extract without opening gold values

```bash
{run_dir / "run_extraction.sh"}
```

The four commands use production `gvf-run` with exact PMID files. Source
recovery and corpus sync are disabled so this is a calibrated comparison over
the frozen source-available cohort; publication is explicitly disabled. The
four gene processes run concurrently, matching the accepted 2026-08-20 test;
each writes a separate `operator_logs/<GENE>.log`, and the launcher fails if any
gene process fails.

## 2. Inspect production completion

Require one successful `RUN_STATUS.json`, final database, and finalized
write-time-verified trace manifest for every gene before proceeding.

## 3. Rebind exact inputs, project, lock, then score

```bash
{run_dir / "lock_and_score.sh"}
```

This second script binds the exact run-local source material and production
trace manifests into `predictions.json`, applies the collaborator-facing trusted
count and identity projection, locks predictions before any gold value is read,
and only then invokes the scorer. A raw (`--trust-mode all --identity-mode all`)
projection may be generated later only as a clearly labeled diagnostic; it must
not replace the locked primary.
"""
    )


def create(args: argparse.Namespace) -> Path:
    manifest = args.paper_manifest.resolve()
    registry = args.registry.resolve()
    contract = cohort_contract(manifest, registry)
    if contract.get("id") != "gold_120":
        raise SetupError(
            f"production setup expects registry tier gold_120, got {contract.get('id')}"
        )
    if int(contract.get("selection_seed")) != args.seed:
        raise SetupError(
            f"seed {args.seed} != registry selection seed {contract.get('selection_seed')}"
        )

    # Capture repository state before ``command_prepare`` creates the run
    # directory. Otherwise the scaffold marks a clean source commit dirty merely
    # because the scaffold itself is untracked.
    state = git_state()

    # Import only after the independent manifest/registry check above.
    from benchmarks.codex_paper_eval.run_eval import command_prepare

    command_prepare(
        SimpleNamespace(
            seed=args.seed,
            per_gene=30,
            minimum_chars=args.minimum_chars,
            corpus_root=args.corpus_root.resolve(),
            gold_root=args.gold_root.resolve(),
            paper_manifest=manifest,
            runs_dir=args.runs_dir.resolve(),
            legacy_source_selection=False,
            run_id=args.run_id,
        )
    )
    run_dir = (args.runs_dir / args.run_id).resolve()
    selection = json.loads((run_dir / "selection.json").read_text())
    selected = [(paper["gene"], str(paper["pmid"])) for paper in selection["papers"]]
    if selected != contract["rows"]:
        raise SetupError(
            "prepared selection does not exactly match the pinned manifest"
        )

    by_gene: dict[str, list[str]] = defaultdict(list)
    for gene, pmid in selected:
        by_gene[gene].append(pmid)
    genes = sorted(by_gene)
    pmid_dir = run_dir / "pmids"
    production_root = run_dir / "production_runs"
    pmid_dir.mkdir()
    production_root.mkdir()
    for gene, pmids in by_gene.items():
        (pmid_dir / f"{gene}.txt").write_text("\n".join(pmids) + "\n")

    # Preserve the virtual-environment entrypoint. Resolving the symlink can
    # select the dependency-free base interpreter and make the generated
    # launcher fail before importing the CLI.
    python = runtime_python()
    if not python.is_file():
        raise SetupError(f"Python interpreter is unavailable: {python}")
    environment_names = configured_environment_names()
    setup = {
        "schema_version": 1,
        "run_id": args.run_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "cohort": {
            **{
                key: contract[key]
                for key in (
                    "id",
                    "sha256",
                    "attempt_count",
                    "unique_pmid_count",
                    "gene_attempt_counts",
                    "selection_seed",
                    "registry",
                )
            },
            "manifest": str(manifest),
        },
        "route": {
            "driver": "python -m cli gvf-run",
            "calibrated_pmid_files": True,
            "source_recovery": False,
            "corpus_sync": False,
            "publish_review": False,
            "full_coverage": False,
            "parallel_gene_processes": 4,
            "primary_projection": "trusted_count_and_identity",
            "reason": (
                "Apples-to-apples current-default recall arm against the accepted "
                "20260820 gold119 production evaluation."
            ),
        },
        "blinding": selection["blinding"],
        "python": str(python),
        "provider_configuration_present": {
            key: key in environment_names for key in PROVIDER_KEYS
        },
        "azure_base_present": "AZURE_AI_API_BASE" in environment_names,
        "repository": state,
    }
    (run_dir / "setup.json").write_text(
        json.dumps(setup, indent=2, sort_keys=True) + "\n"
    )
    write_executable(
        run_dir / "run_extraction.sh",
        make_extraction_script(
            run_dir=run_dir, genes=genes, email=args.email, python=python
        ),
    )
    write_executable(
        run_dir / "lock_and_score.sh",
        make_finalize_script(run_dir=run_dir, python=python),
    )
    write_runbook(run_dir, {**contract, "manifest": str(manifest)}, genes)
    check_run(run_dir)
    return run_dir


def check_run(run_dir: Path) -> None:
    run_dir = run_dir.resolve()
    setup_path = run_dir / "setup.json"
    selection_path = run_dir / "selection.json"
    if not setup_path.is_file() or not selection_path.is_file():
        raise SetupError(f"incomplete setup under {run_dir}")
    setup = json.loads(setup_path.read_text())
    cohort = setup["cohort"]
    contract = cohort_contract(Path(cohort["manifest"]), Path(cohort["registry"]))
    selection = json.loads(selection_path.read_text())
    selected = [(paper["gene"], str(paper["pmid"])) for paper in selection["papers"]]
    if selected != contract["rows"]:
        raise SetupError("selection membership/order drifted from the pinned cohort")
    for gene, expected_count in contract["gene_attempt_counts"].items():
        pmid_path = run_dir / "pmids" / f"{gene}.txt"
        if not pmid_path.is_file():
            raise SetupError(f"missing calibrated PMID file: {pmid_path}")
        actual = [
            line.strip() for line in pmid_path.read_text().splitlines() if line.strip()
        ]
        expected = [pmid for row_gene, pmid in contract["rows"] if row_gene == gene]
        if actual != expected or len(actual) != expected_count:
            raise SetupError(f"{gene} PMID file drifted from the pinned cohort")
    current = runtime_fingerprint()
    prepared = setup["repository"]["runtime_source"]
    if current != prepared:
        raise SetupError(
            "runtime source changed after setup: create a new run scaffold before extraction "
            f"(prepared {prepared}, current {current})"
        )
    extraction = (run_dir / "run_extraction.sh").read_text()
    forbidden = ("run_eval.py score", "--gold-root", "--publish-review")
    if any(token in extraction for token in forbidden):
        raise SetupError("extraction script crosses the gold/publication boundary")
    required = (
        "--pmid-file",
        "--no-source-recovery",
        "--no-corpus-sync",
        "--no-publish-review",
    )
    if any(token not in extraction for token in required):
        raise SetupError("extraction script is missing a calibrated-run safeguard")


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    create_parser = sub.add_parser("create")
    create_parser.add_argument("--run-id", required=True)
    create_parser.add_argument("--email", required=True)
    create_parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    create_parser.add_argument("--minimum-chars", type=int, default=2000)
    create_parser.add_argument("--paper-manifest", type=Path, default=DEFAULT_MANIFEST)
    create_parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    create_parser.add_argument("--corpus-root", type=Path, default=DEFAULT_CORPUS)
    create_parser.add_argument("--gold-root", type=Path, default=DEFAULT_GOLD)
    create_parser.add_argument("--runs-dir", type=Path, default=DEFAULT_RUNS_DIR)
    check_parser = sub.add_parser("check")
    check_parser.add_argument("--run-dir", type=Path, required=True)
    return ap


def main() -> int:
    args = parser().parse_args()
    try:
        if args.command == "create":
            run_dir = create(args)
            print(run_dir)
        else:
            check_run(args.run_dir)
            print(f"setup valid: {args.run_dir.resolve()}")
    except SetupError as exc:
        print(f"setup error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
