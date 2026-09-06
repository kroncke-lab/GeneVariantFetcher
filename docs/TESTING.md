# Installation and runtime checks

Use these checks after installation or maintenance. They establish that the
software builds and runs; they do not establish extraction accuracy or promote
a recall headline. Current metrics and acceptance rules remain in
[RECALL_STATUS.md](RECALL_STATUS.md) and [TASKS.md](../TASKS.md).

## Installation and static checks

Follow [QUICKSTART.md](QUICKSTART.md) for Python 3.11+ and the editable install.
Test both entry points: repository-root tests can pass even when the installed
console script is missing.

```bash
.venv/bin/gvf --help
.venv/bin/gvf gvf-run --help
.venv/bin/python -m cli --help
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
```

Check dependency compatibility with `.venv/bin/python -m pip check`, or
`uv pip check --python .venv/bin/python` in a uv environment without pip.
Do not print API keys while checking configuration. `MODEL_PROVIDER` selects
the provider; a credential alone does not. Set `MODEL_PROVIDER=azure` with the
Azure credentials for this workstation's preferred route. Explicit per-stage
model overrides still win; the shipped unset-provider default is Anthropic.

Compile only tracked Python files, avoiding virtual environments and paper data:

```bash
.venv/bin/python - <<'PY'
import py_compile
import subprocess
import tempfile
from pathlib import Path

paths = [p for p in subprocess.check_output(
    ["git", "ls-files", "-z"], text=True
).split("\0") if p.endswith(".py")]
with tempfile.TemporaryDirectory(prefix="gvf-compile-") as target:
    for index, path in enumerate(paths):
        py_compile.compile(path, cfile=str(Path(target) / f"{index}.pyc"), doraise=True)
print(f"Compiled {len(paths)} tracked Python files")
PY
```

## Offline execution and regression checks

These are the CI suites; no paid extraction is needed:

```bash
.venv/bin/python -m pytest tests/unit -q
.venv/bin/python -m pytest tests/recall/test_bounded_e2e.py -q
.venv/bin/python -m pytest benchmarks/curated_extraction_eval/negative_cases -q
```

The bounded end-to-end fixture exercises SQLite, scoring and recovery-driver
wiring. The unit suite also exercises the full orchestrator with network and
model boundaries replaced by fixtures. Report actual test counts and failures;
there is no fixed historical count to treat as the current target.

The mixed-gold cost checks read tracked predictions and SHA-pinned usage-only
receipts in `tests/fixtures/cost_calibration_usage.json`. Full operator traces
under ignored `results/` are not a CI prerequisite. Re-export those receipts
only when the original calibration traces are available, using
`benchmarks/evaluation_tiers/export_cost_usage_receipts.py`; the exporter checks
source hashes and model totals without changing the historical cost profile.

For import coverage, skip `__main__` modules, which execute CLI argument parsing.
A smoke script must raise or exit nonzero on any import failure; printing
“All modules imported” after catching errors is not a passing check.

## Wheel and script checks

Build and install a wheel into a separate test environment, then run `gvf
--help` from outside the repository. This prevents the checkout from hiding
missing packaged modules or data. `.github/workflows/ci.yml` contains the
current isolated-wheel recipe and reference-resource assertions.

Start local builds with a clean generated `build/` staging directory. A reused
`build/lib` can put previously deleted modules back into an otherwise valid
wheel; archiving that generated directory before rebuilding avoids the stale
copy. Clearing only the wheel output directory does not clear staging. CI checks
each packaged source/data file against the checkout as well as running imports.

Check supported maintenance scripts with `--help` before using them:

```bash
.venv/bin/python scripts/recall_recovery/ingest_clinvar.py --help
.venv/bin/python scripts/recall_recovery/ingest_pubtator.py --help
.venv/bin/python scripts/recover_counts.py --help
.venv/bin/python scripts/extract_figure_variants.py --help
.venv/bin/python scripts/run_recall_suite.py --help
```

## Separately bounded live validation

A live extraction spends API quota and depends on publisher access. Choose a
small explicit PMID manifest, set the provider deliberately, and use a new
output directory. `--max-pmids` limits discovery, not a guaranteed total of model
calls. For a fixed-source reading check, use the registered evaluation harness
and its frozen-source manifest; do not substitute live acquisition midway.
Verify the external `corpus/` link before any corpus job.

Inspect the exit code and the exact run's `RUN_STATUS.json`, active database,
source ledger and trace manifest. An extraction exception after allocating a
run directory records `failed` with exit 3; a missing database records exit 4.
A nonzero exit, absent status, or unresolved source-integrity failure must not
be called a completed run. Healthy execution may still find zero variants or
leave counts unknown when the source lacks evidence.

Use [RECALL_REFRESH_RUNBOOK.md](RECALL_REFRESH_RUNBOOK.md) for measurement.
Old local canonical-baseline DBs and historical figure-reader counts are not
fresh acceptance criteria. A scorer smoke test is not a head-to-head recall
comparison, and historical locks must remain unchanged.
