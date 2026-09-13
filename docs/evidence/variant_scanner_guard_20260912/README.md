# Variant scanner resource guard — 2026-09-12

Prepared on September 12 while BRCA2 shard 3 was still live, using a test-only
import hook and an ignored patch without changing production source. Applied to
`main` on September 13 after Brett authorized landing, the process inventory
showed no grant extraction/recovery jobs, and the run log recorded BRCA2
extraction complete at 10:18 on September 12. The freeze tag and existing run
artifacts are retained. No second checkout or branch was used.

## Behavior

- `SCANNER_MAX_CHARS` defaults to 2,000,000 raw characters. The scanner checks
  length before normalizing, copying the source into IPC, or matching. Oversize
  input is skipped in full; a prefix is never presented as a complete scan.
- `SCANNER_BUDGET_SECONDS` defaults to 30 seconds for the worker. Both settings
  must be positive; the time budget must also be finite. They cannot be disabled.
- Each accepted-size scan runs in a fresh, killable interpreter using the project
  Python. `subprocess.run` kills and reaps the child on timeout. The budget covers
  worker initialization, normalization, matching, attribution and IPC; bounded
  parent serialization/deserialization and audit writes add overhead. No global
  signal handler, uncancellable thread, persistent pool, or new dependency.
- The private worker exposes the existing `utils` package directory without
  executing its eager initializer, which otherwise loads LLM/network clients.
  This does not modify the parent interpreter. Scanner submodules and packaged
  gene/reference data still load normally; source text is only JSON data.
- `ScanResult.status` distinguishes `complete`, `skipped`, `timed_out`, and
  `failed`. Incomplete scans return no candidates/hints. Status, reason, source,
  gene, raw length, limits and elapsed time appear in `variant_scan` extraction
  metadata and decision traces. The paper census also records scanner status.
- Extraction/replay immediately writes incomplete-scan records under
  `pmc_fulltext/variant_scan_audit/*.json`, before the LLM call. Unique filenames
  and atomic replacement support parallel attempts and preserve prior records.
  These records survive LLM failure and a replay regression gate retaining the
  original extraction. Audit I/O errors are logged and surfaced in scan stats.
- The cap protects the pre-scan. Earlier table parsing, later scanner-result
  merging and other extraction stages do not acquire a document timeout from
  this change. Skipping scanner hints is a degraded extraction, not a finding
  that the paper contains no variants.

## Regex changes

Possessive suffix digit runs repair overlapping repetitions in the full,
parenthesized and three-letter protein patterns. Capture groups, match spans and
accepted ordinary tokens are preserved. Possessive whitespace runs also fix
`DELTA_RE`'s two whitespace repetitions separated by an optional hyphen.

The merge helper uses `islice(finditer(...), 100)` instead of materializing all
matches before slicing, retaining the same first 100 mentions.

Before the repair, `p.Arg123fs` + `1` repeated N times + `_` took approximately
0.007 / 0.11 / 1.90 seconds at N=1,000 / 4,000 / 16,000. After repair the three
protein patterns reject N=100,000 in approximately 0.0006–0.0008 seconds each on
this workstation. This demonstrates the synthetic defect; it does not establish
which regex caused the original production stalls.

## Independent advice and dispositions

Brett requested both Grok and Agy CLI reviews and explicitly approved sending
the prepared design prompt containing selected code excerpts and incident
context. Both reviews were read-only; their output is retained in ignored
`tmp/variant_scan_guard/`.

Both recommended a pre-normalization size ceiling, explicit incomplete status,
regex repair, immediate durable skips, and synthetic/offline tests. Both favored
an in-process cooperative timer as a smaller initial change. The implementation
uses a hard worker deadline because a cooperative timer cannot interrupt an
unknown future pathological regex. Their warnings about capture semantics,
threaded callers, IPC fields, child cleanup, and failure-path persistence inform
the regression coverage. Agy's illustrative pattern regrouping was not copied:
it would change capture groups used by the existing normalization code.

Worker overhead measured about 0.06–0.09 seconds on short inputs with offline
built-in gene metadata. A fresh scan with the workstation's explicit
VariantFeatures database took about eight seconds: per-process metadata cache
initialization is a real cost, and the configured budget includes it. No paid
extraction evaluation was performed to calibrate a recall/runtime tradeoff.

## Validation

The dedicated resource-limit tests cover oversize input before normalization or
worker creation, inclusive boundaries, invalid/environment-backed settings,
actual stuck-regex termination from an extraction thread, child reaping, worker
failure, continuation with the next document, audit I/O errors, IPC fidelity,
matching under an external watchdog, differential spans/captures, and audit
persistence across LLM failure and later metadata replacement.

Two scanner cache tests invoke the internal scan engine so their local
monkeypatches still inspect attribution calls. Three parser tests now explicitly
disable their unrelated adjudicator to remain offline. Parser and router fakes
implement the updated scanner metadata interface.

The first applied-source run passed 3,048 tests and exposed one unrelated stale
repository-freshness expectation: current `TASKS.md` already has the committed
"Grant penetrance interpretation and next analysis gate" heading. The expected
heading list now includes it; no checklist content was changed by this fix.

The prepared-source focused suite passed 313 tests; the full offline suite
passed **3,048 tests in 124.72 seconds**, including the delta-whitespace repair
and final metadata-placement regression. Prepared-source Ruff checks and
formatting checks passed. The patch applied to September 13 main with only a
changelog context adjustment around a newer entry.

The final applied-source full offline suite passed **3,049 tests in 121.91
seconds** on September 13, with Ruff lint and formatting checks passing for all
nine changed Python files. The complete log is retained locally in
`tmp/variant_scan_guard/pytest_applied_final_20260913.txt`.

No giant source files are added to Git. This fix does not modify grant-run output,
replay checkpoints, live metrics, scored benchmarks or phenotype figures.
Scientific recall/MAE impact remains unmeasured.

Reproduce validation:

```bash
GVF_DISABLE_LOCAL_DATA=1 LITELLM_LOCAL_MODEL_COST_MAP=True .venv/bin/python -m pytest tests/unit -q
```

The September 12 scratch harness ran that suite with `-p overlay_plugin` and
an explicit `PYTHONPATH` pointing at the test-only import hook. The September 13
validation runs directly against the applied checkout without that hook.
