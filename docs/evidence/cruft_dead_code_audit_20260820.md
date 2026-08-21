# Cruft, stale-code, and deprecation audit — 2026-08-20

## Decision

This audit found real cleanup work, but not a license for a broad dead-code
deletion. GeneVariantFetcher has manual operator entry points, installed-package
compatibility surfaces, dynamic publisher plugins, historical benchmark tools,
and intentionally parked scientific stages that ordinary static analysis does
not understand. Candidates below are therefore labeled by evidence and removal
risk. The evidence-backed queue was subsequently executed; the concrete result
is recorded below.

The documentation-only first patch completed two adjustments: the script
catalog now includes all eight previously omitted top-level scripts, and the
gold-snapshot regeneration example no longer embeds one user's absolute path.
The completed cleanup was removed from the active forward-only `TASKS.md` list.
Grok's independent challenge is recorded in
[`grok46_cruft_review_20260820.md`](grok46_cruft_review_20260820.md).

## Execution result

- Removed all three zero-unique-commit Claude worktrees and their fully merged
  local branches. The two review trees were clean. Before force-removing the
  dirty `review-works-improvements-inefficiencies-1cad7f` tree, its nine paths
  were reconciled individually:
  - `.env.example`, `CHANGELOG.md`, `docs/ARCHITECTURE.md`,
    `docs/QUICKSTART.md`, and `tests/unit/test_env_documentation.py` were ported
    to `main` in updated form.
  - its `TASKS.md`, `docs/PROTOCOL_CHANGELOG.md`,
    `docs/PROTOCOL_COST_EVAL.md`, and
    `tests/unit/test_repository_freshness.py` changes were rejected because they
    were stale or would regress the current single-checklist/freshness contract.
- Removed obsolete embedded command/demo tails from `pipeline/preprocessor.py`,
  `utils/variant_scanner.py`, `utils/variant_normalizer.py`,
  `harvesting/supplement_reference_parser.py`, and
  `gene_literature/pubmind_fetcher.py`. Unique normalizer and Karger supplement
  cases were migrated into unit tests first. That migration found and fixed a
  real boundary bug where the `2` in `KCNH2` was parsed as a variant count.
- Moved offline root tests into `tests/unit`, live tests into
  `tests/integration`, split four offline `SupplementFile` tests away from the
  network-marked suite, moved the live DOI check out of the unit tree, and moved
  the manual Wiley probe to `scripts/check_wiley_api.py`. Its unused inline
  fallback client was deleted.
- Replaced blanket warning suppression with fail-on-project-deprecation plus
  three exact SWIG dependency ignores identified by a full warning inventory.
- Made extraction constants the provenance authority. Provenance schema v3 now
  fingerprints `TEXT_TRUNCATION_MAX_CHARS`,
  `SCANNER_MERGE_MIN_CONFIDENCE`, and `SCANNER_MAX_HINTS` rather than three
  operator settings runtime never read. False settings log a warning when
  explicitly configured; public-looking legacy helpers emit standard
  deprecation warnings for one compatibility window.
- Deleted the proven internal no-caller
  `utils.protein_notation.is_valid_protein_notation`; retained live
  `variants_match`, `create_variant_key`, `intern_model`, dynamic publisher
  plugins, parked scientific gates, recovery code, and reproducibility assets.

## Verification after execution

- Full offline unit suite with project deprecations promoted to errors:
  **2,164 passed** after the follow-up source/readiness/publish hardening.
- Bounded offline end-to-end harness: **6 passed**.
- Curated negative precision guards: **3 passed**.
- Integration default behavior: **4 passed / 43 network tests skipped**; the
  network tests remain explicit opt-in coverage rather than silently collected
  unit tests.
- Repository-wide Ruff: all checks passed; all **429** Python files formatted.
  `git diff --check`, the moved Wiley probe `--help`, the preprocessor import
  smoke, environment-documentation gate, repository-freshness gate, settings
  tests, provenance tests, and parser regressions all pass.
- `git worktree list` contains only the authoritative `main` checkout; no
  `claude/*` branches remain.
- The immediately preceding no-inference/source-scope run is the locked
  **119-paper** gold-value-blind acceptance test: count-eligible variant recall
  **548/633 (86.57%)**, counted-gold PMID precision **97.51%**, and carrier MAE
  **0.425**. Cleanup did not change thresholds or model routing.
- Post-execution cohort/source/trace membership is exact for BRCA1 **50/50**,
  BRCA2 **45/45** after removing the canine PMID, and BMPR2 **50/50**. A later
  fail-closed audit correctly marks all three raw DBs **FAIL** because 346 / 129
  / 56 ambiguous live identities remain outside the allowed novel/cDNA-only
  classes. The Variant Browser trusted projection holds those rows out, but
  public publication remains **HOLD** pending the new cardiac run,
  source-grounded adjudication, and staging/UI verification.

## Labels

| Label | Meaning |
| --- | --- |
| **REMOVE** | Proven local or internal residue with a bounded deletion target. Recheck immediately before deleting. |
| **ADJUST** | Live or potentially live surface whose placement, documentation, naming, or wiring is stale. |
| **DEPRECATE** | No current in-repository caller, but external/operator compatibility is plausible. Warn and inventory before removal. |
| **INVESTIGATE** | Evidence is insufficient or hidden invocation/reproducibility risk remains. |
| **KEEP** | A tempting static false positive that is demonstrably load-bearing or deliberately retained. |

## Original prioritized queue (resolved as recorded above)

| Priority | Label | Target | Evidence | Required action / proof |
| --- | --- | --- | --- | --- |
| P0 | **INVESTIGATE** | `.claude/worktrees/review-works-improvements-inefficiencies-1cad7f` | Branch has zero commits ahead of `main`, but the worktree has 8 modified files (`+289/-2`) plus untracked `tests/unit/test_env_documentation.py`. | Compare the environment-documentation test and useful config/docs changes against current `main`; port only still-valid pieces, then make the tree clean before pruning. |
| P0 | **ADJUST** | Root test topology | `pytest.ini` defaults to `tests/unit`. CI separately names four offline root files; two live-network root suites are outside the documented `tests/integration` path and outside default/CI selection. | Move the four offline files into `tests/unit` and the two network suites into `tests/integration`, apply `requires_network`, update CI, then prove collection and results are unchanged. Do not delete tests to make the tree look tidy. |
| P0 | **ADJUST** | `pytest.ini` warning filters | All `DeprecationWarning` and `PendingDeprecationWarning` messages are globally ignored, hiding the signal needed to retire stale APIs. | Run the suite with warnings restored, inventory first- versus third-party warnings, then replace blanket ignores with narrow dependency-specific filters or fail on project-owned deprecations. |
| P0 | **ADJUST** | `config/settings.py` versus `config/constants.py` | `extraction_max_chars`, `scanner_merge_confidence`, and `scanner_max_hints` are operator-facing settings, but extraction reads `TEXT_TRUNCATION_MAX_CHARS`, `SCANNER_MERGE_MIN_CONFIDENCE`, and `SCANNER_MAX_HINTS` constants instead. | Choose one authority. Wiring settings into extraction is protocol-affecting and requires curated/gold-119 measurement; retiring the knobs requires a deprecation warning first. |
| P1 | **REMOVE** | Two clean Claude worktrees | `code-review-max-e46bda` and `gold-standard-precision-ad637f` are clean; their branches have zero commits ahead of `main`, while `main` is 9 commits ahead. | Last `git status` check, then prune the two exact worktrees and their fully merged branches. Do not include the dirty third tree. |
| P1 | **REMOVE** | `pipeline/preprocessor.py` lines 271-end | The importable `PaperPreprocessor` is live, but its undocumented `__main__` runner defaults to three obsolete `/mnt/temp2/kronckbm/...` paths and has no repository invocation. | Delete the command-line tail only; retain the class. Run preprocessing tests and an import smoke. |
| P1 | **ADJUST — DONE** | `scripts/README.md` | Its cleanup rule requires every script in exactly one section, but eight active scripts were absent. | Catalogued `build_status_dashboard.py`, `collaborator_readiness_audit.py`, `ezproxy_relogin.py`, `fetch_linked_supplements.py`, `final_check_triage_report.py`, `phenotype_value_precision.py`, `publish_curated_review_set.py`, and `smoke_azure_models.py`. |
| P1 | **ADJUST — DONE** | `gene_variant_fetcher_gold_standard/README.md` | Baseline regeneration embedded `/Users/kronckbm/...` paths. | Replaced them with explicit `GVF_REPO` and `VARIANT_BROWSER_REPO` placeholders. |
| P1 | **DEPRECATE** | No-reader settings | Exact repository search finds no runtime reader for `gemini_api_key`, `azure_deployment_gpt5_codex`, `tier1_model`, `enable_tier3`, or `tier1_use_llm`; some remain in provenance snapshots. `intern_model` is a legacy fallback and can still be reached by an empty Tier-2 selection, so it is not safe for immediate deletion. | Emit a startup warning only when a deprecated variable is set; document replacement/no-op behavior; keep for at least one release; then remove fields and update provenance/schema tests. |
| P1 | **DEPRECATE** | Uncalled utility helpers | No exact repository call/import for `utils.html_utils.extract_dois_from_html`, `extract_pmids_from_json_results`; `utils.pubmed_utils.query_pubmed_for_gene`, `validate_pmid`; `utils.protein_notation.is_valid_protein_notation`; or `utils.variant_normalizer.normalize_protein_variant`, `normalize_cdna_variant`, `find_matching_variants`. | Check sibling Variant Browser, notebooks, and installed-use expectations; add deprecation warnings or move to an internal namespace; delete only after a release window. The lack of a local caller alone is not proof. |
| P1 | **INVESTIGATE** | Embedded demonstration runners | Undocumented `__main__` tails exist in `utils/variant_normalizer.py`, `utils/variant_scanner.py`, `harvesting/supplement_reference_parser.py`, and `gene_literature/pubmind_fetcher.py`. The scanner fixture is explicitly copied into unit tests; much of the normalizer print suite is also covered. | Confirm no runbook/manual use. Move any unique assertions into unit tests, then strip demo tails without touching importable production code. |
| P2 | **DEPRECATE** | `scripts/replay_cap_trip_extractions.py`, `scripts/retry_failed_extractions.py` | No callers, tests, or current docs beyond their historical/manual catalog rows. | Add a dated deprecation note and identify the locked artifacts they reproduce. Archive or delete only if the current benchmark harness supersedes the workflows exactly. |
| P2 | **ADJUST** | Historical/manual script status | `discover_recall.py`, `enrich_from_variantfeatures.py`, `targeted_land.py`, and credential probes are intentionally manual/historical but lack a retirement owner/date. | Add `retained_for`, replacement, last-verified date, and removal condition to their catalog rows. Do not infer deadness from a low reference count. |

## Static false positives: keep

| Target | Why it stays |
| --- | --- |
| `pipeline/paper_final_check.py`, `paper_final_check_gate.py`, `paper_final_check_triage.py` | Explicitly parked as a pair because measured cost/latency did not justify default-on use. Tests, reason codes, prompts, and re-enable contract are retained deliberately. |
| `pipeline/count_recovery.py` and `scripts/recover_counts.py` | Default-off pending a clean carrier-only benchmark. It is an experimental, guarded recovery path—not dead code. |
| `harvesting/browser_html/strategies/*.py` | `pkgutil.iter_modules` + `importlib.import_module` dynamically discovers the modules. Registry tests pin required strategies. |
| `scripts/run_priority_extraction.py` | `pipeline/full_coverage.py` shells out to this path; deleting it breaks the full-coverage walk even though an import graph misses the edge. |
| The eight scripts added to `scripts/README.md` | The defect was the catalog, not the scripts. They have current pipeline, runbook, audit, or staging roles. |
| `.pytest_cache/` and Python bytecode caches | They are ignored local tool output, not tracked repository cruft. Do not create a source-code deletion patch for them. |
| Locked benchmark runs and append-only history | Reproducibility evidence is immutable by repository policy; supersede with new evidence or an erratum. |

## Safe execution order

1. Reconcile the dirty worktree; prune only the two clean, fully merged trees.
2. Normalize the test layout and warning policy without changing assertions or
   network opt-in behavior.
3. Remove the obsolete preprocessor runner and migrate proven demo assertions.
4. Add conditional deprecation warnings for settings/utilities/scripts and
   observe one release window.
5. Remove deprecated internals only after sibling/notebook/API checks.
6. Treat any setting-wiring or parked-stage change as a protocol change: run the
   curated extraction benchmark and locked gold cohort before acceptance.

## Audit method and limits

The pass used tracked-file inventory, AST parsing of 224 runtime/script Python
files, import/reference searches, CLI/CI/test discovery review, script catalog
comparison, absolute-path and legacy-marker searches, branch/worktree ancestry
and dirtiness checks, and Grok 4.6 `xhigh` adversarial review. `vulture` was not
installed; no dependency was installed merely to expand the scan. Static
analysis cannot see external notebooks, sibling repositories, operator shell
history, reflection, plugin discovery, or subprocess path edges, so every
deletion candidate still carries a required proof step.
