# Grok 4.6 adversarial cruft review — 2026-08-20

Grok CLI model `grok-4.6` reviewed the evidence at `xhigh`, the highest
reasoning-effort value exposed by the installed CLI. The user requested
"4.6 max"; this CLI has no `max` value. The first read-only run was stopped
after Git repository discovery rejected the repository's
`extensions.relativeworktrees` setting and its fallback scan entered generated
trees. The completed second pass used the verified evidence packet with tools
disabled.

## Corrections and objections

- `.pytest_cache/README.md` is ignored local output, not a tracked file. It is
  not repository debt and must not be included in a removal patch.
- An exact-name static scan cannot prove an operator API, manual script, sibling
  application, dynamic import, or shell-out is unused. The live
  `pipeline/full_coverage.py` to `scripts/run_priority_extraction.py` shell-out
  is the concrete in-repository counterexample.
- Default-off does not mean dead. The final-check pair and count recovery are
  explicitly parked or benchmark-gated, tested, documented, and scientifically
  load-bearing.
- Zero unique commits does not make a dirty worktree disposable. The dirty
  `review-works-improvements-inefficiencies-1cad7f` tree contains 289 lines of
  uncommitted config/documentation/test work plus an untracked test.

## Recommended classification

- **REMOVE:** the obsolete hard-coded `pipeline/preprocessor.py` command-line
  tail only; the two clean Claude worktrees after a last status check.
- **ADJUST:** the incomplete script catalog, hard-coded user paths, split test
  topology, blanket deprecation-warning suppression, and the settings/constants
  split brain.
- **DEPRECATE:** unused-looking operator settings, public-looking utility
  helpers, and historical replay/retry scripts. Warn and inventory external use
  before deletion.
- **INVESTIGATE:** the dirty worktree, the two root live-network test suites,
  and undocumented `__main__` demo blocks other than the obviously obsolete
  preprocessor runner.
- **KEEP:** parked scientific stages, dynamically registered browser publisher
  strategies, the full-coverage priority driver, active manual scripts, and all
  live production bodies surrounding candidate demo tails.

Grok's smallest safe patch was documentation and local-worktree hygiene only:
catalog the eight active omitted scripts, make the gold-builder example
portable, explain or align the CI-only root tests, and prune only the clean
worktrees. It recommended no function, setting, test, or scientific-stage
deletion in the first patch.

## Post-audit execution review

A second read-only Grok CLI run used `grok-4.6`, `xhigh`, no subagents, and no
web search to challenge the exact seven-step implementation plan against the
current checkout, CI, provenance, sibling repositories, and all three
worktrees. It approved the overall shape but required these corrections before
execution:

- split offline `SupplementFile` tests before applying a module-level network
  marker;
- move the live DOI test out of `tests/unit` and the manual Wiley script out of
  pytest collection;
- migrate unique normalizer and Karger parser examples before deleting demo
  runners;
- fingerprint extraction's actual constants instead of false Settings values;
- delete only the internal one-line protein-notation wrapper immediately, while
  preserving public-looking helper names through a deprecation window; and
- keep `intern_model`, `variants_match`, `create_variant_key`, parked scientific
  stages, dynamic publisher strategies, and the subprocessed priority driver.

Those corrections were incorporated. Grok also recommended a separate branch,
but repository policy in `CLAUDE.md` explicitly requires one checkout and
`main` for this handoff, so no fourth worktree or cleanup branch was created.
