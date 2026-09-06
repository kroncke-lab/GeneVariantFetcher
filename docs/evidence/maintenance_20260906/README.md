# Runtime and documentation maintenance — 2026-09-06

The maintenance fixes extraction-summary crashes, incorrect run-directory
selection and successful exit codes after a failed workflow. It restores the
local installed launcher, makes calibration checks portable to fresh CI, and
updates current setup/testing and deprecation guidance. No paper extraction,
gold scoring, source acquisition or headline promotion was performed.

## Reproduced problems and repairs

- `pipeline/steps.py` summed model-authored `total_variants_found`. A string in
  the archived tranche-03 BRCA1 attempt raised `TypeError`. Summary statistics
  now count actual variant-object rows. A network-blocked replay on temporary
  copies succeeds with the original two rows and no model call; all 32 archive
  hashes remain unchanged. See [replay receipt](archived_failure_replay.json).
  Cached variants and metadata are preserved. This count is before aggregation,
  not a count of globally unique variants.
- `gvf-run` used the newest sibling directory after extraction, which could
  select another run or miss an explicit resume outside the output tree. It now
  allocates and passes the exact directory, avoids same-second collisions and
  does not leak a resume environment override. Workflow exceptions and explicit
  failure results write failed status/exit 3 there; no database writes exit 4.
  Failures before allocation have no run-local status. Regression fixtures
  include a newer, successful older run whose bytes must remain unchanged.
- `gvf extract` and `python -m cli.automated_workflow` ignored failure results
  returned as dictionaries. Both now exit nonzero; success/failure CLI fixtures
  cover both entry points.
- The local uv environment had no `gvf` executable despite module tests passing.
  Reinstalling the editable package repairs it. The declared `cli:app` entry
  point and `cli/__main__.py` were already correct. CI now tests console and
  module help inside the isolated wheel environment as well as packaged data.
- Starting-main CI failed two calibration tests because full operator traces
  under ignored `results/` were absent. An exporter validates the original
  hashes and totals and emits 1,809 usage-only receipts for APOE, BRCA1, BRCA2
  and MYBPC3. CI checks the pinned receipt bytes, historical cost profile and
  tracked cardiac predictions. No prompts, responses or credentials are in
  [the fixture](../../../tests/fixtures/cost_calibration_usage.json); no locked
  calibration, registry, prediction or score was rewritten.

## Deprecated code and current documentation

The audit inventoried nine warned public helpers. Their production caller
search found no reason to remove compatibility APIs. HTML adapters own their
source-specific parsers; broad old HTML helpers have no single public drop-in
replacement. The PubMed replacement accepts a full query and returns a list,
where the deprecated wrapper returns a set. Normalization and matching wrappers
also have different interfaces. Docstrings now say so; warning behavior and
return behavior remain compatible. Python 3.11+ can import `Annotated` from
the standard library, so the redundant `typing_extensions` import is removed.

Current documentation now verifies both installed and module entry points,
supports uv environments without pip, labels illustrative timings, explains
that `--no-source-recovery` does not restrict ordinary harvesting to PMC, and
separates installation checks from paid accuracy validation. Stale fixed test
counts, June-only acceptance examples and unsupported blanket browser-access
claims and unsupported publisher-coverage percentages were removed. The PMC
section was checked against its current FAQ and Open Access Subset documentation
(linked from the API guide). Azure setup is explicit for this workstation's preferred
allocation; provider defaults and routing did not change. Historical reports
remain dated evidence.

## Validation and adversarial review

Baseline: 474 Python files compiled; 2,890 unit tests and nine bounded/negative
tests passed locally. Starting-main GitHub CI failed its two workstation-trace
dependencies. The updated source compiles 476 Python files. The full local
suite passed 2,906 tests; a later exporter guard and regression passed in the
42-test focused run. Nine bounded/negative tests and eight current-doc/config
checks also passed. Ruff, formatting and dependency compatibility pass. The
wheel imports all 151 packaged modules outside the checkout; packaged reference
data and all three CLI help checks pass. Thirteen documented command/script
help checks pass. See [verification receipts](verification.json) and
[wheel checks](wheel_smoke.json).

Remote CI runs the complete committed test inventory after push. These are
software-health checks, not a new estimate of extraction accuracy.

Grok, Claude Sonnet and Agy/Gemini 3.1 Pro reviewed bounded source packets, then
received the final diff for adversarial review. Prompts, raw JSON and stderr
are preserved in this directory. Initial findings about the metadata crash,
missing launcher and misleading replacement guidance were adopted. Claims
that `cli:app` was missing, module invocation unsupported, or Azure lacked its
model prefix were rejected after inspecting executable code. Counting actual
rows was selected over coercing model totals, which could retain inflated
statistics. Agy's final non-object JSON warning was adopted in the exporter;
Claude's collision-coverage suggestion produced the same-second isolation test.
Grok's final review added explicit refusal of malformed nested call records;
output-token arithmetic retains the historical calibration convention. The
real export remains byte-identical.
Warnings about missing imports and escaping CLI errors were rejected against
passing executable tests and the existing catch blocks. Detailed dispositions
and costs are in [the review summary](review_summary.json). Review prompt copies
have only trailing whitespace/final-newline normalization from pre-commit.

The two Claude calls report **$0.40454** combined. Grok's completed reviews report **$0.03805** combined; the first final-review
request timed out and was retried with a smaller packet. Timeout billing and
Agy dollar billing are unavailable
and are not counted as zero. No extraction API spend was incurred by maintenance.
The prior improvement campaign's extraction ledger is unchanged.
