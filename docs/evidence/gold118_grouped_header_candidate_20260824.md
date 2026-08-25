# Gold-118 grouped-header and structural candidate — 2026-08-24

## Disposition

The implementation is source-backed and regression-tested, but the metric
projection is a diagnostic rather than a new blind lock. The immutable accepted
result remains `20260824_postfix_gold118`. No accepted extraction JSON, SQLite
database, prediction lock, gold snapshot, or public review surface was mutated.

## General repairs

1. A clinical table with one parent count header followed by exactly two
   opposing phenotype subheaders now maps the children to affected/unaffected
   and sets total carriers to their sum. Footnote markers remain accepted as
   annotations, not digits. Existing attribution, gene, identity, and current-
   study gates still apply.
2. Deterministic markdown rows are available as guarded merge evidence for
   moderate tables, not only very large compact-mode tables. Direct source
   values can therefore correct LLM or verifier drift without inventing rows.
3. SQLite migration preserves strict structural cDNA breakpoints, including
   repeated `c.` endpoints and explicit deleted lengths. Trusted prediction
   export retains structural-only rows only for structural variant classes with
   a canonical structural identity. Distinct breakpoints never collapse merely
   because both normalize to an exon-level alias.

## Live source replay

SCN5A PMID 29709101 contains the flattened header `Number of patients` followed
by `BrS+ | BrS−`. The staged replay at
`benchmarks/codex_paper_eval/runs/20260824_postfix_gold118/production_runs/SCN5A/20260824_075355/refresh_20260824_183305/`
completed without promotion and rebuilt
`SCN5A.refresh_20260824_183305.db`. The repaired rows include:

| Variant | Total | BrS+ | BrS− |
| --- | ---: | ---: | ---: |
| c.1127G>A / p.Arg376His | 12 | 6 | 6 |
| c.2268_2271del / p.Phe756Leufs*8 | 4 | 1 | 3 |
| c.2466G>A / p.Trp822* | 6 | 4 | 2 |
| c.2658T>A / p.His886Gln | 7 | 4 | 3 |
| c.4813+3_4813+6dup | 24 | 13 | 11 |

All twelve count-bearing table variants in the rebuilt DB agree with their
source cells. Grok's adversarial review then narrowed the production change:
post-adjudication repair is existing-only, requires all three nonnegative counts
plus `total = affected + unaffected` and phenotype-labelled provenance, cannot
re-add a filtered identity, and records `deterministic_count_repair.source`
without replacing the model/verifier decision.

## Diagnostic metrics

| Measure | Immutable lock | Candidate projection |
| --- | ---: | ---: |
| TP / FP / FN | 546 / 284 / 86 | 548 / 284 / 84 |
| Recall | 86.39% | 86.71% |
| Raw precision | 65.78% | 65.87% |
| F1 | 74.69% | 74.86% |
| Counted-extra precision | 97.50% | 97.51% |
| Count-bearing-only precision | 93.69% | 93.75% |
| Carrier supplied; absolute error; MAE | 206; 68; 0.330 | 208; 47; 0.226 |
| Affected supplied; absolute error; MAE | 49; 27; 0.551 | 59; 27; 0.458 |
| Unaffected supplied; absolute error; MAE | 18; 9; 0.500 | 28; 9; 0.321 |

The identity gain shown here is the conservative projection of structural-only
rows already present in the current DB. A fresh migration distinguishes an
additional exact structural breakpoint, but it is excluded from this table
until a complete gold-free run is locked and scored.

The diagnostic is machine-reproducible, not a hand-edited table. Its inputs and
summary are:

- `benchmarks/codex_paper_eval/runs/20260824_postfix_gold118/diagnostics/source_backed_overlay_20260824.json`
- `benchmarks/codex_paper_eval/runs/20260824_postfix_gold118/diagnostics/source_backed_report_20260824.json`
- `benchmarks/codex_paper_eval/apply_source_backed_overlay.py`
- `benchmarks/codex_paper_eval/score_source_backed_overlay.py`

Reproduce the compact report with:

```bash
.venv/bin/python benchmarks/codex_paper_eval/score_source_backed_overlay.py \
  --selection benchmarks/codex_paper_eval/runs/20260824_postfix_gold118/selection.json \
  --predictions benchmarks/codex_paper_eval/runs/20260824_postfix_gold118/predictions.json \
  --overlay benchmarks/codex_paper_eval/runs/20260824_postfix_gold118/diagnostics/source_backed_overlay_20260824.json \
  --gold-root gene_variant_fetcher_gold_standard/normalized \
  --summary-only --out /tmp/source_backed_report.json
```

This command is deliberately unblinded and the report says so. It verifies the
immutable baseline prediction digest before applying the source overlay, then
records the selection, overlay, and gold-file digests. It is evidence for the
candidate effect, never a replacement for the next blind lock.

## Answer-key findings

The two dominant KCNH2 affected-count errors are not safe extraction masks.
PMID 16029385 supports the extracted 16-person clinical group while the current
fixture contains 9; PMID 22338672 supports seven K897T torsades carriers while
the fixture contains zero. They remain separately versioned gold-review items.
Changing extraction to mimic those values would improve MAE while degrading
source fidelity.

## Cost and next gate

The locked trace for this paper is a $0.456 public-list-price proxy. Three
targeted replays are $1.37 at normal-effort equivalence, but two used higher
reasoning and the refresh utility failed to persist usage traces. The experiment
therefore charges a conservative $10 reserve against the user's $100 ceiling
and makes no further paid calls. Promotion requires: add refresh tracing; freeze
code and sources; run gold-free; hash predictions and traces; then reveal gold
and require recall/precision/count-coverage non-regression plus lower MAE.

## Independent review

Grok 4.6 ran at its maximum `xhigh` effort on the live diff. Its first verdict
was HOLD: it found overly broad post-adjudication replacement, permissive
structural trust, an exact-vs-generic merge hole, invalid-cDNA promotion order,
moderate-table additive blast radius, and a non-reproducible metric claim. The
implementation, negative tests, and machine-scored overlay above address those
findings. Grok's targeted final re-review returned no unresolved P0-P3 findings
and `MERGE`.

AGY ran Gemini 3.1 Pro at its maximum supported `high` effort, but its first
request timed out after 20 minutes and the saved-session continuation returned
only `The prior inspection did not complete.` It supplies no substantive
approval and is not represented as one.

## Post-lock source-integrity erratum

The acquisition probe used during this work re-folded supplements into four
run-local sources: RYR2 19216760, 22222782, and 34661651, plus SCN5A 29544605.
The ordinary lock command now rejects those paths because their bytes no longer
match `selection.json`. The accepted `LOCK.json`, predictions, and report were
not changed. The diagnostic scorer therefore binds to the immutable prediction
digest and explicit overlay rather than claiming the current source directory
is still selection-valid. Restore the original source snapshots or use a clean
scaffold before a fresh blind run.
