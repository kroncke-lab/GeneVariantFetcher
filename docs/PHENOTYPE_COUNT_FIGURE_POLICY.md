# Phenotype-count figure policy

Last reviewed: 2026-09-04.

This is the authority for phenotype-count evaluation figures. The two figure
families answer different questions and must remain complementary.

## Figure hierarchy

### 1. Canonical progress figure

The canonical per-protocol-change figure is the stratified
automated-versus-reference view built by
`scripts/build_stratified_phenotype_count_recovery.py`.

- It shows affected and unaffected counts for identity-matched variant-paper
  rows. A matched row with no automated count remains `NULL` in stored data but
  is evaluated and displayed as zero.
- Cohorts with different sampling frames or protocol projections share axes but
  remain in separate strata. Never pool or deduplicate them into one
  denominator, bubble, or performance claim.
- The historical cardiac cohort and all opened mixed-gold candidate arms are
  included. Candidate arms come only from append-only tranche consumption logs;
  unopened answer keys are never read.
- The opened-candidate stratum is a cumulative progress/failure map. Its
  disjoint arms may represent successive candidate protocol revisions, so it is
  not a pooled estimate of one protocol's performance.
- Cohort labels and attempt counts must be derived from the selected locked
  runs. Never hard-code a tranche count into the rendered figure.
- The figure has one descriptive title. Cohort and measure names are facet
  labels, not a second title or narrative subtitle. Keep detailed metrics,
  sources, and caveats in the CSV/JSON, changelog, and chat/report text rather
  than filling the plotting area with numbers.

Recipe and outputs:

- recipe: `scripts/build_stratified_phenotype_count_recovery.py`
- shared style/data layer: `scripts/build_phenotype_count_recovery.py`
- opened-candidate combiner:
  `scripts/build_combined_phenotype_count_recovery.py`
- outputs:
  `docs/figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.{svg,png,pdf,csv,json}`
- historical cardiac input artifacts:
  `benchmarks/codex_paper_eval/runs/20260902_false_zero_recovery_gold118/figures/data/`
- manual rebuild:
  `.venv/bin/python scripts/build_stratified_phenotype_count_recovery.py`

### 2. Companion difference figure

The companion is built by `scripts/build_gold_difference_figure.py`. It does
not replace the canonical progress figure.

- It plots `automated count - reference count` against reference count for
  every gold row with an asserted value, with separate panels for carriers,
  affected individuals, and unaffected individuals.
- Identity misses and matched abstentions are both evaluated as zero, so the
  full end-to-end failure surface is visible. Their statuses must remain
  visually distinct from supplied counts and from one another.
- For a candidate arm, the frozen baseline is drawn first on the identical gold
  rows. Paired deltas, registered bounds, and the comparison decision may then
  be shown beneath the panels.
- The candidate-to-baseline relationship is discovered from the tranche
  consumption log and the newest matching `compare*.json`; it is never inferred
  from similarly named directories.
- Per-run outputs live under
  `benchmarks/codex_paper_eval/runs/<run>/figures/gold_difference.*`. Curated
  copies and their README live in `docs/figures/gold_difference/`.

The all-row zero-imputation above is an evaluation rule only. It must not
rewrite an absent automated count to a literal stored zero.

## Shared visual language

- Both families use the palette, typography helpers, and marker-radius function
  from `scripts/build_phenotype_count_recovery.py`.
- Bubble area encodes the number of variant-paper rows at a coordinate. Every
  size key must call the same radius function as the plotted bubbles and show
  the actual pair counts it represents.
- Use direct gene/variant/PMID examples sparingly. Connector lines must land on
  the associated label block and remain clear of neighboring labels and marks.
- Use color plus fill, outline, dash, shape, ordering, or faceting so statuses
  remain legible without color.
- Editable SVG is primary; PNG and PDF are delivery copies. Every figure also
  ships with tidy CSV rows and a JSON contract/provenance record.

## Generation and validation

- `run_eval.py score` writes the per-run agreement artifact and companion
  difference figure for every scored run.
- After every scored candidate arm, `run_eval.py score` also rebuilds the
  opened-candidate agreement figure and the canonical stratified figure.
- `run_eval.py compare` re-renders the candidate companion with the frozen
  baseline, registered bounds, and decision.
- Rendering is diagnostic: a failure is recorded in the report and must not
  alter or block the scientific score.
- Renderers must verify the run lock. The companion must reproduce the locked
  report's TP/FN, supplied-count totals, and conditional MAE; when a comparison
  JSON is used, it must also reproduce its end-to-end MAE and coverage.
- Before handoff, inspect the exported PNG or SVG for truthful denominators,
  stale labels, clipping, collisions, detached connectors, and distinguishable
  marker states.

Every protocol-changing row in `docs/PROTOCOL_CHANGELOG.md` links the canonical
figure and its JSON run manifest. Link the companion too when its all-row or
paired-baseline view is material to the interpretation.
