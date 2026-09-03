# Extraction-blinded paper evaluation

This harness compares paper-reading protocols on a fixed, gold-value-blinded
paper set. `prepare` may use gold only to determine PMID eligibility under the
recorded identity or count rule. Extraction is finalized and locked before
`score` reads any gold identity, value, or row count.

## Fully traced run

```bash
.venv/bin/python benchmarks/codex_paper_eval/run_eval.py prepare \
  --seed 2026072301 \
  --paper-manifest benchmarks/evaluation_tiers/tier1_gold_50.tsv \
  --run-id my_traced_run

.venv/bin/python benchmarks/codex_paper_eval/run_eval.py extract \
  --run-dir benchmarks/codex_paper_eval/runs/my_traced_run \
  --model gpt-5.6-sol

.venv/bin/python benchmarks/codex_paper_eval/run_eval.py lock \
  --run-dir benchmarks/codex_paper_eval/runs/my_traced_run

.venv/bin/python benchmarks/codex_paper_eval/run_eval.py score \
  --run-dir benchmarks/codex_paper_eval/runs/my_traced_run
```

Every paper gets exact route/extraction call traces and explicit route/final
decision events under `llm_traces/<GENE>/<PMID>/`. The trace manifest is inside
the pre-gold lock and is revalidated during scoring. The self-contained
`llm_trace_report.html` is generated automatically and locked alongside it;
open that file directly in a browser to review the run paper by paper. See
[`../../docs/LLM_TRACING.md`](../../docs/LLM_TRACING.md) for the trace contract
and adjudication workflow.

Every eligible `score` invocation also creates a run-scoped phenotype-count
recovery figure under `figures/`. The editable SVG is required; PNG and PDF
copies are emitted when `rsvg-convert` or the workstation's built-in `sips`
converter is available. The adjacent CSV and JSON preserve the plotted
identity-matched variant–paper rows, the missing-count-as-zero evaluation rule,
the marker-size mapping, and the labeled examples. This is an evaluation
artifact only and never changes locked predictions or stored NULL values.

Runs made before trace schema v2 cannot be treated as exact request/response
audits. Their final predictions and rationales remain useful, but rerunning the
fixed manifest is required to produce authentic raw call traces.

External production projections can satisfy the same audit boundary by listing
each finalized `gvf-run` trace manifest under `production_trace_manifests` in
`predictions.json`. `lock` validates the manifest and its call/decision index,
binds its SHA-256 into `LOCK.json`, and `score` repeats the validation before
reading gold. When production stages or refolds richer run-local material after
eligibility selection, run `rebind_production_sources.py` before locking so
`selection.json` hashes the exact FULL_CONTEXT/artifact/PDF/figure inputs used
by extraction; this must not change cohort membership or consult gold values.

The registry-aware production run uses `--trust-mode trusted` and
`--identity-mode trusted`: quarantined count fields are masked and ambiguous
VariantFeatures identity classes are held out, matching the collaborator-facing
projection. Raw `all` projections remain diagnostics and must not replace the
locked primary after scoring.

Production scaffolds also pass `--paper-primary`. `papers[].variants` therefore
contains only rows whose stable origin is `llm_text`, `llm_table`,
`regex_text`, `regex_table`, or `figure`. ClinVar/PubTator-origin rows remain
inside the same hash-locked artifact under `external_linkage_variants`; their
union with paper rows is scored only as the secondary `linkage_assisted`
enrichment diagnostic. Ambiguous `mixed`, legacy, and unknown origins are
locked under `unattributed_variants` and excluded from both scored lanes. A
ClinVar citation can no longer raise the primary claim that the protocol found
a variant in a paper.

The complete cross-gene regression plan is the seeded 29-tranche
[`../evaluation_tiers/mixed_gold/`](../evaluation_tiers/mixed_gold/README.md)
suite. Its manifests use named-variant rather than complete-count eligibility,
so papers with valid identity gold are not lost merely because phenotype counts
are absent. Gold provenance remains stratified in every report. Protocol
changes should be compared as a paired baseline/candidate run on one tranche;
once its score is inspected, move confirmation to the next unopened tranche.
The registry preregisters the delta thresholds and PMID-clustered confidence
rule; apply it with `run_eval.py compare` rather than choosing a pass criterion
after opening the reports.

## Current production gold_120 test

Use the registry-aware setup command for a fresh test of the production pipeline:

```bash
.venv/bin/python benchmarks/codex_paper_eval/setup_production_eval.py create \
  --run-id 20260821_current_changes_gold119 \
  --email brett.kroncke@gmail.com
```

Despite its historical name, `gold_120` is now a pinned 118-attempt / 114-unique-
PMID cohort: KCNH2 has 28 attempts after wrong-paper PMIDs 10086972 and 14642689
were removed; KCNQ1, RYR2, and SCN5A have 30 each. The setup refuses a manifest
digest, count, seed, source-coverage, or code-fingerprint mismatch. It creates per-gene PMID
files plus two executable phases inside the run directory:

1. `run_extraction.sh` runs the current production `gvf-run` route without gold,
   source recovery, corpus mutation, full-coverage add-ons, or publication.
   It passes `--gold-free-run`, which disables `gvf-run`'s otherwise automatic
   per-layer gold discovery/scoring and all file-backed alias maps whose
   provenance may include benchmark gold until the prediction lock is complete.
2. `lock_and_score.sh` rebinds the exact run-local sources, projects the final
   production databases through the maintained `db_to_predictions.py`, locks
   predictions and trace manifests, and only then opens gold for scoring.

The generated `RUNBOOK.md` records the exact absolute commands and rationale.
Create a new scaffold instead of bypassing its source-fingerprint check if code
changes before extraction. The setup refuses multiple completed per-gene runs,
backup DBs, missing or zero-call trace manifests, and any DB other than
`RUN_STATUS.active_db`.

## Standard comparison set: 50 papers (48 cardiac + 2 collaborator-reviewed BRCA2)

**As of 2026-08-11 the active set for comparing paper-processing strategies is
[`../evaluation_tiers/tier1_gold_50.tsv`](../evaluation_tiers/tier1_gold_50.tsv):**
the fixed cardiac 48 plus the two BRCA2 papers with lead-approved Variant
Browser adjudications by Nate (26833046 and 26848529). The dated
`highcarrier48_plus_brca2_collaborator2_20260811.tsv` is an exact source mirror,
pinned by an integrity test. The active BRCA2-only arm is
`brca2_2_collaborator_reviewed_20260811.tsv`.

The dated `highcarrier48_plus_brca2_20260810.tsv` and
`brca2_8_papers_20260810.tsv` files remain frozen historical artifacts. Their
six internally derived BRCA2 papers are not part of the active scored benchmark
or new strategy comparisons. Do not edit or relabel historical runs scored on
them. The live Variant Browser queue is a separate 46-paper operational review
set: Nate's two gold papers plus 44 additional papers awaiting or supporting
review.

Frozen production baselines (comparison evidence, not current defaults):

- Cardiac arm: `runs/20260726_fixed48_production` (full gvf-run pipeline,
  all layers) and `runs/20260726_fixed48_production_paperonly` (paper-derived
  layers only). Single-model references: `runs/20260724_fixed48_sol`,
  `runs/20260724_fixed48_grok`.
- BRCA2 reference extraction: the historical
  `runs/20260810_brca2_8_production` and `_paperonly` sibling. Active
  two-paper metrics are a PMID-filtered projection of those locked predictions;
  rerun the two-paper manifest for a new protocol comparison.

Gold provenance differs by arm, and the harness records it per gene in
`selection.json` and `report.json` under `gold_sources`:

- Cardiac genes score against the manual, human-curated gold standard
  (`gene_variant_fetcher_gold_standard/normalized/`).
- BRCA2 has no manual gold. The active arm selects only the two rows sets in
  `benchmarks/curated_extraction_eval/gold_overrides/BRCA2_recall_input.csv`
  that are traceable to lead-approved Variant Browser adjudications. Report
  BRCA2 separately and never fold it into cardiac headline metrics.

Known scope limitations remain: 26833046 has family-versus-carrier count
semantics and 26848529 is a reviewed-positive subset rather than exhaustive
paper-level gold. Seeded random sampling (`--per-gene` without a
manifest) remains cardiac-only: BRCA2 has too few gold papers to sample and
widening the pool would change historical seeded draws.
