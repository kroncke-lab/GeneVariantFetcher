# Autonomy-at-Scale Roadmap

**Goal.** Run GVF unattended across hundreds of genes and 10^4–10^5 papers, with
**automated quality gates as the primary control**. Human adjudication is an
**exception-only escape hatch** for the rare marginal case — never a per-paper or
per-run step. Must generalize across gene classes (cardiac missense →
BRCA truncating/case-control), without per-gene hand-tuning.

This doc is the forward plan for that goal. Live recall numbers stay in
`docs/RECALL_STATUS.md`; the recall-gap levers stay in `TASKS.md`.

## Landed foundation (PR #140, `trust-validation-hardening`)

Trust/validation hardening + the "fleet-honesty" pass:

- Honest self-measurement: end-to-end count error (misses as zero; median/p95/max);
  adjudicated-precision **calibration** sampler (Wilson CI); adjudication overlay
  consumed by the scorer (opt-in escape hatch).
- Reproducibility: run provenance (git SHA, prompt/extractor + dependency hashes,
  resolved model routing) in `run_manifest.json`; `requirements.lock`.
- Operational trust: per-file `SAVEPOINT` rollback in migration.
- **Fleet honesty**: `gvf-run` returns a non-zero exit (`EXIT_STAGE_WARNINGS`) +
  writes `RUN_STATUS.json` when a completeness stage failed; best-effort stages
  recorded in a separate `stage_warnings` channel; fail-closed regression gate
  (missing baseline / vanished metric = failure) now also gates the e2e count
  error; cold-start harness strips warm-start `GVF_*` env, refuses a non-empty
  isolated corpus, and derives the covered-gene set at runtime.

The foundation **measures and reports** honestly. The per-fact **decision** layer
(the trust gate) landed next — v1 is now merged (see below).

## Keystone — per-fact confidence/trust gate → two-tier DB (v1 LANDED, PR #142)

The piece that actually removes the human from the margin. Sorts every extracted
fact into **trusted** vs **quarantine** using gold-free checks, so downstream
products consume only the trusted tier and the held tier feeds audit/calibration.

**Status — v1 MERGED to main (PR #142):** the two-tier
schema (`penetrance_data.trust_tier / trust_reasons / trust_rule_version`), the
gold-free rule core (`pipeline/trust_gate.py`: `arith_inconsistent`,
`count_is_total`, `population_count`, `paper_outlier`; pure `evaluate_fact` +
soft-quarantine `apply_trust_gate`), default-on wiring in `gvf-run` (Step 3.7),
and `scripts/trust_report.py`. Still to do (checklist below): the role /
evidence-type axis (§3), fail-closed unknown-gene validation (§4), per-stratum
calibration (§5), the fleet acceptance metric (§6), scorer trusted-tier
filtering, and folding the legacy carrier-guard / vf-quarantine into the unified
record.

### 1. One unified per-fact decision record (don't parallel-build)

Today quarantine is fragmented (outlier JSON flags, carrier-guard NULLs, deleted
rows). Replace with a single structure per observation:

```
trust_tier:      trusted | quarantine
trust_reasons:   [structural_rule_ids...]
trust_rule_version: <hash of the rule set>
raw_counts:      always preserved (never NULL/delete the evidence)
```

**Soft-quarantine** (keep the row, exclude it from trusted aggregates) beats
hard-delete — a fleet needs the evidence to debug and to re-tier after a rule fix.

### 2. Gene-class-agnostic structural rules (transfer across cardiac ↔ BRCA)

Prefer rules that work regardless of gene/phenotype:

| Rule class | Example | Reason id |
|---|---|---|
| Internal arithmetic | carriers ≠ affected + unaffected (when all present) | `arith_inconsistent` |
| Scope misuse | count == table grand-total / cohort N / screened N in same paper | `count_is_total` |
| Population contamination | count matches gnomAD AN/AC magnitude or annotation-table pattern (cf. PMID 33013630) | `population_count` |
| Within-paper outlier | count ≫ median × k of the paper's other variants | `paper_outlier` |
| Weak provenance | abstract-only / stub source / no table support / unlocated | `weak_provenance` |

**Avoid** cardiac-shaped heuristics that won't transfer: missense-ish priors
(BRCA is truncation/CNV/case-control heavy), fixed absolute ceilings alone,
cardiac phenotype vocabulary in rules (would false-quarantine oncology).

### 3. Evidence-type / role as a first-class axis

Segregation carrier counts and case-control "cases" must not share one semantic
field. Tag role and score **role consistency**, not just magnitude — otherwise
BRCA case-control numbers look like cardiac penetrance errors to both MAE and the
structural checks. (`pipeline/count_classifier.py` already distinguishes
case/control/cohort; promote role to the record and the rules.)

### 4. Fail-closed validation for unknown genes

`utils/variant_normalizer.py: validate_position` returns `True` for any position
when the gene has no registered length (MLH1, BRCA2 today) — fail **open**. For
the trusted tier, fail **closed**: quarantine with reason `no_gene_length` when
length/transcript is unknown, and resolve length/aliases from an external DB
(ClinVar / RefSeq) rather than hand-built per-gene files. See
`project_brca_generalization_readiness` in memory.

### 5. Calibrate per stratum, not pooled

Calibrate the accept threshold on ≥3 strata with `scripts/precision_sample.py`
**per stratum**, because a pooled Wilson CI hides BRCA failure under cardiac
volume:

1. Cardiac gold (in-distribution).
2. BRCA1/BRCA2 (truncating + case-control) — **must-include**, not cold-start.
3. One cold gene with post-hoc gold (LDLR / MLH1).

Human effort here is O(1) periodic calibration, not O(papers) adjudication.

### 6. Fleet acceptance metric for the gate itself

Gate success ≠ cardiac recall. Require, per release:

- Trusted-tier precision (adjudicated sample) ≥ target on cardiac **and** BRCA.
- Quarantine recall of known FP classes (annotation tables, gnomAD AN, study-wide N).
- Trusted-tier end-to-end count error (misses as zero) non-regressing on gold genes.
- Quarantine rate not exploding on cold genes (else the gate "works" by hiding
  everything).

### 7. Exception review routes through quarantine diffs, not per-paper UI

Variant_Browser should consume **quarantine diffs / calibration samples**, not
"review this run". Keeps human effort on the margin.

## Deferred code-review findings (fold into the next PRs)

From the `/code-review` on PR #140 (correctness/reuse not in the fleet-honesty pass):

- [x] **Overlay is per-row, not per-paper** (`cli/compare_variants.py` ~2303):
      `wrong_paper` / `excluded` drops only the matched row; missed gold variants
      on the excluded paper stay in the recall denominator. Drop **all** rows for
      that PMID. → Done: `apply_adjudication_overlay` now sweeps
      the excluded PMIDs up front and drops every row for them; regression test in
      `test_adjudication_overlay_scorer.py`.
- [x] **Overlay re-encodes the ingest contract** (`compare_variants.py` ~2190,
      ~2299): `_adjudication_variant_key` duplicates
      `ingest_review_adjudications._variant_key`, and the verdict→action branches
      duplicate `VERDICT_TO_ACTION`. Import them — drift silently drops
      adjudications. → Done: `_adjudication_variant_key`
      delegates to ingest's `_variant_key`, and `_overlay_action` resolves rows
      through the imported `VERDICT_TO_ACTION` (lazy imports dodge the load cycle).
- [ ] **End-to-end count error can't see zero-gold over-attribution**
      (`compare_variants.py` ~2487): `ComparisonRow` stores counts as
      `value or None`, so a gold `0` is indistinguishable from missing. Preserve
      raw counts (ties into the trust record's `raw_counts`).
- [ ] **Reuse dedup**: `build_comparison_rows` (a 3rd copy of the comparison
      pipeline), `_int_or_none` (→ `safe_int`), the count field-pairs constant,
      `combine_count_error_end_to_end` vs `combine_mae`, and
      `GVF_APPLY_ADJUDICATIONS` truthy-parse (→ `recall_audit.common.parse_bool`,
      so `=y` isn't silently false).
- [ ] **Latent**: `migrate_to_sqlite` raw `BEGIN` has no guard for a connection
      already in a transaction (safe with current callers). (`step_layers`
      returning a bogus `progression.csv` path on failure was FIXED — merged; it
      returns `None` now.)

## BRCA generalization prerequisites (gate depends on these)

See `project_brca_generalization_readiness` in memory. Do not replicate the 309 KB
per-gene KCNH2 alias file — handle notation **by class**:

- [ ] Register BRCA2 in `PROTEIN_LENGTHS` + `gene_metadata` (trivial; turns its
      position validation back on).
- [ ] Notation-by-class: legacy/BIC bare indels (`185delAG`, `5382insC`), IVS↔c.
      canonicalization (currently KCNH2-only), exon-level CNVs (absent).
- [ ] Silver standard from ClinVar (+ BRCA Exchange / ENIGMA) for BRCA scoring.
- [ ] A **coverage/blind-spot instrument**: count variant-like tokens the parser
      couldn't parse + unregistered-gene flag, so silent drops become visible.
      Run BRCA end-to-end and let it diagnose notation-blindness vs extraction
      error before spending effort.

## Also open (infra)

- [ ] Wire `make regression-gate` into an actual nightly/cron (needs the
      gitignored canonical DBs — self-hosted runner or local cron, not public CI).
- [ ] Consumer Gemini Code Assist sunsets 2026-07-17 — line up a replacement PR
      reviewer.
