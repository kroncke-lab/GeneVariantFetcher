# GVF Handoff Tasks

Last reviewed: 2026-08-13.

This is the only active GVF checklist. Current measurements and caveats live in
[`docs/RECALL_STATUS.md`](docs/RECALL_STATUS.md); completed benchmark history
lives in [`docs/RECALL_HISTORY.md`](docs/RECALL_HISTORY.md); protocol changes
are append-only in
[`docs/PROTOCOL_CHANGELOG.md`](docs/PROTOCOL_CHANGELOG.md). Dated reports and
benchmark run directories are evidence, not competing plans.

Execute the acceptance gates in order. Do not publish a new headline, update the
dashboard, enable a default-off recovery stage, or advance to a larger cohort
until the preceding gate passes and its artifact is recorded.

## 1. Re-establish the scientific baseline

The paired live cardiac rescore and Gate 1 review were recorded on 2026-08-13.
The lead approved advancement from Gate 1 without changing the published
headline. The requested 100--150-paper scale applies to comparison against the
already manually curated cardiac gold standard. Gate 2 is therefore the fixed,
gold-value-blinded `gold_120` sample: 30 source-available, count-eligible papers
per cardiac gene (120 attempts / 116 unique PMIDs; seed 2026081301). It does not
widen the experimental genes.

Gate 2 is complete and **passed**. The patched-system revalidation is locked in
`runs/20260813_gold120_verticalfix`: counted-extra precision is 95.70% raw and
95.87% trusted, above the 77.3% floor; the distinct count-bearing-only
diagnostic is 86.05% / 85.80%. Variant recall remains 84.09%, and carrier MAE is
0.308 raw / 0.299 trusted, still about half the 0.614 canonical all-paper
baseline. The immediately preceding stochastic run was lower at 0.266 / 0.243,
so the revalidation is a precision improvement, not a claim that every
conditional error metric improved on that one sample. The vertical-table fix
retained all 42 classification-table identities and correctly left their
carrier counts null.

The revalidation used 527 calls / 2.351M tokens and a $9.774 public-price
proxy. All four production traces are write-time verified and bound into the
prediction lock. The approved BMPR2 50 / BRCA1 50 / BRCA2 46 launch began only
after this lock and score completed.

- [ ] Finish Gate 3, `reviewer_546`: 546 attempts across 507
      unique PMIDs in the 11 populated reviewer workspaces. The first approved
      experimental strata are the existing BMPR2 50-, BRCA1 50-, and BRCA2
      46-paper queues; do not expand those to 100 papers per gene. Start with
      those three queues after the patched gold-120 revalidation. The fixed
      50/50/46 run is complete and cost/QC-locked: 972 calls, 4.261M tokens,
      $23.664 proxy, 45/50 + 50/50 + 45/46 full-text source integrity, and
      write-time-verified traces. The three fixed collaborator queues were
      refreshed in Variant Browser staging on 2026-08-13 under dataset label
      `collaborator_reextract_current_system_20260813`; pinned paper membership
      and order remained exactly 50/50/46. All 111 historical BRCA2 calls were
      preserved and moved to re-review, with zero stale calls eligible for the
      default adjudication/gold export. Public annotations were not published.
- [ ] After an accepted rescore, update `docs/RECALL_STATUS.md`, append
      `docs/RECALL_HISTORY.md` and `docs/PROTOCOL_CHANGELOG.md`, and regenerate
      the dashboard. Until then, the public dashboard remains an archived
      pre-correction snapshot.

The three manifests in
[`benchmarks/evaluation_tiers/`](benchmarks/evaluation_tiers/README.md) are the
only rollout cohorts. The curated extraction fixture and dated benchmark runs
remain specialized or historical evidence, not extra gates.

## 2. Recover missing source before adding inference

- [ ] Complete the one-time EZproxy profile bootstrap with
      `.venv/bin/python scripts/ezproxy_relogin.py --bootstrap`, then rerun
      paywall recovery and `refresh_run_db.py` for the BMPR2 stub backlog.
- [ ] Validate automatic headless EZproxy refresh against the next real session
      expiry; unit tests currently cover only the simulated browser path.
- [ ] Rebuild the two source cards that still lack decisive evidence: RYR2 PMID
      28237968 `c.13352del` and KCNH2 PMID 10862094 `c.526C>T`. Review the other
      seven untested variants in RYR2 PMID 33606749.
- [ ] Work the current source-audit ranking, beginning with blocked SCN5A
      supplements such as PMID 29325976. Use manual PDF acquisition only after
      the current audit confirms the artifact is still missing, then ingest the
      real PDF through `harvesting/format_converters.py`.
- [ ] Re-run the Karger and Sage/Liebert access probes only when credentials or
      publisher access change; the last recorded attempts remained blocked.

Source acquisition is still the dominant recall opportunity. Do not replace a
missing paper or table with speculative count arithmetic.

## 3. Measure count semantics and recovery

- [ ] Move compact count-semantics verification after all extraction and
      recovery layers so it evaluates the final winning observation rather than
      only the early `ExpertExtractor` result.
- [ ] Route only count-bearing multi-cohort or large-count ambiguities. Keep
      broad missing-slot recovery off; its 162-gap probe grounded no additions.
- [ ] Recreate clean DB copies from
      `benchmarks/codex_paper_eval/runs/20260726_fixed48_production`, run the
      hardened carrier-only count-recovery path through the current
      `run_eval.py` validation/lock workflow, and inspect rejections by gene.
- [ ] Treat RYR2 as its own acceptance stratum. Resolve curator questions for
      PMID 29622001 `c.453delC` (24 versus gold 0) and `p.L552S` (73 versus 74)
      before classifying them as recovery defects.
- [ ] Enable `COUNT_RECOVERY_ENABLED` only if the clean measurement improves
      count recall without unacceptable MAE, attribution, or RYR2 regression.
      Keep `COUNT_RECOVERY_FIELDS=carriers`; affected/unaffected widening needs
      a separate benchmark.
- [ ] Resolve the two active BRCA2 gold-scope limitations: family-versus-carrier
      semantics for PMID 26833046 and whether PMID 26848529 is exhaustive
      paper-level gold.
- [ ] Remove the known PMID 33013630 population-annotation contamination from
      canonical data through a reversible, backed-up refresh; its gnomAD and
      predictor columns are not patient carriers.

## 4. Finish the trust and fleet boundary

- [ ] Move report, publish, and other downstream count consumers to the trusted
      projection while keeping raw observations visible for audit. Report
      trusted-tier precision and quarantine rates.
- [ ] Promote count role/evidence type to a first-class persisted axis and have
      the trust gate distinguish patient, cohort-total, control, and population
      observations.
- [ ] Calibrate the trust gate separately on cardiac, BRCA, and one cold-gene
      stratum with `scripts/precision_sample.py`; do not use a pooled result to
      hide a weak stratum.
- [ ] Fail closed for unknown-gene position validation, expand BRCA notation by
      class, and build a BRCA silver standard before claiming generalization.
- [ ] Fold destructive legacy guards into the trust record so raw values remain
      recoverable, then define and schedule a self-hosted nightly regression
      gate with per-stratum acceptance metrics.
- [ ] Create or import a source-reconciled KCNE1 per-PMID gold input before
      making KCNE1 recall claims.

## 5. Engineering handoff follow-ups

- [ ] Enable run-scoped LLM tracing for standalone extraction, discovery,
      targeted landing, figure extraction, and recall-audit entry points.
- [ ] Give `scripts/smoke_azure_models.py` a normalized, scoped connectivity
      decision event.
- [ ] Move reusable metadata-backfill, dashboard trust-reader, and EZproxy
      self-heal code out of `scripts.*` so installed-package behavior matches a
      source checkout.
- [ ] Add accepted-link tests for clinical triage, extraction-priority triage,
      claim verification, final check, source-grounded summary, and synonym
      relevance.
- [ ] Perform desktop and mobile interaction/visual QA on a completed production
      trace report in an approved browser surface.
- [ ] Extend matcher coverage for indels, isoform offsets, compound variants,
      frameshifts, splice/IVS notation, and structural/CNV events; keep structural
      rows out of precision claims until they are matchable or explicitly
      excluded.
- [ ] Investigate the residual `regex_table` count-role omissions and column
      confusion without reintroducing source-free phenotype zeros.
- [ ] Establish the Friday recall rerun/compare cadence only after the new
      explicit-zero baseline is accepted.

## Deliberate decisions and non-goals

- **One checkout, one branch.** `main` is the handoff branch. Do not create local
  feature branches or experimental worktrees unless the repository owner
  explicitly changes this policy.
- **No `unassessed_count` schema.** The partial experiment was reviewed and
  retired during the 2026-08-12 consolidation. GVF records observed carrier,
  affected, unaffected, and uncertain counts; an unobserved phenotype partition
  remains SQL NULL. It must not be inferred as a residual. A future proposal
  would need end-to-end aggregation, migration, trust/final-check/guard support,
  adjudication and Variant_Browser compatibility, provenance versioning, and a
  paired benchmark before changing the schema.
- **Count recovery stays default off** until the clean carrier-only evaluation
  passes.
- **The per-paper final check and composer stay default off**; their measured
  cost/latency did not justify default-on use.
- **Historical artifacts are immutable.** Correct them with an erratum or a new
  run; do not rewrite locked predictions, historical manifests, or append-only
  recall/protocol history.
